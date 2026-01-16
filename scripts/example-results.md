# Example Benchmark Results

This document shows example benchmark results for different file sizes, operations, and access patterns to help you understand what performance to expect.

## Test Setup

- **Location:** EC2 c5.4xlarge in us-east-1
- **S3 Bucket:** Same region (us-east-1)
- **Network:** 10 Gbps
- **Cache:** 100GB NVMe local storage

## Test 1: Medium BAM File (12 GB)

**File:** `s3://bucket/sample-12gb.bam`
**Operation:** `samtools flagstat` (sequential read, full file scan)

| Method | Duration | Throughput | Speedup | Notes |
|--------|----------|------------|---------|-------|
| copy-then-process | 145.2s | 0.08 GB/s | 1.0x (baseline) | Download: 120s, Process: 25s |
| fuse-mountpoint | 52.1s | 0.23 GB/s | 2.8x | Cache warm: ~35s on second run |
| fuse-goofys | 68.5s | 0.18 GB/s | 2.1x | Slower than mountpoint |
| fuse-s3fs | 112.3s | 0.11 GB/s | 1.3x | Not recommended |
| streaming | 48.9s | 0.25 GB/s | **3.0x** | **Winner** |
| direct-s3 | N/A | N/A | N/A | samtools not compiled with S3 |

**Recommendation:** Use streaming for this operation.

**Key insight:** Download dominates time. Streaming processes while downloading, saving ~100s.

## Test 2: Small BAM File (500 MB)

**File:** `s3://bucket/sample-500mb.bam`
**Operation:** `samtools flagstat`

| Method | Duration | Throughput | Speedup | Notes |
|--------|----------|------------|---------|-------|
| copy-then-process | 6.2s | 0.08 GB/s | 1.0x | Download: 5.5s, Process: 0.7s |
| fuse-mountpoint | 5.8s | 0.09 GB/s | 1.1x | Mount overhead ~0.5s |
| streaming | 5.5s | 0.09 GB/s | **1.1x** | Marginal benefit |

**Recommendation:** Copy-then-process is fine for small files.

**Key insight:** For small files, overhead dominates and benefits are minimal.

## Test 3: Large BAM File (150 GB)

**File:** `s3://bucket/sample-150gb.bam`
**Operation:** `samtools flagstat`

| Method | Duration | Throughput | Speedup | Notes |
|--------|----------|------------|---------|-------|
| copy-then-process | 1847s (31min) | 0.08 GB/s | 1.0x | Download: 1800s, Process: 47s |
| fuse-mountpoint | 682s (11min) | 0.22 GB/s | 2.7x | First access |
| streaming | 621s (10min) | 0.24 GB/s | **3.0x** | **Winner** |

**Recommendation:** Always use direct S3 access for large files!

**Key insight:** Saved 20 minutes. For 100 samples = 33 hours saved.

## Test 4: Random Access with Index

**File:** `s3://bucket/sample-12gb.bam` + `.bai`
**Operation:** `samtools view chr1:1000000-2000000` (1 MB region)

| Method | Duration | Throughput | Speedup | Notes |
|--------|----------|------------|---------|-------|
| copy-then-process | 122.5s | N/A | 1.0x | Must download entire 12GB! |
| fuse-mountpoint | 1.8s | N/A | **68x** | **Winner** - Only fetches needed chunks |
| streaming | N/A | N/A | N/A | Can't seek in stream |
| direct-s3 | 2.1s | N/A | 58x | If tool supports (samtools with S3) |

**Recommendation:** FUSE or direct S3 essential for indexed queries.

**Key insight:** Avoid downloading entire file for small region queries!

## Test 5: Multiple Region Queries

**File:** `s3://bucket/sample-12gb.bam`
**Operation:** 100 random region queries (each ~1MB)

| Method | Duration | Speedup | Notes |
|--------|----------|---------|-------|
| copy-then-process | 125s | 1.0x | Download once: 120s, Query: 5s |
| fuse-mountpoint (cold cache) | 185s | 0.7x | Each query fetches from S3 |
| fuse-mountpoint (warm cache) | 8s | **16x** | **Winner** - All data cached |
| direct-s3 | 92s | 1.4x | Repeated S3 requests |

**Recommendation:**
- One-time queries: Copy-then-process
- Repeated queries: FUSE with cache (best after first access)

**Key insight:** Caching matters! First query warms cache, subsequent queries are fast.

## Test 6: FASTQ Quality Control

**File:** `s3://bucket/sample-50gb_R1.fastq.gz`
**Operation:** `fastqc` (sequential read)

| Method | Duration | Throughput | Speedup | Notes |
|--------|----------|------------|---------|-------|
| copy-then-process | 612s (10min) | 0.08 GB/s | 1.0x | Download: 600s, Process: 12s |
| fuse-mountpoint | 198s (3min) | 0.25 GB/s | 3.1x | |
| streaming | 185s (3min) | 0.27 GB/s | **3.3x** | **Winner** |

**Recommendation:** Stream FASTQ for QC - huge time savings.

## Test 7: VCF Query (Columnar Format)

**File:** `s3://bucket/variants.parquet` (10 GB, 100M variants)
**Operation:** Query chromosome 1 only (10% of data)
**Tool:** Python with pyarrow

| Method | Duration | Data Read | Speedup | Notes |
|--------|----------|-----------|---------|-------|
| copy-then-process | 128s | 10 GB | 1.0x | Download: 120s, Query: 8s |
| direct-s3 (row format) | 125s | 10 GB | 1.0x | Must read entire file |
| direct-s3 (parquet) | 18s | 1 GB | **7.1x** | **Winner** - Columnar + predicate pushdown |

**Recommendation:** Use columnar formats (Parquet) for selective queries.

**Key insight:** Format matters! Parquet reads only needed data.

## Test 8: Repeated Access Pattern

**Workflow:** Analyze same 10 BAM files (5GB each), 5 different operations each

| Method | Total Time | Notes |
|--------|------------|-------|
| copy-then-process | 3200s (53min) | Download once: 600s, Process: 2600s |
| fuse-mountpoint | 2100s (35min) | First access slow, cached after |
| streaming (each time) | 2850s (48min) | Re-downloads for each operation |

**Recommendation:**
- Multiple operations on same files: FUSE with cache or copy once
- Single operation per file: Streaming

## Summary: When to Use Each Method

### Copy-Then-Process
**Use when:**
- ✅ Files are small (<1GB)
- ✅ You'll access the file many times
- ✅ You have cheap/fast local storage
- ✅ Other methods don't work

**Avoid when:**
- ❌ Files are large (>10GB)
- ❌ One-time sequential processing
- ❌ Storage is limited/expensive

### FUSE Mounting
**Use when:**
- ✅ Random access needed (indexed queries)
- ✅ Tools don't support streaming
- ✅ Repeated access to same files
- ✅ Want transparent S3 access

**Avoid when:**
- ❌ Only sequential access
- ❌ One-time processing
- ❌ Can't install FUSE tools

**Best with:**
- Large cache on fast storage
- Files accessed multiple times
- Operations that benefit from caching

### Streaming
**Use when:**
- ✅ Sequential read only
- ✅ Large files (>5GB)
- ✅ One-time processing
- ✅ Storage is limited

**Avoid when:**
- ❌ Need random access / seeking
- ❌ Tool doesn't support stdin
- ❌ Multiple passes over data

**Best for:**
- FASTQ quality control
- Format conversion
- Filtering pipelines
- One-off analyses

### Direct S3 (Native Tool Support)
**Use when:**
- ✅ Tool supports s3:// URIs
- ✅ Any file size
- ✅ Any access pattern

**Currently available:**
- samtools/bcftools (if compiled with S3)
- Many Python tools (pandas, pyarrow, etc.)
- Growing support in other tools

**Long-term goal:**
- This is the ideal - tools should just support S3 natively!

## Cost Considerations

Beyond time savings, consider S3 costs:

### S3 Pricing (us-east-1, approximate)
- **Storage:** $0.023/GB/month
- **GET requests:** $0.0004 per 1,000 requests
- **Data transfer OUT:** $0.09/GB (to internet), FREE (within region)

### Cost Examples

**Scenario: Process 100 files (10GB each) once**

| Method | S3 Requests | Data Transfer | S3 Cost |
|--------|-------------|---------------|---------|
| Copy-then-process | ~10,000 | 1000 GB | $90.04 |
| FUSE (cold cache) | ~100,000 | 1000 GB | $90.40 |
| Streaming | ~10,000 | 1000 GB | $90.04 |

**Insight:** S3 data transfer is the major cost. Direct access doesn't add significant cost.

**Within AWS (EC2 in same region):** Data transfer is FREE! Direct S3 access costs almost nothing.

## Recommendations by Use Case

### Use Case 1: Daily Pipeline (100s of samples)
- **Use:** FUSE mounting with large cache
- **Why:** Amortize download over multiple operations
- **Expected:** 2-3x speedup after cache warming

### Use Case 2: QC Before Analysis Decision
- **Use:** Streaming
- **Why:** Need quick QC to decide if full analysis is worth it
- **Expected:** 3x speedup, avoid full download

### Use Case 3: Interactive Analysis
- **Use:** FUSE mounting
- **Why:** Random access, repeated queries
- **Expected:** 10-100x speedup for region queries

### Use Case 4: One-Time Big Analysis
- **Use:** Copy-then-process or streaming
- **Why:** Will access data many times, copy once
- **Expected:** Time saved on processing offsets copy time

## Your Turn

Run benchmarks on your own data:

```bash
cd scripts
./quick-benchmark.sh s3://your-bucket/your-file.bam samtools flagstat
```

Results will vary based on:
- File size
- Network speed
- Operation type
- Access pattern
- Cache configuration

Document your results and choose the best method for your workflows!
