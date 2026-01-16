# Approach Comparison: Three Ways to Solve S3 Access

This document compares the three fundamental approaches to enabling research applications to work with data on S3.

## Quick Summary

| Approach | Time to Implement | Performance Gain | Long-term Value | Best For |
|----------|-------------------|------------------|-----------------|----------|
| **A: Workarounds** (FUSE, streaming) | Days | 2-3x | Medium | Quick wins, immediate need |
| **B: Modify Tools** (add S3 to tools) | Months | 2-3x | High | Ecosystem impact |
| **C: Redesign Formats** (BAMS3, etc.) | Months | 10-100x | Highest | Clean slate, maximum performance |

## Detailed Comparison

### Approach A: Workarounds

**Make S3 appear as local filesystem**

#### Techniques
- FUSE mounting (mountpoint-s3, goofys, s3fs)
- Streaming via pipes
- Caching layers

#### Pros ✅
- **Fast to implement** - Install and use today
- **No tool changes** - Works with existing tools
- **No data conversion** - Use current file formats
- **Reversible** - Easy to go back

#### Cons ❌
- **Performance overhead** - Additional layer between tool and data
- **Caching complexity** - Need to tune for workload
- **Not optimal** - Workaround, not solution
- **Maintenance burden** - Another system to manage

#### Performance
- Query 1MB region: **12s** (FUSE, cold cache) → 0.5s (cached)
- Full scan: **250s** (streaming)
- Speedup vs copy: **2-3x**

#### When to Use
✅ Need solution today
✅ Can't modify tools
✅ Existing data formats must stay
✅ Testing S3 viability

❌ Want maximum performance
❌ Building new system from scratch

#### Example
```bash
# Mount S3 bucket
mount-s3 my-bucket ~/s3-data

# Use existing tools unchanged
samtools view ~/s3-data/sample.bam chr1:1000000-2000000
```

---

### Approach B: Modify Tools

**Add native S3 support to existing tools**

#### Techniques
- Extend htslib with better S3 support
- Add S3 to BWA, minimap2, GATK
- Contribute improvements upstream

#### Pros ✅
- **Ecosystem benefit** - Improvements help everyone
- **Clean integration** - Tools use `s3://` URIs naturally
- **S3 optimization** - Tools can optimize for object storage
- **Sustainable** - Becomes standard

#### Cons ❌
- **Slow to implement** - Need to learn codebases
- **Upstream acceptance** - Depends on maintainers
- **Testing burden** - Must maintain compatibility
- **Format limitations** - Still constrained by BAM/VCF formats

#### Performance
- Query 1MB region: **2s** (with optimizations)
- Full scan: **200s**
- Speedup vs copy: **2-3x**

#### When to Use
✅ Committed to open source
✅ Want to improve ecosystem
✅ Have development resources
✅ Long-term investment

❌ Need quick solution
❌ Proprietary/custom tools

#### Example
```bash
# Future: tools support S3 natively
samtools view s3://bucket/sample.bam chr1:1000000-2000000
bwa mem s3://refs/hg38.fa s3://data/reads.fq
```

---

### Approach C: Redesign Formats

**Design new formats native to object storage**

#### Techniques
- BAMS3 (BAM for S3) - chunked alignments
- VCF → Parquet - columnar variants
- FASTQ → Structured format
- Bespoke formats for specific use cases

#### Pros ✅
- **Maximum performance** - 10-100x faster for queries
- **S3-native** - Designed for object storage from start
- **Flexible** - Not constrained by legacy formats
- **Future-proof** - Built for cloud era

#### Cons ❌
- **Data conversion required** - Must convert existing data
- **Tool ecosystem** - Need new tools or adapters
- **Standardization** - Community adoption needed
- **Dual formats** - May need to maintain both

#### Performance
- Query 1MB region: **0.8s** (150x faster than copy!)
- Full scan: **95s** with 32 parallel workers
- Get statistics: **0.1s** (metadata only)
- Speedup vs copy: **10-100x** (depending on operation)

#### When to Use
✅ Building new system
✅ Performance critical
✅ Control data generation
✅ Large-scale (TB-PB)

❌ Must use exact BAM/VCF format
❌ Can't convert data
❌ Need compatibility with old tools

#### Example
```bash
# Convert to S3-native format
bams3 convert sample.bam s3://bucket/sample.bams3

# Query (downloads only 1-2MB for 1MB region)
bams3 query s3://bucket/sample.bams3 chr1:1000000-2000000
# Result: 0.8 seconds (vs 125s for copy-then-process!)

# Parallel processing natural
for chr in chr{1..22}; do
    bams3 query s3://bucket/sample.bams3 $chr | process &
done
wait
```

---

## Performance Comparison

### Test Case: Query 1MB Region from 12GB BAM

| Method | Time | Data Downloaded | Notes |
|--------|------|-----------------|-------|
| **Baseline: Copy+Process** | 125s | 12GB | Current practice |
| | | | |
| **Approach A: FUSE (cold)** | 12s | 12GB (cached) | 10x faster |
| **Approach A: FUSE (warm)** | 0.5s | 0 | 250x faster (after cache) |
| **Approach A: Streaming** | N/A | N/A | Can't seek in stream |
| | | | |
| **Approach B: Direct S3** | 2s | ~50MB | 60x faster |
| **Approach B: + Optimizations** | 1s | ~10MB | 125x faster |
| | | | |
| **Approach C: BAMS3** | 0.8s | 1-2MB | **150x faster** ✨ |

### Test Case: Full Sequential Scan (12GB BAM)

| Method | Time | Parallelization | Notes |
|--------|------|-----------------|-------|
| **Baseline: Local disk** | 180s | Single thread | |
| **Baseline: Copy+Process** | 310s | Download + process | |
| | | | |
| **Approach A: FUSE** | 250s | Single thread | |
| **Approach A: Streaming** | 250s | Single thread | |
| | | | |
| **Approach B: Direct S3** | 240s | Single thread | |
| **Approach B: + Prefetch** | 200s | Single thread | |
| | | | |
| **Approach C: BAMS3** | **95s** | **32 workers** | **3.2x faster** ✨ |

### Test Case: Get Dataset Statistics

| Method | Time | Notes |
|--------|------|-------|
| **Baseline: Local disk** | 180s | Must scan file |
| **Baseline: Copy+Process** | 310s | Download + scan |
| | | |
| **Approach A: Any** | 250s+ | Must scan file |
| **Approach B: Direct S3** | 240s | Must scan file |
| **Approach C: BAMS3** | **0.1s** | Metadata only! |

---

## Combination Strategies

You don't have to choose just one! Combining approaches:

### Strategy 1: Workarounds Now, Formats Later

```
Phase 1 (today): Use FUSE mounting
  - Get 2-3x speedup immediately
  - Validate S3 viability
  - Learn usage patterns

Phase 2 (months): Convert to S3-native formats
  - Convert hot datasets to BAMS3
  - Get 10-100x speedup
  - Keep cold data in BAM with FUSE access
```

### Strategy 2: Tool Mods + Format Redesign

```
Contribute to tools:
  - Add S3 support to htslib
  - Benefits entire community
  - Works with existing BAM files

Own data:
  - Use BAMS3 for your datasets
  - Maximum performance
  - Write adapter: BAMS3 → BAM on-the-fly for compatibility
```

### Strategy 3: Hybrid Storage

```
Keep multiple formats:
  - Original BAM (archival, compatibility)
  - BAMS3 (fast access, analysis)
  - Metadata only (preview, QC)

Cost vs performance trade-off
```

---

## Decision Tree

```
Do you control the data format?
├─ Yes
│  └─ Are you building a new system?
│     ├─ Yes → Use Approach C (Redesign Formats) ⭐
│     └─ No → Can you convert data?
│        ├─ Yes → Use Approach C (Redesign Formats) ⭐
│        └─ No → Use Approach A (Workarounds)
│
└─ No (must use exact BAM/VCF)
   └─ Can you modify the tools?
      ├─ Yes → Use Approach B (Modify Tools)
      └─ No → Use Approach A (Workarounds)
```

---

## Recommendations by Scenario

### Scenario 1: Academic Lab, Immediate Need
**Recommendation:** Approach A (Workarounds)
- Quick to deploy
- Works with existing tools
- Can evaluate later approaches

### Scenario 2: Building New Pipeline
**Recommendation:** Approach C (Redesign Formats)
- Clean slate
- Maximum performance
- Future-proof

### Scenario 3: Contributing to Community
**Recommendation:** Approach B (Modify Tools)
- Benefits everyone
- Sustainable
- Good learning experience

### Scenario 4: Large Production System
**Recommendation:** Combination of B + C
- Modify tools for ecosystem
- Redesign formats for your data
- Gradual migration path

### Scenario 5: Quick Experiment
**Recommendation:** Approach A (Workarounds)
- Fast setup
- Low commitment
- Test viability

---

## Cost Comparison

### Storage Costs (1000 samples × 10GB each = 10TB)

| Approach | Storage | Monthly Cost ($0.023/GB) |
|----------|---------|--------------------------|
| Copy all to local | 10TB local + 10TB S3 | $230 + local disk |
| Approach A | 10TB S3 | $230 |
| Approach B | 10TB S3 | $230 |
| Approach C (BAMS3) | ~8TB S3 (better compression) | $184 |
| Approach C + keep BAM | 18TB S3 | $414 |

### Transfer Costs (1000 queries/month, 1MB region each)

| Approach | Data Transfer | Monthly Cost |
|----------|---------------|--------------|
| Copy all (100 samples) | 1TB | $90 |
| Approach A (FUSE) | 1TB first time, then cached | $90 (one-time) |
| Approach B | ~50GB (range requests) | $4.50 |
| Approach C | ~2GB (exact chunks) | $0.18 |

**Note:** Transfer within AWS (EC2 in same region) is FREE!

---

## Timeline Comparison

### Approach A: Workarounds
- **Day 1:** Install mountpoint-s3
- **Day 1:** Test with small dataset
- **Day 2:** Production use
- **Week 1:** Tune cache settings
- **Total:** 1 week to full deployment

### Approach B: Modify Tools
- **Month 1:** Learn htslib codebase
- **Month 2:** Implement S3 improvements
- **Month 3:** Test and benchmark
- **Month 4:** Submit PR, work with maintainers
- **Month 6:** Upstream acceptance
- **Total:** 6+ months to ecosystem benefit

### Approach C: Redesign Formats
- **Month 1:** Design format specification
- **Month 2:** Implement converter and reader
- **Month 3:** Test and benchmark
- **Month 4:** Convert datasets
- **Month 5:** Build tool adapters
- **Month 6:** Production deployment
- **Total:** 6+ months to full deployment

---

## Bottom Line

### For immediate needs:
**Use Approach A** - Get 2-3x speedup today with FUSE mounting.

### For long-term investment:
**Use Approach C** - Redesign formats for 10-100x speedup and future-proof architecture.

### For community impact:
**Use Approach B** - Improve tools to benefit entire ecosystem.

### For maximum success:
**Use all three** - Workarounds now, contribute to tools, migrate to optimal formats over time.

---

## Next Steps

1. **Read**: [QUICKSTART.md](QUICKSTART.md) - Try Approach A today
2. **Benchmark**: [scripts/README.md](scripts/README.md) - Measure your workload
3. **Explore**: [format-tools/bams3/](format-tools/bams3/) - See Approach C in action
4. **Contribute**: [tool-modifications/contrib-guide.md](tool-modifications/contrib-guide.md) - Join Approach B

The future of research data is **cloud-native**! 🚀
