# AWS Direct S3 Access for Research Applications

> **⚡ NEW:** BAMS3 format is now **working and tested!** See [test-data/TEST_RESULTS.md](test-data/TEST_RESULTS.md) for validation results showing **10-150x speedup** for queries. Try it: `cd test-data && cat README.md`

## The Problem

Research data increasingly lives in S3 (genomics sequencing data, imaging data, etc.), but most scientific tools are POSIX-only applications that expect local filesystem paths. This creates a painful workflow:

1. **Data generated** → lands in S3 (e.g., sequencer output)
2. **Users copy data locally** → `aws s3 cp` for hundreds of GB or TB
3. **Run tools** → BWA, GATK, samtools, etc.
4. **Upload results** → back to S3
5. **Delete local copies** → to free expensive local storage

**Problems with this approach:**
- Massive time spent copying data
- Requires expensive local storage (NVMe, large EBS volumes)
- Data duplication and version confusion
- Users often don't know there are alternatives

### The Origin Story: Why This Project Exists

**The specific problem that sparked BAMS3:**

We were pulling genomics data from AWS Open Data (us-east-1) → copying to our own S3 (us-west-2) → then copying again to local compute. **Two completely unnecessary copies, costing $4+ per 100GB file.**

**The realization:** Instead of copying data around, just launch a $1.53/hour EC2 instance in the same region as the data and stream what you need.

**The math:**
- Cross-region copy: $2.00 + 30 minutes
- Download to local: $9.00 + 30 minutes
- **Total waste:** $11.00 and 60 minutes before analysis even starts

**BAMS3 approach:**
- Launch c5.9xlarge in source region: $1.53/hour
- Stream only what you need: $0.0045 per gene
- Start analysis immediately: 0 wait time
- **Actual cost:** ~$0.26 for most queries (10 min × $1.53/hour)

**The insight:** The EC2 instance costs less than the transfer you avoided, and you get your answer 10× faster. See [docs/COST_OPTIMIZATION.md](docs/COST_OPTIMIZATION.md) for detailed cost analysis.

## Use Case: Genomics Pipeline Example

```bash
# Current workflow (SLOW and WASTEFUL)
aws s3 cp s3://my-bucket/sample1_R1.fastq.gz ./data/     # 50GB, takes 10+ minutes
aws s3 cp s3://my-bucket/sample1_R2.fastq.gz ./data/     # 50GB
bwa mem -t 16 reference.fa data/sample1_R1.fastq.gz data/sample1_R2.fastq.gz > aligned.sam
samtools view -b aligned.sam > aligned.bam
aws s3 cp aligned.bam s3://my-bucket/results/            # Upload back
rm -rf data/                                              # Clean up
```

## Quick Start

**Want to test immediately?** See [QUICKSTART.md](QUICKSTART.md) for a 5-minute setup guide.

**Want to benchmark your workflows?** Use our benchmarking tools to compare approaches:

```bash
# Quick benchmark: compare all methods
cd scripts
./quick-benchmark.sh s3://your-bucket/sample.bam samtools flagstat

# Example output:
# copy-then-process:    145.23s  (baseline)
# fuse-mountpoint:       52.18s  (2.78x faster) ✓
# streaming:             48.91s  (2.97x faster) ✓
#
# Recommendation: Direct S3 access provides significant benefit!
```

See [scripts/README.md](scripts/README.md) for detailed benchmarking guide.

## Three Fundamental Approaches

There are three ways to solve the S3 access problem:

1. **Adapt S3 to look like POSIX** (workarounds: FUSE, streaming)
2. **Modify tools to support S3** (add S3 support to tools)
3. **Design formats for S3** (rethink data formats entirely) ⭐

### Approach A: Workarounds (Quick Start)

**Make S3 appear as local filesystem**

#### 1. FUSE Filesystem Mounting
Mount S3 as a POSIX filesystem, tools see it as local files.

**Options:**
- **s3fs-fuse** - Mature, widely used, but slower
- **goofys** - Faster than s3fs-fuse, optimized for read
- **mountpoint-for-amazon-s3** - AWS official, high performance

**Pros:** Zero code changes to applications
**Cons:** Performance varies, caching strategies matter

#### 2. Streaming/Pipes
For tools that accept stdin/stdout, stream directly from S3.

```bash
# Example: stream data instead of copying
aws s3 cp s3://bucket/file.fastq.gz - | gunzip | tool -
```

**Pros:** No local storage needed
**Cons:** Only works for sequential read tools

### Approach B: Modify Tools (Long-term Solution)

**Add native S3 support to existing tools**

Modify open source genomics tools to natively support `s3://` URIs.

**Key targets:**
- **htslib** - Foundation library for samtools, bcftools, etc. Already has partial S3 support via plugins
  - Enhance existing S3 support
  - Add range request optimization
  - Improve error handling and retry logic

- **BWA/minimap2** - Read aligners
  - Accept S3 URIs for reference genomes
  - Stream FASTQ input from S3

- **GATK** - Already Java-based, could leverage AWS SDK
  - Add native S3 support alongside current Google Cloud support

**Benefits:**
- Changes propagate to entire ecosystem
- Can be contributed upstream
- Users just use `s3://` URIs instead of paths
- Tools can optimize for S3 (range requests, prefetching)

**Example target usage:**
```bash
# Future ideal workflow - no copying!
samtools view s3://my-bucket/sample1.bam chr1:1000000-2000000
bwa mem s3://refs/hg38.fa s3://data/sample1_R1.fastq.gz s3://data/sample1_R2.fastq.gz | \
  samtools view -b - > aligned.bam
aws s3 cp aligned.bam s3://my-bucket/results/
```

See [docs/tool-modification-guide.md](docs/tool-modification-guide.md) for details.

### Approach C: Redesign Formats ⭐ (Most Transformative)

**Design new formats native to object storage**

Instead of "how do we make BAM work on S3?", ask "what would we design if S3 was the primary storage?"

#### Example: BAMS3 (BAM for S3)

**Traditional BAM:**
```
sample.bam          (single 500GB file)
sample.bam.bai      (index file)
```

**BAMS3:**
```
sample.bams3/
  _metadata.json              # Dataset info, statistics
  _index/spatial.idx          # Position → chunk mapping
  data/
    chr1/
      00000000-01000000.chunk  # 1MB independent chunks
      01000000-02000000.chunk
      ...
    chr2/
      ...
```

**Key benefits:**
- **Object-per-chunk**: Each chunk is independently accessible
- **Parallel processing**: Process 32 chromosomes simultaneously
- **Minimal data transfer**: Query 1MB region → download 1-2MB (not 500GB!)
- **No external index**: Metadata embedded in dataset
- **Modern compression**: zstd/lz4 instead of BGZF

**Performance:**

| Operation | BAM (copy+use) | BAM (FUSE) | BAMS3 |
|-----------|----------------|------------|-------|
| Query 1MB region | 125s | 12s | **0.8s** |
| Full scan | 180s | 250s | **95s** (parallel) |
| Get statistics | 180s (scan) | 180s | **0.1s** (metadata) |

**Usage:**
```bash
# Convert BAM → BAMS3
bams3 convert input.bam s3://bucket/output.bams3

# Query (only downloads needed chunks)
bams3 query s3://bucket/sample.bams3 chr1:1000000-2000000

# Statistics (instant - reads metadata only)
bams3 stats s3://bucket/sample.bams3
```

**More format redesigns:**
- **VCF → Columnar** (Parquet): Query by column, predicate pushdown
- **FASTQ → Structured** (Arrow): Random access, metadata queries
- **Multi-sample cohorts**: One dataset for 1000 samples by region

See:
- [format-tools/bams3-spec.md](format-tools/bams3-spec.md) - Complete BAMS3 specification
- [format-tools/bams3/README.md](format-tools/bams3/README.md) - Proof-of-concept implementation
- [docs/format-design-for-s3.md](docs/format-design-for-s3.md) - General format design principles

**This is the future**: Formats designed for cloud from the ground up, not retrofitted.

## Project Structure

```
.
├── README.md                          # This file
├── docs/
│   ├── problem-analysis.md           # Detailed problem breakdown
│   ├── fuse-comparison.md            # FUSE options comparison
│   ├── tool-modification-guide.md    # How to add S3 support to tools
│   ├── htslib-s3-status.md          # Current state of htslib S3 support
│   ├── format-optimization.md        # Cloud-optimized format strategies
│   ├── genomics-tools.md             # Tool-specific notes
│   └── benchmarks.md                 # Performance comparisons
├── examples/
│   ├── fuse-mounting/                # FUSE mount examples
│   ├── streaming/                    # Streaming approaches
│   ├── hybrid/                       # Combined strategies
│   ├── genomics-pipeline/            # Real pipeline examples
│   └── format-conversion/            # BAM → cloud-optimized formats
├── tool-modifications/
│   ├── htslib-patches/               # Patches for htslib S3 improvements
│   ├── bwa-s3/                       # BWA with S3 support
│   ├── minimap2-s3/                  # minimap2 with S3 support
│   └── contrib-guide.md              # How to contribute upstream
├── format-tools/
│   ├── bam-chunker/                  # Split BAMs for parallel access
│   ├── vcf-to-parquet/               # Convert VCF to columnar format
│   ├── cloud-optimize/               # General optimization tool
│   └── validators/                   # Check if files are cloud-ready
└── scripts/
    ├── setup-mountpoint.sh           # Install and configure
    ├── benchmark.py                  # Performance testing
    └── cache-analysis.py             # Cache hit rate analysis
```

## Goals

### Short-term (Workarounds)
1. **Document** real-world trade-offs of FUSE and streaming approaches
2. **Benchmark** performance with genomics-sized data
3. **Provide** practical examples for naive users to avoid unnecessary copying

### Long-term (Sustainable Solutions)
4. **Modify** key genomics tools to natively support S3
5. **Contribute** improvements upstream to htslib, BWA, etc.
6. **Develop** format optimization tools and best practices
7. **Create** cloud-optimized data layouts for common genomics workflows
8. **Build** conversion and validation tools for cloud-ready formats

## Success Metrics

- **User experience:** `tool s3://bucket/file.bam` just works
- **Performance:** S3-native access competitive with or better than copy-then-process
- **Community impact:** Patches accepted into upstream projects
- **Adoption:** Researchers stop copying data unnecessarily

## Getting Started

See individual approach directories for setup instructions and examples.
