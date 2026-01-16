# BAMS3 Format Specification

**Status:** Production Ready (v0.2.0)
**Date:** 2026-01-15
**Authors:** Scott Friedman
**License:** Apache 2.0

## Overview

BAMS3 (BAM for S3) is a cloud-native alignment format designed specifically for object storage systems. Unlike BAM, which is optimized for POSIX filesystems, BAMS3 uses an object-per-chunk architecture that enables selective data access, parallel processing, and instant metadata queries.

## Version History

### v0.2.0 - Binary Format (Production) ✅ **RECOMMENDED**
**Status:** Production Ready
**Date:** 2026-01-15
**Format:** Compact binary encoding with zstd compression

**Key Features:**
- 10% smaller than BAM files
- 3.7x faster queries than JSON format
- 92x faster than traditional BAM workflows
- All SAM fields preserved including tags
- Power-of-2 chunk sizes (256K, 512K, 1M, 2M, 4M, 8M)
- Backward compatible reader (auto-detects format)

**Documentation:** See [BAMS3_BINARY_FORMAT_v0.2.md](./BAMS3_BINARY_FORMAT_v0.2.md)
**Validation:** See [BINARY_FORMAT_VALIDATION.md](./bams3-go/BINARY_FORMAT_VALIDATION.md)

### v0.1.0 - JSON Format (Proof of Concept)
**Status:** Deprecated (Legacy Support Only)
**Date:** 2026-01-15
**Format:** JSON encoding with optional zstd compression

**Key Features:**
- Human-readable format for debugging
- 2% smaller than BAM with compression
- 25x faster queries than traditional workflows
- Simple implementation for prototyping

**Documentation:** See [BAMS3_SPECIFICATION_v0.1.md](./BAMS3_SPECIFICATION_v0.1.md)

**⚠️ Note:** JSON format is deprecated for production use. Use binary format (v0.2.0) instead.

## Quick Start

### Converting BAM to BAMS3

```bash
# Production (binary format, recommended)
bams3 convert sample.bam sample.bams3

# With custom chunk size
bams3 convert --chunk-size 512K sample.bam sample.bams3

# Legacy JSON format
bams3 convert --format json sample.bam sample.bams3
```

### Querying BAMS3

```bash
# Query a region
bams3 query sample.bams3 chr1:1000000-2000000

# Count reads only
bams3 query sample.bams3 chr1:1000000-2000000 --count

# Get dataset statistics (instant!)
bams3 stats sample.bams3
```

### Converting Back to BAM

```bash
bams3 to-bam sample.bams3 reconstructed.bam
```

## Design Principles

### 1. Object-Native Architecture

BAMS3 stores alignment data as independent objects rather than a single monolithic file. Each genomic region (chunk) is a separate object that can be accessed independently.

**Benefits:**
- **Selective access:** Download only needed regions
- **Parallel processing:** Process multiple chunks simultaneously
- **Cloud-optimized:** Natural fit for S3, GCS, Azure Blob
- **Cost efficient:** Pay only for data you access

**Example:** Query 1MB region from 50GB genome
- Traditional: Download all 50GB ($4.50)
- BAMS3: Download 16MB chunk ($0.001)
- **Savings:** 99.98% reduction

### 2. Metadata Separation

Dataset metadata, statistics, and indexes are stored separately from alignment data. This enables instant queries without reading alignment data.

**Benefits:**
- **Instant statistics:** No file scanning (0.006s vs 120s)
- **Fast query planning:** Know which chunks to access
- **Efficient browsing:** Understand data without downloading
- **Version detection:** Auto-detect format version

**Metadata files:**
- `_metadata.json` - Statistics, chunk manifest, version
- `_header.json` - SAM header (references, read groups, etc.)
- `_index/spatial.json` - Position → chunk mapping

### 3. Format Versioning

BAMS3 datasets include explicit version information for forward compatibility. Readers automatically detect and handle both v0.1 (JSON) and v0.2 (binary) formats.

**Version detection:**
- Check metadata file for version string ("0.1.0" or "0.2.0")
- Check chunk format (binary magic number or JSON)
- Automatic fallback to appropriate parser

## Directory Structure

A BAMS3 dataset is a directory (local filesystem) or prefix (object storage) containing:

```
dataset.bams3/
├── _metadata.json              # Required: Dataset metadata and manifest
├── _header.json                # Required: SAM header
├── _index/                     # Optional: Indexes
│   └── spatial.json            # Position → chunk mapping
└── data/                       # Required: Alignment data
    ├── chr1/                   # Per-reference directories
    │   ├── 000000000-001048576.chunk    # Binary format (1M chunks)
    │   ├── 001048576-002097152.chunk
    │   └── ...
    ├── chr2/
    │   └── ...
    └── unmapped.chunk          # Unmapped reads (if any)
```

**Chunk naming conventions:**
- Power-of-2 boundaries: `000000000-001048576.chunk` (1M = 2^20 bytes)
- Start position is inclusive, end is exclusive
- 9-digit zero-padded coordinates
- Extension: `.chunk` for both binary and JSON formats

## Chunk Sizes

### Power-of-2 Sizes (v0.2.0)

BAMS3 v0.2.0 uses power-of-2 chunk sizes optimized for S3 access patterns:

| Size | Bytes | Use Case | Latency | Cost |
|------|-------|----------|---------|------|
| 256K | 262,144 | Ultra-low latency | Lowest | Higher |
| 512K | 524,288 | Low latency queries | Low | Medium |
| **1M** | **1,048,576** | **Default (recommended)** | **Medium** | **Optimal** |
| 2M | 2,097,152 | Larger queries | Medium | Lower |
| 4M | 4,194,304 | Archival | Higher | Lowest |
| 8M | 8,388,608 | Cold storage | Highest | Lowest |

**Tradeoff:** Smaller chunks mean more S3 GET requests (higher cost, lower latency). Larger chunks mean fewer requests (lower cost, higher latency).

**Default:** 1M provides optimal balance for most use cases.

### Legacy Sizes (v0.1.0)

JSON format used base-10 sizes (1,000,000 bp by default).

## Format Comparison

| Feature | v0.1.0 (JSON) | v0.2.0 (Binary) |
|---------|---------------|-----------------|
| **Status** | Deprecated | ✅ Production |
| **Encoding** | JSON text | Compact binary |
| **Size vs BAM** | 0.98x (2% smaller) | 0.90x (10% smaller) |
| **Compression** | Optional zstd | zstd by default |
| **Query Speed** | 25x faster | 92x faster |
| **Parse Speed** | Slow (JSON) | Fast (binary) |
| **Human Readable** | Yes | No |
| **Tags Preserved** | Limited | Full |
| **Chunk Sizes** | Base-10 | Power-of-2 |
| **Magic Number** | None | 0x42414D33 |

## Performance

### Real Data Results (GIAB chr22: 11.9M reads, 591 MB)

| Format | Size | Conversion | Query (100KB) | Statistics |
|--------|------|------------|---------------|------------|
| Original BAM | 591 MB | - | 120s (download) | 120s (scan) |
| JSON v0.1 | 577 MB | 231s | 4.8s | 0.006s |
| **Binary v0.2** | **531 MB** | **316s** | **1.3s** | **0.006s** |

**Key Achievements:**
- ✅ 10% smaller than BAM
- ✅ 92x faster queries
- ✅ 20,000x faster statistics
- ✅ 99.7% cost savings on data transfer

### Conversion Throughput

| Format | Reads/sec | Notes |
|--------|-----------|-------|
| JSON v0.1 | 51,600 | Simple JSON encoding |
| Binary v0.2 | 37,700 | More complex binary encoding |
| **Target v0.3** | **100,000+** | With parallel compression |

## Use Cases

### When to Use BAMS3

✅ **Cloud-native workflows**
- Data stored in S3, GCS, Azure Blob
- Need selective region access
- Cost-sensitive workloads
- Distributed processing

✅ **Frequent queries**
- >10 queries per dataset
- Interactive analysis
- Real-time variant calling
- Coverage visualization

✅ **Large-scale analysis**
- Population studies (1000+ samples)
- Multi-sample cohorts
- Continuous data accumulation
- Long-term archival

### When to Use Traditional BAM

❌ **Local-only processing**
- Never stored in cloud
- Full-file processing always needed
- Legacy tool compatibility required

❌ **One-time use**
- Single analysis run
- Immediate deletion after use
- No repeated access needed

## Implementation

### Current Implementations

1. **bams3-go** (Official)
   - Language: Go
   - Status: Production ready
   - CLI: convert, query, stats, to-bam
   - Format: Binary v0.2.0 (default), JSON v0.1.0 (legacy)

2. **Python POC** (Reference)
   - Language: Python
   - Status: Proof of concept only
   - Format: JSON v0.1.0 only
   - Use: Educational and testing

### Future Implementations

- Rust library (high-performance, zero-copy)
- C/C++ library (legacy tool integration)
- Python library (production, wraps Go/Rust)
- Java library (Spark/Hadoop integration)

## Roadmap

### v0.3.0 (Planned)
- Parallel compression (3-5x faster conversion)
- Direct S3 reading (no local download)
- Streaming queries (process without loading)
- Memory optimization

### v0.4.0 (Planned)
- Columnar storage option
- Multi-level indexing
- Pre-computed statistics per chunk
- Adaptive chunk sizing

### v1.0.0 (Future)
- Multi-sample datasets
- Variant calling integration
- Expression quantification format
- Community adoption

## Tools and Libraries

### Official Tools

```bash
# Install bams3-go CLI
go install github.com/scttfrdmn/bams3-go/cmd/bams3@latest

# Or build from source
git clone https://github.com/scttfrdmn/bams3-go.git
cd bams3-go
go build ./cmd/bams3
```

### Usage Examples

```bash
# Convert BAM to BAMS3 (binary format, 1M chunks, zstd compression)
bams3 convert sample.bam sample.bams3

# Query specific region
bams3 query sample.bams3 chr1:1000000-2000000

# Get instant statistics
bams3 stats sample.bams3

# Convert back to BAM
bams3 to-bam sample.bams3 reconstructed.bam

# Custom chunk size for low-latency queries
bams3 convert --chunk-size 512K sample.bam sample_512k.bams3

# No compression (faster conversion, larger size)
bams3 convert --compression none sample.bam sample_uncompressed.bams3

# Legacy JSON format
bams3 convert --format json sample.bam sample_json.bams3
```

## Cost Analysis

### Storage Costs (S3 Standard: $0.023/GB/month)

| Dataset | BAM | BAMS3 v0.2 | Monthly Savings |
|---------|-----|------------|-----------------|
| 591 MB sample | $0.014 | $0.012 | $0.002 |
| 50 GB WGS | $1.15 | $1.04 | $0.11 |
| 1 TB cohort | $23.00 | $20.71 | $2.29 |
| 100 TB project | $2,300 | $2,071 | $229 |

**Annual savings (100 TB):** $2,748

### Data Transfer Costs (First 10TB: $0.09/GB)

**Scenario:** Query 10 different 1MB regions from 50GB genome, 100 times

| Method | Data Transfer | Cost | Savings |
|--------|---------------|------|---------|
| Traditional | 500 TB (50GB × 10,000) | $45,000 | - |
| FUSE (cached) | 5 TB (50GB × 100) | $450 | 90% |
| **BAMS3** | **16 GB** (160MB × 100) | **$1.44** | **99.997%** |

**Annual savings (typical research lab):** $44,998.56

## Migration Guide

### From BAM to BAMS3

1. **Convert your BAM files:**
   ```bash
   bams3 convert --chunk-size 1M sample.bam sample.bams3
   ```

2. **Upload to S3:**
   ```bash
   aws s3 sync sample.bams3 s3://my-bucket/datasets/sample.bams3/
   ```

3. **Update your workflows:**
   - Replace `samtools view` with `bams3 query`
   - Replace `samtools stats` with `bams3 stats`
   - Keep BAM files for tools that require them

4. **Gradual adoption:**
   - BAMS3 for cloud storage and queries
   - Convert back to BAM for legacy tools

### Backward Compatibility

BAMS3 v0.2.0 reader automatically detects and reads both formats:
- JSON v0.1.0 datasets (legacy)
- Binary v0.2.0 datasets (current)

No user intervention required - format detection is automatic.

## Community and Support

### Getting Help

- Documentation: [github.com/scttfrdmn/bams3-go](https://github.com/scttfrdmn/bams3-go)
- Issues: [github.com/scttfrdmn/bams3-go/issues](https://github.com/scttfrdmn/bams3-go/issues)
- Discussions: [github.com/scttfrdmn/bams3-go/discussions](https://github.com/scttfrdmn/bams3-go/discussions)

### Contributing

Contributions welcome! Areas of interest:
- Additional implementations (Rust, C++, Java)
- Cloud provider integrations
- Performance optimizations
- Documentation improvements
- Real-world use cases

## References

- **SAM/BAM Specification:** https://samtools.github.io/hts-specs/SAMv1.pdf
- **zstd Compression:** https://facebook.github.io/zstd/
- **AWS Open Data Registry:** https://registry.opendata.aws/

## License

Apache License 2.0

Copyright 2026 Scott Friedman

See [LICENSE](./LICENSE) file for details.

---

**Specification Version:** 0.2.0
**Last Updated:** 2026-01-15
**Status:** Production Ready ✅
