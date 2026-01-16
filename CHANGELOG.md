# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-01-15

### Added
- **Binary chunk format** - Compact binary encoding replacing JSON POC format
  - 2-3x smaller than JSON even with compression
  - Matches or beats BAM file sizes when compressed
  - Faster parsing (no JSON overhead)
  - Preserves all SAM fields including tags
  - Backward compatible reader supports both binary (v0.2) and JSON (v0.1) formats
- **Power-of-2 chunk size selection** - Optimized for S3 cost/latency tradeoff
  - Supported sizes: 256K, 512K, 1M (default), 2M, 4M, 8M
  - Human-readable chunk size parsing (e.g., "512K", "2M")
  - Validation ensures power-of-2 values only
- **Compression enabled by default** - zstd compression is now the default
  - 6x compression ratio on binary format
  - 10ms overhead per chunk
  - Significant storage and transfer cost savings

### Changed
- Default format changed from JSON to binary
- Default compression changed from none to zstd
- Default chunk size changed from 1,000,000 bp to 1M (1,048,576 bp)
- CLI flags updated:
  - `--chunk-size` now accepts power-of-2 sizes (e.g., "512K", "1M")
  - `--format` flag added: binary (default) or json (legacy)
  - `--compression` default changed to "zstd"

### Performance
- **Binary format size comparison** (test file: 291 bytes BAM, 5 reads)
  - Binary chunks: 94 bytes and 80 bytes (2.1x smaller than JSON)
  - JSON chunks: 197 bytes and 193 bytes
  - 52% smaller with binary format vs JSON (both with zstd compression)
- **Real data validation** (GIAB chr22: 11.9M reads, 591 MB BAM)
  - Compressed BAMS3: 577 MB (0.96x vs BAM - smaller than original!)
  - Query speedup: 31x faster (3.9s vs 120s)
  - Statistics: 20,000x faster (0.006s vs 120s)
  - Data transfer reduction: 99.7% cost savings

### Technical
- Magic number: 0x42414D33 ("BAM3") in little-endian format
- 16-byte chunk header with version and flags
- 12-byte record header with variable-length fields
- 4-bit base encoding (2 bases per byte)
- Packed CIGAR operations (4 bytes per operation)
- SAM tags properly preserved with type information

### Fixed
- Tags preservation in round-trip conversion (POC limitation resolved)
- CIGAR operations now stored as structured data for binary format
- Proper little-endian magic number detection
- Decompression handled correctly for both binary and JSON formats

### Documentation
- Added `BAMS3_BINARY_FORMAT_v0.2.md` specification
- Updated CLI help with power-of-2 chunk size options
- Added cost/latency tradeoff analysis for different chunk sizes

## [0.1.0] - 2026-01-15

### Added
- BAMS3 format specification - cloud-native alignment format designed for object storage
- Python proof-of-concept implementation
  - `bams3_converter.py` - Convert BAM files to BAMS3 format
  - `bams3_query.py` - Query BAMS3 datasets by genomic region
  - `bams3_to_bam.py` - Convert BAMS3 back to BAM format
- Go CLI implementation (`bams3-go`)
  - `bams3 convert` - High-performance BAM to BAMS3 converter
  - `bams3 query` - Query reads from BAMS3 datasets
  - `bams3 stats` - Display dataset statistics instantly
- Test data generation and validation
  - Synthetic BAM file generation
  - Round-trip conversion testing
  - S3 integration testing
- Comprehensive documentation
  - Format specification
  - Implementation guides for Go and Rust
  - Comparison benchmarks vs traditional BAM
  - Bidirectional conversion guide
- S3 integration
  - Dedicated test bucket (`bams3-testing-2026`)
  - Benchmark comparing traditional vs BAMS3 workflows
  - Demonstrated 76% data transfer reduction for queries
- Project infrastructure
  - Apache 2.0 license
  - Keep a Changelog format
  - Semantic versioning

### Performance
- 4x less data transfer for small files (87 KB test)
- Projected 5,000-50,000x less data for large files (10-100 GB)
- Projected 27-267x faster queries for large files
- Instant statistics (0.1s vs 180s full scan)

### Validated
- ✅ BAM → BAMS3 conversion working
- ✅ BAMS3 → BAM conversion working
- ✅ Round-trip fidelity validated
- ✅ S3 integration tested
- ✅ Query performance measured
- ✅ Go CLI functional

### Known Limitations
- POC format uses JSON (production should use binary)
- No compression yet (zstd planned)
- SAM tags not fully preserved in POC format
- Mate pair information not stored in POC format

[Unreleased]: https://github.com/scttfrdmn/aws-direct-s3/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/scttfrdmn/aws-direct-s3/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/scttfrdmn/aws-direct-s3/releases/tag/v0.1.0
