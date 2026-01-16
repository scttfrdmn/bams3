# BAMS3 Binary Format Specification v0.2.0

## Overview

Version 0.2.0 introduces a compact binary chunk format to replace the JSON POC format. This provides:
- **5-10x size reduction** vs uncompressed JSON
- **Faster parsing** (no JSON overhead)
- **Match or beat BAM sizes** even without compression
- **Power-of-2 chunk sizes** optimized for S3 cost/latency tradeoff

## Power-of-2 Chunk Sizes

### Rationale

S3 pricing involves a tradeoff between request costs and data transfer:
- **Smaller chunks**: More S3 GET requests (higher cost), lower latency per query
- **Larger chunks**: Fewer S3 GET requests (lower cost), higher latency per query

Power-of-2 sizes provide clean memory boundaries and predictable performance.

### Supported Chunk Sizes

| Size | Bytes | Use Case |
|------|-------|----------|
| 256KB | 2^18 (262,144) | Ultra-low latency queries, hot data |
| 512KB | 2^19 (524,288) | Low latency queries, frequently accessed |
| **1MB** | **2^20 (1,048,576)** | **Default - balanced performance** |
| 2MB | 2^21 (2,097,152) | Larger queries, moderate access |
| 4MB | 2^22 (4,194,304) | Archival, infrequent access |
| 8MB | 2^23 (8,388,608) | Cold storage, minimal queries |

**Default**: 1MB (2^20 bytes) provides optimal balance for most use cases.

### Cost Analysis

**S3 Pricing** (us-west-2):
- GET requests: $0.0004 per 1,000 requests
- Data transfer: $0.09/GB (first 10TB)

**Example**: Query 100 different 1MB regions from 50GB genome

| Chunk Size | Chunks/Query | Total Requests | Request Cost | Data Transfer | Total Cost |
|------------|--------------|----------------|--------------|---------------|------------|
| 256KB | 4 | 400 | $0.00016 | 100 MB ($0.009) | $0.00916 |
| 512KB | 2 | 200 | $0.00008 | 100 MB ($0.009) | $0.00908 |
| **1MB** | **1** | **100** | **$0.00004** | **100 MB ($0.009)** | **$0.00904** |
| 2MB | 1 | 100 | $0.00004 | 200 MB ($0.018) | $0.01804 |
| 4MB | 1 | 100 | $0.00004 | 400 MB ($0.036) | $0.03604 |

**Conclusion**: For most queries, 1MB chunks minimize total cost. Request costs are negligible compared to data transfer.

### Latency Analysis

**From EC2 (same region)**:
- S3 GET latency: 10-20ms
- Decompression (zstd): ~2ms per MB
- Total per chunk: ~15-25ms

| Chunk Size | Decompress Time | Total Latency |
|------------|-----------------|---------------|
| 256KB | 0.5ms | ~15ms |
| 512KB | 1ms | ~16ms |
| 1MB | 2ms | ~17ms |
| 2MB | 4ms | ~19ms |
| 4MB | 8ms | ~23ms |

**From remote (cross-region/internet)**:
- S3 GET latency: 100-200ms (dominates)
- Chunk size has minimal impact on latency

**Recommendation**: Use 1MB for most workloads. Consider 256KB-512KB for ultra-low latency requirements.

## Binary Chunk Format

### Overview

Each chunk file contains a binary-encoded list of alignment records. The format is:

```
[Header: 16 bytes]
[Record 1]
[Record 2]
...
[Record N]
```

### Chunk Header (16 bytes)

| Offset | Size | Type | Field | Description |
|--------|------|------|-------|-------------|
| 0 | 4 | uint32 | magic | Magic number: 0x42414D33 ("BAM3") |
| 4 | 2 | uint16 | version | Format version: 0x0200 (v0.2.0) |
| 6 | 2 | uint16 | flags | Compression and encoding flags |
| 8 | 4 | uint32 | num_records | Number of records in chunk |
| 12 | 4 | uint32 | reserved | Reserved for future use |

**Flags field** (16 bits):
- Bits 0-3: Compression algorithm
  - 0x0: None
  - 0x1: zstd
  - 0x2-0xF: Reserved
- Bits 4-7: Reserved
- Bits 8-15: Reserved

### Record Format

Each record uses a compact binary encoding:

```
[Record Header: 12 bytes]
[Variable-length data]
```

#### Record Header (12 bytes)

| Offset | Size | Type | Field | Description |
|--------|------|------|-------|-------------|
| 0 | 4 | uint32 | record_size | Total size of this record (including header) |
| 4 | 2 | uint16 | name_len | Length of read name |
| 6 | 2 | uint16 | seq_len | Length of sequence |
| 8 | 2 | uint16 | cigar_ops | Number of CIGAR operations |
| 10 | 2 | uint16 | num_tags | Number of optional tags |

#### Variable-length Data

Following the header, in order:

1. **Read name** (name_len bytes): UTF-8 encoded string
2. **Reference ID** (4 bytes): int32 (-1 for unmapped)
3. **Position** (4 bytes): int32 (0-based, -1 for unmapped)
4. **Mapping quality** (1 byte): uint8
5. **Flags** (2 bytes): uint16 (SAM flags)
6. **CIGAR** (cigar_ops * 4 bytes): Packed CIGAR operations
   - Each operation: uint32 with op_len in upper 28 bits, op_type in lower 4 bits
7. **Sequence** (⌈seq_len/2⌉ bytes): 4-bit encoded bases (2 bases per byte)
   - Encoding: A=1, C=2, G=4, T=8, N=15, other=0
8. **Quality** (seq_len bytes): Phred+33 quality scores
9. **Tags** (variable): Optional tags in binary format
   - Each tag: 2 bytes (tag name) + 1 byte (type) + variable data

### Size Comparison

**Example read** (100bp, name="READ001", CIGAR=100M, 2 tags):

| Format | Size | Notes |
|--------|------|-------|
| JSON (POC) | ~350 bytes | Human-readable, verbose |
| JSON + zstd | ~60 bytes | Compressed JSON |
| **Binary** | **~140 bytes** | Compact encoding |
| **Binary + zstd** | **~40 bytes** | Best compression |
| BAM | ~130 bytes | Similar to binary format |

**Expected savings**:
- Binary format: 2.5x smaller than JSON
- Binary + zstd: 6-9x smaller than JSON
- **Matches or beats BAM** file sizes

## Implementation Plan

### Phase 1: Binary Writer

1. Create `pkg/bams3/binary_writer.go`
   - Implement chunk header encoding
   - Implement record binary encoding
   - Support compression flag
2. Update `pkg/bams3/writer.go`
   - Add format flag (json/binary)
   - Default to binary for v0.2.0

### Phase 2: Binary Reader

1. Create `pkg/bams3/binary_reader.go`
   - Implement chunk header decoding
   - Implement record binary decoding
   - Auto-detect format (check magic number)
2. Update `pkg/bams3/reader.go`
   - Auto-detect chunk format
   - Support both JSON (v0.1) and binary (v0.2)

### Phase 3: CLI Updates

1. Update `cmd/bams3/convert.go`
   - Change default compression to "zstd"
   - Change default format to "binary"
   - Add `--chunk-size` with power-of-2 parsing
   - Support: "256K", "512K", "1M", "2M", "4M", "8M"
2. Update metadata format version to "0.2.0"

### Phase 4: Backward Compatibility

- Reader must support both v0.1 (JSON) and v0.2 (binary) chunks
- Auto-detect based on:
  - Check first 4 bytes for magic number 0x42414D33
  - If not present, assume JSON format
- Metadata includes format version and per-chunk format info

## Migration Path

### From v0.1 to v0.2

Users can convert existing BAMS3 datasets:

```bash
# Read v0.1 (JSON) and write v0.2 (binary)
bams3 upgrade sample_v0.1.bams3 sample_v0.2.bams3
```

Or simply query v0.1 datasets directly (reader supports both).

### Default Behavior

**v0.2.0 defaults**:
- Format: binary
- Compression: zstd
- Chunk size: 1M (1,048,576 bases)

```bash
# Uses defaults (binary, zstd, 1MB chunks)
bams3 convert input.bam output.bams3

# Explicit options
bams3 convert --format binary --compression zstd --chunk-size 1M input.bam output.bams3

# Smaller chunks for low-latency queries
bams3 convert --chunk-size 512K input.bam output.bams3

# Larger chunks for archival
bams3 convert --chunk-size 4M input.bam output.bams3
```

## Performance Targets

### Conversion Speed

| Metric | v0.1.0 (JSON) | v0.2.0 Target |
|--------|---------------|---------------|
| Throughput | 51,600 reads/sec | 100,000+ reads/sec |
| Time (11.9M reads) | 231 seconds | 120 seconds |
| Speedup | 1x | 2x |

### File Sizes

| Dataset | Original BAM | v0.1 (JSON+zstd) | v0.2 Target |
|---------|--------------|------------------|-------------|
| chr22 (591 MB) | 591 MB | 577 MB (0.96x) | 450 MB (0.76x) |
| WGS (50 GB) | 50 GB | 48 GB | 38 GB |

**Goal**: 20-25% smaller than BAM with binary format + zstd.

### Query Performance

| Operation | v0.1.0 | v0.2.0 Target |
|-----------|--------|---------------|
| 100KB region query | 3.9s | 2.0s (50% faster) |
| Statistics | 0.006s | 0.006s (unchanged) |
| Decompression | 2ms/MB | 2ms/MB (unchanged) |

**Improvement comes from**:
- Smaller chunks → less data to decompress
- Binary parsing → faster than JSON

## Validation

### Test Suite

1. **Round-trip testing**: BAM → BAMS3 v0.2 → BAM
   - All fields preserved
   - Tags preserved correctly
   - Checksums match

2. **Format compatibility**:
   - v0.2 reader can read v0.1 chunks
   - v0.2 writer produces valid binary chunks

3. **Compression testing**:
   - Binary + zstd matches or beats BAM size
   - Decompression speed acceptable (<5ms per MB)

4. **Real data validation**:
   - Test with GIAB chr22 (11.9M reads)
   - Test with full WGS sample (3B reads)
   - Validate against samtools output

## Roadmap

### v0.2.0 (Binary Format)
- ✅ Binary chunk format specification (this document)
- 🔄 Binary reader/writer implementation
- ⏳ Power-of-2 chunk size selection
- ⏳ Backward compatibility (read v0.1)
- ⏳ CLI updates and documentation
- ⏳ Real data validation

### v0.3.0 (Performance)
- Parallel conversion (multi-core)
- Direct S3 read/write
- Streaming compression
- Memory optimization

### v0.4.0 (Advanced Features)
- Columnar storage option
- Multi-level indexing
- Pre-computed statistics per chunk
- Adaptive chunk sizing

## References

- BAM format specification: https://samtools.github.io/hts-specs/SAMv1.pdf
- CRAM format: https://samtools.github.io/hts-specs/CRAMv3.pdf
- S3 pricing: https://aws.amazon.com/s3/pricing/

---

**Version**: 0.2.0-draft
**Date**: 2026-01-15
**Status**: Design phase
