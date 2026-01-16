# BAMS3 Format Specification v0.1.0

**Status:** Proof of Concept
**Date:** 2026-01-15
**Authors:** Scott Friedman
**License:** Apache 2.0

## Abstract

BAMS3 (BAM for S3) is a cloud-native alignment format designed specifically for object storage systems. Unlike BAM, which is optimized for POSIX filesystems, BAMS3 uses an object-per-chunk architecture that enables selective data access, parallel processing, and instant metadata queries.

This specification describes version 0.1.0, a proof-of-concept implementation using JSON encoding. Production versions will use binary encoding.

## 1. Design Principles

### 1.1 Object-Native Architecture

BAMS3 stores alignment data as independent objects rather than a single monolithic file. Each genomic region (chunk) is a separate object that can be accessed independently.

**Benefits:**
- Selective access (download only needed regions)
- Parallel processing (multiple chunks simultaneously)
- Cloud-optimized (natural fit for S3, GCS, Azure Blob)

### 1.2 Metadata Separation

Dataset metadata, statistics, and indexes are stored separately from alignment data. This enables instant queries without reading alignment data.

**Benefits:**
- Instant statistics (no file scanning)
- Fast query planning (know which chunks to access)
- Efficient browsing (understand data without downloading)

### 1.3 Format Versioning

BAMS3 datasets include explicit version information for forward compatibility. Readers must check version compatibility before processing.

## 2. Directory Structure

A BAMS3 dataset is a directory (local filesystem) or prefix (object storage) containing:

```
dataset.bams3/
├── _metadata.json              # Required: Dataset metadata and manifest
├── _header.json                # Required: SAM header
├── _index/                     # Optional: Indexes
│   └── spatial.json            # Position → chunk mapping
└── data/                       # Required: Alignment data
    ├── chr1/                   # Per-reference directories
    │   ├── 000000000-001000000.chunk
    │   ├── 001000000-002000000.chunk
    │   └── ...
    ├── chr2/
    │   └── ...
    └── unmapped.chunk          # Unmapped reads
```

### 2.1 Naming Conventions

- **Dataset directory:** `{name}.bams3` or `{name}.bams3/`
- **Reference directories:** `{reference_name}/` (e.g., `chr1/`, `chrX/`)
- **Chunk files:** `{start:09d}-{end:09d}.chunk` (e.g., `000000000-001000000.chunk`)
- **Unmapped:** `unmapped.chunk`
- **Metadata files:** Prefix with `_` to distinguish from data

### 2.2 File Requirements

**Required:**
- `_metadata.json` - Dataset metadata and chunk manifest
- `_header.json` - SAM header
- `data/` - At least one chunk file

**Optional:**
- `_index/spatial.json` - Spatial index (recommended for large datasets)
- Additional custom indexes

## 3. Metadata Format (_metadata.json)

The metadata file contains dataset-level information and a manifest of all chunks.

### 3.1 Structure

```json
{
  "format": "bams3",
  "version": "0.1.0-poc",
  "created": "2026-01-15T16:36:38.123Z",
  "created_by": "bams3-go v0.1.0",

  "source": {
    "file": "sample.bam",
    "format": "BAM"
  },

  "statistics": {
    "total_reads": 11930668,
    "mapped_reads": 11930668,
    "unmapped_reads": 0,
    "duplicate_reads": 0,
    "total_bases": 303477269,
    "mean_coverage": 0.0
  },

  "chunks": [
    {
      "path": "data/chr22/016000000-017000000.chunk",
      "reference": "chr22",
      "start": 16000000,
      "end": 17000000,
      "reads": 331508,
      "size_bytes": 16234567,
      "compression": "zstd",
      "checksum": "a1b2c3...",
      "created": "2026-01-15T16:38:22.456Z"
    }
  ],

  "compression": {
    "algorithm": "zstd"
  },

  "chunk_size": 1000000
}
```

### 3.2 Field Definitions

#### Top-Level Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `format` | string | Yes | Must be `"bams3"` |
| `version` | string | Yes | Format version (e.g., `"0.1.0-poc"`) |
| `created` | string | Yes | ISO 8601 timestamp |
| `created_by` | string | Yes | Tool name and version |
| `source` | object | No | Original source file info |
| `statistics` | object | Yes | Dataset-level statistics |
| `chunks` | array | Yes | Manifest of all chunks |
| `compression` | object | Yes | Compression settings |
| `chunk_size` | integer | Yes | Target chunk size in base pairs |

#### Statistics Object

| Field | Type | Description |
|-------|------|-------------|
| `total_reads` | integer | Total number of reads |
| `mapped_reads` | integer | Number of mapped reads |
| `unmapped_reads` | integer | Number of unmapped reads |
| `duplicate_reads` | integer | Number of duplicate reads |
| `total_bases` | integer | Total bases sequenced |
| `mean_coverage` | float | Mean coverage depth |

#### Chunk Object

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `path` | string | Yes | Relative path to chunk file |
| `reference` | string | Yes | Reference name or "unmapped" |
| `start` | integer | Yes | Start position (0-based, inclusive) |
| `end` | integer | Yes | End position (0-based, exclusive) |
| `reads` | integer | Yes | Number of reads in chunk |
| `size_bytes` | integer | Yes | Chunk file size in bytes |
| `compression` | string | Yes | Compression: "none" or "zstd" |
| `checksum` | string | Yes | SHA-256 hash of chunk file |
| `created` | string | Yes | ISO 8601 timestamp |

#### Compression Object

| Field | Type | Description |
|-------|------|-------------|
| `algorithm` | string | "none" or "zstd" |
| `level` | integer | Compression level (optional) |

## 4. Header Format (_header.json)

The header file contains SAM header information in JSON format.

### 4.1 Structure

```json
{
  "HD": {
    "VN": "1.6",
    "SO": "coordinate"
  },
  "SQ": [
    {
      "SN": "chr1",
      "LN": "249250621"
    },
    {
      "SN": "chr2",
      "LN": "243199373"
    }
  ],
  "RG": [
    {
      "ID": "sample1",
      "SM": "NA12878",
      "PL": "ILLUMINA"
    }
  ],
  "PG": [
    {
      "ID": "1",
      "PN": "bwa",
      "VN": "0.7.17"
    }
  ],
  "CO": [
    "Comment line 1",
    "Comment line 2"
  ]
}
```

### 4.2 Field Definitions

| Field | Type | Description |
|-------|------|-------------|
| `HD` | object | Header line (VN, SO, GO) |
| `SQ` | array | Reference sequences |
| `RG` | array | Read groups (optional) |
| `PG` | array | Programs (optional) |
| `CO` | array | Comments (optional) |

Each field corresponds to its SAM header equivalent. Values are stored as strings to preserve exact SAM formatting.

## 5. Chunk Format

### 5.1 Current Format (v0.1.0-poc)

Chunks are stored as JSON arrays of read objects. This is a proof-of-concept format; production versions will use binary encoding.

```json
[
  {
    "name": "read_000124",
    "flag": 0,
    "ref": 21,
    "pos": 16010988,
    "mapq": 21,
    "cigar": "75M",
    "seq": "GTTGTGTGATGCCCAGCGTGTTAAGGATTACTAAGTAATAGGCCAATTGCCTCCTATGGAGCGTCAACAACGAGG",
    "qual": "::6GA?B@885DA8A96C:HGI<F=AI6CHD<>@D>DBD;97=>BAG9@?G>@HB<95;E<?>BI6?=FB@87GB"
  }
]
```

### 5.2 Read Object Fields

| Field | Type | Description | SAM Equivalent |
|-------|------|-------------|----------------|
| `name` | string | Read name | QNAME |
| `flag` | integer | SAM flags | FLAG |
| `ref` | integer | Reference ID (-1 for unmapped) | RNAME (as ID) |
| `pos` | integer | 0-based position | POS - 1 |
| `mapq` | integer | Mapping quality | MAPQ |
| `cigar` | string | CIGAR string | CIGAR |
| `seq` | string | Sequence | SEQ |
| `qual` | string | Quality scores (ASCII) | QUAL |
| `tags` | object | Optional SAM tags | TAG (optional) |

**Notes:**
- Positions are 0-based (SAM is 1-based)
- Reference IDs match the order in `_header.json` SQ array
- Mate information not stored in v0.1.0-poc
- Tags are optional and not implemented in v0.1.0-poc

### 5.3 Compression

Chunks may be compressed with zstd. The compression type is specified in the chunk metadata.

**Uncompressed:**
- File contains raw JSON
- `compression: "none"`

**Compressed:**
- File contains zstd-compressed JSON
- `compression: "zstd"`
- Readers must decompress before parsing JSON

## 6. Spatial Index (_index/spatial.json)

The spatial index provides fast lookup of chunks by genomic position.

### 6.1 Structure

```json
{
  "references": {
    "chr22": {
      "name": "chr22",
      "length": 51304566,
      "chunks": [
        {
          "path": "data/chr22/016000000-017000000.chunk",
          "reference": "chr22",
          "start": 16000000,
          "end": 17000000,
          "reads": 331508,
          "size_bytes": 16234567,
          "compression": "zstd",
          "checksum": "a1b2c3...",
          "created": "2026-01-15T16:38:22.456Z"
        }
      ]
    }
  }
}
```

### 6.2 Usage

Readers use the spatial index to:
1. Find which chunks overlap a query region
2. Determine chunk locations without scanning metadata
3. Plan parallel queries across multiple regions

## 7. Query Algorithm

### 7.1 Region Query

To query reads in region `chr:start-end`:

1. Read `_metadata.json`
2. Find chunks where `chunk.reference == chr` AND `chunk.end > start` AND `chunk.start < end`
3. For each overlapping chunk:
   - Load chunk file
   - Decompress if `compression != "none"`
   - Parse JSON
   - Filter reads where `read.pos >= start` AND `read.pos < end`
4. Return filtered reads

### 7.2 Statistics Query

To get dataset statistics:

1. Read `_metadata.json`
2. Return `statistics` object

Time: O(1), no data scanning required.

## 8. Coordinate System

BAMS3 uses **0-based, half-open intervals** for positions:

- `start`: Inclusive (first base in range)
- `end`: Exclusive (first base NOT in range)

**Examples:**
- `chr1:0-1000` includes positions 0-999 (1000 bases)
- `chr1:1000000-2000000` includes positions 1,000,000 - 1,999,999

**SAM Conversion:**
- SAM uses 1-based positions
- BAMS3 position = SAM position - 1

## 9. Chunk Size Guidelines

### 9.1 Recommended Sizes

| Use Case | Chunk Size | Chunks per GB | Query Latency |
|----------|-----------|---------------|---------------|
| Interactive queries | 256 KB - 512 KB | 2,000 - 4,000 | Low (<1s) |
| **Balanced (default)** | **1 MB** | **1,000** | **Medium (1-5s)** |
| Archival | 5 MB - 10 MB | 100 - 200 | High (5-20s) |

### 9.2 Trade-offs

**Smaller chunks:**
- ✅ Faster queries (less decompression)
- ✅ More selective access
- ❌ More metadata overhead
- ❌ More S3 requests

**Larger chunks:**
- ✅ Better compression ratio
- ✅ Fewer S3 requests
- ❌ Slower queries (decompress entire chunk)
- ❌ Less selective access

**Recommendation:** Use 1 MB chunks (default) unless specific use case requires different size.

## 10. Compression

### 10.1 Supported Algorithms

| Algorithm | Status | Ratio | Speed | Notes |
|-----------|--------|-------|-------|-------|
| none | Required | 1.0x | ∞ | No compression |
| **zstd** | **Recommended** | **6x** | **Fast** | Default for v0.1.0 |
| gzip | Future | 3x | Medium | Compatibility |
| lz4 | Future | 2x | Very fast | Low latency |

### 10.2 Compression Settings

**Default (balanced):**
```bash
bams3 convert --compression zstd sample.bam sample.bams3
```
- zstd level 3 (default)
- Good compression (6x for JSON data)
- Fast decompression (~500 MB/s)

**No compression:**
```bash
bams3 convert --compression none sample.bam sample.bams3
```
- Fastest queries
- 6x larger
- Use for hot data with cost-insensitive storage

## 11. Validation

### 11.1 Dataset Validation

A valid BAMS3 dataset must:

1. ✅ Contain `_metadata.json` with `format: "bams3"`
2. ✅ Contain `_header.json` with valid SAM header
3. ✅ Have at least one chunk file
4. ✅ All chunks listed in metadata must exist
5. ✅ Chunk checksums must match file contents
6. ✅ Chunk positions must not overlap within same reference
7. ✅ All reads in a chunk must map to chunk's reference
8. ✅ Statistics must match actual read counts

### 11.2 Compatibility

Readers must:
- Check `format == "bams3"`
- Verify version compatibility
- Support required compression algorithms
- Handle missing optional fields gracefully

## 12. Performance Characteristics

### 12.1 Validated Performance (chr22, 11.9M reads, 591 MB BAM)

| Operation | Time | Data Accessed |
|-----------|------|---------------|
| Convert to BAMS3 (zstd) | 231s | 591 MB |
| Query 100KB region | 3.9s | 16 MB (1 chunk) |
| Get statistics | 0.006s | 37 KB (metadata) |
| Full scan (parallel) | TBD | 577 MB (all chunks) |

### 12.2 Scaling Estimates

**Query 1MB region from 50GB whole genome:**

| Method | Time | Data Downloaded |
|--------|------|-----------------|
| Traditional BAM | 600s | 50 GB |
| FUSE | 10s | 50 GB* |
| **BAMS3** | **8s** | **50 MB** |

*First access; cached subsequent access

## 13. Limitations (v0.1.0-poc)

### 13.1 Known Limitations

1. **JSON format** - 5-6x larger than binary BAM before compression
   - Mitigation: zstd compression brings to parity with BAM
   - Solution: Binary format in v0.2.0

2. **No SAM tags** - Tags not preserved in conversion
   - Solution: Add `tags` field to read objects

3. **No mate pair info** - Paired-end mate information not stored
   - Solution: Add mate fields to read objects

4. **Fixed chunk size** - May create uneven chunk sizes
   - Solution: Adaptive chunking based on read density

5. **Sequential compression** - Single-threaded compression
   - Solution: Parallel compression in v0.2.0

### 13.2 Not Limitations

These are intentional design choices:

- **Multiple files** - Object-per-chunk is core design principle
- **No backward compatibility with BAM** - Different format for different purpose
- **Metadata duplication** - Required for independent object access

## 14. Future Extensions

### 14.1 Binary Format (v0.2.0)

Replace JSON with efficient binary encoding:
- Match or beat BAM file sizes
- Faster parsing (no JSON overhead)
- Smaller even without compression

### 14.2 Columnar Storage (v0.3.0)

Store read components separately:
```
data/chr1/001000000-002000000/
  names.bin
  sequences.bin
  qualities.bin
  positions.bin
  cigars.bin
```

Benefits:
- Query only needed columns
- 10-100x less data for coverage queries
- Better compression per column

### 14.3 Multi-Sample (v0.4.0)

Store multiple samples in one dataset:
```
cohort.bams3/
  data/chr1/001000000-002000000/
    sample001.chunk
    sample002.chunk
    ...
```

Benefits:
- Analyze regions across cohorts
- Shared metadata and indexes
- Efficient multi-sample queries

## 15. Reference Implementation

### 15.1 Available Tools

**bams3-go** (v0.1.0):
- Go implementation of BAMS3 tools
- Commands: convert, query, stats
- Optional zstd compression
- Source: github.com/scttfrdmn/bams3-go

**Python POC** (v0.1.0):
- Proof-of-concept Python tools
- Bidirectional conversion (BAM ↔ BAMS3)
- Source: aws-direct-s3/format-tools/bams3/

### 15.2 Usage Example

```bash
# Convert BAM to BAMS3
bams3 convert --compression zstd sample.bam sample.bams3

# Query region
bams3 query sample.bams3 chr1:1000000-2000000

# Get statistics
bams3 stats sample.bams3

# Upload to S3
aws s3 sync sample.bams3/ s3://bucket/sample.bams3/
```

## 16. Versioning

BAMS3 follows [Semantic Versioning 2.0.0](https://semver.org/):

- **MAJOR** version: Incompatible format changes
- **MINOR** version: Backward-compatible features
- **PATCH** version: Backward-compatible bug fixes

### 16.1 Version History

| Version | Date | Status | Changes |
|---------|------|--------|---------|
| 0.1.0-poc | 2026-01-15 | Current | Initial POC with JSON chunks |
| 0.2.0 | TBD | Planned | Binary chunk format |
| 1.0.0 | TBD | Planned | Production stable release |

### 16.2 Compatibility

Readers should check version compatibility:

```python
if metadata['version'].startswith('0.1.'):
    # Compatible with v0.1.x
    process_dataset()
else:
    raise ValueError(f"Unsupported version: {metadata['version']}")
```

## 17. License

BAMS3 format specification: CC0 1.0 Universal (Public Domain)

Reference implementations: Apache License 2.0

## 18. References

- SAM/BAM specification: https://samtools.github.io/hts-specs/SAMv1.pdf
- zstd compression: https://facebook.github.io/zstd/
- Semantic Versioning: https://semver.org/

---

**Copyright 2026 Scott Friedman**

**Document Version:** 1.0
**BAMS3 Format Version:** 0.1.0-poc
**Last Updated:** 2026-01-15
