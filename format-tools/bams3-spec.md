# BAMS3 Format Specification

**BAMS3:** BAM for S3 - A cloud-native alignment format

## Overview

BAMS3 is a replacement for BAM designed specifically for object storage like S3. Instead of a single monolithic file, BAMS3 is a **collection of objects** organized for efficient cloud access.

**Key differences from BAM:**
- ❌ Single file → ✅ Multiple objects (one per chunk)
- ❌ External index → ✅ Embedded metadata
- ❌ BGZF compression → ✅ Modern compression (zstd, lz4)
- ❌ Sequential scan → ✅ Parallel access
- ❌ 64KB blocks → ✅ Larger chunks (1-10MB)

## Directory Structure

```
sample.bams3/
├── _metadata.json              # Global metadata, schema, manifest
├── _header.json                # SAM header (references, read groups, etc.)
├── _index/
│   ├── spatial.idx             # Position → chunk mapping
│   └── stats.json              # Pre-computed statistics
└── data/
    ├── chr1/
    │   ├── 00000000-01000000.chunk     # 1M base pair chunks
    │   ├── 01000000-02000000.chunk
    │   └── ...
    ├── chr2/
    │   └── ...
    └── unmapped.chunk
```

## File Specifications

### _metadata.json

```json
{
  "format": "bams3",
  "version": "1.0.0",
  "created": "2024-01-15T10:30:00Z",
  "created_by": "bam-to-bams3 v1.0",

  "sample": {
    "name": "NA12878",
    "library": "Illumina-HiSeq",
    "platform": "ILLUMINA",
    "center": "BI"
  },

  "reference": {
    "name": "hg38",
    "accession": "GCA_000001405.15"
  },

  "statistics": {
    "total_reads": 450000000,
    "mapped_reads": 445000000,
    "unmapped_reads": 5000000,
    "duplicate_reads": 45000000,
    "total_bases": 67500000000,
    "mean_coverage": 30.2
  },

  "chunks": [
    {
      "path": "data/chr1/00000000-01000000.chunk",
      "reference": "chr1",
      "start": 0,
      "end": 1000000,
      "reads": 12534,
      "size_bytes": 1048576,
      "compression": "zstd",
      "checksum": "sha256:abc123...",
      "created": "2024-01-15T10:30:15Z"
    },
    {
      "path": "data/chr1/01000000-02000000.chunk",
      "reference": "chr1",
      "start": 1000000,
      "end": 2000000,
      "reads": 12891,
      "size_bytes": 1102345,
      "compression": "zstd",
      "checksum": "sha256:def456...",
      "created": "2024-01-15T10:30:22Z"
    }
    // ... more chunks
  ],

  "compression": {
    "algorithm": "zstd",
    "level": 3
  },

  "chunk_size": 1000000,  // Base pairs per chunk

  "compatibility": {
    "can_convert_to_bam": true,
    "bam_source": "s3://original-data/sample.bam",
    "conversion_time": "2024-01-15T10:30:00Z"
  }
}
```

### _header.json

Standard SAM header in JSON format:

```json
{
  "HD": {
    "VN": "1.6",
    "SO": "coordinate"
  },
  "SQ": [
    {
      "SN": "chr1",
      "LN": 248956422
    },
    {
      "SN": "chr2",
      "LN": 242193529
    }
    // ... all references
  ],
  "RG": [
    {
      "ID": "rg1",
      "SM": "NA12878",
      "LB": "lib1",
      "PL": "ILLUMINA"
    }
  ],
  "PG": [
    {
      "ID": "bwa",
      "VN": "0.7.17",
      "CL": "bwa mem ref.fa reads.fq"
    }
  ]
}
```

### Chunk Format (.chunk files)

Binary format for each chunk:

```
CHUNK FILE STRUCTURE
====================

[HEADER: 128 bytes fixed]
  Magic: "BAS3" (4 bytes)
  Version: uint16 (2 bytes)
  Flags: uint16 (2 bytes)
  Reference ID: uint32 (4 bytes)
  Start position: uint32 (4 bytes)
  End position: uint32 (4 bytes)
  Read count: uint32 (4 bytes)
  Compressed size: uint64 (8 bytes)
  Uncompressed size: uint64 (8 bytes)
  Compression: uint8 (1 byte) - 0=none, 1=zstd, 2=lz4
  Checksum algorithm: uint8 (1 byte)
  Checksum: (32 bytes)
  Reserved: (64 bytes)

[INDEX SECTION: variable]
  Index offset: uint64 (where in file the index starts)
  Index size: uint32
  Index data: [array of read offsets within DATA section]
    - Allows random access to individual reads within chunk

[DATA SECTION: variable, compressed]
  Compressed read records

  Each read record (uncompressed format):
    Record length: uint32
    Read name length: uint8
    Read name: variable
    Flag: uint16
    Reference ID: uint32
    Position: uint32
    Mapping quality: uint8
    CIGAR length: uint16
    CIGAR: variable (packed format)
    Mate reference ID: uint32
    Mate position: uint32
    Template length: int32
    Sequence length: uint32
    Sequence: variable (4-bit encoding: A=0,C=1,G=2,T=3,N=15)
    Quality: variable (8-bit per base)
    Tags: variable (SAM tag format)
```

### _index/spatial.idx

Quick lookup: position → chunk

```
SPATIAL INDEX FORMAT
====================

[HEADER]
  Magic: "IDX1"
  Reference count: uint32

[PER-REFERENCE SECTION]
  Reference ID: uint32
  Bin count: uint32

  [BINS]
    For each bin (1MB regions):
      Start position: uint32
      End position: uint32
      Chunk path: string (offset into string table)
      Read count: uint32
```

## Usage Examples

### Query Region

```python
import bams3

# Open BAMS3 dataset
dataset = bams3.open('s3://bucket/sample.bams3')

# Query specific region
reads = dataset.query('chr1', 1000000, 2000000)

# Efficient: Only downloads relevant chunks
# - Reads _metadata.json (small)
# - Reads spatial.idx (small)
# - Downloads data/chr1/01000000-02000000.chunk only

for read in reads:
    print(read.name, read.position, read.sequence)
```

### Parallel Processing

```python
import bams3
from concurrent.futures import ProcessPoolExecutor

dataset = bams3.open('s3://bucket/sample.bams3')

# Get all chunks for parallel processing
chunks = dataset.get_chunks('chr1')

# Process in parallel
def process_chunk(chunk_path):
    reads = bams3.read_chunk(chunk_path)
    return calculate_coverage(reads)

with ProcessPoolExecutor(max_workers=32) as executor:
    results = executor.map(process_chunk, chunks)

# 32 workers downloading different chunks in parallel!
```

### Stream All Reads

```python
# Sequential scan of entire dataset
dataset = bams3.open('s3://bucket/sample.bams3')

for read in dataset.iter_reads():
    process(read)

# Internally: streams chunks one at a time
```

### Convert to/from BAM

```python
import bams3

# BAM → BAMS3
bams3.convert(
    input='sample.bam',
    output='s3://bucket/sample.bams3',
    chunk_size=1000000,  # 1Mbp per chunk
    compression='zstd',
    compression_level=3,
    parallel=True
)

# BAMS3 → BAM
bams3.to_bam(
    input='s3://bucket/sample.bams3',
    output='sample.bam'
)
```

## Implementation Strategy

### Phase 1: Core Library (Python)

```python
# bams3/__init__.py

class BAMS3Dataset:
    """Main interface for BAMS3 datasets."""

    def __init__(self, uri):
        """Open BAMS3 dataset from S3 URI."""
        self.uri = uri
        self.metadata = self._read_metadata()
        self.header = self._read_header()
        self.index = self._read_index()

    def query(self, reference, start, end):
        """Query reads in genomic region."""
        # 1. Use index to find relevant chunks
        chunks = self.index.get_chunks(reference, start, end)

        # 2. Download and decompress chunks
        for chunk_info in chunks:
            chunk = self._read_chunk(chunk_info['path'])

            # 3. Filter reads within exact region
            for read in chunk.reads:
                if start <= read.position < end:
                    yield read

    def iter_reads(self):
        """Iterate all reads sequentially."""
        for chunk_info in self.metadata['chunks']:
            chunk = self._read_chunk(chunk_info['path'])
            yield from chunk.reads

    def get_coverage(self, reference, start, end):
        """Get coverage in region (uses pre-computed if available)."""
        pass

    def statistics(self):
        """Get dataset statistics."""
        return self.metadata['statistics']
```

### Phase 2: Conversion Tools

```python
# bams3/convert.py

def bam_to_bams3(
    input_bam: str,
    output_uri: str,
    chunk_size: int = 1000000,
    compression: str = 'zstd',
    parallel: bool = True
):
    """
    Convert BAM to BAMS3 format.

    Process:
    1. Read BAM header → create _header.json
    2. Partition reads by position into chunks
    3. Compress each chunk
    4. Upload to S3
    5. Create metadata and index
    """

    import pysam
    import boto3

    # Open input BAM
    bam = pysam.AlignmentFile(input_bam, 'rb')

    # Initialize output structure
    s3 = boto3.client('s3')
    bucket, prefix = parse_s3_uri(output_uri)

    # Create header
    header = convert_sam_header_to_json(bam.header)
    upload_json(s3, bucket, f"{prefix}/_header.json", header)

    # Process reads into chunks
    chunks = []
    current_chunk = ChunkWriter(
        reference='chr1',
        start=0,
        end=chunk_size,
        compression=compression
    )

    for read in bam:
        if read.reference_start >= current_chunk.end:
            # Finish current chunk
            chunk_data = current_chunk.finalize()
            chunk_path = f"data/{read.reference_name}/{current_chunk.start:09d}-{current_chunk.end:09d}.chunk"
            upload_chunk(s3, bucket, f"{prefix}/{chunk_path}", chunk_data)

            chunks.append({
                'path': chunk_path,
                'reference': read.reference_name,
                'start': current_chunk.start,
                'end': current_chunk.end,
                'reads': len(current_chunk),
                # ... other metadata
            })

            # Start next chunk
            current_chunk = ChunkWriter(...)

        current_chunk.add_read(read)

    # Create metadata
    metadata = create_metadata(chunks, header, statistics)
    upload_json(s3, bucket, f"{prefix}/_metadata.json", metadata)

    # Create index
    index = create_spatial_index(chunks)
    upload_index(s3, bucket, f"{prefix}/_index/spatial.idx", index)
```

### Phase 3: CLI Tools

```bash
# Install
pip install bams3

# Convert BAM → BAMS3
bams3 convert \
    --input sample.bam \
    --output s3://bucket/sample.bams3 \
    --chunk-size 1000000 \
    --compression zstd \
    --parallel

# Query BAMS3
bams3 view \
    s3://bucket/sample.bams3 \
    chr1:1000000-2000000

# Statistics
bams3 stats s3://bucket/sample.bams3

# Convert back to BAM
bams3 to-bam \
    --input s3://bucket/sample.bams3 \
    --output sample.bam
```

### Phase 4: Integration with Existing Tools

```python
# Adapter: Make BAMS3 appear as BAM to existing tools

class BAMS3ToBamAdapter:
    """
    Presents BAMS3 as a BAM file to existing tools.

    Uses FUSE or similar to create virtual BAM file.
    Fetches BAMS3 chunks on-demand.
    """

    def __init__(self, bams3_uri, cache_dir):
        self.dataset = bams3.open(bams3_uri)
        self.cache = DiskCache(cache_dir)

    def read(self, offset, length):
        """
        Read bytes at offset.

        Transparently converts BAMS3 chunks to BAM format.
        """
        # Map offset to BAMS3 chunk
        chunk = self._offset_to_chunk(offset)

        # Convert to BAM format on-the-fly
        bam_data = self._chunk_to_bam(chunk)

        return bam_data[offset:offset+length]
```

## Performance Expectations

### Query Single Region (1Mbp on chr1)

| Method | Time | Data Downloaded | Notes |
|--------|------|-----------------|-------|
| BAM on S3 (FUSE) | 12s | 12GB (entire file) | First access |
| BAM on S3 (copy) | 125s | 12GB | Download + query |
| BAMS3 | 0.8s | 1-2MB (one chunk) | **15x faster** |

### Full Sequential Scan

| Method | Time | Data Downloaded | Notes |
|--------|------|-----------------|-------|
| BAM local | 180s | N/A | Baseline |
| BAM on S3 (stream) | 250s | 12GB | |
| BAMS3 (parallel) | 95s | 12GB | **2.6x faster** with 32 workers |

### 100 Random Region Queries

| Method | Time | Data Downloaded | Notes |
|--------|------|-----------------|-------|
| BAM on S3 (FUSE, cold) | 1200s | 12GB × 100 | Terrible! |
| BAM on S3 (FUSE, cached) | 15s | 12GB once | After first query |
| BAM local | 10s | N/A | |
| BAMS3 | 12s | ~200MB | Only relevant chunks |

## Migration Path

### For New Data

**Just use BAMS3 from the start:**

```bash
# After alignment
bwa mem ref.fa reads.fq | \
    samtools view -b | \
    bams3 convert --output s3://bucket/sample.bams3
```

### For Existing Data

**Gradual migration:**

```bash
# Convert high-value datasets
bams3 convert \
    --input s3://old-bucket/*.bam \
    --output s3://new-bucket/{sample}.bams3 \
    --parallel

# Keep BAMs initially for compatibility
# After validation, delete BAMs
```

### For Tool Compatibility

**Provide both formats:**

```
s3://bucket/
  bam/
    sample.bam       # For legacy tools
  bams3/
    sample.bams3/    # For modern tools
```

## Next Steps

1. **Implement core library** (Python first)
2. **Create conversion tools** (BAM ↔ BAMS3)
3. **Build CLI tools** (view, stats, query)
4. **Benchmark** with real data
5. **Refine format** based on results
6. **Document** and release
7. **Build community** around format

## Format Extensions

Future versions could add:

- **BAMS3-QC**: Embedded quality metrics
- **BAMS3-VCF**: Variants alongside alignments
- **BAMS3-Multi**: Multiple samples in one dataset
- **BAMS3-Sparse**: For low-coverage data (long reads)

## Open Questions

1. **Chunk size**: 1Mbp? 10Mbp? Configurable?
2. **Compression**: zstd, lz4, or both?
3. **Index format**: Current simple index vs. more sophisticated?
4. **Compatibility**: How much BAM compatibility to maintain?
5. **Standardization**: Try to standardize format or keep experimental?

## Get Involved

See `format-tools/bams3/` for:
- Reference implementation
- Conversion tools
- Examples and benchmarks
- Discussion and development

**The future is object-native!** 🚀
