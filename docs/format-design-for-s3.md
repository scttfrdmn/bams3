# Designing Data Formats for Object Storage

## The Big Idea

Instead of adapting S3 to work with POSIX formats, **design formats that are native to object storage**.

**Traditional thinking:** "How do we make BAM work on S3?"

**Better thinking:** "What format would we design if we started with S3 as the primary storage?"

## Key Principles for S3-Native Formats

### 1. Object-Per-Chunk Architecture

**POSIX format:**
```
huge_file.bam   (single 500GB file)
```

**S3-native format:**
```
dataset/
  _metadata           (manifest, schema)
  chunk_0000          (independent unit)
  chunk_0001
  chunk_0002
  ...
  chunk_9999
```

**Why:** Each chunk is:
- Independently readable (parallel access)
- Right-sized for network requests (1-100MB)
- Cacheable unit
- Processing unit

### 2. Self-Describing Chunks

Each chunk contains:
- Schema/format version
- Compression codec
- Byte offsets/indexes
- Checksums

**Why:** Can process chunk without reading others.

### 3. Metadata Co-location

**Bad:**
```
Must read file to know what's in it
```

**Good:**
```
_metadata file describes all chunks:
- What regions they contain
- Statistics (min/max, counts)
- File locations
- Provenance
```

**Why:** Query metadata before fetching data.

### 4. Columnar When Appropriate

**Row-oriented (traditional):**
```
Record 1: [field1, field2, field3, ...]
Record 2: [field1, field2, field3, ...]
```

**Column-oriented (S3-optimized):**
```
Column: field1 [val1, val2, val3, ...]
Column: field2 [val1, val2, val3, ...]
Column: field3 [val1, val2, val3, ...]
```

**Why:** Query only columns you need.

### 5. Indexes as First-Class Objects

**Traditional:** Index is separate file (.bai, .tbi)

**S3-native:**
```
dataset/
  _index/
    spatial_index     (chromosome, position)
    sample_index      (by sample ID)
    quality_index     (by QC metrics)
  data/
    chunk_0000
    ...
```

**Why:**
- Multiple indexes for different query patterns
- Indexes are objects too
- Can prefetch index with data

## Real Examples

## Example 1: Genomic Reads - From BAM to "S3-Optimized Alignments"

### Current: BAM Format
```
sample.bam (single file)
├─ Header
├─ Read1 (chr1:1000)
├─ Read2 (chr1:1050)
├─ Read3 (chr2:5000)
├─ ... (millions more, interleaved)
└─ [Compressed with BGZF]

Needs: sample.bam.bai (index file)
```

**Problems for S3:**
- Single monolithic file
- Must download all or use FUSE
- Index is separate object
- Difficult to parallelize

### S3-Native: Partitioned Alignment Format (PAF)

```
sample.paf/
  _metadata.json              # Global metadata
  _index/
    chr1.idx                  # Spatial index
    chr2.idx
    ...
  chr1/
    chunk_00000-01000000.bin  # 1MB chunks
    chunk_01000001-02000000.bin
    ...
  chr2/
    chunk_00000-01000000.bin
    ...
  chrX/
    ...
```

**Format spec:**

**_metadata.json:**
```json
{
  "format": "paf",
  "version": "1.0",
  "sample": "NA12878",
  "reference": "hg38",
  "chromosomes": ["chr1", "chr2", ..., "chrY"],
  "chunk_size": 1000000,
  "compression": "zstd",
  "created": "2024-01-15T10:30:00Z",
  "read_count": 450000000,
  "chunks": [
    {
      "path": "chr1/chunk_00000-01000000.bin",
      "reads": 12534,
      "size_bytes": 1048576,
      "checksum": "sha256:abc123..."
    },
    ...
  ]
}
```

**Chunk format (binary):**
```
[HEADER: 64 bytes]
  - Magic number: "PAF1"
  - Version: uint16
  - Compression: uint8
  - Chromosome: uint8
  - Start position: uint32
  - End position: uint32
  - Read count: uint32
  - Checksum: 32 bytes

[INDEX: variable]
  - Per-read offsets for random access within chunk

[DATA: variable]
  - Compressed read records
  - Each record: [name, sequence, quality, position, CIGAR, flags]
```

**Benefits:**
1. **Parallel processing:** Process each chromosome independently
2. **Range queries:** Fetch only chr1:1M-2M chunks
3. **No index needed:** Metadata tells you which chunks to fetch
4. **Streaming friendly:** Download chunks as needed
5. **Cacheable:** Cache hot chromosomes/regions

**Tools needed:**
- `bam-to-paf` - Convert BAM → PAF
- `paf-query` - Query PAF by region
- `paf-to-bam` - Convert back if needed

### Example 2: Variants - From VCF to Columnar Store

### Current: VCF/BCF
```
variants.vcf.gz
  ##header
  #CHROM  POS     ID      REF ALT QUAL FILTER INFO FORMAT Sample1 Sample2
  chr1    100     rs123   A   G   99   PASS   DP=50  GT:DP  0/1:25  0/0:25
  chr1    200     rs124   C   T   99   PASS   DP=48  GT:DP  1/1:24  0/1:24
  ...
```

**Problems for S3:**
- Row-oriented (must read entire record)
- Compressed as whole (can't seek easily)
- Queries like "all variants with DP>30" require full scan

### S3-Native: Variant Columnar Store (VCS)

```
variants.vcs/
  _metadata.json
  _schema.json
  columns/
    CHROM/
      chunk_000.parquet
      chunk_001.parquet
    POS/
      chunk_000.parquet
      chunk_001.parquet
    REF/
      chunk_000.parquet
    ALT/
      chunk_000.parquet
    QUAL/
      chunk_000.parquet
    INFO_DP/
      chunk_000.parquet
    samples/
      Sample1_GT/
        chunk_000.parquet
      Sample1_DP/
        chunk_000.parquet
```

**Using Apache Parquet with optimizations:**

```python
import pyarrow as pa
import pyarrow.parquet as pq

# Write VCF as columnar
schema = pa.schema([
    ('CHROM', pa.string()),
    ('POS', pa.int64()),
    ('REF', pa.string()),
    ('ALT', pa.string()),
    ('QUAL', pa.float32()),
    ('DP', pa.int32()),
    # ... more fields
])

# Create row groups by chromosome
pq.write_table(
    table,
    'variants.vcs/data.parquet',
    row_group_size=100000,  # 100k variants per chunk
    compression='zstd',
    use_dictionary=True,    # Compress repeated values (chromosomes)
    write_statistics=True   # Enable predicate pushdown
)
```

**Query examples:**

```python
import pyarrow.parquet as pq
import pyarrow.dataset as ds

# Open dataset on S3
dataset = ds.dataset('s3://bucket/variants.vcs/', format='parquet')

# Query 1: Get only CHROM and POS for chr1 (selective columns + predicate)
table = dataset.to_table(
    columns=['CHROM', 'POS'],
    filter=ds.field('CHROM') == 'chr1'
)
# Only reads CHROM and POS columns, only chr1 row groups!

# Query 2: High quality variants
table = dataset.to_table(
    filter=ds.field('QUAL') > 30
)
# Pushes filter to storage - skips entire chunks where max(QUAL) < 30

# Query 3: Complex query
table = dataset.to_table(
    columns=['CHROM', 'POS', 'ALT'],
    filter=(
        (ds.field('CHROM') == 'chr1') &
        (ds.field('POS') > 1000000) &
        (ds.field('POS') < 2000000) &
        (ds.field('QUAL') > 20)
    )
)
# Massive speedup - only fetches relevant data!
```

**Performance:**
- **Traditional VCF:** Read entire 10GB file
- **Columnar VCS:** Read only 100MB for selective query (100x less data!)

### Example 3: FASTQ - Structured Read Store

### Current: FASTQ
```
sample_R1.fastq.gz  (single compressed file)
@READ1
ACGTACGTACGT
+
IIIIIIIIIIII
@READ2
GCTAGCTAGCTA
+
JJJJJJJJJJJJ
...
```

**Problems for S3:**
- Sequential access only
- GZIP not seekable
- No metadata (must scan to count reads)
- No indexing

### S3-Native: Structured Read Format (SRF)

```
sample.srf/
  _metadata.json
  chunks/
    0000000000-0000999999.srf  # First 1M reads
    0001000000-0001999999.srf  # Next 1M reads
    ...
  index/
    by_quality.idx              # Index by quality scores
    by_length.idx               # Index by read length
```

**Chunk format (using Apache Arrow IPC):**

```python
import pyarrow as pa

# Schema for reads
schema = pa.schema([
    ('read_id', pa.string()),
    ('sequence', pa.binary()),          # Raw bytes (ACGT)
    ('quality', pa.list_(pa.int8())),   # Quality scores
    ('length', pa.int32()),
    ('avg_quality', pa.float32()),      # Pre-computed
])

# Write chunk
chunk = pa.Table.from_arrays([...], schema=schema)
with pa.OSFile('chunk_0000.srf', 'wb') as f:
    with pa.RecordBatchFileWriter(f, schema) as writer:
        writer.write_table(chunk)
```

**Benefits:**

1. **Random access:** Jump to read 1,000,000 without scanning
2. **Selective queries:**
   ```python
   # Get only high-quality reads
   reads = read_chunk('chunk_0000.srf')
   high_qual = reads.filter(reads['avg_quality'] > 30)
   ```

3. **Parallel processing:**
   ```bash
   # Process chunks in parallel
   for chunk in chunks/*.srf; do
       process_chunk $chunk &
   done
   ```

4. **Metadata available:**
   ```json
   {
     "total_reads": 450000000,
     "total_bases": 67500000000,
     "avg_quality": 35.2,
     "read_length": 150,
     "chunks": 450
   }
   ```

## Example 4: Bespoke Format - Alignment Matrix

### Use Case: Multi-sample genomic region analysis

**Problem:** Want to analyze chr1:1M-2M across 1000 samples.

**Traditional approach:**
```
Download 1000 BAM files (100GB each)
Extract region from each
= 100TB download for ~1GB of actual data needed!
```

**S3-Native: Regional Alignment Matrix (RAM)**

```
cohort-1000samples.ram/
  _metadata.json
  regions/
    chr1_00000000-01000000/
      _region_metadata.json
      alignments.arrow          # All 1000 samples, this region
      coverage.arrow            # Pre-computed coverage
      variants.arrow            # Pre-called variants
    chr1_01000000-02000000/
      ...
```

**Format design:**

```python
# Each region is a columnar table with all samples
schema = pa.schema([
    ('sample_id', pa.string()),
    ('read_name', pa.string()),
    ('position', pa.int32()),
    ('sequence', pa.binary()),
    ('quality', pa.list_(pa.int8())),
    ('mapping_quality', pa.int8()),
    ('flags', pa.int16()),
])

# All reads from all samples in this region
# Organized by position for cache efficiency
```

**Usage:**

```python
# Analyze one region across all samples
region_data = read_ram_region('chr1_01000000-02000000')

# Extract just one sample
sample_reads = region_data.filter(
    region_data['sample_id'] == 'NA12878'
)

# Compute coverage across all samples
coverage = compute_coverage_by_sample(region_data)

# Download size: 100MB instead of 100TB!
```

**When to use:**
- Multi-sample cohort analysis
- Region-focused studies
- Population genomics

**Creation:**
```bash
# Pre-process 1000 BAMs into RAM format
bam-to-ram \
    --input sample_*.bam \
    --output cohort.ram \
    --region-size 1000000 \
    --parallel 32
```

## Design Patterns

### Pattern 1: Partition by Query Dimension

**Principle:** Organize by how you'll query.

**Example:**
- Query by chromosome? → Partition by chromosome
- Query by sample? → Partition by sample
- Query by time? → Partition by date

### Pattern 2: Hierarchical Chunking

```
dataset/
  level0/          # Coarse chunks (100MB)
    chunk_000/
      level1/      # Fine chunks (1MB)
        chunk_000
        chunk_001
```

**Use case:**
- Coarse for full scans
- Fine for targeted queries

### Pattern 3: Multiple Views

Store data in multiple formats optimized for different access patterns:

```
dataset/
  raw/             # Original data
  by_chromosome/   # For regional queries
  by_sample/       # For sample-level analysis
  columnar/        # For analytical queries
```

**Trade-off:** Storage space vs query speed

### Pattern 4: Immutable Chunks + Manifest

**Chunks never change** (immutable)
**Manifest tracks current dataset state**

```
data/chunks/
  abc123.chunk    # Immutable
  def456.chunk    # Immutable

manifests/
  v1.json         # Points to chunks [abc123]
  v2.json         # Points to chunks [abc123, def456]
  v3.json         # Points to chunks [def456, ghi789]  (removed abc123)
```

**Benefits:**
- Versioning for free
- Concurrent writes (new chunks)
- Rollback (use old manifest)
- Deduplication (same chunk in multiple versions)

## Conversion Strategy

### Phase 1: Assessment

1. **Analyze access patterns:**
   - Sequential vs random?
   - Full file vs regions?
   - Which fields queried?
   - Read vs write frequency?

2. **Measure current performance:**
   - Time to access
   - Data volume transferred
   - Storage costs

### Phase 2: Format Design

1. **Choose chunking strategy:**
   - Spatial (genomic regions)?
   - Temporal (time ranges)?
   - Logical (samples, experiments)?

2. **Select encoding:**
   - Columnar (Parquet, Arrow)?
   - Row-oriented (custom binary)?
   - Hybrid?

3. **Design metadata:**
   - What should be queryable without reading chunks?
   - What indexes are needed?

### Phase 3: Implementation

1. **Build converters:**
   - Original format → S3-native format
   - S3-native format → Original format (for compatibility)

2. **Create access libraries:**
   - Python: Read/write format
   - Command-line tools
   - Integration with existing tools

### Phase 4: Migration

1. **Convert a dataset:**
2. **Benchmark:**
   - Compare old vs new format performance
   - Measure storage efficiency
   - Test typical queries

3. **Iterate:**
   - Adjust chunk size
   - Refine indexes
   - Optimize encoding

4. **Scale:**
   - Convert remaining datasets
   - Update pipelines
   - Document best practices

## Tools and Libraries

### Existing Formats
- **Parquet** - Columnar, analytics
- **Arrow IPC** - Columnar, interop
- **Zarr** - Chunked arrays
- **HDF5 (cloud-optimized)** - Scientific data
- **TileDB** - Multi-dimensional arrays

### Build Your Own
- **Python:** `pyarrow`, `struct`, `pickle`
- **Rust:** `serde`, `bincode`
- **C++:** `flatbuffers`, `protobuf`

### Design Checklist
- [ ] Chunks are independently readable
- [ ] Chunk size appropriate (1-100MB)
- [ ] Metadata separate from data
- [ ] Schema versioning
- [ ] Compression appropriate for data
- [ ] Indexes for common queries
- [ ] Checksums for integrity
- [ ] Documentation and examples

## Success Metrics

**Before (POSIX format on S3):**
- Query time: 120s (download 10GB)
- Data transferred: 10GB
- Parallelization: Difficult

**After (S3-native format):**
- Query time: 5s (download 50MB)
- Data transferred: 50MB (200x less!)
- Parallelization: Trivial (process chunks independently)

## Next Steps

See `format-tools/` directory for:
- Conversion utilities
- Format specifications
- Example implementations
- Benchmarks

The future of genomics data is **object-native formats** designed for cloud storage from the ground up!
