# Format Optimization for S3

## Core Principle

**Local filesystem assumption:** Files are cheap to open, sequential reads are fast, seeking is nearly free.

**S3 reality:** Opening objects has latency, sequential reads need large chunks, seeking requires new requests.

**Solution:** Design formats that minimize requests and maximize throughput.

## Key Optimizations

### 1. Self-Describing Files

**Problem:** Need to download file just to check what's in it.

**Solution:** Embed metadata in predictable locations.

#### Example: BAM Header
```
BAM structure:
[Header with metadata] [Index] [Data blocks...]
```

**Good:** Header is at the start, can read with a small range request.

**Better:** Add extended metadata in a standard location:
- Sample ID
- Sequencing platform
- Quality metrics
- Processing provenance

```bash
# Can query metadata without downloading file
aws s3api head-object --bucket my-bucket --key sample1.bam \
    --query 'Metadata'
```

### 2. Chunking and Partitioning

**Problem:** Must download entire file to access small portion.

**Solution:** Split into smaller, independently accessible chunks.

#### Spatial Partitioning (Genomics)
```
# Instead of:
huge_dataset.bam  (500 GB)

# Use:
dataset/
  chr1.bam
  chr2.bam
  ...
  chr22.bam
  chrX.bam
  chrY.bam
  chrM.bam
```

**Benefits:**
- Process chromosomes in parallel
- Only download chromosomes you need
- Better S3 request distribution

#### Temporal/Sample Partitioning
```
cohort/
  batch_2024_01/
  batch_2024_02/
  ...
```

**Partitioning tool pseudocode:**
```python
def partition_bam_by_chromosome(input_bam, output_prefix):
    """Split BAM into per-chromosome files."""
    chromosomes = ['chr1', 'chr2', ..., 'chrY']
    writers = {chr: open(f"{output_prefix}/{chr}.bam", 'w')
               for chr in chromosomes}

    for read in read_bam(input_bam):
        chr = read.reference_name
        writers[chr].write(read)

    for writer in writers.values():
        writer.close()
        # Upload to S3
```

### 3. Columnar Storage

**Problem:** Need one column from file, must download all columns.

**Solution:** Use columnar formats that allow column selection.

#### VCF → Parquet Conversion

**VCF (row-oriented):**
```
CHROM  POS     ID      REF ALT QUAL FILTER INFO
chr1   100000  rs123   A   G   99   PASS   DP=50;AF=0.5
chr1   100100  rs124   C   T   99   PASS   DP=48;AF=0.3
```

To get just `CHROM` and `POS`, must read entire file.

**Parquet (columnar):**
```
Column: CHROM [chr1, chr1, chr1, ...]
Column: POS   [100000, 100100, ...]
Column: ALT   [G, T, ...]
...
```

Can read just `CHROM` and `POS` columns!

**Conversion tool:**
```python
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

def vcf_to_parquet(vcf_file, output_parquet):
    # Parse VCF
    df = parse_vcf_to_dataframe(vcf_file)

    # Write as Parquet with optimal settings
    pq.write_table(
        pa.Table.from_pandas(df),
        output_parquet,
        compression='snappy',
        use_dictionary=True,  # Good for categorical data
        row_group_size=100000  # Tune for access patterns
    )
```

**S3 access:**
```python
import pyarrow.parquet as pq
import s3fs

fs = s3fs.S3FileSystem()

# Read only specific columns (efficient!)
table = pq.read_table(
    's3://bucket/variants.parquet',
    columns=['CHROM', 'POS', 'ALT'],
    filesystem=fs
)

# Apply predicates (even more efficient!)
table = pq.read_table(
    's3://bucket/variants.parquet',
    columns=['CHROM', 'POS', 'ALT'],
    filters=[('CHROM', '=', 'chr1')],
    filesystem=fs
)
```

### 4. Index Co-location

**Problem:** Need index file to efficiently access data file. Requires two requests.

**Solution:** Store index predictably, prefetch with data.

#### BAM + BAI Pattern

**Current practice:**
```
s3://bucket/sample.bam
s3://bucket/sample.bam.bai
```

**Optimization 1: Predictable naming**
Tools should automatically look for `.bai` when opening `.bam`.

**Optimization 2: Embed index** (format change)
```
Modified BAM:
[Header] [Embedded Index] [Data]
```

Can read index with single range request to header.

**Optimization 3: Index in S3 metadata**
For very small indexes, store in object metadata:
```bash
aws s3api put-object \
    --bucket my-bucket \
    --key sample.bam \
    --metadata "index-locations=chr1:0-1000000,chr2:1000001-2000000"
```

### 5. Compression Trade-offs

**Local files:** Aggressive compression saves disk space.
**S3:** Different considerations.

#### Compression Factors

| Format | Compression | Seekable | S3 Efficiency |
|--------|-------------|----------|---------------|
| Uncompressed | 1x | Yes | ⭐⭐⭐⭐⭐ |
| GZIP | 3-5x | No | ⭐ |
| BGZF (BAM) | 3-5x | Yes | ⭐⭐⭐⭐ |
| CRAM | 5-10x | Yes | ⭐⭐⭐⭐ |
| LZ4 | 2-3x | Yes | ⭐⭐⭐⭐⭐ |

**GZIP problem:**
```
compressed.fastq.gz (50GB)
Want reads 1000-2000.
Must decompress from start → download entire 50GB!
```

**BGZF solution (block GZIP):**
```
file.bam uses BGZF
├─ Block 1 (65KB compressed)
├─ Block 2 (65KB compressed)
...
└─ Block N

Index says: reads 1000-2000 are in blocks 50-55
→ Request only those blocks!
```

**Recommendation for new formats:**
- Use block-based compression (BGZF, LZ4)
- Make blocks ~64KB - 1MB
- Allow random access to blocks

### 6. Cloud-Optimized Variants

#### Cloud-Optimized BAM (COBAM)

Hypothetical optimized format:
```
Structure:
[Header + Metadata]
[Block Index - embedded]
[Data blocks - 1MB each]
[Footer with statistics]

Features:
- Larger block size (1MB vs 64KB)
- Embedded index
- Richer metadata
- Optional: column store within blocks (reads vs qualities separate)
```

#### Cloud-Optimized FASTQ

```
FASTQ is inherently bad for S3:
- Must read sequentially
- GZIP makes it worse
- No index

Better alternatives:
1. Convert to BAM/CRAM (indexed, compressed)
2. Chunk into smaller files:
   sample_R1.part001.fastq.gz (1GB each)
   sample_R1.part002.fastq.gz
   ...
3. Use structured format:
   - Parquet with (read_id, sequence, quality)
   - SQLite database
```

### 7. Zarr for Array Data

**Use case:** Multidimensional array data (imaging, tensor data).

**Zarr advantages:**
- Chunked N-dimensional arrays
- Each chunk is separate S3 object
- Access only chunks you need
- Built for cloud

**Structure:**
```
dataset.zarr/
  .zarray              # Metadata
  0.0.0                # Chunk [0,0,0]
  0.0.1                # Chunk [0,0,1]
  1.0.0                # Chunk [1,0,0]
  ...
```

**Python usage:**
```python
import zarr
import s3fs

fs = s3fs.S3FileSystem()
store = s3fs.S3Map('s3://bucket/data.zarr', s3=fs)
z = zarr.open(store, mode='r')

# Access subset - only relevant chunks downloaded
subset = z[0:100, 0:100, 0:100]
```

## Format Conversion Strategy

### Assessment

For each dataset, ask:
1. **Access pattern:** Sequential or random?
2. **Query pattern:** Whole file or subsets?
3. **Update frequency:** Read-only or append?
4. **File size:** Small (<1GB) or large?

### Decision Matrix

| Use Case | Recommended Format | Rationale |
|----------|-------------------|-----------|
| Large BAM, random access by region | Partitioned BAM (by chr) | Parallel access, smaller files |
| VCF analysis (column queries) | Parquet | Columnar, predicate pushdown |
| FASTQ (must keep) | Chunked FASTQ or convert to BAM | Enable parallel processing |
| Reference genomes | BGZF FASTA with .fai | Existing tools support it |
| Array data (imaging) | Zarr | Purpose-built for chunked access |

### Conversion Pipeline Example

```bash
#!/bin/bash
# Optimize a BAM file for S3

INPUT_BAM="s3://raw-data/sample.bam"
OUTPUT_PREFIX="s3://optimized-data/sample"

# 1. Validate input
samtools quickcheck $INPUT_BAM

# 2. Split by chromosome
for chr in chr1 chr2 chr3 ... chrY; do
    samtools view -b $INPUT_BAM $chr > ${chr}.bam &
done
wait

# 3. Create indexes
for bam in chr*.bam; do
    samtools index $bam &
done
wait

# 4. Upload to S3 with metadata
for bam in chr*.bam; do
    chr=$(basename $bam .bam)
    aws s3 cp $bam ${OUTPUT_PREFIX}/${chr}.bam \
        --metadata "chromosome=$chr,source=$INPUT_BAM,optimized=true"
    aws s3 cp ${bam}.bai ${OUTPUT_PREFIX}/${chr}.bam.bai
done

# 5. Create manifest
cat > manifest.json <<EOF
{
  "sample": "sample1",
  "chromosomes": ["chr1", ..., "chrY"],
  "format": "BAM",
  "optimized": true,
  "created": "$(date -Iseconds)"
}
EOF
aws s3 cp manifest.json ${OUTPUT_PREFIX}/manifest.json
```

## Validation

### Is a File Cloud-Optimized?

Create a validator tool:

```python
def validate_cloud_optimization(file_path):
    checks = []

    # Check 1: Is it chunked/partitioned?
    if is_single_large_file(file_path):
        checks.append("❌ Large monolithic file - consider partitioning")
    else:
        checks.append("✅ Partitioned or reasonably sized")

    # Check 2: Does index exist?
    if requires_index(file_path) and not index_exists(file_path):
        checks.append("❌ Missing index file")
    else:
        checks.append("✅ Index present or not needed")

    # Check 3: Compression type
    compression = detect_compression(file_path)
    if compression == "GZIP" and file_size(file_path) > 1GB:
        checks.append("❌ GZIP on large file - not seekable")
    elif compression in ["BGZF", "LZ4"]:
        checks.append("✅ Block compression - seekable")

    # Check 4: Metadata present?
    metadata = get_s3_metadata(file_path)
    if not metadata:
        checks.append("⚠️  No S3 metadata - consider adding")

    return checks
```

## Measuring Success

### Metrics

1. **Request count:** Fewer requests = better
2. **Data transfer:** Less data = better
3. **Latency:** Time to first byte
4. **Throughput:** Overall processing speed

### Benchmark

```python
import time
import boto3

def benchmark_access(file_path, access_pattern):
    """Compare optimized vs non-optimized."""

    s3 = boto3.client('s3')

    start = time.time()

    if access_pattern == "sequential":
        # Read entire file
        obj = s3.get_object(Bucket=bucket, Key=key)
        data = obj['Body'].read()

    elif access_pattern == "random":
        # Read 10 random 1MB chunks
        for _ in range(10):
            offset = random.randint(0, file_size - 1MB)
            range_header = f"bytes={offset}-{offset+1MB}"
            obj = s3.get_object(Bucket=bucket, Key=key, Range=range_header)
            chunk = obj['Body'].read()

    elapsed = time.time() - start
    print(f"{file_path}: {elapsed:.2f}s")
```

## Next Steps

See `format-tools/` directory for:
- Conversion utilities
- Validation scripts
- Benchmark tools
- Example pipelines
