# Streaming Examples

## Overview

For tools that accept stdin/stdout, you can stream data directly from S3 without downloading entire files first.

**Advantages:**
- No local storage needed
- Start processing immediately
- Works with unmodified tools

**Limitations:**
- Sequential access only (no seeking)
- No random access
- Tool must support stdin/stdout

## Simple Examples

### Using AWS CLI

```bash
# Stream and decompress FASTQ
aws s3 cp s3://bucket/sample.fastq.gz - | gunzip | head -n 1000

# Stream BAM, convert to SAM
aws s3 cp s3://bucket/sample.bam - | samtools view - | less

# Chain multiple operations
aws s3 cp s3://bucket/reads.fastq.gz - | \
    gunzip | \
    fastq-stats - | \
    tee stats.txt
```

### Using Custom Streaming Tool

The `s3-stream.py` utility provides better control:

```bash
# Make executable
chmod +x s3-stream.py

# Stream FASTQ
./s3-stream.py s3://bucket/sample.fastq.gz | gunzip | head

# With custom chunk size (for faster streaming)
./s3-stream.py --chunk-size 8388608 s3://bucket/large.bam | samtools view -
```

## Real Genomics Workflows

### FASTQ Quality Check

```bash
#!/bin/bash
# Stream FASTQ from S3 and run FastQC

S3_FASTQ="s3://my-bucket/sample_R1.fastq.gz"
OUTPUT_DIR="./qc-results"

mkdir -p $OUTPUT_DIR

# Stream to FastQC (accepts stdin)
aws s3 cp $S3_FASTQ - | fastqc stdin -o $OUTPUT_DIR

# Or if tool needs a filename, use process substitution
fastqc <(aws s3 cp $S3_FASTQ -) -o $OUTPUT_DIR
```

### Read Filtering

```bash
#!/bin/bash
# Stream, filter low-quality reads, upload result

INPUT="s3://bucket/raw.fastq.gz"
OUTPUT="s3://bucket/filtered.fastq.gz"

aws s3 cp $INPUT - | \
    gunzip | \
    fastp --stdin --stdout \
        --qualified_quality_phred 20 \
        --length_required 50 | \
    gzip | \
    aws s3 cp - $OUTPUT
```

### BAM Region Extraction

```bash
#!/bin/bash
# Extract specific chromosome without downloading entire BAM

BAM="s3://bucket/large.bam"
REGION="chr1:1000000-2000000"

# Note: This only works if tool supports seeking
# For full BAM, need FUSE mount or tool modification

# Stream entire BAM (inefficient but works)
aws s3 cp $BAM - | samtools view -b - $REGION > region.bam
```

## Python Streaming Examples

### Stream and Process Line-by-Line

```python
import boto3
import gzip
from io import BytesIO

def process_fastq_from_s3(bucket, key):
    """Stream FASTQ from S3 and process line by line."""
    s3 = boto3.client('s3')

    # Get object
    response = s3.get_object(Bucket=bucket, Key=key)

    # Wrap in gzip if needed
    if key.endswith('.gz'):
        body = gzip.GzipFile(fileobj=response['Body'])
    else:
        body = response['Body']

    # Process line by line
    read_count = 0
    for line in body:
        line = line.decode('utf-8').strip()

        # FASTQ format: 4 lines per read
        # @SEQID
        # SEQUENCE
        # +
        # QUALITY

        if line.startswith('@'):
            read_count += 1

            # Process every 1000th read
            if read_count % 1000 == 0:
                print(f"Processed {read_count} reads...")

    print(f"Total reads: {read_count}")


# Usage
process_fastq_from_s3('my-bucket', 'sample.fastq.gz')
```

### Stream with Progress Bar

```python
import boto3
from tqdm import tqdm

def stream_with_progress(bucket, key, output_file):
    """Stream S3 object with progress bar."""
    s3 = boto3.client('s3')

    # Get object size
    head = s3.head_object(Bucket=bucket, Key=key)
    size = head['ContentLength']

    # Stream with progress
    response = s3.get_object(Bucket=bucket, Key=key)

    with open(output_file, 'wb') as f:
        with tqdm(total=size, unit='B', unit_scale=True) as pbar:
            for chunk in response['Body'].iter_chunks(chunk_size=1024*1024):
                f.write(chunk)
                pbar.update(len(chunk))


# Usage
stream_with_progress('my-bucket', 'large-file.bam', 'local-copy.bam')
```

## Advanced: Parallel Streaming

For tools that can process chunks independently:

```python
import boto3
from concurrent.futures import ThreadPoolExecutor
import math

def process_chunk(bucket, key, start, end):
    """Process a byte range of an S3 object."""
    s3 = boto3.client('s3')

    # Request specific byte range
    response = s3.get_object(
        Bucket=bucket,
        Key=key,
        Range=f'bytes={start}-{end}'
    )

    # Process chunk
    data = response['Body'].read()
    # ... do something with data ...

    return len(data)


def parallel_stream(bucket, key, chunk_size=10*1024*1024):
    """Stream S3 object in parallel chunks."""
    s3 = boto3.client('s3')

    # Get object size
    head = s3.head_object(Bucket=bucket, Key=key)
    size = head['ContentLength']

    # Calculate chunks
    num_chunks = math.ceil(size / chunk_size)
    chunks = [
        (i * chunk_size, min((i + 1) * chunk_size - 1, size - 1))
        for i in range(num_chunks)
    ]

    # Process in parallel
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = [
            executor.submit(process_chunk, bucket, key, start, end)
            for start, end in chunks
        ]

        total_bytes = sum(f.result() for f in futures)

    print(f"Processed {total_bytes} bytes")


# Usage
parallel_stream('my-bucket', 'large-file.txt')
```

## When to Use Streaming

### Good Use Cases ✅
- **Sequential processing** - Tools that read start to finish
- **Filtering** - Remove unwanted data before storing locally
- **Format conversion** - Convert on the fly
- **Quick inspection** - Check first N lines/records
- **Preprocessing** - Quality filtering, adapter trimming

### Bad Use Cases ❌
- **Random access** - Tools that seek around file
- **Multiple passes** - Tools that read file multiple times
- **Indexed access** - BAM region queries (need index)
- **Very large files** with poor network

## Tools That Work Well with Streaming

### ✅ Works Great
- `fastqc` - Quality control for FASTQ
- `fastp` - FASTQ filtering
- `cutadapt` - Adapter trimming
- `gzip/gunzip` - Compression
- `head/tail/grep` - Text processing
- `samtools view` (without region) - BAM to SAM

### ⚠️ Works with Limitations
- `samtools view` with region - Needs entire file or index
- `bwa mem` - Can stream FASTQ but not reference
- `bcftools` - Can stream but slow without index

### ❌ Doesn't Work
- `samtools sort` - Needs random access
- `picard MarkDuplicates` - Multiple passes
- IGV/genome browsers - Random access

## Optimization Tips

### 1. Use Larger Buffers

```python
# Bad: default buffer (small)
for line in response['Body']:
    process(line)

# Good: larger chunks
for chunk in response['Body'].iter_chunks(chunk_size=10*1024*1024):
    process(chunk)
```

### 2. Process as You Stream

```python
# Bad: accumulate then process
data = response['Body'].read()  # Downloads everything
process(data)

# Good: process incrementally
for chunk in response['Body'].iter_chunks():
    process(chunk)  # Process as data arrives
```

### 3. Use S3 Transfer Acceleration

```python
# Enable for faster downloads
s3 = boto3.client(
    's3',
    config=Config(s3={'use_accelerate_endpoint': True})
)
```

## Next Steps

1. Identify tools in your workflow that support stdin/stdout
2. Test streaming with small datasets
3. Measure time vs copying locally
4. Automate with scripts
5. For tools that don't support streaming, consider:
   - FUSE mounting (see `../fuse-mounting/`)
   - Tool modification (see `../../tool-modifications/`)

See `../genomics-pipeline/` for complete workflow examples using streaming.
