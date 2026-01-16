# BAMS3 S3 Workflows (v0.3.0)

Zero-copy operations for cloud-native genomics workflows.

## Overview

BAMS3 v0.3.0+ supports direct S3 operations - no need to download entire datasets locally!

**Key Features:**
- ✅ Convert BAM → BAMS3 directly to S3
- ✅ Query BAMS3 datasets directly from S3
- ✅ Convert BAMS3 → BAM with S3 input
- ✅ Zero-copy operations (stream from/to S3)
- ✅ Works with AWS credentials (IAM roles, environment variables, profiles)

## Prerequisites

### AWS Credentials Setup

```bash
# Option 1: Environment variables
export AWS_ACCESS_KEY_ID="your-access-key"
export AWS_SECRET_ACCESS_KEY="your-secret-key"
export AWS_REGION="us-west-2"

# Option 2: AWS CLI profile
aws configure
# Enter your credentials when prompted

# Option 3: IAM role (for EC2, Lambda, etc.)
# No configuration needed - uses instance role automatically
```

### S3 Bucket Permissions

Your IAM user/role needs these permissions:

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:GetObject",
        "s3:PutObject",
        "s3:ListBucket"
      ],
      "Resource": [
        "arn:aws:s3:::your-bucket-name/*",
        "arn:aws:s3:::your-bucket-name"
      ]
    }
  ]
}
```

## Basic S3 Operations

### 1. Convert BAM to BAMS3 in S3

Upload while converting - true zero-copy!

```bash
# Convert local BAM → S3 BAMS3
bams3 convert input.bam s3://my-bucket/datasets/sample.bams3

# Convert with custom chunk size
bams3 convert --chunk-size 2M input.bam s3://my-bucket/datasets/sample.bams3

# Parallel conversion with 8 workers
bams3 convert --workers 8 input.bam s3://my-bucket/datasets/sample.bams3
```

**What happens:**
- Reads BAM file locally
- Compresses chunks
- Uploads chunks to S3 as they're created
- No temporary local storage needed!

### 2. Query BAMS3 Datasets from S3

Query without downloading the entire dataset!

```bash
# Query specific region
bams3 query s3://my-bucket/datasets/sample.bams3 chr1:1000000-2000000

# Count reads only (even faster!)
bams3 query s3://my-bucket/datasets/sample.bams3 chr1:1000000-2000000 --count

# Show first 10 reads
bams3 query s3://my-bucket/datasets/sample.bams3 chr20:10M-11M --show 10
```

**What happens:**
- Downloads only relevant chunks (not entire dataset!)
- For 1M chunk size, queries are 90%+ faster than downloading full dataset

**Cost Optimization:**
```bash
# Query chr1 region (1 MB)
# Only downloads ~2 chunks (2 GET requests, ~40 KB)
# Cost: $0.0000008 (2 × $0.0004 per 1000 requests)
bams3 query s3://my-bucket/datasets/sample.bams3 chr1:10000000-11000000 --count

# Compare to downloading full 50 GB BAM:
# Cost: ~$0.45 (data transfer)
# Time: 5-10 minutes vs 1-2 seconds
```

### 3. Get Statistics from S3

Instant statistics - no download needed!

```bash
# Get complete dataset statistics
bams3 stats s3://my-bucket/datasets/sample.bams3
```

**Output:**
```
===========================================
BAMS3 Dataset Statistics
===========================================

Format: bams3 v0.3.0
Source: s3://my-bucket/datasets/sample.bams3

Statistics:
  Total reads: 11,930,668
  Mapped reads: 11,930,668 (100.00%)
  Total bases: 303,477,269
  Mean coverage: 30.5x

Structure:
  Chunk size: 1048576 bp
  Total chunks: 245
  Compression: zstd

Storage:
  Total size: 531 MB
  Compression ratio: 6.1x
```

**What happens:**
- Downloads only metadata file (~15 KB)
- 3 GET requests total
- Instant results!

### 4. Convert BAMS3 from S3 to BAM

Download and convert in one step:

```bash
# S3 BAMS3 → local BAM
bams3 to-bam s3://my-bucket/datasets/sample.bams3 output.bam

# Extract specific region only
bams3 to-bam s3://my-bucket/datasets/sample.bams3 output.bam --region chr20:1-10000000
```

## Advanced S3 Workflows

### Population-Scale Analysis

Process 1000 samples without downloading!

```bash
# List all samples in S3
aws s3 ls s3://cohort-bucket/samples/ | grep .bams3 > samples.txt

# Query all samples for specific region (parallel)
cat samples.txt | parallel -j 50 \
  'bams3 query s3://cohort-bucket/samples/{} chr20:10M-11M --count > {/.}_count.txt'

# Aggregate results
awk '{sum+=$1} END {print "Total reads:", sum}' *_count.txt
```

**Benefits:**
- No local storage needed
- 50 parallel queries = 50x faster
- Each query downloads only ~2 chunks (~40 KB)
- Total data transfer: ~2 MB (vs 50 TB for full downloads!)

### Distributed Processing with Spark on EMR

```python
from pyspark import SparkContext, SparkConf
import subprocess

conf = SparkConf().setAppName("BAMS3 S3 Population Analysis")
sc = SparkContext(conf=conf)

# List of S3 paths (1000 samples)
samples = [f"s3://cohort/sample_{i}.bams3" for i in range(1, 1001)]

# Parallelize across cluster
samples_rdd = sc.parallelize(samples, numSlices=1000)

def query_sample(s3_path):
    """Query directly from S3 - no download needed!"""
    result = subprocess.check_output([
        "bams3", "query", s3_path, "chr20:10M-11M", "--count"
    ])
    return (s3_path, int(result.strip()))

# Process all samples in parallel
results = samples_rdd.map(query_sample).collect()

# Aggregate
total_reads = sum(count for _, count in results)
print(f"Total reads across 1000 samples: {total_reads:,}")
```

**Expected Performance:**
- 100-node cluster: processes 1000 samples in ~2 minutes
- Each node queries 10 samples
- Total data transfer: ~40 MB per node
- Compare to downloading all BAMs: 50 TB transfer, hours of time!

### AWS Batch Array Jobs

Massively parallel queries across cohorts:

```bash
# job.json - Batch job definition
{
  "jobDefinitionName": "bams3-s3-query",
  "type": "container",
  "containerProperties": {
    "image": "genomics/bams3:latest",
    "vcpus": 2,
    "memory": 4096,
    "command": [
      "bams3", "query",
      "s3://cohort/Ref::sample",
      "Ref::region",
      "--count"
    ],
    "jobRoleArn": "arn:aws:iam::ACCOUNT:role/BatchJobRole"
  }
}

# Submit array job (1000 samples)
aws batch submit-job \
  --job-name cohort-variant-calling \
  --job-queue genomics-queue \
  --job-definition bams3-s3-query \
  --array-properties size=1000 \
  --parameters \
    sample="sample_${AWS_BATCH_JOB_ARRAY_INDEX}.bams3" \
    region="chr20:10M-11M"
```

**Benefits:**
- Auto-scales to 1000 concurrent jobs
- Pay only for compute time (~$0.08 per sample-hour)
- No data transfer between S3 and compute (same region)
- Results written back to S3

### Lambda Integration

Serverless queries for real-time applications:

```python
# lambda_function.py
import boto3
import subprocess
import json

def lambda_handler(event, context):
    """
    Serverless BAMS3 query from S3
    """
    s3_path = event['s3_path']
    region = event['region']

    # Query directly from S3
    result = subprocess.run([
        "bams3", "query", s3_path, region, "--count"
    ], capture_output=True, text=True)

    count = int(result.stdout.strip())

    return {
        'statusCode': 200,
        'body': json.dumps({
            'sample': s3_path,
            'region': region,
            'read_count': count
        })
    }
```

**Invocation:**
```bash
aws lambda invoke \
  --function-name bams3-query \
  --payload '{"s3_path": "s3://bucket/sample.bams3", "region": "chr1:1M-2M"}' \
  response.json

cat response.json
# {"sample": "s3://bucket/sample.bams3", "region": "chr1:1M-2M", "read_count": 1523}
```

## S3 Performance Optimization

### 1. Choose Optimal Chunk Size

Trade-off between parallelism and request costs:

```bash
# Small chunks (256K) - better parallelism, more requests
bams3 convert --chunk-size 256K input.bam s3://bucket/sample.bams3
# Queries: Download 4-8 chunks per region
# Cost: Higher S3 request costs

# Large chunks (4M) - fewer requests, less parallelism
bams3 convert --chunk-size 4M input.bam s3://bucket/sample.bams3
# Queries: Download 1-2 chunks per region
# Cost: Lower S3 request costs

# Recommended: 1M chunks (default)
# Good balance for most workflows
```

### 2. Use S3 Transfer Acceleration

For large uploads from distant regions:

```bash
# Enable transfer acceleration on bucket
aws s3api put-bucket-accelerate-configuration \
  --bucket my-bucket \
  --accelerate-configuration Status=Enabled

# Use accelerated endpoint
export AWS_S3_USE_ACCELERATE=true
bams3 convert input.bam s3://my-bucket/sample.bams3
# Up to 50% faster uploads from distant regions!
```

### 3. Same-Region Processing

Avoid data transfer costs:

```bash
# Run bams3 on EC2 instance in same region as S3 bucket
# Example: bucket in us-west-2, EC2 in us-west-2

# Check your region
aws configure get region  # Should match bucket region

# Query without data transfer charges
bams3 query s3://us-west-2-bucket/sample.bams3 chr1:1M-2M
# Free data transfer within same region!
```

### 4. Parallel Queries with Chunked Regions

Split large regions for parallel processing:

```bash
# Query entire chr1 (249 Mb) - slow
bams3 query s3://bucket/sample.bams3 chr1:1-249250621 --count

# Split into 10 chunks - 10x faster!
for i in {0..9}; do
  start=$((i * 24925062 + 1))
  end=$(((i + 1) * 24925062))
  bams3 query s3://bucket/sample.bams3 chr1:$start-$end --count &
done
wait
```

## Cost Analysis

### Storage Costs

S3 Standard storage: $0.023 per GB/month

```bash
# 50 GB BAM file
# → 8.3 GB BAMS3 (6x compression)
# Monthly cost: $0.19

# 1000 samples × 50 GB = 50 TB
# → 8.3 TB BAMS3
# Monthly cost: $190 (vs $1,150 for uncompressed BAM)
```

### Request Costs

S3 GET requests: $0.0004 per 1,000 requests

```bash
# Query 1 region (downloads 2 chunks)
# Cost: $0.0000008

# Query 1000 samples × 100 regions = 100,000 queries
# 100,000 × 2 = 200,000 GET requests
# Cost: $0.08

# Compare to downloading all BAMs:
# 1000 × 50 GB = 50 TB transfer
# Cost: $4,500 (data transfer) + 10 hours
```

### Intelligent-Tiering

For infrequently accessed datasets:

```bash
# Enable Intelligent-Tiering
aws s3api put-bucket-intelligent-tiering-configuration \
  --bucket my-bucket \
  --id default-config \
  --intelligent-tiering-configuration '{
    "Id": "default-config",
    "Status": "Enabled",
    "Tierings": [
      {
        "Days": 90,
        "AccessTier": "ARCHIVE_ACCESS"
      }
    ]
  }'

# Datasets not accessed for 90 days automatically move to cheaper storage
# Savings: Up to 70% on storage costs
```

## Multi-Region Replication

Replicate datasets for global access:

```bash
# Enable versioning (required for replication)
aws s3api put-bucket-versioning \
  --bucket source-bucket \
  --versioning-configuration Status=Enabled

# Create replication rule
aws s3api put-bucket-replication \
  --bucket source-bucket \
  --replication-configuration file://replication.json

# replication.json
{
  "Role": "arn:aws:iam::ACCOUNT:role/ReplicationRole",
  "Rules": [
    {
      "Status": "Enabled",
      "Priority": 1,
      "Destination": {
        "Bucket": "arn:aws:s3:::destination-bucket",
        "ReplicationTime": {
          "Status": "Enabled",
          "Time": {
            "Minutes": 15
          }
        }
      }
    }
  ]
}
```

## Troubleshooting

### Error: "Access Denied"

Check your IAM permissions:

```bash
# Test S3 access
aws s3 ls s3://your-bucket/

# Verify IAM user
aws sts get-caller-identity

# Check bucket policy
aws s3api get-bucket-policy --bucket your-bucket
```

### Error: "No credentials found"

Set up AWS credentials:

```bash
# Option 1: AWS CLI
aws configure

# Option 2: Environment variables
export AWS_ACCESS_KEY_ID="..."
export AWS_SECRET_ACCESS_KEY="..."

# Option 3: Instance profile (EC2/ECS/Lambda)
# Automatically detected - no configuration needed
```

### Slow S3 Queries

Optimize for your workload:

```bash
# 1. Use same region for processing
# 2. Enable S3 Transfer Acceleration
# 3. Use larger chunk sizes (fewer requests)
# 4. Query smaller regions when possible
```

## Best Practices

### 1. Dataset Organization

```
s3://genomics-bucket/
├── projects/
│   ├── project-001/
│   │   ├── sample-A.bams3/
│   │   ├── sample-B.bams3/
│   │   └── sample-C.bams3/
│   └── project-002/
│       └── samples/
├── reference/
│   └── hg38.bams3/
└── analysis-results/
    └── variant-calls/
```

### 2. Metadata Tagging

```bash
# Tag datasets with metadata
aws s3api put-object-tagging \
  --bucket my-bucket \
  --key datasets/sample.bams3/_metadata.json \
  --tagging 'TagSet=[
    {Key=project,Value=genomics-001},
    {Key=sample-id,Value=NA12878},
    {Key=tissue,Value=blood},
    {Key=coverage,Value=30x}
  ]'
```

### 3. Lifecycle Policies

```bash
# Archive old datasets after 1 year
aws s3api put-bucket-lifecycle-configuration \
  --bucket my-bucket \
  --lifecycle-configuration file://lifecycle.json

# lifecycle.json
{
  "Rules": [
    {
      "Status": "Enabled",
      "Transitions": [
        {
          "Days": 365,
          "StorageClass": "GLACIER"
        }
      ]
    }
  ]
}
```

## Summary

**BAMS3 + S3 enables:**

✅ **Zero-copy queries** - download only relevant chunks
✅ **Massive cost savings** - 90%+ reduction in data transfer
✅ **Instant scalability** - query 1000s of samples in parallel
✅ **Cloud-native workflows** - no local storage needed

**Perfect for:**
- Population-scale studies
- Real-time variant calling
- Serverless genomics pipelines
- Multi-region collaborations

---

**Next Steps:**
- [Parallel Workflows](./parallel_workflows.md) - Distributed processing
- [Variant Calling](./variant_calling.md) - Pipeline integration
- [Cost Optimization](./cost_optimization.md) - Minimize cloud costs
