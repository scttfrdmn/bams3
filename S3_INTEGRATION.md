# S3 Direct Integration

## Overview

BAMS3 provides native S3 integration, allowing you to stream alignment data directly to S3 without local storage. This enables true cloud-native genomics workflows with zero local disk usage.

## Architecture

### Storage Abstraction

BAMS3 uses a pluggable storage backend that automatically detects S3 vs local filesystem:

```go
type Storage interface {
    ReadFile(path string) ([]byte, error)
    WriteFile(path string, data []byte) error
    List(prefix string) ([]string, error)
    Exists(path string) (bool, error)
    // ...
}
```

- **Local paths** (`/path/to/file.bams3`) → `LocalStorage`
- **S3 URIs** (`s3://bucket/prefix`) → `S3Storage`

No code changes needed - just use `s3://` prefix in your paths!

## Setup

### 1. AWS Credentials

Configure AWS credentials using one of these methods:

**Environment Variables:**
```bash
export AWS_ACCESS_KEY_ID=your_access_key
export AWS_SECRET_ACCESS_KEY=your_secret_key
export AWS_REGION=us-east-1
```

**AWS CLI Configuration:**
```bash
aws configure
```

**IAM Role (EC2/ECS):**
Attach IAM role with S3 permissions to your compute instance.

### 2. Required IAM Permissions

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:PutObject",
        "s3:GetObject",
        "s3:ListBucket",
        "s3:HeadObject"
      ],
      "Resource": [
        "arn:aws:s3:::your-bucket/*",
        "arn:aws:s3:::your-bucket"
      ]
    }
  ]
}
```

## Usage

### Streaming Conversion to S3

**Zero-copy pipeline:**
```bash
bwa mem ref.fa reads_R1.fq reads_R2.fq | \
    bams3 convert --stdin s3://my-bucket/samples/sample001.bams3
```

**With custom settings:**
```bash
bwa mem -t 32 ref.fa reads_R1.fq reads_R2.fq | \
    bams3 convert --stdin s3://my-bucket/wgs/sample001.bams3 \
        --workers 32 \
        --sort-buffer 64G \
        --chunk-size 1M
```

### Querying from S3

**Region query (selective chunk download):**
```bash
bams3 query s3://my-bucket/samples/sample001.bams3 chr1:1000000-2000000
```

**Statistics:**
```bash
bams3 stats s3://my-bucket/samples/sample001.bams3
```

### Cross-Region Usage

```bash
# Specify region explicitly if bucket is in different region
AWS_REGION=us-west-2 bams3 convert --stdin s3://west-coast-bucket/sample.bams3
```

## Performance Optimization

### Same-Region Deployment

**Critical for cost optimization:**
- Deploy compute in same region as S3 bucket
- Data transfer within same region is **FREE**
- Cross-region transfer: ~$0.02/GB

**Example:**
```bash
# Bucket in us-east-1, EC2 in us-east-1 → Free transfer
# Bucket in us-east-1, EC2 in us-west-2 → $0.02/GB
```

### Multipart Upload

BAMS3 automatically uses S3 multipart upload for chunks:
- Part size: 10MB (optimized for chunk size)
- Concurrent uploads: 3 parts at a time
- Automatic retry on failure

### Selective Downloads

Range requests minimize data transfer:
```bash
# Only downloads chunks overlapping chr1:1M-2M
bams3 query s3://bucket/sample.bams3 chr1:1000000-2000000

# Traditional BAM would download entire file first
```

## Cost Analysis

### Storage Costs

**S3 Standard:**
- $0.023/GB/month
- First 50 TB

**Example:** 30x WGS (50GB BAM → 8GB BAMS3)
- Traditional: 50GB × $0.023 = **$1.15/month**
- BAMS3: 8GB × $0.023 = **$0.18/month**
- **Savings: 84%**

### Data Transfer

**Within Same Region:**
- Upload (PUT): **Free**
- Download (GET): **Free**
- Total: **$0**

**Cross-Region:**
- Upload: Free
- Download: $0.02/GB
- Avoid if possible!

### Request Costs

**S3 Standard:**
- PUT requests: **Free** (included with Standard tier)
- GET requests: $0.0004 per 1,000 requests

**Example query (5 chunks):**
- 5 GET requests × $0.0004/1000 = **$0.000002**
- Effectively free

### Total Cost Comparison

**Traditional BAM Workflow (50GB file):**
```
Download BAM:     50GB × $0.09/GB      = $4.50
Storage:          50GB × $0.023/mo     = $1.15/month
Per-query:        50GB × $0.09/GB      = $4.50
First analysis:   $4.50 + $1.15/month
Each subsequent:  $4.50
```

**BAMS3 Workflow (8GB BAMS3):**
```
Upload BAMS3:     Free (same region)
Storage:          8GB × $0.023/mo      = $0.18/month
Per-query:        ~50MB × $0.09/GB     = $0.0045 (selective)
First analysis:   $0 + $0.18/month
Each subsequent:  $0.0045
```

**Savings:**
- First analysis: **91% reduction**
- Storage: **84% reduction**
- Per-query: **99.9% reduction**

## S3 Structure

BAMS3 creates this structure in S3:

```
s3://my-bucket/samples/sample001.bams3/
├── _metadata.json           # Dataset metadata
├── _header.json             # SAM header
├── _index/
│   └── spatial.json         # Spatial index for queries
└── data/
    ├── chr1/
    │   ├── 000000000-001048576.chunk
    │   ├── 001048576-002097152.chunk
    │   └── ...
    ├── chr2/
    │   └── ...
    └── unmapped/
        └── 000000000-001048576.chunk
```

Each chunk is independently queryable without downloading others.

## Production Workflows

### AWS Batch

```bash
#!/bin/bash
# Task array for parallel sample processing

SAMPLE_ID=$AWS_BATCH_JOB_ARRAY_INDEX
BUCKET=s3://genomics-data

# Download FASTQ from S3
aws s3 cp ${BUCKET}/fastq/${SAMPLE_ID}_R1.fq.gz - | gunzip > R1.fq &
aws s3 cp ${BUCKET}/fastq/${SAMPLE_ID}_R2.fq.gz - | gunzip > R2.fq &
wait

# Align and stream directly to S3
bwa mem -t ${AWS_BATCH_JOB_NUM_CPUS} \
    /ref/GRCh38.fa \
    R1.fq R2.fq | \
    bams3 convert --stdin ${BUCKET}/aligned/${SAMPLE_ID}.bams3 \
        --workers ${AWS_BATCH_JOB_NUM_CPUS} \
        --sort-buffer 32G

# No cleanup needed - no local BAM files created!
```

### Nextflow

```nextflow
process ALIGN {
    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2)
    path(reference)

    output:
    val("s3://bucket/${sample_id}.bams3")

    script:
    """
    bwa mem -t ${task.cpus} ${reference} ${reads_r1} ${reads_r2} | \
        bams3 convert --stdin s3://genomics-bucket/aligned/${sample_id}.bams3 \
            --workers ${task.cpus} \
            --sort-buffer ${task.memory.toGiga()}G
    """
}

process VARIANT_CALL {
    input:
    val(bams3_uri)
    val(region)

    output:
    path("${region}.vcf")

    script:
    """
    bams3 to-bam ${bams3_uri} - --region ${region} | \
        gatk HaplotypeCaller \
            -I /dev/stdin \
            -R /ref/GRCh38.fa \
            -O ${region}.vcf \
            -L ${region}
    """
}
```

## Troubleshooting

### Authentication Issues

```bash
# Test AWS credentials
aws sts get-caller-identity

# Test S3 access
aws s3 ls s3://your-bucket/
```

### Permission Denied

Ensure IAM role/user has required permissions:
- `s3:PutObject` - Upload chunks
- `s3:GetObject` - Download for queries
- `s3:ListBucket` - List chunks
- `s3:HeadObject` - Check existence

### Slow Uploads

Check if you're in the same region:
```bash
aws s3api get-bucket-location --bucket your-bucket
aws ec2 describe-availability-zones --query 'AvailabilityZones[0].Region'
```

If different regions, either:
1. Create bucket in same region as compute
2. Accept cross-region transfer costs

### Connection Timeouts

For large uploads, increase timeout:
```bash
export AWS_DEFAULT_CONNECTION_TIMEOUT=300
export AWS_DEFAULT_READ_TIMEOUT=300
```

## Best Practices

1. **Same-Region Everything**: Compute, S3 bucket, and downstream tools in same region
2. **Lifecycle Policies**: Archive old BAMS3 files to S3 Glacier for long-term storage
3. **Bucket Versioning**: Enable for protection against accidental deletion
4. **Server-Side Encryption**: Enable SSE-S3 or SSE-KMS for compliance
5. **VPC Endpoints**: Use S3 VPC endpoints to avoid internet gateway costs
6. **Monitoring**: Enable S3 access logs and CloudWatch metrics

## Advanced Features

### S3 VPC Endpoint (No NAT Gateway)

```bash
# Create VPC endpoint for S3
aws ec2 create-vpc-endpoint \
    --vpc-id vpc-xxx \
    --service-name com.amazonaws.us-east-1.s3 \
    --route-table-ids rtb-xxx

# All S3 traffic now stays within AWS network (free)
```

### Server-Side Encryption

```bash
# KMS encryption (automatic)
bams3 convert --stdin s3://bucket/sample.bams3
# BAMS3 uploads with standard S3 encryption
```

### S3 Intelligent-Tiering

```bash
# Automatic cost optimization for variable access patterns
aws s3api put-bucket-intelligent-tiering-configuration \
    --bucket your-bucket \
    --id bams3-tiering \
    --intelligent-tiering-configuration '{
        "Id": "bams3-tiering",
        "Status": "Enabled",
        "Tierings": [{
            "Days": 90,
            "AccessTier": "ARCHIVE_ACCESS"
        }]
    }'
```

## References

- [AWS S3 Pricing](https://aws.amazon.com/s3/pricing/)
- [S3 Transfer Acceleration](https://docs.aws.amazon.com/AmazonS3/latest/userguide/transfer-acceleration.html)
- [S3 VPC Endpoints](https://docs.aws.amazon.com/vpc/latest/privatelink/vpc-endpoints-s3.html)
- [AWS Batch Best Practices](https://docs.aws.amazon.com/batch/latest/userguide/best-practices.html)
