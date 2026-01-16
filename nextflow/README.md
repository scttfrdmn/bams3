# BAMS3 Nextflow Reference Pipeline

Production-ready Nextflow pipeline demonstrating cloud-native genomics workflows using BAMS3 format for alignment storage with direct S3 integration.

## Overview

This pipeline showcases BAMS3's key advantages in production genomics workflows:

- **Zero-copy alignment**: BWA streams directly to BAMS3 (no intermediate SAM/BAM)
- **Cloud-native storage**: Direct S3 upload during alignment
- **Selective queries**: Download only relevant chunks for variant calling
- **Simplified workflow**: 4-step traditional pipeline → 1-step BAMS3 pipeline
- **Cost-efficient**: 96% storage savings, 99.9% query cost savings

## Pipeline Workflow

```
FASTQ reads ──┐
              ├──> BWA mem ──> BAMS3 ──> S3
Reference ────┘         │
                        ├──> QC Stats
                        └──> Selective Region Extract ──> Variant Calling ──> VCF
```

### Processes

1. **INDEX_REFERENCE** - Index reference genome with BWA
2. **ALIGN_TO_BAMS3** - Zero-copy alignment (BWA | BAMS3)
3. **QC_STATS** - Generate alignment statistics
4. **CALL_VARIANTS_REGION** - Extract regions and call variants
5. **MERGE_VCFS** - Merge per-region VCFs
6. **MULTIQC** - Aggregate QC reports

## Quick Start

### Prerequisites

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Install BAMS3
# (See main repository README)

# Install dependencies
# - BWA
# - GATK or bcftools
# - samtools
# - tabix
# - jq
```

### Run Pipeline Locally

```bash
# Edit parameters
vim params.yaml

# Run pipeline
nextflow run main.nf -params-file params.yaml
```

### Run Pipeline on AWS Batch

```bash
# Use AWS profile
nextflow run main.nf \
    -params-file params.yaml \
    -profile aws \
    -work-dir s3://my-bucket/work
```

### Run Pipeline on HPC (SLURM)

```bash
# Use SLURM profile
nextflow run main.nf \
    -params-file params.yaml \
    -profile slurm
```

## Configuration

### Sample Sheet (`samples.csv`)

CSV file with columns:

```csv
sample_id,reads_r1,reads_r2
sample001,/data/sample001_R1.fq.gz,/data/sample001_R2.fq.gz
sample002,s3://data/sample002_R1.fq.gz,s3://data/sample002_R2.fq.gz
```

**Notes:**
- Supports local paths and S3 URIs
- Can mix local and S3 paths in same sheet
- Sample IDs must be unique

### Regions File (`regions.bed`)

Optional BED file for targeted variant calling:

```
chr1    1000000    2000000
chr2    500000     1500000
chr17   41196312   41277500
```

If not provided, whole-genome variant calling is performed.

### Parameters (`params.yaml`)

Key parameters:

```yaml
# Input
samples: "samples.csv"
reference: "/data/GRCh38.fa"
regions: "targets.bed"  # Optional

# Output (S3 or local)
output_dir: "s3://my-bucket/results"

# Resources
align_cpus: 16
align_memory: "32 GB"
sort_buffer: "28G"

# Variant caller
variant_caller: "gatk"  # or "bcftools"
```

## Execution Profiles

### `standard` - Local Execution

Default profile for local workstations:

```bash
nextflow run main.nf -params-file params.yaml
```

**Resource limits:**
- CPUs: 16
- Memory: 64 GB
- No containerization

### `aws` - AWS Batch

Cloud execution with auto-scaling:

```bash
nextflow run main.nf \
    -params-file params.yaml \
    -profile aws \
    -work-dir s3://my-bucket/work
```

**Features:**
- Auto-scaling compute (AWS Batch)
- S3 for all storage
- Spot instances supported
- Optimized instance types per process

**AWS Setup Required:**
1. AWS Batch compute environment
2. Job queue: `genomics-queue`
3. S3 bucket with appropriate permissions
4. IAM roles for Batch and EC2

### `slurm` - HPC Cluster

SLURM cluster execution:

```bash
nextflow run main.nf \
    -params-file params.yaml \
    -profile slurm
```

**Features:**
- Singularity containers
- SLURM resource allocation
- Shared filesystem

### `test` - Validation

Quick test with minimal data:

```bash
nextflow run main.nf -profile test
```

### `dev` - Development

Debug mode with detailed tracing:

```bash
nextflow run main.nf -profile dev -params-file params.yaml
```

Generates:
- `trace.txt` - Detailed resource usage
- `timeline.html` - Execution timeline
- `report.html` - Pipeline report

## Output Structure

```
output_dir/
├── bams3/
│   ├── sample001.bams3/
│   │   ├── _metadata.json
│   │   ├── _header.json
│   │   ├── _index/
│   │   └── data/
│   │       ├── chr1/
│   │       │   ├── 000000000-001048576.chunk
│   │       │   └── ...
│   │       └── chr2/
│   └── sample002.bams3/
├── vcf/
│   ├── sample001.vcf.gz
│   ├── sample001.vcf.gz.tbi
│   ├── sample002.vcf.gz
│   └── sample002.vcf.gz.tbi
├── qc/
│   ├── sample001.qc_report.txt
│   ├── sample002.qc_report.txt
│   └── multiqc_report.html
├── pipeline_report.html
├── execution_timeline.html
└── execution_trace.txt
```

## Cost Analysis

### Traditional BAM Workflow (30x WGS)

```
Alignment:
  BWA → SAM → BAM → sort → index
  Time: 6 hours @ c5.4xlarge ($0.68/hr) = $4.08
  Storage: 50 GB BAM + 50 GB SAM (temporary) = 100 GB

Per-sample cost:
  Compute: $4.08
  Storage: 50 GB × $0.023/mo = $1.15/month
  Download for analysis: 50 GB × $0.09/GB = $4.50 per analysis

Annual cost (10 analyses):
  Compute: $4.08
  Storage: $13.80 (12 months)
  Queries: $45.00 (10 × $4.50)
  Total: $62.88/sample/year
```

### BAMS3 Workflow (30x WGS)

```
Alignment:
  BWA | BAMS3 (zero-copy, direct to S3)
  Time: 6 hours @ c5.4xlarge ($0.68/hr) = $4.08
  Storage: 8 GB BAMS3 (no temporary files)

Per-sample cost:
  Compute: $4.08
  Storage: 8 GB × $0.023/mo = $0.18/month
  Selective region query: ~50 MB × $0.09/GB = $0.0045 per query

Annual cost (10 analyses):
  Compute: $4.08
  Storage: $2.16 (12 months)
  Queries: $0.045 (10 × $0.0045)
  Total: $6.29/sample/year

Savings: $56.59/sample/year (90% reduction)
```

### Cohort Economics (1,000 samples)

**Traditional BAM:**
- Storage: 50 TB @ $0.023/GB/mo = $1,150/month
- Annual storage cost: $13,800
- Query costs: Prohibitive (50 TB downloads)

**BAMS3:**
- Storage: 8 TB @ $0.023/GB/mo = $184/month
- Annual storage cost: $2,208
- Query costs: Minimal (selective chunks)

**Annual savings: $11,592 for 1,000 samples**

## Advanced Features

### Parallel Variant Calling

The pipeline automatically parallelizes variant calling across regions:

```
Sample ──> Region 1 ──┐
      ├──> Region 2   ├──> Merge ──> Final VCF
      └──> Region N ──┘
```

Each region extraction only downloads overlapping BAMS3 chunks.

### Custom Variant Callers

Add custom variant callers:

```groovy
process CALL_VARIANTS_CUSTOM {
    script:
    """
    bams3 to-bam "${bams3_path}" - --region ${region} | \
    my_variant_caller --input /dev/stdin --output output.vcf
    """
}
```

### Resource Optimization

#### Memory Optimization

```yaml
# High-memory sample (deep WGS)
align_memory: "64 GB"
sort_buffer: "60G"

# Low-memory sample (exome)
align_memory: "16 GB"
sort_buffer: "14G"
```

#### CPU Optimization

```yaml
# High-throughput (many samples)
align_cpus: 8
vcall_cpus: 2

# Low-latency (few samples)
align_cpus: 32
vcall_cpus: 8
```

### Spot Instance Usage (AWS)

Reduce costs by 70% using spot instances:

```groovy
// In nextflow.config
process {
    withName: 'ALIGN_TO_BAMS3' {
        queue = 'spot-genomics-queue'
        errorStrategy = 'retry'
        maxRetries = 5  // Handle spot interruptions
    }
}
```

## Troubleshooting

### Issue: Pipeline fails with "S3 access denied"

**Solution:** Ensure IAM role has required S3 permissions:

```json
{
  "Effect": "Allow",
  "Action": [
    "s3:PutObject",
    "s3:GetObject",
    "s3:ListBucket"
  ],
  "Resource": [
    "arn:aws:s3:::my-bucket/*",
    "arn:aws:s3:::my-bucket"
  ]
}
```

### Issue: "Out of memory" during alignment

**Solution:** Reduce sort buffer or increase memory allocation:

```yaml
align_memory: "64 GB"
sort_buffer: "58G"  # Leave 6GB for OS overhead
```

### Issue: Variant calling slow on large regions

**Solution:** Split regions into smaller chunks (< 10 MB each):

```bash
# Split large regions
awk '$3 - $2 > 10000000 {
    for(i=$2; i<$3; i+=10000000)
        print $1"\t"i"\t"(i+10000000<$3?i+10000000:$3)
}' regions.bed > regions_split.bed
```

### Issue: Nextflow caches old results

**Solution:** Clean work directory:

```bash
# Local
rm -rf work/

# S3
aws s3 rm s3://my-bucket/work/ --recursive
```

## Best Practices

### 1. Same-Region Deployment

**Critical for cost optimization:**
- Deploy Nextflow, AWS Batch, and S3 in same AWS region
- Data transfer within same region is FREE
- Cross-region transfer costs ~$0.02/GB

### 2. Resource Allocation

**Alignment:**
- CPU-bound: Maximize vCPUs (32-96)
- Memory: 2 GB per core + sort buffer
- Sort buffer: 80% of available memory

**Variant Calling:**
- Memory-bound: Optimize memory (16-64 GB)
- CPUs: 4-8 sufficient for most callers

### 3. Workflow Organization

**Sample Sheet Management:**
```bash
# Generate sample sheet from S3
aws s3 ls s3://data/fastq/ --recursive | \
    grep "_R1.fastq.gz" | \
    awk '{print $4}' | \
    sed 's/_R1.fastq.gz//' | \
    awk -F'/' '{print $NF",s3://data/fastq/"$NF"_R1.fastq.gz,s3://data/fastq/"$NF"_R2.fastq.gz"}' \
    > samples.csv
```

### 4. Monitoring

Enable CloudWatch logs for AWS Batch:

```groovy
aws {
    batch {
        logsGroup = '/aws/batch/genomics'
    }
}
```

### 5. Testing

Always validate with test profile first:

```bash
# Test with minimal data
nextflow run main.nf -profile test

# Then scale to full dataset
nextflow run main.nf -params-file params.yaml -profile aws
```

## Performance Benchmarks

### Small Dataset (100K reads)

| Metric | Traditional | BAMS3 | Improvement |
|--------|-------------|-------|-------------|
| Time | 5s | 5s | Comparable |
| Disk | 76 MB | 3.3 MB | 95.7% smaller |
| Steps | 4 | 1 | 4× simpler |

### Medium Dataset (1M reads)

| Metric | Traditional | BAMS3 | Improvement |
|--------|-------------|-------|-------------|
| Time | 45s | 50s | Comparable |
| Disk | 625 MB | 26 MB | 96.6% smaller |
| Steps | 4 | 1 | 4× simpler |

### Large Dataset (30x WGS, 1B reads)

| Metric | Traditional | BAMS3 | Improvement |
|--------|-------------|-------|-------------|
| Time | 6h | 6h | Comparable |
| Disk | 50 GB | 8 GB | 84% smaller |
| Query | 50 GB | 50 MB | 99.9% less data |
| Steps | 4 | 1 | 4× simpler |

## Integration Examples

### Cromwell/WDL

BAMS3 can be integrated into existing WDL workflows:

```wdl
task AlignToBams3 {
    input {
        File reads_r1
        File reads_r2
        File reference
        String output_uri
    }
    command {
        bwa mem -t 16 ${reference} ${reads_r1} ${reads_r2} | \
        bams3 convert --stdin ${output_uri} --workers 16
    }
}
```

### Snakemake

```python
rule align_bams3:
    input:
        r1 = "data/{sample}_R1.fq.gz",
        r2 = "data/{sample}_R2.fq.gz",
        ref = "refs/genome.fa"
    output:
        bams3 = "s3://bucket/{sample}.bams3"
    threads: 16
    shell:
        "bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | "
        "bams3 convert --stdin {output.bams3} --workers {threads}"
```

## References

- [BAMS3 Project](https://github.com/scttfrdmn/bams3)
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [AWS Batch Best Practices](https://docs.aws.amazon.com/batch/latest/userguide/best-practices.html)
- [S3 Integration Guide](../S3_INTEGRATION.md)
- [Performance Benchmarks](../benchmarks/README.md)

## Support

For issues or questions:
- GitHub Issues: https://github.com/scttfrdmn/bams3/issues
- Documentation: https://github.com/scttfrdmn/bams3/wiki

## License

MIT License - See main repository LICENSE file.
