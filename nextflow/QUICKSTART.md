# BAMS3 Nextflow Pipeline - Quick Start Guide

Get the pipeline running in 5 minutes.

## Prerequisites

Install required tools:

```bash
# Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# BAMS3
# See: https://github.com/scttfrdmn/bams3

# Dependencies
# - BWA: brew install bwa  (or apt-get install bwa)
# - samtools: brew install samtools
# - bcftools: brew install bcftools
# - GATK: Download from https://github.com/broadinstitute/gatk/releases
# - jq: brew install jq
```

## Option 1: Test Run (1 minute)

Validate pipeline with minimal test data:

```bash
cd nextflow

# Generate test data
cd test_data
bash generate_test_data.sh
cd ..

# Run test pipeline
nextflow run main.nf -profile test
```

**Expected output:**
```
✓ Reference indexed
✓ Aligned test_sample (2000 reads)
✓ Generated QC report
✓ Called variants (2 regions)

Results in: test_results/
```

## Option 2: Local Run (5 minutes)

Run with your own data:

### Step 1: Prepare Sample Sheet

Create `my_samples.csv`:

```csv
sample_id,reads_r1,reads_r2
sample001,/data/fastq/sample001_R1.fq.gz,/data/fastq/sample001_R2.fq.gz
sample002,/data/fastq/sample002_R1.fq.gz,/data/fastq/sample002_R2.fq.gz
```

### Step 2: Configure Parameters

Create `my_params.yaml`:

```yaml
samples: "my_samples.csv"
reference: "/data/references/GRCh38.fa"
output_dir: "results"

# Optional: target regions
regions: "targets.bed"

# Resources
align_cpus: 8
align_memory: "16 GB"
sort_buffer: "14G"
```

### Step 3: Run Pipeline

```bash
nextflow run main.nf -params-file my_params.yaml
```

### Step 4: Check Results

```bash
ls -lh results/
├── bams3/
│   ├── sample001.bams3/
│   └── sample002.bams3/
├── vcf/
│   ├── sample001.vcf.gz
│   └── sample002.vcf.gz
└── qc/
    ├── sample001.qc_report.txt
    ├── sample002.qc_report.txt
    └── multiqc_report.html
```

## Option 3: AWS Batch (10 minutes setup)

Run on AWS with auto-scaling:

### Step 1: AWS Setup

```bash
# Configure AWS credentials
aws configure

# Create S3 bucket
aws s3 mb s3://my-genomics-bucket

# Upload reference
aws s3 cp reference.fa s3://my-genomics-bucket/reference/
```

### Step 2: Configure for AWS

Create `aws_params.yaml`:

```yaml
samples: "s3://my-genomics-bucket/samples.csv"
reference: "s3://my-genomics-bucket/reference/GRCh38.fa"
output_dir: "s3://my-genomics-bucket/results"

align_cpus: 32
align_memory: "64 GB"
sort_buffer: "58G"
```

### Step 3: Run on AWS Batch

```bash
nextflow run main.nf \
    -params-file aws_params.yaml \
    -profile aws \
    -work-dir s3://my-genomics-bucket/work
```

### Step 4: Monitor Progress

```bash
# Check Nextflow output
tail -f .nextflow.log

# Check AWS Batch console
# https://console.aws.amazon.com/batch/
```

## Common Scenarios

### Exome Sequencing

```yaml
samples: "exome_samples.csv"
reference: "GRCh38.fa"
regions: "exome_targets.bed"
output_dir: "s3://bucket/exome_results"

align_cpus: 16
align_memory: "32 GB"
variant_caller: "gatk"
```

### Whole Genome Sequencing

```yaml
samples: "wgs_samples.csv"
reference: "GRCh38.fa"
regions: null  # Whole genome
output_dir: "s3://bucket/wgs_results"

align_cpus: 32
align_memory: "64 GB"
sort_buffer: "58G"
variant_caller: "gatk"
```

### Targeted Panel

```yaml
samples: "panel_samples.csv"
reference: "GRCh38.fa"
regions: "gene_panel_targets.bed"
output_dir: "results/panel"

align_cpus: 8
align_memory: "16 GB"
variant_caller: "bcftools"  # Faster for small regions
```

## Troubleshooting

### Pipeline hangs at alignment

**Check:** CPU/memory allocation

```bash
# Monitor resources
top  # or htop

# Reduce if needed
align_cpus: 4
align_memory: "8 GB"
sort_buffer: "6G"
```

### S3 upload fails

**Check:** AWS credentials and permissions

```bash
# Test access
aws s3 ls s3://my-bucket/

# Check IAM permissions
aws iam get-user
```

### Out of memory

**Solution:** Increase memory or reduce sort buffer

```yaml
align_memory: "32 GB"
sort_buffer: "28G"  # Leave 4GB for OS
```

### VCF files empty

**Check:** Regions overlap with aligned data

```bash
# Check alignment statistics
bams3 stats results/bams3/sample001.bams3

# Verify regions are valid
head regions.bed
```

## Next Steps

- **Production Setup**: See [README.md](README.md) for full documentation
- **AWS Optimization**: Configure Batch queues and spot instances
- **Custom Analysis**: Add custom processes to pipeline
- **Cost Analysis**: Review [S3_INTEGRATION.md](../S3_INTEGRATION.md)

## Examples

### Resume Failed Run

```bash
nextflow run main.nf -params-file params.yaml -resume
```

### Run Specific Samples

```bash
# Create filtered sample sheet
head -n 1 samples.csv > subset_samples.csv
grep "sample001\|sample005" samples.csv >> subset_samples.csv

# Run with subset
nextflow run main.nf \
    --samples subset_samples.csv \
    --reference reference.fa
```

### Custom Resource Allocation

```bash
nextflow run main.nf \
    --samples samples.csv \
    --reference reference.fa \
    --align_cpus 24 \
    --align_memory "48 GB" \
    --sort_buffer "44G"
```

## Getting Help

- **Documentation**: [README.md](README.md)
- **Issues**: https://github.com/scttfrdmn/bams3/issues
- **Examples**: See `test_data/` directory

## Performance Tips

1. **Use same region**: Deploy compute, Nextflow, and S3 in same AWS region
2. **Optimize resources**: CPU for alignment, memory for variant calling
3. **Split regions**: Break large regions into smaller chunks (< 10 MB)
4. **Use spot instances**: Save 70% on AWS compute costs
5. **Resume on failure**: Use `-resume` to continue from last successful step

## Cost Estimates

### Local Workstation
- Free (uses your hardware)
- Storage: Local disk

### AWS Batch (per 30x WGS sample)
- Compute: ~$4-6 (6 hours @ c5.4xlarge)
- Storage: $0.18/month (8 GB BAMS3)
- Queries: $0.0045 per analysis (selective chunks)

**Total: ~$6/sample (vs ~$63 traditional BAM)**
