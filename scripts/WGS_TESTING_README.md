# BAMS3 WGS Production Testing Scripts

Complete AWS-based testing suite for BAMS3 validation with real whole-genome sequencing data.

## Overview

These scripts automate the provisioning, testing, and reporting of BAMS3 performance with production-scale genomics data (NA12878, 30x WGS coverage, ~1 billion reads).

**Total Test Time**: ~6 hours
**Total Cost**: ~$11 (EC2 + S3 + data transfer)
**Region**: us-west-2
**AWS Profile**: `aws`

## Prerequisites

### Local Requirements
- AWS CLI configured with profile named `aws`
- SSH key pair in us-west-2
- Security group allowing SSH (port 22)
- IAM role: `BAMS3-Testing-Role` with permissions:
  - S3 read from `s3://1000genomes/*`
  - S3 read/write to `s3://bams3-testing-*/*`
  - CloudWatch logs write

### AWS Resources
The test suite will create:
- S3 bucket: `s3://bams3-testing-${USER}`
- EC2 instance: c5.9xlarge (36 vCPU, 72GB RAM, 500GB disk)
- Estimated cost: ~$11 for full test suite

## Scripts

| Script | Purpose | Time | Where to Run |
|--------|---------|------|--------------|
| `aws-setup.sh` | Create S3 bucket and launch EC2 | 5 min | **Local** |
| `ec2-setup-userdata.sh` | Install dependencies, download reference | 90 min | **EC2** |
| `run-chr22-test.sh` | Chr22 validation test (~5M reads) | 30 min | **EC2** |
| `run-full-wgs-test.sh` | Full WGS test (~1B reads) | 4 hours | **EC2** |
| `run-gatk-integration-test.sh` | GATK variant calling test | 30 min | **EC2** |
| `generate-test-report.sh` | Generate comprehensive markdown report | 5 min | **EC2** |

## Execution Workflow

### Phase 1: AWS Infrastructure Setup (Local)

```bash
cd scripts

# Launch EC2 instance and create S3 bucket
./aws-setup.sh <your-key-name> <your-security-group-id>

# Example:
# ./aws-setup.sh my-ec2-key sg-0123456789abcdef0

# Note the instance IP from output
```

**Output**: Instance ID, public IP, S3 bucket name saved to `instance-info.txt`

### Phase 2: EC2 Environment Setup (On EC2 Instance)

```bash
# From local machine, copy setup script
scp -i ~/.ssh/<key-name>.pem scripts/ec2-setup-userdata.sh ec2-user@<instance-ip>:~/

# SSH to instance
ssh -i ~/.ssh/<key-name>.pem ec2-user@<instance-ip>

# Run setup (installs Go, BWA, samtools, BAMS3, downloads GRCh38)
chmod +x ec2-setup-userdata.sh
sudo ./ec2-setup-userdata.sh

# This takes ~90 minutes (reference indexing is slow)
```

### Phase 3: Copy Test Scripts to EC2

```bash
# From local machine
scp -i ~/.ssh/<key-name>.pem scripts/run-*.sh ec2-user@<instance-ip>:/data/
scp -i ~/.ssh/<key-name>.pem scripts/generate-test-report.sh ec2-user@<instance-ip>:/data/

# SSH back to instance
ssh -i ~/.ssh/<key-name>.pem ec2-user@<instance-ip>
cd /data
chmod +x *.sh
```

### Phase 4: Run Tests (On EC2 Instance)

#### Step 1: Chr22 Validation Test (30 minutes)

```bash
./run-chr22-test.sh
```

**What it does**:
- Downloads chr22 FASTQ subset (~5M reads)
- Runs BWA → BAMS3 pipeline
- Uploads to S3
- Tests region extraction (BRCA1)
- Validates success criteria

**Success criteria**:
- ✓ Time < 30 minutes
- ✓ No errors in logs
- ✓ Region extraction works
- ✓ S3 upload succeeds

#### Step 2: Full WGS Test (4 hours)

```bash
./run-full-wgs-test.sh
```

**What it does**:
- Downloads full NA12878 FASTQ (~100GB, ~1B reads)
- Runs complete BWA → BAMS3 pipeline
- Uploads BAMS3 file to S3
- Collects comprehensive metrics
- Uploads results to S3

**Success criteria**:
- ✓ Completes in 2-4 hours
- ✓ ~1 billion reads processed
- ✓ ~8-12GB BAMS3 output (vs ~100GB BAM)
- ✓ Memory < 72GB
- ✓ No disk spill

#### Step 3: GATK Integration Test (30 minutes)

```bash
./run-gatk-integration-test.sh
```

**What it does**:
- Installs GATK (if needed)
- Streams BRCA1, BRCA2, TP53 regions from S3 BAMS3
- Runs GATK HaplotypeCaller
- Validates VCF outputs
- Uploads results to S3

**Success criteria**:
- ✓ All VCFs validate
- ✓ Variant counts reasonable
- ✓ Streaming works correctly

#### Step 4: Generate Report (5 minutes)

```bash
./generate-test-report.sh
```

**What it does**:
- Consolidates all test results
- Generates comprehensive markdown report
- Uploads to S3
- Commits to GitHub (if credentials available)

**Output**: `WGS_TEST_RESULTS.md` with full analysis and cost projections

## Expected Results

### Storage Efficiency
- Traditional BAM: ~100 GB
- BAMS3: ~10 GB
- **Savings**: 90%

### Cost Efficiency (1000 queries)
- Traditional: $900 (download full BAM each time)
- BAMS3: $4.50 (selective region access)
- **Savings**: 99.5%

### Performance
- Full WGS conversion: 4 hours
- Region extraction: <30 seconds (vs 32 minutes for traditional)
- Memory usage: <65GB (no disk spill)

## Viewing Results

### On EC2 Instance

```bash
# Chr22 validation
cat /data/results/chr22_test/chr22_results.txt

# Full WGS
cat /data/results/full_wgs/full_wgs_results.txt

# GATK integration
cat /data/results/gatk_integration/gatk_integration_results.txt

# Final report
cat /data/results/WGS_TEST_RESULTS.md
```

### From S3

```bash
# List all results
aws --profile aws s3 ls s3://bams3-testing-${USER}/results/

# Download final report
aws --profile aws s3 cp s3://bams3-testing-${USER}/results/WGS_TEST_RESULTS.md .

# Download logs
aws --profile aws s3 sync s3://bams3-testing-${USER}/results/logs/ ./logs/
```

## Troubleshooting

### EC2 Setup Fails

**Problem**: Reference download or indexing fails

**Solution**:
```bash
# Check disk space
df -h

# Restart setup from reference download
cd /data/references
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip *.fna.gz
mv *.fna GRCh38.fa
bwa index GRCh38.fa
samtools faidx GRCh38.fa
```

### BAMS3 Conversion Fails

**Problem**: Out of memory or disk space

**Solution**:
```bash
# Check memory
free -h

# Check disk
df -h

# Reduce sort buffer (in test script)
--sort-buffer 40G  # Instead of 60G

# Clear temporary files
rm -rf /tmp/*
```

### S3 Upload Fails

**Problem**: AWS credentials or permissions

**Solution**:
```bash
# Verify AWS profile
aws --profile aws sts get-caller-identity

# Check S3 bucket access
aws --profile aws s3 ls s3://bams3-testing-${USER}/

# Manually upload results
aws --profile aws s3 cp /data/results/ s3://bams3-testing-${USER}/results/ --recursive
```

## Cleanup

### Keep Results, Terminate EC2

```bash
# From local machine
aws --profile aws ec2 terminate-instances --instance-ids <instance-id> --region us-west-2

# S3 bucket and results will remain
```

### Full Cleanup

```bash
# Delete S3 bucket
aws --profile aws s3 rb s3://bams3-testing-${USER} --force --region us-west-2

# Terminate instance
aws --profile aws ec2 terminate-instances --instance-ids <instance-id> --region us-west-2
```

## Cost Breakdown

| Item | Cost | Notes |
|------|------|-------|
| EC2 (6 hours) | $9.18 | c5.9xlarge @ $1.53/hour |
| S3 storage | $0.23/month | ~10GB @ $0.023/GB/month |
| Data transfer | $2.00 | Cross-region (us-east-1 → us-west-2) |
| **Total** | **$11.20** | One-time cost |

## Next Steps

After successful testing:

1. **Document Results**: Review `WGS_TEST_RESULTS.md`
2. **Share Findings**: Upload report to GitHub
3. **Iterate**: Use results to optimize BAMS3 parameters
4. **Expand**: Implement VCFS3 and migration tools

## Support

For issues or questions:
- GitHub Issues: https://github.com/scttfrdmn/aws-direct-s3/issues
- Documentation: `/docs` directory
- Email: scott@friedman.com
