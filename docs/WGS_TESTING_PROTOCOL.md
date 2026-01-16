# BAMS3 Real WGS Testing Protocol

Testing BAMS3 with production-scale whole genome sequencing data to validate performance, costs, and tool compatibility.

## Test Dataset

### Primary: 1000 Genomes NA12878 (Genome in a Bottle)

**Why NA12878:**
- Well-characterized reference sample
- ~30x coverage (realistic for clinical/research)
- Multiple sequencing platforms available
- Gold standard for validation
- Available on AWS S3 (no download costs!)

**Dataset Location:**
```bash
# 1000 Genomes on AWS
s3://1000genomes/phase3/data/NA12878/sequence_read/

# Genome in a Bottle (NIST)
s3://giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/
```

**Dataset Specifications:**
- Sample: NA12878 (CEPH individual)
- Coverage: ~30x
- Read length: 150bp paired-end
- Platform: Illumina HiSeq
- Expected reads: ~1 billion
- Expected BAM size: ~100GB
- Reference: GRCh38

### Alternative: Synthetic Data (Controlled Testing)

For specific stress tests, generate synthetic data:

```bash
# Install wgsim
brew install bwa samtools

# Generate synthetic WGS (30x, chr1 only for quick test)
samtools faidx GRCh38.fa chr1 > chr1.fa
wgsim -N 50000000 -1 150 -2 150 chr1.fa reads_R1.fq reads_R2.fq

# This creates ~50M reads (30x coverage of chr1)
```

---

## Test Environment Setup

### AWS Configuration

**EC2 Instance:**
- Instance type: `c5.9xlarge` or `c5.12xlarge`
  - 36-48 vCPUs (for parallel alignment)
  - 72-96 GB RAM (for sort buffer)
  - EBS: 500GB gp3 (for temporary files)
- AMI: Amazon Linux 2 or Ubuntu 22.04
- Region: `us-east-1` (same as 1000 Genomes bucket)

**S3 Bucket:**
- Create bucket in `us-east-1` (same region = free data transfer!)
- Name: `bams3-testing-<your-name>`
- Enable versioning (for safety)
- Lifecycle: Delete after 30 days (to control costs)

```bash
aws s3 mb s3://bams3-testing-yourusername --region us-east-1
```

**IAM Permissions:**
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
        "arn:aws:s3:::1000genomes/*",
        "arn:aws:s3:::bams3-testing-*/*",
        "arn:aws:s3:::bams3-testing-*"
      ]
    }
  ]
}
```

### Software Installation

```bash
# Update system
sudo yum update -y  # Amazon Linux
# or: sudo apt update && sudo apt upgrade -y  # Ubuntu

# Install dependencies
sudo yum install -y gcc make git wget bzip2

# Install Go (1.21+)
wget https://go.dev/dl/go1.21.5.linux-amd64.tar.gz
sudo tar -C /usr/local -xzf go1.21.5.linux-amd64.tar.gz
export PATH=$PATH:/usr/local/go/bin
echo 'export PATH=$PATH:/usr/local/go/bin' >> ~/.bashrc

# Install BWA
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17 && make
sudo cp bwa /usr/local/bin/

# Install samtools
wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
tar -xjf samtools-1.18.tar.bz2
cd samtools-1.18 && ./configure && make
sudo make install

# Install BAMS3
git clone https://github.com/scttfrdmn/bams3.git
cd bams3/bams3-go
go build -o bams3 ./cmd/bams3
sudo cp bams3 /usr/local/bin/

# Verify installations
bwa 2>&1 | head -3
samtools --version
bams3 --version
```

### Download Reference Genome

```bash
# Download GRCh38 (no-alt analysis set recommended)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Unzip
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38.fa

# Index for BWA (takes ~1 hour)
bwa index GRCh38.fa

# Index for samtools
samtools faidx GRCh38.fa
```

---

## Test Scenarios

### Test 1: Full WGS Alignment and Conversion

**Goal**: Validate BAMS3 handles 1B reads, measures actual performance

**Commands**:
```bash
# Set variables
REFERENCE="GRCh38.fa"
SAMPLE="NA12878"
FASTQ_R1="s3://1000genomes/phase3/data/NA12878/sequence_read/SRR622461_1.filt.fastq.gz"
FASTQ_R2="s3://1000genomes/phase3/data/NA12878/sequence_read/SRR622461_2.filt.fastq.gz"
OUTPUT_BAMS3="s3://bams3-testing-yourusername/${SAMPLE}.bams3"
WORKERS=32
SORT_BUFFER="60G"

# Stream FASTQ from S3, align, convert to BAMS3
time (
  aws s3 cp $FASTQ_R1 - | gunzip > /tmp/R1.fq &
  aws s3 cp $FASTQ_R2 - | gunzip > /tmp/R2.fq &
  wait

  bwa mem -t $WORKERS $REFERENCE /tmp/R1.fq /tmp/R2.fq 2> bwa.log | \
    bams3 convert --stdin $OUTPUT_BAMS3 \
      --workers $WORKERS \
      --sort-buffer $SORT_BUFFER \
      2>&1 | tee bams3_convert.log
)
```

**Metrics to Record**:
```bash
# From logs
grep "Total reads processed" bams3_convert.log
grep "Throughput" bams3_convert.log
grep "Elapsed time" bams3_convert.log

# Memory usage (run in separate terminal during conversion)
watch -n 5 "ps aux | grep bams3 | head -1"

# Disk usage
df -h
du -sh /tmp/*

# S3 costs (after completion)
aws s3 ls s3://bams3-testing-yourusername/${SAMPLE}.bams3/ --recursive --human-readable --summarize
```

**Expected Results**:
- Conversion time: 2-4 hours (depends on instance type)
- Peak memory: ~65GB (60GB sort buffer + 5GB overhead)
- Disk spill: Should NOT trigger with 60GB buffer
- Output size: ~8-12GB (vs ~100GB BAM)
- Reads: ~1 billion (validate with stats)

---

### Test 2: Region Extraction and GATK Integration

**Goal**: Validate selective access and tool compatibility

**Test Regions**:
- Small region (100kb): BRCA1 gene
- Medium region (1Mb): Multi-gene locus
- Large region (10Mb): Entire chromosome arm
- Whole chromosome: chr22 (smallest)

**Commands**:
```bash
# Test BRCA1 region (chr17:41196312-41277500, 81kb)
echo "Testing BRCA1 region extraction..."
time bams3 to-bam $OUTPUT_BAMS3 - --region chr17:41196312-41277500 > brca1.bam
samtools view -c brca1.bam  # Count reads
samtools index brca1.bam

# Variant calling with GATK
time gatk HaplotypeCaller \
  -I brca1.bam \
  -R $REFERENCE \
  -O brca1_variants.vcf \
  -L chr17:41196312-41277500

# Test streaming (no intermediate file)
echo "Testing streaming to GATK..."
time (
  bams3 to-bam $OUTPUT_BAMS3 - --region chr17:41196312-41277500 | \
  gatk HaplotypeCaller \
    -I /dev/stdin \
    -R $REFERENCE \
    -O brca1_stream.vcf \
    -L chr17:41196312-41277500
)

# Compare VCFs (should be identical)
diff brca1_variants.vcf brca1_stream.vcf
```

**S3 Data Transfer Metrics**:
```bash
# Enable S3 access logging
aws s3api put-bucket-logging \
  --bucket bams3-testing-yourusername \
  --bucket-logging-status file://logging.json

# After queries, analyze logs
aws s3 sync s3://bams3-testing-yourusername/logs/ ./s3_logs/
python analyze_s3_costs.py s3_logs/
```

**Expected Results**:
- BRCA1 extraction: <5 seconds, ~5MB download
- Variant calling: Works identically to traditional BAM
- VCFs match: Streaming vs file-based identical
- Data transfer: 99%+ reduction vs full BAM download

---

### Test 3: Nextflow Pipeline Validation

**Goal**: Validate production workflow

**Setup**:
```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Create sample sheet
cat > samples.csv << EOF
sample_id,reads_r1,reads_r2
NA12878,s3://1000genomes/.../_1.filt.fastq.gz,s3://1000genomes/.../_2.filt.fastq.gz
EOF

# Create regions file (for targeted calling)
cat > regions.bed << EOF
chr17	41196312	41277500
chr13	32889611	32973805
chr19	1220373	1234567
EOF
```

**Run Pipeline**:
```bash
cd bams3/nextflow

nextflow run main.nf \
  --samples ../../samples.csv \
  --reference ../../GRCh38.fa \
  --regions ../../regions.bed \
  --output_dir s3://bams3-testing-yourusername/results \
  --bams3_dir s3://bams3-testing-yourusername/bams3 \
  --align_cpus 32 \
  --align_memory "64 GB" \
  --sort_buffer "60G" \
  -profile standard \
  -with-report report.html \
  -with-timeline timeline.html
```

**Expected Results**:
- Pipeline completes without errors
- BAMS3 created in S3
- VCFs generated for each region
- MultiQC report successful
- Resource usage within limits

---

### Test 4: Cost Analysis

**Goal**: Measure actual AWS costs, validate projections

**Cost Components to Track**:

1. **Compute (EC2)**:
   ```bash
   # c5.9xlarge in us-east-1: $1.53/hour
   # Expected: 2-4 hours for WGS = $3-6
   ```

2. **Storage (S3)**:
   ```bash
   # BAMS3 size: ~10GB
   # Cost: 10GB × $0.023/GB/month = $0.23/month

   # Compare to traditional: 100GB BAM
   # Cost: 100GB × $0.023/GB/month = $2.30/month
   # Savings: 90%
   ```

3. **Data Transfer**:
   ```bash
   # Within us-east-1: FREE
   # Upload to S3: FREE
   # Download from 1000 Genomes: FREE (public dataset)

   # Query costs (region extraction):
   # Download 5MB chunk vs 100GB full BAM
   # 5MB × $0.09/GB = $0.00045 (vs $9.00 for full BAM)
   # Savings: 99.995%
   ```

4. **S3 Requests**:
   ```bash
   # PUT: Free with S3 Standard
   # GET: $0.0004 per 1000 requests
   # Chunk queries: ~10 GETs = $0.000004
   # Effectively free
   ```

**Total Cost for Test**:
- EC2: $3-6 (one-time)
- S3 storage: $0.23/month
- Queries: <$0.01
- **Total: ~$6 for complete validation**

**Compare to Traditional**:
- EC2: $3-6 (same)
- S3 storage: $2.30/month (10x more)
- Queries: $9 per full download
- **10 queries: $96 vs <$1 with BAMS3**

---

## Test Execution Checklist

### Pre-Test
- [ ] EC2 instance launched (c5.9xlarge or larger)
- [ ] S3 bucket created in us-east-1
- [ ] IAM permissions configured
- [ ] Software installed (BWA, samtools, BAMS3, GATK)
- [ ] Reference genome downloaded and indexed
- [ ] Monitoring configured (CloudWatch, S3 logging)

### During Test
- [ ] Monitor memory usage (`htop`, `watch ps aux`)
- [ ] Monitor disk usage (`df -h`, `du -sh /tmp`)
- [ ] Capture timing for each step
- [ ] Log all commands and outputs
- [ ] Take CloudWatch snapshots

### Post-Test
- [ ] Validate read counts (BAMS3 vs traditional)
- [ ] Verify GATK compatibility (VCFs valid)
- [ ] Calculate actual costs from AWS billing
- [ ] Document performance metrics
- [ ] Update documentation with real numbers
- [ ] Clean up resources (terminate EC2, delete test files)

---

## Metrics Collection Template

Create `test_results.md`:

```markdown
# BAMS3 WGS Test Results

**Date**: YYYY-MM-DD
**Instance**: c5.9xlarge (36 vCPU, 72 GB RAM)
**Dataset**: NA12878, ~30x WGS
**Reference**: GRCh38

## Conversion Performance

| Metric | Value |
|--------|-------|
| Total reads | X billion |
| Conversion time | X hours |
| Throughput | X K reads/sec |
| Peak memory | X GB |
| Disk spill | Yes/No, X GB |
| Output size (BAMS3) | X GB |
| Output size (BAM est.) | X GB |
| Compression ratio | X% |

## Region Extraction

| Region | Size | Extraction Time | Data Downloaded | Reads Extracted |
|--------|------|-----------------|-----------------|-----------------|
| BRCA1 | 81kb | X sec | X MB | X reads |
| 1Mb region | 1Mb | X sec | X MB | X reads |
| chr22 | 51Mb | X sec | X MB | X reads |

## Tool Compatibility

| Tool | Test | Result | Notes |
|------|------|--------|-------|
| GATK HC | BRCA1 calling | Pass/Fail | VCF valid |
| Streaming | GATK stdin | Pass/Fail | Works |
| samtools | flagstat | Pass/Fail | Stats match |

## Cost Analysis

| Component | Cost | Notes |
|-----------|------|-------|
| EC2 compute | $X | X hours |
| S3 storage | $X/month | X GB |
| Data transfer | $X | S3 GETs |
| **Total** | **$X** | One-time + monthly |

**Comparison to Traditional**:
- Storage savings: X%
- Query cost savings: X%
- Total savings: $X/year

## Issues Encountered

- [ ] None OR list issues

## Recommendations

- Update documentation with actual numbers
- Adjust default parameters if needed
- Known limitations to document
```

---

## Next Steps After Testing

1. **Update Documentation**:
   - Replace projected numbers with actual measurements
   - Add case study to repository
   - Update README with real benchmarks

2. **Blog Post**:
   - Write detailed case study
   - Include commands, results, costs
   - Share on bioinformatics forums

3. **GitHub Release**:
   - Create v1.0.0 with confidence
   - Include test results in release notes

4. **Community**:
   - Share results for feedback
   - Identify early adopters
   - Gather feature requests

---

## Troubleshooting

### Out of Memory
```bash
# Reduce sort buffer
--sort-buffer 40G  # instead of 60G

# Or use larger instance
c5.12xlarge  # 96 GB RAM
```

### Disk Spill Triggered
```bash
# This is OK! It means the safety net works
# Check spill logs:
grep "spill" bams3_convert.log

# Validate output still correct:
bams3 stats output.bams3
```

### S3 Upload Slow
```bash
# Check you're in same region:
aws s3api get-bucket-location --bucket bams3-testing-yourusername
# Should show: us-east-1

# If not, recreate bucket in us-east-1
```

### GATK Errors
```bash
# Check BAM is valid:
samtools quickcheck region.bam && echo "OK" || echo "FAIL"

# Validate header:
samtools view -H region.bam

# Check GATK version:
gatk --version  # Should be 4.0+
```

---

This protocol provides comprehensive validation of BAMS3 at production scale. Expected total time: 1 day (mostly waiting for conversion).
