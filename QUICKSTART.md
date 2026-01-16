# Quick Start Guide

Get started with S3-native genomics workflows in 5 minutes.

## Prerequisites

- AWS account with S3 access
- AWS credentials configured (`aws configure`)
- Some genomics data in S3

## Immediate Solutions (No Code Changes)

### Option 1: Mount S3 as Filesystem (Easiest)

**Best for:** Existing tools, no modifications needed

```bash
# Install mountpoint-for-amazon-s3
# macOS:
brew install mountpoint-s3

# Linux:
wget https://s3.amazonaws.com/mountpoint-s3-release/latest/x86_64/mount-s3.deb
sudo apt install ./mount-s3.deb

# Mount your bucket
mkdir ~/s3-data
mount-s3 YOUR-BUCKET-NAME ~/s3-data

# Use tools normally
samtools view ~/s3-data/sample.bam | head
bcftools query -f '%CHROM\t%POS\n' ~/s3-data/variants.vcf.gz | head

# Unmount when done
umount ~/s3-data
```

**Performance tips:**
```bash
# Add caching for better performance
mount-s3 YOUR-BUCKET \
    ~/s3-data \
    --cache /tmp/s3-cache \
    --max-cache-size 50000  # 50GB cache
```

### Option 2: Stream for Sequential Processing

**Best for:** FASTQ processing, filtering, format conversion

```bash
# Quality check FASTQ without downloading
aws s3 cp s3://bucket/sample.fastq.gz - | fastqc stdin -o ./qc-results

# Filter and upload result
aws s3 cp s3://bucket/raw.fastq.gz - | \
    gunzip | \
    fastp --stdin --stdout --qualified_quality_phred 20 | \
    gzip | \
    aws s3 cp - s3://bucket/filtered.fastq.gz
```

## Test Your Setup

### 1. Check Credentials

```bash
aws s3 ls
# Should list your buckets
```

### 2. Upload Test Data

```bash
# Create a small test file
echo "@read1
ACGTACGTACGT
+
IIIIIIIIIIII" > test.fastq

# Upload to S3
aws s3 cp test.fastq s3://YOUR-BUCKET/test/test.fastq
```

### 3. Test FUSE Mount

```bash
mkdir ~/s3-mount
mount-s3 YOUR-BUCKET ~/s3-mount

# Read the file through mount
cat ~/s3-mount/test/test.fastq

# Should see your test data
```

### 4. Test Streaming

```bash
# Stream from S3
aws s3 cp s3://YOUR-BUCKET/test/test.fastq - | wc -l

# Should output: 4 (lines in FASTQ)
```

## Real Example: Complete Workflow

### Traditional (Slow) Approach

```bash
# ❌ DON'T DO THIS - wastes time and storage

# Download data (10+ minutes for large files)
aws s3 cp s3://bucket/sample1_R1.fastq.gz ./
aws s3 cp s3://bucket/sample1_R2.fastq.gz ./
aws s3 cp s3://bucket/reference.fa ./

# Process
bwa mem reference.fa sample1_R1.fastq.gz sample1_R2.fastq.gz > aligned.sam
samtools view -b aligned.sam > aligned.bam

# Upload results
aws s3 cp aligned.bam s3://bucket/results/

# Clean up
rm -f sample1_R1.fastq.gz sample1_R2.fastq.gz reference.fa aligned.sam aligned.bam
```

### Modern (Fast) Approach

```bash
# ✅ DO THIS - instant start, no wasted storage

# Mount S3
mount-s3 YOUR-BUCKET ~/s3-data

# Process directly (reference cached on first access)
bwa mem \
    ~/s3-data/reference.fa \
    ~/s3-data/sample1_R1.fastq.gz \
    ~/s3-data/sample1_R2.fastq.gz | \
    samtools view -b - > /tmp/aligned.bam

# Upload result
aws s3 cp /tmp/aligned.bam s3://bucket/results/

# Cleanup (only local result, inputs stay in S3)
rm /tmp/aligned.bam

# Unmount
umount ~/s3-data
```

**Time saved:** 10-20 minutes per sample
**Storage saved:** Hundreds of GB

## Common Patterns

### Pattern 1: Read-Only Analysis

```bash
# Mount S3
mount-s3 --read-only data-bucket ~/data

# Run analyses
for bam in ~/data/*.bam; do
    samtools flagstat $bam > $(basename $bam .bam).stats
done

# Upload results
aws s3 sync . s3://results-bucket/stats/ --include "*.stats"
```

### Pattern 2: Streaming Pipeline

```bash
#!/bin/bash
# stream-pipeline.sh

INPUT="s3://input-bucket/sample.fastq.gz"
OUTPUT="s3://output-bucket/sample.bam"

aws s3 cp $INPUT - | \
    gunzip | \
    bwa mem reference.fa - | \
    samtools view -b - | \
    samtools sort - | \
    aws s3 cp - $OUTPUT
```

### Pattern 3: Hybrid (Cache Hot Data)

```bash
# Download frequently accessed data (e.g., reference genome)
aws s3 cp s3://refs/hg38.fa ./reference/

# Mount sample data
mount-s3 samples-bucket ~/samples

# Process with local reference, remote samples
for sample in ~/samples/*.fastq.gz; do
    bwa mem ./reference/hg38.fa $sample > $(basename $sample .fastq.gz).sam
done
```

## Benchmarking Your Workflow

Test to see which approach is fastest for your data:

```bash
#!/bin/bash
# benchmark.sh

SAMPLE="s3://bucket/sample.bam"

echo "Test 1: Copy then process"
time {
    aws s3 cp $SAMPLE ./sample.bam
    samtools flagstat ./sample.bam > stats.txt
    rm sample.bam
}

echo "Test 2: Mount and process"
time {
    mkdir -p /tmp/mnt
    mount-s3 bucket /tmp/mnt
    samtools flagstat /tmp/mnt/sample.bam > stats.txt
    umount /tmp/mnt
}

echo "Test 3: Stream (if supported)"
time {
    aws s3 cp $SAMPLE - | samtools flagstat - > stats.txt
}
```

Run this with your real data to see which works best.

## When You're Ready for More

### Short-term improvements:
1. **Optimize mount settings** - Tune cache size, thread count
2. **Reorganize data** - Partition by chromosome, sample, etc.
3. **Use streaming where possible** - Avoid unnecessary downloads

### Long-term improvements:
4. **Modify tools** - Add native S3 support to your key tools
5. **Optimize formats** - Convert to cloud-optimized formats
6. **Contribute upstream** - Share improvements with community

See the full documentation for details:
- [FUSE Mounting](examples/fuse-mounting/README.md)
- [Streaming](examples/streaming/README.md)
- [Tool Modification](docs/tool-modification-guide.md)
- [Format Optimization](docs/format-optimization.md)

## Troubleshooting

### "Permission denied" mounting S3
```bash
# Check credentials
aws sts get-caller-identity

# Verify bucket access
aws s3 ls s3://YOUR-BUCKET/
```

### "Transport endpoint not connected"
```bash
# Unmount forcefully
fusermount -u ~/s3-data  # Linux
umount -f ~/s3-data      # macOS

# Remount
mount-s3 YOUR-BUCKET ~/s3-data
```

### Slow performance
```bash
# Add caching
mount-s3 YOUR-BUCKET ~/s3-data \
    --cache /path/to/fast/storage \
    --max-cache-size 100000 \
    --thread-count 16
```

### Tool doesn't work with mounted S3
- Check if tool needs random access (BAM with regions)
- Verify index files exist (.bai, .tbi, etc.)
- Consider streaming approach or tool modification instead

## Get Help

- Read the docs in `docs/`
- Check examples in `examples/`
- Open an issue in this repository
- AWS Mountpoint docs: https://github.com/awslabs/mountpoint-s3

## Next Steps

1. ✅ Get FUSE mounting working
2. ✅ Test with your real workflows
3. ✅ Measure performance improvements
4. ⏭️ Explore format optimization
5. ⏭️ Consider tool modifications for long-term solution
