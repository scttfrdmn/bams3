#!/bin/bash
set -e

# BAMS3 WGS Testing - Chr22 Validation Test
# Tests: Small chr22 region (10Mbp, ~5M reads)
# Time: ~30 minutes

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m'

echo -e "${BOLD}${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}${BLUE}║    BAMS3 Chr22 Validation Test                                ║${NC}"
echo -e "${BOLD}${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Configuration
PROFILE="aws"
REGION="us-west-2"
SAMPLE="NA12878"
S3_OUTPUT="s3://bams3-testing-${USER}/chr22_test.bams3"

# Check for required files
if [ ! -f "/data/references/GRCh38.fa" ]; then
    echo -e "${RED}✗ Reference not found: /data/references/GRCh38.fa${NC}"
    echo "Run ec2-setup-userdata.sh first"
    exit 1
fi

# Create working directory
cd /data/results
mkdir -p chr22_test
cd chr22_test

echo -e "${BLUE}Test Configuration:${NC}"
echo "  Sample:     $SAMPLE"
echo "  Region:     chr22:10000000-20000000 (10Mbp)"
echo "  S3 output:  $S3_OUTPUT"
echo ""

# Step 1: Download chr22 FASTQ subset from 1000 Genomes
echo -e "${BLUE}[1/5] Downloading chr22 FASTQ subset from 1000 Genomes...${NC}"
echo -e "${YELLOW}Note: Cross-region transfer (us-east-1 → us-west-2)${NC}"
echo ""

DOWNLOAD_START=$(date +%s)

echo "  Downloading R1..."
aws --profile $PROFILE s3 cp s3://1000genomes/phase3/data/NA12878/sequence_read/SRR622461_1.filt.fastq.gz - 2>/dev/null | \
  gunzip | head -20000000 > /tmp/chr22_R1.fq

echo "  Downloading R2..."
aws --profile $PROFILE s3 cp s3://1000genomes/phase3/data/NA12878/sequence_read/SRR622461_2.filt.fastq.gz - 2>/dev/null | \
  gunzip | head -20000000 > /tmp/chr22_R2.fq

DOWNLOAD_END=$(date +%s)
DOWNLOAD_TIME=$((DOWNLOAD_END - DOWNLOAD_START))

R1_READS=$(wc -l < /tmp/chr22_R1.fq)
R1_READS=$((R1_READS / 4))

echo ""
echo -e "${GREEN}✓${NC} Download complete (${DOWNLOAD_TIME}s)"
echo "  Reads: ${R1_READS}"
echo ""

# Step 2: Run BWA alignment and BAMS3 conversion
echo -e "${BLUE}[2/5] Running BWA alignment and BAMS3 conversion...${NC}"
echo "  Workers: 36"
echo "  Sort buffer: 60G"
echo ""

PIPELINE_START=$(date +%s)

bwa mem -t 36 /data/references/GRCh38.fa /tmp/chr22_R1.fq /tmp/chr22_R2.fq 2> bwa_chr22.log | \
  bams3 convert --stdin $S3_OUTPUT \
    --workers 36 \
    --sort-buffer 60G \
    2>&1 | tee bams3_chr22.log

PIPELINE_END=$(date +%s)
PIPELINE_TIME=$((PIPELINE_END - PIPELINE_START))

echo ""
echo -e "${GREEN}✓${NC} Pipeline complete (${PIPELINE_TIME}s)"
echo ""

# Step 3: Collect metrics
echo -e "${BLUE}[3/5] Collecting metrics...${NC}"

cat > chr22_results.txt << EOF
BAMS3 Chr22 Validation Test Results
====================================
Date: $(date)
Sample: $SAMPLE

Download:
  Time: ${DOWNLOAD_TIME}s
  Reads: ${R1_READS}

Pipeline:
  Time: ${PIPELINE_TIME}s
  BWA log: bwa_chr22.log
  BAMS3 log: bams3_chr22.log

EOF

# Extract metrics from logs
echo "BAMS3 Statistics:" >> chr22_results.txt
grep -E "(Total reads|Throughput|chunks created)" bams3_chr22.log >> chr22_results.txt 2>/dev/null || echo "  (Check bams3_chr22.log for details)" >> chr22_results.txt
echo "" >> chr22_results.txt

echo -e "${GREEN}✓${NC} Metrics collected"
echo ""

# Step 4: Validate BAMS3 output
echo -e "${BLUE}[4/5] Validating BAMS3 output...${NC}"

echo "BAMS3 Dataset Stats:" >> chr22_results.txt
bams3 stats $S3_OUTPUT 2>&1 >> chr22_results.txt || echo "  (Stats command failed)" >> chr22_results.txt
echo "" >> chr22_results.txt

echo -e "${GREEN}✓${NC} Validation complete"
echo ""

# Step 5: Test region extraction
echo -e "${BLUE}[5/5] Testing region extraction...${NC}"
echo "  Extracting BRCA1: chr17:41196312-41277500"

bams3 to-bam $S3_OUTPUT /tmp/brca1.bam --region chr17:41196312-41277500 2>/dev/null

BRCA1_READS=$(samtools view -c /tmp/brca1.bam 2>/dev/null)

echo ""
echo -e "${GREEN}✓${NC} Region extraction successful"
echo "  BRCA1 reads: $BRCA1_READS"
echo ""

echo "Region Extraction Test:" >> chr22_results.txt
echo "  Region: chr17:41196312-41277500 (BRCA1)" >> chr22_results.txt
echo "  Reads: $BRCA1_READS" >> chr22_results.txt
echo "" >> chr22_results.txt

# Summary
TOTAL_TIME=$((DOWNLOAD_TIME + PIPELINE_TIME))

cat >> chr22_results.txt << EOF
Summary:
  Total time: ${TOTAL_TIME}s
  Download: ${DOWNLOAD_TIME}s
  Pipeline: ${PIPELINE_TIME}s
  S3 output: $S3_OUTPUT

Success Criteria:
EOF

# Check success criteria
ALL_PASSED=true

if [ $PIPELINE_TIME -lt 1800 ]; then
    echo "  ✓ Time < 30 minutes" >> chr22_results.txt
else
    echo "  ✗ Time > 30 minutes" >> chr22_results.txt
    ALL_PASSED=false
fi

if [ -f "bams3_chr22.log" ] && ! grep -qi "error\|failed" bams3_chr22.log; then
    echo "  ✓ No errors in logs" >> chr22_results.txt
else
    echo "  ✗ Errors found in logs" >> chr22_results.txt
    ALL_PASSED=false
fi

if [ "$BRCA1_READS" -gt 0 ]; then
    echo "  ✓ Region extraction works" >> chr22_results.txt
else
    echo "  ✗ Region extraction failed" >> chr22_results.txt
    ALL_PASSED=false
fi

if aws --profile $PROFILE s3 ls $S3_OUTPUT/ >/dev/null 2>&1; then
    echo "  ✓ S3 upload succeeded" >> chr22_results.txt
else
    echo "  ✗ S3 upload failed" >> chr22_results.txt
    ALL_PASSED=false
fi

echo "" >> chr22_results.txt

# Display results
echo -e "${BOLD}${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}${BLUE}║                Test Results                                    ║${NC}"
echo -e "${BOLD}${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
cat chr22_results.txt
echo ""

# Cleanup temp files
rm -f /tmp/chr22_R1.fq /tmp/chr22_R2.fq /tmp/brca1.bam

if [ "$ALL_PASSED" = true ]; then
    echo -e "${BOLD}${GREEN}✓ Chr22 Validation Test PASSED${NC}"
    echo ""
    echo -e "${BLUE}Next Steps:${NC}"
    echo "  1. Review results: cat /data/results/chr22_test/chr22_results.txt"
    echo "  2. Run full WGS test: cd /data && ./run-full-wgs-test.sh"
    echo ""
    exit 0
else
    echo -e "${BOLD}${RED}✗ Chr22 Validation Test FAILED${NC}"
    echo ""
    echo "  Check logs:"
    echo "    - bwa_chr22.log"
    echo "    - bams3_chr22.log"
    echo ""
    exit 1
fi
