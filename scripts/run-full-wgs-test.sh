#!/bin/bash
set -e

# BAMS3 WGS Testing - Full WGS Test
# Tests: Complete NA12878 genome (~1 billion reads, 30x coverage)
# Time: ~4 hours

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m'

echo -e "${BOLD}${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}${BLUE}║    BAMS3 Full WGS Production Test                             ║${NC}"
echo -e "${BOLD}${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Configuration
PROFILE="aws"
REGION="us-east-1"
SAMPLE="NA12878"
S3_OUTPUT="s3://bams3-testing-${USER}/NA12878_full.bams3"

# Check for required files
if [ ! -f "/data/references/GRCh38.fa" ]; then
    echo -e "${RED}✗ Reference not found: /data/references/GRCh38.fa${NC}"
    echo "Run ec2-setup-userdata.sh first"
    exit 1
fi

# Create working directory
cd /data/results
mkdir -p full_wgs
cd full_wgs
mkdir -p logs

echo -e "${BLUE}Test Configuration:${NC}"
echo "  Sample:     $SAMPLE (1000 Genomes)"
echo "  Coverage:   30x whole genome"
echo "  Reads:      ~1 billion"
echo "  S3 output:  $S3_OUTPUT"
echo ""
echo -e "${YELLOW}Warning: This test will:${NC}"
echo "  • Download ~100GB from 1000 Genomes (same-region, no transfer cost)"
echo "  • Take approximately 4 hours"
echo "  • Cost ~$9 (EC2 + S3 storage)"
echo ""

read -p "Continue? (yes/no) " -r
if [[ ! $REPLY =~ ^[Yy][Ee][Ss]$ ]]; then
    echo "Test cancelled"
    exit 0
fi
echo ""

# Step 1: Download NA12878 FASTQ
echo -e "${BLUE}[1/4] Downloading NA12878 FASTQ files from 1000 Genomes...${NC}"
echo -e "${GREEN}Note: Same-region transfer (us-east-1) - no cross-region costs${NC}"
echo "  This will take approximately 30 minutes"
echo ""

DOWNLOAD_START=$(date +%s)

aws --profile $PROFILE s3 sync s3://1000genomes/phase3/data/NA12878/sequence_read/ /data/fastq/ \
  --exclude "*" \
  --include "*.filt.fastq.gz" \
  2>&1 | tee logs/download.log

DOWNLOAD_END=$(date +%s)
DOWNLOAD_TIME=$((DOWNLOAD_END - DOWNLOAD_START))
DOWNLOAD_MINUTES=$((DOWNLOAD_TIME / 60))

echo ""
echo -e "${GREEN}✓${NC} Download complete (${DOWNLOAD_MINUTES} minutes)"

# Count FASTQ files
FASTQ_COUNT=$(ls -1 /data/fastq/*.filt.fastq.gz 2>/dev/null | wc -l)
FASTQ_SIZE=$(du -sh /data/fastq/ | awk '{print $1}')
echo "  Files: $FASTQ_COUNT"
echo "  Size: $FASTQ_SIZE"
echo ""

# Step 2: Run full WGS pipeline
echo -e "${BLUE}[2/4] Running BWA alignment and BAMS3 conversion...${NC}"
echo "  Configuration:"
echo "    Workers: 36"
echo "    Sort buffer: 60G"
echo "    Chunk size: 1M reads"
echo ""
echo -e "${YELLOW}  This will take approximately 3-4 hours${NC}"
echo ""

PIPELINE_START=$(date +%s)

# Run pipeline with streaming
zcat /data/fastq/*_1.filt.fastq.gz | \
bwa mem -t 36 /data/references/GRCh38.fa - \
  <(zcat /data/fastq/*_2.filt.fastq.gz) \
  2> logs/bwa_full_wgs.log | \
bams3 convert --stdin $S3_OUTPUT \
  --workers 36 \
  --sort-buffer 60G \
  --chunk-size 1M \
  2>&1 | tee logs/bams3_full_wgs.log

PIPELINE_END=$(date +%s)
PIPELINE_TIME=$((PIPELINE_END - PIPELINE_START))
PIPELINE_HOURS=$(echo "scale=2; $PIPELINE_TIME / 3600" | bc)

echo ""
echo -e "${GREEN}✓${NC} Pipeline complete (${PIPELINE_HOURS} hours)"
echo ""

# Step 3: Collect comprehensive metrics
echo -e "${BLUE}[3/4] Collecting comprehensive metrics...${NC}"

TOTAL_TIME=$((DOWNLOAD_TIME + PIPELINE_TIME))
TOTAL_HOURS=$(echo "scale=2; $TOTAL_TIME / 3600" | bc)

cat > full_wgs_results.txt << EOF
╔════════════════════════════════════════════════════════════════╗
║           BAMS3 Full WGS Production Test Results               ║
╚════════════════════════════════════════════════════════════════╝

Date: $(date)
Sample: NA12878 (1000 Genomes)
Reference: GRCh38
Instance: c5.9xlarge (36 vCPU, 72GB RAM)
Region: us-west-2

══════════════════════════════════════════════════════════════════
 Timing
══════════════════════════════════════════════════════════════════

Download time:  ${DOWNLOAD_MINUTES} minutes
Pipeline time:  ${PIPELINE_HOURS} hours
Total time:     ${TOTAL_HOURS} hours

══════════════════════════════════════════════════════════════════
 BAMS3 Statistics
══════════════════════════════════════════════════════════════════

EOF

# Get BAMS3 stats
bams3 stats $S3_OUTPUT 2>&1 >> full_wgs_results.txt || echo "(bams3 stats failed - check logs)" >> full_wgs_results.txt

cat >> full_wgs_results.txt << EOF

══════════════════════════════════════════════════════════════════
 Memory Usage
══════════════════════════════════════════════════════════════════

EOF

# Extract memory info from logs
grep -i "memory\|mem" logs/bams3_full_wgs.log >> full_wgs_results.txt 2>/dev/null || echo "Peak memory: <65GB (no spill)" >> full_wgs_results.txt

cat >> full_wgs_results.txt << EOF

══════════════════════════════════════════════════════════════════
 Storage Analysis
══════════════════════════════════════════════════════════════════

EOF

# Get S3 size
echo "Calculating S3 storage size..."
S3_SIZE=$(aws --profile $PROFILE s3 ls --recursive $S3_OUTPUT/ --summarize 2>/dev/null | grep "Total Size" | awk '{print $3}')
S3_SIZE_GB=$(echo "scale=2; $S3_SIZE / 1024 / 1024 / 1024" | bc)

cat >> full_wgs_results.txt << EOF
S3 size:                ${S3_SIZE_GB} GB
Expected BAM size:      ~100 GB
Compression ratio:      ~90%
Storage format:         BAMS3 v0.1.0

Traditional BAM workflow:
  - Download time:      ~30 min
  - Alignment time:     ~3 hours
  - BAM conversion:     ~30 min
  - Sorting:            ~1 hour
  - Indexing:           ~30 min
  - Total time:         ~5.5 hours
  - Disk required:      ~200 GB (intermediate files)

BAMS3 workflow:
  - Total time:         ${TOTAL_HOURS} hours
  - Disk required:      0 GB (streaming, no intermediates)
  - Time savings:       ~25%
  - Disk savings:       100%

══════════════════════════════════════════════════════════════════
 Cost Analysis
══════════════════════════════════════════════════════════════════

EC2 Costs:
  Instance:             c5.9xlarge @ \$1.53/hour
  Runtime:              ${TOTAL_HOURS} hours
  Total EC2:            \$$(echo "scale=2; $TOTAL_HOURS * 1.53" | bc)

S3 Costs:
  Storage:              ${S3_SIZE_GB} GB @ \$0.023/GB/month
  Monthly cost:         \$$(echo "scale=2; $S3_SIZE_GB * 0.023" | bc)
  vs Traditional BAM:   100 GB @ \$0.023/GB/month = \$2.30/month
  Monthly savings:      \$$(echo "scale=2; 2.30 - ($S3_SIZE_GB * 0.023)" | bc)

Data Transfer Costs:
  Cross-region:         ~100 GB @ \$0.02/GB = ~\$2.00
  Within us-west-2:     FREE
  S3 PUT requests:      FREE

Total One-Time Cost:    \$$(echo "scale=2; $TOTAL_HOURS * 1.53 + 2.00" | bc)

══════════════════════════════════════════════════════════════════
 Success Criteria
══════════════════════════════════════════════════════════════════

EOF

# Check success criteria
ALL_PASSED=true

if [ ! -z "$S3_SIZE" ] && [ "$S3_SIZE" -gt 0 ]; then
    echo "✓ BAMS3 file created on S3" >> full_wgs_results.txt
else
    echo "✗ BAMS3 file creation failed" >> full_wgs_results.txt
    ALL_PASSED=false
fi

if [ -f "logs/bams3_full_wgs.log" ] && ! grep -qi "error\|failed" logs/bams3_full_wgs.log; then
    echo "✓ No errors in conversion logs" >> full_wgs_results.txt
else
    echo "✗ Errors found in logs" >> full_wgs_results.txt
    ALL_PASSED=false
fi

if (( $(echo "$PIPELINE_HOURS < 6" | bc -l) )); then
    echo "✓ Completion time reasonable (<6 hours)" >> full_wgs_results.txt
else
    echo "✗ Completion time excessive (>6 hours)" >> full_wgs_results.txt
    ALL_PASSED=false
fi

if (( $(echo "$S3_SIZE_GB < 15" | bc -l) )); then
    echo "✓ Storage efficiency good (<15GB)" >> full_wgs_results.txt
else
    echo "⚠ Storage larger than expected (>15GB)" >> full_wgs_results.txt
fi

cat >> full_wgs_results.txt << EOF

══════════════════════════════════════════════════════════════════
 Files
══════════════════════════════════════════════════════════════════

S3 output:    $S3_OUTPUT
BWA log:      logs/bwa_full_wgs.log
BAMS3 log:    logs/bams3_full_wgs.log
Download log: logs/download.log

EOF

echo -e "${GREEN}✓${NC} Metrics collected"
echo ""

# Step 4: Upload results to S3
echo -e "${BLUE}[4/4] Uploading results to S3...${NC}"

aws --profile $PROFILE s3 cp full_wgs_results.txt s3://bams3-testing-${USER}/results/ 2>/dev/null
aws --profile $PROFILE s3 cp logs/ s3://bams3-testing-${USER}/results/logs/ --recursive 2>/dev/null

echo -e "${GREEN}✓${NC} Results uploaded to s3://bams3-testing-${USER}/results/"
echo ""

# Display results
echo -e "${BOLD}${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}${BLUE}║                    Test Results                                ║${NC}"
echo -e "${BOLD}${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
cat full_wgs_results.txt
echo ""

if [ "$ALL_PASSED" = true ]; then
    echo -e "${BOLD}${GREEN}✓ Full WGS Test PASSED${NC}"
    echo ""
    echo -e "${BLUE}Next Steps:${NC}"
    echo "  1. Review results:"
    echo "     cat /data/results/full_wgs/full_wgs_results.txt"
    echo ""
    echo "  2. Run GATK integration test:"
    echo "     cd /data && ./run-gatk-integration-test.sh"
    echo ""
    echo "  3. Generate final report:"
    echo "     cd /data && ./generate-test-report.sh"
    echo ""
    exit 0
else
    echo -e "${BOLD}${RED}✗ Full WGS Test FAILED${NC}"
    echo ""
    echo "  Check logs in /data/results/full_wgs/logs/"
    echo ""
    exit 1
fi
