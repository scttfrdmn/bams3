#!/bin/bash
set -e

# WGS Testing Environment Setup
# Based on docs/WGS_TESTING_PROTOCOL.md

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== BAMS3 WGS Testing Environment Setup ===${NC}\n"

# Load testing configuration
if [ -f .testing.config ]; then
    source .testing.config
    echo -e "${GREEN}✓${NC} Loaded testing configuration"
    echo "  Test data directory: $TEST_DATA_DIR"
else
    echo -e "${YELLOW}!${NC} No .testing.config found, using default"
    TEST_DATA_DIR="/Volumes/External HD/bams3-testing"
fi

# Create test directory
mkdir -p "$TEST_DATA_DIR"
cd "$TEST_DATA_DIR"
echo -e "${GREEN}✓${NC} Created test directory: $(pwd)\n"

# Check required software
echo -e "${BLUE}Checking required software...${NC}"

# BWA
if command -v bwa &> /dev/null; then
    BWA_VERSION=$(bwa 2>&1 | grep "Version" | awk '{print $2}')
    echo -e "${GREEN}✓${NC} BWA $BWA_VERSION installed"
else
    echo -e "${YELLOW}!${NC} BWA not found. Install with: brew install bwa"
    exit 1
fi

# samtools
if command -v samtools &> /dev/null; then
    SAMTOOLS_VERSION=$(samtools --version | head -1 | awk '{print $2}')
    echo -e "${GREEN}✓${NC} samtools $SAMTOOLS_VERSION installed"
else
    echo -e "${YELLOW}!${NC} samtools not found. Install with: brew install samtools"
    exit 1
fi

# Go
if command -v go &> /dev/null; then
    GO_VERSION=$(go version | awk '{print $3}')
    echo -e "${GREEN}✓${NC} Go $GO_VERSION installed"
else
    echo -e "${YELLOW}!${NC} Go not found. Install from: https://go.dev/dl/"
    exit 1
fi

# BAMS3
BAMS3_BIN="$OLDPWD/bams3-go/bams3"
if [ -f "$BAMS3_BIN" ]; then
    BAMS3_VERSION=$($BAMS3_BIN version | head -1 | awk '{print $2}')
    echo -e "${GREEN}✓${NC} BAMS3 $BAMS3_VERSION built"
    # Create symlink in test directory
    ln -sf "$BAMS3_BIN" bams3
else
    echo -e "${YELLOW}!${NC} BAMS3 not built. Building..."
    cd "$OLDPWD/bams3-go"
    go build -o bams3 ./cmd/bams3
    cd "$TEST_DATA_DIR"
    ln -sf "$BAMS3_BIN" bams3
    echo -e "${GREEN}✓${NC} BAMS3 built successfully"
fi

# GATK (optional but recommended)
if command -v gatk &> /dev/null; then
    GATK_VERSION=$(gatk --version 2>&1 | grep "GATK" | awk '{print $6}')
    echo -e "${GREEN}✓${NC} GATK $GATK_VERSION installed"
else
    echo -e "${YELLOW}!${NC} GATK not found (optional)"
    echo "  Install from: https://github.com/broadinstitute/gatk/releases"
    echo "  Or with conda: conda install -c bioconda gatk4"
fi

# AWS CLI
if command -v aws &> /dev/null; then
    AWS_VERSION=$(aws --version 2>&1 | awk '{print $1}' | cut -d/ -f2)
    echo -e "${GREEN}✓${NC} AWS CLI $AWS_VERSION installed"
else
    echo -e "${YELLOW}!${NC} AWS CLI not found (required for 1000 Genomes data)"
    echo "  Install with: brew install awscli"
    exit 1
fi

echo ""

# Check AWS credentials
echo -e "${BLUE}Checking AWS configuration...${NC}"
if aws sts get-caller-identity &> /dev/null; then
    AWS_ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
    AWS_REGION=$(aws configure get region || echo "us-east-1")
    echo -e "${GREEN}✓${NC} AWS credentials configured"
    echo "  Account: $AWS_ACCOUNT"
    echo "  Region: $AWS_REGION"
else
    echo -e "${YELLOW}!${NC} AWS credentials not configured"
    echo "  Run: aws configure"
    echo "  You'll need access to 1000 Genomes public dataset"
fi

echo ""

# Download test reference (small for validation)
echo -e "${BLUE}Setting up test reference genome...${NC}"
if [ ! -f chr22.fa ]; then
    echo "Downloading chromosome 22 (smallest, for quick testing)..."
    # We'll download chr22 from UCSC
    curl -L -o chr22.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
    gunzip chr22.fa.gz
    echo -e "${GREEN}✓${NC} Downloaded chr22.fa (49 MB)"

    echo "Indexing reference..."
    samtools faidx chr22.fa
    bwa index chr22.fa
    echo -e "${GREEN}✓${NC} Indexed chr22.fa"
else
    echo -e "${GREEN}✓${NC} chr22.fa already present"
fi

echo ""

# Check disk space
echo -e "${BLUE}Checking disk space...${NC}"
AVAIL_GB=$(df -g "$TEST_DATA_DIR" | tail -1 | awk '{print $4}')
echo "  Available space: ${AVAIL_GB} GB"

if [ "$AVAIL_GB" -lt 50 ]; then
    echo -e "${YELLOW}!${NC} Warning: Less than 50GB available"
    echo "  Recommended for WGS testing: 200GB"
    echo "  You have enough for small/medium tests"
else
    echo -e "${GREEN}✓${NC} Sufficient disk space for testing"
fi

echo ""

# Create directory structure
mkdir -p {references,fastq,bams3,results,logs}
echo -e "${GREEN}✓${NC} Created directory structure:"
echo "  references/ - Reference genomes"
echo "  fastq/      - FASTQ input files"
echo "  bams3/      - BAMS3 output files"
echo "  results/    - Test results and reports"
echo "  logs/       - Execution logs"

echo ""
echo -e "${GREEN}=== Setup Complete ===${NC}\n"
echo "Test directory: $TEST_DATA_DIR"
echo ""
echo "Next steps:"
echo "  1. For small test (chr22, quick validation):"
echo "     cd $TEST_DATA_DIR"
echo "     ./quick-test.sh"
echo ""
echo "  2. For full WGS test (1000 Genomes NA12878):"
echo "     Follow: docs/WGS_TESTING_PROTOCOL.md"
echo ""
echo "  3. Check available test data:"
echo "     aws s3 ls s3://1000genomes/phase3/data/NA12878/sequence_read/ --no-sign-request"
