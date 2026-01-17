#!/bin/bash
set -e

# BAMS3 WGS Testing - EC2 Instance Setup
# Installs: Go, BWA, samtools, BAMS3, GRCh38 reference
# Time: ~90 minutes (reference indexing takes ~60 minutes)

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m'

echo -e "${BOLD}${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}${BLUE}║    BAMS3 WGS Testing - EC2 Setup                              ║${NC}"
echo -e "${BOLD}${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

START_TIME=$(date +%s)

# Step 1: Update system
echo -e "${BLUE}[1/8] Updating system packages...${NC}"
yum update -y
yum install -y gcc make git wget bzip2 zlib-devel ncurses-devel bzip2-devel xz-devel
echo -e "${GREEN}✓${NC} System updated"
echo ""

# Step 2: Install Go
echo -e "${BLUE}[2/8] Installing Go 1.21.5...${NC}"
cd /tmp
wget -q https://go.dev/dl/go1.21.5.linux-amd64.tar.gz
tar -C /usr/local -xzf go1.21.5.linux-amd64.tar.gz
export PATH=$PATH:/usr/local/go/bin
echo 'export PATH=$PATH:/usr/local/go/bin' >> /etc/profile.d/go.sh
echo -e "${GREEN}✓${NC} Go installed: $(go version)"
echo ""

# Step 3: Install BWA
echo -e "${BLUE}[3/8] Installing BWA 0.7.17...${NC}"
cd /tmp
wget -q https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make > /dev/null 2>&1
cp bwa /usr/local/bin/
cd /tmp
rm -rf bwa-0.7.17*
echo -e "${GREEN}✓${NC} BWA installed: $(bwa 2>&1 | grep Version | head -1)"
echo ""

# Step 4: Install samtools
echo -e "${BLUE}[4/8] Installing samtools 1.18...${NC}"
cd /tmp
wget -q https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
tar -xjf samtools-1.18.tar.bz2
cd samtools-1.18
./configure --prefix=/usr/local > /dev/null 2>&1
make > /dev/null 2>&1
make install > /dev/null 2>&1
cd /tmp
rm -rf samtools-1.18*
echo -e "${GREEN}✓${NC} samtools installed: $(samtools --version | head -1)"
echo ""

# Step 5: Install BAMS3
echo -e "${BLUE}[5/8] Installing BAMS3...${NC}"
cd /opt
git clone https://github.com/scttfrdmn/aws-direct-s3.git bams3 2>&1 | grep -v "Cloning" || true
cd /opt/bams3/bams3-go
go build -o bams3 ./cmd/bams3 2>&1 | grep -v "go: downloading" || true
cp bams3 /usr/local/bin/
echo -e "${GREEN}✓${NC} BAMS3 installed: $(bams3 version 2>&1 || echo 'v0.1.0')"
echo ""

# Step 6: Create directories
echo -e "${BLUE}[6/8] Creating data directories...${NC}"
mkdir -p /data/references
mkdir -p /data/results/logs
mkdir -p /data/fastq
echo -e "${GREEN}✓${NC} Directories created"
echo ""

# Step 7: Download GRCh38 reference
echo -e "${BLUE}[7/8] Downloading GRCh38 reference (~3GB)...${NC}"
cd /data/references
wget -q --show-progress \
  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
echo ""
echo "Decompressing..."
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38.fa
echo -e "${GREEN}✓${NC} Reference downloaded: $(ls -lh GRCh38.fa | awk '{print $5}')"
echo ""

# Step 8: Index reference (takes ~60 minutes)
echo -e "${BLUE}[8/8] Indexing reference (this will take ~60 minutes)...${NC}"
echo ""

echo -e "${YELLOW}Starting BWA indexing...${NC}"
INDEX_START=$(date +%s)
bwa index GRCh38.fa
INDEX_END=$(date +%s)
INDEX_TIME=$((INDEX_END - INDEX_START))
echo -e "${GREEN}✓${NC} BWA index complete (${INDEX_TIME}s)"
echo ""

echo -e "${YELLOW}Starting samtools indexing...${NC}"
samtools faidx GRCh38.fa
echo -e "${GREEN}✓${NC} samtools index complete"
echo ""

END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))
MINUTES=$((TOTAL_TIME / 60))

echo -e "${BOLD}${GREEN}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}${GREEN}║                Setup Complete!                                 ║${NC}"
echo -e "${BOLD}${GREEN}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${BLUE}Total setup time: ${MINUTES} minutes${NC}"
echo ""
echo -e "${BLUE}Installed Components:${NC}"
echo "  • Go:       $(go version | awk '{print $3}')"
echo "  • BWA:      $(bwa 2>&1 | grep Version | head -1 | awk '{print $2}')"
echo "  • samtools: $(samtools --version | head -1 | awk '{print $2}')"
echo "  • BAMS3:    v0.1.0"
echo ""
echo -e "${BLUE}Data Directories:${NC}"
echo "  • References: /data/references/"
echo "  • Results:    /data/results/"
echo "  • FASTQ:      /data/fastq/"
echo ""
echo -e "${BLUE}Reference:${NC}"
echo "  • GRCh38.fa:      $(ls -lh /data/references/GRCh38.fa | awk '{print $5}')"
echo "  • BWA index:      $(ls /data/references/GRCh38.fa.bwt 2>&1 >/dev/null && echo '✓' || echo '✗')"
echo "  • samtools index: $(ls /data/references/GRCh38.fa.fai 2>&1 >/dev/null && echo '✓' || echo '✗')"
echo ""
echo -e "${BOLD}${YELLOW}Ready to run tests!${NC}"
echo ""
echo -e "${BLUE}Next Steps:${NC}"
echo "  1. Copy test scripts to instance:"
echo "     scp scripts/run-chr22-test.sh ec2-user@<instance-ip>:/data/"
echo "     scp scripts/run-full-wgs-test.sh ec2-user@<instance-ip>:/data/"
echo "     scp scripts/run-gatk-integration-test.sh ec2-user@<instance-ip>:/data/"
echo ""
echo "  2. Run chr22 validation test:"
echo "     cd /data && chmod +x run-chr22-test.sh && ./run-chr22-test.sh"
echo ""
echo "  3. Run full WGS test (if chr22 passes):"
echo "     cd /data && chmod +x run-full-wgs-test.sh && ./run-full-wgs-test.sh"
echo ""
echo -e "${YELLOW}Estimated test times:${NC}"
echo "  • Chr22 test:     30 minutes"
echo "  • Full WGS test:  4 hours"
echo "  • GATK test:      30 minutes"
echo ""
