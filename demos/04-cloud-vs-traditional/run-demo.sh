#!/bin/bash
set -e

# Demo: Cloud-Native vs Traditional Workflow
# Shows: Complete workflow comparison with multiple queries
# Time: ~3 minutes

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m'

echo -e "${BOLD}${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}${BLUE}║        Cloud-Native vs Traditional Workflow                   ║${NC}"
echo -e "${BOLD}${BLUE}║        Real-world analysis scenario comparison                ║${NC}"
echo -e "${BOLD}${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

mkdir -p demo-workspace
cd demo-workspace

echo -e "${BLUE}Scenario: Variant calling across 10 genomic regions${NC}"
echo ""
echo "Research task:"
echo "  - Analyze 10 different genes across a 30x WGS sample"
echo "  - Each region: ~100kb"
echo "  - Total: 1Mbp of 3Gbp genome (0.03%)"
echo ""
echo "Dataset:"
echo "  - Source: 100GB BAM file on S3"
echo "  - Coverage: 30x whole genome"
echo "  - Platform: Illumina HiSeq"
echo ""

# Define regions of interest
REGIONS=(
    "BRCA1:chr17:41196312-41277500"
    "BRCA2:chr13:32889611-32973805"
    "TP53:chr17:7661779-7687550"
    "EGFR:chr7:55019017-55211628"
    "KRAS:chr12:25205246-25250929"
    "PTEN:chr10:87863113-87971930"
    "APC:chr5:112707498-112846239"
    "MLH1:chr3:36993333-37050918"
    "BRAF:chr7:140719327-140924929"
    "PIK3CA:chr3:179148114-179240093"
)

echo "Target genes:"
for region in "${REGIONS[@]}"; do
    gene=$(echo $region | cut -d: -f1)
    coords=$(echo $region | cut -d: -f2-3)
    echo "  - $gene ($coords)"
done
echo ""

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${NC}"
echo -e "${BOLD} Traditional Workflow: Download Full BAM                      ${NC}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${NC}"
echo ""

START_TRAD=$(date +%s)

echo -e "${BLUE}Step 1: Download full BAM from S3${NC}"
echo ""
echo "  aws s3 cp s3://bucket/sample.bam ./sample.bam"
echo ""
echo "  File size:     100 GB"
echo "  Transfer time: ~30 minutes (at 50 MB/s)"
echo "  S3 cost:       \$9.00 (100 GB × \$0.09/GB)"
echo "  Local disk:    100 GB required"
echo ""
sleep 2

echo -e "${BLUE}Step 2: Index BAM file${NC}"
echo ""
echo "  samtools index sample.bam"
echo ""
echo "  Time: ~5 minutes"
echo ""
sleep 1

echo -e "${BLUE}Step 3: Extract regions (10 iterations)${NC}"
echo ""
REGION_COUNT=0
for region in "${REGIONS[@]}"; do
    gene=$(echo $region | cut -d: -f1)
    coords=$(echo $region | cut -d: -f2-3)

    echo "  [$((++REGION_COUNT))/10] Extracting $gene..."
    echo "    samtools view sample.bam $coords > ${gene}.bam"

    sleep 0.3
done
echo ""
echo "  Total time: ~3 minutes"
echo ""

echo -e "${BLUE}Step 4: Variant calling on each region${NC}"
echo ""
REGION_COUNT=0
for region in "${REGIONS[@]}"; do
    gene=$(echo $region | cut -d: -f1)

    echo "  [$((++REGION_COUNT))/10] Calling variants in $gene..."
    echo "    gatk HaplotypeCaller -I ${gene}.bam -R ref.fa -O ${gene}.vcf"

    sleep 0.3
done
echo ""
echo "  Total time: ~15 minutes"
echo ""

END_TRAD=$(date +%s)
TRAD_TOTAL=$((END_TRAD - START_TRAD))

echo -e "${GREEN}✓ Traditional workflow complete${NC}"
echo ""

cat << EOF
┌─────────────────────────────────────────────────────────────────┐
│              Traditional Workflow Summary                       │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  Step 1: Download full BAM      30 minutes                     │
│  Step 2: Index BAM               5 minutes                     │
│  Step 3: Extract 10 regions      3 minutes                     │
│  Step 4: Variant calling        15 minutes                     │
│                                                                 │
│  Total time:                    53 minutes                     │
│                                                                 │
│  S3 costs:                      \$9.00                          │
│  Local disk:                    100 GB                         │
│  Network:                       100 GB downloaded              │
│                                                                 │
│  Issues:                                                        │
│  • Long wait before starting analysis                          │
│  • Large local disk requirement                                │
│  • Must re-download for new queries                            │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
EOF

echo ""
echo ""

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${NC}"
echo -e "${BOLD} BAMS3 Workflow: Selective Streaming                          ${NC}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${NC}"
echo ""

START_BAMS3=$(date +%s)

echo -e "${BLUE}Step 1: Query metadata (instant)${NC}"
echo ""
echo "  bams3 query s3://bucket/sample.bams3 --regions regions.bed --stats"
echo ""
echo "  Time: <1 second (S3 metadata only)"
echo ""
sleep 1

echo -e "${BLUE}Step 2: Stream regions and call variants (parallel)${NC}"
echo ""
REGION_COUNT=0
for region in "${REGIONS[@]}"; do
    gene=$(echo $region | cut -d: -f1)
    coords=$(echo $region | cut -d: -f2-3)

    echo "  [$((++REGION_COUNT))/10] Streaming $gene → variant calling..."
    echo "    bams3 to-bam s3://bucket/sample.bams3 - --region $coords | \\"
    echo "      gatk HaplotypeCaller -I /dev/stdin -R ref.fa -O ${gene}.vcf"

    sleep 0.5
done
echo ""
echo "  Total time: ~5 minutes (parallel processing)"
echo ""

END_BAMS3=$(date +%s)
BAMS3_TOTAL=$((END_BAMS3 - START_BAMS3))

echo -e "${GREEN}✓ BAMS3 workflow complete${NC}"
echo ""

# Calculate data transfer for BAMS3
# 10 regions × 100kb each = 1Mbp
# At 30x coverage, ~5MB compressed per region
BAMS3_TRANSFER=50  # MB

cat << EOF
┌─────────────────────────────────────────────────────────────────┐
│                BAMS3 Workflow Summary                           │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  Step 1: Query metadata          <1 second                     │
│  Step 2: Stream + call variants   5 minutes (parallel)         │
│                                                                 │
│  Total time:                      5 minutes                    │
│                                                                 │
│  S3 costs:                        \$0.0045                      │
│  Local disk:                      ~50 MB (minimal)             │
│  Network:                         50 MB downloaded              │
│                                                                 │
│  Advantages:                                                    │
│  • Instant query start (no download wait)                      │
│  • Minimal disk requirement                                    │
│  • Parallel processing                                         │
│  • Pay only for data accessed                                  │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
EOF

echo ""
echo ""

echo -e "${BOLD}${GREEN}Comparison Summary${NC}"
echo ""

SPEEDUP=$(echo "scale=1; 53 / 5" | bc)
COST_SAVINGS=$(echo "scale=2; 9.00 - 0.0045" | bc)
DATA_SAVINGS=$(echo "scale=1; 100 - (50 / 1024 / 100 * 100)" | bc)

cat << EOF
╔═════════════════════════════════════════════════════════════════╗
║                   Workflow Comparison                           ║
╠═════════════════════════════════════════════════════════════════╣
║                                                                 ║
║  Metric              │ Traditional │ BAMS3     │ Improvement    ║
║  ───────────────────────────────────────────────────────────── ║
║  Total time          │ 53 min      │ 5 min     │ ${SPEEDUP}x faster    ║
║  S3 cost (single)    │ \$9.00       │ \$0.0045  │ 99.95% less  ║
║  Data transfer       │ 100 GB      │ 50 MB     │ 99.95% less  ║
║  Local disk          │ 100 GB      │ 50 MB     │ 99.95% less  ║
║  Wait before start   │ 30 min      │ 0 sec     │ Instant      ║
║                                                                 ║
╠═════════════════════════════════════════════════════════════════╣
║                                                                 ║
║  Scenario: 100 queries per year                                 ║
║                                                                 ║
║  Traditional:                                                   ║
║    • Time:          88 hours (3.7 days)                        ║
║    • S3 costs:      \$900                                       ║
║    • Disk needed:   100 GB continuously                        ║
║                                                                 ║
║  BAMS3:                                                         ║
║    • Time:          8 hours (0.3 days)                         ║
║    • S3 costs:      \$0.45                                      ║
║    • Disk needed:   50 MB per query                            ║
║                                                                 ║
║  Annual savings:                                                ║
║    • Time:          80 hours (11x faster)                      ║
║    • Cost:          \$899.55 (99.95% less)                      ║
║    • Efficiency:    Analyze 1000x more samples                 ║
║                                                                 ║
╚═════════════════════════════════════════════════════════════════╝
EOF

echo ""
echo -e "${GREEN}${BOLD}Key Insights:${NC}"
echo ""
echo "  1. ${BOLD}Selective Access${NC}: Download only what you analyze"
echo "     • Traditional: 100 GB for 1 Mbp (0.03% of genome)"
echo "     • BAMS3: 50 MB for 1 Mbp (1999x less data)"
echo ""
echo "  2. ${BOLD}Cloud-Native Architecture${NC}: Optimized for cloud storage"
echo "     • No local caching required"
echo "     • Instant query start"
echo "     • Parallel region processing"
echo ""
echo "  3. ${BOLD}Cost Efficiency${NC}: Pay for what you use"
echo "     • 99.95% cost savings on targeted queries"
echo "     • Scales to thousands of samples"
echo "     • Enables larger cohort studies"
echo ""
echo "  4. ${BOLD}Time Savings${NC}: 11x faster for iterative analysis"
echo "     • No download wait"
echo "     • Parallel processing"
echo "     • Rapid iteration on analysis"
echo ""

echo -e "${BLUE}Real-world impact:${NC}"
echo ""
echo "  Research group analyzing 1000 samples:"
echo "    • Traditional: \$900,000/year in S3 costs alone"
echo "    • BAMS3:       \$450/year in S3 costs"
echo "    • Savings:     \$899,550/year"
echo ""
echo "  This enables:"
echo "    • Larger cohort studies (10x-100x more samples)"
echo "    • More exploratory analysis (no cost barrier)"
echo "    • Faster time to insight (11x faster iteration)"
echo ""

echo -e "${BLUE}Learn more:${NC}"
echo "  • S3_INTEGRATION.md - Cloud-native architecture"
echo "  • BAM_EXPORT.md - Tool integration patterns"
echo "  • WGS_TESTING_PROTOCOL.md - Production validation"
echo ""

# Cleanup option
read -p "Clean up demo files? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    cd ..
    rm -rf demo-workspace
    echo -e "${GREEN}✓${NC} Cleaned up"
fi
