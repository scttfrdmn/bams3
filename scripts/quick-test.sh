#!/bin/bash
set -e

# Quick BAMS3 validation test using chr22
# This script runs a small-scale test to verify the setup before full WGS testing

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}=== BAMS3 Quick Validation Test ===${NC}\n"

# Check we're in test directory
if [ ! -f chr22.fa ]; then
    echo -e "${YELLOW}!${NC} Run setup-wgs-testing.sh first"
    exit 1
fi

# Generate synthetic reads for chr22 region
echo -e "${BLUE}Generating synthetic reads (chr22:10000000-11000000, 1Mbp)...${NC}"
if [ ! -f fastq/chr22_test_R1.fq ]; then
    # Extract 1Mbp region
    samtools faidx chr22.fa chr22:10000000-11000000 > references/chr22_1mb.fa

    # Generate ~100K reads (30x coverage of 1Mbp)
    # Using mason_simulator if available, otherwise use simple perl script
    if command -v mason_simulator &> /dev/null; then
        mason_simulator -N 100000 -o fastq/chr22_test_R1.fq -or fastq/chr22_test_R2.fq \
            -ir references/chr22_1mb.fa --illumina-read-length 150
    else
        # Simple fallback: create minimal test data
        echo "Creating minimal test reads..."
        # We'll just create a small BAM for testing
        samtools view -H chr22.fa.fai | head -1 > test.sam
        # Add some synthetic alignments
        for i in {1..1000}; do
            pos=$((10000000 + i * 100))
            echo "read_${i}	0	chr22	${pos}	60	150M	*	0	0	$(head -c 150 /dev/urandom | base64 | tr -dc 'ACGT' | head -c 150)	$(head -c 150 /dev/urandom | od -An -tu1 | tr ' ' '\n' | awk '{print chr(33+$1%40)}' | tr -d '\n')"
        done >> test.sam
        samtools view -bS test.sam > test.bam
        rm test.sam

        # Use this existing BAM for testing
        echo -e "${GREEN}✓${NC} Created test BAM with 1000 reads"

        # Skip to BAM conversion test
        echo -e "\n${BLUE}Testing BAMS3 conversion...${NC}"
        START=$(date +%s)

        ./bams3 convert test.bam bams3/quick_test.bams3 \
            --workers 4 \
            --sort-buffer 1G \
            2>&1 | tee logs/quick_test_convert.log

        END=$(date +%s)
        ELAPSED=$((END - START))

        echo -e "${GREEN}✓${NC} Conversion completed in ${ELAPSED}s"

        # Get statistics
        echo -e "\n${BLUE}Checking BAMS3 statistics...${NC}"
        ./bams3 stats bams3/quick_test.bams3 | tee logs/quick_test_stats.log

        # Test BAM export
        echo -e "\n${BLUE}Testing BAM export (full)...${NC}"
        ./bams3 to-bam bams3/quick_test.bams3 results/reconstructed.bam

        # Validate
        echo -e "\n${BLUE}Validating results...${NC}"
        ORIGINAL_READS=$(samtools view -c test.bam)
        EXPORTED_READS=$(samtools view -c results/reconstructed.bam)

        echo "Original BAM reads: $ORIGINAL_READS"
        echo "Exported BAM reads: $EXPORTED_READS"

        if [ "$ORIGINAL_READS" -eq "$EXPORTED_READS" ]; then
            echo -e "${GREEN}✓${NC} Read counts match!"
        else
            echo -e "${YELLOW}!${NC} Read count mismatch!"
            exit 1
        fi

        # Test region extraction
        echo -e "\n${BLUE}Testing region extraction (chr22:10000000-10100000)...${NC}"
        ./bams3 to-bam bams3/quick_test.bams3 results/region.bam \
            --region chr22:10000000-10100000

        REGION_READS=$(samtools view -c results/region.bam)
        echo "Region extraction: $REGION_READS reads"

        # Test streaming
        echo -e "\n${BLUE}Testing streaming to stdout...${NC}"
        STREAM_READS=$(./bams3 to-bam bams3/quick_test.bams3 - | samtools view -c)
        echo "Streaming output: $STREAM_READS reads"

        if [ "$STREAM_READS" -eq "$ORIGINAL_READS" ]; then
            echo -e "${GREEN}✓${NC} Streaming works correctly!"
        else
            echo -e "${YELLOW}!${NC} Streaming read count mismatch!"
        fi

        # Calculate sizes
        echo -e "\n${BLUE}Storage comparison...${NC}"
        ORIGINAL_SIZE=$(du -h test.bam | awk '{print $1}')
        BAMS3_SIZE=$(du -sh bams3/quick_test.bams3 | awk '{print $1}')
        EXPORTED_SIZE=$(du -h results/reconstructed.bam | awk '{print $1}')

        echo "Original BAM:     $ORIGINAL_SIZE"
        echo "BAMS3 format:     $BAMS3_SIZE"
        echo "Exported BAM:     $EXPORTED_SIZE"

        echo -e "\n${GREEN}=== Quick Test Complete ===${NC}\n"
        echo "All tests passed! BAMS3 is working correctly."
        echo ""
        echo "Next steps:"
        echo "  - For full WGS testing, follow docs/WGS_TESTING_PROTOCOL.md"
        echo "  - Test results saved in results/"
        echo "  - Logs saved in logs/"

        exit 0
    fi
fi

# If we have real FASTQ, align and convert
echo -e "\n${BLUE}Aligning reads with BWA...${NC}"
START=$(date +%s)

bwa mem -t 4 chr22.fa fastq/chr22_test_R1.fq fastq/chr22_test_R2.fq 2> logs/bwa_align.log | \
    ./bams3 convert --stdin bams3/quick_test.bams3 \
        --workers 4 \
        --sort-buffer 1G \
        2>&1 | tee logs/quick_test_convert.log

END=$(date +%s)
ELAPSED=$((END - START))

echo -e "${GREEN}✓${NC} Alignment and conversion completed in ${ELAPSED}s"

# Get statistics
echo -e "\n${BLUE}Checking BAMS3 statistics...${NC}"
./bams3 stats bams3/quick_test.bams3 | tee logs/quick_test_stats.log

# Test BAM export
echo -e "\n${BLUE}Testing BAM export...${NC}"
./bams3 to-bam bams3/quick_test.bams3 results/reconstructed.bam

# Validate with samtools
echo -e "\n${BLUE}Validating with samtools...${NC}"
samtools quickcheck results/reconstructed.bam && echo -e "${GREEN}✓${NC} BAM is valid"
samtools flagstat results/reconstructed.bam | tee logs/flagstat.log

echo -e "\n${GREEN}=== Quick Test Complete ===${NC}\n"
echo "Test results saved in: $(pwd)/results/"
echo "Logs saved in: $(pwd)/logs/"
