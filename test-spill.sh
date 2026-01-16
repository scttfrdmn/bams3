#!/bin/bash
# Test script to verify disk spill functionality

set -e

cd bams3-go

echo "Building bams3..."
go build -o bams3 ./cmd/bams3

echo ""
echo "Generating test SAM data..."

# Generate SAM header
cat > /tmp/test-spill.sam << 'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248956422
@SQ	SN:chr2	LN:242193529
@SQ	SN:chr3	LN:198295559
@PG	ID:test	PN:test	VN:1.0
EOF

# Generate many reads across chromosomes to create memory pressure
# This will create reads that can't be incrementally flushed due to frontier
echo "Generating 100K reads per chromosome (300K total)..."

# Generate reads for chr1
for i in $(seq 1 100000); do
    pos=$((1000000 + i * 100))
    echo "read_${i}_chr1	0	chr1	${pos}	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	*"
done >> /tmp/test-spill.sam

# Generate reads for chr2
for i in $(seq 1 100000); do
    pos=$((1000000 + i * 100))
    echo "read_${i}_chr2	0	chr2	${pos}	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	*"
done >> /tmp/test-spill.sam

# Generate reads for chr3
for i in $(seq 1 100000); do
    pos=$((1000000 + i * 100))
    echo "read_${i}_chr3	0	chr3	${pos}	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	*"
done >> /tmp/test-spill.sam

echo "Generated SAM file: $(wc -l /tmp/test-spill.sam | awk '{print $1}') lines"
echo ""

echo "Testing conversion with small sort buffer (to force spilling)..."
# Use a very small sort buffer (50MB) to force spilling
cat /tmp/test-spill.sam | ./bams3 convert --stdin /tmp/test-output.bams3 \
    --sort-buffer 50M \
    --workers 2

echo ""
echo "Conversion complete!"
echo ""

# Verify output exists
if [ -d "/tmp/test-output.bams3" ]; then
    echo "Output directory created successfully"
    du -sh /tmp/test-output.bams3
    echo ""
    echo "Contents:"
    find /tmp/test-output.bams3 -type f | head -20
else
    echo "ERROR: Output directory not created"
    exit 1
fi

# Cleanup
echo ""
echo "Cleaning up test files..."
rm -f /tmp/test-spill.sam
# rm -rf /tmp/test-output.bams3

echo "Test complete!"
