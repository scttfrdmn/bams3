#!/bin/bash
# Test script to force disk spill by generating reads that prevent incremental flushing

set -e

cd bams3-go

echo "Building bams3..."
go build -o bams3 ./cmd/bams3

echo ""
echo "Generating challenging test SAM data..."

# Generate SAM header
cat > /tmp/test-spill-hard.sam << 'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248956422
@PG	ID:test	PN:test	VN:1.0
EOF

# Generate reads in a pattern that prevents incremental flushing:
# - All reads are on chr1
# - Reads alternate between positions at the start and end of the range
# - This keeps the frontier advancing slowly while accumulating many chunks
# that can't be flushed yet

echo "Generating 200K reads with interleaved positions (forces memory pressure)..."

for i in $(seq 1 200000); do
    # Alternate between early and late positions
    if [ $((i % 2)) -eq 0 ]; then
        # Early position
        pos=$((1000000 + (i / 2) * 50))
    else
        # Late position
        pos=$((50000000 + (i / 2) * 50))
    fi
    echo "read_${i}	0	chr1	${pos}	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	*"
done >> /tmp/test-spill-hard.sam

echo "Generated SAM file: $(wc -l /tmp/test-spill-hard.sam | awk '{print $1}') lines"
echo ""

echo "Testing conversion with tiny sort buffer (10MB) to force spilling..."
# Use a very small sort buffer (10MB) to force spilling
cat /tmp/test-spill-hard.sam | ./bams3 convert --stdin /tmp/test-output-hard.bams3 \
    --sort-buffer 10M \
    --workers 2

echo ""
echo "Conversion complete!"
echo ""

# Verify output exists
if [ -d "/tmp/test-output-hard.bams3" ]; then
    echo "Output directory created successfully"
    du -sh /tmp/test-output-hard.bams3
    echo ""

    # Count chunks
    num_chunks=$(find /tmp/test-output-hard.bams3/data -name "*.chunk" 2>/dev/null | wc -l)
    echo "Total chunks created: $num_chunks"

    # Check metadata
    if [ -f "/tmp/test-output-hard.bams3/_metadata.json" ]; then
        echo ""
        echo "Metadata preview:"
        head -20 /tmp/test-output-hard.bams3/_metadata.json
    fi
else
    echo "ERROR: Output directory not created"
    exit 1
fi

# Cleanup
echo ""
echo "Cleaning up test files..."
rm -f /tmp/test-spill-hard.sam
rm -rf /tmp/test-output-hard.bams3

echo "Test complete!"
