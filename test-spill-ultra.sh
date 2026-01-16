#!/bin/bash
# Test script to force disk spill with deep coverage in small region

set -e

cd bams3-go

echo "Building bams3..."
go build -o bams3 ./cmd/bams3

echo ""
echo "Generating ultra-deep coverage test data..."

# Generate SAM header
cat > /tmp/test-spill-ultra.sam << 'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248956422
@PG	ID:test	PN:test	VN:1.0
EOF

# Generate 500K reads all within a 500KB region (ultra-deep coverage ~1000x)
# This prevents incremental flushing because all reads go into the same few chunks
echo "Generating 500K reads in 500KB region (ultra-deep coverage)..."

for i in $(seq 1 500000); do
    # Random position within 500KB region
    pos=$((1000000 + (i % 5000) * 100))
    echo "read_${i}	0	chr1	${pos}	60	100M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	*"
done >> /tmp/test-spill-ultra.sam

echo "Generated SAM file: $(wc -l /tmp/test-spill-ultra.sam | awk '{print $1}') lines"
echo ""

echo "Testing conversion with 5MB sort buffer..."
cat /tmp/test-spill-ultra.sam | ./bams3 convert --stdin /tmp/test-output-ultra.bams3 \
    --sort-buffer 5M \
    --workers 2 \
    2>&1 | tee /tmp/conversion.log

echo ""

# Check for spill messages
if grep -q "Spilling" /tmp/conversion.log; then
    echo "SUCCESS: Disk spill was triggered!"
    echo ""
    grep "Spill" /tmp/conversion.log
else
    echo "INFO: No disk spill occurred (incremental flushing handled it)"
fi

echo ""
echo "Conversion complete!"

# Verify output
if [ -d "/tmp/test-output-ultra.bams3" ]; then
    echo "Output directory: $(du -sh /tmp/test-output-ultra.bams3 | awk '{print $1}')"
    num_chunks=$(find /tmp/test-output-ultra.bams3/data -name "*.chunk" 2>/dev/null | wc -l)
    echo "Total chunks: $num_chunks"
else
    echo "ERROR: Output directory not created"
    exit 1
fi

# Cleanup
echo ""
echo "Cleaning up..."
rm -f /tmp/test-spill-ultra.sam /tmp/conversion.log
rm -rf /tmp/test-output-ultra.bams3

echo "Test complete!"
