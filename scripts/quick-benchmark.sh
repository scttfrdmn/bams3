#!/bin/bash
#
# Quick benchmark script to compare copy-then-process vs FUSE mounting
#
# Usage:
#   ./quick-benchmark.sh s3://bucket/file.bam samtools flagstat
#   ./quick-benchmark.sh s3://bucket/file.vcf.gz bcftools view

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
S3_FILE="${1:-}"
TOOL="${2:-samtools}"
OPERATION="${3:-flagstat}"
CACHE_DIR="${BENCHMARK_CACHE:-/tmp/benchmark-cache}"
MOUNT_POINT="$CACHE_DIR/mount"

usage() {
    echo "Usage: $0 <s3-file> [tool] [operation]"
    echo ""
    echo "Examples:"
    echo "  $0 s3://bucket/sample.bam samtools flagstat"
    echo "  $0 s3://bucket/variants.vcf.gz bcftools view"
    echo "  $0 s3://bucket/reads.fastq.gz head -n 1000"
    echo ""
    echo "Environment variables:"
    echo "  BENCHMARK_CACHE - Cache directory (default: /tmp/benchmark-cache)"
    exit 1
}

if [ -z "$S3_FILE" ]; then
    usage
fi

if [[ ! "$S3_FILE" =~ ^s3:// ]]; then
    echo -e "${RED}Error: First argument must be an S3 URI (s3://bucket/key)${NC}"
    exit 1
fi

# Parse S3 URI
S3_PATH="${S3_FILE#s3://}"
BUCKET="${S3_PATH%%/*}"
KEY="${S3_PATH#*/}"
FILENAME="$(basename "$KEY")"

echo -e "${BLUE}======================================${NC}"
echo -e "${BLUE}S3 Access Method Benchmark${NC}"
echo -e "${BLUE}======================================${NC}"
echo ""
echo "File:      $S3_FILE"
echo "Tool:      $TOOL $OPERATION"
echo "Cache dir: $CACHE_DIR"
echo ""

# Get file size
echo -e "${YELLOW}Getting file size...${NC}"
FILE_SIZE_BYTES=$(aws s3api head-object --bucket "$BUCKET" --key "$KEY" --query 'ContentLength' --output text 2>/dev/null || echo "0")
FILE_SIZE_GB=$(echo "scale=2; $FILE_SIZE_BYTES / 1024 / 1024 / 1024" | bc)
echo "Size: $FILE_SIZE_GB GB"
echo ""

# Create cache directory
mkdir -p "$CACHE_DIR"

# Function to run and time a command
time_command() {
    local method="$1"
    shift
    local cmd="$@"

    echo -e "${YELLOW}Testing: $method${NC}"
    echo "Command: $cmd"

    local start=$(date +%s.%N)
    if eval "$cmd" > /dev/null 2>&1; then
        local end=$(date +%s.%N)
        local duration=$(echo "$end - $start" | bc)
        local throughput=$(echo "scale=1; $FILE_SIZE_GB / $duration" | bc)
        echo -e "${GREEN}✓ Success${NC}"
        echo "Duration: ${duration}s"
        echo "Throughput: ${throughput} GB/s"
        echo "$duration" > "$CACHE_DIR/${method}.time"
    else
        local end=$(date +%s.%N)
        local duration=$(echo "$end - $start" | bc)
        echo -e "${RED}✗ Failed${NC}"
        echo "Duration: ${duration}s (failed)"
        echo "999999" > "$CACHE_DIR/${method}.time"  # Large number so it sorts last
    fi
    echo ""
}

# Cleanup function
cleanup() {
    echo -e "${YELLOW}Cleaning up...${NC}"

    # Unmount if mounted
    if mountpoint -q "$MOUNT_POINT" 2>/dev/null; then
        umount "$MOUNT_POINT" 2>/dev/null || fusermount -u "$MOUNT_POINT" 2>/dev/null || true
    fi

    # Remove local copy if exists
    rm -f "$CACHE_DIR/$FILENAME"

    echo "Done"
}
trap cleanup EXIT

echo -e "${BLUE}======================================${NC}"
echo -e "${BLUE}Running Benchmarks${NC}"
echo -e "${BLUE}======================================${NC}"
echo ""

# Benchmark 1: Copy then process (baseline)
time_command "copy-then-process" "
    aws s3 cp '$S3_FILE' '$CACHE_DIR/$FILENAME' && \
    $TOOL $OPERATION '$CACHE_DIR/$FILENAME'
"

# Benchmark 2: FUSE mount (try mountpoint-s3 first, then others)
FUSE_MOUNTED=false

if command -v mount-s3 &> /dev/null; then
    mkdir -p "$MOUNT_POINT"
    mkdir -p "$CACHE_DIR/s3-cache"

    # Mount
    if mount-s3 "$BUCKET" "$MOUNT_POINT" --cache "$CACHE_DIR/s3-cache" --read-only 2>/dev/null; then
        FUSE_MOUNTED=true
        time_command "fuse-mountpoint" "$TOOL $OPERATION '$MOUNT_POINT/$KEY'"
        umount "$MOUNT_POINT" 2>/dev/null || fusermount -u "$MOUNT_POINT" 2>/dev/null || true
    fi
elif command -v goofys &> /dev/null; then
    mkdir -p "$MOUNT_POINT"

    if goofys "$BUCKET" "$MOUNT_POINT" 2>/dev/null; then
        FUSE_MOUNTED=true
        time_command "fuse-goofys" "$TOOL $OPERATION '$MOUNT_POINT/$KEY'"
        umount "$MOUNT_POINT" 2>/dev/null || fusermount -u "$MOUNT_POINT" 2>/dev/null || true
    fi
elif command -v s3fs &> /dev/null; then
    mkdir -p "$MOUNT_POINT"

    if s3fs "$BUCKET" "$MOUNT_POINT" -o use_cache="$CACHE_DIR/s3-cache" -o ro 2>/dev/null; then
        FUSE_MOUNTED=true
        time_command "fuse-s3fs" "$TOOL $OPERATION '$MOUNT_POINT/$KEY'"
        umount "$MOUNT_POINT" 2>/dev/null || fusermount -u "$MOUNT_POINT" 2>/dev/null || true
    fi
fi

if [ "$FUSE_MOUNTED" = false ]; then
    echo -e "${YELLOW}No FUSE tools available (install mount-s3, goofys, or s3fs)${NC}"
    echo ""
fi

# Benchmark 3: Streaming (if supported)
STREAMING_SUPPORTED=false

if [[ "$TOOL" == "samtools" ]] && [[ "$OPERATION" == "flagstat" ]]; then
    STREAMING_SUPPORTED=true
    time_command "streaming" "aws s3 cp '$S3_FILE' - | $TOOL $OPERATION -"
elif [[ "$TOOL" == "samtools" ]] && [[ "$OPERATION" == "view" ]]; then
    STREAMING_SUPPORTED=true
    time_command "streaming" "aws s3 cp '$S3_FILE' - | $TOOL $OPERATION -"
elif [[ "$TOOL" == "head" ]] || [[ "$TOOL" == "wc" ]] || [[ "$TOOL" == "grep" ]]; then
    STREAMING_SUPPORTED=true
    time_command "streaming" "aws s3 cp '$S3_FILE' - | $TOOL $OPERATION"
fi

if [ "$STREAMING_SUPPORTED" = false ]; then
    echo -e "${YELLOW}Streaming not supported for: $TOOL $OPERATION${NC}"
    echo ""
fi

# Benchmark 4: Direct S3 (if supported)
DIRECT_S3_SUPPORTED=false

if [[ "$TOOL" == "samtools" ]] || [[ "$TOOL" == "bcftools" ]]; then
    # Check if htslib has S3 support
    if $TOOL --version 2>&1 | grep -qi "s3"; then
        DIRECT_S3_SUPPORTED=true
        time_command "direct-s3" "$TOOL $OPERATION '$S3_FILE'"
    fi
fi

if [ "$DIRECT_S3_SUPPORTED" = false ]; then
    echo -e "${YELLOW}Direct S3 not supported (tool not compiled with S3 support)${NC}"
    echo ""
fi

# Summary
echo -e "${BLUE}======================================${NC}"
echo -e "${BLUE}Summary${NC}"
echo -e "${BLUE}======================================${NC}"
echo ""

declare -A times
for method in copy-then-process fuse-mountpoint fuse-goofys fuse-s3fs streaming direct-s3; do
    if [ -f "$CACHE_DIR/${method}.time" ]; then
        times[$method]=$(cat "$CACHE_DIR/${method}.time")
    fi
done

# Sort methods by time
sorted_methods=$(for method in "${!times[@]}"; do
    echo "${times[$method]} $method"
done | sort -n | awk '{print $2}')

baseline_time=${times[copy-then-process]:-0}

printf "%-25s %12s %15s %10s\n" "Method" "Duration" "Speedup" "Status"
echo "----------------------------------------------------------------------"

for method in $sorted_methods; do
    time=${times[$method]}

    if [ "$time" = "999999" ]; then
        printf "%-25s %12s %15s %10s\n" "$method" "FAILED" "N/A" "✗"
    else
        duration=$(printf "%.2fs" "$time")

        if [ "$method" = "copy-then-process" ]; then
            speedup="1.00x (baseline)"
        else
            speedup_val=$(echo "scale=2; $baseline_time / $time" | bc)
            speedup="${speedup_val}x"
        fi

        printf "%-25s %12s %15s %10s\n" "$method" "$duration" "$speedup" "✓"
    fi
done

echo ""

# Find fastest method
fastest_method=$(echo "$sorted_methods" | head -n 1)
fastest_time=${times[$fastest_method]}

if [ "$fastest_method" != "copy-then-process" ] && [ "$fastest_time" != "999999" ]; then
    speedup=$(echo "scale=2; $baseline_time / $fastest_time" | bc)
    echo -e "${GREEN}Best method: $fastest_method (${speedup}x faster than copy-then-process)${NC}"
else
    echo -e "${YELLOW}Copy-then-process was fastest (or only working method)${NC}"
fi

echo ""
echo "Recommendation:"
if (( $(echo "$baseline_time / $fastest_time > 1.5" | bc -l) )); then
    echo -e "${GREEN}✓ Direct S3 access provides significant benefit!${NC}"
    echo "  Consider using $fastest_method for your workflows."
elif (( $(echo "$baseline_time / $fastest_time > 1.1" | bc -l) )); then
    echo -e "${YELLOW}⚠ Modest improvement with direct S3 access${NC}"
    echo "  May be worth using for large datasets or repeated access."
else
    echo -e "${RED}✗ No significant benefit from direct S3 access${NC}"
    echo "  Copy-then-process may be adequate for your use case."
    echo "  Consider file size, access patterns, and network speed."
fi

echo ""
