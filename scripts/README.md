# Benchmarking Scripts

Tools to benchmark different S3 access methods and compare performance.

## Quick Start

### Simple Shell Script (Recommended for Quick Tests)

```bash
# Make executable
chmod +x quick-benchmark.sh

# Run benchmark
./quick-benchmark.sh s3://my-bucket/sample.bam samtools flagstat

# Or with other tools
./quick-benchmark.sh s3://my-bucket/variants.vcf.gz bcftools view
./quick-benchmark.sh s3://my-bucket/reads.fastq.gz head -n 1000
```

**What it tests:**
1. ✅ Copy-then-process (baseline)
2. ✅ FUSE mounting (auto-detects available tools)
3. ✅ Streaming (if operation supports it)
4. ✅ Direct S3 (if tool has S3 support)

**Output:**
- Duration for each method
- Throughput (GB/s)
- Speedup vs baseline
- Recommendation

### Python Script (More Detailed)

```bash
# Install dependencies
pip install pyyaml boto3

# Run benchmark
python benchmark.py \
    --s3-file s3://my-bucket/sample.bam \
    --operation flagstat \
    --output results.json

# Results saved to results.json
```

**Features:**
- More detailed metrics
- Multiple operations
- JSON output for analysis
- Memory usage tracking (TODO)

## Example Output

```
======================================
S3 Access Method Benchmark
======================================

File:      s3://my-bucket/sample.bam
Tool:      samtools flagstat
Cache dir: /tmp/benchmark-cache

Getting file size...
Size: 12.45 GB

======================================
Running Benchmarks
======================================

Testing: copy-then-process
Command: aws s3 cp ... && samtools flagstat ...
✓ Success
Duration: 145.23s
Throughput: 0.09 GB/s

Testing: fuse-mountpoint
Command: samtools flagstat /mnt/sample.bam
✓ Success
Duration: 52.18s
Throughput: 0.24 GB/s

Testing: streaming
Command: aws s3 cp - | samtools flagstat -
✓ Success
Duration: 48.91s
Throughput: 0.25 GB/s

Testing: direct-s3
Command: samtools flagstat s3://my-bucket/sample.bam
✗ Failed (tool not compiled with S3 support)

======================================
Summary
======================================

Method                    Duration        Speedup       Status
----------------------------------------------------------------------
streaming                 48.91s          2.97x         ✓
fuse-mountpoint           52.18s          2.78x         ✓
copy-then-process         145.23s         1.00x (baseline) ✓
direct-s3                 FAILED          N/A           ✗

Best method: streaming (2.97x faster than copy-then-process)

Recommendation:
✓ Direct S3 access provides significant benefit!
  Consider using streaming for your workflows.
```

## Supported Operations

### BAM/SAM Files
- `flagstat` - Count flags in BAM file
- `view-head` - View first records
- `idxstats` - Index statistics
- `count-reads` - Count total reads
- `region-query` - Query specific region (requires index)

### VCF/BCF Files
- `view` - View variants
- `query` - Query specific fields

### Generic
- `head` - First N lines
- `wc` - Line count
- `grep` - Search patterns

## Understanding Results

### Duration
Total time including:
- **Copy-then-process**: Download time + processing time
- **FUSE mount**: Mount setup + processing (first access may be slower due to cache warming)
- **Streaming**: Download + process simultaneously
- **Direct S3**: Tool's native S3 access

### Throughput
- Calculated as: File Size / Duration
- Higher is better
- Network limited? Check if < network bandwidth

### Speedup
- Relative to copy-then-process baseline
- **>2x**: Significant improvement, strongly recommend
- **1.5-2x**: Good improvement, recommend for large datasets
- **1.1-1.5x**: Modest improvement, consider for repeated access
- **<1.1x**: Minimal benefit, copy-then-process may be adequate

## Factors Affecting Performance

### 1. File Size
- **Small files (<1GB)**: Copy overhead is minimal
- **Medium files (1-10GB)**: Significant benefit from direct access
- **Large files (>10GB)**: Maximum benefit, avoid copying

### 2. Access Pattern
- **Sequential full read**: Streaming is optimal
- **Random access**: FUSE with good cache or direct S3
- **Multiple small queries**: Direct S3 with indexes
- **One-time access**: Copy may be fine
- **Repeated access**: FUSE with cache or copy once

### 3. Network
- **Fast network (10+ Gbps)**: Streaming competitive with local
- **Moderate network (1-10 Gbps)**: Still better than copy-then-process
- **Slow network (<1 Gbps)**: Benefit reduced but still present

### 4. Operation Type
- **CPU-bound** (compression, decompression): Streaming hides download
- **I/O-bound** (simple reads): Direct access faster
- **Seeks required** (region queries): Need FUSE or direct S3 with index

## Advanced Usage

### Benchmark Multiple Files

```bash
#!/bin/bash
# benchmark-batch.sh

for file in s3://my-bucket/samples/*.bam; do
    echo "Benchmarking $file"
    ./quick-benchmark.sh "$file" samtools flagstat
    echo ""
done
```

### Compare Different Operations

```bash
# Test different operations on same file
for op in flagstat idxstats count-reads; do
    echo "Operation: $op"
    python benchmark.py \
        --s3-file s3://bucket/sample.bam \
        --operation $op \
        --output "results-${op}.json"
done
```

### Benchmark with Different Cache Sizes

```bash
# Test FUSE mount with different cache sizes
for cache_size in 10000 50000 100000; do
    echo "Cache size: ${cache_size}MB"

    mount-s3 my-bucket /mnt/s3 \
        --cache /tmp/cache \
        --max-cache-size $cache_size

    time samtools flagstat /mnt/s3/sample.bam

    umount /mnt/s3
done
```

## Interpreting Results for Your Workflow

### Scenario 1: One-Time Analysis
**If you process each file once:**
- Small benefit from direct S3 (unless file is huge)
- Copy-then-process acceptable
- Consider streaming to avoid storage

### Scenario 2: Repeated Access
**If you access files multiple times:**
- FUSE with cache: Excellent performance after first access
- Direct S3: Good for indexed queries
- Copy once: Also acceptable

### Scenario 3: Large-Scale Pipeline
**Processing 100s-1000s of files:**
- Direct S3 essential (avoid massive storage)
- FUSE + cache for hot data
- Streaming for sequential processing
- **Don't copy everything locally!**

### Scenario 4: Interactive Analysis
**Genome browsers, manual inspection:**
- FUSE mount: Most convenient
- Direct S3 if tool supports (IGV can use S3)
- Caching critical for responsiveness

## Troubleshooting

### "Command not found: mount-s3"
```bash
# Install mountpoint-for-amazon-s3
# macOS:
brew install mountpoint-s3

# Linux:
# See examples/fuse-mounting/README.md for installation
```

### "Permission denied"
```bash
# Check AWS credentials
aws sts get-caller-identity

# Configure if needed
aws configure
```

### "Streaming failed"
Some operations can't be streamed:
- Random access (region queries with index)
- Multiple passes over data
- Operations requiring file size upfront

### Slow FUSE performance
```bash
# Increase cache size
mount-s3 bucket /mnt/s3 \
    --cache /fast-disk/cache \
    --max-cache-size 200000  # 200GB

# Use faster storage for cache (NVMe)
# Increase thread count
mount-s3 bucket /mnt/s3 --thread-count 32
```

### Results don't match expectations
- Check network speed: `wget http://speedtest.tele2.net/100MB.zip`
- Check S3 transfer rates: `aws s3 cp s3://bucket/large-file /dev/null`
- Verify no other processes using bandwidth
- Try from EC2 in same region as S3 bucket

## Next Steps

1. Run benchmarks on your real data
2. Identify which method works best for your workflow
3. Document your findings
4. Update your pipelines to use optimal method
5. Consider:
   - Data reorganization (partitioning)
   - Format optimization
   - Tool modifications for native S3 support

## Related Documentation

- [FUSE Mounting Guide](../examples/fuse-mounting/README.md)
- [Streaming Guide](../examples/streaming/README.md)
- [Format Optimization](../docs/format-optimization.md)
