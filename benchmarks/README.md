# BAMS3 Performance Benchmarks

Comprehensive test harness for comparing BAMS3 against traditional BAM workflows.

## Overview

This benchmark suite provides objective, reproducible comparisons across multiple dimensions:

1. **End-to-End Performance**: Pipeline execution time
2. **Resource Usage**: Memory, CPU, disk I/O
3. **Storage Efficiency**: Disk space, compression ratios
4. **Data Integrity**: Read counts, sorting correctness
5. **Scalability**: Performance across dataset sizes
6. **Cloud Economics**: AWS cost analysis

## Quick Start

### Basic Comparison

```bash
cd benchmarks
chmod +x compare-workflows.sh
./compare-workflows.sh
```

### Different Dataset Sizes

```bash
# Small dataset (10K reads, ~2 minutes)
DATASET_SIZE=small ./compare-workflows.sh

# Medium dataset (100K reads, ~10 minutes)
DATASET_SIZE=medium ./compare-workflows.sh

# Large dataset (1M reads, ~30 minutes)
DATASET_SIZE=large ./compare-workflows.sh
```

### Custom Configuration

```bash
# Specify workers and memory
DATASET_SIZE=medium WORKERS=16 MEMORY_LIMIT=32G ./compare-workflows.sh
```

## Benchmark Scripts

### 1. `compare-workflows.sh`

**Purpose:** Complete pipeline comparison

**Measures:**
- Pipeline execution time (total and per-step)
- Memory usage (peak RSS)
- Disk usage (intermediate files, final output)
- Read counts and data integrity
- Compression efficiency

**Output:**
- Detailed metrics in `results/comparison_*.txt`
- Markdown summary in `results/summary_*.md`
- Individual step timings in `results/*_time.txt`

**Traditional Workflow Steps:**
1. BWA alignment → SAM
2. SAM → BAM conversion
3. BAM sorting
4. BAM indexing

**BAMS3 Workflow Steps:**
1. BWA alignment → BAMS3 (zero-copy, coordinate sorted)

### 2. `query-benchmark.sh` (Coming Soon)

**Purpose:** Compare query performance

**Tests:**
- Small region queries (1KB - 1MB)
- Large region queries (1MB - 100MB)
- Whole chromosome queries
- Cold vs warm cache performance

**Metrics:**
- Query latency
- Data transferred
- Memory usage

### 3. `scalability-test.sh` (Coming Soon)

**Purpose:** Test scaling characteristics

**Tests:**
- Worker count scaling (1, 2, 4, 8, 16, 32 workers)
- Memory scaling (different buffer sizes)
- Dataset size scaling (1K to 10M reads)

### 4. `cloud-cost-analysis.sh` (Coming Soon)

**Purpose:** AWS cost calculations

**Analyzes:**
- Storage costs (S3 Standard, Glacier)
- Data transfer costs (same-region, cross-region)
- Compute costs (EC2, AWS Batch)
- Request costs (GET, PUT)

## Understanding Results

### Performance Metrics

**Speedup Factor:**
- `> 2.0x`: Significantly faster
- `1.5-2.0x`: Moderately faster
- `1.0-1.5x`: Comparable performance
- `< 1.0x`: Investigate (may need tuning)

**Expected Results:**
- Small datasets (< 100K reads): 1.5-2x speedup
- Medium datasets (100K-1M reads): 2-3x speedup
- Large datasets (> 1M reads): 3-5x speedup

*Speedup increases with scale due to zero intermediate files*

### Storage Efficiency

**Compression Ratios:**
- **Target:** 60-80% of traditional BAM size
- Better compression for:
  - High coverage depth (more repetitive data)
  - Longer reads (better base-level compression)
  - Sorted data (BAMS3 is always sorted)

**Disk Space Savings:**
- Eliminates intermediate SAM (often 3-5x BAM size)
- Single format (no separate index file)
- Compressed metadata

### Data Integrity

**Validation Checks:**
- ✓ Read counts must match exactly
- ✓ Mapped vs unmapped counts must match
- ✓ Reference distribution should match
- ✓ Coordinate sorting verified

**If counts differ:**
- Check for BWA version differences
- Verify no data truncation
- Review warning messages

## Results Interpretation

### Example Output

```
Performance Metrics:
===================

Pipeline Time:
  Traditional: 120s
  BAMS3:       45s
  Speedup: 2.67x faster

Data Integrity:
  Traditional reads: 100000
  BAMS3 reads:       100000
  ✓ Read counts match

Disk Usage:
  Traditional: 256000 KB (includes SAM, BAM, sorted BAM, index)
  BAMS3:       51200 KB (single format)
  Savings: 80% reduction
  Compression ratio: 20%

Key Advantages:
===============
  ✓ Zero intermediate files (no SAM, no unsorted BAM)
  ✓ Single-step pipeline (vs 4-step traditional)
  ✓ Coordinate sorted during alignment
  ✓ Queryable chunks (no full-file downloads)
  ✓ Cloud-ready (direct S3 upload capability)
```

### What This Means

1. **2.67x Speedup**: BAMS3 completes in 1/3 the time
2. **80% Space Savings**: Uses 1/5 the disk space
3. **Data Integrity**: All reads preserved correctly
4. **Simplified Workflow**: Single command vs 4 commands

## Advanced Usage

### Profiling with perf

```bash
# Profile traditional workflow
perf record -g ./compare-workflows.sh

# Analyze
perf report
```

### Memory Profiling

```bash
# Use valgrind for detailed memory analysis
valgrind --tool=massif bams3 convert --stdin output.bams3 < input.sam
```

### I/O Profiling

```bash
# Monitor disk I/O
iostat -x 1 > io_stats.txt &
./compare-workflows.sh
pkill iostat
```

## Automated Regression Testing

Run benchmarks in CI/CD:

```yaml
# .github/workflows/benchmark.yml
name: Performance Benchmark

on: [pull_request]

jobs:
  benchmark:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Run benchmarks
        run: |
          cd benchmarks
          DATASET_SIZE=small ./compare-workflows.sh
      - name: Upload results
        uses: actions/upload-artifact@v2
        with:
          name: benchmark-results
          path: benchmarks/results/
```

## Reproducibility

### Environment

Record environment for reproducibility:

```bash
# System info
uname -a > results/system_info.txt
lscpu >> results/system_info.txt
free -h >> results/system_info.txt

# Tool versions
bwa 2>&1 | head -n 3 >> results/versions.txt
samtools --version >> results/versions.txt
./bams3 --version >> results/versions.txt
```

### Random Seed

For synthetic data generation:

```bash
# Set random seed for reproducible data
export RANDOM_SEED=42
./compare-workflows.sh
```

## Common Issues

### Out of Memory

```bash
# Reduce memory usage
MEMORY_LIMIT=4G WORKERS=2 ./compare-workflows.sh
```

### Disk Space

```bash
# Cleanup previous results
rm -rf results/*.sam results/*.bam
```

### Performance Variance

Run multiple iterations and average:

```bash
for i in {1..5}; do
  DATASET_SIZE=small ./compare-workflows.sh
done

# Analyze variance in results/
```

## Contributing

To add new benchmarks:

1. Create script in `benchmarks/`
2. Follow naming convention: `*-benchmark.sh`
3. Output results to `results/` directory
4. Document in this README
5. Add to regression test suite

## References

- [Issue #6: Performance Benchmarking](https://github.com/scttfrdmn/bams3/issues/6)
- [Traditional BAM Workflow](https://www.htslib.org/)
- [BWA Documentation](http://bio-bwa.sourceforge.net/)
- [SAMtools Documentation](http://www.htslib.org/doc/samtools.html)

## Results Archive

Benchmark results are stored in `results/` with timestamps:

```
results/
├── comparison_small_20260116_120000.txt    # Raw metrics
├── summary_20260116_120000.md               # Markdown report
├── traditional_1_align_time.txt             # Step timings
├── traditional_2_sam_to_bam_time.txt
├── traditional_3_sort_time.txt
├── traditional_4_index_time.txt
└── bams3_1_align_and_convert_time.txt
```

Keep historical results to track performance over time and detect regressions.
