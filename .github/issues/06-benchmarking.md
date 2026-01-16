# Performance benchmarking and optimization

**Labels:** performance, documentation, priority-high
**Priority:** High

## Goal
Systematically measure and document performance characteristics to validate cost/speed improvements.

## Benchmark Categories

### 1. Conversion Performance

**Traditional Workflow:**
```bash
time bwa mem ref.fa reads.fq > output.sam
time samtools view -bS output.sam > output.bam
time samtools sort output.bam -o sorted.bam
time samtools index sorted.bam
du -h output.sam output.bam sorted.bam
```

**BAMS3 Workflow:**
```bash
time bwa mem ref.fa reads.fq | bams3 convert --stdin output.bams3
du -h output.bams3
```

**Metrics to Measure:**
- [ ] Total time (end-to-end)
- [ ] Peak memory usage
- [ ] Disk space required (including intermediates)
- [ ] CPU utilization
- [ ] Throughput (reads/sec, MB/sec)

### 2. Query Performance

**Test Cases:**
```bash
# Small region (1KB)
time bams3 query sample.bams3 chr1:1000000-1001000 --count
time samtools view -c sample.bam chr1:1000000-1001000

# Medium region (1MB)
time bams3 query sample.bams3 chr1:1000000-2000000 --count
time samtools view -c sample.bam chr1:1000000-2000000

# Large region (10MB)
time bams3 query sample.bams3 chr1:1000000-11000000 --count
time samtools view -c sample.bam chr1:1000000-11000000

# Whole chromosome
time bams3 query sample.bams3 chr22 --count
time samtools view -c sample.bam chr22
```

**Metrics:**
- [ ] Query latency (cold vs warm)
- [ ] Data transferred (selective chunk download)
- [ ] Memory usage
- [ ] Cache effectiveness

### 3. Storage Efficiency

**Test Files:**
- WGS 30x coverage (~50GB BAM)
- Exome sequencing (~5GB BAM)
- RNA-Seq (~10GB BAM)
- Single chromosome

**Metrics:**
```bash
# Size comparison
du -sh original.bam
du -sh output.bams3
# Calculate compression ratio

# Format breakdown
du -sh output.bams3/_metadata.json
du -sh output.bams3/_header.json
du -sh output.bams3/data/
```

### 4. Scalability Testing

**Worker Count Scaling:**
```bash
for workers in 1 2 4 8 16 32; do
  time bams3 convert --stdin test.bams3 \
    --workers $workers \
    < input.sam
done
```

**Memory Scaling:**
```bash
for buffer in 512M 1G 2G 4G 8G 16G; do
  time bams3 convert --stdin test.bams3 \
    --sort-buffer $buffer \
    < input.sam
done
```

**Dataset Size Scaling:**
- 1M reads
- 10M reads
- 100M reads
- 1B reads (if feasible)

### 5. Cloud Cost Analysis

**AWS Cost Breakdown:**

**Traditional Approach:**
```
Download BAM: 50GB × $0.09/GB = $4.50
Compute: r5.4xlarge × 2 hours × $1.008/hr = $2.02
Upload BAM: 50GB × $0.09/GB = $4.50
Storage: 50GB × $0.023/GB/month = $1.15/month
Per-query download: 50GB × $0.09/GB = $4.50

First analysis: $11.02 + $1.15/month
Each subsequent: $4.50
```

**BAMS3 Approach:**
```
Stream in-region: 50GB × $0.00/GB (same region) = $0.00
Compute: r5.4xlarge × 1 hour × $1.008/hr = $1.01
Upload BAMS3: 8GB × $0.00/GB (same region) = $0.00
Storage: 8GB × $0.023/GB/month = $0.18/month
Per-query: 50MB × $0.09/GB = $0.0045 (selective chunks)

First analysis: $1.01 + $0.18/month
Each subsequent: $0.0045
```

**Calculate:**
- [ ] Cost per sample (first run)
- [ ] Cost per sample (subsequent analyses)
- [ ] Storage cost over 1 year
- [ ] Total cost for 100 samples
- [ ] Break-even point

## Test Infrastructure

### Benchmark Suite
**Path:** `benchmarks/` (new directory)

```
benchmarks/
├── conversion_test.sh       # Conversion benchmarks
├── query_test.sh           # Query benchmarks
├── storage_test.sh         # Storage benchmarks
├── scalability_test.sh     # Scaling benchmarks
├── cloud_cost.sh           # Cost calculations
├── results/                # Benchmark results
│   ├── conversion.json
│   ├── query.json
│   └── ...
└── analysis/               # Analysis scripts
    ├── plot_performance.py
    └── generate_report.py
```

### Automated Reporting
- [ ] Generate performance graphs
- [ ] Create comparison tables
- [ ] Export to markdown/PDF
- [ ] CI integration (track performance over time)

## Deliverables

### Documentation
- [ ] `docs/performance.md` - Detailed benchmark results
- [ ] `docs/cloud-costs.md` - Cost analysis and comparisons
- [ ] `README.md` - Performance summary section
- [ ] `benchmarks/README.md` - How to run benchmarks

### Visualizations
- [ ] Conversion time vs dataset size
- [ ] Query latency by region size
- [ ] Memory usage patterns
- [ ] Cost comparison charts
- [ ] Scalability curves

## Acceptance Criteria
- [ ] Comprehensive benchmark suite implemented
- [ ] Results documented with graphs and tables
- [ ] Cost analysis shows 75%+ savings
- [ ] Performance shows 10x+ speedup for queries
- [ ] Conversion speed competitive with traditional tools
- [ ] Results reproducible

## Related
- Hyperfine for benchmarking: https://github.com/sharkdp/hyperfine
- AWS Cost Calculator: https://calculator.aws/

## Notes
Critical for demonstrating value proposition and validating claims.
