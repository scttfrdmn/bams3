# BAMS3/VCFS3 Demonstrations

Interactive demos showcasing key features of cloud-native genomics formats.

## Overview

This directory contains self-contained demos that can be run quickly to demonstrate specific features. Each demo is designed to:
- Run in <5 minutes
- Use small synthetic datasets (no large downloads)
- Show clear before/after comparisons
- Include visual output (charts, tables)
- Be presentable to stakeholders

## Available Demos

### 1. Basic Workflow Demo (`01-basic-workflow/`)
**Time**: 2 minutes
**Shows**: Complete BAMS3 workflow from alignment to variant calling

```bash
cd demos/01-basic-workflow
./run-demo.sh
```

**Demonstrates**:
- BWA alignment streaming directly to BAMS3
- Zero intermediate BAM files
- Region extraction
- GATK integration
- Storage savings (96% reduction)

**Output**: Side-by-side comparison table

---

### 2. S3 Sharding Performance Demo (`02-sharding-performance/`)
**Time**: 3 minutes
**Shows**: Throughput difference with/without sharding

```bash
cd demos/02-sharding-performance
./run-demo.sh
```

**Demonstrates**:
- Directory structure comparison (sequential vs sharded)
- Simulated concurrent queries (1, 10, 100, 1000 queries)
- Request rate measurements
- Visual chart showing bottleneck without sharding
- 100× throughput improvement

**Output**:
- Directory tree visualization
- Performance chart (ASCII or HTML)
- Summary table

---

### 3. VCFS3 Selective Access Demo (`03-vcfs3-selective-access/`)
**Time**: 3 minutes
**Shows**: Cost savings from selective sample/region queries

```bash
cd demos/03-vcfs3-selective-access
./run-demo.sh
```

**Demonstrates**:
- Create small multi-sample VCF (10 samples, chr22 subset)
- Convert to VCFS3 with 2D chunking
- Single sample extraction (downloads 1/10th of data)
- Region extraction (downloads 1/1000th of data)
- Sample + region (downloads 1/10,000th of data)
- Cost comparison

**Output**:
- Data transfer bar chart
- Cost savings table
- Query time comparison

---

### 4. Cloud vs Traditional Workflow Demo (`04-cloud-vs-traditional/`)
**Time**: 3 minutes
**Shows**: Complete workflow comparison

```bash
cd demos/04-cloud-vs-traditional
./run-demo.sh
```

**Demonstrates**:
- Traditional: Download full BAM (100GB) → extract region → analyze
- BAMS3: Stream region (5MB) → analyze
- Time and cost comparison
- Multiple query scenario (10 regions)

**Output**:
- Workflow diagram (ASCII art)
- Time/cost comparison table
- Cumulative savings chart

---

### 5. Multi-Format Demo (`05-multi-format/`)
**Time**: 4 minutes
**Shows**: BAMS3, VCFS3, and core library reuse

```bash
cd demos/05-multi-format
./run-demo.sh
```

**Demonstrates**:
- Same data processed through multiple formats
- Shared core library components
- Consistent query patterns across formats
- Format-specific optimizations

**Output**:
- Format comparison table
- Code reuse visualization
- Performance comparison

---

### 6. Production Pipeline Demo (`06-production-pipeline/`)
**Time**: 5 minutes
**Shows**: End-to-end Nextflow pipeline

```bash
cd demos/06-production-pipeline
./run-demo.sh
```

**Demonstrates**:
- FASTQ → BAMS3 → Variant calling
- Parallel sample processing
- S3 integration throughout
- Resource usage monitoring
- MultiQC report generation

**Output**:
- Pipeline execution report
- Resource usage graphs
- QC metrics dashboard

---

## Demo Requirements

### Software
- Docker (recommended) or local installation:
  - BWA 0.7.17+
  - samtools 1.18+
  - GATK 4.0+ (optional for some demos)
  - BAMS3/VCFS3 tools
  - Python 3.8+ (for visualizations)

### Hardware
- Minimal: 4 CPU cores, 8GB RAM, 5GB disk
- All demos use synthetic or small public datasets

### Cloud (Optional)
- AWS account for S3 demos
- Free tier sufficient for all demos
- Estimated cost: <$0.50 per complete demo run

---

## Running All Demos

```bash
# Run complete demo suite
./run-all-demos.sh

# Generates consolidated report: demos/DEMO_REPORT.html
```

**Output includes**:
- Performance metrics
- Cost analysis
- Visual comparisons
- Executive summary

---

## Demo Data

All demos use:
1. **Synthetic data**: Generated on-the-fly (small, fast)
2. **Public data**: Small subsets from 1000 Genomes
3. **Cached data**: Pre-generated for consistency

No large downloads required. Total demo data: <100MB.

---

## Creating New Demos

See `DEMO_TEMPLATE/` for demo structure:

```
demos/XX-demo-name/
├── README.md           # Demo description
├── run-demo.sh         # Main demo script
├── setup.sh            # One-time setup
├── data/               # Small test datasets
├── scripts/            # Demo-specific scripts
└── expected-output/    # Expected results
```

**Guidelines**:
- Keep demos <5 minutes
- Use clear visual output
- Include cleanup script
- Document every step
- Test on clean system

---

## Presentation Mode

For live demonstrations:

```bash
# Run with commentary
./run-demo.sh --verbose

# Slow down output for visibility
./run-demo.sh --slow

# Pause at key points
./run-demo.sh --interactive
```

---

## CI/CD Integration

All demos run in CI to ensure they stay working:

```yaml
# .github/workflows/demos.yml
- name: Run Demos
  run: |
    cd demos
    ./run-all-demos.sh --ci-mode
```

Demos are tested on every commit to catch regressions.

---

## Video Recordings

Pre-recorded demo videos available at:
- YouTube: [BAMS3 Demos Playlist](#)
- Website: [Project demos page](#)

Each demo includes:
- Narrated walkthrough (3-5 min)
- Slides explaining the feature
- Live terminal session
- Results visualization

---

## Demo Datasets

### Synthetic Chromosome 22
- 10,000 reads
- Generated with `wgsim`
- ~5MB SAM/BAM
- Included in repo

### 1000 Genomes Subset
- NA12878, chr22:10000000-11000000 (1Mbp)
- ~100K reads
- Downloaded from S3 (public, no auth)
- Cached locally

### Multi-sample VCF
- 10 samples
- 1000 variants
- chr22 only
- Generated from 1000 Genomes

---

## Metrics Collected

Each demo tracks:
- Execution time
- Memory usage (peak)
- Disk usage
- Network transfer (S3)
- File sizes
- Cost estimates

Results stored in `demos/results/` for comparison over time.

---

## Troubleshooting

**Demo fails to start**:
```bash
# Check requirements
./check-requirements.sh

# Install missing tools
./install-demo-tools.sh
```

**Slow performance**:
- Demos use small datasets; should be fast
- Check available CPU/RAM
- Try Docker version (pre-configured)

**S3 demos fail**:
- Check AWS credentials: `aws configure`
- Verify S3 bucket access
- Try local-only demos first

---

## Demo Schedule

**Quick overview** (10 minutes):
1. Basic workflow (2 min)
2. S3 sharding (3 min)
3. VCFS3 selective access (3 min)
4. Cost comparison (2 min)

**Technical deep-dive** (30 minutes):
1. All demos (20 min)
2. Q&A (10 min)

**Executive presentation** (15 minutes):
1. Cloud vs traditional (3 min)
2. Cost analysis (5 min)
3. Production pipeline (5 min)
4. Executive summary (2 min)

---

## Future Demos

Planned:
- [ ] Migration tools demo (existing BAM → BAMS3)
- [ ] Multi-cloud demo (S3, Azure, GCS)
- [ ] Real-time streaming demo
- [ ] Large cohort demo (1000 Genomes full dataset)
- [ ] Integration with Terra/DNAnexus
- [ ] Performance tuning demo
- [ ] Debugging and troubleshooting demo
