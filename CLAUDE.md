# BAMS3 Project Vision & Context

## Mission

**Fix how researchers use the cloud for bioinformatics.**

Researchers currently use cloud infrastructure inefficiently - downloading terabytes to analyze megabytes, storing massive intermediate files, running sequential workflows on expensive parallel hardware, and paying 100-1000x more than necessary.

BAMS3 demonstrates the correct approach: cloud-optimized data formats that enable fast, efficient, cost-conscious workflows without materially altering scientific results.

## Core Principles

1. **Fast** - 10-100x faster than traditional approaches
2. **Efficient** - Minimal data transfer, no wasteful intermediates
3. **Cost-conscious** - 75-99% cost reduction
4. **Scientifically identical** - Bit-for-bit identical results
5. **Easy adoption** - Works with existing tools and workflows
6. **Smart defaults** - Non-technical users get optimal performance
7. **Expert configurability** - Power users can tune everything

## The Problem

### Traditional Cloud Bioinformatics (Wrong)

```
Sequencer → FASTQ (250 GB in S3)
  ↓ Download entire file to compute (250 GB transfer, $22.50, 30 min)
  ↓ bwa mem → 150 GB SAM file on disk
  ↓ samtools sort → 50 GB BAM file on disk (delete SAM)
  ↓ samtools markdup → 50 GB BAM file on disk
  ↓ Upload BAM to S3 (50 GB transfer, $4.50, 10 min)
  ↓ Later: Download BAM for analysis (50 GB transfer, $4.50, 10 min)
  ↓ GATK HaplotypeCaller (sequential, per-chromosome)

Cost per sample: ~$50 first run, ~$10 per subsequent analysis
Time: 8-12 hours first run, 2-3 hours per analysis
Waste: 200+ GB temporary files, 100+ GB unnecessary transfers
```

### Cloud-Optimized Bioinformatics (Correct)

```
Sequencer → FASTQ (250 GB in S3)
  ↓ Stream directly to compute (no download, same region = free)
  ↓ bwa mem | bams3 convert --stdin s3://bucket/sample.bams3
      (no temp files, direct upload while converting, 8 GB final)
  ↓ Parallel variant calling (24 chromosomes simultaneously)
      bams3 to-bam s3://bucket/sample.bams3 - --region chr${N} | GATK
  ↓ Merge VCFs

Cost per sample: ~$12 first run, ~$0.01 per subsequent analysis
Time: 1-2 hours first run, 5-10 minutes per analysis
Waste: 0 GB temporary files, 0 GB unnecessary transfers
Savings: 75% cost, 85% time on first run; 99.9% cost on subsequent analyses
```

## Solution Architecture

### BAMS3 Format Design

**Core insight:** Store alignment data as independent, self-contained chunks in object storage.

```
sample.bams3/
├── _metadata.json          # Dataset statistics, chunk inventory (15 KB)
├── _header.json           # SAM header (2 KB)
├── _index/
│   └── spatial.json       # Spatial index for O(log n) chunk lookup (50 KB)
└── data/
    ├── chr1/
    │   ├── 000000000-001048576.chunk    # 1 MB genomic region
    │   ├── 001048576-002097152.chunk
    │   └── ...
    ├── chr2/
    └── ...
```

**Why this works:**
- Each chunk is independent (no locks, no coordination)
- Selective downloads (query chr1 → download only chr1 chunks)
- Parallel processing (process all chunks simultaneously)
- Cloud-native (optimized for object storage access patterns)
- Compressed (zstd 6x ratio, smaller than BAM)

### Key Technical Components

1. **Streaming Conversion** - Zero-copy pipeline from aligner to S3
2. **S3 Integration** - Direct read/write, no local storage needed
3. **Spatial Indexing** - Fast chunk lookup for queries
4. **Parallel Processing** - Linear scaling with worker count
5. **Smart Defaults** - Auto-configure for optimal performance
6. **Tool Integration** - Works with existing bioinformatics tools

## Demonstration Scope

### Target: Complete Reference Workflows

Demonstrate cloud-optimized approaches across multiple domains:

#### 1. Whole Genome Sequencing (Primary)
```
FASTQ → Alignment → BAMS3 → Variant Calling → VCF
Tools: BWA, BAMS3, GATK
Benchmark: 30x WGS, compare traditional vs cloud-optimized
Show: 75% cost reduction, 85% time reduction, identical variants
```

#### 2. RNA-Seq (Secondary)
```
FASTQ → Alignment → BAMS3 → Gene Counting → Differential Expression
Tools: STAR, BAMS3, HTSeq/featureCounts, DESeq2
Benchmark: 100 samples, repeated analyses
Show: 99.9% transfer cost reduction on repeated analyses
```

#### 3. Microbiome/Metagenomics (Secondary)
```
Reads → Taxonomic Classification → BAMS3 → Abundance/Diversity Analysis
Tools: Kraken2, BAMS3, custom analysis
Benchmark: Multiple samples, per-organism queries
Show: Query specific taxa without downloading entire dataset
```

#### 4. Proteomics (Future)
```
Similar pattern: chunk by scan range, query specific m/z windows
Needs: Cloud-optimized mzML-equivalent format
```

### Deliverable: Research Publication

**Title:** *Cloud-Optimized Data Formats for Bioinformatics: A 100× Cost Reduction Without Compromising Scientific Validity*

**Thesis:** Traditional bioinformatics formats (BAM, SAM, CRAM) were designed for local storage. Cloud-optimized alternatives deliver identical scientific results at 10-100x lower cost and time.

**Evidence:**
- Complete workflows: FASTQ → final analysis
- Multiple domains: genomics, transcriptomics, metagenomics
- Quantitative metrics: cost, time, storage, transfer
- Scientific validation: bit-for-bit identical results
- Reproducibility: public Docker images, versioned workflows

**Impact:**
- Researchers save 75-99% on cloud costs
- Analyses run 10-100x faster
- Enables population-scale studies (1000s of samples)
- Sets pattern for other bioinformatics domains

## Technical Roadmap

### Phase 1: Complete Genomics Workflow (4-6 weeks)
- [ ] Streaming conversion (stdin → BAMS3 → S3, zero-copy)
- [ ] Region-based BAM export (BAMS3 → BAM for GATK)
- [ ] BED file queries (for gene counting)
- [ ] Nextflow reference pipeline (complete WGS workflow)
- [ ] Traditional pipeline (for comparison)
- [ ] Benchmark infrastructure (cost/time tracking)
- [ ] Initial results (NA12878 30x WGS)

### Phase 2: Tool Integration (4-6 weeks)
- [ ] htslib plugin (native samtools support)
- [ ] Python library (bams3-py, pysam-compatible API)
- [ ] R package (Bioconductor integration)
- [ ] GATK wrapper (gatk-bams3)
- [ ] Performance optimization (streaming compression, prefetching)

### Phase 3: Multi-Domain Workflows (4-6 weeks)
- [ ] RNA-Seq workflow (STAR → BAMS3 → counts)
- [ ] Microbiome workflow (Kraken2 → BAMS3 → abundance)
- [ ] Benchmarks across all domains
- [ ] Cost/performance analysis

### Phase 4: Publication & Dissemination (4-6 weeks)
- [ ] Write paper (methods, results, discussion)
- [ ] Generate figures (cost comparisons, performance scaling)
- [ ] Document best practices (migration guide)
- [ ] Public release (GitHub, Docker, documentation)
- [ ] Present to stakeholders (AWS, UCLA, collaborators)

## Success Metrics

### Technical Metrics
- **Query speed:** 100x faster than traditional BAM (1-2s vs 120s)
- **Storage:** 80% smaller (6x compression + no redundancy)
- **Cost:** 75-99% reduction depending on workflow
- **Scalability:** Linear scaling to 1000+ parallel workers
- **Compatibility:** Works with 95%+ of existing tools (via conversion/streaming)

### Adoption Metrics
- **Ease of use:** Single command conversion, auto-configuration
- **Scientific validity:** Bit-for-bit identical results
- **Integration:** Works with existing pipelines (Nextflow, Snakemake)
- **Documentation:** Complete workflows, migration guides

### Impact Metrics
- **Research enabled:** Population-scale studies (1000s of samples)
- **Cost savings:** $5,000-$50,000 per large study
- **Time savings:** Hours instead of days for analysis
- **Cloud efficiency:** 10-100x better resource utilization

## Design Decisions

### Smart Defaults Philosophy

Every configuration parameter must have an intelligent default:

```bash
# This should "just work" optimally for most users
bams3 convert --stdin output.bams3 < input.sam

# Auto-detects:
# - Workers: runtime.NumCPU()
# - Buffer: 25% of available RAM
# - Chunk size: 1MB (optimal for most queries)
# - Compression: zstd level 3 (best speed/ratio)
# - Sort buffer: 8GB or available RAM, whichever is smaller
```

Power users can override everything:
```bash
bams3 convert --stdin output.bams3 \
  --workers 32 \
  --buffer 16G \
  --sort-buffer 32G \
  --chunk-size 2M \
  --compression-level 5 \
  < input.sam
```

### Graceful Degradation

System must work under constrained resources:
- Low memory → spill sort to disk
- Slow network → reduce parallel uploads
- Disk full → fail fast with clear error
- S3 errors → retry with exponential backoff

### Error Messages Must Be Actionable

Bad:
```
Error: out of memory
```

Good:
```
Error: sort buffer full (8 GB)
Options:
  1. Reduce buffer: --sort-buffer 4G
  2. Enable disk spill: --spill-to-disk (default: enabled)
  3. Increase memory: use instance with more RAM
Current memory usage: 7.8 GB / 8.0 GB
```

## Integration Strategy

### Tier 1: Essential (Must Work)
- samtools (via htslib plugin)
- GATK (via region-based conversion)
- BWA/STAR/Bowtie2 (via streaming conversion)
- IGV (via local proxy or plugin)

### Tier 2: High Value (Should Work)
- Python (native library: bams3-py)
- R (native package: BAMS3)
- bedtools (via streaming)
- deepTools (native coverage calculation)

### Tier 3: Nice to Have (Can Work)
- All other BAM consumers (via conversion)
- Custom tools (via Python/R libraries)
- Workflow managers (Nextflow, Snakemake, Cromwell)

### Integration Pattern

```
Level 1: Pipe-through (works today)
  bams3 query sample.bams3 chr1:1M-2M | tool --stdin

Level 2: Conversion (works today)
  bams3 to-bam sample.bams3 temp.bam --region chr20
  tool temp.bam

Level 3: Native integration (future)
  tool sample.bams3 --region chr20:1M-2M
```

## Non-Goals (For Now)

- **Not** replacing CRAM for archival storage
- **Not** supporting every edge case (focus on common workflows)
- **Not** building a startup/commercial product (research project)
- **Not** working through AWS HealthOmics (internal constraints)
- **Not** waiting for standards bodies (move fast, standardize later)

## Development Practices

### Code Quality
- Clear, documented code (this may be reference implementation)
- Comprehensive error handling (fail gracefully)
- Performance instrumentation (measure everything)
- Testing (unit tests, integration tests, benchmarks)

### Documentation
- Vision (this document)
- Architecture (design decisions)
- API reference (godoc)
- User guides (examples, tutorials)
- Migration guides (traditional → cloud-optimized)

### Project Management
- Use GitHub Projects for tracking (NOT markdown documents)
- Milestones for phases
- Issues for specific tasks
- Labels for categorization (enhancement, bug, documentation, etc.)

### Benchmarking
- Every claim must be measurable
- Reproduce all benchmarks (Docker, versioned workflows)
- Track: time, cost, storage, transfer, compute resources
- Compare: traditional vs cloud-optimized (apples-to-apples)

## Target Audience

### Primary: Bioinformatics Researchers
- Need: Fast, cheap analyses
- Pain: Cloud costs are exploding
- Solution: Cloud-optimized formats save 75-99%

### Secondary: Cloud Providers
- Need: Reference architectures for genomics
- Pain: Customers complain about costs
- Solution: Show how to use cloud correctly

### Tertiary: Tool Developers
- Need: Patterns for cloud-native tools
- Pain: Unclear how to adapt tools for cloud
- Solution: BAMS3 as reference implementation

## Long-Term Vision

**Near-term (6-12 months):**
- BAMS3 proven for genomics workflows
- Complete benchmarks published
- Early adopters (UCLA, AWS customers)

**Mid-term (1-2 years):**
- BAMS3 adopted by major biobanks
- Other formats follow pattern (proteomics, imaging)
- Tool ecosystem supports cloud-optimized formats
- GA4GH consideration as standard

**Long-term (3-5 years):**
- Cloud-optimized formats are default for new projects
- Legacy formats (BAM) used only for archives
- 10-100x cost reduction is table stakes
- Enables billion-genome databases

## Why This Matters

**Scientific Impact:**
- Democratizes large-scale genomics (smaller labs can afford it)
- Enables real-time analysis (minutes not hours)
- Accelerates discovery (more analyses, faster iteration)

**Economic Impact:**
- Saves researchers $5,000-$500,000 per study
- Reduces cloud provider transfer revenue (but increases compute/storage)
- Makes cloud genomics economically viable

**Technical Impact:**
- Shows correct cloud architecture patterns
- Reference implementation for other domains
- Pushes tool ecosystem toward cloud-native designs

---

**This is not just a file format. It's a demonstration of how bioinformatics should be done in the cloud.**

**Fast. Efficient. Cost-conscious. Correct.**
