# GitHub Issues and Milestones Plan

This document outlines the complete issue/milestone structure for the project.
Use this to create issues on GitHub.

## Milestones

### Milestone 1: BAMS3 Production Validation (v1.0.0)
**Timeline**: 2 weeks
**Goal**: Validate BAMS3 with real production data and finalize v1.0 release

### Milestone 2: VCFS3 Prototype (v0.1.0)
**Timeline**: 4 weeks
**Goal**: Prove the pattern works for a second format (VCF)

### Milestone 3: Core Library Extraction
**Timeline**: 2 weeks
**Goal**: Extract shared components for reuse across formats

### Milestone 4: Format Expansion (TBD)
**Timeline**: TBD
**Goal**: Research video, medical imaging, or other formats

---

## Labels

Create these labels for organization:

- `type: enhancement` - New features
- `type: bug` - Bug fixes
- `type: documentation` - Documentation improvements
- `type: testing` - Testing and validation
- `type: research` - Research and design
- `priority: critical` - Must have for release
- `priority: high` - Important but not blocking
- `priority: medium` - Nice to have
- `priority: low` - Future consideration
- `format: bams3` - BAMS3 specific
- `format: vcfs3` - VCFS3 specific
- `format: research-video` - Research video specific
- `component: core` - Shared core library
- `component: docs` - Documentation
- `component: benchmarks` - Performance testing
- `good first issue` - Good for new contributors
- `help wanted` - Looking for contributors

---

## Issues for Milestone 1: BAMS3 Production Validation

### Issue #7: Test BAMS3 with Real WGS Data (30x Coverage)
**Labels**: `type: testing`, `priority: critical`, `format: bams3`
**Milestone**: BAMS3 Production Validation (v1.0.0)

**Description**:
Validate BAMS3 with production-scale whole genome sequencing data to ensure it handles real-world workloads.

**Test Dataset**:
- 1000 Genomes Project sample (NA12878 - Genome in a Bottle)
- ~30x coverage, ~1 billion reads
- ~100GB BAM file equivalent

**Test Scenarios**:
1. Full WGS alignment and conversion
   ```bash
   bwa mem -t 32 GRCh38.fa NA12878_R1.fq NA12878_R2.fq | \
     bams3 convert --stdin s3://test-bucket/NA12878.bams3 \
       --workers 32 --sort-buffer 60G
   ```

2. Region extraction and GATK integration
   ```bash
   bams3 to-bam s3://test-bucket/NA12878.bams3 - \
     --region chr17:41196312-41277500 | \
     gatk HaplotypeCaller -I /dev/stdin -R ref.fa -O output.vcf
   ```

3. Nextflow pipeline validation
   ```bash
   nextflow run bams3/nextflow/main.nf \
     --samples samples.csv \
     --reference GRCh38.fa \
     -profile aws
   ```

**Metrics to Collect**:
- [ ] Total conversion time
- [ ] Peak memory usage
- [ ] Disk spill triggered? (how much data)
- [ ] Actual S3 storage cost
- [ ] S3 data transfer costs
- [ ] Query performance (various region sizes)
- [ ] Read count validation (vs traditional BAM)
- [ ] GATK compatibility verification

**Success Criteria**:
- Completes without errors
- Memory usage within expected bounds
- Read counts match traditional workflow
- Cost savings match projections (90%+)
- GATK pipeline produces valid VCF

**Deliverables**:
- Test report with metrics
- Updated documentation with real-world numbers
- Case study write-up

---

### Issue #8: Create GitHub Release v1.0.0
**Labels**: `type: enhancement`, `priority: critical`, `component: docs`
**Milestone**: BAMS3 Production Validation (v1.0.0)

**Description**:
Create official v1.0.0 release to mark production-ready status.

**Tasks**:
- [ ] Create git tag `v1.0.0`
- [ ] Write comprehensive release notes
- [ ] Build binaries for major platforms (Linux, macOS, maybe Windows)
- [ ] Create GitHub release with binaries attached
- [ ] Update version numbers in code
- [ ] Add release badges to README

**Release Notes Should Include**:
- Summary of features
- Performance benchmarks
- Cost savings analysis
- Breaking changes (if any)
- Upgrade guide
- Known limitations
- Links to documentation

---

### Issue #9: Update Main README for v1.0.0 Release
**Labels**: `type: documentation`, `priority: high`, `component: docs`
**Milestone**: BAMS3 Production Validation (v1.0.0)

**Description**:
Update the main README to reflect production-ready status and guide new users.

**Sections to Add/Update**:
- [ ] Add "Production Ready" badge at top
- [ ] Feature highlights with cost savings
- [ ] Quick start guide (5-minute getting started)
- [ ] Architecture diagram
- [ ] Comparison table (BAMS3 vs traditional BAM)
- [ ] Links to all documentation
- [ ] Contributing guidelines
- [ ] Citation information
- [ ] License information
- [ ] Acknowledgments

**Example Quick Start**:
```bash
# Install
go install github.com/scttfrdmn/bams3/cmd/bams3@latest

# Convert BAM to BAMS3
bams3 convert input.bam output.bams3

# Stream from BWA to S3
bwa mem ref.fa R1.fq R2.fq | \
  bams3 convert --stdin s3://bucket/sample.bams3

# Extract region and pipe to GATK
bams3 to-bam s3://bucket/sample.bams3 - --region chr17 | \
  gatk HaplotypeCaller -I /dev/stdin -R ref.fa -O output.vcf
```

---

### Issue #10: Document Real-World Performance Numbers
**Labels**: `type: documentation`, `priority: high`, `component: docs`
**Milestone**: BAMS3 Production Validation (v1.0.0)

**Description**:
Update all documentation with real-world performance numbers from Issue #7 testing.

**Files to Update**:
- [ ] README.md - Performance summary
- [ ] BAM_EXPORT.md - Query performance numbers
- [ ] S3_INTEGRATION.md - Actual S3 costs
- [ ] benchmarks/README.md - Real WGS benchmarks
- [ ] IMPLEMENTATION_SUMMARY.md - Update with v1.0 status

**Include**:
- Actual AWS costs (not projections)
- Performance with 1B reads
- Memory usage patterns
- Query latency for various region sizes
- Comparison with traditional workflow

---

### Issue #11: Write Case Study: 30x WGS Analysis
**Labels**: `type: documentation`, `priority: medium`, `component: docs`
**Milestone**: BAMS3 Production Validation (v1.0.0)

**Description**:
Create a detailed case study showing complete workflow from FASTQ to VCF using BAMS3.

**Content**:
- Step-by-step workflow
- Commands used
- Time/cost breakdown
- Comparison with traditional approach
- Lessons learned
- Best practices

**Format**: Markdown document with code snippets and metrics

**Target audience**: Bioinformaticians considering BAMS3 adoption

---

## Issues for Milestone 2: VCFS3 Prototype

### Issue #12: Design VCFS3 Architecture
**Labels**: `type: research`, `priority: critical`, `format: vcfs3`
**Milestone**: VCFS3 Prototype (v0.1.0)

**Description**:
Design cloud-native VCF format with selective sample and region access.

**Design Goals**:
1. **Sample-centric queries**: Extract single sample without downloading cohort
2. **Region-centric queries**: Query specific genomic region
3. **Combined queries**: Sample X, region Y
4. **Tool compatibility**: Export to standard VCF format

**Key Design Questions**:
- [ ] 1D chunking (genomic position only) or 2D (position × sample)?
- [ ] Sparse genotype encoding (most match reference)
- [ ] How to handle multi-sample metadata?
- [ ] Sample manifest structure
- [ ] Variant annotation storage
- [ ] Index structure for fast queries

**Proposed Architecture**:
```
sample_cohort.vcfs3/
├── _metadata.json          # Cohort metadata
├── _samples.json           # Sample manifest
├── _header.json            # VCF header
├── _index/
│   ├── by_position.idx    # Genomic position index
│   ├── by_sample.idx      # Sample index
│   └── by_variant.idx     # Variant ID index
└── chunks/
    ├── chr1/
    │   ├── 00000000-01000000/
    │   │   ├── sample_001-050.chunk    # Samples 1-50
    │   │   ├── sample_051-100.chunk    # Samples 51-100
    │   │   └── ...
```

**Deliverable**: Design document (VCFS3_DESIGN.md)

---

### Issue #13: Implement VCFS3 Converter (VCF → VCFS3)
**Labels**: `type: enhancement`, `priority: critical`, `format: vcfs3`
**Milestone**: VCFS3 Prototype (v0.1.0)

**Description**:
Implement VCF to VCFS3 converter supporting multi-sample VCFs.

**Requirements**:
- [ ] Parse VCF format (handle BCF binary too?)
- [ ] Extract sample genotypes
- [ ] Chunk by genomic position and sample groups
- [ ] Compress chunks (zstd)
- [ ] Generate metadata and indexes
- [ ] Support streaming input (for large VCFs)
- [ ] Progress reporting

**API**:
```bash
vcfs3 convert input.vcf.gz output.vcfs3 --workers 16
```

**Test with**:
- Small multi-sample VCF (10 samples, chr22 only)
- 1000 Genomes VCF subset (100 samples, single chromosome)

---

### Issue #14: Implement VCFS3 Query (Sample + Region Extraction)
**Labels**: `type: enhancement`, `priority: critical`, `format: vcfs3`
**Milestone**: VCFS3 Prototype (v0.1.0)

**Description**:
Implement selective query functionality for VCFS3.

**Query Types**:
1. **Sample-only**: Extract all variants for sample X
2. **Region-only**: Extract all samples for region Y
3. **Sample + Region**: Extract sample X, region Y
4. **Variant ID**: Find specific variant across all samples

**API**:
```bash
# Single sample, all variants
vcfs3 query cohort.vcfs3 --sample NA12878 -o output.vcf

# Single region, all samples
vcfs3 query cohort.vcfs3 --region chr17:41196312-41277500 -o output.vcf

# Sample + region
vcfs3 query cohort.vcfs3 --sample NA12878 --region chr17 -o output.vcf

# Multiple samples
vcfs3 query cohort.vcfs3 --samples samples.txt --region chr17 -o output.vcf
```

**Performance target**: Download only relevant chunks (>99% reduction)

---

### Issue #15: Implement VCFS3 to VCF Export
**Labels**: `type: enhancement`, `priority: high`, `format: vcfs3`
**Milestone**: VCFS3 Prototype (v0.1.0)

**Description**:
Export VCFS3 back to standard VCF format for tool compatibility.

**Requirements**:
- [ ] Reconstruct VCF header
- [ ] Export selected samples/regions
- [ ] Stream to stdout (for piping)
- [ ] Output to file
- [ ] Support BCF output format
- [ ] Maintain VCF specification compliance

**API**:
```bash
# Export to VCF
vcfs3 to-vcf cohort.vcfs3 output.vcf --sample NA12878

# Stream to bcftools
vcfs3 to-vcf cohort.vcfs3 - --region chr17 | \
  bcftools view -i 'QUAL>30'
```

---

### Issue #16: VCFS3 Benchmark: 1000 Genomes Cohort
**Labels**: `type: testing`, `priority: high`, `format: vcfs3`
**Milestone**: VCFS3 Prototype (v0.1.0)

**Description**:
Benchmark VCFS3 with real 1000 Genomes VCF data.

**Test Dataset**:
- 1000 Genomes Phase 3 VCF
- 2,504 samples
- ~80M variants
- Original VCF: ~500GB (all chromosomes)

**Benchmark Scenarios**:
1. **Conversion**: VCF → VCFS3
2. **Single sample extraction**: 1 sample, all variants
3. **Region query**: All samples, single gene region
4. **Sample + region**: 1 sample, 1 region
5. **Cohort subset**: 100 samples, chromosome 17

**Metrics**:
- Storage size (VCFS3 vs original VCF)
- Conversion time
- Query time for each scenario
- Data transfer (S3 GET size)
- Cost comparison

**Success Criteria**:
- 90%+ reduction in data transfer for selective queries
- Storage size ≤ original VCF
- Query time < 1 minute for single sample

---

### Issue #17: VCFS3 Documentation
**Labels**: `type: documentation`, `priority: high`, `format: vcfs3`
**Milestone**: VCFS3 Prototype (v0.1.0)

**Description**:
Create comprehensive documentation for VCFS3 format and tools.

**Documents to Create**:
- [ ] VCFS3_DESIGN.md - Format specification
- [ ] VCFS3_USAGE.md - User guide with examples
- [ ] VCFS3_BENCHMARKS.md - Performance results
- [ ] VCFS3_TOOLS.md - bcftools/GATK integration

**Include**:
- Format specification
- Usage examples
- Tool integration patterns
- Performance benchmarks
- Cost analysis
- Migration guide (VCF → VCFS3)

---

### Issue #18: S3 Parallel Prefix Optimization for High Throughput
**Labels**: `type: enhancement`, `priority: medium`, `format: vcfs3`, `component: core`
**Milestone**: VCFS3 Prototype (v0.1.0)

**Description**:
Implement hash-based prefix sharding for S3 to enable higher aggregate throughput for large cohorts and concurrent queries.

**Problem**:
Current design uses sequential prefixes (e.g., `chunks/chr1/...`). This limits S3 request rates to ~5,500 requests/second per prefix. For large cohorts with concurrent queries, this becomes a bottleneck:
- VCFS3 with 2,504 samples: 26 chunks per genomic position
- 100 concurrent queries: 2,600 S3 GET requests
- Could hit rate limits with >200 concurrent queries

**Solution - Hash-Based Sharding**:
Distribute chunks across multiple S3 prefixes using hash-based sharding:

**Before** (sequential):
```
chunks/chr1/00000000-01000000/samples_000-099.chunk.zst
chunks/chr1/00000000-01000000/samples_100-199.chunk.zst
```

**After** (hash-sharded):
```
chunks/8a/chr1/00000000-01000000/samples_000-099.chunk.zst
chunks/f2/chr1/00000000-01000000/samples_100-199.chunk.zst
chunks/1d/chr1/01000000-02000000/samples_000-099.chunk.zst
```

**Benefits**:
- 100x higher aggregate throughput (256 prefixes → ~1.4M requests/sec)
- Supports 50,000+ concurrent queries
- Reduced latency under high load
- Better scaling for large cohorts

**Implementation**:
- [ ] Add `ShardingStrategy` interface to core library
- [ ] Implement `HashSharding` (SHA-256 based, configurable bits: 4, 6, 8)
- [ ] Implement `NoSharding` (default, backward compatible)
- [ ] Add `--enable-sharding` flag to converters
- [ ] Update manifest to track sharding configuration
- [ ] Modify query engine to scan multiple prefixes when sharding enabled
- [ ] Document when to enable sharding (cohorts >1000 samples, high query volume)

**Configuration**:
```bash
# Enable for high-throughput scenarios
vcfs3 convert cohort.vcf s3://bucket/cohort.vcfs3 \
    --enable-sharding \
    --sharding-bits 8  # 256 prefixes

# Queries work transparently
vcfs3 query s3://bucket/cohort.vcfs3 --sample NA12878
```

**Trade-offs**:
- Pros: Higher throughput, better scaling, reduced latency
- Cons: More complex directory structure, listing chunks requires scanning multiple prefixes
- Recommendation: Default disabled, enable for large cohorts (>1000 samples)

**Benchmarking**:
Compare query performance with/without sharding:
- Single query: Should be similar (no benefit)
- 100 concurrent queries: Should see improvement
- 1000 concurrent queries: Should see 10x+ improvement

**Success Criteria**:
- Supports 1000+ concurrent queries without S3 rate limiting
- No performance regression for single queries
- Documentation clearly explains when to enable

**References**:
- AWS S3 Performance Guidelines: https://docs.aws.amazon.com/AmazonS3/latest/userguide/optimizing-performance.html
- S3 request rate limits: 5,500 GET/HEAD per second per prefix

---

## Issues for Milestone 3: Core Library Extraction

### Issue #19: Design Shared Core Library Architecture
**Labels**: `type: research`, `priority: high`, `component: core`
**Milestone**: Core Library Extraction

**Description**:
Design shared library that can be used by BAMS3, VCFS3, and future formats.

**Components to Extract**:
1. **Chunking framework**
   - Generic chunking algorithms
   - Chunk size optimization
   - Overlap handling

2. **S3 client**
   - Multipart upload
   - Range requests
   - Credential management
   - Progress tracking

3. **Compression**
   - zstd, gzip wrappers
   - Streaming compression
   - Compression level selection

4. **Indexing**
   - Generic index structures
   - B-tree, R-tree implementations
   - Query optimization

5. **Metadata**
   - JSON schema validation
   - Metadata storage patterns
   - Version management

**Package Structure**:
```
core/
├── chunking/
│   ├── chunker.go          # Generic chunker interface
│   ├── algorithms.go       # Chunking algorithms
│   └── optimizer.go        # Chunk size optimization
├── storage/
│   ├── storage.go          # Storage interface
│   ├── s3.go              # S3 implementation
│   ├── local.go           # Local filesystem
│   └── multipart.go       # Multipart upload
├── compression/
│   ├── compressor.go      # Compression interface
│   └── zstd.go            # zstd implementation
├── index/
│   ├── index.go           # Index interface
│   └── btree.go           # B-tree implementation
└── metadata/
    ├── metadata.go        # Metadata structures
    └── versioning.go      # Version management
```

**Deliverable**: Architecture document (CORE_LIBRARY_DESIGN.md)

---

### Issue #20: Extract S3 Client to Core Library
**Labels**: `type: enhancement`, `priority: high`, `component: core`
**Milestone**: Core Library Extraction

**Description**:
Extract S3 client code from BAMS3 into shared core library.

**Tasks**:
- [ ] Create `core/storage` package
- [ ] Move S3 code from `bams3/pkg/bams3/storage.go`
- [ ] Generalize interfaces (remove BAM-specific code)
- [ ] Add unit tests
- [ ] Update BAMS3 to use core library
- [ ] Update VCFS3 to use core library

**Interface Design**:
```go
type Storage interface {
    ReadFile(path string) ([]byte, error)
    WriteFile(path string, data []byte) error
    ReadRange(path string, offset, length int64) ([]byte, error)
    MultipartUpload(path string, parts [][]byte) error
    List(prefix string) ([]string, error)
    Delete(path string) error
}
```

---

### Issue #21: Extract Compression to Core Library
**Labels**: `type: enhancement`, `priority: medium`, `component: core`
**Milestone**: Core Library Extraction

**Description**:
Extract compression utilities into shared core library.

**Requirements**:
- Support multiple algorithms (zstd, gzip, lz4)
- Streaming compression/decompression
- Configurable compression levels
- Benchmarking utilities

---

### Issue #22: Create Format Design Guide
**Labels**: `type: documentation`, `priority: medium`, `component: docs`
**Milestone**: Core Library Extraction

**Description**:
Document best practices for designing new cloud-native formats.

**Content**:
1. **Design Principles**
   - Chunking strategies
   - Index design
   - Metadata organization
   - Compression selection

2. **Implementation Checklist**
   - Required components
   - Testing requirements
   - Documentation standards
   - Performance targets

3. **Examples**
   - Walk through BAMS3 design decisions
   - Walk through VCFS3 design decisions
   - Common patterns

4. **Tools**
   - Code templates
   - Testing frameworks
   - Benchmark harnesses

**Deliverable**: FORMAT_DESIGN_GUIDE.md

---

## Issues for Milestone 4: Format Expansion (Research)

### Issue #23: Research Video Format - Requirements Gathering
**Labels**: `type: research`, `priority: medium`, `format: research-video`
**Milestone**: Format Expansion

**Description**:
Gather requirements from researchers who work with annotated video data.

**Target Users**:
- Psychology researchers (behavioral coding)
- Linguistics researchers (conversation analysis)
- Education researchers (classroom observation)
- Medical researchers (surgical video analysis)

**Questions to Answer**:
- [ ] What tools do you currently use?
- [ ] What are the pain points?
- [ ] How large are typical datasets?
- [ ] What kinds of queries do you need?
- [ ] What annotation schemas do you use?
- [ ] Do you need cloud storage?
- [ ] Would you pay for a solution?

**Deliverable**: Requirements document with user interviews

---

### Issue #24: Research Video Format - Design Specification
**Labels**: `type: research`, `priority: medium`, `format: research-video`
**Milestone**: Format Expansion

**Description**:
Design cloud-native research video format with integrated annotations.

**Design Goals**:
1. Co-locate video, audio, transcripts, and annotations
2. Time-based indexing for fast queries
3. Annotation-based queries ("find all moments of laughter")
4. Support multiple annotation layers
5. Export to standard formats (MP4, WebM, SRT, CSV)

**Format Structure**:
```
study.rvid/
├── _metadata.json          # Study metadata
├── _schema.json            # Annotation schema
├── _index/
│   ├── by_time.idx        # Time-based index
│   ├── by_annotation.idx   # Annotation index
│   └── by_speaker.idx     # Speaker index
└── chunks/
    ├── 00000-00030/       # 0-30 seconds
    │   ├── video.mp4      # Video chunk
    │   ├── audio.aac      # Audio chunk
    │   ├── transcript.json # Speech-to-text
    │   ├── frames.json     # Frame annotations
    │   └── behaviors.json  # Behavioral codes
```

**Deliverable**: Format specification document

---

### Issue #25: Imaging Formats - Requirements Survey
**Labels**: `type: research`, `priority: medium`, `format: imaging`
**Milestone**: Format Expansion

**Description**:
Survey imaging community needs across medical, microscopy, and research imaging.

**Target Domains**:

**A) Medical Imaging (DICOM)**
- Radiology (CT, MRI, PET)
- Pathology (whole slide imaging)
- Cardiology (echo, angiography)

**B) Microscopy**
- CryoEM (electron microscopy for structure)
- Super-resolution (STORM, PALM, SIM)
- Light sheet (developmental biology)
- High-content screening (drug discovery)
- Live-cell imaging (time-lapse)

**C) Research Imaging**
- Behavioral video (as discussed)
- Satellite/drone (extends COG)

**Key Questions**:
- [ ] What are current storage formats? (DICOM, MRC, TIFF, HDF5, Zarr, OME-TIFF)
- [ ] Dataset sizes? (GB to TB per study/experiment)
- [ ] Query patterns? (slices, channels, timepoints, ROIs)
- [ ] Cloud adoption blockers?
- [ ] Tool ecosystems? (ImageJ, Napari, CellProfiler, RELION)
- [ ] Compliance needs? (HIPAA for medical)

**Pain Points by Domain**:

**CryoEM**:
- Movie stacks: 1000s of frames, need subsets for motion correction
- Micrographs: TB datasets, need random access for particle picking
- Current: Download entire dataset to local cluster
- Tools: RELION, cryoSPARC (not cloud-native)

**Super-Resolution**:
- Raw data: 10,000s of frames per acquisition
- Need: Specific time windows for reconstruction
- Current: TIFF stacks (no selective access)

**Live-Cell Imaging**:
- 4D data: X, Y, Z, time, multiple channels
- Need: Slice by time, channel, z-plane
- Current: OME-TIFF (better than TIFF but still monolithic)

**Deliverable**: Requirements report prioritizing which imaging domain to tackle first

---

### Issue #26: CryoEM Format Design - Prototype
**Labels**: `type: research`, `priority: medium`, `format: cryoem`
**Milestone**: Format Expansion

**Description**:
Design cloud-native format for CryoEM movie stacks and micrographs.

**Current Problems**:
- Movie stacks: 1000s of frames (50-100GB), need subsets for motion correction
- Micrographs: TB datasets, must download all for particle picking
- No cloud-native processing (RELION, cryoSPARC require local data)

**Design Goals**:
1. Chunk movie stacks for selective frame access
2. Tile micrographs for region-based processing
3. Co-locate metadata (CTF, particle picks)
4. Enable cloud-based processing pipelines
5. Export to MRC/TIFF for tool compatibility

**Format Structure** (draft):
```
cryoem_session.crys3/
├── _metadata.json          # Session metadata
├── _ctf/                   # CTF parameters
├── _particles/             # Particle picks
└── data/
    ├── movies/             # Movie stacks
    │   ├── movie_0001/
    │   │   ├── frames_000-099.chunk
    │   │   ├── frames_100-199.chunk
    │   │   └── metadata.json
    └── micrographs/        # Processed micrographs
        ├── tile_000_000.chunk
        └── ...
```

**Target Users**:
- Structural biologists
- CryoEM facilities
- Cloud-based processing services

**Deliverable**: CryoEM format specification + prototype converter

---

### Issue #27: ML Embeddings Format - Feasibility Study
**Labels**: `type: research`, `priority: low`, `format: embeddings`
**Milestone**: Format Expansion

**Description**:
Evaluate whether chunked embeddings format would be valuable.

**Research Questions**:
- Do vector DBs fully solve this problem?
- Are there use cases for S3-native embeddings?
- What about offline/research scenarios?
- Would clustering-based chunking actually help?

**Deliverable**: Feasibility report with recommendation

---

## Priority Order for Issue Creation

**Create immediately** (blocking v1.0.0):
1. Issue #7 - Test BAMS3 with Real WGS Data
2. Issue #8 - Create GitHub Release v1.0.0
3. Issue #9 - Update Main README
4. Issue #10 - Document Real-World Performance

**Create after Issue #7 completes**:
5. Issue #11 - Write Case Study

**Create to start VCFS3 work**:
6. Issue #12 - Design VCFS3 Architecture
7. Issue #13 - Implement VCFS3 Converter
8. Issue #14 - Implement VCFS3 Query
9. Issue #15 - Implement VCFS3 Export
10. Issue #16 - VCFS3 Benchmark
11. Issue #17 - VCFS3 Documentation

**Create for core library**:
12. Issue #18 - Design Core Library
13. Issue #19 - Extract S3 Client
14. Issue #20 - Extract Compression
15. Issue #21 - Format Design Guide

**Create for research** (imaging focus after genomics):
16. Issue #24 - Imaging Formats Requirements Survey
17. Issue #25 - CryoEM Format Design
18. Issue #22 - Research Video Requirements
19. Issue #23 - Research Video Design
20. Issue #26 - ML Embeddings Feasibility (low priority)

---

## Additional Issues to Consider

### Issue #28: Add Telemetry/Analytics (Optional)
**Labels**: `type: enhancement`, `priority: low`

**Description**:
Add optional telemetry to understand usage patterns (with user consent).

**Metrics to Track**:
- Command usage frequency
- Dataset sizes processed
- Error rates
- Performance metrics

**Privacy**: Fully anonymized, opt-in only

---

### Issue #29: Create Homebrew Formula
**Labels**: `type: enhancement`, `priority: medium`

**Description**:
Create Homebrew formula for easy installation on macOS.

```bash
brew install bams3
```

---

### Issue #30: Create Bioconda Package
**Labels**: `type: enhancement`, `priority: high`

**Description**:
Submit to Bioconda for bioinformatics community distribution.

---

### Issue #31: Write Announcement Blog Post
**Labels**: `type: documentation`, `priority: high`

**Description**:
Write comprehensive blog post announcing BAMS3 v1.0.0.

**Target Audience**: Bioinformatics community

**Content**:
- Problem statement
- Solution overview
- Performance benchmarks
- Cost savings
- Getting started guide
- Roadmap

**Distribution**:
- Personal blog
- Biostars
- Reddit r/bioinformatics
- Twitter/X

---

### Issue #32: Submit to bioRxiv (Future)
**Labels**: `type: documentation`, `priority: low`

**Description**:
Write academic preprint describing BAMS3 format and performance.

**Include**:
- Format specification
- Benchmarks with real data
- Cost analysis
- Tool compatibility study
- Adoption recommendations

---

## Strategic Direction

**Phase 1** (Current - 3 months): **Genomics Excellence**
- BAMS3 production validation
- VCFS3 prototype and validation
- Core library extraction
- Community building

**Phase 2** (3-6 months): **Imaging Expansion**
- Priority: CryoEM (structural biology has resources and need)
- Secondary: Medical DICOM (huge market but regulatory complexity)
- Tertiary: Research video (niche but underserved)

**Rationale for Imaging Focus**:
1. ✅ Similar technical challenges (large files, selective access)
2. ✅ Proven pattern from genomics (90%+ savings)
3. ✅ Underserved markets (no good cloud-native solutions)
4. ✅ High-value users (research institutions, pharma)
5. ✅ Core library reuse validates multi-domain approach

**CryoEM specifically**:
- TB-scale datasets (bigger than genomics!)
- Well-funded domain (pharma, NIH)
- Conservative tooling (RELION, cryoSPARC) ripe for disruption
- Cloud adoption starting but tools lacking

---

This provides a complete roadmap from v1.0.0 release through multi-domain expansion.
