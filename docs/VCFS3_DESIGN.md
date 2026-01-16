# VCFS3 Architecture Design

Cloud-native VCF format for efficient sample and region queries in cohort genomics.

## Problem Statement

**Current VCF limitations:**

1. **Row-oriented storage**: Variants stored by genomic position
   - Query "sample X, all variants" → Must scan entire file
   - Query "variant Y, all samples" → Must scan entire file

2. **Monolithic files**: Single file per cohort
   - 1000 samples × 80M variants = 500GB VCF
   - Extract 1 sample → Download all 500GB
   - No selective access

3. **Cloud unfriendly**: Designed for local filesystem
   - Tabix indexing requires full file download first
   - No HTTP range request optimization
   - Column queries impossible

**Real-world impact:**
```
Cohort: 10,000 samples, 100M variants = 5TB VCF

Query: "Give me patient X's variants"
Traditional: Download 5TB, extract 1 sample = hours, $450 transfer cost
Needed: Download patient X only = seconds, <$1

Query: "Who has variant rs123?"
Traditional: Download 5TB, scan all samples = hours
Needed: Query variant index = seconds
```

---

## Design Goals

1. **Sample-centric queries**: Extract single sample without cohort download
2. **Region-centric queries**: Extract genomic region across samples
3. **Combined queries**: Sample X, region Y (most common in practice)
4. **Variant-level queries**: "Which samples have this variant?"
5. **Tool compatibility**: Export to standard VCF format
6. **Cloud-native**: Optimized for S3 with selective chunk access
7. **Compression**: Match or beat gzipped VCF
8. **Scalability**: Handle 100K+ sample cohorts

---

## Core Architecture

### 1. Two-Dimensional Chunking

**Key innovation**: Chunk by **genomic position × sample groups**

```
Traditional VCF (row-oriented):
Position | Sample1 | Sample2 | ... | Sample10000
chr1:100 | 0/1     | 0/0     | ... | 1/1
chr1:200 | 0/0     | 0/1     | ... | 0/1
...

VCFS3 (2D chunked):
Chunk[chr1:0-1Mb][samples 1-100]
Chunk[chr1:0-1Mb][samples 101-200]
...
Chunk[chr1:1-2Mb][samples 1-100]
...
```

**Benefits**:
- Sample query: Download only sample's chunks
- Region query: Download only region's chunks
- Sample + region: Download single chunk
- Scales independently in both dimensions

### 2. Chunk Dimensions

**Genomic dimension** (rows):
- Chunk size: 1 Mbp (megabase pairs)
- Rationale: Balance query granularity vs overhead
  - Too small (100kb): Many chunks, high overhead
  - Too large (10Mb): Large downloads for small regions
  - 1Mb: Sweet spot for most analyses

**Sample dimension** (columns):
- Chunk size: 100 samples per chunk
- Rationale: Balance flexibility vs efficiency
  - Too small (10 samples): Many chunks for cohort queries
  - Too large (1000 samples): Large downloads for single sample
  - 100 samples: Reasonable overhead (~1% for 10K cohort)

**Example: 1000 Genomes (2,504 samples)**:
```
Chromosomes: 22 autosomes + X + Y + MT = 25 refs
Total genomic size: ~3,000 Mbp = 3,000 chunks (genomic)
Sample chunks: 2,504 samples / 100 = 26 chunks (sample)
Total chunks: 3,000 × 26 = 78,000 chunks

Average chunk size: 500GB / 78,000 = 6.4 MB per chunk
```

### 3. Sparse Genotype Encoding

**Observation**: Most genotypes match reference (0/0)

**Strategy**: Store only non-reference genotypes
```
Dense storage (all genotypes):
chr1:100 | 0/0 | 0/0 | 0/1 | 0/0 | 0/0 | ... | 1/1

Sparse storage (non-ref only):
chr1:100 | samples: [3:0/1, 10000:1/1]
         | 99.98% are reference (not stored)
```

**Compression ratio**: 10-50× for common variants

**Implementation**:
```json
{
  "position": 100,
  "ref": "A",
  "alt": ["G"],
  "genotypes": [
    {"sample_idx": 3, "gt": "0/1", "gq": 99},
    {"sample_idx": 10000, "gt": "1/1", "gq": 99}
  ]
}
```

### 4. Sample Manifest

**Purpose**: Map sample names to chunk indices

```json
{
  "samples": [
    {"id": "NA12878", "chunk_index": 0, "offset_in_chunk": 0},
    {"id": "NA12891", "chunk_index": 0, "offset_in_chunk": 1},
    ...
    {"id": "HG00101", "chunk_index": 1, "offset_in_chunk": 0}
  ],
  "chunk_size": 100
}
```

**Fast lookup**: Sample name → chunk index → download only needed chunks

---

## Format Specification

### Directory Structure

```
cohort.vcfs3/
├── _metadata.json              # Dataset metadata
├── _samples.json               # Sample manifest
├── _header.vcf                 # VCF header (for export)
├── _index/
│   ├── position.idx            # Genomic position index
│   ├── sample.idx              # Sample name index
│   └── variant.idx             # Variant ID index (rsID, etc.)
└── chunks/
    ├── chr1/
    │   ├── 00000000-01000000/
    │   │   ├── samples_000-099.chunk.zst
    │   │   ├── samples_100-199.chunk.zst
    │   │   ├── samples_200-299.chunk.zst
    │   │   └── ...
    │   ├── 01000000-02000000/
    │   │   └── ...
    │   └── ...
    ├── chr2/
    │   └── ...
    └── ...
```

### Metadata Format (_metadata.json)

```json
{
  "format": "vcfs3",
  "version": "0.1.0",
  "created": "2026-01-16T12:00:00Z",
  "created_by": "vcfs3-go",
  "source": {
    "file": "1000genomes_phase3.vcf.gz",
    "format": "VCF",
    "version": "4.2"
  },
  "statistics": {
    "total_samples": 2504,
    "total_variants": 81271745,
    "references": ["chr1", "chr2", ..., "chrM"],
    "variant_types": {
      "snp": 77818198,
      "indel": 3453547
    }
  },
  "chunking": {
    "genomic_chunk_size": 1000000,
    "sample_chunk_size": 100,
    "total_chunks": 78012,
    "compression": "zstd"
  }
}
```

### Sample Manifest (_samples.json)

```json
{
  "sample_count": 2504,
  "sample_chunk_size": 100,
  "samples": [
    {
      "id": "HG00096",
      "population": "GBR",
      "superpopulation": "EUR",
      "sex": "male",
      "chunk_index": 0,
      "offset": 0
    },
    {
      "id": "HG00097",
      "population": "GBR",
      "superpopulation": "EUR",
      "sex": "female",
      "chunk_index": 0,
      "offset": 1
    },
    ...
  ]
}
```

### Chunk Format (Binary)

**Each chunk contains**:
- Chunk header (position range, sample range)
- Variant index (position → offset in chunk)
- Variant records (sparse genotypes)

**Binary format** (protobuf or custom):
```
ChunkHeader:
  - reference: string (chr1)
  - start_pos: int64 (0)
  - end_pos: int64 (1000000)
  - sample_start: int32 (0)
  - sample_end: int32 (99)
  - variant_count: int32

VariantIndex:
  - positions: []int64 (sorted)
  - offsets: []int64 (byte offset in chunk)

VariantRecords:
  - For each variant:
    - position: int64
    - ref: string
    - alt: []string
    - info: map[string]string
    - format: []string
    - sparse_genotypes: []SparseGT

SparseGT:
  - sample_idx: int32 (relative to chunk)
  - gt: string (0/1, 1/1, etc.)
  - fields: map[string]string (GQ, DP, etc.)
```

### Position Index (_index/position.idx)

**Purpose**: Fast lookup of which chunks contain a genomic region

```json
{
  "chr1": [
    {"start": 0, "end": 1000000, "chunks": ["chunks/chr1/00000000-01000000/"]},
    {"start": 1000000, "end": 2000000, "chunks": ["chunks/chr1/01000000-02000000/"]},
    ...
  ],
  "chr2": [...],
  ...
}
```

### Sample Index (_index/sample.idx)

**Purpose**: Fast lookup of which chunks contain a sample

```json
{
  "HG00096": {
    "chunk_index": 0,
    "chunks": [
      "chunks/chr1/00000000-01000000/samples_000-099.chunk.zst",
      "chunks/chr1/01000000-02000000/samples_000-099.chunk.zst",
      ...
    ]
  },
  ...
}
```

---

## Query Algorithms

### Query 1: Single Sample, All Variants

**Example**: "Give me all variants for patient HG00096"

**Algorithm**:
```python
1. Read _samples.json → find HG00096 → chunk_index = 0
2. Read _index/sample.idx → get all chunks for sample 0
3. For each chunk:
   - Download chunk file
   - Extract genotypes for sample offset 0
   - Convert to VCF records
4. Output VCF
```

**Data transfer**:
- Traditional: 500GB full VCF
- VCFS3: ~200 MB (0.04% of data!)
- **Savings: 99.96%**

### Query 2: Single Region, All Samples

**Example**: "Give me all samples for chr17:41196312-41277500 (BRCA1)"

**Algorithm**:
```python
1. Read _index/position.idx → find chunks overlapping region
   → chr17/41000000-42000000/
2. For each sample chunk in that genomic chunk:
   - Download chunk file
   - Extract variants in region
3. Merge genotypes across sample chunks
4. Output VCF
```

**Data transfer**:
- Traditional: 500GB full VCF
- VCFS3: ~150 MB (26 sample chunks × 6 MB each)
- **Savings: 99.97%**

### Query 3: Single Sample + Single Region

**Example**: "Give me variants for HG00096 in BRCA1"

**Algorithm**:
```python
1. Find sample chunk: _samples.json → chunk_index = 0
2. Find genomic chunk: _index/position.idx → chr17/41000000-42000000/
3. Download single chunk: chr17/41000000-42000000/samples_000-099.chunk.zst
4. Extract sample 0, region BRCA1
5. Output VCF
```

**Data transfer**:
- Traditional: 500GB full VCF
- VCFS3: ~6 MB (single chunk)
- **Savings: 99.9988%**

### Query 4: Variant Lookup

**Example**: "Which samples have rs123456?"

**Algorithm**:
```python
1. Read _index/variant.idx → find position of rs123456
2. Identify chunks containing that position
3. Download relevant chunks
4. For each chunk:
   - Find variant at position
   - Extract non-reference samples
5. Return list of samples with variant
```

**Use case**: GWAS, case-control studies

---

## Compression Strategy

### 1. Sparse Encoding

**Before compression**: Only store non-ref genotypes
- Typical common variant: 1-5% ALT frequency
- Store: ~250 samples (1-5% of 10K) instead of 10K
- **10-20× reduction before compression**

### 2. Delta Encoding

**Genomic positions**: Store as deltas
```
Positions: [100, 150, 175, 200]
Deltas:    [100, +50, +25, +25]
```
**Saves**: ~50% on position storage (smaller integers compress better)

### 3. Dictionary Compression

**Genotypes**: Limited vocabulary (0/0, 0/1, 1/1, ./.)
- Build dictionary: 0→0, 1→0/1, 2→1/1, 3→./.
- Store as bytes instead of strings
- **5-10× reduction**

### 4. zstd Compression

**Apply zstd to chunks**:
- Level 3 (fast) for real-time conversion
- Level 19 (max) for archival
- **Additional 2-3× reduction**

**Combined**: 10× (sparse) × 2× (dict) × 2.5× (zstd) = **50× total compression**

**Reality check**: Expect 20-30× vs uncompressed, 5-10× vs gzipped VCF

---

## Implementation Plan

### Phase 1: Converter (VCF → VCFS3)

**Input**: Standard VCF/BCF file
**Output**: VCFS3 directory structure

**Algorithm**:
```python
1. Parse VCF header → extract samples, contigs
2. Create _metadata.json, _samples.json
3. Initialize chunk writers (one per genomic×sample chunk)
4. Stream VCF records:
   - For each variant:
     - Determine genomic chunk
     - Extract genotypes
     - Encode as sparse
     - Write to appropriate chunk files
5. Finalize chunks (compress, write indexes)
6. Write position/sample/variant indexes
```

**CLI**:
```bash
vcfs3 convert input.vcf.gz output.vcfs3 \
  --workers 16 \
  --genomic-chunk-size 1000000 \
  --sample-chunk-size 100 \
  --compression zstd
```

### Phase 2: Query Tool

**CLI**:
```bash
# Single sample
vcfs3 query cohort.vcfs3 --sample HG00096 -o output.vcf

# Single region
vcfs3 query cohort.vcfs3 --region chr17:41196312-41277500 -o output.vcf

# Sample + region
vcfs3 query cohort.vcfs3 \
  --sample HG00096 \
  --region chr17:41196312-41277500 \
  -o output.vcf

# Multiple samples
vcfs3 query cohort.vcfs3 \
  --samples samples.txt \
  --region chr17 \
  -o output.vcf

# Variant lookup
vcfs3 query cohort.vcfs3 \
  --variant rs123456 \
  --output-samples-only
```

### Phase 3: VCF Export

**CLI**:
```bash
# Export to VCF
vcfs3 to-vcf cohort.vcfs3 output.vcf.gz \
  --sample HG00096 \
  --region chr17

# Stream to bcftools
vcfs3 to-vcf cohort.vcfs3 - --region chr17 | \
  bcftools view -i 'QUAL>30'

# Export to BCF (binary VCF)
vcfs3 to-vcf cohort.vcfs3 output.bcf \
  --format bcf \
  --sample HG00096
```

### Phase 4: Tool Integration

**bcftools plugin**:
```bash
bcftools +vcfs3 cohort.vcfs3 \
  --sample HG00096 \
  --region chr17 | \
  bcftools stats
```

**GATK SelectVariants**:
```bash
# VCFS3 exports on-the-fly
gatk SelectVariants \
  -V <(vcfs3 to-vcf cohort.vcfs3 - --sample HG00096) \
  -L chr17 \
  -O output.vcf
```

---

## Performance Projections

### 1000 Genomes (2,504 samples, 81M variants)

**Storage**:
- Original VCF.gz: 500 GB
- VCFS3: 50-100 GB (5-10× compression)
- **Savings: 80-90%**

**Query Performance**:

| Query Type | Traditional | VCFS3 | Speedup |
|------------|-------------|-------|---------|
| 1 sample, all variants | 500GB download, 20 min | 200MB download, 10 sec | 120× |
| 1 region, all samples | 500GB download, 20 min | 150MB download, 8 sec | 150× |
| 1 sample, 1 region | 500GB download, 20 min | 6MB download, 1 sec | 1200× |

**Cost (S3 data transfer, us-east-1)**:

| Query Type | Traditional | VCFS3 | Savings |
|------------|-------------|-------|---------|
| 1 sample, all variants | $45 | $0.018 | 99.96% |
| 1 region, all samples | $45 | $0.014 | 99.97% |
| 1 sample, 1 region | $45 | $0.0005 | 99.999% |

---

## Validation Strategy

### Test Suite

1. **Correctness Tests**:
   - Convert VCF → VCFS3 → VCF
   - Diff original vs round-trip (should be identical)
   - Test with: 1KG, gnomAD, TOPMed samples

2. **Performance Tests**:
   - Benchmark queries on 1KG data
   - Measure compression ratios
   - Profile memory usage

3. **Scalability Tests**:
   - Test with 10K, 50K, 100K samples
   - Validate chunk distribution
   - Ensure query time stays constant

4. **Integration Tests**:
   - Export to bcftools → validate output
   - Export to GATK → run SelectVariants
   - Export to Plink → run association tests

---

## Open Questions

### 1. Annotation Storage

**Question**: Where to store INFO fields (population frequencies, functional annotations)?

**Options**:
- A) Store with variants (in each chunk)
  - Pro: Self-contained chunks
  - Con: Duplication across sample chunks
- B) Separate annotation file
  - Pro: No duplication
  - Con: Extra download for annotated queries

**Recommendation**: Option A for now (simpler), optimize later if needed

### 2. Phasing Information

**Question**: How to handle phased genotypes (0|1 vs 0/1)?

**Options**:
- Store phasing bit separately
- Store as string (current VCF way)

**Recommendation**: Store as string initially (simpler, compatible)

### 3. Multi-Allelic Variants

**Question**: Handle multi-allelic variants (>2 alleles)?

**Recommendation**: Store as-is (match VCF spec), no special handling needed

### 4. Metadata Updates

**Question**: Allow updating metadata (e.g., adding samples) without rewrite?

**Recommendation**: v1.0 - No updates (immutable). v2.0 - Consider append-only

---

## Success Metrics

**VCFS3 is successful if**:

1. **Storage**: ≤ original gzipped VCF size
2. **Sample query**: 99%+ reduction in data transfer
3. **Region query**: 99%+ reduction in data transfer
4. **Combined query**: 99.9%+ reduction in data transfer
5. **Conversion**: Processes 1KG in <1 hour on standard hardware
6. **Compatibility**: Exports valid VCF, works with bcftools/GATK
7. **Scalability**: Handles 100K+ samples without performance degradation

---

## Next Steps

1. **Prototype** (Week 1-2):
   - Implement basic converter (VCF → VCFS3)
   - Test with small VCF (chr22 only)
   - Validate chunk structure

2. **Query Engine** (Week 3-4):
   - Implement sample query
   - Implement region query
   - Implement combined query

3. **Benchmarking** (Week 5-6):
   - Test with 1000 Genomes
   - Measure compression ratios
   - Validate cost savings

4. **Tool Integration** (Week 7-8):
   - VCF export
   - bcftools plugin
   - GATK integration

5. **Documentation** (Week 8):
   - Format specification
   - User guide
   - Migration guide

**Total timeline**: 8 weeks to v0.1.0 prototype

---

This design provides a solid foundation for cloud-native cohort genomics. The 2D chunking strategy is the key innovation enabling efficient sample and region queries.
