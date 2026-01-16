# BAMS3 Usage Examples

This directory contains practical examples of using BAMS3 for genomics workflows.

## Basic Examples

### 1. Converting BAM to BAMS3

```bash
# Default settings (binary format, zstd compression, 1M chunks, auto-detect workers)
bams3 convert sample.bam sample.bams3

# Custom settings
bams3 convert \
  --chunk-size 512K \
  --compression zstd \
  --format binary \
  --workers 8 \
  input.bam output.bams3
```

### 2. Querying Regions

```bash
# Query specific region
bams3 query sample.bams3 chr1:1000000-2000000

# Count reads only
bams3 query sample.bams3 chr1:1000000-2000000 --count

# Show only first 5 reads
bams3 query sample.bams3 chr1:1000000-2000000 --show 5
```

### 3. Getting Statistics

```bash
# Instant statistics (no file scan needed!)
bams3 stats sample.bams3
```

### 4. Converting Back to BAM

```bash
# Round-trip conversion
bams3 to-bam sample.bams3 reconstructed.bam

# Verify with samtools
samtools flagstat reconstructed.bam
```

## Parallel Processing Examples

See [parallel_workflows.md](./parallel_workflows.md) for detailed parallel processing examples.

## Cloud Storage Examples

See [s3_workflows.md](./s3_workflows.md) for S3 integration examples (v0.3.0+).

## Advanced Examples

See individual files for specific use cases:
- [variant_calling.md](./variant_calling.md) - Variant calling workflows
- [coverage_analysis.md](./coverage_analysis.md) - Coverage calculation
- [population_studies.md](./population_studies.md) - Multi-sample analysis
- [distributed_processing.md](./distributed_processing.md) - Spark/Hadoop integration
