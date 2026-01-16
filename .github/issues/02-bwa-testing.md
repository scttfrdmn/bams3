# End-to-end testing with BWA aligner

**Labels:** testing, priority-high
**Priority:** High

## Goal
Validate complete zero-copy pipeline from sequencer to BAMS3 format.

## Requirements
- BWA aligner installed
- Reference genome (e.g., GRCh38)
- FASTQ test data (or subset)

## Test Scenarios

### 1. Basic Pipeline
```bash
bwa mem ref.fa reads.fq | bams3 convert --stdin sample.bams3
```
**Verify:**
- [ ] No errors during conversion
- [ ] All reads processed (count matches input)
- [ ] Chunks created correctly
- [ ] Output queryable

### 2. Direct S3 Upload (when S3 integration ready)
```bash
bwa mem ref.fa reads.fq | bams3 convert --stdin s3://bucket/sample.bams3
```
**Verify:**
- [ ] Streams directly to S3
- [ ] No local intermediate files
- [ ] Correct S3 object structure

### 3. Large Dataset
Test with full chromosome or whole genome:
```bash
bwa mem -t 16 ref.fa large_reads.fq | \
  bams3 convert --stdin output.bams3 \
    --workers 16 \
    --sort-buffer 16G
```
**Verify:**
- [ ] Handles multi-gigabyte input
- [ ] Memory usage stays within bounds
- [ ] Throughput matches expectations
- [ ] Incremental flushing triggers correctly

### 4. Quality Control
```bash
# Statistics comparison
samtools flagstat input.bam
bams3 stats output.bams3

# Read count verification
samtools view -c input.bam
bams3 query output.bams3 --count
```
**Verify:**
- [ ] Read counts match exactly
- [ ] Mapped/unmapped counts match
- [ ] No data loss in conversion

## Acceptance Criteria
- [ ] Successfully processes BWA output without errors
- [ ] All reads preserved (count verification)
- [ ] Output format correct and queryable
- [ ] Performance documented (reads/sec, memory usage)
- [ ] Instructions documented in README

## Blockers
- BWA not currently installed on test system
- Need reference genome and test data
