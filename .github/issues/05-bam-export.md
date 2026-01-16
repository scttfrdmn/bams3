# Implement BAM export for tool integration

**Labels:** enhancement, integration, priority-high
**Priority:** High

## Goal
Enable integration with existing bioinformatics tools (GATK, samtools, etc.) by exporting BAMS3 to BAM format.

## Command Interface
```bash
# Export entire file
bams3 to-bam sample.bams3 output.bam

# Export specific region (selective chunk download)
bams3 to-bam sample.bams3 output.bam --region chr1:1000000-2000000

# Stream to stdout for piping
bams3 to-bam sample.bams3 - --region chr1 | samtools view
bams3 to-bam s3://bucket/sample.bams3 - --region chr22 | gatk HaplotypeCaller -I /dev/stdin
```

## Core Functionality

**Path:** `pkg/bams3/bam_export.go` (new file)

```go
type BAMExporter struct {
    reader *Reader
    writer *bam.Writer
}

func (e *BAMExporter) Export(region *Region) error {
    // 1. Read header from BAMS3
    // 2. Write BAM header
    // 3. Query chunks for region (or all if region == nil)
    // 4. Convert reads to BAM format
    // 5. Write BAM records
    // 6. Close writer (finalizes BAM)
}
```

## Features
- [ ] Full file export
- [ ] Region-based export (selective chunk download)
- [ ] Streaming to stdout
- [ ] Progress reporting for large exports
- [ ] Coordinate sorting (already sorted in BAMS3)
- [ ] Index generation (optional)

## BAM Format Requirements
- [ ] Correct header (@HD, @SQ, @RG, @PG)
- [ ] BGZF compression
- [ ] Proper record encoding
- [ ] Compatible with samtools/GATK/Picard

## CLI Integration

**New subcommand:** `bams3 to-bam`

**Flags:**
```go
--region string     // Genomic region to export (chr:start-end)
--index            // Generate .bai index file (default: true)
--compression int  // BGZF compression level (default: 6)
--threads int      // Parallel decompression (default: auto)
```

## Performance Considerations
- [ ] Only download needed chunks (region query)
- [ ] Parallel decompression of chunks
- [ ] Streaming output (no temp files)
- [ ] BGZF multi-threaded compression

## Testing

### Compatibility Tests
```bash
# Verify samtools compatibility
bams3 to-bam sample.bams3 test.bam
samtools view test.bam | head
samtools flagstat test.bam
samtools index test.bam

# Verify GATK compatibility
bams3 to-bam sample.bams3 - --region chr22 | \
  gatk ValidateSamFile -I /dev/stdin

# Verify Picard compatibility
java -jar picard.jar ValidateSamFile I=test.bam
```

### Correctness Tests
```bash
# Round-trip test
samtools view -h original.bam | bams3 convert --stdin test.bams3
bams3 to-bam test.bams3 roundtrip.bam
diff <(samtools view original.bam) <(samtools view roundtrip.bam)
```

## Acceptance Criteria
- [ ] Produces valid BAM files (passes samtools/GATK validation)
- [ ] Bit-for-bit identical to original BAM (round-trip test)
- [ ] Region queries only download necessary chunks
- [ ] Streaming output works for large files
- [ ] Performance within 2x of direct BAM access
- [ ] Documentation with integration examples

## Dependencies
- `github.com/biogo/hts` - BAM format handling
- Or implement BAM writer directly for better control

## Notes
Required for GATK integration and tool ecosystem compatibility.
