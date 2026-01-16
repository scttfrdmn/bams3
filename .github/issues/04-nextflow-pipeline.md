# Create Nextflow reference pipeline

**Labels:** pipeline, documentation, priority-high
**Priority:** High

## Goal
Demonstrate complete WGS workflow using cloud-optimized BAMS3 format in production-grade pipeline.

## Pipeline Stages

### 1. Alignment
```nextflow
process ALIGN {
    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2)
    path(reference)

    output:
    path("${sample_id}.bams3")

    script:
    """
    bwa mem -t ${task.cpus} ${reference} ${reads_r1} ${reads_r2} | \
      bams3 convert --stdin ${sample_id}.bams3 \
        --workers ${task.cpus} \
        --sort-buffer ${task.memory.toGiga()}G
    """
}
```

### 2. Quality Control
```nextflow
process QC {
    input:
    path(bams3)

    output:
    path("${bams3.baseName}_stats.json")

    script:
    """
    bams3 stats ${bams3} --json > ${bams3.baseName}_stats.json
    """
}
```

### 3. Parallel Variant Calling
```nextflow
process VARIANT_CALL {
    input:
    path(bams3)
    val(region)
    path(reference)

    output:
    path("${region}.vcf")

    script:
    """
    bams3 to-bam ${bams3} - --region ${region} | \
      gatk HaplotypeCaller \
        -I /dev/stdin \
        -R ${reference} \
        -O ${region}.vcf \
        -L ${region}
    """
}
```

### 4. Merge VCFs
```nextflow
process MERGE_VCFS {
    input:
    path(vcfs)

    output:
    path("merged.vcf.gz")

    script:
    """
    bcftools concat ${vcfs} -Oz -o merged.vcf.gz
    bcftools index merged.vcf.gz
    """
}
```

## Features to Demonstrate

### Cloud-Native Architecture
- [ ] Direct S3 read/write (no local storage)
- [ ] Parallel processing (all chromosomes simultaneously)
- [ ] Minimal data transfer (selective region access)
- [ ] No intermediate BAM files

### Performance Optimizations
- [ ] Resource-adaptive (auto-configure based on instance)
- [ ] Parallel compression during alignment
- [ ] Chunked processing for variant calling

### Cost Optimization
- [ ] Same-region data transfer (free)
- [ ] Spot instance support
- [ ] Efficient resource allocation

## Deliverables
- [ ] `workflows/wgs.nf` - Main workflow
- [ ] `nextflow.config` - Configuration
- [ ] `README.md` - Usage instructions
- [ ] `docs/cloud-costs.md` - Cost analysis
- [ ] `docs/performance.md` - Performance benchmarks

## Acceptance Criteria
- [ ] Pipeline runs successfully on AWS Batch
- [ ] Generates correct VCF output
- [ ] Demonstrates 75%+ cost reduction vs traditional
- [ ] Documentation includes cost/performance comparisons
- [ ] Can process 30x WGS in <4 hours

## Related
- Nextflow documentation: https://www.nextflow.io
- nf-core best practices: https://nf-co.re/docs

## Notes
Key demonstration of real-world value proposition.
