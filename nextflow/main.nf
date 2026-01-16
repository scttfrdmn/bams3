#!/usr/bin/env nextflow

/*
 * BAMS3 Reference Pipeline
 *
 * Demonstrates cloud-native genomics workflow using BAMS3 format
 * for alignment storage with direct S3 integration and selective
 * region queries for downstream analysis.
 *
 * Usage:
 *   nextflow run main.nf -c nextflow.config -params-file params.yaml
 */

nextflow.enable.dsl=2

/*
 * Pipeline parameters with defaults
 */
params.samples = null              // CSV file with sample metadata
params.reference = null            // Reference genome FASTA
params.output_dir = "s3://genomics-data/results"
params.bams3_dir = "${params.output_dir}/bams3"
params.vcf_dir = "${params.output_dir}/vcf"
params.qc_dir = "${params.output_dir}/qc"

// Analysis parameters
params.regions = null              // BED file with regions for variant calling
params.variant_caller = "gatk"    // "gatk" or "bcftools"

// Resource parameters
params.align_cpus = 16
params.align_memory = "32 GB"
params.sort_buffer = "28G"
params.vcall_cpus = 4
params.vcall_memory = "16 GB"

log.info """
╔═══════════════════════════════════════════════════════════╗
║         BAMS3 CLOUD-NATIVE GENOMICS PIPELINE              ║
╚═══════════════════════════════════════════════════════════╝

Reference:     ${params.reference}
Samples:       ${params.samples}
Output:        ${params.output_dir}
Regions:       ${params.regions ?: 'Whole genome'}
Variant caller: ${params.variant_caller}

Resources:
  Alignment:   ${params.align_cpus} CPUs, ${params.align_memory}
  Calling:     ${params.vcall_cpus} CPUs, ${params.vcall_memory}

═══════════════════════════════════════════════════════════
"""

/*
 * Process 1: Index reference genome (if needed)
 */
process INDEX_REFERENCE {
    tag "${reference.name}"

    input:
    path reference

    output:
    tuple path(reference), path("${reference}.*"), emit: indexed

    script:
    """
    # Check if already indexed
    if [ ! -f ${reference}.bwt ]; then
        echo "Indexing reference genome..."
        bwa index ${reference}
    else
        echo "Reference already indexed"
        touch ${reference}.bwt
    fi
    """
}

/*
 * Process 2: Align reads and convert to BAMS3 (zero-copy pipeline)
 *
 * This is the core BAMS3 advantage: BWA streams directly to BAMS3
 * with coordinate sorting, no intermediate SAM/BAM files.
 */
process ALIGN_TO_BAMS3 {
    tag "${sample_id}"
    publishDir "${params.bams3_dir}", mode: 'copy', enabled: params.bams3_dir.startsWith('s3://')

    cpus params.align_cpus
    memory params.align_memory

    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2)
    tuple path(reference), path(ref_indices)

    output:
    tuple val(sample_id), path("${sample_id}.bams3"), emit: bams3
    tuple val(sample_id), path("${sample_id}.align_stats.txt"), emit: stats

    script:
    def bams3_output = params.bams3_dir.startsWith('s3://') ?
        "${params.bams3_dir}/${sample_id}.bams3" :
        "${sample_id}.bams3"

    """
    echo "Starting alignment: ${sample_id}"
    echo "Output: ${bams3_output}"

    # Zero-copy pipeline: BWA → BAMS3 (single step!)
    bwa mem -t ${task.cpus} \
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' \
        ${reference} \
        ${reads_r1} ${reads_r2} \
        2> bwa.log | \
    bams3 convert --stdin "${bams3_output}" \
        --workers ${task.cpus} \
        --sort-buffer ${params.sort_buffer}

    # Extract alignment statistics
    if [[ "${bams3_output}" == s3://* ]]; then
        # Download metadata for stats
        aws s3 cp "${bams3_output}/_metadata.json" metadata.json
        jq -r '.statistics | "Total reads: \\(.total_reads)\\nMapped reads: \\(.mapped_reads)\\nUnmapped reads: \\(.unmapped_reads)\\nMean coverage: \\(.mean_coverage)"' \
            metadata.json > ${sample_id}.align_stats.txt
    else
        jq -r '.statistics | "Total reads: \\(.total_reads)\\nMapped reads: \\(.mapped_reads)\\nUnmapped reads: \\(.unmapped_reads)\\nMean coverage: \\(.mean_coverage)"' \
            ${sample_id}.bams3/_metadata.json > ${sample_id}.align_stats.txt
    fi

    echo "Alignment complete for ${sample_id}"
    """
}

/*
 * Process 3: Quality control statistics
 */
process QC_STATS {
    tag "${sample_id}"
    publishDir "${params.qc_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bams3_dir)

    output:
    tuple val(sample_id), path("${sample_id}.qc_report.txt"), emit: report

    script:
    def bams3_path = params.bams3_dir.startsWith('s3://') ?
        "${params.bams3_dir}/${sample_id}.bams3" :
        bams3_dir

    """
    echo "Generating QC report for ${sample_id}"

    # Get comprehensive statistics from BAMS3
    bams3 stats "${bams3_path}" > ${sample_id}.qc_report.txt

    echo "QC complete for ${sample_id}"
    """
}

/*
 * Process 4: Variant calling on specific regions
 *
 * Demonstrates selective region extraction from BAMS3 - only downloads
 * chunks overlapping the target region instead of the entire file.
 */
process CALL_VARIANTS_REGION {
    tag "${sample_id}-${region_name}"
    publishDir "${params.vcf_dir}/${sample_id}", mode: 'copy'

    cpus params.vcall_cpus
    memory params.vcall_memory

    input:
    tuple val(sample_id), path(bams3_dir), val(region_name), val(region)
    tuple path(reference), path(ref_indices)

    output:
    tuple val(sample_id), val(region_name), path("${sample_id}.${region_name}.vcf.gz"), emit: vcf

    script:
    def bams3_path = params.bams3_dir.startsWith('s3://') ?
        "${params.bams3_dir}/${sample_id}.bams3" :
        bams3_dir

    if (params.variant_caller == "gatk")
        """
        echo "Calling variants on ${region} for ${sample_id}"

        # Extract region from BAMS3 to BAM (streams only relevant chunks!)
        bams3 to-bam "${bams3_path}" - --region ${region} | \
        gatk HaplotypeCaller \
            -I /dev/stdin \
            -R ${reference} \
            -O ${sample_id}.${region_name}.vcf.gz \
            -L ${region} \
            --native-pair-hmm-threads ${task.cpus}

        echo "Variant calling complete for ${sample_id} ${region_name}"
        """
    else if (params.variant_caller == "bcftools")
        """
        echo "Calling variants on ${region} for ${sample_id}"

        # Extract region from BAMS3 to BAM
        bams3 to-bam "${bams3_path}" ${sample_id}.${region_name}.bam --region ${region}
        samtools index ${sample_id}.${region_name}.bam

        # Call variants with bcftools
        bcftools mpileup -Ou -f ${reference} \
            --regions ${region} \
            ${sample_id}.${region_name}.bam | \
        bcftools call -mv -Oz -o ${sample_id}.${region_name}.vcf.gz

        # Cleanup temporary BAM
        rm ${sample_id}.${region_name}.bam ${sample_id}.${region_name}.bam.bai

        echo "Variant calling complete for ${sample_id} ${region_name}"
        """
}

/*
 * Process 5: Merge per-region VCFs into whole-genome VCF
 */
process MERGE_VCFS {
    tag "${sample_id}"
    publishDir "${params.vcf_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: merged_vcf

    script:
    """
    echo "Merging VCFs for ${sample_id}"

    # Index all VCFs
    for vcf in ${vcfs}; do
        tabix -p vcf \$vcf
    done

    # Merge all region VCFs
    bcftools concat -Oz -o ${sample_id}.vcf.gz ${vcfs}
    tabix -p vcf ${sample_id}.vcf.gz

    echo "VCF merging complete for ${sample_id}"
    """
}

/*
 * Process 6: Generate MultiQC report
 */
process MULTIQC {
    publishDir "${params.qc_dir}", mode: 'copy'

    input:
    path(qc_files)

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc . -n multiqc_report.html
    """
}

/*
 * Workflow definition
 */
workflow {
    // Validate required parameters
    if (!params.samples || !params.reference) {
        error "Missing required parameters. Use: --samples <csv> --reference <fasta>"
    }

    // Parse sample sheet
    Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample_id, file(row.reads_r1), file(row.reads_r2))
        }
        .set { samples_ch }

    // Reference genome
    reference_ch = Channel.fromPath(params.reference)

    // Index reference genome
    INDEX_REFERENCE(reference_ch)

    // Align reads to BAMS3 (zero-copy!)
    ALIGN_TO_BAMS3(samples_ch, INDEX_REFERENCE.out.indexed)

    // Generate QC statistics
    QC_STATS(ALIGN_TO_BAMS3.out.bams3)

    // Variant calling workflow
    if (params.regions) {
        // Parse regions BED file
        regions_ch = Channel
            .fromPath(params.regions)
            .splitCsv(sep: '\t')
            .map { row ->
                def region_name = "${row[0]}_${row[1]}_${row[2]}"
                def region = "${row[0]}:${row[1]}-${row[2]}"
                tuple(region_name, region)
            }

        // Create combinations of samples and regions
        samples_regions_ch = ALIGN_TO_BAMS3.out.bams3
            .combine(regions_ch)
            .map { sample_id, bams3, region_name, region ->
                tuple(sample_id, bams3, region_name, region)
            }

        // Call variants per region
        CALL_VARIANTS_REGION(
            samples_regions_ch,
            INDEX_REFERENCE.out.indexed
        )

        // Group VCFs by sample and merge
        CALL_VARIANTS_REGION.out.vcf
            .groupTuple(by: 0)
            .map { sample_id, region_names, vcfs ->
                tuple(sample_id, vcfs)
            }
            .set { vcfs_by_sample_ch }

        MERGE_VCFS(vcfs_by_sample_ch)
    }

    // Generate MultiQC report
    MULTIQC(QC_STATS.out.report.collect())
}

workflow.onComplete {
    log.info """
    ═══════════════════════════════════════════════════════════
    Pipeline completed!
    Status:     ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration:   ${workflow.duration}

    Results:
      BAMS3:    ${params.bams3_dir}
      VCFs:     ${params.vcf_dir}
      QC:       ${params.qc_dir}
    ═══════════════════════════════════════════════════════════
    """
}

workflow.onError {
    log.error """
    ═══════════════════════════════════════════════════════════
    Pipeline failed!
    Error:      ${workflow.errorMessage}
    ═══════════════════════════════════════════════════════════
    """
}
