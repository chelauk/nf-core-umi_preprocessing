// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../../../nf-core/software/functions'
params.options = [:]
def options    = initOptions(params.options)

process MARK_DUPLICATES {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::gatk4=4.1.9" : null)
    container "quay.io/biocontainers/gatk4:4.1.9--0"

    input:
    tuple val(meta), file(bam)

    output:
    tuple meta, path("*md.ba{m,i}"), emit: bam
    // set idPatient, idSample into tsv_bam_duplicates_marked
    tuple meta, path("*.metrics"), emit: mark_dup_report

    script:
    markdup_java_options = task.memory.toGiga() < 8 ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
    metrics = "-M ${meta.patient}.${meta.sample}.bam.metrics"
    if (params.no_gatk_spark)
    """
    gatk --java-options ${meta.patient}.${meta.sample}.bam.metrics" \\
        MarkDuplicates \\
        --MAX_RECORDS_IN_RAM 50000 \\
        --INPUT $bam \\
        --METRICS_FILE  \\
        --TMP_DIR . \\
        --ASSUME_SORT_ORDER coordinate \\
        --CREATE_INDEX true \\
        --OUTPUT ${meta.patient}.${meta.sample}.md.bam
    """
    else
    """
    gatk --java-options ${markdup_java_options} \\
        MarkDuplicatesSpark \\
        -I $bam \\
        -O ${meta.patient}.${meta.sample}.md.bam \\
        ${metrics} \\
        --tmp-dir . \\
        --create-output-bam-index true \\
        --spark-master local[${task.cpus}]
    """
}