// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
params.options = [:]
def options    = initOptions(params.options)

process MARK_DUPLICATES {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::gatk4=4.1.9" : null)
    //container "quay.io/biocontainers/gatk4:4.1.9--0"

    input:
    tuple val(meta), file(bam)

    output:
    tuple val(meta), path("*bam"), emit: md_bam
    tuple val(meta), path("*bai"), emit: md_bai
    val meta,                          emit: tsv
    tuple val(meta), path("*.metrics"),     emit: report

    script:
    markdup_java_options = task.memory.toGiga() < 8 ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
    metrics = "-M ${meta.patient}_${meta.sample}.bam.metrics"

    if (params.no_gatk_spark)
    """
    gatk --java-options \\
        MarkDuplicates \\
        --MAX_RECORDS_IN_RAM 50000 \\
        --INPUT $bam \\
        --METRICS_FILE ${meta.patient}_${meta.sample}.bam.metrics \\
        --TMP_DIR . \\
        --ASSUME_SORT_ORDER coordinate \\
        --CREATE_INDEX true \\
        --OUTPUT ${meta.patient}_${meta.sample}.md.bam
        mv ${meta.patient}_${meta.sample}.md.bai ${meta.patient}_${meta.sample}.md.bam.bai
    """
    else
    """
    gatk --java-options ${markdup_java_options} \\
        MarkDuplicatesSpark \\
        -I $bam \\
        -O ${meta.patient}_${meta.sample}.md.bam \\
        ${metrics} \\
        --tmp-dir . \\
        --create-output-bam-index true \\
        --spark-master local[${task.cpus}]
    """
    stub:
    """
    touch ${meta.patient}_${meta.sample}.bam
    touch ${meta.patient}_${meta.sample}.bai
    touch ${meta.patient}_${meta.sample}.bam.metrics
    """

}