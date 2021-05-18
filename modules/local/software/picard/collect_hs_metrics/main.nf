// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PICARD_COLLECT_HS_METRICS {
    tag "${meta.id}"

    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::picard=2.23.8" : null)
    //container "quay.io/biocontainers/picard:2.23.8--0"

    input:
    tuple val(meta), file(bam)
    //path interval_list
    path fasta
    path fai
    path dict
    path interval_list

    output:
    tuple val(meta), file("*metrics.txt" ), emit: hs_metrics
    path  "*.version.txt"                 , emit: version

    script:
    def software  = getSoftwareName(task.process)
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def avail_mem = 3
        if (!task.memory) {
        log.info '[Picard CollectHsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard -Xmx${avail_mem}g CollectHsMetrics \\
    R=${fasta} \\
    ${options.args} \\
    I=${bam} \\
    O=${meta.patient}_${meta.id}_hs_metrics_2.txt \\
    BAIT_INTERVALS=${interval_list} \\
    TARGET_INTERVALS=${interval_list} 
    echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d: > ${software}.version.txt
    """
    stub:
    """
    touch ${meta.id}_hs_metrics.txt
    touch software.version.txt
    """
    }