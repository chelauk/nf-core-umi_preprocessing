// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process COLLECT_DUPLEX_SEQ_METRICS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::fgbio=1.3.0" : null)
    //container "quay.io/biocontainers/fgbio:1.3.0--0"

    input:
    tuple val(meta), path(bam)
    path interval_list

    output:
    tuple val(meta), file("*.txt"), emit: metrics

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g CollectDuplexSeqMetrics \\
    --input ${bam} \\
    --intervals ${interval_list} \\
    --duplex-umi-counts true \\
    --output ${meta.id}
    """
    stub:
    """
    touch ${meta.id}.txt
    """

}