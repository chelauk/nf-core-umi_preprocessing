// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
params.options = [:]
def options    = initOptions(params.options)

process MARK_ILLUMINA_ADAPTERS {
    tag "{$meta.id}"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::picard=2.23.8" : null)
    //container "quay.io/biocontainers/picard:2.23.8--0"

    input:
    tuple val(meta), file(bam)

    output:
    tuple val(meta), file("*bam"), emit : bam
    //path "*_mark_adapter.metrics", emit: mark_adaptor_log

    script:
    """
    picard -Xmx${task.memory.toGiga()}g  MarkIlluminaAdapters \\
    MAX_RECORDS_IN_RAM=4000000 \\
    INPUT=$bam \\
    OUTPUT="${meta.id}_unaln_umi_marked.bam" \\
    M="${meta.patient}_${meta.id}_mark_adapter.metrics"
    """
    }