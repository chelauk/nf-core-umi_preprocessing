// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
params.options = [:]
def options    = initOptions(params.options)

process MARK_ILLUMINA_ADAPTERS {
    tag "{$meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::picard=2.23.8" : null)
    //container "quay.io/biocontainers/picard:2.23.8--0"

    input:
    tuple val(meta), file(bam)

    output:
    tuple val(meta), file("*bam"), emit : bam
    tuple val(meta), file("*_mark_adapter.metrics"), emit: mark_adaptor_log

    script:
    def max_records = task.memory.toGiga() * 100000
    def prefix   = params.stage == "two" ? "${meta.id}" : "${meta.id}_${meta.run}"
    """
    picard -Xmx${task.memory.toGiga()}g  MarkIlluminaAdapters \\
    MAX_RECORDS_IN_RAM=${max_records} \\
    INPUT=$bam \\
    OUTPUT="${prefix}_unaln_umi_marked.bam" \\
    METRICS="${prefix}_mark_adapter.metrics"
    """
    stub:
    def prefix   = params.stage == "two" ? "${meta.id}" : "${meta.id}_${meta.run}"
    """
    touch ${prefix}_unaln_umi_marked.bam
    touch ${prefix}_mark_adapter.metrics
	"""
    }

