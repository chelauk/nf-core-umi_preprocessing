// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:] 
def options    = initOptions(params.options)

process PICARD_UMI_AWARE_MARKDUPLICATES_WITH_MATE_CIGAR {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::picard=2.26.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/picard:2.26.2--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/YOUR-TOOL-HERE"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.txt"), emit: umi_aware_metrics
    path "*.version.txt"          , emit: version
    
    script:
    def software = getSoftwareName(task.process)
    def avail_mem = 3
    if (!task.memory) {
        log.info '[picard EstimateLibraryComplexity] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    picard -Xmx${avail_mem}g \\
    UmiAwareMarkDuplicatesWithMateCigar \\
    --INPUT  $bam \\
    --OUTPUT ${prefix}_umi_aware_md.bam \\
    --METRICS ${prefix}_duplicate_metrics.txt \\
    --UMI_METRICS ${prefix}_umi_metrics.txt
    echo \$(picard UmiAwareMarkDuplicatesWithMateCigar --version 2>&1) | sed 's/^.*Version://' > ${software}.version.txt
    """
    
    stub:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    touch ${prefix}_umi_aware_md.bam
    touch ${prefix}_duplicate_metrics.txt
    touch ${prefix}_umi_metrics.txt

    echo 2.6.2 > ${software}.version.txt
    """
}
