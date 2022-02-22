// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:] 
def options    = initOptions(params.options)

process PICARD_UMI_AWARE_MARKDUPLICATES_WITH_MATE_CIGAR {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::picard=2.26.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/picard:2.26.2--hdfd78af_0"
    } else {
        container "quay.io/picard:2.26.2--hdfd78af_0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*bam"),      emit: md_bam
    tuple val(meta), path("*bai"),      emit: md_bai
    val (meta),                         emit: tsv
    tuple val(meta), path("*.metrics"), emit: report 

    script:
    def software = getSoftwareName(task.process)
    def avail_mem = 3
    if (!task.memory) {
        log.info '[picard UmiAwareMarkDuplicatesWithMateCigar] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    picard -Xmx${avail_mem}g \\
    UmiAwareMarkDuplicatesWithMateCigar \\
    --INPUT  $bam \\
    --OUTPUT ${prefix}_umi_aware_md.bam \\
    --ASSUME_SORT_ORDER  coordinate \\
    --METRICS_FILE ${prefix}_duplicate.metrics \\
    --UMI_METRICS_FILE ${prefix}_umi.metrics \\
    --CREATE_INDEX true
    echo \$(picard UmiAwareMarkDuplicatesWithMateCigar --version 2>&1) | sed 's/^.*Version://' > ${software}.version.txt
    """
    
    stub:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    touch ${prefix}_umi_aware_md.bam
    touch ${prefix}_umi_aware_md.bam.bai
    touch ${prefix}_duplicate.metrics
    touch ${prefix}_umi.metrics

    echo 2.6.2 > ${software}.version.txt
    """
}
