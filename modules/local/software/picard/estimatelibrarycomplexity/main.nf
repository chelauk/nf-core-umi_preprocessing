// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:] 
def options    = initOptions(params.options)

process PICARD_ESTIMATELIBRARYCOMPLEXITY {
    tag "$meta.id"
    label 'process_medium'
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
    tuple val(meta), path("*.library_complexity.txt"), emit: metrics
    path "*.version.txt"          , emit: version
    
    script:
    def software = getSoftwareName(task.process)
    def avail_mem = 3
    if (!task.memory) {
        log.info '[picard EstimateLibaryComplexity] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def prefix   = params.stage == "two" ? "${meta.id}" : "${meta.id}_${meta.run}"
    def max_records = task.memory.toGiga() * 100000
    """
    picard -Xmx${avail_mem}g \\
    EstimateLibraryComplexity \\
    --MAX_RECORDS_IN_RAM ${max_records} \\
    --INPUT  $bam \\
    --BARCODE_TAG RX \\
    --OUTPUT ${prefix}.library_complexity.txt 
    echo \$(picard EstimateLibraryComplexity --version 2>&1) | sed 's/^.*Version://' > ${software}.version.txt
    """
    
    stub:
    def software = getSoftwareName(task.process)
    def prefix   = params.stage == "two" ? "${meta.id}" : "${meta.id}_${meta.run}"
    """
    touch ${prefix}.library_complexity.txt 
    echo 2.6.2 > ${software}.version.txt
    """
}
