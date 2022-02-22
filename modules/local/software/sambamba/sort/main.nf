// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../../../../nf-core/software/functions'
params.options = [:]
def options    = initOptions(params.options)
process SAMBAMBA_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sambamba=0.8.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sambamba:0.8.1--hadffe2f_1' :
        'quay.io/biocontainers/sambamba:0.8.1--hadffe2f_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*sorted.bam"), emit: bam
    path  "*.version.txt"               , emit: version

    script:
    def software = getSoftwareName(task.process)
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    sambamba sort  \
        $bam \
        --memory-limit=${task.memory.toGiga()}G \
        --tmpdir=./temp \
        --nthreads=${task.cpus}

    cat <<-END_VERSIONS > ${software}.version.txt
        sambamba: \$(echo \$(sambamba --version 2>&1) | sed 's/^.*sambamba //; s/ by.*//')
    END_VERSIONS
    """
    stub:
    def software = getSoftwareName(task.process)
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${bam}.sorted.bam
    touch ${software}.version.txt
    """
}