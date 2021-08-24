// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../../../../nf-core/software/functions'

params.options = [:]
def options    = initOptions(params.options)

process FASTQ_TO_BAM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::fgbio=1.3.0" : null)
//    container "quay.io/biocontainers/fgbio:1.3.0--0"

    input:
    tuple val(meta), path(reads)
    val library
    val rstructure 

    output:
    tuple val(meta), file("*bam"), emit: bam

    script:
    """  
    mkdir temp
    fgbio -Xmx${task.memory.toGiga()}g --tmp-dir=./temp FastqToBam \\
    --input $reads \\
    --output ${meta.id}_unaln.bam \\
    --read-structures $rstructure \\
    --umi-tag RX \\
    --sample ${meta.id} \\
    --library $library
    """
    stub:
    """
    touch ${meta.id}_unaln.bam
    """
}
