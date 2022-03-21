include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

environment = params.enable_conda ? "bioconda::samtools=1.10" : null
//container = "quay.io/biocontainers/samtools:1.10--h2e538c0_3"

process SAMTOOLS_MERGE_BAM {
    echo true
    label 'process_high'
    tag "${meta.id}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${meta.id}.bam"), emit: bam

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    samtools merge --threads ${task.cpus} temp.bam ${bam}
    samtools sort  --threads ${task.cpus} -o ${prefix}.bam  temp.bam
    """
    stub:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
