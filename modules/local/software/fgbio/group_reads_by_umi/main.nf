// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GROUP_READS_BY_UMI {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::fgbio=1.3.0" : null)
    container "quay.io/biocontainers/fgbio:1.3.0--0"

    input:
    tuple val(meta), file(bam)

    output:
    tuple val(meta), file("*umi_group.bam"), emit: bam
    tuple val(meta), file("*.metrics"),      emit: group_metrics

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g GroupReadsByUmi \\
    -i ${bam} \\
    -f ${meta.patient}_${meta.sample}_aln_merged_umi.metrics \\
    -s adjacency -m 30 -t RX -T MI --min-umi-length 9 \\
    -o ${meta.patient}_${meta.sample}_umi_group.bam 
   """
}