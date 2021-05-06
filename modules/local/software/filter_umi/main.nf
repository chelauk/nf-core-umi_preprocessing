// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process FILTER_UMIS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(${prefix}_1_val_1.fq.gz,${prefix}_2_val_2.fq.gz,${prefix}_3_val_3.fq.gz)    , emit: reads

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    [ ! -f  ${prefix}_1_val_1.fq.gz ] && ln -s ${reads[0]} ${prefix}_1_val_1.fq.gz
    [ ! -f  ${prefix}_3_val_2.fq.gz ] && ln -s ${reads[2]} ${prefix}_3_val_3.fq.gz
    filter_umis.py -u ${prefix}_2.fq.gz -v ${prefix}_1_val_1.fq.gz
    """

}