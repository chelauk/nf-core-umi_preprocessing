// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../../../../nf-core/software/functions'
params.options = [:]
def options    = initOptions(params.options)

process PICARD_SORT_BAM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::picard=2.23.8" : null)
    container "quay.io/biocontainers/picard:2.23.8--0"

    input:
    tuple val(meta), file(bam)
    path fasta
    path dict

    output:
    tuple val(meta), file("*sort*.bam")

    script:
    picard_opts = params.second_file ? "mv ${meta.id}_sort.bam ${meta.id}_sort_2.bam" : ""
    """
    picard -Xmx${task.memory.toGiga()}g SortSam \\
    MAX_RECORDS_IN_RAM=4000000 \\
    SORT_ORDER=queryname \\
    INPUT=$bam \\
    OUTPUT=${meta.id}_sort.bam \\
    VALIDATION_STRINGENCY=LENIENT
    ${picard_opts}
    """ 
    }