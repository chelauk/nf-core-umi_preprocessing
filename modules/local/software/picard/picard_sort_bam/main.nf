/*
 * A rule of thumb for reads of ~100bp is to set MAX_RECORDS_IN_RAM to be 250,000 
 * reads per each GB given to the -Xmx parameter for SortSam. 
 * Thanks to Keiran Raine for performing the experiments to arrive at these numbers.
 * I am going for 100,000
 */
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../../../../nf-core/software/functions'
params.options = [:]
def options    = initOptions(params.options)

process PICARD_SORT_BAM {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::picard=2.23.8" : null)
    //container "quay.io/biocontainers/picard:2.23.8--0"

    input:
    tuple val(meta), file(bam)
    path fasta
    path dict

    output:
    tuple val(meta), file("*sort*.bam")

    script:
    def picard_opts = params.second_file ? "mv ${meta.id}_sort.bam ${meta.id}_sort_2.bam" : ""
    def max_records = task.memory.toGiga() * 100000
    """
    mkdir tmpdir
    picard -Xmx${task.memory.toGiga()}g SortSam \\
    MAX_RECORDS_IN_RAM=${max_records} \\
    SORT_ORDER=queryname \\
    TMP_DIR=./tmpdir \\
    INPUT=${bam[0]} \\
    OUTPUT=${meta.id}_sort.bam \\
    VALIDATION_STRINGENCY=LENIENT
    ${picard_opts}
    """
    stub:
    picard_opts = params.second_file ? "mv ${meta.id}_sort.bam ${meta.id}_sort_2.bam" : ""
    """
    touch ${meta.id}_sort.bam
    ${picard_opts}
    """
    }
