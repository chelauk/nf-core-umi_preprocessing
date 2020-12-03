// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../../../../nf-core/software/functions'

params.options = [:]
def options    = initOptions(params.options)

process BAM_TO_FASTQ {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::picard=2.23.8" : null)
    container "quay.io/biocontainers/picard:2.23.8--0"

    input:
    tuple val(meta), file(bam)

    output:
    tuple val(meta), file("*fastq"), emit: fastq

    script:
    """
    picard SamToFastq \\
    MAX_RECORDS_IN_RAM=4000000 \\
    INPUT=$bam \\
    FASTQ="${meta.id}.fastq" \\
    CLIPPING_ATTRIBUTE=XT \\
    CLIPPING_ACTION=2 \\
    INTERLEAVE=true \\
    NON_PF=true 
    """
    }