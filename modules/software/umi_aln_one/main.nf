// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../../nf-core/software/functions'

params.options = [:]
def options    = initOptions(params.options)

process UMI_ALN_ONE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::bwa=0.7.17" : null)
    container "quay.io/biocontainers/0.7.17--hed695b0_6"

    input:
    tuple val(meta), file(reads) 
    path bwa
    path fasta
    path fasta_fai

    output:
    tuple meta, file("*bam"), emit: bam

    script:
    CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    """
    bwa mem \\
    -M -p -k 50 \\
    -t ${task.cpus} \\
    ${fasta} ${reads} | \
    samtools sort  -n -@$task.cpus -o ${meta.id}.bam -
    echo \$(bwa version 2>&1) > bwa.version.txt
    """
    }