// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PICARD_MERGE_BAMS {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::picard=2.23.8" : null)
    //container "quay.io/biocontainers/picard:2.23.8--0"

    input:
    tuple val(meta), file(aligned_unmarked_sam), file(unaligned_marked_bam)
    path fasta
    path dict

    output:
    tuple val(meta), file("*merged.bam"), emit: merged_bam
    path  "*.version.txt"               , emit: version

    script:
    def max_records = task.memory.toGiga() * 100000
    def software  = getSoftwareName(task.process)
    def prefix   = params.stage == "two" ? "${meta.id}" : "${meta.id}_${meta.run}"
    def avail_mem = 3
        if (!task.memory) {
        log.info '[Picard MergeSamFiles] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard -Xmx${avail_mem}g MergeBamAlignment \\
    MAX_RECORDS_IN_RAM=${max_records} \\
    R=$fasta \\
    ${options.args} \\
    UNMAPPED_BAM=${unaligned_marked_bam} \\
    ALIGNED_BAM=${aligned_unmarked_sam} \\
    O="${prefix}.merged.bam" \\
    ADD_MATE_CIGAR=true CLIP_ADAPTERS=false \\
    CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true \\
    MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant \\
    ATTRIBUTES_TO_RETAIN=XS
    echo \$(picard MergeSamFiles --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d: > ${software}.version.txt
    """
    stub:
    def prefix   = params.stage == "two" ? "${meta.id}" : "${meta.id}_${meta.run}"
    """
    touch ${prefix}.merged.bam
    touch software.version.txt
    """
    }
