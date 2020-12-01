// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../../nf-core/software/functions'

params.options = [:]
def options    = initOptions(params.options)

process MERGE_BAM_ALIGNMENT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::picard=2.23.8" : null)
    container "quay.io/biocontainers/picard:2.23.8--0"

    input:
    tuple val(meta), file(aligned_unmarked_sam), file(unaligned_marked_bam)
    path fasta
    path dict

    output:
    tuple val(meta), file("*aligned_marked.bam")

    script:
    """
    picard -Xmx${task.memory.toGiga()}g MergeBamAlignment \\
    MAX_RECORDS_IN_RAM=4000000 \\
    R=$fasta \\
    UNMAPPED_BAM=$unaligned_marked_bam \\
    ALIGNED_BAM=$aligned_unmarked_sam \\
    O="${meta.id}.aligned_marked.bam" \\
    ADD_MATE_CIGAR=true CLIP_ADAPTERS=false \\
    CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true \\
    MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant \\
    ATTRIBUTES_TO_RETAIN=XS
    """ 
    }