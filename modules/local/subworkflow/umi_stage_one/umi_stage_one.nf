/*
===================================
          UMI stage one
===================================
*/

include { FASTQ_TO_BAM }           from '../../software/fgbio/fastq_to_bam/main' 
include { MARK_ILLUMINA_ADAPTERS } from '../../software/picard/mark_illumina_adapters/main' 
include { BAM_TO_FASTQ }           from '../../software/picard/bam_to_fastq/main' 
include { BWA_ALN }                from '../../software/bwa/bwa_aln/main'               addParams(options: params.bwamem1_mem_options)
include { PICARD_MERGE_BAMS }      from '../../software/picard/picard_merge_bams/main' 
include { MERGE_RUNS}              from '../../subworkflow/merge_runs/main'
include { GROUP_READS_BY_UMI }     from '../../software/fgbio/group_reads_by_umi/main' 
include { FGBIO_SORT_BAM }         from '../../software/fgbio/fgbio_sort_bam/main'
include { CALL_CONSENSUS }         from '../../software/fgbio/call_consensus/main'
include { FILTER_CONSENSUS }       from '../../software/fgbio/filter_consensus/main'

workflow UMI_STAGE_ONE {
    take:
    input_samples
    read_structure
    bwa_index
    fasta
    fasta_fai
    dict
    min_reads

    main:
    FASTQ_TO_BAM(input_samples, read_structure)
    MARK_ILLUMINA_ADAPTERS(FASTQ_TO_BAM.out.bam)
    BAM_TO_FASTQ(MARK_ILLUMINA_ADAPTERS.out.bam)
    BWA_ALN(BAM_TO_FASTQ.out.fastq,bwa_index,fasta,fasta_fai)
    PICARD_MERGE_BAMS(BWA_ALN.out.bam.join(MARK_ILLUMINA_ADAPTERS.out.bam),fasta,dict)
    MERGE_RUNS(PICARD_MERGE_BAMS.out.merged_bam)
    GROUP_READS_BY_UMI(MERGE_RUNS.out)
    FGBIO_SORT_BAM(GROUP_READS_BY_UMI.out.bam)
    CALL_CONSENSUS(FGBIO_SORT_BAM.out)
    FILTER_CONSENSUS(CALL_CONSENSUS.out,fasta, min_reads)

    emit:
    FILTER_CONSENSUS.out
}