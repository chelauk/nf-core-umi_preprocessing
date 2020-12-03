params.second_sort = true

include { BAM_TO_FASTQ }                    from '../../software/picard/bam_to_fastq/main'
include { BWA_ALN }                         from '../../software/bwa/bwa_aln/main'             addParams(options: params.bwamem1_mem_options)
include { PICARD_SORT_BAM }                 from '../../software/picard/picard_sort_bam/main'
include { PICARD_SORT_BAM as SORT_BAM_TWO } from '../../software/picard/picard_sort_bam/main'   addParams(second_file: true )
include { PICARD_MERGE_BAMS }               from '../../software/picard/picard_merge_bams/main' 
include { MARK_DUPLICATES }                 from '../../software/gatk/mark_duplicates/main'

workflow UMI_STAGE_TWO {

    take:
    filtered_bam
    bwa_index
    fasta
    fasta_fai
    dict

    main:
    BAM_TO_FASTQ(filtered_bam)
    BWA_ALN(BAM_TO_FASTQ.out,bwa_index,fasta,fasta_fai)
    PICARD_SORT_BAM(BWA_ALN.out.bam,fasta,dict)
    SORT_BAM_TWO(filtered_bam,fasta,dict)
    PICARD_MERGE_BAMS(PICARD_SORT_BAM.out.join(SORT_BAM_TWO.out),fasta,dict)
    MARK_DUPLICATES(PICARD_MERGE_BAMS.out.merged_bam)    
}