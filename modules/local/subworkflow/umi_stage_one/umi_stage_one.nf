/*
===================================
          UMI stage one
===================================
*/

include { BED_TO_INTERVAL_LIST }       from '../../software/picard/bed_to_interval_list/main'   addParams(options: params.bed_to_intervals_options)
include { FASTQ_TO_BAM }               from '../../software/fgbio/fastq_to_bam/main'            addParams(options: params.fastq_to_bam_options)
include { MARK_ILLUMINA_ADAPTERS }     from '../../software/picard/mark_illumina_adapters/main' addParams(options: params.mark_adapters_options)
include { BAM_TO_FASTQ }               from '../../software/picard/bam_to_fastq/main'           addParams(options: params.bam_to_fastq_options)
include { BWA_ALN }                    from '../../software/bwa/bwa_aln/main'                   addParams(options: params.bwamem1_mem_options)
include { PICARD_MERGE_BAMS }          from '../../software/picard/picard_merge_bams/main'      addParams(options: params.picard_merge_bams_options)
include { MERGE_RUNS}                  from '../../subworkflow/merge_runs/main'                 addParams(options: params.merge_runs_mapping_options)
include { PICARD_COLLECT_HS_METRICS }  from '../../software/picard/collect_hs_metrics/main'     addParams(options: params.collect_hs_metrics_options)
include { ERRORRATE_BY_READ_POSITION } from '../../software/fgbio/error_rate/main'              addParams(options: params.error_rate_options)
include { GROUP_READS_BY_UMI }         from '../../software/fgbio/group_reads_by_umi/main'      addParams(options: params.group_reads_by_umi_options)
include { FGBIO_SORT_BAM }             from '../../software/fgbio/fgbio_sort_bam/main'          addParams(options: params.fgbio_sort_mapping_options)
include { CALL_CONSENSUS }             from '../../software/fgbio/call_consensus/main'          addParams(options: params.fgbio_call_consensus_mapping_options)
include { FILTER_CONSENSUS }           from '../../software/fgbio/filter_consensus/main'        addParams(options: params.fgbio_filter_mapping_options)

workflow UMI_STAGE_ONE {
    take:
    input_samples
    library
    read_structure
    bwa_index
    fasta
    fasta_fai
    dict
    min_reads
    target_bed
    dbsnp
    dbsnp_index

    main:
    BED_TO_INTERVAL_LIST(target_bed, dict)
    FASTQ_TO_BAM(input_samples, library, read_structure)
    MARK_ILLUMINA_ADAPTERS(FASTQ_TO_BAM.out.bam)
    BAM_TO_FASTQ(MARK_ILLUMINA_ADAPTERS.out.bam)
    BWA_ALN(BAM_TO_FASTQ.out.fastq,bwa_index,fasta,fasta_fai)
    PICARD_MERGE_BAMS(BWA_ALN.out.bam.join(MARK_ILLUMINA_ADAPTERS.out.bam),fasta,dict)
    MERGE_RUNS(PICARD_MERGE_BAMS.out.merged_bam) // note this is the MERGE_RUNS subworkflow
    PICARD_COLLECT_HS_METRICS(MERGE_RUNS.out.bam_mapped, fasta, fasta_fai, dict, BED_TO_INTERVAL_LIST.out.interval_list)
    ERRORRATE_BY_READ_POSITION(MERGE_RUNS.out.bam_mapped,fasta,dict,dbsnp,dbsnp_index,BED_TO_INTERVAL_LIST.out.interval_list)
    GROUP_READS_BY_UMI(MERGE_RUNS.out)
    FGBIO_SORT_BAM(GROUP_READS_BY_UMI.out.bam)
    CALL_CONSENSUS(FGBIO_SORT_BAM.out)
    FILTER_CONSENSUS(CALL_CONSENSUS.out,fasta, min_reads)

    emit:
    filtered_bam  = FILTER_CONSENSUS.out
    hs_metrics_1    = PICARD_COLLECT_HS_METRICS.out.hs_metrics
    error_rate   = ERRORRATE_BY_READ_POSITION.out.error_rate
    group_metrics = GROUP_READS_BY_UMI.out.group_metrics
    iv_list       = BED_TO_INTERVAL_LIST.out.interval_list
    
}
