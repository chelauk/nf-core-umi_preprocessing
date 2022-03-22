/*
===================================
          UMI stage one
===================================
*/

include { BED_TO_INTERVAL_LIST }       from '../../software/picard/bed_to_interval_list/main'        addParams(options: params.bed_to_intervals_options)
include { PICARD_COLLECT_HS_METRICS }  from '../../software/picard/collect_hs_metrics/main'     addParams(options: params.collect_hs_metrics_options)
include { ERRORRATE_BY_READ_POSITION } from '../../software/fgbio/error_rate/main'              addParams(options: params.error_rate_options)
include { GROUP_READS_BY_UMI }         from '../../software/fgbio/group_reads_by_umi/main'      addParams(options: params.group_reads_mapping_options)
include { FGBIO_SORT_BAM }             from '../../software/fgbio/fgbio_sort_bam/main'          addParams(options: params.fgbio_sort_mapping_options)
include { CALL_CONSENSUS }             from '../../software/fgbio/call_consensus/main'          addParams(options: params.fgbio_call_consensus_mapping_options)
include { FILTER_CONSENSUS }           from '../../software/fgbio/filter_consensus/main'        addParams(options: params.fgbio_filter_mapping_options)

workflow POST_MERGE {
    take:
    input_samples
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
    PICARD_COLLECT_HS_METRICS(input_samples, fasta, fasta_fai, dict, BED_TO_INTERVAL_LIST.out.interval_list)
    ERRORRATE_BY_READ_POSITION(input_samples,fasta,dict,dbsnp,dbsnp_index,BED_TO_INTERVAL_LIST.out.interval_list)
    GROUP_READS_BY_UMI(input_samples)
    FGBIO_SORT_BAM(GROUP_READS_BY_UMI.out.bam)
    CALL_CONSENSUS(FGBIO_SORT_BAM.out)
    FILTER_CONSENSUS(CALL_CONSENSUS.out,fasta, min_reads)

    emit:
    filtered_bam        = FILTER_CONSENSUS.out
    hs_metrics          = PICARD_COLLECT_HS_METRICS.out.hs_metrics
    error_rate          = ERRORRATE_BY_READ_POSITION.out.error_rate
    group_metrics       = GROUP_READS_BY_UMI.out.group_metrics
    iv_list             = BED_TO_INTERVAL_LIST.out.interval_list
    
}
