/*
===================================
        TRIMGALORE
===================================
*/

include { TRIMGALORE }                 from '../../../nf-core/software/trimgalore/main'   addParams(options: params.trimgalore_options)
include { FILTER_UMIS }                from '../../software/filter_umi/main.nf'           addParams(options: params.trim_umi_options)

workflow TRIMGALORE_WF {
    take:
    input_samples

    main:
    TRIMGALORE(input_samples)
    FILTER_UMIS(TRIMGALORE.out.reads)

    emit:
    trimmed_samples   = FILTER_UMIS.out.reads
    trim_qc           = TRIMGALORE.out.zip
}