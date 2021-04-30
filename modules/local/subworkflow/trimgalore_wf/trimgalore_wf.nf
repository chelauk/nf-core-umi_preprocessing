/*
===================================
        TRIMGALORE
===================================
*/

include { TRIMGALORE }                 from '../../../nf-core/software/trimgalore/main'   addParams(options: params.trimgalore_options)

workflow TRIMGALORE_WF {
    take:
    input_samples

    main:
    TRIMGALORE(input_samples)

    emit:
    trimmed_samples   = TRIMGALORE.out.reads
    trim_qc           = TRIMGALORE.out.zip
}