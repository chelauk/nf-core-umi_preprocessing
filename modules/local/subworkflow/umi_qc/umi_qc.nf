/*
===================================
          UMI QC
===================================
*/

include { FASTQC }           from '../../software/fastqc/fastqc/main'            addParams(options: params.fastqc_options)
include { MULTIQC }          from '../../software/multiqc/multiqc'               addParams(options: params.multiqc_options)

qc_reports          = Channel.empty()

workflow UMI_QC {
    take:
    input_samples
//    trimmed_samples
    multiqc_config
    multiqc_custom_config
    workflow_summary
    hs_metrics
    md_hs_metrics
    error_rate
    group_metrics
    md_report
    error_rate_2

    main:
    FASTQC(input_samples)
    MULTIQC(multiqc_config,
            multiqc_custom_config.ifEmpty([]),
            workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            FASTQC.out.zip.collect{it[1]}.ifEmpty([]),
//            trimmed_samples.collect{it[1]}.ifEmpty([]),
            hs_metrics.collect{it[1]}.ifEmpty([]),
            md_hs_metrics.collect{it[1]}.ifEmpty([]),
            group_metrics.collect{it[1]}.ifEmpty([]),
            md_report.collect{it[1]}.ifEmpty([]),
            error_rate.collect{it[1]}.ifEmpty([]),
            error_rate_2.collect{it[1]}.ifEmpty([]))
    
}
