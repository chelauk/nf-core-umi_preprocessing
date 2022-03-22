/*
===================================
          UMI QC
===================================
*/

include { MULTIQC }          from '../../software/multiqc/multiqc'               addParams(options: params.multiqc_options)

qc_reports          = Channel.empty()

workflow UMI_QC_MERGED {
    take:
    multiqc_config
    multiqc_custom_config
    workflow_summary
    hs_metrics
    md_hs_metrics
    error_rate
    group_metrics
    md_report
    error_rate_2
    bamqc_out

    main:
    MULTIQC(multiqc_config,
            multiqc_custom_config.ifEmpty([]),
            workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            hs_metrics.collect{it[1]}.ifEmpty([]),
            md_hs_metrics.collect{it[1]}.ifEmpty([]),
            error_rate.collect{it[1]}.ifEmpty([]),
            group_metrics.collect{it[1]}.ifEmpty([]),
            md_report.collect{it[1]}.ifEmpty([]),
            error_rate_2.collect{it[1]}.ifEmpty([]),
            bamqc_out.collect{it[1]}.ifEmpty([]))
    
}
