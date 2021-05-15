/*
===================================
          UMI QC
===================================
*/

include { MULTIQC_2 }          from '../../software/multiqc_2/multiqc_2'

qc_reports          = Channel.empty()

workflow UMI_QC_2 {
    take:
    multiqc_config
    multiqc_custom_config
    md_report
    error_rate_2

    main:
    MULTIQC_2(multiqc_config,
            multiqc_custom_config.ifEmpty([]),
            md_report.collect{it[1]}.ifEmpty([]),
            error_rate_2.collect{it[1]}.ifEmpty([]))
    
}
