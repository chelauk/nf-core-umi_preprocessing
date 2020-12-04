/*
===================================
          UMI QC
===================================
*/

include { FASTQC }           from '../../software/fastqc/fastqc/main'            addParams(options: params.fastqc_options)
include { MULTIQC }          from '../../software/multiqc/multiqc'           

qc_reports          = Channel.empty()

workflow UMI_QC {
    take:
    input_samples
    multiqc_config
    multiqc_custom_config
    workflow_summary


    main:
    FASTQC(input_samples)
    MULTIQC(multiqc_config,
            multiqc_custom_config.ifEmpty([]),
            workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    
}