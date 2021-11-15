// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MULTIQC {
    echo true
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "bioconda::multiqc=1.9" : null)
    //container "quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"

    input:
    path multiqc_config
    path multiqc_custom_config
    path  workflow_summary
    path ('fastqc/*')
    path ('hs_metrics/*')
    path ('md_hs_metrics/*')
    path ('error_rate/*')
    path ('group_metrics/*')
	path ('duplex_seq_metrics/*')
    path ('md_metrics/*')
    path ('error_rate_2/*')
    path ('hs_metrics_2/*')
    path ('bamqc_out/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots

    script:
    //title = custom_runName ? "--title \"${custom_runName}\"" : ''
    //filename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config = params.multiqc_config ? "--config ${multiqc_custom_config}" : ''
    """
    multiqc -f $options.args . --config $multiqc_config
    """
    stub:
    """
    touch stub_multiqc_report.html
    touch stub_data
    """
}

