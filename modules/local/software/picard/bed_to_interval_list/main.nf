// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BED_TO_INTERVAL_LIST {
    tag "intervals"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "bioconda::picard=2.23.8" : null)
    //container "quay.io/biocontainers/picard:2.23.8--0"

    input:
    path target_bed
    path dict

    output:
    path "interval.list",   emit: interval_list
    path "*.version.txt",   emit: version

    script:
    def software  = getSoftwareName(task.process)
    """
    picard -Xmx${task.memory.toGiga()}g BedToIntervalList \\
    I=${target_bed} \\
    SD=${dict} \\
    O=interval.list
    echo \$(picard  BedToIntervalList --version 2>&1) > ${software}.version.txt
    """
    stub:
    """
    touch interval.list
    touch software.version.txt
    """
    }
