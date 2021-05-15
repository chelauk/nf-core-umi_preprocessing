// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process ERRORRATE_BY_READ_POSITION {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::fgbio=1.3.0" : null)
    //container "quay.io/biocontainers/fgbio:1.3.0--0"

    input:
    tuple val(meta), path(bam)
    path fasta
    path dict 
    path dbsnp
    path dbsnp_index 
    path interval_list

    output:
    tuple val(meta), file("*.error_rate_by_read_position.txt"), emit: error_rate

    script:
    output_options = params.second_file ? "--output ${meta.patient}_${meta.sample}_st2_qc" : "--output ${meta.patient}_${meta.sample}_st1_qc"
    """
    fgbio -Xmx${task.memory.toGiga()}g ErrorRateByReadPosition \\
    --input ${bam} \\
    ${output_options} \\
    --intervals ${interval_list} \\
    --ref ${fasta} \\
    --variants ${dbsnp}
    """
    stub:
    output_options = params.second_file ? "--output ${meta.patient}_${meta.sample}_st2_qc" : "--output ${meta.patient}_${meta.sample}_st1_qc"
    """
    output="${output_options}"
    touch \${output:8}.error_rate_by_read_position.txt
    """


}