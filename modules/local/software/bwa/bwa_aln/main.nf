// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BWA_ALN {
    label 'process_high'

    tag "${meta.id}"
    
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.10" : null )
    //container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:eabfac3657eda5818bae4090db989e3d41b01542-0"


    input:
    tuple val(meta), file(reads) 
    path index
    path fasta
    path fai

    output:
    tuple meta, file("*bam"), emit: bam
    path  "*.version.txt"   , emit: version

    script:
    def software   = getSoftwareName(task.process)
    CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    """
    bwa mem \\
    ${options.args} \\
    -t ${task.cpus} \\
    ${fasta} ${reads} | \\
    samtools sort --threads ${task.cpus} -m 2G -o ${meta.id}.bam
    echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//' > ${software}.version.txt
    """
    }