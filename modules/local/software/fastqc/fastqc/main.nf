// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process FASTQC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
    //container "quay.io/biocontainers/fastqc:0.11.9--0"
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "*.version.txt"          , emit: version

    script:
    // Add soft-links to original FastQs for consistent naming in pipeline
    def software = getSoftwareName(task.process)
	def prefix   = options.suffix ? "${meta.patient}_${meta.id}.${options.suffix}" : "${meta.patient}_${meta.id}_${meta.run}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        fastqc $options.args --threads $task.cpus ${prefix}.fastq.gz
        fastqc --version | sed -e "s/FastQC v//g" > ${software}.version.txt
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        [ ! -f  ${prefix}_3.fastq.gz ] && ln -s ${reads[2]} ${prefix}_3.fastq.gz
        fastqc $options.args --threads $task.cpus ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${prefix}_3.fastq.gz
        fastqc --version | sed -e "s/FastQC v//g" > ${software}.version.txt
        """
    }

    stub:
	def prefix   = options.suffix ? "${meta.patient}_${meta.id}.${options.suffix}" : "${meta.patient}_${meta.id}_${meta.run}"
    """
    touch ${prefix}.html
    touch ${prefix}.zip
    touch fastq.version.txt
    """
}
