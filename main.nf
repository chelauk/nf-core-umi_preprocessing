nextflow.preview.dsl = 2

params.outdir = "./results"
params.input = "input.tsv"
params.publish_dir_mode ="copy"
params.read_structure = "+T +M +T"
ch_read_structure = params.read_structure ? Channel.value(params.read_structure) : "null"

/* 
 * include umi variant calling functions
 */

include {
   extract_fastq;
   has_extension
} from './modules/functions'

// Handle input
tsv_path = null
if (params.input && (has_extension(params.input, "tsv") )) tsv_path = params.input
if (tsv_path) {
    tsv_file = file(tsv_path)
    input_samples = extract_fastq(tsv_file)
}

input_samples = input_samples.dump(tag: "reads")

process fastqc {
    echo true
    tag "FASTQC $meta.id"
    publishDir params.outdir, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    path "fastqc_${meta.patient}-${meta.sample}-${meta.run}_logs"

    script:
    """
    fastqc.sh "${meta.patient}-${meta.sample}-${meta.run}" "$reads"
    """ 
}

process trim_galore {
    label 'process_high'

    tag "TRIM $meta.id"

    publishDir "${params.outdir}/trim_galore/${meta.patient}/${meta.id}", mode: params.publish_dir_mode,
      saveAs: {filename ->
        if (filename.indexOf("_fastqc") > 0) "fastqc/$filename"
        else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
//        else if (params.save_trimmed) filename
        else null
      }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_1.fq"), path("${reads[1]}"), path("*_3.fq"), emit: reads
    path "*.html" ,                                        emit: html optional true
    path "*.txt" ,                                         emit: log
//  path "*.version.txt",                                  emit: version
    path "*.zip" ,                                         emit: zip optional true

    script:
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
      cores = (task.cpus as int) - 4
      if (cores < 1) cores = 1
      if (cores > 4) cores = 4
      }

    // Added soft-links to original fastqs for consistent naming in MultiQC
    // prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    prefix = "${meta.id}"

    //nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''
    """
    [ ! -f  ${prefix}_R1.fastq ] && ln -s ${reads[0]} ${prefix}_R1.fastq
    [ ! -f  ${prefix}_R3.fastq ] && ln -s ${reads[2]} ${prefix}_R3.fastq
    trim_galore \\
        --cores ${cores} \\
        --paired \\
        --fastqc \\
        ${prefix}_R1.fastq ${prefix}_R3.fastq
    mv ${prefix}_R1_val_1.fq ${prefix}_R1_val_1.fq
    mv ${prefix}_R3_val_2.fq ${prefix}_R3_val_3.fq 
    """
}

process fix_trim{
    tag "fix trim $meta.id"
    echo true

    input:
    tuple val(meta), file(reads1), file(reads2), file(reads3)

    output:
    tuple val(meta), path("*R1*paired*.fq"), path("*R2*paired*.fq"), path("*R3*paired*.fq"), emit: reads

    script:
    """
    fastq_pair $reads1 $reads2 
    fastq_pair $reads3 $reads2 
    """
}

process fastq_to_bam {
    
    label 'process_high'

    tag "fastq to bam ${meta.id}"

    input:
    tuple val(meta), file(reads1), file(reads2), file(reads3)
    val(rstructure)

    output:
    tuple val(meta) file("*bam")

    script:
    """  
    mkdir temp
    fgbio --tmp-dir=./temp FastqToBam \\
    --input ${reads1} ${reads2} ${reads3} \\
    --output ${meta.id}_unaln_umi.bam \\
    --read-structures $rstructure \\
    --sort true \\
    --umi-tag RX \\
    --sample ${meta.id} \\
    --library "test" \\
    --read-group-id ${meta.patient}-${meta.id}
    """
    }

process mark_illumina_adapters {

    tag "mark adapters ${meta.id}"

    label 'process_high'

    publishDir "$params.outdir/mark_illumina_adapters/${meta.patient}/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".metrics") > 0) filename
                      else null
                }

    input:
    tuple val(meta), file(bam)

    output:
    tuple val(meta), file("*bam"), emit : bam
    path "*_mark_adaptor.metrics", emit: mark_adaptor_log

    script:
    """
    picard MarkIlluminaAdapters \\
    MAX_RECORDS_IN_RAM=4000000 \\
    INPUT=$bam \\
    OUTPUT="${meta.id}_unaln_umi_marked.bam \\
    M="${meta.patient}_${meta.id}"_mark_adaptor.metrics
    """
    }

workflow {
    fastqc(input_samples)
    trim_galore(input_samples)
    fix_trim(trim_galore.out.reads)
    fastq_to_bam(fix_trim.out.reads, ch_read_structure)
    //mark_illumina_adapters
}