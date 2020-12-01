nextflow.preview.dsl = 2
/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each params in params.genomes, catch the command line first if it was defined
params.bwa                     = params.genome ? params.genomes[params.genome].bwa                     ?: false : false
params.dbsnp                   = params.genome ? params.genomes[params.genome].dbsnp                   ?: false : false
params.dbsnp_index             = params.genome ? params.genomes[params.genome].dbsnp_index             ?: false : false
params.dict                    = params.genome ? params.genomes[params.genome].dict                    ?: false : false
params.fasta                   = params.genome ? params.genomes[params.genome].fasta                   ?: false : false
params.fasta_fai               = params.genome ? params.genomes[params.genome].fasta_fai               ?: false : false
params.germline_resource       = params.genome ? params.genomes[params.genome].germline_resource       ?: false : false
params.germline_resource_index = params.genome ? params.genomes[params.genome].germline_resource_index ?: false : false
params.intervals               = params.genome ? params.genomes[params.genome].intervals               ?: false : false
params.known_indels            = params.genome ? params.genomes[params.genome].known_indels            ?: false : false
params.known_indels_index      = params.genome ? params.genomes[params.genome].known_indels_index      ?: false : false
params.mappability             = params.genome ? params.genomes[params.genome].mappability             ?: false : false

file("${params.outdir}/no_file").text = "no_file\n"
// Initialize files channels based on params, not defined within the params.genomes[params.genome] scope
target_bed        = params.target_bed        ? file(params.target_bed)       : file("${params.outdir}/no_file")

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
bwa               = params.bwa               ? file(params.bwa)               : file("${params.outdir}/no_file")
dbsnp             = params.dbsnp             ? file(params.dbsnp)             : file("${params.outdir}/no_file")
dict              = params.dict              ? file(params.dict)              : file("${params.outdir}/no_file")
fasta             = params.fasta             ? file(params.fasta)             : file("${params.outdir}/no_file")
fasta_fai         = params.fasta_fai         ? file(params.fasta_fai)         : file("${params.outdir}/no_file")
germline_resource = params.germline_resource ? file(params.germline_resource) : file("${params.outdir}/no_file")
known_indels      = params.known_indels      ? file(params.known_indels)      : file("${params.outdir}/no_file")
target_bed        = params.target_bed        ? file(params.target_bed)        : file("${params.outdir}/no_file")

// Initialize value channels based on params, not defined within the params.genomes[params.genome] scope
min_reads         = params.min_reads         ?: Channel.empty()
read_structure    = params.read_structure    ?: Channel.empty()

/* 
 * include umi variant calling functions
 */

include {
   extract_fastq;
   has_extension
} from './modules/functions'

include { FASTQC }                 from './modules/nf-core/software/fastqc/main' addParams( options: [:] )
include { FASTQ_TO_BAM }           from './modules/software/fastq_to_bam/main' addParams( options: [:] )
include { MARK_ILLUMINA_ADAPTERS } from './modules/software/mark_illumina_adapters/main' addParams( options: [:] )
include { BAM_TO_FASTQ }           from './modules/software/bam_to_fastq/main' addParams( options: [:] )
include { UMI_ALN_ONE }            from './modules/software/umi_aln_one/main' addParams( options: [:] )
include { MERGE_BAM_ALIGNMENT }    from './modules/software/merge_bam_alignment/main' addParams( options: [:] )
include { MERGE_RUNS}              from './modules/local/subworkflow/merging/main' addParams( options : [:])
include { GROUP_READS_BY_UMI }     from './modules/software/group_reads_by_umi/main' addParams( options : [:])



// Handle input
tsv_path = null
if (params.input && (has_extension(params.input, "tsv") )) tsv_path = params.input
if (tsv_path) {
    tsv_file = file(tsv_path)
    input_samples = extract_fastq(tsv_file)
}

process sort_bam_one{

    tag "sort ${meta.sample}"

    input:
    tuple val(meta), file(bam_in)
    
    output:
    tuple val(meta), file("*sort1.bam")

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g SortBam \\
    -i $bam_in \\
    -o ${meta.sample}_sort1.bam \\
    --sort-order TemplateCoordinate --max-records-in-ram 4000000
    """
    }

process call_consensus_reads{

    label 'MEMORY_MAX'

    tag "${meta.sample}"

    input:
    tuple val(meta), file(sort1_bam)

    output:
    tuple val(meta), file("*consensus.bam")

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g -XX:+AggressiveOpts -XX:+AggressiveHeap CallMolecularConsensusReads \\
    -i $sort1_bam \\
    -o ${meta.sample}_consensus.bam \\
    --min-reads 1 --min-input-base-quality 30 --tag MI
    """
    }

process filter_consensus_reads {

    tag "${meta.sample}"

    input:
    tuple val(meta), file(consensus_bam) 
    path fasta
    val min_reads

    output:
    tuple meta, file("*_filt.bam")

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g FilterConsensusReads \\
    -i ${consensus_bam} \\
    -o ${meta.sample}_cons_filt.bam \\
    -r ${fasta} \\
    --min-reads $min_reads \\
    --max-read-error-rate 0.05 \\
    --min-base-quality 30 \\
    --max-base-error-rate 0.1 \\
    --max-no-call-fraction 0.1 \\
    --reverse-per-base-tags true
    """
    }

process sam_to_fastq_two {

    tag "${meta.sample}"

    label 'MAX_MEMORY'

    input:
    tuple val(meta), file(bam)

    output:
    tuple val(meta), file("*.fastq")

    script:
    """
    picard -Xmx${task.memory.toGiga()}g SamToFastq \\
    MAX_RECORDS_IN_RAM=4000000 \\
    INPUT=$bam \\
    FASTQ=${meta.sample}.filt.fastq \\
    CLIPPING_ATTRIBUTE=XT \\
    CLIPPING_ACTION=2 \\
    INTERLEAVE=true \\
    NON_PF=true
    """
}

process umi_aln_two {

    tag "bwa ${meta.sample}"

    input:
    tuple val(meta), file(fastq) 
    path bwa
    path fasta
    path fasta_fai

    output:
    tuple meta, file("*sam"), emit: sam

    script:
    """
    bwa mem -M -c 1 -t $task.cpus -k 50 -p $fasta $fastq > ${meta.id}.sam
    """
    }

process sort_bam_two {

    tag "from aln two: ${meta.sample}"

    label 'MEMORY_MAX'

    input:
    tuple val(meta), file(sam)

    output:
    tuple val(meta), file("*sort.sam")

    script:
    """
    picard -Xmx${task.memory.toGiga()}g SortSam \\
    MAX_RECORDS_IN_RAM=4000000 \\
    SORT_ORDER=queryname \\
    INPUT=$sam \\
    OUTPUT=${meta.sample}_aln_merged_umi_group_sort_consensus_filtered_sort.sam \\
    CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
    """
    }

process sort_bam_three {

    tag "from filter consensus ${meta.sample}"

    label 'MEMORY_MAX'

    input:
    tuple val(meta), file(bam)

    output:
    tuple(meta), file("*sort_three.bam")

    script:
    """
    picard -Xmx${task.memory.toGiga()}g SortSam \\
    MAX_RECORDS_IN_RAM=4000000 \\
    SORT_ORDER=queryname \\
    INPUT=$bam \\
    OUTPUT=${meta.sample}.sort_three.bam \\
    CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
    """

    }

process merge_two{
    
    echo true

    tag "second merge: ${meta.sample}"

    label 'MEMORY_MAX'

    input:
    tuple val(meta), file(sort_three_bam), file(sort_two_bam)
    path fasta
    path dict

    output:
    tuple val(meta), file("*bam"), file("*bai")

    script:
    """
    picard -Xmx${task.memory.toGiga()}g MergeBamAlignment \\
    MAX_RECORDS_IN_RAM=4000000 \\
    R=${fasta} \\
    UNMAPPED_BAM=${sort_three_bam} \\
    ALIGNED_BAM=${sort_two_bam} \\
    O=${meta.sample}_aln_umis_final.bam \\
    CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=true CLIP_OVERLAPPING_READS=true \\
    INCLUDE_SECONDARY_ALIGNMENTS=true SO=coordinate MAX_INSERTIONS_OR_DELETIONS=-1 \\
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS
    """
    }


process mark_duplicates {

    tag "mark duplicates ${meta.patient} ${meta.sample}"

    label 'MEMORY_MAX'

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${meta.patient}.${meta.sample}.bam.metrics") "Reports/MarkDuplicates/${meta.patient}/${meta.sample}/${it}"
            else "Preprocessing/bam/duplicates_marked//${meta.patient}/${meta.sample}/${it}"
        }

    input:
    tuple val(meta), file(bam), file(bai)

    output:
    tuple val(meta), file("*md.bam"), file("*md.bai"), emit: md_out
    tuple val(meta), file ("*.bam.metrics"), emit: md_metrics

    script:
    markdup_java_options = task.memory.toGiga() < 8 ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
    metrics = "-M ${meta.patient}.${meta.sample}.bam.metrics"
    //if (params.no_gatk_spark)
    """
    gatk --java-options ${markdup_java_options} \\
        MarkDuplicates \\
        --MAX_RECORDS_IN_RAM 500000 \\
        --INPUT $bam \\
        --METRICS_FILE ${meta.patient}.${meta.sample}.bam.metrics \\
        --TMP_DIR . \\
        --ASSUME_SORT_ORDER coordinate \\
        --CREATE_INDEX true \\
        --OUTPUT ${meta.patient}.${meta.sample}.md.bam
    """
    /*
    else
    """
    gatk --java-options ${markdup_java_options} \\
        MarkDuplicatesSpark \\
        -I ${idPatient}.${idSample}.bam \\
        -O ${idPatient}.${idSample}.md.bam \\
        ${metrics} \\
        --tmp-dir . \\
        --create-output-bam-index true \\
        --spark-master local[${task.cpus}]
    
    """
    */
}

workflow {
    FASTQC(input_samples)
    FASTQ_TO_BAM(input_samples, read_structure)
    MARK_ILLUMINA_ADAPTERS(FASTQ_TO_BAM.out.bam)
    BAM_TO_FASTQ(MARK_ILLUMINA_ADAPTERS.out.bam)
    UMI_ALN_ONE(BAM_TO_FASTQ.out.fastq,bwa,fasta,fasta_fai)
    MERGE_BAM_ALIGNMENT(UMI_ALN_ONE.out.bam.join(MARK_ILLUMINA_ADAPTERS.out.bam),fasta,dict)
    MERGE_RUNS(MERGE_BAM_ALIGNMENT.out)
    GROUP_READS_BY_UMI(MERGE_RUNS.out)
    /*
    sort_bam_one(group_reads_by_umi.out.group_by_umi)
    call_consensus_reads(sort_bam_one.out)
    filter_consensus_reads(call_consensus_reads.out,fasta, min_reads)
    sam_to_fastq_two(filter_consensus_reads.out)
    umi_aln_two(sam_to_fastq_two.out,bwa,fasta,fasta_fai)
    sort_bam_two(umi_aln_two.out)
    sort_bam_three(filter_consensus_reads.out)
    merge_two(sort_bam_three.out.join(sort_bam_two.out), fasta, dict)
    mark_duplicates(merge_two.out)
    */
}