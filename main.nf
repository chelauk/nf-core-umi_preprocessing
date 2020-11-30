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

// Handle input
tsv_path = null
if (params.input && (has_extension(params.input, "tsv") )) tsv_path = params.input
if (tsv_path) {
    tsv_file = file(tsv_path)
    input_samples = extract_fastq(tsv_file)
}

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

process fastq_to_bam {
    tag "fastq to bam: ${meta.id}"

    input:
    tuple val(meta), file(reads)
    val rstructure 

    output:
    tuple val(meta), file("*bam"), emit: bam

    script:
    """  
    mkdir temp
    fgbio --tmp-dir=./temp FastqToBam \\
    --input $reads \\
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
    OUTPUT="${meta.id}_unaln_umi_marked.bam" \\
    M="${meta.patient}_${meta.id}_mark_adaptor.metrics"
    """
    }

    process  sam_to_fastq {

    tag "sam to fastq ${meta.id}"

    label 'process_high'

    input:
    tuple val(meta), file(bam)

    output:
    tuple val(meta), file("*fastq"), emit: fastq

    script:
    """
    picard SamToFastq \\
    MAX_RECORDS_IN_RAM=4000000 \\
    INPUT=$bam \\
    FASTQ="${meta.id}.fastq" \\
    CLIPPING_ATTRIBUTE=XT \\
    CLIPPING_ACTION=2 \\
    INTERLEAVE=true \\
    NON_PF=true 
    """
    }

process umi_aln {

    tag "bwa ${meta.id}"

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

process merge_one {

    echo true

    tag "merge ${meta.id}"
    
    input:
    tuple val(meta), file(aligned_unmarked_sam), file(unaligned_marked_bam)
    path fasta
    path dict

    output:
    tuple val(meta), file("*aligned_marked.ba{m,i}")
    
    script:
    """
    picard MergeBamAlignment \\
    MAX_RECORDS_IN_RAM=4000000 \\
    R=$fasta \\
    UNMAPPED_BAM=$unaligned_marked_bam \\
    ALIGNED_BAM=$aligned_unmarked_sam \\
    O="${meta.id}.aligned_marked.bam" \\
    CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false \\
    CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true \\
    MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant \\
    ATTRIBUTES_TO_RETAIN=XS
    """ 
}

process MERGE_BAM {
    label 'CPUS_MAX'

    tag "${meta.id}"


    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${meta.sample}.ba{m,i}"), emit: bam

    script:
    """
    samtools merge --threads ${task.cpus} ${meta.sample}.bam ${bam}
    samtools index ${meta.sample}.bam
    """
}

process group_reads_by_umi{

    tag "group reads ${meta.sample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${meta.sample}_aln_merged_umi.metrics") "Reports/group_reads_by_umi/${meta.patient}/${meta.sample}/${it}"
            else null
        }

    input:
    tuple val(meta), file(bam_merged)

    output:
    tuple val(meta), file("*bam"), emit: group_by_umi
    tuple val(meta.patient), val(meta.sample), file("*.metrics"), emit: group_by_umi_metrics

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g GroupReadsByUmi \\
    -i ${bam_merged[1]} \\
    -f ${meta.sample}_aln_merged_umi.metrics \\
    -s adjacency -m 30 -t RX -T MI --min-umi-length 9 \\
    -o ${meta.sample}_aln_merged_umi_group.bam 
   """
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
    fastqc(input_samples)
    fastq_to_bam(input_samples, read_structure)
    mark_illumina_adapters(fastq_to_bam.out.bam)
    sam_to_fastq(mark_illumina_adapters.out.bam)
    umi_aln(sam_to_fastq.out.fastq,bwa,fasta,fasta_fai)
    merge_one(umi_aln.out.sam.join(mark_illumina_adapters.out.bam),fasta,dict)
    merge_one_bam = merge_one.out
    
    merge_one_bam.map{ meta, bam ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            [patient, sample, gender, status, bam]
        }.groupTuple(by: [0,1])
            .branch{
                single:   it[4].size() == 1
                multiple: it[4].size() > 1
            }.set{ merge_one_bam_to_sort }

    merge_one_single = merge_one_bam_to_sort.single.map {
            patient, sample, gender, status, bam ->

            def meta = [:]
            meta.patient = patient
            meta.sample = sample
            meta.gender = gender[0]
            meta.status = status[0]
            meta.id = sample

            [meta, bam[0]]
        }
    
    bam_bwa_multiple = merge_one_bam_to_sort.multiple.map {
            patient, sample, gender, status, bam ->

            def meta = [:]
            meta.patient = patient
            meta.sample = sample
            meta.gender = gender[0]
            meta.status = status[0]
            meta.id = sample

            [meta, bam]
        }
    
    MERGE_BAM(merge_one_bam_to_sort.multiple)
    group_reads_by_umi(merge_one_single)
    sort_bam_one(group_reads_by_umi.out.group_by_umi)
    call_consensus_reads(sort_bam_one.out)
    filter_consensus_reads(call_consensus_reads.out,fasta, min_reads)
    sam_to_fastq_two(filter_consensus_reads.out)
    umi_aln_two(sam_to_fastq_two.out,bwa,fasta,fasta_fai)
    sort_bam_two(umi_aln_two.out)
    sort_bam_three(filter_consensus_reads.out)
    merge_two(sort_bam_three.out.join(sort_bam_two.out), fasta, dict)
    mark_duplicates(merge_two.out)
}