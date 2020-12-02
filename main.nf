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
bwa_index         = params.bwa               ? file(params.bwa)               : file("${params.outdir}/no_file")
dbsnp             = params.dbsnp             ? file(params.dbsnp)             : file("${params.outdir}/no_file")
dict              = params.dict              ? file(params.dict)              : file("${params.outdir}/no_file")
fasta             = params.fasta             ? file(params.fasta)             : file("${params.outdir}/no_file")
fasta_fai         = params.fasta_fai         ? file(params.fasta_fai)         : file("${params.outdir}/no_file")
germline_resource = params.germline_resource ? file(params.germline_resource) : file("${params.outdir}/no_file")
known_indels      = params.known_indels      ? file(params.known_indels)      : file("${params.outdir}/no_file")
target_bed        = params.target_bed        ? file(params.target_bed)        : file("${params.outdir}/no_file")

// Initialize value channels based on params, not defined within the params.genomes[params.genome] scope
min_reads         = params.min_reads         ?: Channel.empty()  // this the minimum reads parameter passed to FilterConsensusReads
read_structure    = params.read_structure    ?: Channel.empty()

//  include functions
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

include { FASTQC }        from './modules/nf-core/software/fastqc/main' 
include { UMI_STAGE_ONE } from './modules/local/subworkflow/umi_stage_one/umi_stage_one'
include { UMI_STAGE_TWO } from './modules/local/subworkflow/umi_stage_two/umi_stage_two'

workflow {   
    FASTQC(input_samples)
    UMI_STAGE_ONE(input_samples, read_structure, bwa_index, fasta, fasta_fai, dict, min_reads)
    UMI_STAGE_TWO(UMI_STAGE_ONE.out, bwa_index, fasta, fasta_fai,dict)
}