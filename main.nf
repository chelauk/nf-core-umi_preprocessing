nextflow.enable.dsl=2

/*
=================================
          PRINT HELP
=================================
*/

def json_schema = "$baseDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/rnaseq --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

/*
=================================
       PARAMETER SUMMARY
=================================
*/

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

/*
=================================
       PARAMETER CHECKS
=================================
*/

Checks.aws_batch(workflow, params) // Check AWS batch settings
Checks.hostname(workflow, params, log)  // Check the hostnames against configured profiles

// MultiQC - Stage config files

multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)
/*
================================================================================
                     UPDATE MODULES OPTIONS BASED ON PARAMS
================================================================================
*/
modules = params.modules
/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each params in params.genomes, catch the command line first if it was defined
params.bwa                     = params.genome ? params.genomes[params.genome].bwa         ?: false : false
params.dbsnp                   = params.genome ? params.genomes[params.genome].dbsnp       ?: false : false
params.dbsnp_index             = params.genome ? params.genomes[params.genome].dbsnp_index ?: false : false
params.dict                    = params.genome ? params.genomes[params.genome].dict        ?: false : false
params.fasta                   = params.genome ? params.genomes[params.genome].fasta       ?: false : false
params.fasta_fai               = params.genome ? params.genomes[params.genome].fasta_fai   ?: false : false
file("${params.outdir}/no_file").text = "no_file\n"

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
bwa_index         = params.bwa               ? file(params.bwa)               : file("${params.outdir}/no_file")
dbsnp             = params.dbsnp             ? file(params.dbsnp)             : file("${params.outdir}/no_file")
dbsnp_index       = params.dbsnp_index       ? file(params.dbsnp_index)       : file("${params.outdir}/no_file")
dict              = params.dict              ? file(params.dict)              : file("${params.outdir}/no_file")
fasta             = params.fasta             ? file(params.fasta)             : file("${params.outdir}/no_file")
fasta_fai         = params.fasta_fai         ? file(params.fasta_fai)         : file("${params.outdir}/no_file")
target_bed        = params.target_bed        ? file(params.target_bed)        : file("${params.outdir}/no_file")

// Initialize value channels based on params, not defined within the params.genomes[params.genome] scope
min_reads         = params.min_reads         ?: Channel.empty()  // this the minimum reads parameter passed to FilterConsensusReads
read_structure    = params.read_structure    ?: Channel.empty()
library           = params.library           ?: Channel.empty()
params.enable_conda = false
params.second_file  = false

//  include functions
include {
    extract_fastq;
    extract_bam;
    has_extension
} from './modules/functions'

// Handle input
tsv_path = null
if (params.input && (has_extension(params.input, "tsv")) && (params.stage != "two") && (params.stage != "merged") ) tsv_path = params.input
if (tsv_path) {
    tsv_file = file(tsv_path)
    input_samples = extract_fastq(tsv_file)
    }

// handle bam input for stage 2
bam_tsv_path = null
if (params.input && (has_extension(params.input, "tsv")) && (params.stage == "two") ) bam_tsv_path = params.input
if (bam_tsv_path) {
    bam_tsv_file = file(bam_tsv_path)
    input_samples = extract_bam(bam_tsv_file)
    }
// handle bam input for stage 2
merged_tsv_path = null
if (params.input && (has_extension(params.input, "tsv")) && (params.stage == "merged") ) merged_tsv_path = params.input
if (merged_tsv_path) {
    merged_tsv_file = file(merged_tsv_path)
    input_samples = extract_bam(merged_tsv_file)
    }

// params summary for MultiQC
workflow_summary = Schema.params_summary_multiqc(workflow, summary_params)
workflow_summary = Channel.value(workflow_summary)

//include { TRIMGALORE_WF } from './modules/local/subworkflow/trimgalore_wf/trimgalore_wf' addParams (
//    trimgalore_options:                   modules['trimgalore'],
//    trim_umi_options:                     modules['trim_umi']
//	)

include { UMI_STAGE_ONE } from './modules/local/subworkflow/umi_stage_one/umi_stage_one' addParams(
    bed_to_intervals_options:             modules['bed_to_intervals'],
    bwamem1_mem_options:                  modules['bwa_mem1_mem'],
    fastq_to_bam_options:                 modules['fastq_to_bam_mapping'],
    estimate_complexity_options:          modules['estimate_complexity'],
    mark_adapters_options:                modules['mark_illumina_adapters_mapping'],
    bam_to_fastq_options:                 modules['bam_to_fastq_mapping'],
    picard_merge_bams_options:            modules['picard_merge_bams_mapping'],
    merge_runs_mapping_options:           modules['merge_runs_mapping'],
    collect_hs_metrics_options:           modules['picard_hs_metrics'],
    error_rate_options:                   modules['fgbio_error_rate'],
    group_reads_mapping_options:          modules['group_reads_mapping'],
    fgbio_sort_mapping_options:           modules['fgbio_sort_mapping'],
    fgbio_call_consensus_mapping_options: modules['fgbio_call_consensus_mapping'],
    fgbio_filter_mapping_options:         modules['fgbio_filter_mapping']
)

include { POST_MERGE } from './modules/local/subworkflow/post_merge/post_merge' addParams(
    bed_to_intervals_options:             modules['bed_to_intervals'],
    collect_hs_metrics_options:           modules['picard_hs_metrics'],
    error_rate_options:                   modules['fgbio_error_rate'],
    group_reads_mapping_options:          modules['group_reads_mapping'],
    fgbio_sort_mapping_options:           modules['fgbio_sort_mapping'],
    fgbio_call_consensus_mapping_options: modules['fgbio_call_consensus_mapping'],
    fgbio_filter_mapping_options:         modules['fgbio_filter_mapping']
)

include { UMI_STAGE_TWO } from './modules/local/subworkflow/umi_stage_two/umi_stage_two' addParams(
    bed_to_intervals_options:             modules['bed_to_intervals'],
    bwamem1_mem_options:                  modules['bwa_mem1_mem'],
    bam_to_fastq_options:                 modules['bam_to_fastq_mapping'],
    picard_sort_mapping_options:          modules['picard_sort_bams_mapping'],
    picard_merge_bams_options:            modules['picard_merge_bams_mapping'],
    sambamba_sort_options:                modules['sambamba_sort'],
    collect_umi_metrics_options:          modules['umi_aware_mark_duplicates'],
    qualimap_bamqc_mapping_options:       modules['qualimap_bamqc_mapping'],
    collect_hs_metrics_options:           modules['picard_hs_metrics'],
    error_rate_options:                   modules['fgbio_error_rate']
)

include { UMI_QC }       from './modules/local/subworkflow/umi_qc/umi_qc'                addParams(
    fastqc_options:                       modules['fastqc'],
	multiqc_options:                      modules['multiqc']
)

include { UMI_QC_2 }       from './modules/local/subworkflow/umi_qc_2/umi_qc_2'

include { UMI_QC_MERGED}   from './modules/local/subworkflow/umi_qc_merged/umi_qc.nf'

workflow {
    if (params.stage == 'merged') {
        POST_MERGE (input_samples, bwa_index, fasta, fasta_fai, dict, min_reads, target_bed, dbsnp, dbsnp_index)
        filtered_bam = POST_MERGE.out.filtered_bam
    }
    if ( params.stage != 'two' && params.stage != 'merged' ) {
        UMI_STAGE_ONE(input_samples, library, read_structure, bwa_index, fasta, fasta_fai, dict, min_reads, target_bed, dbsnp, dbsnp_index)
        filtered_bam = UMI_STAGE_ONE.out.filtered_bam
        }
    if ( params.stage == 'two' ) { filtered_bam = input_samples }
    UMI_STAGE_TWO(filtered_bam, bwa_index, fasta, fasta_fai,target_bed, dict, dbsnp, dbsnp_index)
    if ( params.stage != 'two' && params.stage != "merged") {
        UMI_QC(
            input_samples,
            multiqc_config,
            multiqc_custom_config,
            workflow_summary,
            UMI_STAGE_ONE.out.hs_metrics,
            UMI_STAGE_TWO.out.md_hs_metrics,
            UMI_STAGE_ONE.out.error_rate,
            UMI_STAGE_ONE.out.group_metrics,
            UMI_STAGE_TWO.out.md_report,
            UMI_STAGE_TWO.out.error_rate_2,
            UMI_STAGE_TWO.out.bamqc_out
            )
        }
    if ( params.stage == "merged") {
        UMI_QC_MERGED(
            multiqc_config,
            multiqc_custom_config,
            workflow_summary,
            POST_MERGE.out.hs_metrics,
            UMI_STAGE_TWO.out.md_hs_metrics,
            POST_MERGE.out.error_rate,
            POST_MERGE.out.group_metrics,
            UMI_STAGE_TWO.out.md_report,
            UMI_STAGE_TWO.out.error_rate_2,
            UMI_STAGE_TWO.out.bamqc_out
        )
    }
    if (params.stage == 'two') {
            UMI_QC_2(
            multiqc_config,
            multiqc_custom_config,
            UMI_STAGE_TWO.out.md_report,
            UMI_STAGE_TWO.out.error_rate_2
            )
        }
}
