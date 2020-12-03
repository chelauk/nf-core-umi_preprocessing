/*
================================================================================
                                MERGING
================================================================================
*/

include { SAMTOOLS_MERGE_BAM }              from '../../software/samtools/merge_bam/merge_bam'     

workflow MERGE_RUNS {
    take:
       take: bam_bwa

    main:
    //    bam_bwa = MERGE_BAM_ALIGNMENT.out.bam
        bam_bwa.map{ meta, bam ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            [patient, sample, gender, status, bam]
        }.groupTuple(by: [0,1])
            .branch{
                single:   it[4].size() == 1
                multiple: it[4].size() > 1
            }.set{ bam_bwa_to_sort }

        bam_bwa_single = bam_bwa_to_sort.single.map {
            patient, sample, gender, status, bam ->

            def meta = [:]
            meta.patient = patient
            meta.sample = sample
            meta.gender = gender[0]
            meta.status = status[0]
            meta.id = sample

            [meta, bam[0]]
        }

        bam_bwa_multiple = bam_bwa_to_sort.multiple.map {
            patient, sample, gender, status, bam ->

            def meta = [:]
            meta.patient = patient
            meta.sample = sample
            meta.gender = gender[0]
            meta.status = status[0]
            meta.id = sample

            [meta, bam]
        }

        // STEP 1.5: MERGING AND INDEXING BAM FROM MULTIPLE LANES 
        SAMTOOLS_MERGE_BAM(bam_bwa_multiple)
        bam_mapped       = bam_bwa_single.mix(SAMTOOLS_MERGE_BAM.out.bam)

    emit:
        bam = bam_mapped
}