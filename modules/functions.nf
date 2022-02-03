/*
 * This file holds several functions used to perform operation in umi
 */
 
// Check if a row has the expected number of item
def check_number_of_item(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
}

// Check parameter existence
def check_parameter_existence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

// Compare each parameter with a list of parameters
def check_parameter_list(list, realList) {
    return list.every{ check_parameter_existence(it, realList) }
}


// Channeling the TSV file containing BAM.
// Format is: "patient gender status sample bam bai"
def extract_bam(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            check_number_of_item(row, 5)
            def meta = [:]
            meta.patient = row[0]
            meta.gender  = row[1]
            meta.status  = return_status(row[2].toInteger())
            meta.sample  = row[3]
            meta.id      = "${meta.patient}_${meta.sample}"
            def bam      = return_file(row[4])
            def bai      = return_file(row[5])
            if (!has_extension(bam, "bam")) exit 1, "File: ${bam} has the wrong extension. See --help for more information"
            if (!has_extension(bai, "bai")) exit 1, "File: ${bai} has the wrong extension. See --help for more information"
            return [meta, bam, bai]
        }
}

// Channeling the TSV file containing FASTQ or BAM
// Format is: "patient gender status sample lane fastq1 fastq2"
// or: "patient gender status sample lane bam"
def extract_fastq(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            def meta = [:]
            meta.patient = row[0]
            meta.gender  = row[1]
            meta.status  = return_status(row[2].toInteger())
            meta.sample  = row[3]
            meta.run     = row[4]
            meta.id      = "${meta.patient}_${meta.sample}"
            def read1    = return_file(row[5])
            def read2    = "null"
            def read3    = "null"
            if (has_extension(read1, "fastq.gz") || has_extension(read1, "fq.gz") || has_extension(read1, "fastq") || has_extension(read1, "fq")) {
                check_number_of_item(row, 8)
                read2 = return_file(row[6])
                read3 = return_file(row[7])
            if (!has_extension(read2, "fastq.gz") && !has_extension(read2, "fq.gz")  && !has_extension(read2, "fastq") && !has_extension(read2, "fq")) exit 1, "File: ${file2} has the wrong extension. See --help for more information"
        }
        else if (has_extension(read1, "bam")) check_number_of_item(row, 6)
        else exit 1, "No recognisable extention for input file: ${read1}"

        return [meta, [read1, read2, read3]]
    }
}


// Channeling the TSV file containing Recalibration Tables.
// Format is: "patient gender status sample bam bai recalTable"
def extract_recal(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            check_number_of_item(row, 7)
            def meta = [:]

            meta.patient = row[0]
            meta.gender  = row[1]
            meta.status  = return_status(row[2].toInteger())
            meta.sample  = row[3]
            meta.id      = "${meta.patient}_${meta.sample}"
            def bam      = return_file(row[4])
            def bai      = return_file(row[5])
            def table    = return_file(row[6])

            if (!has_extension(bam, "bam")) exit 1, "File: ${bam} has the wrong extension. See --help for more information"
            if (!has_extension(bai, "bai")) exit 1, "File: ${bai} has the wrong extension. See --help for more information"
            if (!has_extension(table, "recal.table")) exit 1, "File: ${table} has the wrong extension. See --help for more information"

            return [meta, bam, bai, table]
        }
}



// Check file extension
def has_extension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Return file if it exists
def return_file(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

// Remove .ann .gz and .vcf extension from a VCF file
def reduce_vcf(file) {
    return file.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")
}

// Return status [0,1]
// 0 == Normal, 1 == Tumor
def return_status(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
}
