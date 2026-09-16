// subworkflow/report.nf — join everything on sample_id
include { GENERATE_VALIDATION_REPORT } from '../modules/report.nf'
include { AGGREGATE_REPORT           } from '../modules/report.nf'

workflow VALIDATION_REPORT {
    take:
    background_truth    // tuple(sample_id, fastq) or tuple(sample_id, r1, r2) -- pre-spike, from SAMPLES_SETUP
    reference_truth     // single fastq -- from REFERENCE_PARSING (broadcast with .first())
    isolate_out         // tuple(sample_id, json, fastq) -- from SINGLE_SAMPLE_REMOVAL / PAIRED_SAMPLE_REMOVAL
    depleted_out        // tuple(sample_id, json, fastq) -- from SINGLE_REFERENCE_REMOVAL / PAIRED_REFERENCE_REMOVAL
    read_type           // val, 'single' or 'paired'

    main:

    def bg_norm  = background_truth.map { row -> tuple(row[0], row[1..-1].flatten()) }
    def iso_norm = isolate_out.map      { row -> tuple(row[0], row[1], row[2..-1]) }
    def dep_norm = depleted_out.map     { row -> tuple(row[0], row[1], row[2..-1]) }
    def ref_norm = reference_truth.map { ref -> [ref instanceof List ? ref : [ref]] }

    def joined = bg_norm
        .join(iso_norm)
        .join(dep_norm)
        .combine(ref_norm)
        .map { sample_id, bg_fqs, iso_json, iso_fqs, dep_json, dep_fqs, ref_fqs ->
            tuple(sample_id, read_type, ref_fqs, bg_fqs, iso_fqs, dep_fqs, iso_json)
        }

    GENERATE_VALIDATION_REPORT(joined)

    def grouped = GENERATE_VALIDATION_REPORT.out.confusion.collect().map { jsons -> tuple(read_type, jsons) }
    
    AGGREGATE_REPORT(grouped)

    emit:
    report     = AGGREGATE_REPORT.out.report
    report_pdf = AGGREGATE_REPORT.out.report_pdf
}