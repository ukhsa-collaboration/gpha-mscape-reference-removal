#!/usr/bin/env nextflow
include {REFERENCE_VALIDATION} from './workflows/reference_validation'
include {REFERENCE_REMOVAL} from './workflows/reference_removal'


workflow {
    // Handle validation user input
    if (!(params.validation instanceof Boolean)) {
        exit 1, "ERROR: --validation must be true or false (got '${params.validation}')"
    }

    // Handle validation process if requested:
    if (params.validation) {
        if (!params.fasta) {
            exit(1, "Validation requires --fasta (synthetic reference reads must be generated from sequence).")
        }
        def valid_read_type = ['paired', 'single', 'both']
        if (!(params.read_type in valid_read_type)) {
            exit 1, "ERROR: --read_type must be one of ${valid_read_type.join(', ')} (got '${params.read_type}')"
        }
        if ((!params.samplesheet && !params.sample_data_dir) ||
            (params.samplesheet && params.sample_data_dir)) {
            exit(1, "Validation requires exactly one of --samplesheet or --sample_data_dir")
        }
        log.info "Running reference validation with ${params.fasta}"
        REFERENCE_VALIDATION(params.fasta, params.idx, params.samplesheet, params.sample_data_dir, params.read_type)

    } else { // Handle standard reference removal
        if ((!params.fasta && !params.idx) || (params.fasta && params.idx)) {
            exit(1, "Reference removal requires exactly one of --fasta or --idx")
        }
        if ((!params.samplesheet && !params.sample_data_dir) ||
            (params.samplesheet && params.sample_data_dir)) {
            exit(1, "Reference removal requires exactly one of --samplesheet or --sample_data_dir")
        }
        log.info "Running reference removal with ${params.fasta ?: params.idx}"
        log.info "Note: --read_type is not used for plain reference removal -- single-end and paired-end samples are both auto-detected and processed regardless of what it's set to."
        REFERENCE_REMOVAL(params.fasta, params.idx, params.samplesheet, params.sample_data_dir)
    }
}