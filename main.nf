#!/usr/bin/env nextflow
include {REFERENCE_VALIDATION} from './workflows/reference_validation'
include {REFERENCE_REMOVAL} from './workflows/reference_removal'


workflow {
    // 1. Check User inputs
    // Check Reference has been supplied
    if (!params.fasta){
        exit(1, "Please specify --fasta FASTA file to use.")
    }

    // Performing Validation?
    if (!(params.validation instanceof Boolean)){
        exit 1, "ERROR: --validation must be one of ${validation_options.join(', ')} (Got '${params.validation}')"
    }

    if (params.validation){
        // Validation always requires fasta
        if (!params.fasta){
            exit(1, "Please specify --reference FASTA file to use.")
        }
        log.info "Running reference validation with ${params.fasta}"
        if (params.idx){
            log.info "Testing with included index: ${params.idx}"
        }

        def valid_read_type = ['paired', 'single', 'both']
        if (!(params.read_type in valid_read_type)){
            exit 1, "ERROR: --read_type must be one of ${valid_read_type.join(', ')} (got '${params.read_type}')"
        }
        log.info "Generating ${params.read_type} read synthetic data. Designated by --read_type, default='single', options='paired','single','both'"

        // Validation requires either samplesheet or test-directory
        if ((!params.samplesheet && !params.sample_data_dir) ||
            (params.samplesheet && params.sample_data_dir)) {
            exit(1, "Validation requires exactly one of --samplesheet or --sample_data_dir")
        }
        REFERENCE_VALIDATION(params.fasta, params.idx, params.samplesheet, params.sample_data_dir, params.read_type)

    }
    else {
        log.info "Running reference removal with ${params.fasta}"
        if(!params.samplesheet){
            exit(1, "When running reference removal please specify --background_samplesheet.")
        }
        else{
            REFERENCE_REMOVAL(params.fasta, params.samplesheet)
        }
        
    }
    // Handle no background specified, download background dataset
    // if (!params.background){
    // }

    // 2. Check FASTA Files
    // 3. Generate synthetic reads
    // 4. Run Filtering process
    // 5. Generate Summary Table outputs
}