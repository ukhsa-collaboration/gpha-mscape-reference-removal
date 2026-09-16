#!/usr/bin/env nextflow

include { GENERATE_IDX     } from '../modules/generate_idx'
include { SINGLE_SYNTH_READS } from '../modules/synthesise_reads'
include { PAIRED_SYNTH_READS } from '../modules/synthesise_reads'


workflow REFERENCE_PARSING{
    take:
    fasta_fp
    read_type

    main:
    // Check if reference is index file or fasta
    if (fasta_fp.endsWith('.fasta') || fasta_fp.endsWith('.fa') || fasta_fp.endsWith('.fna')) {
        log.info "Generating Deacon index file"
        ref_id = file(fasta_fp).baseName
        log.info "Reference ID: ${ref_id}"
        // Create index
        GENERATE_IDX(fasta_fp, ref_id)
        ref_idx = GENERATE_IDX.out.ref_idx
        
        ref_idx.subscribe { idx ->
            idx_simp = file(idx).baseName
            log.info "Generated Deacon index file: ${idx_simp}"
        }
    } 

    def single_reads = ['single', 'both']
    def paired_reads = ['paired', 'both']
    // Generate synthetic reads for index
    ref_single_synth = null
    ref_paired_synth = null

    if (read_type  in  single_reads){
        SINGLE_SYNTH_READS(fasta_fp, ref_id)
        ref_single_synth = SINGLE_SYNTH_READS.out.ref_single_synth

        ref_single_synth.subscribe { single_ref ->
            log.info "Generated synthetic reference single reads: ${single_ref}"
        }
    }
    if (read_type  in  paired_reads){
        PAIRED_SYNTH_READS(fasta_fp, ref_id)
        ref_paired_synth = PAIRED_SYNTH_READS.out.ref_paired_synth

        ref_paired_synth.subscribe { paired_ref ->
            log.info "Generated synthetic reference paired reads: ${paired_ref.join(', ')}"
        }
    }
    
    

    emit:
    ref_id
    ref_idx
    ref_single_synth
    ref_paired_synth
}