#!/usr/bin/env nextflow
include { GENERATE_IDX                           } from '../modules/generate_idx'
include { SINGLE_SYNTH_READS; PAIRED_SYNTH_READS } from '../modules/synthesise_reads'

workflow REFERENCE_PARSING{
    take:
    fasta_fp   // always required for validation -- needed to synthesise reference reads
    idx_fp     // optional for validation -- used directly, skips re-indexing
    read_type

    main:
    // fasta_fp is always required (needed to synthesise reference reads
    // regardless of idx_fp), so ref_id is always derived from it.
    ref_id = file(fasta_fp).baseName
    log.info "Reference ID: ${ref_id}"

    // idx_fp is optional -- use it directly if supplied, to skip re-indexing
    if (idx_fp) {
        log.info "Using supplied Deacon index directly: ${idx_fp}"
        ref_idx = Channel.fromPath(idx_fp)
    } else {
        log.info "Generating Deacon index file from FASTA: ${fasta_fp}"
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