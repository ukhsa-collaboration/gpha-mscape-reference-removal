#!/usr/bin/env nextflow
include { GENERATE_IDX } from '../modules/generate_idx'

workflow RESOLVE_REFERENCE {
    take:
    fasta_fp   // always required for validation
    idx_fp     // optional -- used preferentially if supplied, skips re-indexing

    main:
    def supplied_idx = idx_fp ?: (fasta_fp?.toString()?.endsWith('.idx') ? fasta_fp : null)

    if (supplied_idx) {
        log.info "Using supplied Deacon index directly: ${supplied_idx}"
        ref_id  = file(supplied_idx).baseName
        ref_idx = Channel.fromPath(supplied_idx)
    } else if (fasta_fp) {
        log.info "Generating Deacon index from FASTA: ${fasta_fp}"
        ref_id = file(fasta_fp).baseName
        GENERATE_IDX(fasta_fp, ref_id)
        ref_idx = GENERATE_IDX.out.ref_idx
    } else {
        error "Please supply either --fasta (a FASTA file) or --idx (a prebuilt Deacon index)."
    }

    emit:
    ref_id
    ref_idx
}