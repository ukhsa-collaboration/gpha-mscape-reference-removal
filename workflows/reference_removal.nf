#!/usr/bin/env nextflow
include { GENERATE_IDX } from '../modules/generate_idx.nf'

workflow REFERENCE_REMOVAL{
    take:
    reference
    background_samplesheet

    // Check reference input
    // main:
    // Check if reference is index file or fasta

}