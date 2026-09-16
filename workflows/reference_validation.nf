#!/usr/bin/env nextflow
include { REFERENCE_PARSING } from '../subworkflow/reference_parsing'
include { SAMPLES_SETUP     } from '../subworkflow/samples_parsing'
include { FILTER_READS      } from '../subworkflow/filter_reads'
include { VALIDATION_REPORT } from '../subworkflow/report'

def checkBackgroundPresent(ch, String label, boolean required) {
    ch.count().subscribe { n ->
        if (n == 0) {
            def msg = "No ${label} background samples found to spike the synthetic reference into."
            if (required) {
                error "${msg} Cannot proceed with read_type='${params.read_type}'. Supply matching background data, or change read_type."
            } else {
                log.warn "${msg} Skipping ${label} validation; continuing with what is available."
            }
        }
    }
    return ch
}

process SPIKE_SINGLE_READS {
    label 'process_low'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(background_fq)
    path(ref_synth_fq)

    output:
    tuple val(sample_id), path("${sample_id}.spiked.single.fq.gz"), emit: spiked

    script:
    """
    cat ${background_fq} ${ref_synth_fq} > ${sample_id}.spiked.single.fq.gz
    """
}

process SPIKE_PAIRED_READS {
    label 'process_low'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(background_reads)
    tuple path(ref_r1), path(ref_r2)

    output:
    tuple val(sample_id), path("${sample_id}.spiked.R1.fq.gz"), path("${sample_id}.spiked.R2.fq.gz"), emit: spiked

    script:
    """
    cat ${background_reads[0]} ${ref_r1} > ${sample_id}.spiked.R1.fq.gz
    cat ${background_reads[1]} ${ref_r2} > ${sample_id}.spiked.R2.fq.gz
    """
}

workflow REFERENCE_VALIDATION{
    take:
    fasta
    idx
    background_samplesheet_fp
    background_data_dir
    read_type

    main:
    def single_reads  = ['single', 'both']
    def paired_reads = ['paired', 'both']
    def strict_both = !(params.allow_partial_validation ?: false)

    // Generate index files for reference
    REFERENCE_PARSING(fasta, read_type)
    // Set up background data
    SAMPLES_SETUP(background_samplesheet_fp, background_data_dir, read_type)

    
    def spiked_long_ch  = Channel.empty()
    def spiked_short_ch = Channel.empty()

    // Spike in syntheised reference reads into test sample(s)
    if (read_type in paired_reads) {
        // 'paired' alone → missing paired background makes the whole run pointless → always required
        // 'both'        → required unless the user opted into partial validation
        def required = (read_type == 'paired') || (read_type == 'both' && strict_both)
        def paired_ch = checkBackgroundPresent(SAMPLES_SETUP.out.paired_end, 'paired-end (short-read)', required)
        SPIKE_PAIRED_READS(paired_ch, REFERENCE_PARSING.out.ref_paired_synth.first())
        spiked_paired_ch = SPIKE_PAIRED_READS.out.spiked
    }

    if (read_type in single_reads) {
        def required = (read_type == 'single') || (read_type == 'both' && strict_both)
        def single_ch = checkBackgroundPresent(SAMPLES_SETUP.out.single_end, 'single-end (long-read)', required)
        SPIKE_SINGLE_READS(single_ch, REFERENCE_PARSING.out.ref_single_synth.first())
        spiked_single_ch = SPIKE_SINGLE_READS.out.spiked
    }

    // Run Reference Removal on Spiked samples
    FILTER_READS(spiked_single_ch, spiked_paired_ch, REFERENCE_PARSING.out.ref_idx)
   
    // --- Reporting: one VALIDATION_REPORT call per read type, since the ---
    // --- underlying (background_truth, isolate_out, depleted_out) shapes  ---
    // --- differ between long (single fastq) and paired (R1+R2 pair)        ---
    if (read_type in single_reads) {
        VALIDATION_REPORT(
            SAMPLES_SETUP.out.single_end,           // background_truth, pre-spike
            REFERENCE_PARSING.out.ref_single_synth.first(),
            FILTER_READS.out.single_ref_only,
            FILTER_READS.out.single_depleted,
            'single',
        )
    }

    if (read_type in paired_reads) {
        VALIDATION_REPORT(
            SAMPLES_SETUP.out.paired_end,
            REFERENCE_PARSING.out.ref_paired_synth.first(),
            FILTER_READS.out.paired_ref_only,
            FILTER_READS.out.paired_depleted,
            'paired',
        )
    }
    

}