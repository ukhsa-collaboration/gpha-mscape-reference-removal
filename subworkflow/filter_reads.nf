// subworkflow/filter_reads.nf
#!/usr/bin/env nextflow

include { SINGLE_REFERENCE_REMOVAL; SINGLE_SAMPLE_REMOVAL     } from '../modules/filter'
include { PAIRED_REFERENCE_REMOVAL; PAIRED_SAMPLE_REMOVAL } from '../modules/filter'

workflow FILTER_READS {
    take:
    spiked_single   // tuple(sample_id, fastq)         from SPIKE_LONG_READS  (or Channel.empty())
    spiked_paired  // tuple(sample_id, r1, r2)         from SPIKE_PAIRED_READS (or Channel.empty())
    ref_idx       // single path, from REFERENCE_PARSING.out.ref_idx

    main:
    // .first() — same fix as the AIBLAST queue-vs-value issue: broadcasts
    // the one reference index to every sample rather than being consumed once
    def idx_ch = ref_idx.first()

    SINGLE_REFERENCE_REMOVAL(spiked_single, idx_ch)
    SINGLE_SAMPLE_REMOVAL(spiked_single, idx_ch)
    PAIRED_REFERENCE_REMOVAL(spiked_paired, idx_ch)
    PAIRED_SAMPLE_REMOVAL(spiked_paired, idx_ch)

    emit:
    single_depleted = SINGLE_REFERENCE_REMOVAL.out.single_depleted   // depleted background — should be reference-free
    single_ref_only = SINGLE_SAMPLE_REMOVAL.out.single_ref_only   // isolated matches — should be ~exactly the spiked reference reads
    paired_depleted  = PAIRED_REFERENCE_REMOVAL.out.paired_depleted
    paired_ref_only  = PAIRED_SAMPLE_REMOVAL.out.paired_ref_only
}