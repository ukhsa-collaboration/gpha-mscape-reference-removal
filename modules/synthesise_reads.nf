#!/usr/bin/env nextflow

process SINGLE_SYNTH_READS {
    /*
        This process takes a FASTA input file and generates Synthetic single reads using
        PBSim

        Inputs:
            - FASTA filepath
        Outputs:
            - FASTQ file

    */

    container 'community.wave.seqera.io/library/pbsim3:3.0.5--86541aa3eccd4c3c'
    label 'process_low'
    maxForks 10

    input:
    path(fasta_fp)
    val(sample_id)

    output:
    path("${sample_id}.pbsim.fq.gz"), emit: ref_single_synth

    script:
    """
    pbsim \
        --prefix ${sample_id} \
        --id-prefix ${sample_id}__ \
        --seed 1 \
        --strategy wgs \
        --genome ${fasta_fp} \
        --depth 10 \
        --method errhmm \
        --errhmm ${params.pbsim_model} \
        --accuracy-mean 0.98;

    cat ${sample_id}_*.fq.gz > ${sample_id}.pbsim.fq.gz
    """

}


process PAIRED_SYNTH_READS {
    /*
        This process takes a FASTA input file and generates synthetic paired reads using
        DWGSIM

        Inputs:
            - FASTA filepath
        Outputs:
            - FASTQ file

    */

    container 'community.wave.seqera.io/library/dwgsim:1.1.14--b4033839f1e4b148'
    label 'process_low'
    maxForks 10

    input:
    path(fasta_fp)
    val(sample_id)

    output:
    tuple path("${sample_id}.bwa.R1.fastq.gz"), path("${sample_id}.bwa.R2.fastq.gz"), emit: ref_paired_synth

    script:
    """
        dwgsim \
        -C 10 \
        -1 150 \
        -2 150 \
        -y 0.0 \
        -o 1 \
        -z 1 \
        -F 0.0 \
        -r 0.0 \
        -e 0.01 \
        -E 0.01 "$fasta_fp" "$sample_id";
    """

}