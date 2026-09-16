process SINGLE_REFERENCE_REMOVAL {
        /*
        This process takes single FASTQ input file and a reference index file and removes references reads using Deacon

        Inputs:
            - Sample ID
            - FASTQ
            - Index

        Outputs:
            - FASTQ file
            - Reference Report** Not sure what this looks like currently

    */

    container 'community.wave.seqera.io/library/deacon:0.17.0--43cd5289edd1686c'
    label 'process_medium'
    maxForks 10

    input:
    tuple val(sample_id), path(fastq_fp)
    path(ref_idx)

    output:
    tuple val(sample_id), path("${sample_id}.sample_summary.json"), path("${sample_id}.sample_reads.fq.gz"), emit: single_depleted

    script:
    """
    deacon \
        filter \
        -d ${ref_idx} \
        ${fastq_fp} \
        -s ${sample_id}.sample_summary.json \
        -o ${sample_id}.sample_reads.fq.gz \
        -t ${task.cpus}
    """
}

process SINGLE_SAMPLE_REMOVAL {
        /*
        This process takes single FASTQ input file and a reference index file and removes sample reads using Deacon

        Inputs:
            - Sample ID
            - FASTQ
            - Index

        Outputs:
            - FASTQ file
            - Reference Report** Not sure what this looks like currently

    */

    container 'community.wave.seqera.io/library/deacon:0.17.0--43cd5289edd1686c'
    label 'process_medium'
    maxForks 10

    input:
    tuple val(sample_id), path(fastq_fp)
    path(ref_idx)

    output:
    tuple val (sample_id), path("${sample_id}.reference_summary.json"), path("${sample_id}.reference_reads.fq.gz"), emit: single_ref_only

    script:
    """
    deacon \
        filter \
        ${ref_idx} \
        ${fastq_fp} \
        -s ${sample_id}.reference_summary.json \
        -o ${sample_id}.reference_reads.fq.gz \
        -t ${task.cpus}
    """
}

process PAIRED_REFERENCE_REMOVAL {
        /*
        This process takes paired FASTQ input files and a reference index file and removes references reads using Deacon

        Inputs:
            - SAMPLE ID
            - FASTQ R1
            - FASTQ R2
            - Index


        Outputs:
            - FASTQ file

    */

    container 'community.wave.seqera.io/library/deacon:0.17.0--43cd5289edd1686c'
    label 'process_medium'
    maxForks 10

    input:
    tuple val(sample_id), path(fastq_r1_fp), path(fastq_r2_fp)
    path(ref_idx)

    output:
    tuple val(sample_id), path("${sample_id}.sample_summary.json"), path("${sample_id}.sample_reads.R1.fq.gz"), path("${sample_id}.sample_reads.R2.fq.gz"), emit: paired_depleted


    script:
    """
    deacon \
        filter \
        -d ${ref_idx} \
        ${fastq_r1_fp} \
        ${fastq_r2_fp} \
        -s ${sample_id}.sample_summary.json \
        -o ${sample_id}.sample_reads.R1.fq.gz \
        -o ${sample_id}.sample_reads.R2.fq.gz \
        -t ${task.cpus}
    """
}

process PAIRED_SAMPLE_REMOVAL {
        /*
        This process takes paired FASTQ input files and a reference index file and removes sample reads using Deacon

        Inputs:
            - SAMPLE ID
            - FASTQ R1
            - FASTQ R2
            - Index
        Outputs:
            - FASTQ file

    */

    container 'community.wave.seqera.io/library/deacon:0.17.0--43cd5289edd1686c'
    label 'process_medium'
    maxForks 10

    input:
    tuple val(sample_id), path(fastq_r1_fp), path(fastq_r2_fp)
    path(ref_idx)

    output:
    tuple val(sample_id), path("${sample_id}.reference_summary.json"), path("${sample_id}.reference_reads.R1.fq.gz"), path("${sample_id}.reference_reads.R2.fq.gz"), emit: paired_ref_only

    script:
    """
    deacon \
        filter \
        ${ref_idx} \
        ${fastq_r1_fp} \
        ${fastq_r2_fp} \
        -s ${sample_id}.reference_summary.json \
        -o ${sample_id}.reference_reads.R1.fq.gz \
        -o ${sample_id}.reference_reads.R2.fq.gz \
        -t ${task.cpus}
    """
}