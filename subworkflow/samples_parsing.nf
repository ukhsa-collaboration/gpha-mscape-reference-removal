#!/usr/bin/env nextflow
workflow SAMPLES_SETUP{
    take:
    samplesheet_fp
    background_data_dir
    read_length

    main:
    // Declare channels up front so they're visible outside the if/else
    single_end_ch      = Channel.empty()
    paired_end_ch      = Channel.empty()
    reference_genomes_ch = Channel.empty()

    // Check if background data samplesheet has been supplied
    if (samplesheet_fp){
        // Handle either tsv or csv
        def sep_char = samplesheet_fp.toString().toLowerCase().endsWith('.tsv') ? '\t' : ','


        def parsed_ch = Channel
            .fromPath(samplesheet_fp)
            .splitCsv(header: false, skip: 1, sep: sep_char)  // drop header row; column names unknown/irrelevant
            .map { row ->
                if (row.size() < 2) {
                    exit 1, "ERROR: Malformed row in samplesheet (${samplesheet_fp}): ${row}"
                }

                def sample_id = row[0]
                def read1_fp  = file(row[1])
                // Treat a missing/blank third column as single-end
                def read2_raw = (row.size() >= 3) ? row[2]?.trim() : null
                def read2_fp  = (read2_raw) ? file(read2_raw) : null

                if (!read1_fp.exists()) {
                    exit 1, "ERROR: read1 file not found for sample '${sample_id}': ${row[1]}"
                }
                if (read2_fp && !read2_fp.exists()) {
                    exit 1, "ERROR: read2 file not found for sample '${sample_id}': ${read2_raw}"
                }

                tuple(sample_id, read1_fp, read2_fp)
            }

        def branched = parsed_ch.branch { sample_id, read1_fp, read2_fp ->
            paired: read2_fp != null
            single: read2_fp == null
        }

        single_end_ch = branched.single.map { sample_id, read1_fp, read2_fp -> tuple(sample_id, read1_fp) }
        paired_end_ch = branched.paired.map { sample_id, read1_fp, read2_fp -> tuple(sample_id, [read1_fp, read2_fp]) }

        single_end_ch.subscribe { sample_id, read1 -> log.info "Single-end sample: ${sample_id} -> ${read1}" }
        paired_end_ch.subscribe { sample_id, reads -> log.info "Paired-end sample: ${sample_id} -> ${reads.join(', ')}" }   

    } else if(background_data_dir){
        log.info "Specified test data directory ${background_data_dir}. Auto-discovering background FASTQ samples."
        // Generate a channel for each of the FASTQ files in the directory
        log.info "Specified test data directory ${background_data_dir}. Auto-discovering background FASTQ samples."

        def grouped_ch = Channel
            .fromPath("${background_data_dir}/**.{fastq,fq,fastq.gz,fq.gz}")
            .map { fq ->
                // Strip common mate-pair suffixes to get a sample-level grouping key
                def sample_id = fq.getName().replaceAll(/(?i)(?:[._-](?:read)?r?[12](?:_001)?)?\.(fastq|fq)(\.gz)?$/, '')
                tuple(sample_id, fq)
            }
        .groupTuple()

        def branched = grouped_ch.branch { sample_id, fqs ->
            paired: fqs.size() == 2
            single: fqs.size() == 1
        }

        single_end_ch = branched.single.map { sample_id, fqs -> tuple(sample_id, fqs[0]) }
        paired_end_ch = branched.paired.map { sample_id, fqs -> tuple(sample_id, fqs.sort()) } // sort → R1 before R2

        single_end_ch.subscribe { id, fq  -> log.info "Background single-end sample: ${id} -> ${fq}" } 
        paired_end_ch.subscribe { id, fqs -> log.info "Background paired-end sample: ${id} -> ${fqs.join(', ')}" }
    }
    emit:
    single_end          = single_end_ch
    paired_end          = paired_end_ch
}