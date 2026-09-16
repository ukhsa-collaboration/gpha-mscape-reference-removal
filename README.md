# Reference  Removal

An application that takes in reference genome and looks to remove it from fastq datasets. Statistical sumamry will be provided for the reference sequence that has been removed from the sample.

An additional component can be used which will test how effectively a reference is removed from either a supplied dataset or a downloaded dataset.

## Validation
## Input data
- Reference FASTA (required)
- Reference IDX. Index created by Deacon (Optional)
- Samplesheet (Optional)
    - Row one should be a header with column names
    - Column 1: Sample ID; Column 2: Read 1; Column 3: Read 2 (for paired-end only)
- Test FASTQ or Samplesheet (required for validation)


## To Do
- [x] Set up README    
- [x] Set up Nextflow config  
- [x] Read in reference FASTA    
- [x] Read in FASTA | FASTQ to remove reference from  
- [x] Provide test FASTQ to remove reference from
- [x] Generate synthetic reads for reference FASTA [Optional]  
- [x] Spike reference reads into into background FASTQ [Optional]  
- [x] Generate index for Reference FASTA
- [x] Remove reference from background FASTQ
- [x] Generate summary statistics
- [x] Generate summary report
- [ ] Complete process for only removing reference i.e. not validation