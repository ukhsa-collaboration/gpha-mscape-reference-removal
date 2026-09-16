# Running 
```
nextflow run \
    -latest \
    -r fea/background_data_setup \
    nfellaby/reference-removal \
    -profile docker 
```
- Created test folder
- Added GCA_000854365.1 (TMV reference fasta)

- Create test samplesheet. Includes sample_id, and fasta location
