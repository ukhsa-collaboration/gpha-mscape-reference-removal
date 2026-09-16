# Reference Removal

## Description
The project looks to evaluate methods of removing references from metagenomic datasets. Specifically, internal controls that have been spiked-in i.e. TMV, lambda phage etc.

## De-Hosting Tools
- [Hostile](https://github.com/bede/hostile): Hostile removes host sequences from short and long read (meta)genomes, consuming single or paired FASTQ from files or stdin.
- [Deacon](https://github.com/bede/deacon): Deacon filters DNA sequences in FASTA/Q files and streams using SIMD-accelerated minimizer comparison with query sequence(s), emitting either matching sequences (search mode), or sequences without matches (deplete mode).
- [NoHuman](https://github.com/mbhall88/nohuman/): nohuman removes human reads from sequencing reads by classifying them with kraken2 against a custom database built from all of the genomes in the Human Pangenome Reference Consortium's ( HPRC) second release. It can take any type of sequencing technology.
- [detaxizer](https://github.com/nf-core/detaxizer/tree/dev): A pipeline to identify (and remove) certain sequences from raw genomic data. Default taxon to identify (and remove) is Homo sapiens. Removal is optional.

## Reference Datasets
The following datasets were canabilised from the Deacon pre-pub [paper](https://www.biorxiv.org/content/10.1101/2025.06.09.658732v1.full.pdf):
- [Viral Genomes](https://zenodo.org/records/15411280): 14,861 complete NCBI RefSeq virus sequences downloaded on 2026-07-31
- [Bacterial genomes](https://zenodo.org/records/15424142): all 1,428 complete [FDA-ARGOS bacterial reference genomes]((https://www.nature.com/articles/s41467-019-11306-6)) including plasmids, downloaded on 2026-07-31.
- Fungal Genomes: all 5,488 complete NCBI RefSeq, downloaded on 2026-07-31.
- Phage Genomes: Curated list of 23,122 genomes down-sampled from all Caudoviricetes on genbank that pass basic QC like <10% Ns etc. (via AK)

## Synthetic Read Generation Tools
- Short Read, [DWGSIM](https://github.com/nh13/DWGSIM)
- Long Read, [PBSIM](https://github.com/yukiteruono/pbsim3)

## Reference Organisms for Internal Control Spike-in
Viral:
- Tobacco Mosaic Virus ([Adela Alcolea-Medina et al., 2025](https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(25)00102-8/fulltext))
- Hazara Virus ([Kuiama Lewandowski et al., 2019](https://pubmed.ncbi.nlm.nih.gov/31666364/))
- Lambda Phage - NC_001416.1 ([Jiayi Duan et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12633245/))
- T1 Phage ([Zhangfan Fu et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10714923/))
- Thermus thermophilus ([Zhangfan Fu et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10714923/))

Bacterial FASTA References 
- [ZymoBIOMICs references](https://www.bioscience.co.uk/cpl/zymobiomics-dna-kits?gad_source=1&gad_campaignid=21301459676&gbraid=0AAAAAD_TCRJyensW6-q0ajdbGar2qCICn&gclid=CjwKCAjwj7HTBhBiEiwA8s35OmX-N4EqACiysya4QZFb7_7QZMjsWO357qzxvQ40A8Ae7OlWGuWCbhoC0j4QAvD_BwE):
  - Imtechella halotolerans
  - Allobacillus halotolerans
  - Truepera radiovictrix

[Predefined by mSCAPE](https://github.com/artic-network/scylla/tree/main/resources/spike_ins):  
spike-in:
- ERCC-RNA_4456740
- bacillus_ms2phage
- ms2-phage
- phix
- tobacco_mosaic_virus
- zymo_D6320
- zymo_D6321

Internal Controls:
- [NIBSC_11/242](https://nibsc.org/products/brm_product_catalogue/detail_page.aspx?catid=11/242): Mock community containing 25 human pathogenic viruses
- [NIBSC_20/170](https://nibsc.org/about_us/latest_news/winter_multiplex.aspx): Flu A (H1N1, H3N2), Flu B, RSV A, RSV B, SARS-CoV-2
- bacillus_ms2phage
- zepto_rp2.1: Adenovirus 1, 3, 31; C. pneumoniae; Influenza A 2009 H1N1pdm, H3N2; Metapneumovirus 8; M. pneumoniae; Parainfluenza Type 1, 4; Rhinovirus 1A; SARS-CoV-2; B. parapertussis; B. pertussis; Coronavirus 229E, HKU-1, BL63, OC43; Influenza AH1, Influenza B, Parainfluenza Type 2, 3; RSV A
- [zymo-mc_D6300](https://zymoresearch.eu/collections/zymobiomics-microbial-community-standards/products/zymobiomics-microbial-community-standard): Listeria monocytogenes - 12%, Pseudomonas aeruginosa - 12%, Bacillus subtilis - 12%, Escherichia coli - 12%, Salmonella enterica - 12%, Lactobacillus fermentum - 12%, Enterococcus faecalis - 12%, Staphylococcus aureus - 12%, Saccharomyces cerevisiae - 2%, and Cryptococcus neoformans - 2%.



## Creating Reference Datasets
Given the references will be in the refseq/argos downloads, need to filter out references from those datasets. 
### 1. Remove References from Dataset
1. Run Kraken2 and de-hosting tools for reference
    - Evaluate number of hits
2. Remove matched reads 
3. Create version for each reference without reference, i.e. NCBI_refseq_visues_20260731.noTMV.fasta
4. Re-Run to see if any additional reads are annotated.

### 2.1 Generate Synthetic reads for viral datasets and references
Again, canabalised from Deacon pre-pub [paper](https://www.biorxiv.org/content/10.1101/2025.06.09.658732v1.full.pdf)
1. Long Read
   - Model: ERRHMM-ONT-HQ
   - Depth: 10x
   - Mean read length: 1,000bp
   - Max read length: 10,000bp
   - Mean accuracy: 0.98
   - Random seed: 1

Command used:
```
 for fasta in rsviruses17900/*.fa; do
    acc=$(basename "$fasta" .fa) \ 
    pbsim \
        --seed 1 \
        --strategy wgs \
        --method errhmm \
        --errhmm pbsim3/data/ERRHMM-ONT-HQ.model \
        --depth 10 \
        --genome ${fasta} \
        --prefix ${acc} \
        --id-prefix ${acc}__ \
        --length-mean 1000 \
        --length-max 10000 \
        --accuracy-mean 0.98; \
    cat ${acc}*.fastq | pigz > ${acc}.fastq.gz 
done;
```
2. Short Read
   - Read length: 2x150bp (paired)
   - Depth: 10x
   - Random read probability (-y): 0
   - Error rate (-e and -E): 0.01
   - Mutation rate (-r): 0.0 (of which low frequency somatic mutations (-F): 0.0)
   - Random seed (-z): 1

Command used:
```
for fasta in rsviruses17900/*.fa; do \
    acc=$(basename "$fasta" .fa) \
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
        -E 0.01 "$fasta" "$acc";
done;
```

### 2.2 Generate Synthetic reads for bacterial datasets and references
Again, canabalised from Deacon pre-pub [paper](https://www.biorxiv.org/content/10.1101/2025.06.09.658732v1.full.pdf)
1. Long Read
   - Model: ERRHMM-ONT-HQ
   - Depth: 10x
   - Mean read length: 5,000bp
   - Max read length: 50,000bp
   - Mean accuracy: 0.98
   - Random seed: 1

Command used:
```
 for fasta in bacteria/*.fa; do
    acc=$(basename "$fasta" .fa) \ 
    pbsim \
        --seed 1 \
        --strategy wgs \
        --method errhmm \
        --errhmm pbsim3/data/ERRHMM-ONT-HQ.model \
        --depth 10 \
        --genome ${fasta} \
        --prefix ${acc} \
        --id-prefix ${acc}__ \
        --length-mean 5000 \
        --length-max 50000 \
        --accuracy-mean 0.98; \
    cat ${acc}*.fastq | pigz > ${acc}.fastq.gz 
done;
```
2. Short Read
   - Read length: 2x150bp (paired)
   - Depth: 10x
   - Random read probability (-y): 0
   - Error rate (-e and -E): 0.01
   - Mutation rate (-r): 0.0 (of which low frequency somatic mutations (-F): 0.0)
   - Random seed (-z): 1

Command used:
```
for fasta in bacteria/*.fa; do \
    acc=$(basename "$fasta" .fa) \
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
        -E 0.01 "$fasta" "$acc";
done;
```

### Expected outputs:
All Viral Genomes
- NCBI_refseq_viruses_20260731.fasta
Viral FASTA References
- NCBI_refseq_viruses_20260731.noTMV.fasta
   - Remove GCA_000854365.1
- NCBI_refseq_viruses_20260731.noHazara.fasta
   - GCA_002831085.1
- NCBI_refseq_viruses_20260731.noLambda.fasta
   - GCA_000840245.1 (E.coli phage lambda)
   - Note: GCA_000840825.1 for Lambdapapillomavirus 2
- NCBI_refseq_viruses_20260731.noT1.fasta
   - GCA_000845005.1 (Escherichia phage T1)

Bacterial FASTA References
All bacterial genomes
- FDA-ARGOS_bacteria_20260731.fasta

Fungal FASTA References
- NCBI_refseq_fungi_20260731.fasta

### Set-up
1. Generate synthetic reads per FASTA
2. Combine Viral and Bacterial synth reads for background datasets 
   - i.e. not including references
3. Test for presence of reference
   - Filter reads mapping to reference if required (note # of reads)
4. Spike in reference reads
5. Test de-hosting process

### Starting with long read content only
- Transferred viral and bacterial datasets across to HPC working directory
- Uncompress directories for bacteria and viruses. Note: Leaving fungi for the moment.
```
mkdir FDA-ARGOS_bacterial_20260731/
unzip FDA-ARGOS_bacterial_20260731.zip -d FDA-ARGOS_bacterial_20260731/
mkdir NCBI_refseq_viruses_20260731
unzip NCBI_refseq_viruses_20260731.zip -d NCBI_refseq_viruses_20260731
```

- Generate synthetic long reads for viruses:
   - Create PBSim3 environment:
   ````
   mm create -n pbsim3
   mm activate  pbsim3
   mm install -c bioconda pbsim3
   ````
```
mkdir -p synth-reads/viral
cd $_


for fasta in /home/phe.gov.uk/nicholas.ellaby/scratch/NCBI_refseq_viruses_20260731/ncbi_dataset/data/*/*.fna; do
    acc=$(basename "$fasta" .fa)
    pbsim --seed 1 --strategy wgs --method errhmm --errhmm /home/phe.gov.uk/nicholas.ellaby/micromamba/envs/pbsim3/data/ERRHMM-ONT-HQ.model --depth 10 --genome ${fasta} --prefix ${acc} --id-prefix ${acc}__ --length-mean 1000 --length-max 10000 --accuracy-mean 0.98; cat ${acc}*.fastq | pigz > ${acc}.fastq.gz
done
```
- Generate synthetic long reads for bacteria

```
for fasta in /home/phe.gov.uk/nicholas.ellaby/scratch/FDA-ARGOS_bacterial_20260731/data/*/*fna; do
    acc=$(basename "$fasta" .fa)
    pbsim --seed 1 --strategy wgs --method errhmm --errhmm /home/phe.gov.uk/nicholas.ellaby/micromamba/envs/pbsim3/data/ERRHMM-ONT-HQ.model --depth 10 --genome ${fasta} --prefix ${acc} --id-prefix ${acc}__ --length-mean 5000 --length-max 50000 --accuracy-mean 0.98; cat ${acc}*.fastq | pigz > ${acc}.fastq.gz;
done;
```

- Minor issue, output files are in .fq.gz, I wanted to consolidate them into single files:
```
for fasta in /home/phe.gov.uk/nicholas.ellaby/scratch/FDA-ARGOS_bacterial_20260731/data/*/*fna; do     acc=$(basename "$fasta" .fa); echo ${acc}; cat ${acc}_*.fq.gz > ${acc}.combined.fq.gz; done;
```

- Do the same for viruses:
```
for fasta in /home/phe.gov.uk/nicholas.ellaby/scratch/NCBI_refseq_viruses_20260731/ncbi_dataset/data/*/*.fna; do     acc=$(basename "$fasta" .fa); echo ${acc}; cat ${acc}_*.fq.gz > ${acc}.combined.fq.gz; done;
```

- Move reference fastq's to a separate directory:
```
mkdir reference_synth_reads
mv viral/GCF_000854365.1_ViralProj15071_genomic.fna.combined.fq.gz reference_synth_reads/
```
- Combine all viral and bacterial datasets, without TMV
```
cat viral/*.fna.combined.fq.gz bacteria/*.fna.combined.fq.gz >>All_viral_bacterial.minusTMV.combined.fq.gz
```
- Copy TMV back into original directory
```
cp reference_synth_reads/GCF_000854365.1_ViralProj15071_genomic.fna.combined.fq.gz  viral/
```
- Some minor clean-up of temp files

### Running tools
#### 1. Deacon
- Running on TMV:
   - Build index: `deacon index build ../NCBI_refseq_viruses_20260731/ncbi_dataset/data/GCF_000854365.1/GCF_000854365.1_ViralProj15071_genomic.fna >GCF_000854365.deacon.idx`
   - Check how many reads it finds associated with TMV in each fa:  
      - Running against dataset without reference:
        ```
        deacon filter -d GCF_000854365.deacon.idx ../synth-reads/All_viral_bacterial.minusTMV.combined.fq.gz -o All_viral_bacterial.minusTMV.combined.deacon.TMV_filt.fq.gz -s All_viral_bacterial.minusTMV.combined.deacon.TMV_filt.summary.json

        Deacon v0.15.0; mode: deplete; input: single; options: abs_threshold=2, rel_threshold=0.01, threads=8(4f+4c)
        Loaded index (k=31, w=15) in 4.35ms
        Retained 24835274/24835287 sequences (100.000%), 127882883577/127882914955 bp (100.000%) in 859.51s. 28895 seqs/s (148.8 Mbp/s)

        "seqs_removed": 13,
        "bp_removed": 31378,
        ```
        To work out which sequences have been removed:
        ```
        # Extract sorted, unique read IDs from each file
        seqkit seq -n ../synth-reads/All_viral_bacterial.minusTMV.combined.fq.gz | sort -u > All_viral_bacterial.minusTMV.combined.ids.txt
        seqkit seq -n All_viral_bacterial.minusTMV.combined.deacon.TMV_filt.fq.gz | sort -u > All_viral_bacterial.minusTMV.combined.deacon.TMV_filt.ids.txt

        # Reads present in All_viral_bacterial.minusTMV.combined but not All_viral_bacterial.minusTMV.combined.deacon.TMV_filt
        comm -23 All_viral_bacterial.minusTMV.combined.ids.txt All_viral_bacterial.minusTMV.combined.deacon.TMV_filt.ids.txt > ids_only_in_All_viral_bacterial.minusTMV.combined.ids.txt

        # Reads present in All_viral_bacterial.minusTMV.combined.deacon.TMV_filt but not All_viral_bacterial.minusTMV.combined
        comm -13 All_viral_bacterial.minusTMV.combined.ids.txt All_viral_bacterial.minusTMV.combined.deacon.TMV_filt.ids.txt > only_in_All_viral_bacterial.minusTMV.combined.deacon.TMV_filt.ids.txt

        ```
        - Sequences present in `All_viral_bacterial.minusTMV.combined` but not in `All_viral_bacterial.minusTMV.combined.deacon.TMV_filt` (no sequences were unique in the opposite direction):
        - So two references generated synthetic reads:  GCF_000870525.1 ([Rehmannia mosaic virus](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000870525.1/), 12 seqs) and GCF_000911995.1([Tomato mottle mosaic virus](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000911995.1/), 1 seq)



      - Running against dataset with reference:  
        ```
        deacon filter -d GCF_000854365.deacon.idx ../synth-reads/All_viral_bacterial.combined.fq.gz -o All_viral_bacterial.combined.deacon.TMV_filt.fq.gz -s All_viral_bacterial.combined.deacon.TMV_filt.summary.json
        

        Deacon v0.15.0; mode: deplete; input: single; options: abs_threshold=2, rel_threshold=0.01, threads=8(4f+4c)
        Loaded index (k=31, w=15) in 11.63ms
        Retained 24835274/24835326 sequences (100.000%), 127882883577/127882978905 bp (100.000%) in 947.57s. 26210 seqs/s (135.0 Mbp/s)
        Filter summary saved to "All_viral_bacterial.combined.deacon.TMV_filt.summary.json"

        "seqs_removed": 52
        "bp_removed": 95328
        ```

   - 13 reads (31,378 bp) removed from dataset without TMV spiked in (false positives)
   - 52 reads (95,328 bp) removed from dataset with TMV spiked in

   - Number of synthetic reads generated for TMV: 39 (63,950 bp):  
   `seqkit stats ../synth-reads/viral/GCF_000854365.1_ViralProj15071_genomic.fna.combined.fq.gz`
       - This suggests that the synthetic reads generated for TMV have been successfully removed from the full dataset (read # from dataset - TMV) - (read # from dataset + TMV)
       - I suspect that the 13 reads present in the total dataset are TMV, but have been included in a reference genome. Will review IDs that match to these sequences, understand which samples they are retrieved from.

   - Validate that those reads are infact from the TMV syntetic reference reads:
      - Get read ids from TMV spiked datasets filtered by deacon:  
      ```
      # Extract sorted, unique read IDs from each file
      seqkit seq -n ../synth-reads/All_viral_bacterial.combined.fq.gz | sort -u > All_viral_bacterial.combined.ids.txt
      
      seqkit seq -n All_viral_bacterial.combined.deacon.TMV_filt.fq.gz | sort -u  > All_viral_bacterial.combined.deacon.TMV_filt.ids.txt

      # Reads present in `All_viral_bacterial.combined.ids.txt` but not in `All_viral_bacterial.combined.deacon.TMV_filt.ids.txt`
      comm -23 All_viral_bacterial.combined.ids.txt All_viral_bacterial.combined.deacon.TMV_filt.ids.txt > ids_only_in_All_viral_bacterial.combined.ids.txt

      comm -13 All_viral_bacterial.combined.ids.txt All_viral_bacterial.combined.deacon.TMV_filt.ids.txt > ids_only_All_viral_bacterial.combined.deacon.TMV_filt.ids.txt
      ```
      - Of the 52 unique reads filtered by Deacon for TMV from `All_viral_bacterial.combined.ids.txt`:
         - GCF_000854365.1: 39 ([Tobacco mosaic virus](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000854365.1/))
         - GCF_000870525.1: 12 ([Rehmannia mosaic virus](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000870525.1/))
         - GCF_000911995.1: 1 ([Tomato mottle mosaic virus](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000911995.1/))

      - With respect to GCF_000870525.1, there are a total of 39 sequences (63,950 bases). Therefore Deacon mistakenly filtered out ~31% of Rehmannia mosaic virus synthetic reads. 
      - With respect to GCF_000911995.1, there are a total of 36 sequences (63,980 bases). Therefore Deacon mistakenly filtered out ~2.8% of Tomato mottle mosaic virus reads.

      - What is Rehmannia mosaic virus:
         - Rehmannia mosaic virus (ReMV) is a plant-infecting virus
      
      - Potential for Further investigation:
         - Map reads to TMV reference, are they exact matches? Identical regions of genomes?
         - If they are exact matches this would limit the utility of those reads?

         - Extrac mis-assigned reads from dataset:
         `seqkit grep -n -f misannotated_tmv_reads.txt All_viral_bacterial.minusTMV.combined -o misannotated_tmv_reads.fq.gz`
         - Map to reference:
         - Remove shared minimizers with:
             - GCF_000854365:TMV
             - GCF_000870525.1: Rehmannia mosaic virus
             - GCF_000911995.1: Tomato mottle mosaic virus
         `deacon index diff 1.idx 2.idx > 1-2.idx`  
         `deacon index build ../NCBI_refseq_viruses_20260731/ncbi_dataset/data/GCF_000870525.1/GCF_000870525.1_ViralProj18885_genomic.fna > GCF_000870525.1.deacon.idx` - Rehmannia mosaic virus   
         `deacon index build ../NCBI_refseq_viruses_20260731/ncbi_dataset/data/GCF_000911995.1/GCF_000911995.1_ViralProj217881_genomic.fna > GCF_000911995.1.deacon.idx` - Tomato mottle mosaic virus
         - Identify and remove shared minimizers between GCF_000854365 (TMV) and GCF_000870525.1 (Rehmannia mosaic virus)
         `deacon index diff GCF_000854365.deacon.idx GCF_000870525.1.deacon.idx> GCF_000854365-GCF_000870525.1.idx`
         ```
         First index: loaded 815 minimizers
         Second index: loaded 793 minimizers
         Removed 13 minimizers, 802 remaining
         Completed diff operation in 4.65ms
         ```
         - Identify and remove shared minimizers between GCF_000854365 (TMV) and GCF_000911995.1 (Tomato mottle mosaic virus)
         `deacon index diff GCF_000854365.deacon.idx GCF_000911995.1.deacon.idx > GCF_000854365-GCF_000911995.1.deacon.idx`   
         ```
         First index: loaded 815 minimizers
         Second index: loaded 815 minimizers
         Removed 3 minimizers, 812 remaining
         Completed diff operation in 1.87ms
         ```
         - Remove common index from all three?
         `deacon index diff GCF_000854365-GCF_000870525.1.idx GCF_000911995.1.deacon.idx > GCF_000854365-GCF_000870525.1-GCF_000911995.1.idx`
         ```
         First index: loaded 802 minimizers
         Second index: loaded 815 minimizers
         Removed 2 minimizers, 800 remaining
         Completed diff operation in 2.37ms
         ```
         - Re-run deacon with new combined index (just on full dataset):   
         `deacon filter -d GCF_000854365-GCF_000870525.1-GCF_000911995.1.idx ../synth-reads/All_viral_bacterial.combined.fq.gz -o All_viral_bacterial.combined.deacon.TMV_filt.combined_index.fq.gz -s All_viral_bacterial.combined.deacon.TMV_filt.combined_index.summary.json`
         - 39 sequences removed, matches the number of TMV reads synthesized and added to the dataset, need to check though

- Running with Hazara Virus (GCA_002831085.1):
   - Move Hazara Virus (GCA_002831085.1) into reference folder:  
   `mv ../synth-reads/viral/GCF_002831085.1_ASM283108v1_genomic.fna.combined.fq.gz ../synth-reads/reference_synth_reads/`
   - Combine all bacterial and viral (minus GCA_002831085.1):
   `cat viral/*.fna.combined.fq.gz bacteria/*.fna.combined.fq.gz >>All_viral_bacterial.minusHazara.combined.fq.gz`
   - Generate Deacon reference:  
   `deacon index build ../NCBI_refseq_viruses_20260731/ncbi_dataset/data/GCF_002831085.1/GCF_002831085.1_ASM283108v1_genomic.fna >GCF_002831085.1.deacon.idx` 
   - Run Deacon with reference on completed dataset:
   `deacon filter -d GCF_002831085.1.deacon.idx ../synth-reads/All_viral_bacterial.combined.fq.gz -o All_viral_bacterial.combined.deacon.Hazara_filt.fq.gz -s All_viral_bacterial.combined.deacon.Hazara_filt.summary.json`
   - Run Deacon with reference on dataset minus GCA_002831085.1:  
   `deacon filter -d GCF_002831085.1.deacon.idx ../synth-reads/All_viral_bacterial.minusHazara.combined.fq.gz -o All_viral_bacterial.minusHazara.combined.deacon.Hazara_filt.fq.gz -s All_viral_bacterial.minusHazara.combined.deacon.Hazara_filt.summary.json`
   - Get Seqs IDs for `All_viral_bacterial.combined.deacon.Hazara_filt.fq.gz`:   
   `seqkit seq -n All_viral_bacterial.combined.deacon.Hazara_filt.fq.gz | sort -u > All_viral_bacterial.combined.deacon.Hazara_filt.ids.txt`
   - Get Seqs IDs for `All_viral_bacterial.minusHazara.combined.fq.gz`:
   `seqkit seq -n ../synth-reads/All_viral_bacterial.minusHazara.combined.fq.gz | sort -u > All_viral_bacterial.minusHazara.combined.ids.txt`
   - Get Seqs IDs for `All_viral_bacterial.minusHazara.combined.deacon.Hazara_filt.fq.gz`:   
   `seqkit seq -n All_viral_bacterial.minusHazara.combined.deacon.Hazara_filt.fq.gz | sort -u > All_viral_bacterial.minusHazara.combined.deacon.Hazara_filt.ids.txt`
   - What reads were present in All bacterial viral not present in the Hazara filtered  
   `comm -23 All_viral_bacterial.combined.ids.txt All_viral_bacterial.combined.deacon.Hazara_filt.ids.txt > only_in_All_viral_bacterial.combined.deacon.Hazara_filt.ids.txt`
      - 109 GCF_002831085.1 (Hazara virus) reads
      - No reads were only in All_viral_bacterial.combined.deacon.Hazara_filt.ids.txt (sanity check)
   - What reads were present in All bacterial minus Hazara vs Hazara filtered:   
   `comm -23 All_viral_bacterial.minusHazara.combined.ids.txt All_viral_bacterial.minusHazara.combined.deacon.Hazara_filt.combined.ids.txt > only_in_All_viral_bacterial.minusHazara.combined.ids.txt`
      - No reads present only in All_viral_bacterial.minusHazara vs All_viral_bacterial.minusHazara.combined.deacon.Hazara_filt
      - No reads present only in All_viral_bacterial.minusHazara.combined.deacon.Hazara_filt vs All_viral_bacterial.minusHazara (sanity check)
   
- Running with Lambda Phage (NC_001416.1): 
   - Setting up data, AK sent all sequences combined into single fastas. Seperate them out and edit names:   
   `seqkit split -i 2731619_genbank_qc_phage.fasta --out-dir individual_seqs`
   `cd individual_seqs`   
   `for f in *2731619_genbank_qc_phage.part_*; do mv "${f}" "${f//2731619_genbank_qc_phage.part_/}"; done`   
   - Generate synthetic data for each sequence
   - Looking at the sequences 6,554 were <3Kb in size. Removed from dataset
   `mkdir small_seqs`
   `while read f; do mv individual_seqs/${f} small_seqs/; done<sequences_under_3kb.txt`
   - Generate syntetic data for phage now:
      ```
      for fasta in ../../phage/individual_seqs/*.fasta; do \
      acc=$(basename "$fasta" .fasta) \
      pbsim --seed 1 --strategy wgs --method errhmm --errhmm /home/phe.gov.uk/nicholas.ellaby/micromamba/envs/pbsim3/data/ERRHMM-ONT-HQ.model --depth 10 --genome ${fasta} --prefix ${acc} --id-prefix ${acc}__ --accuracy-mean 0.98; cat ${acc}*.fq.gz | pigz > ${acc}.fq.gz
      done;
      ```
   - Combine phages, virus and bacterial into a single fq:   
      `cat phages/*.fq.gz All_viral_bacterial.combined.fq.gz >All_viral_bacterial_phage.combined.fq.gz`
   - Run deacon with phage references:
      - Lambda: GCF_000840245.1, not present in original dataset. Download and generate synth reads   
         `pbsim --seed 1 --strategy wgs --method errhmm --errhmm /home/phe.gov.uk/nicholas.ellaby/micromamba/envs/pbsim3/data/ERRHMM-ONT-HQ.model --depth 10 --genome ../phage/individual_seqs/GCF_000840245.1_ViralProj14204_genomic.fna --prefix GCF_000840245.1 --id-prefix GCF_000840245.1__ --accuracy-mean 0.98`   
         - Combine with total dataset:
         `cat ../synth-reads/All_viral_bacterial_phage.minus_lambda.combined.fq.gz ../synth-reads/phages/GCF_000840245.1_0001.fq.gz >>../synth-reads/All_viral_bacterial_phage.combined.fq.gz`
         - Run Deacon to create index:   
         `deacon index build ../phage/individual_seqs/GCF_000840245.1_ViralProj14204_genomic.fna >GCF_000840245.1_ViralProj14204.deacon.idx`
         - Run Deacon to remove reads from `All_viral_bacterial_phage.combined.fq.gz`:
         `deacon filter -d GCF_000840245.1_ViralProj14204.deacon.idx ../synth-reads/All_viral_bacterial_phage.combined.fq.gz -o All_viral_bacterial_phage.combined.deacon.Lambda_filt.fq.gz -s All_viral_bacterial_phage.combined.deacon.Lambda_filt.summary.json`
            Deacon v0.15.0; mode: search; input: single; options: abs_threshold=2, rel_threshold=0.01, threads=8(4f+4c) Loaded index (k=31, w=15)
            Retained 25827000 / 25868704 (99.83%), 138585775731 / 138916500929 (99.76%)


         - Run Deacon to retain reads from `All_viral_bacterial_phage.combined.fq.gz`:   
         `deacon filter GCF_000840245.1_ViralProj14204.deacon.idx ../synth-reads/All_viral_bacterial_phage.combined.fq.gz -o All_viral_bacterial_phage.combined.deacon.Lambda.fq.gz -s All_viral_bacterial_phage.combined.deacon.Lambda.summary.json`  
            Deacon v0.15.0; mode: search; input: single; options: abs_threshold=2, rel_threshold=0.01, threads=8(4f+4c) Loaded index (k=31, w=15) in 5.32ms  
            Retained 41704/25868704 sequences (0.161%), 330725198/138916500929 bp (0.238%) in 366.72s. 70543 seqs/s (378.8 Mbp/s)

- Running with phage MS2 (NC_001417.2):
   
#### 2. Hostile
- Install:  
`mm create -y -n hostile -c conda-forge -c bioconda hostile`  
`mm activate hostile`
- Generate index, minimapper used for long-reads:   
`minimap2 -d GCF_000854365.1_ViralProj15071_genomic.mni ../NCBI_refseq_viruses_20260731/ncbi_dataset/data/GCF_000854365.1/GCF_000854365.1_ViralProj15071_genomic.fna`

- Run clean using index with Dataset that doesn't contain TMV:  
`hostile clean --fastq1 ../synth-reads/All_viral_bacterial.minusTMV.combined.fq.gz --index GCF_000854365.1_ViralProj15071_genomic.mni -t 12 -o All_viral_bacterial.minusTMV.combined.hostile.TMV_filt.fq.gz`
    - Results: 153 reads removed
    - Get sequence IDs for those reads that have passed filtering (requires seqkit install btw):  
    `seqkit seq -n All_viral_bacterial.minusTMV.combined.hostile.TMV_filt.fq.gz/All_viral_bacterial.minusTMV.combined.clean.fastq.gz | sort -u > All_viral_bacterial.minusTMV.combined.hostile.TMV_filt.ids.txt`


- Run clean using index with dataset that contains TMV:  
`hostile clean --fastq1 ../synth-reads/All_viral_bacterial.combine
d.fq.gz --index GCF_000854365.1_ViralProj15071_genomic.mni -t 12 -o All_viral_bacterial.combined.hostile.TMV_filt.fq.gz`
    - Results: 192 reads removed (+39 reads)
    - Get sequence IDs for those reads that have passed filtering:
    `seqkit seq -n All_viral_bacterial.combined.hostile.TMV_filt.fq.gz/All_viral_bacterial.combined.clean.fastq.gz | sort -u > All_viral_bacterial.combined.hostile.TMV_filt.ids.txt`
    `seqkit seq -n All_viral_bacterial.minusTMV.combined.hostile.TMV_filt.fq.gz/All_viral_bacterial.minusTMV.combined.hostile.TMV_filt.fq.gz | sort -u >All_viral_bacterial.minusTMV.combined.hostile.TMV_filt.ids.txt`
    - What sequence ids are present only in All_viral_bacterial.combined.ids.txt vs All_viral_bacterial.combined.hostile.TMV_filt.ids.txt
    `comm -23 ../deacon/All_viral_bacterial.combined.ids.txt All_viral_bacterial.combined.hostile.TMV_filt.ids.txt > only_../deacon/All_viral_bacterial.combined.ids.txt`
       - Hostile returned many Taxa: 39 (reads) GCF_000854365.1 (TMV); 23 GCF_001461485.1(Tomato brown rugose fruit); 23 GCF_000870525.1 (Rehmannia mosaic); 21 GCF_000853705.1 (Tomato mosaic); 17 GCF_000911995.1 (Tomato mottle mosaic); 17 GCF_000873425.1 (Bell pepper mottle); 16 GCF_000859645.1 (Pepper mild mottle); 9 GCF_000879495.1 (Brugmansia mild mottle); 7 GCF_001654245.1 (Tropical soda apple mosaic); 5 GCF_000847545.1 (Tobacco mild green mosaic); 4 GCF_002145505.1 (Hoya chlorotic spot); 4 GCF_000912435.1 (Yellow tailflower mild mottle); 3 GCF_000853745.1 (Paprika mild mottle); 2 GCF_000859905.1 (Odontoglossum ringspot); 1 GCF_000870125.1 (Streptocarpus flower break); 1 GCF_000852265.1 (Obuda pepper)
          - All plant viruses

- NB: There is an option to invert the Hostile process, keeping the mapped ("filtered") reads, rather than the non-mapped reads. However, given this isn't the proposed use case, will run as standard. I don't expect there to be a difference between inverted reads and those that are removed in the standard clean process.

#### 2. Detaxizer
- This is an nf-core nextflow process, requires the installation of nextflow
- Install:    
`mm create -n detaxizser`  
`mm activate detaxizer`  
`mm install -c bioconda nextflow`  
   - Note: It is not recommended to install nextflow in this manner, but makes sense on this HPC
- Required reducing strictness of syntax: `export NXF_SYNTAX_PARSER=v1`
- Test install:   
`nextflow run nf-core/detaxizer -profile test,apptainer --outdir test_profile`
- Issues running on SMED



