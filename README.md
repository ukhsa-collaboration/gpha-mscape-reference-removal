# Reference Removal
 
A Nextflow (DSL2) pipeline for removing reference/spike-in sequences from FASTQ data using [Deacon](https://github.com/bede/deacon), with a companion **validation mode** that measures how accurately a given reference is actually removed before you trust the pipeline against real samples.
 
Built for internal-control depletion in metagenomic sequencing — e.g. confirming that a spiked-in control (Tobacco Mosaic Virus, Lambda phage, ERCC, etc.) is cleanly separated from the biological signal in a sample.

## Contents
 
- [How it works](#how-it-works)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Reference removal](#reference-removal)
- [Reference validation](#reference-validation)
- [Input formats](#input-formats)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [Testing](#testing)
- [Background & research notes](#background--research-notes)

## How it works
 
The pipeline has two modes, selected with `--validation true|false`:
 
- **Reference removal**: the production path. Give it a reference (FASTA or a prebuilt Deacon index) and some samples; it runs `deacon   filter` and hands back the depleted reads, the matched (removed) reads, and per-sample summary stats.  
- **Reference validation**: a self-test. Given only a reference FASTA, it generates synthetic reads *from that reference* (via [pbsim3](https://github.com/yukiteruono/pbsim3) for long reads and [dwgsim](https://github.com/nh13/DWGSIM) for short reads), spikes them into your background sample(s), runs the same `deacon filter` step, and then — because the identity of every spiked-in read is known — computes a proper confusion matrix (true/false positive/negative, recall, precision, specificity) and renders it as an HTML + PDF report. This is how you'd confirm a reference or index actually depletes cleanly before using it for real.  

Both modes accept either **paired-end (short-read)** or **single-end (long-read)** samples, auto-detected — you don't need to tell the pipeline which is which.  

## Requirements
 
- [Nextflow](https://www.nextflow.io/) (DSL2)
- One of: Docker, Singularity, or Conda/Mamba (all containers/environments
  are pinned per-process — no manual tool installation needed)  

## Quick start
 
**Remove a reference from some samples:**
```bash
nextflow run nfellaby/reference-removal -profile docker \
  --validation false \
  --fasta reference.fna \
  --samplesheet samples.csv \
  --outdir results/
```
 
**Validate how well a reference depletes:**
```bash
nextflow run nfellaby/reference-removal -profile docker \
  --validation true \
  --fasta reference.fna \
  --read_type both \
  --sample_data_dir background_fastqs/ \
  --outdir results/
```

## Reference removal
 
`--validation false` (the default).
 
**Reference**: supply exactly one of `--fasta` or `--idx`. If you give a FASTA, a Deacon index is built for you (`deacon index build`); if you already have a prebuilt index, pass it with `--idx` to skip that step.
 
**Samples**: supply exactly one of `--samplesheet` or `--sample_data_dir`. Every sample,  single-end or paired-end, is auto-detected and processed; `--read_type` has no effect in this mode (single- and paired-end sample found in your input are both handled regardless of what it's set to).
 
**What runs**: for every sample, `deacon filter` is run against the index twice, once to produce the *depleted* reads (reference removed, what you'd keep), and once to isolate the *matched* reads (what was removed) so you get both the clean output and a record of exactly what was taken out.

## Reference validation
 
`--validation true`.
 
**Reference**: `--fasta` is always required (synthetic reads have to be generated from real sequence). `--idx` is optional, if supplied it's used directly for the `deacon filter` step (skipping index rebuild); the FASTA is still used to generate the synthetic spike-in reads either way.
 
**Background samples**: supply exactly one of `--samplesheet` or `--sample_data_dir`, this is the "clean" background data the synthetic reference reads get spiked into.
 
**`--read_type`**: `single`, `paired`, or `both` (default). Controls which kind of synthetic reads get generated and spiked in, and which kind of background data is required. With `both`, background data for *each* read type is required by default, if your background data only covers one, either restrict `--read_type` to match, or set `--allow_partial_validation true` to downgrade the missing type to a warning instead of a hard failure.
 
**What runs**: synthetic reference reads → spiked into background sample(s) → `deacon filter` (depleted + matched, same as removal mode) → per-sample confusion matrix (comparing known spiked-read IDs against what `deacon` actually caught) → aggregated HTML + PDF report.

## Input formats
 
### Samplesheet (`--samplesheet`)
 
CSV or TSV (detected by file extension), **first row is a header and is always skipped**, its exact column names don't matter, only position:
 
| Column | Contents |
|---|---|
| 1 | Sample ID |
| 2 | Read 1 FASTQ path |
| 3 | Read 2 FASTQ path — leave blank for single-end samples |
 
Single-end and paired-end rows can be freely mixed in the same samplesheet; each row is classified independently.

### Directory of FASTQs (`--sample_data_dir`)
 
Point at a directory and every `.fastq`/`.fq`/`.fastq.gz`/`.fq.gz` file anywhere in it (including subdirectories) is picked up automatically. Mate pairs are detected by filename, case-insensitively, using any of these conventions immediately before the extension: `_1`/`_2`, `_R1`/`_R2`, `_read1`/`_read2`, `.R1`/`.R2`, with or without a trailing `_001` (standard Illumina naming). Two files that reduce to the same sample ID after stripping the mate marker are treated as one paired-end sample; anything left unmatched is treated as single-end.

## Parameters
 
| Parameter | Default | Description |
|---|---|---|
| `--validation` | `false` | `true` runs the validation workflow, `false` runs plain removal |
| `--fasta` | `null` | Reference FASTA. Always required for validation; one of `--fasta`/`--idx` required for removal |
| `--idx` | `null` | Prebuilt Deacon index. Optional in both modes — used in place of rebuilding from `--fasta` when supplied |
| `--samplesheet` | `null` | Samplesheet path — exactly one of this or `--sample_data_dir` required |
| `--sample_data_dir` | `""` | Directory of FASTQs to auto-discover — exactly one of this or `--samplesheet` required |
| `--read_type` | `both` | `single`, `paired`, or `both`. Validation only — ignored for plain removal |
| `--allow_partial_validation` | `false` | With `--read_type both`, allow validation to proceed with only one read type present (warns instead of failing) |
| `--pbsim_model` | *(site-specific path)* | PBSIM3 error model file, used for long-read synthesis during validation |
| `--outdir` | `output` | Where results are published |
| `--publish_dir_mode` | `copy` | How outputs are published (`copy`, `symlink`, etc. — any valid Nextflow `publishDir` mode) |
| `--max_memory` / `--max_cpus` / `--max_time` | `10.GB` / `4` / `240.h` | Resource ceilings applied across all processes |

## Outputs
 
```
<outdir>/
├── data/
│   └── <sample_id>_single_reads/     # or _paired_reads/
│       ├── <sample_id>.sample_summary.json      # depleted-read stats
│       ├── <sample_id>.sample_reads.fq.gz        # depleted reads (R1/R2 for paired)
│       └── <sample_id>.reference_reads.fq.gz     # matched/removed reads (R1/R2 for paired)
└── report/                            # validation mode only
    ├── single.validation_report.html
    ├── single.validation_report.pdf
    ├── paired.validation_report.html
    └── paired.validation_report.pdf
```

The validation report shows, per sample: a confusion matrix (TP/FP/FN/TN), recall/precision/specificity/F1, and flagged caveats e.g. recall below a configurable threshold, any false positives, or a mismatch between `deacon`'s own read counts and the known truth set (a sign something upstream doesn't match what's actually being compared).

## Testing
 
`tests/` has two layers — see [`tests/README.md`](tests/README.md) for details on running either:
 
- **Python** (`tests/bin/`) — pytest, covers the report-generation logic (confusion-matrix math, caveat triggers, CLI contract) directly, no Nextflow required.
- **Nextflow** ([nf-test](https://www.nf-test.com/)) — module and subworkflow tests covering index resolution (`--fasta` vs `--idx` preference), sample auto-discovery, and output shape correctness.

## Background & research notes
 
[`docs/Workflow.md`](docs/Workflow.md) has the exploratory research behind this pipeline, reference dataset provenance, comparisons against other de-hosting tools (Hostile, NoHuman, detaxizer), and the synthetic-read parameters used during early tool evaluation. It's a lab notebook, not pipeline documentation, this README reflects what the pipeline actually
does today.