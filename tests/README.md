# Test suite

## Python (`tests/bin/`) — logic in `bin/*.py`

Run directly with pytest, no Nextflow/deacon/containers needed:

```bash
pip install pytest pyfastx pandas plotly matplotlib reportlab
pytest tests/bin/ -v
```

These were run against the actual `bin/generate_validation_report.py` and
`bin/build_report.py` while writing this suite. `build_report.py` passed
cleanly. `generate_validation_report.py`

## Nextflow (`tests/modules/`, `tests/subworkflow/`, `tests/main.nf.test`)

Written against [nf-test](https://www.nf-test.com/) syntax and checked
for correctness against the current module/subworkflow signatures.

Run for real with:

```bash
# https://www.nf-test.com/installation/
nf-test test tests/ --tag modules,subworkflow,pipeline
```

What each layer is specifically guarding against:

- **`tests/modules/generate_idx.nf.test`, `filter.nf.test`** — asserts
  the declared `output:` filenames match what the `deacon` invocation
  actually writes. This exact mismatch (e.g. declaring
  `*.sample_summary.json` while the script wrote `*.ref_removed.json`)
  broke every filter process at one point and would only have surfaced
  as a runtime "Missing output file(s)" error without a test like this.
- **`tests/subworkflow/reference_parsing.nf.test`,
  `resolve_reference.nf.test`** — asserts `GENERATE_IDX` does **not**
  run when a prebuilt `--idx` is supplied, and **does** run when it
  isn't. This is the regression test for `idx_fp` being accepted as a
  parameter but silently never read.
- **`tests/subworkflow/samples_parsing.nf.test`** — asserts correct
  single/paired branching from both a samplesheet and a directory of
  FASTQs, and that a malformed samplesheet row fails loudly rather than
  silently dropping data.
- **`tests/main.nf.test`** — pipeline-level parameter validation: the
  fasta/idx and samplesheet/sample_data_dir XOR-style checks in
  `main.nf` are exactly the logic that drifted out of sync with the
  subworkflows underneath more than once during development.

## Fixtures (`tests/test-data/`)

Deliberately tiny (a 300bp repeat sequence, ~2 reads per file) so the
whole suite runs in seconds rather than minutes — correctness of the
logic doesn't depend on realistic sequence content, only on the shapes
(single vs paired, matching vs non-matching read IDs) being right.

## Not yet covered

- End-to-end `workflow REFERENCE_VALIDATION` / `REFERENCE_REMOVAL` tests
  (`tests/workflows/`) — these need a full run through `deacon`, `pbsim`,
  and `dwgsim` in containers, which is why they're not included yet
  given this environment can't verify them. Worth adding once you can
  run them against real containers to confirm expected output shapes.
- A CI workflow (`.github/workflows/test.yml`) wiring both pytest and
  `nf-test` into GitHub Actions on every PR — happy to draft this next.
