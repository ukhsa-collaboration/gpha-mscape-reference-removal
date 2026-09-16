"""
Unit tests for bin/generate_validation_report.py

Run with: pytest tests/bin/test_generate_validation_report.py -v

These test the pure logic (ID trimming, set-based confusion-matrix math)
without needing deacon, Nextflow, or real sequencing data -- small
in-memory FASTQ fixtures are enough to exercise every code path.
"""

import json
import sys
from pathlib import Path

import pytest

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN_DIR))

import generate_validation_report as gvr

# --------------------------------------------------------------------------
# trim_read_id
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("read1/1", "read1"),
        ("read1/2", "read1"),
        ("read1", "read1"),
        ("read1/12", "read1/12"),  # doesn't end in exactly "/1" or "/2" -- untouched
    ],
)
def test_trim_read_id(raw, expected):
    assert gvr.trim_read_id(raw) == expected


# --------------------------------------------------------------------------
# get_ids -- needs real FASTQ files since it calls pyfastx
# --------------------------------------------------------------------------


def _write_fastq(path, read_ids):
    with open(path, "w") as fh:
        fh.writelines(f"@{rid}\nACGT\n+\nIIII\n" for rid in read_ids)


def test_get_ids_pools_across_multiple_files_and_trims_mates(tmp_path):
    f1 = tmp_path / "a.fastq"
    f2 = tmp_path / "b.fastq"
    _write_fastq(f1, ["r1/1", "r2/1"])
    _write_fastq(f2, ["r1/2", "r3"])
    ids = gvr.get_ids([str(f1), str(f2)])
    # r1/1 and r1/2 collapse to the same trimmed id "r1"
    assert ids == {"r1", "r2", "r3"}


# --------------------------------------------------------------------------
# Confusion-matrix logic, exercised via a small end-to-end fixture set.
# This is a regression test for the exact TP/FP/FN/TN definitions --
# if someone flips a set operation, this catches it immediately.
# --------------------------------------------------------------------------


@pytest.fixture
def fixture_dir(tmp_path):
    """
    Ground truth:
      reference reads: ref1, ref2, ref3      (3 total)
      background reads: bg1, bg2, bg3, bg4    (4 total)

    Deacon's isolate output (predicted positive) contains: ref1, ref2, bg1
      -> TP = {ref1, ref2}          (2)
      -> FN = {ref3}                (1, escaped depletion)
      -> FP = {bg1}                 (1, background wrongly caught)
    Deacon's depleted output (kept) contains: ref3, bg2, bg3, bg4
      -> TN = background_truth - predicted_pos = {bg2, bg3, bg4}  (3)
    """
    ref = tmp_path / "reference.fastq"
    bg = tmp_path / "background.fastq"
    iso = tmp_path / "isolate.fastq"
    dep = tmp_path / "depleted.fastq"
    deacon_json = tmp_path / "deacon.json"

    _write_fastq(ref, ["ref1", "ref2", "ref3"])
    _write_fastq(bg, ["bg1", "bg2", "bg3", "bg4"])
    _write_fastq(iso, ["ref1", "ref2", "bg1"])
    _write_fastq(dep, ["ref3", "bg2", "bg3", "bg4"])

    deacon_json.write_text(json.dumps({"seqs_in": 7, "seqs_out": 3}))

    return {
        "reference": str(ref),
        "background": str(bg),
        "isolate": str(iso),
        "depleted": str(dep),
        "deacon_json": str(deacon_json),
        "outprefix": str(tmp_path / "sample1"),
    }


def _run_main(args, fixture_dir):
    argv = [
        "generate_validation_report.py",
        "--sample-id",
        "sample1",
        "--read-type",
        "single",
        "--reference-fastq",
        fixture_dir["reference"],
        "--background-fastq",
        fixture_dir["background"],
        "--isolate-fastq",
        fixture_dir["isolate"],
        "--depleted-fastq",
        fixture_dir["depleted"],
        "--deacon-json",
        fixture_dir["deacon_json"],
        "-o",
        fixture_dir["outprefix"],
    ] + args
    old_argv = sys.argv
    sys.argv = argv
    try:
        gvr.main()
    finally:
        sys.argv = old_argv


def test_confusion_matrix_counts_are_correct(fixture_dir):
    _run_main([], fixture_dir)
    result = json.loads(Path(fixture_dir["outprefix"] + ".confusion.json").read_text())

    assert result["tp"] == 2
    assert result["fp"] == 1
    assert result["fn"] == 1
    assert result["tn"] == 3
    assert result["ref_truth_total"] == 3
    assert result["background_truth_total"] == 4
    # seqs_in (7) == ref_truth_total (3) + background_truth_total (4)
    assert result["count_consistency_ok"] is True
    assert result["untraceable_reads"] == 0

    assert result["precision"] == pytest.approx(2 / 3)  # TP / (TP+FP)
    assert result["recall"] == pytest.approx(2 / 3)  # TP / (TP+FN)
    assert result["specificity"] == pytest.approx(3 / 4)  # TN / (TN+FP)


def test_count_consistency_flag_catches_mismatched_deacon_json(fixture_dir):
    # Deliberately corrupt the deacon stats so seqs_in no longer matches
    # ref_truth_total + background_truth_total -- this is the exact
    # caveat that would otherwise let a mismatched idx/reads pairing
    # slip through silently.
    Path(fixture_dir["deacon_json"]).write_text(
        json.dumps({"seqs_in": 999, "seqs_out": 3})
    )
    _run_main([], fixture_dir)
    result = json.loads(Path(fixture_dir["outprefix"] + ".confusion.json").read_text())
    assert result["count_consistency_ok"] is False


def test_untraceable_reads_detected(tmp_path):
    """A read in the isolate output that belongs to neither truth set
    signals an ID-normalisation mismatch (e.g. unstripped mate suffixes)."""
    ref = tmp_path / "reference.fastq"
    bg = tmp_path / "background.fastq"
    iso = tmp_path / "isolate.fastq"
    dep = tmp_path / "depleted.fastq"
    deacon_json = tmp_path / "deacon.json"

    _write_fastq(ref, ["ref1"])
    _write_fastq(bg, ["bg1"])
    _write_fastq(iso, ["ref1", "mystery_read"])
    _write_fastq(dep, ["bg1"])
    deacon_json.write_text(json.dumps({"seqs_in": 2, "seqs_out": 2}))

    fixture = {
        "reference": str(ref),
        "background": str(bg),
        "isolate": str(iso),
        "depleted": str(dep),
        "deacon_json": str(deacon_json),
        "outprefix": str(tmp_path / "sample2"),
    }
    _run_main([], fixture)
    result = json.loads(Path(fixture["outprefix"] + ".confusion.json").read_text())
    assert result["untraceable_reads"] == 1


# --------------------------------------------------------------------------
# CLI contract tests -- these are the ones that catch vocabulary drift
# between this script and the rest of the pipeline (e.g. the
# long/short -> single/paired rename that never reached this file).
# --------------------------------------------------------------------------


@pytest.mark.parametrize("read_type", ["single", "paired"])
def test_read_type_accepts_current_pipeline_vocabulary(fixture_dir, read_type):
    """
    Regression test: subworkflow/report.nf calls this script with
    --read-type set to 'single' or 'paired' (see subworkflow/report.nf,
    'single' or 'paired' comment on the read_type take: parameter).
    If this script's argparse --choices haven't been updated to match,
    this test fails with SystemExit(2) from argparse.
    """
    try:
        _run_main(["--read-type", read_type], fixture_dir)
    except SystemExit as e:
        pytest.fail(
            f"generate_validation_report.py rejected --read-type {read_type!r} "
            f"(exit code {e.code}) -- its argparse --choices are out of sync "
            "with the 'single'/'paired' vocabulary used elsewhere in the pipeline."
        )
