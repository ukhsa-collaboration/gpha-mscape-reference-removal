"""
Unit tests for bin/build_report.py

Run with: pytest tests/bin/test_build_report.py -v

Covers the pure logic (caveat triggers, formatting, record loading) plus
a full CLI smoke test. Chart rendering (matplotlib/reportlab) is exercised
via the smoke test rather than pixel-inspected -- correctness of what's
*drawn* is a visual-review concern, not a unit-test one; what this suite
guards is that the report reliably gets built at all, and that the
caveat logic that flags a bad run is correct.
"""

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN_DIR))

import build_report as br


def _confusion_record(**overrides):
    base = {
        "sample_id": "sample1",
        "read_type": "single",
        "tp": 95,
        "fp": 0,
        "fn": 5,
        "tn": 900,
        "precision": 1.0,
        "recall": 0.95,
        "specificity": 1.0,
        "f1": 0.974,
        "untraceable_reads": 0,
        "ref_truth_total": 100,
        "background_truth_total": 900,
        "deacon_seqs_in": 1000,
        "deacon_seqs_out": 95,
        "count_consistency_ok": True,
    }
    base.update(overrides)
    return base


# --------------------------------------------------------------------------
# load_records
# --------------------------------------------------------------------------


def test_load_records_reads_multiple_jsons(tmp_path):
    p1 = tmp_path / "a.confusion.json"
    p2 = tmp_path / "b.confusion.json"
    p1.write_text(json.dumps(_confusion_record(sample_id="a")))
    p2.write_text(json.dumps(_confusion_record(sample_id="b")))
    df = br.load_records([str(p1), str(p2)])
    assert set(df["sample_id"]) == {"a", "b"}


def test_load_records_rejects_missing_field(tmp_path):
    p1 = tmp_path / "bad.confusion.json"
    bad = _confusion_record()
    del bad["recall"]
    p1.write_text(json.dumps(bad))
    with pytest.raises(ValueError, match="missing expected field"):
        br.load_records([str(p1)])


def test_load_records_rejects_empty_list():
    with pytest.raises(ValueError, match="No confusion JSONs"):
        br.load_records([])


# --------------------------------------------------------------------------
# fmt_pct
# --------------------------------------------------------------------------


def test_fmt_pct_formats_and_handles_nan():
    assert br.fmt_pct(0.9541) == "95.41%"
    assert br.fmt_pct(float("nan")) == "n/a"


# --------------------------------------------------------------------------
# caveats_table -- one test per trigger, and one confirming a clean
# sample produces no caveats. This is the part of the report users
# actually rely on to know when NOT to trust the results.
# --------------------------------------------------------------------------


def _df(records):
    return pd.DataFrame(records)


def test_caveats_clean_sample_has_none():
    df = _df([_confusion_record()])
    caveats = br.caveats_table(df)
    assert caveats.iloc[0]["Caveats"] == "None detected"


def test_caveats_flags_count_inconsistency():
    df = _df([_confusion_record(count_consistency_ok=False, deacon_seqs_in=500)])
    caveats = br.caveats_table(df)
    assert "does not match" in caveats.iloc[0]["Caveats"]


def test_caveats_flags_untraceable_reads():
    df = _df([_confusion_record(untraceable_reads=3)])
    caveats = br.caveats_table(df)
    assert "3 read ID(s)" in caveats.iloc[0]["Caveats"]


def test_caveats_flags_low_recall():
    df = _df([_confusion_record(recall=0.80, fn=20, tp=80)])
    caveats = br.caveats_table(df)
    assert "Recall is 80.00%" in caveats.iloc[0]["Caveats"]


def test_caveats_flags_false_positives():
    df = _df([_confusion_record(fp=4)])
    caveats = br.caveats_table(df)
    assert "4 background read(s)" in caveats.iloc[0]["Caveats"]


def test_caveats_recall_threshold_is_configurable():
    """Regression guard: --recall-threshold must actually change the
    trigger point, not just the displayed reference line."""
    df = _df([_confusion_record(recall=0.97)])
    br.RECALL_REFERENCE_LINE = 0.99
    try:
        caveats = br.caveats_table(df)
        assert "Recall is 97.00%" in caveats.iloc[0]["Caveats"]
    finally:
        br.RECALL_REFERENCE_LINE = 0.95  # restore module-level default


# --------------------------------------------------------------------------
# Full CLI smoke test -- both HTML and PDF get written and are non-empty.
# This is exactly the check that would have caught the earlier
# build_html_report.py / build_report.py filename mismatch between the
# script's own name and what modules/report.nf actually calls.
# --------------------------------------------------------------------------


def test_cli_produces_html_and_pdf(tmp_path):
    json_path = tmp_path / "sample1.confusion.json"
    json_path.write_text(json.dumps(_confusion_record()))

    html_out = tmp_path / "report.html"
    result = subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "build_report.py"),
            "--jsons",
            str(json_path),
            "-o",
            str(html_out),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    pdf_out = html_out.with_suffix(".pdf")
    assert html_out.exists() and html_out.stat().st_size > 0
    assert pdf_out.exists() and pdf_out.stat().st_size > 0


def test_cli_respects_recall_threshold_flag(tmp_path):
    json_path = tmp_path / "sample1.confusion.json"
    json_path.write_text(json.dumps(_confusion_record(recall=0.97)))
    html_out = tmp_path / "report.html"

    result = subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "build_report.py"),
            "--jsons",
            str(json_path),
            "-o",
            str(html_out),
            "--recall-threshold",
            "0.99",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    html = html_out.read_text()
    assert "97.00%" in html  # summary table renders the actual value
    assert "below the 99% reference line" in html  # caveat uses the flag's threshold
