#!/usr/bin/env python3
"""generate_validation_report.py — cross-reference read IDs to compute
TP/FP/FN/TN for a spiked reference-removal validation run."""

import argparse
import json

import pyfastx


def trim_read_id(read_id):
    if read_id.endswith(("/1", "/2")):
        read_id = read_id[:-2]
    return read_id


def get_ids(fastq_paths):
    ids = set()
    for path in fastq_paths:
        ids |= {
            trim_read_id(name) for name, *_ in pyfastx.Fastq(path, build_index=False)
        }
    return ids


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample-id", required=True)
    p.add_argument("--read-type", required=True, choices=["single", "paired"])
    p.add_argument(
        "--reference-fastq",
        required=True,
        nargs="+",
        help="pre-spike synthetic reference reads",
    )
    p.add_argument(
        "--background-fastq",
        required=True,
        nargs="+",
        help="pre-spike background reads",
    )
    p.add_argument(
        "--isolate-fastq",
        required=True,
        nargs="+",
        help="deacon SAMPLE_REMOVAL output (matched reads)",
    )
    p.add_argument(
        "--depleted-fastq",
        required=True,
        nargs="+",
        help="deacon REFERENCE_REMOVAL output (kept/background)",
    )
    p.add_argument(
        "--deacon-json",
        required=True,
        help="deacon stats JSON, for a sanity-check cross-reference",
    )
    p.add_argument("-o", "--outprefix", required=True)
    args = p.parse_args()

    ref_truth = get_ids(args.reference_fastq)
    background_truth = get_ids(args.background_fastq)
    predicted_pos = get_ids(args.isolate_fastq)
    kept_ids = get_ids(args.depleted_fastq)  # should approximate TN ∪ FN

    tp = predicted_pos & ref_truth
    fp = predicted_pos - ref_truth
    fn = ref_truth - predicted_pos
    tn = background_truth - predicted_pos

    # Consistency check: flag anything deacon output that isn't traceable
    # to either truth set at all -- signals an ID-normalisation mismatch
    # (e.g. mate-pair suffixing) rather than a real classification result.
    untraceable = predicted_pos - (ref_truth | background_truth)

    precision = len(tp) / len(tp | fp) if (tp or fp) else float("nan")
    recall = (
        len(tp) / len(tp | fn) if (tp or fn) else float("nan")
    )  # sensitivity -- the safety-critical number
    specificity = len(tn) / len(tn | fp) if (tn or fp) else float("nan")
    f1 = (
        (2 * precision * recall / (precision + recall))
        if (precision + recall)
        else float("nan")
    )

    with open(args.deacon_json) as fh:
        deacon_stats = json.load(fh)

    result = {
        "sample_id": args.sample_id,
        "read_type": args.read_type,
        "tp": len(tp),
        "fp": len(fp),
        "fn": len(fn),
        "tn": len(tn),
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "f1": f1,
        "untraceable_reads": len(untraceable),
        "ref_truth_total": len(ref_truth),
        "background_truth_total": len(background_truth),
        "deacon_seqs_in": deacon_stats.get("seqs_in"),
        "deacon_seqs_out": deacon_stats.get("seqs_out"),
        "count_consistency_ok": deacon_stats.get("seqs_in")
        == len(ref_truth) + len(background_truth),
    }

    with open(f"{args.outprefix}.confusion.json", "w") as fh:
        json.dump(result, fh, indent=2)


if __name__ == "__main__":
    main()
