#!/usr/bin/env python3
"""
build_html_report.py

Aggregates per-sample confusion-matrix JSONs (produced by
generate_validation_report.py) into:
  1. a self-contained interactive HTML report (Plotly), and
  2. a static PDF report (matplotlib + reportlab)

...summarising reference-removal validation performance: recall/sensitivity,
precision, specificity, per-sample confusion matrices, and data-integrity
caveats.

The PDF path deliberately does NOT depend on Plotly/Kaleido, which requires
a Chrome install to rasterise figures -- too heavy a dependency for a
lightweight per-tool Nextflow container. matplotlib produces the static
chart images instead; reportlab assembles the PDF.

Usage:
    build_html_report.py --jsons sample1.confusion.json sample2.confusion.json ... \
        -o validation_report.html
    (writes validation_report.html and validation_report.pdf by default;
     use --no-pdf / --no-html to skip one, or --pdf-output to name it explicitly)
"""

import argparse
import json
import tempfile
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    Image,
    PageBreak,
)


REQUIRED_FIELDS = [
    "sample_id",
    "read_type",
    "tp",
    "fp",
    "fn",
    "tn",
    "precision",
    "recall",
    "specificity",
    "f1",
    "untraceable_reads",
    "ref_truth_total",
    "background_truth_total",
    "deacon_seqs_in",
    "deacon_seqs_out",
    "count_consistency_ok",
]

RECALL_REFERENCE_LINE = 0.95
BLUE = "#2f6fed"
LIGHT_BLUE = "#f4f9ff"


# --------------------------------------------------------------------------
# Data loading / shared logic
# --------------------------------------------------------------------------


def load_records(json_paths):
    records = []
    for path in json_paths:
        with open(path) as fh:
            rec = json.load(fh)
        missing = [f for f in REQUIRED_FIELDS if f not in rec]
        if missing:
            raise ValueError(f"{path} is missing expected field(s): {missing}")
        records.append(rec)
    if not records:
        raise ValueError("No confusion JSONs provided.")
    return pd.DataFrame(records)


def fmt_pct(x):
    return "n/a" if pd.isna(x) else f"{x * 100:.2f}%"


def caveats_table(df):
    rows = []
    for _, row in df.iterrows():
        issues = []
        if not row.count_consistency_ok:
            issues.append(
                f"deacon seqs_in ({row.deacon_seqs_in}) does not match "
                f"ref_truth + background_truth ({row.ref_truth_total + row.background_truth_total}) "
                "— the input FASTQ passed to deacon may not match the truth "
                "sets used for this comparison."
            )
        if row.untraceable_reads > 0:
            issues.append(
                f"{row.untraceable_reads} read ID(s) in the isolated-match output "
                "could not be traced to either truth set — check for read-ID "
                "normalisation mismatches (e.g. mate-pair /1 /2 suffixes)."
            )
        if row.recall < RECALL_REFERENCE_LINE:
            issues.append(
                f"Recall is {fmt_pct(row.recall)} — below the "
                f"{RECALL_REFERENCE_LINE * 100:.0f}% reference line. "
                f"{row.fn} true reference read(s) escaped depletion."
            )
        if row.fp > 0:
            issues.append(
                f"{row.fp} background read(s) were incorrectly classified as "
                "reference matches (false positives)."
            )
        rows.append(
            {
                "Sample": row.sample_id,
                "Read type": row.read_type,
                "Caveats": "; ".join(issues) if issues else "None detected",
            }
        )
    return pd.DataFrame(rows)


def summary_table_df(df):
    t = df[
        [
            "sample_id",
            "read_type",
            "tp",
            "fp",
            "fn",
            "tn",
            "recall",
            "specificity",
            "precision",
            "f1",
        ]
    ].copy()
    for col in ["recall", "specificity", "precision", "f1"]:
        t[col] = t[col].map(fmt_pct)
    t.columns = [
        "Sample",
        "Read type",
        "TP",
        "FP",
        "FN",
        "TN",
        "Recall",
        "Specificity",
        "Precision",
        "F1",
    ]
    return t


# --------------------------------------------------------------------------
# HTML report (Plotly, interactive)
# --------------------------------------------------------------------------


def metrics_bar_figure_plotly(df):
    metrics = ["recall", "specificity", "precision", "f1"]
    labels = {
        "recall": "Recall (sensitivity)",
        "specificity": "Specificity",
        "precision": "Precision",
        "f1": "F1",
    }
    fig = go.Figure()
    for m in metrics:
        fig.add_trace(
            go.Bar(
                name=labels[m],
                x=df["sample_id"],
                y=df[m] * 100,
                text=[fmt_pct(v) for v in df[m]],
                textposition="outside",
            )
        )
    fig.update_layout(
        barmode="group",
        title="Recall, specificity, precision and F1 by sample",
        yaxis_title="%",
        yaxis_range=[0, 105],
        legend_title_text="",
        height=450,
        margin=dict(t=60, b=80),
    )
    fig.add_hline(
        y=RECALL_REFERENCE_LINE * 100,
        line_dash="dot",
        line_color="crimson",
        annotation_text=f"{RECALL_REFERENCE_LINE * 100:.0f}% recall reference line",
        annotation_position="top left",
    )
    return fig


def confusion_matrix_grid_plotly(df):
    n = len(df)
    ncols = 2
    nrows = -(-n // ncols)
    fig = make_subplots(
        rows=nrows,
        cols=ncols,
        subplot_titles=[f"{r.sample_id} ({r.read_type})" for r in df.itertuples()],
        horizontal_spacing=0.15,
        vertical_spacing=0.12,
    )
    for i, row in enumerate(df.itertuples()):
        r, c = divmod(i, ncols)
        z = [[row.tp, row.fn], [row.fp, row.tn]]
        text = [[f"TP: {row.tp}", f"FN: {row.fn}"], [f"FP: {row.fp}", f"TN: {row.tn}"]]
        fig.add_trace(
            go.Heatmap(
                z=z,
                text=text,
                texttemplate="%{text}",
                textfont={"size": 13},
                x=["Matched", "Not matched"],
                y=["Ref truth", "Background truth"],
                colorscale=[[0, LIGHT_BLUE], [1, BLUE]],
                showscale=False,
                xgap=3,
                ygap=3,
            ),
            row=r + 1,
            col=c + 1,
        )
    fig.update_yaxes(autorange="reversed")
    fig.update_layout(
        height=320 * nrows, title="Per-sample confusion matrices", margin=dict(t=80)
    )
    return fig


def df_to_html_table(d, row_colors=None):
    header_cells = "".join(f"<th>{c}</th>" for c in d.columns)
    body_rows = []
    for i, (_, r) in enumerate(d.iterrows()):
        color = f' style="background:{row_colors[i]}"' if row_colors else ""
        cells = "".join(f"<td>{v}</td>" for v in r)
        body_rows.append(f"<tr{color}>{cells}</tr>")
    return f"<table><thead><tr>{header_cells}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def build_html_report(df, caveats_df, outpath):
    n = len(df)
    any_caveats = (caveats_df["Caveats"] != "None detected").any()
    summary_fig = metrics_bar_figure_plotly(df)
    matrix_fig = confusion_matrix_grid_plotly(df)
    summary_table = summary_table_df(df)
    caveat_row_colors = [
        "#fff3f3" if c != "None detected" else "#f3fff5" for c in caveats_df["Caveats"]
    ]

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Reference removal validation report</title>
<style>
  body {{ font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
          margin: 0; padding: 0 0 60px 0; background: #fafbfc; color: #1a1a1a; }}
  header {{ background: #1f2d3d; color: white; padding: 28px 40px; }}
  header h1 {{ margin: 0 0 6px 0; font-size: 22px; }}
  header p {{ margin: 0; color: #b9c4d1; font-size: 14px; }}
  main {{ max-width: 1100px; margin: 0 auto; padding: 30px 40px; }}
  section {{ margin-bottom: 40px; }}
  h2 {{ font-size: 17px; border-bottom: 1px solid #e1e5ea; padding-bottom: 8px; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 13px; margin-top: 10px; }}
  th, td {{ text-align: left; padding: 8px 10px; border-bottom: 1px solid #e8ebee; }}
  th {{ background: #f0f2f5; font-weight: 600; }}
  .banner {{ padding: 14px 18px; border-radius: 6px; font-size: 14px; margin-bottom: 20px; }}
  .banner.warn {{ background: #fff3f3; border: 1px solid #f4b8b8; color: #7a1f1f; }}
  .banner.ok {{ background: #f3fff5; border: 1px solid #b8e6c2; color: #1f6b34; }}
  .legend {{ font-size: 13px; color: #555; margin-top: 8px; }}
</style>
</head>
<body>
<header>
  <h1>Reference removal validation report</h1>
  <p>{n} sample(s) &middot; generated from spiked-sample deacon output vs. known reference/background read IDs</p>
</header>
<main>

<section>
  <div class="banner {"warn" if any_caveats else "ok"}">
    {"⚠ One or more samples have caveats flagged below — review before treating results as clean." if any_caveats else "✓ No data-integrity caveats detected across samples."}
  </div>
  <h2>Summary</h2>
  {df_to_html_table(summary_table)}
  <p class="legend">
    <b>Recall (sensitivity)</b> — proportion of true reference reads that were removed. This is the
    safety-critical metric: low recall means reference material is escaping depletion. &nbsp;
    <b>Precision</b> — proportion of reads flagged as reference matches that were actually reference reads.
    &nbsp; <b>Specificity</b> — proportion of true background reads correctly left alone.
  </p>
</section>

<section>
  <h2>Recall / specificity / precision / F1 by sample</h2>
  {summary_fig.to_html(full_html=False, include_plotlyjs="cdn")}
</section>

<section>
  <h2>Per-sample confusion matrices</h2>
  {matrix_fig.to_html(full_html=False, include_plotlyjs=False)}
</section>

<section>
  <h2>Caveats</h2>
  {df_to_html_table(caveats_df, row_colors=caveat_row_colors)}
</section>

</main>
</body>
</html>
"""
    with open(outpath, "w") as fh:
        fh.write(html)


# --------------------------------------------------------------------------
# PDF report (matplotlib for static charts, reportlab for layout)
# --------------------------------------------------------------------------


def mpl_metrics_bar_image(df, path):
    metrics = ["recall", "specificity", "precision", "f1"]
    labels = ["Recall", "Specificity", "Precision", "F1"]
    x = range(len(df))
    width = 0.2

    fig, ax = plt.subplots(figsize=(9, 4.2), dpi=200)
    for i, (m, lab) in enumerate(zip(metrics, labels)):
        offsets = [xi + (i - 1.5) * width for xi in x]
        vals = (df[m] * 100).tolist()
        ax.bar(offsets, vals, width=width, label=lab)

    ax.axhline(
        RECALL_REFERENCE_LINE * 100, color="crimson", linestyle=":", linewidth=1.2
    )
    ax.text(
        len(df) - 0.5,
        RECALL_REFERENCE_LINE * 100 + 1.5,
        f"{RECALL_REFERENCE_LINE * 100:.0f}% recall reference",
        color="crimson",
        fontsize=8,
        ha="right",
    )

    ax.set_xticks(list(x))
    ax.set_xticklabels(
        [_wrap_label(s, width=20) for s in df["sample_id"]],
        rotation=20,
        ha="right",
        fontsize=8,
    )
    ax.set_ylabel("%")
    ax.set_ylim(0, 108)
    ax.set_title("Recall, specificity, precision and F1 by sample")
    ax.legend(loc="lower right", fontsize=8, ncol=4)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _wrap_label(text, width=28):
    """Break a long sample_id onto multiple lines so matplotlib titles
    don't clip/overflow a fixed-size figure."""
    import textwrap

    return "\n".join(textwrap.wrap(text, width=width)) or text


def mpl_confusion_matrix_image(row, path):
    z = [[row.tp, row.fn], [row.fp, row.tn]]
    fig, ax = plt.subplots(figsize=(3.6, 3.4), dpi=200)
    ax.imshow(z, cmap="Blues")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Matched", "Not matched"], fontsize=8)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Ref truth", "Background truth"], fontsize=8)
    labels = [["TP", "FN"], ["FP", "TN"]]
    vmax = max(max(r) for r in z) or 1
    for i in range(2):
        for j in range(2):
            color = "white" if z[i][j] > vmax * 0.6 else "black"
            ax.text(
                j,
                i,
                f"{labels[i][j]}\n{z[i][j]}",
                ha="center",
                va="center",
                fontsize=9,
                color=color,
            )
    ax.set_title(f"{_wrap_label(row.sample_id)}\n({row.read_type})", fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def reportlab_table(
    data, header_bg=colors.HexColor("#f0f2f5"), row_bg_flags=None, col_widths=None
):
    t = Table(data, colWidths=col_widths, repeatRows=1)
    style = [
        ("BACKGROUND", (0, 0), (-1, 0), header_bg),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#e1e5ea")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]
    if row_bg_flags:
        for i, flagged in enumerate(row_bg_flags):
            if flagged:
                style.append(
                    ("BACKGROUND", (0, i + 1), (-1, i + 1), colors.HexColor("#fff3f3"))
                )
    t.setStyle(TableStyle(style))
    return t


def build_pdf_report(df, caveats_df, outpath):
    any_caveats = (caveats_df["Caveats"] != "None detected").any()
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "TitleBig", parent=styles["Title"], fontSize=18, spaceAfter=4
    )
    subtitle_style = ParagraphStyle(
        "Subtitle", parent=styles["Normal"], textColor=colors.HexColor("#555")
    )
    h2_style = ParagraphStyle(
        "H2", parent=styles["Heading2"], spaceBefore=14, spaceAfter=6
    )
    body_style = ParagraphStyle("Body", parent=styles["Normal"], fontSize=9, leading=13)
    banner_style = ParagraphStyle(
        "Banner",
        parent=styles["Normal"],
        fontSize=10,
        backColor=colors.HexColor("#fff3f3")
        if any_caveats
        else colors.HexColor("#f3fff5"),
        borderColor=colors.HexColor("#f4b8b8")
        if any_caveats
        else colors.HexColor("#b8e6c2"),
        borderWidth=1,
        borderPadding=8,
        textColor=colors.HexColor("#7a1f1f")
        if any_caveats
        else colors.HexColor("#1f6b34"),
    )

    story = []
    story.append(Paragraph("Reference removal validation report", title_style))
    story.append(
        Paragraph(
            f"{len(df)} sample(s) &middot; generated from spiked-sample deacon output vs. known "
            "reference/background read IDs",
            subtitle_style,
        )
    )
    story.append(Spacer(1, 12))

    banner_text = (
        "WARNING: One or more samples have caveats flagged below — review before treating results as clean."
        if any_caveats
        else "OK: No data-integrity caveats detected across samples."
    )
    story.append(Paragraph(banner_text, banner_style))
    story.append(Spacer(1, 10))

    # Summary table
    cell_style = ParagraphStyle("Cell", parent=styles["Normal"], fontSize=8, leading=10)

    story.append(Paragraph("Summary", h2_style))
    summary_table = summary_table_df(df)
    table_data = [list(summary_table.columns)] + [
        [Paragraph(str(v), cell_style) for v in row]
        for row in summary_table.values.tolist()
    ]
    story.append(
        reportlab_table(
            table_data,
            col_widths=[
                4 * cm,
                2 * cm,
                1.3 * cm,
                1.3 * cm,
                1.3 * cm,
                1.5 * cm,
                1.8 * cm,
                2 * cm,
                1.8 * cm,
                1.3 * cm,
            ],
        )
    )
    story.append(Spacer(1, 8))
    story.append(
        Paragraph(
            "<b>Recall (sensitivity)</b> — proportion of true reference reads that were removed. "
            "This is the safety-critical metric: low recall means reference material is escaping "
            "depletion. <b>Precision</b> — proportion of reads flagged as reference matches that were "
            "actually reference reads. <b>Specificity</b> — proportion of true background reads "
            "correctly left alone.",
            body_style,
        )
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Metrics bar chart
        bar_path = tmpdir / "bar.png"
        mpl_metrics_bar_image(df, bar_path)
        story.append(
            Paragraph("Recall / specificity / precision / F1 by sample", h2_style)
        )
        story.append(Image(str(bar_path), width=16 * cm, height=16 * cm * (4.2 / 9)))

        # Confusion matrices, two per row
        story.append(PageBreak())
        story.append(Paragraph("Per-sample confusion matrices", h2_style))
        cm_paths = []
        for i, row in enumerate(df.itertuples()):
            p = tmpdir / f"cm_{i}.png"
            mpl_confusion_matrix_image(row, p)
            cm_paths.append(p)

        img_flowables = [
            Image(str(p), width=7.5 * cm, height=7.5 * cm * (3.0 / 3.6))
            for p in cm_paths
        ]
        rows = [img_flowables[i : i + 2] for i in range(0, len(img_flowables), 2)]
        if rows and len(rows[-1]) == 1:
            rows[-1].append("")
        story.append(Table(rows, colWidths=[8 * cm, 8 * cm]))

        # Caveats
        story.append(PageBreak())
        story.append(Paragraph("Caveats", h2_style))
        caveat_data = [list(caveats_df.columns)]
        for _, r in caveats_df.iterrows():
            caveat_data.append(
                [
                    Paragraph(r["Sample"], cell_style),
                    Paragraph(r["Read type"], cell_style),
                    Paragraph(r["Caveats"], cell_style),
                ]
            )
        flagged = [c != "None detected" for c in caveats_df["Caveats"]]
        story.append(
            reportlab_table(
                caveat_data,
                row_bg_flags=flagged,
                col_widths=[3.5 * cm, 2.5 * cm, 10 * cm],
            )
        )

        doc = SimpleDocTemplate(
            str(outpath),
            pagesize=A4,
            leftMargin=1.8 * cm,
            rightMargin=1.8 * cm,
            topMargin=1.8 * cm,
            bottomMargin=1.8 * cm,
        )
        doc.build(story)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def main():
    global RECALL_REFERENCE_LINE
    p = argparse.ArgumentParser()
    p.add_argument(
        "--jsons", nargs="+", required=True, help="per-sample confusion JSON files"
    )
    p.add_argument("-o", "--output", required=True, help="HTML output path")
    p.add_argument(
        "--pdf-output",
        default=None,
        help="PDF output path (default: same basename as --output, .pdf extension)",
    )
    p.add_argument("--no-html", action="store_true", help="skip HTML generation")
    p.add_argument("--no-pdf", action="store_true", help="skip PDF generation")
    p.add_argument(
        "--recall-threshold",
        type=float,
        default=RECALL_REFERENCE_LINE,
        help="recall reference line / caveat threshold (default: 0.95)",
    )
    args = p.parse_args()

    RECALL_REFERENCE_LINE = args.recall_threshold

    df = load_records(args.jsons)
    df = df.sort_values(["read_type", "sample_id"]).reset_index(drop=True)
    caveats_df = caveats_table(df)

    if not args.no_html:
        build_html_report(df, caveats_df, args.output)

    if not args.no_pdf:
        pdf_path = args.pdf_output or str(Path(args.output).with_suffix(".pdf"))
        build_pdf_report(df, caveats_df, pdf_path)


if __name__ == "__main__":
    main()
