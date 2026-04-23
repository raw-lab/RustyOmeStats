#!/usr/bin/env python3
"""Render seaborn plots from the CSV files written by `rustyomestats`, and
bundle them into a single self-contained HTML report.

Inputs: CSV files produced by `rustyomestats genome` and/or `rustyomestats u50`.
Outputs: one PNG per figure, plus `report.html` with every PNG embedded
(base64) so the HTML can be shared as a single file.

Stack: polars (read CSV), seaborn + matplotlib (plots).
"""
from __future__ import annotations

import argparse
import base64
import html
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns

FRAME_ORDER = ["+1", "+2", "+3", "-1", "-2", "-3"]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _save(fig, outdir: Path, name: str) -> Path:
    path = outdir / f"{name}.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# per-sequence: length distribution, GC distribution, GC vs length
# ---------------------------------------------------------------------------
def plot_length_hist(outdir: Path) -> Path | None:
    p = outdir / "per_sequence.csv"
    if not p.exists():
        return None
    df = pl.read_csv(p)
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.histplot(
        df["length"].to_numpy(),
        bins=50,
        log_scale=(True, False),
        color="steelblue",
        ax=ax,
    )
    ax.set(title="Sequence length distribution",
           xlabel="length (bp, log10)", ylabel="# sequences")
    sns.despine()
    fig.tight_layout()
    return _save(fig, outdir, "plot_length_histogram")


def plot_gc_distribution(outdir: Path) -> Path | None:
    p = outdir / "per_sequence.csv"
    if not p.exists():
        return None
    df = pl.read_csv(p)
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.histplot(df["gc"].to_numpy(), bins=40, kde=True, color="seagreen", ax=ax)
    ax.set(title="GC% per sequence", xlabel="GC %", ylabel="count")
    sns.despine()
    fig.tight_layout()
    return _save(fig, outdir, "plot_gc_distribution")


def plot_gc_vs_length(outdir: Path) -> Path | None:
    p = outdir / "per_sequence.csv"
    if not p.exists():
        return None
    df = pl.read_csv(p)
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(
        x=df["length"].to_numpy(),
        y=df["gc"].to_numpy(),
        s=14, alpha=0.5, ax=ax,
    )
    ax.set_xscale("log")
    ax.set(title="GC% vs sequence length", xlabel="length (bp)", ylabel="GC %")
    sns.despine()
    fig.tight_layout()
    return _save(fig, outdir, "plot_gc_vs_length")


# ---------------------------------------------------------------------------
# codon density
# ---------------------------------------------------------------------------
def plot_codon_usage_bar(outdir: Path) -> Path | None:
    abs_p  = outdir / "codon_absolute_aggregate.csv"
    pred_p = outdir / "codon_predicted_aggregate.csv"
    parts: list[pl.DataFrame] = []
    if abs_p.exists():
        parts.append(pl.read_csv(abs_p).with_columns(
            pl.lit("6-frame (absolute)").alias("source")))
    if pred_p.exists():
        parts.append(pl.read_csv(pred_p).with_columns(
            pl.lit("FragGeneScan (predicted)").alias("source")))
    if not parts:
        return None
    pdf = pl.concat(parts, how="vertical").to_pandas()

    fig, ax = plt.subplots(figsize=(14, 5))
    sns.barplot(data=pdf, x="codon", y="density", hue="source", ax=ax)
    ax.set(title="Codon density: 6-frame vs FragGeneScan",
           ylabel="density (fraction)")
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    sns.despine()
    fig.tight_layout()
    return _save(fig, outdir, "plot_codon_usage_bar")


def plot_codon_heatmap_by_frame(outdir: Path) -> Path | None:
    p = outdir / "codon_absolute.csv"
    if not p.exists():
        return None
    agg = (
        pl.read_csv(p)
          .group_by(["frame", "codon"])
          .agg(pl.col("count").sum())
          .with_columns(
              (pl.col("count") / pl.col("count").sum().over("frame"))
              .alias("density")
          )
    )
    pivot = (
        agg.to_pandas()
           .pivot(index="frame", columns="codon", values="density")
           .reindex(FRAME_ORDER)
    )
    fig, ax = plt.subplots(figsize=(18, 4))
    sns.heatmap(pivot, cmap="viridis", ax=ax, cbar_kws={"label": "density"})
    ax.set(title="6-frame codon density heatmap",
           xlabel="codon", ylabel="frame")
    ax.tick_params(axis="x", labelsize=7)
    fig.tight_layout()
    return _save(fig, outdir, "plot_codon_heatmap_by_frame")


def plot_codon_enrichment(outdir: Path) -> Path | None:
    p = outdir / "codon_comparison.csv"
    if not p.exists():
        return None
    pdf = (
        pl.read_csv(p)
          .sort("enrichment_pred_over_abs", descending=True)
          .to_pandas()
    )
    fig, ax = plt.subplots(figsize=(14, 5))
    sns.barplot(data=pdf, x="codon", y="enrichment_pred_over_abs",
                color="firebrick", ax=ax)
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1)
    ax.set(title="Codon enrichment (predicted / absolute)",
           ylabel="enrichment")
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    sns.despine()
    fig.tight_layout()
    return _save(fig, outdir, "plot_codon_enrichment")


# ---------------------------------------------------------------------------
# U50 plots
# ---------------------------------------------------------------------------
def _read_u50_summary(outdir: Path) -> dict | None:
    p = outdir / "u50_summary.csv"
    if not p.exists():
        return None
    df = pl.read_csv(p, schema_overrides={"metric": pl.Utf8, "value": pl.Utf8})
    return {row["metric"]: row["value"] for row in df.iter_rows(named=True)}


def plot_u50_contig_lengths(outdir: Path) -> Path | None:
    p = outdir / "u50_contigs.csv"
    if not p.exists():
        return None
    df = pl.read_csv(p)
    if df.height == 0:
        return None
    long = pl.concat(
        [
            df.select(pl.col("orig_length").alias("length"),
                      pl.lit("original").alias("kind")),
            df.select(pl.col("unique_bp").alias("length"),
                      pl.lit("unique").alias("kind")),
        ],
        how="vertical",
    ).filter(pl.col("length") > 0)

    pdf = long.to_pandas()
    fig, ax = plt.subplots(figsize=(9, 5))
    sns.histplot(data=pdf, x="length", hue="kind", bins=40,
                 log_scale=(True, False), multiple="layer",
                 palette={"original": "steelblue", "unique": "firebrick"},
                 alpha=0.55, ax=ax)
    ax.set(title="Contig lengths: original vs unique-bp (post-masking)",
           xlabel="length (bp, log10)", ylabel="# contigs")
    sns.despine()
    fig.tight_layout()
    return _save(fig, outdir, "plot_u50_contig_lengths")


def plot_u50_summary_bars(outdir: Path) -> Path | None:
    summ = _read_u50_summary(outdir)
    if summ is None:
        return None
    metric_order = ["N50", "NG50", "U50", "UG50"]
    try:
        lengths = [int(summ[m]) for m in metric_order]
    except (KeyError, ValueError):
        return None
    ref_len  = int(summ.get("ref_length", "0"))
    ug50_pct = float(summ.get("UG50_percent", "0") or 0)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    pal = ["#4c72b0", "#55a868", "#c44e52", "#8172b2"]
    sns.barplot(x=metric_order, y=lengths, hue=metric_order,
                palette=pal, ax=ax, legend=False)
    for i, v in enumerate(lengths):
        ax.text(i, v, f"{v}", ha="center", va="bottom", fontsize=9)
    if ref_len > 0:
        ax.axhline(ref_len, color="grey", linestyle=":", linewidth=1,
                   label=f"ref length = {ref_len}")
        ax.legend(loc="upper right")
    ax.set(title=f"Castro assembly metrics (UG50% = {ug50_pct:.2f}%)",
           ylabel="length (bp)", xlabel="")
    sns.despine()
    fig.tight_layout()
    return _save(fig, outdir, "plot_u50_summary")


def plot_u50_coverage(outdir: Path, bins: int = 500) -> Path | None:
    summ = _read_u50_summary(outdir)
    gi = outdir / "u50_gap_intervals.csv"
    if summ is None or not gi.exists():
        return None
    try:
        ref_len = int(summ["ref_length"])
    except (KeyError, ValueError):
        return None
    if ref_len == 0:
        return None

    mask = np.ones(ref_len, dtype=np.uint8)
    gaps = pl.read_csv(gi)
    if gaps.height > 0:
        for s, e in zip(gaps["start"].to_numpy(), gaps["end"].to_numpy()):
            mask[int(s):int(e)] = 0

    bin_size = max(1, ref_len // bins)
    n_bins = ref_len // bin_size
    if n_bins == 0:
        return None
    trimmed = mask[: n_bins * bin_size].reshape(n_bins, bin_size)
    cov = trimmed.mean(axis=1)
    xs = np.arange(n_bins) * bin_size

    fig, ax = plt.subplots(figsize=(12, 3))
    ax.fill_between(xs, 0, cov, step="mid", color="steelblue", alpha=0.8)
    ax.set(title="Reference coverage by mapped contigs (gaps carved)",
           xlabel=f"reference position (bins of {bin_size} bp)",
           ylabel="fraction covered", ylim=(0, 1.05))
    sns.despine()
    fig.tight_layout()
    return _save(fig, outdir, "plot_u50_coverage")


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------
HTML_HEAD = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>rustyomestats report</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         max-width: 1100px; margin: 2em auto; padding: 0 1em; color: #222; }}
  h1 {{ border-bottom: 2px solid #444; padding-bottom: .3em; }}
  h2 {{ margin-top: 2em; color: #2a4b7c; }}
  .meta {{ color: #666; font-size: .9em; }}
  figure {{ margin: 1em 0 2em 0; text-align: center; }}
  figure img {{ max-width: 100%; height: auto;
                border: 1px solid #ddd; border-radius: 4px; }}
  figcaption {{ font-size: .9em; color: #555; margin-top: .4em; }}
  table {{ border-collapse: collapse; margin: 1em 0; }}
  th, td {{ border: 1px solid #ccc; padding: .4em .8em; text-align: left; }}
  th {{ background: #f5f5f5; }}
  code {{ background: #f5f5f5; padding: .1em .3em; border-radius: 3px; }}
</style>
</head>
<body>
<h1>rustyomestats report</h1>
<p class="meta">Generated {when} from <code>{outdir}</code>.</p>
"""

HTML_TAIL = """</body>
</html>
"""


def _embed(png_path: Path) -> str:
    data = base64.b64encode(png_path.read_bytes()).decode("ascii")
    return f'<img src="data:image/png;base64,{data}" alt="{png_path.stem}">'


def _render_summary_table(csv_path: Path, title: str) -> str:
    if not csv_path.exists():
        return ""
    df = pl.read_csv(csv_path,
                     schema_overrides={"metric": pl.Utf8, "value": pl.Utf8})
    rows = "".join(
        f"<tr><td><code>{html.escape(r['metric'])}</code></td>"
        f"<td>{html.escape(r['value'])}</td></tr>"
        for r in df.iter_rows(named=True)
    )
    return (
        f"<h2>{html.escape(title)}</h2>\n"
        f"<table><thead><tr><th>metric</th><th>value</th></tr></thead>"
        f"<tbody>{rows}</tbody></table>\n"
    )


def write_html_report(outdir: Path,
                      figures: list[tuple[str, Path | None, str]]) -> Path:
    """
    figures is a list of (section_title, png_path_or_None, caption).
    None entries are skipped silently (that section wasn't produced by Rust).
    """
    pieces = [
        HTML_HEAD.format(
            when=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            outdir=html.escape(str(outdir.resolve())),
        )
    ]

    pieces.append(_render_summary_table(
        outdir / "summary_stats.csv",
        "Genome summary stats (summary_stats.csv)",
    ))
    pieces.append(_render_summary_table(
        outdir / "u50_summary.csv",
        "U50 assembly metrics (u50_summary.csv)",
    ))

    for title, png, caption in figures:
        if png is None or not png.exists():
            continue
        pieces.append(
            f"<h2>{html.escape(title)}</h2>\n"
            f"<figure>{_embed(png)}"
            f"<figcaption>{html.escape(caption)}</figcaption></figure>\n"
        )

    pieces.append(HTML_TAIL)
    out_html = outdir / "report.html"
    out_html.write_text("".join(pieces), encoding="utf-8")
    return out_html


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Plot rustyomestats CSV outputs with seaborn; "
                    "emit PNGs and a self-contained HTML report.",
    )
    ap.add_argument(
        "-d", "--dir",
        type=Path,
        default=Path("rustyomestats_out"),
        help="rustyomestats output directory (default: rustyomestats_out)",
    )
    ap.add_argument(
        "--no-html",
        action="store_true",
        help="skip report.html generation (PNGs still written)",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    if not args.dir.is_dir():
        sys.exit(f"error: {args.dir} is not a directory")

    sns.set_theme(style="whitegrid", context="notebook")

    figures = [
        ("Length histogram",
         plot_length_hist(args.dir),
         "Distribution of sequence lengths (log10 x-axis)."),
        ("GC distribution",
         plot_gc_distribution(args.dir),
         "Per-sequence GC% (ambiguous bases excluded from the denominator)."),
        ("GC vs length",
         plot_gc_vs_length(args.dir),
         "Per-sequence GC% plotted against log-length."),
        ("Codon usage (6-frame vs FragGeneScan)",
         plot_codon_usage_bar(args.dir),
         "64-codon density. Compare the 6-frame background (absolute) to "
         "predicted ORFs (FGS)."),
        ("6-frame codon heatmap",
         plot_codon_heatmap_by_frame(args.dir),
         "Codon density across all 6 reading frames."),
        ("Codon enrichment (predicted / absolute)",
         plot_codon_enrichment(args.dir),
         "Codons overrepresented in predicted ORFs relative to the "
         "6-frame background (dashed line = no enrichment)."),
        ("U50: contig length distributions",
         plot_u50_contig_lengths(args.dir),
         "Original contig lengths vs unique-bp lengths after greedy masking."),
        ("U50: assembly metric summary",
         plot_u50_summary_bars(args.dir),
         "N50 / NG50 / U50 / UG50. Dotted line = reference length."),
        ("U50: reference coverage",
         plot_u50_coverage(args.dir),
         "Binned fraction of the reference covered by at least one contig."),
    ]

    produced = [png for _, png, _ in figures if png is not None]
    print(f"[rustyomestats] wrote {len(produced)} PNG(s) to {args.dir}")
    for p in produced:
        print(f"                {p.name}")

    if not args.no_html:
        html_path = write_html_report(args.dir, figures)
        print(f"[rustyomestats] wrote {html_path.name}")


if __name__ == "__main__":
    main()
