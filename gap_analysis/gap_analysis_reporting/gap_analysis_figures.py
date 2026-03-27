#!/usr/bin/env python3
"""
UKBOL Gap Analysis Figure Generator
====================================
Produces publication-quality figures from gap analysis CSV files.

Usage:
    python gap_analysis_figures.py input.csv
    python gap_analysis_figures.py input.csv --output-dir figures/
    python gap_analysis_figures.py input.csv --prefix dataset1 --format pdf png svg
    python gap_analysis_figures.py input.csv --top-n 25 --dpi 600

Requirements:
    pip install pandas matplotlib numpy

Expects a CSV with at minimum these columns:
    - species_status:  GREEN | AMBER | BLUE | ORANGE | RED | BLACK
    - bags_grade:      A | B | C | D | E | F
    - phylum_division: Taxonomic phylum
    - class:           Taxonomic class
    - order:           Taxonomic order
    - family:          Taxonomic family
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit these if your column names or categories differ
# ═══════════════════════════════════════════════════════════════════════════════

STATUS_ORDER = ["GREEN", "AMBER", "BLUE", "ORANGE", "RED", "BLACK"]
STATUS_COLORS = {
    "GREEN":  "#2E7D32",
    "AMBER":  "#F57F17",
    "BLUE":   "#1565C0",
    "ORANGE": "#E65100",
    "RED":    "#C62828",
    "BLACK":  "#333333",
}
STATUS_LABELS = {
    "GREEN":  "Green – Only valid name has records",
    "AMBER":  "Amber – Valid name & synonym(s) have records",
    "BLUE":   "Blue – Only synonym(s), valid name absent",
    "ORANGE": "Orange – BIN shared only with non-Linnaean names (interim ID conflict)",
    "RED":    "Red – Taxonomic conflict (Linnaean name shares BIN)",
    "BLACK":  "Black – No records found",
}

BAGS_ORDER = ["A", "B", "C", "D", "E", "F"]
BAGS_COLORS = {
    "A": "#2E7D32", "B": "#43A047", "C": "#F57F17",
    "D": "#FF8F00", "E": "#C62828", "F": "#333333",
}
BAGS_LABELS = {
    "A": "A – Single BIN, >10 records",
    "B": "B – Single BIN, 3–10 records",
    "C": "C – Multiple BINs (potential split)",
    "D": "D – Single BIN, 1–2 records",
    "E": "E – Taxonomic conflict",
    "F": "F – No records / no BIN",
}

# Columns in the input CSV
COL_STATUS = "species_status"
COL_BAGS   = "bags_grade"
COL_PHYLUM = "phylum_division"
COL_CLASS  = "class"
COL_ORDER  = "order"
COL_FAMILY = "family"

TAXONOMIC_LEVELS = [
    ("phylum", COL_PHYLUM),
    ("class",  COL_CLASS),
    ("order",  COL_ORDER),
    ("family", COL_FAMILY),
]


# ═══════════════════════════════════════════════════════════════════════════════
# STYLE
# ═══════════════════════════════════════════════════════════════════════════════

def apply_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "Georgia"],
        "font.size": 9,
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "axes.edgecolor": "#cccccc",
        "axes.linewidth": 0.5,
        "xtick.color": "#444444",
        "ytick.color": "#444444",
        "text.color": "#222222",
        "axes.labelcolor": "#222222",
        "grid.color": "#e8e8e8",
        "grid.linewidth": 0.4,
    })


# ═══════════════════════════════════════════════════════════════════════════════
# PLOTTING HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def make_stacked_horizontal(ax, ct, title, proportional=False, show_n=True):
    """Horizontal stacked bar chart on a given axes."""
    ct = ct.reindex(columns=STATUS_ORDER, fill_value=0)
    ct = ct.loc[ct.sum(axis=1).sort_values(ascending=True).index]
    totals = ct.sum(axis=1)

    if proportional:
        ct_plot = ct.div(totals, axis=0) * 100
    else:
        ct_plot = ct.copy()

    y_pos = np.arange(len(ct_plot))
    left = np.zeros(len(ct_plot))

    for status in STATUS_ORDER:
        if status in ct_plot.columns:
            vals = ct_plot[status].values
            ax.barh(y_pos, vals, left=left, height=0.7,
                    color=STATUS_COLORS[status], edgecolor="white",
                    linewidth=0.3, label=status, zorder=3)
            left += vals

    ax.set_yticks(y_pos)
    if show_n:
        labels = [f"{name}  (n={int(totals[name])})" for name in ct_plot.index]
    else:
        labels = list(ct_plot.index)
    ax.set_yticklabels(labels, fontsize=8, fontstyle="italic")
    ax.set_title(title, fontsize=10, fontweight="bold", pad=10, loc="left")

    if proportional:
        ax.set_xlim(0, 100)
        ax.set_xlabel("Proportion of species (%)", fontsize=8)
        ax.xaxis.set_major_formatter(mticker.PercentFormatter())
    else:
        ax.set_xlabel("Number of species", fontsize=8)

    ax.grid(axis="x", alpha=0.3, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def coverage_color(pct):
    """Return a colour based on coverage percentage thresholds."""
    if pct >= 80:
        return "#2E7D32"
    elif pct >= 50:
        return "#F57F17"
    elif pct >= 20:
        return "#E65100"
    else:
        return "#C62828"


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE GENERATORS
# ═══════════════════════════════════════════════════════════════════════════════

def fig_summary(df):
    """Fig 1: Summary donut (status) + BAGS grade bars."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5),
                             gridspec_kw={"width_ratios": [1, 1.3]})
    total = len(df)

    # ── Panel A: donut ──
    ax = axes[0]
    status_counts = df[COL_STATUS].value_counts()
    pie_data = [status_counts.get(s, 0) for s in STATUS_ORDER]
    pie_colors = [STATUS_COLORS[s] for s in STATUS_ORDER]

    wedges, _, _ = ax.pie(
        pie_data, labels=None, colors=pie_colors, autopct="",
        startangle=90, counterclock=False,
        wedgeprops={"edgecolor": "white", "linewidth": 1.5},
        pctdistance=0.75,
    )
    for i, (w, s) in enumerate(zip(wedges, STATUS_ORDER)):
        ang = (w.theta2 - w.theta1) / 2.0 + w.theta1
        x, y = np.cos(np.deg2rad(ang)), np.sin(np.deg2rad(ang))
        count = pie_data[i]
        pct = count / sum(pie_data) * 100
        if pct > 3:
            ax.annotate(
                f"{count}\n({pct:.1f}%)",
                xy=(0.75 * x, 0.75 * y), fontsize=7, ha="center", va="center",
                fontweight="bold", color="#333" if s in ("AMBER", "ORANGE") else "white",
            )

    centre = plt.Circle((0, 0), 0.45, fc="white", ec="none", zorder=5)
    ax.add_artist(centre)
    ax.text(0, 0.06, f"{total:,}", ha="center", va="center",
            fontsize=16, fontweight="bold", color="#222", zorder=6)
    ax.text(0, -0.1, "species", ha="center", va="center",
            fontsize=8, color="#666", zorder=6)
    ax.set_title("(a) Species by taxonomic concordance status",
                 fontsize=10, fontweight="bold", pad=12, loc="left")

    # ── Panel B: BAGS grades ──
    ax = axes[1]
    bags_counts = df[COL_BAGS].value_counts()
    bags_data = [bags_counts.get(g, 0) for g in BAGS_ORDER]
    y_pos = np.arange(len(BAGS_ORDER))

    ax.barh(y_pos, bags_data, height=0.6,
            color=[BAGS_COLORS[g] for g in BAGS_ORDER],
            edgecolor="white", linewidth=0.5, zorder=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([BAGS_LABELS[g] for g in BAGS_ORDER], fontsize=8)
    ax.invert_yaxis()
    for i, (v, g) in enumerate(zip(bags_data, BAGS_ORDER)):
        pct = v / total * 100
        ax.text(v + total * 0.01, i, f"{v:,} ({pct:.1f}%)",
                va="center", fontsize=7.5, color="#444")
    ax.set_xlabel("Number of species", fontsize=8)
    ax.set_title("(b) Species by BAGS grade",
                 fontsize=10, fontweight="bold", pad=12, loc="left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", alpha=0.3, zorder=0)
    ax.set_xlim(0, max(bags_data) * 1.25)

    fig.tight_layout(w_pad=3)

    # Shared legend
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=STATUS_COLORS[s], edgecolor="none")
               for s in STATUS_ORDER]
    fig.legend(handles, [STATUS_LABELS[s] for s in STATUS_ORDER],
               loc="lower center", ncol=2, fontsize=7, frameon=False,
               bbox_to_anchor=(0.35, -0.12),
               handlelength=1.2, handletextpad=0.5, columnspacing=1.5)

    return fig


def fig_stacked_bars(df, col, level_name, top_n=None):
    """Fig: Paired absolute + proportional stacked bars for a taxonomic level."""
    ct = pd.crosstab(df[col], df[COL_STATUS])
    ct = ct.reindex(columns=STATUS_ORDER, fill_value=0)

    if top_n is not None:
        top_taxa = ct.sum(axis=1).nlargest(top_n).index
        ct = ct.loc[top_taxa]

    n_taxa = len(ct)
    height = max(5, n_taxa * 0.32 + 1.5)

    fig, axes = plt.subplots(1, 2, figsize=(12, height),
                             gridspec_kw={"width_ratios": [1.1, 1]})
    make_stacked_horizontal(axes[0], ct,
                            f"(a) Absolute counts by {level_name}")
    make_stacked_horizontal(axes[1], ct,
                            f"(b) Proportional by {level_name}",
                            proportional=True, show_n=False)
    fig.tight_layout(w_pad=4)
    return fig


def fig_coverage_heatmap(df, top_n=30):
    """Fig: Horizontal bars showing % of species with any records, by order."""
    ct = pd.crosstab(df[COL_ORDER], df[COL_STATUS])
    ct = ct.reindex(columns=STATUS_ORDER, fill_value=0)
    totals = ct.sum(axis=1)

    has_records = ct.drop(columns="BLACK", errors="ignore").sum(axis=1)
    pct_with_records = has_records / totals * 100

    # Top N by species richness, then sort ascending for horizontal bar
    top_idx = totals.nlargest(top_n).index
    pct_plot = pct_with_records[top_idx].sort_values(ascending=False)
    # Reverse so largest-count is at top
    sort_by_count = totals[top_idx].sort_values(ascending=True)
    pct_plot = pct_with_records[sort_by_count.index]
    totals_plot = totals[sort_by_count.index]

    height = max(6, len(pct_plot) * 0.28 + 2)
    fig, ax = plt.subplots(figsize=(8, height))

    y_pos = np.arange(len(pct_plot))
    colors = [coverage_color(p) for p in pct_plot.values]

    ax.barh(y_pos, pct_plot.values, height=0.7, color=colors,
            edgecolor="white", linewidth=0.3, zorder=3, alpha=0.85)
    ax.set_yticks(y_pos)
    labels = [f"{name}  (n={int(totals_plot[name])})" for name in pct_plot.index]
    ax.set_yticklabels(labels, fontsize=7.5, fontstyle="italic")
    ax.set_xlim(0, 108)
    ax.set_xlabel("Species with records (%)", fontsize=9)
    ax.set_title(
        f"Proportion of species with BIN records by order\n"
        f"(top {len(pct_plot)} orders ranked by species richness)",
        fontsize=10, fontweight="bold", pad=12, loc="left",
    )

    for i, (v, name) in enumerate(zip(pct_plot.values, pct_plot.index)):
        ax.text(v + 1.2, i, f"{v:.0f}%", va="center", fontsize=7, color="#444")

    ax.axvline(x=50, color="#aaa", linestyle="--", linewidth=0.6, alpha=0.6, zorder=1)
    ax.axvline(x=80, color="#aaa", linestyle="--", linewidth=0.6, alpha=0.6, zorder=1)
    ax.text(50, len(pct_plot) - 0.3, "50%", fontsize=7, color="#999", ha="center")
    ax.text(80, len(pct_plot) - 0.3, "80%", fontsize=7, color="#999", ha="center")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", alpha=0.2, zorder=0)

    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate publication-quality gap analysis figures.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("input_csv", type=Path,
                        help="Path to the gap analysis CSV file.")
    parser.add_argument("--output-dir", "-o", type=Path, default=None,
                        help="Output directory (default: same directory as input).")
    parser.add_argument("--prefix", "-p", type=str, default=None,
                        help="Filename prefix (default: derived from input filename).")
    parser.add_argument("--format", "-f", nargs="+", default=["pdf", "png"],
                        choices=["pdf", "png", "svg", "eps", "tiff"],
                        help="Output formats (default: pdf png).")
    parser.add_argument("--dpi", type=int, default=300,
                        help="Resolution for raster formats (default: 300).")
    parser.add_argument("--top-n", type=int, default=None,
                        help="Max taxa to show per level. Default: all for phylum, "
                             "15 for class, 20 for order, 30 for family.")
    parser.add_argument("--skip", nargs="*", default=[],
                        choices=["summary", "phylum", "class", "order",
                                 "family", "coverage"],
                        help="Figure(s) to skip.")
    return parser.parse_args()


def save_figure(fig, output_dir, name, formats, dpi):
    """Save a figure in all requested formats."""
    paths = []
    for fmt in formats:
        path = output_dir / f"{name}.{fmt}"
        fig.savefig(path, bbox_inches="tight", dpi=dpi, facecolor="white")
        paths.append(path)
    plt.close(fig)
    return paths


def main():
    args = parse_args()

    # ── Validate input ──
    if not args.input_csv.exists():
        print(f"Error: {args.input_csv} not found.", file=sys.stderr)
        sys.exit(1)

    # ── Resolve output directory and prefix ──
    output_dir = args.output_dir or args.input_csv.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or args.input_csv.stem

    # ── Load data ──
    print(f"Loading {args.input_csv} ...")
    df = pd.read_csv(args.input_csv, low_memory=False)

    # Validate required columns
    required = [COL_STATUS, COL_BAGS, COL_PHYLUM, COL_CLASS, COL_ORDER, COL_FAMILY]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"Error: Missing columns: {missing}", file=sys.stderr)
        print(f"Available columns: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    print(f"  {len(df):,} species across {df[COL_PHYLUM].nunique()} phyla, "
          f"{df[COL_CLASS].nunique()} classes, {df[COL_ORDER].nunique()} orders, "
          f"{df[COL_FAMILY].nunique()} families")
    print(f"  Status: {df[COL_STATUS].value_counts().to_dict()}")
    print(f"  Output: {output_dir}/  prefix={prefix}  formats={args.format}")
    print()

    apply_style()

    # ── Default top-N per level ──
    default_top_n = {
        "phylum": None,   # show all
        "class":  15,
        "order":  20,
        "family": 30,
    }

    figures_saved = []

    # ── Fig 1: Summary ──
    if "summary" not in args.skip:
        print("Generating Fig 1: Summary ...")
        fig = fig_summary(df)
        paths = save_figure(fig, output_dir, f"{prefix}_fig1_summary",
                            args.format, args.dpi)
        figures_saved.extend(paths)

    # ── Figs 2–5: Stacked bars by taxonomic level ──
    for i, (level_name, col) in enumerate(TAXONOMIC_LEVELS, start=2):
        if level_name in args.skip:
            continue
        top_n = args.top_n if args.top_n is not None else default_top_n[level_name]
        label = level_name
        if top_n:
            label += f" (top {top_n})"
        print(f"Generating Fig {i}: {label} ...")
        fig = fig_stacked_bars(df, col, level_name, top_n=top_n)
        paths = save_figure(fig, output_dir, f"{prefix}_fig{i}_{level_name}",
                            args.format, args.dpi)
        figures_saved.extend(paths)

    # ── Fig 6: Coverage heatmap ──
    if "coverage" not in args.skip:
        cov_n = args.top_n or 30
        print(f"Generating Fig 6: Coverage heatmap (top {cov_n} orders) ...")
        fig = fig_coverage_heatmap(df, top_n=cov_n)
        paths = save_figure(fig, output_dir, f"{prefix}_fig6_coverage",
                            args.format, args.dpi)
        figures_saved.extend(paths)

    # ── Done ──
    print(f"\nDone — {len(figures_saved)} files saved to {output_dir}/")
    for p in sorted(set(str(f.parent) for f in figures_saved)):
        files_in_dir = [f.name for f in figures_saved if str(f.parent) == p]
        for fn in sorted(files_in_dir):
            print(f"  {fn}")


if __name__ == "__main__":
    main()
