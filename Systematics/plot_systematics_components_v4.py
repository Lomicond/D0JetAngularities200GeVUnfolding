#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_systematics_components_v3.py

Plot individual systematic-component contributions from the wide systematics table.

Input:
  Systematics/systematics_components_wide_v3.tsv

Output by default:
  Systematics/Plots/SystematicsComponents/all_systematics_components_abs.pdf
  Systematics/Plots/SystematicsComponents/png/*.png

The script intentionally skips all *_max_abs_pct summary columns.  It plots the
individual variations/components only, e.g. sWeight Cheby2, D0 PID up/down,
JetsReco DCA, prior-shape +/- and binning RMS.  In addition to the
original all-components plots, it also creates one set of plots per systematic
source category, following the summary slide structure.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


BASE_COLUMNS = {
    "job_label", "observable", "observable_pretty", "hist_name",
    "centrality", "centrality_label", "method", "iteration_display",
    "iteration_root", "bin", "bin_low", "bin_high", "nominal_value",
    "nominal_stat_abs", "nominal_stat_pct",
}

GROUP_ORDER = [
    "iteration",
    "binning",
    "priorshape",
    "sweight",
    "jetsreco",
    "d0meson",
    "stat",
    "other",
]

GROUP_PREFIXES = {
    "binning_": "binning",
    "priorshape_": "priorshape",
    "first_variable_": "priorshape",
    "second_variable_": "priorshape",
    "sweight_": "sweight",
    "jetsreco_": "jetsreco",
    "d0meson_": "d0meson",
    "iteration_": "iteration",
}

PRETTY_GROUP = {
    "iteration": "Iteration",
    "binning": "Binning",
    "priorshape": "Prior shape",
    "sweight": "sWeight",
    "jetsreco": "Jets reco",
    "d0meson": "D0 meson",
    "stat": "Stat.",
    "other": "Other",
}

# Compact labels for common components.
PRETTY_LABELS = {
    "iteration_minus_1_pct": "Iter. -1",
    "iteration_plus_1_pct": "Iter. +1",
    "binning_binning_first_variable_rms_pct": "Binning: jet pT RMS",
    "binning_binning_second_variable_rms_pct": "Binning: 2nd var. RMS",
    "first_variable_plus20_pct": "Prior: jet pT +20%",
    "first_variable_minus20_pct": "Prior: jet pT -20%",
    "second_variable_plus20_pct": "Prior: 2nd var. +20%",
    "second_variable_minus20_pct": "Prior: 2nd var. -20%",
    "sweight_cheby2_bkg_pct": "sW: Cheby2 bkg",
    "sweight_double_gauss_sig_pct": "sW: double Gauss",
    "sweight_student_t_sig_pct": "sW: Student t",
    "sweight_narrow_fit_range_pct": "sW: narrow range",
    "sweight_wide_fit_range_pct": "sW: wide range",
    "sweight_zeroed_negative_bins_pct": "sW: zero neg. bins",
    "jetsreco_jet_rec_efficiency_pct": "Jet reco: track eff.",
    "jetsreco_jet_nhitsfit13_pct": "Jet reco: nHitsFit13",
    "jetsreco_jet_nhitsfit17_pct": "Jet reco: nHitsFit17",
    "jetsreco_jet_ktdrop_pct": "Jet reco: kT dropped",
    "jetsreco_jet_dca2_8_pct": "Jet reco: DCA 2.8",
    "jetsreco_jet_dca3_2_pct": "Jet reco: DCA 3.2",
    "jetsreco_jet_hadroniccorr_pct": "Jet reco: hadr. corr.",
    "d0meson_tpc_track_up_pct": "D0: TPC track up",
    "d0meson_tpc_track_down_pct": "D0: TPC track down",
    "d0meson_pid_up_pct": "D0: PID up",
    "d0meson_pid_down_pct": "D0: PID down",
    "d0meson_single_track_pt_up_pct": "D0: single track up",
    "d0meson_single_track_pt_down_pct": "D0: single track down",
    "d0meson_topo_eff_150pct_pct": "D0: topo eff. 150%",
    "d0meson_topo_eff_50pct_pct": "D0: topo eff. 50%",
    "d0meson_double_counting_up_pct": "D0: double count up",
    "d0meson_double_counting_down_pct": "D0: double count down",
    "d0meson_vertex_correction_up_pct": "D0: vertex corr. up",
    "d0meson_vertex_correction_down_pct": "D0: vertex corr. down",
    "d0meson_secondary_track_up_pct": "D0: secondary track up",
    "d0meson_secondary_track_down_pct": "D0: secondary track down",
    "nominal_stat_pct": "Stat.",
}



# Large palette + independent markers/linestyles.  With the default matplotlib
# cycle (10 colors) many components get identical colors; this keeps colors
# unique for the usual ~30--40 systematic components and also distinguishes
# categories by line style.
COLOR_MAPS = ("tab20", "tab20b", "tab20c")
MARKERS = ("o", "s", "^", "v", "D", "P", "X", "h", "<", ">", "p", "*")
LINESTYLE_BY_GROUP = {
    "iteration": "-",
    "binning": "--",
    "priorshape": "-.",
    "sweight": ":",
    "jetsreco": (0, (5, 1)),
    "d0meson": (0, (3, 1, 1, 1)),
    "stat": "-",
    "other": "-",
}

# Component sets used for the extra category plots.
# "Set" here follows the overview slide structure:
#   Unfolding, sPlot, D0 meson, Jets reconstruction.
# The individual lines inside each set remain separate columns/variations
# (e.g. up/down, +/-20%, different fit models, etc.).
# The JetsReco ICS-vs-SEAM component is intentionally excluded.
GLOBAL_SET_EXCLUDE_PATTERNS = [
    r"seam",
    r"ics[_ -]*vs[_ -]*seam",
    r"seam[_ -]*vs[_ -]*ics",
]

COMPONENT_SETS = [
    {
        "key": "unfolding",
        "title": "Unfolding",
        "patterns": [
            r"^binning_.*_pct$",
            r"^iteration_.*_pct$",
            r"^first[_-]?variable_.*_pct$",
            r"^second[_-]?variable_.*_pct$",
            r"^priorshape_.*_pct$",
        ],
    },
    {
        "key": "splot",
        "title": "sPlot",
        "patterns": [r"^sweight_.*_pct$"],
    },
    {
        "key": "d0meson",
        "title": r"D0 meson",
        "patterns": [r"^d0meson_.*_pct$"],
    },
    {
        "key": "jetsreco",
        "title": "Jets reconstruction",
        "patterns": [r"^jetsreco_.*_pct$"],
    },
]

# Optional aliases, so both short and slide-like names work with --sets.
COMPONENT_SET_ALIASES = {
    "unfold": "unfolding",
    "unfolding": "unfolding",
    "splot": "splot",
    "s_plot": "splot",
    "d0": "d0meson",
    "d0meson": "d0meson",
    "d0_meson": "d0meson",
    "jets": "jetsreco",
    "jetsreco": "jetsreco",
    "jets_reco": "jetsreco",
    "jets_reconstruction": "jetsreco",
}

COMPONENT_SET_BY_KEY = {s["key"]: s for s in COMPONENT_SETS}


def distinct_colors(n: int):
    """Return at least n visually distinct matplotlib colors."""
    colors = []
    for cmap_name in COLOR_MAPS:
        cmap = plt.get_cmap(cmap_name)
        if hasattr(cmap, "colors"):
            colors.extend(list(cmap.colors))
        else:
            colors.extend([cmap(i / max(1, cmap.N - 1)) for i in range(cmap.N)])

    # tab20-like palettes contain neighbouring light/dark pairs.  This reorder
    # avoids placing the most similar shades next to each other in the legend.
    colors = colors[::2] + colors[1::2]

    if not colors:
        colors = [None]
    return [colors[i % len(colors)] for i in range(n)]

def parse_float(x: object) -> float:
    try:
        s = str(x).strip()
        if not s or s.lower() in {"nan", "none", "<na>"}:
            return float("nan")
        return float(s)
    except Exception:
        return float("nan")


def finite(x: float) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def sanitize_filename(text: str) -> str:
    text = str(text)
    text = text.replace("%", "pct")
    text = re.sub(r"[^A-Za-z0-9_.+-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "plot"


def component_group(col: str) -> str:
    if col == "nominal_stat_pct":
        return "stat"
    for prefix, group in GROUP_PREFIXES.items():
        if col.startswith(prefix):
            return group
    return "other"


def pretty_label(col: str) -> str:
    if col in PRETTY_LABELS:
        return PRETTY_LABELS[col]
    label = col
    label = re.sub(r"_pct$", "", label)
    label = label.replace("priorshape_", "Prior: ")
    label = label.replace("sweight_", "sW: ")
    label = label.replace("jetsreco_", "Jet reco: ")
    label = label.replace("d0meson_", "D0: ")
    label = label.replace("binning_", "Binning: ")
    label = label.replace("_", " ")
    label = label.replace(" plus20", " +20%")
    label = label.replace(" minus20", " -20%")
    return label


def is_component_column(col: str, include_stat: bool) -> bool:
    if not col.endswith("_pct"):
        return False
    if col in BASE_COLUMNS:
        return include_stat and col == "nominal_stat_pct"
    if "max_abs_pct" in col:
        return False
    if col.endswith("_stat_pct"):
        return include_stat and col == "nominal_stat_pct"
    return True


def read_wide_table(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        header = reader.fieldnames or []
    return header, rows


def select_component_columns(header: List[str], include_stat: bool, groups: Iterable[str], regex: str) -> List[str]:
    allowed_groups = {g.strip().lower() for g in groups if g.strip()}
    pattern = re.compile(regex) if regex else None
    cols = []
    for col in header:
        if not is_component_column(col, include_stat=include_stat):
            continue
        group = component_group(col)
        if allowed_groups and group not in allowed_groups:
            continue
        if pattern and not pattern.search(col):
            continue
        cols.append(col)

    def sort_key(c: str):
        g = component_group(c)
        gi = GROUP_ORDER.index(g) if g in GROUP_ORDER else len(GROUP_ORDER)
        return gi, c

    return sorted(cols, key=sort_key)


def regex_matches_any(text: str, patterns: Iterable[str]) -> bool:
    return any(re.search(pattern, text, flags=re.IGNORECASE) for pattern in patterns)


def component_set_columns(component_cols: List[str], set_def: Dict[str, object]) -> List[str]:
    patterns = [str(p) for p in set_def.get("patterns", [])]
    exclude = [str(p) for p in GLOBAL_SET_EXCLUDE_PATTERNS]
    exclude.extend(str(p) for p in set_def.get("exclude", []))
    out = []
    for col in component_cols:
        if patterns and not regex_matches_any(col, patterns):
            continue
        if exclude and regex_matches_any(col, exclude):
            continue
        out.append(col)
    return out


def selected_component_sets(raw: str) -> List[Dict[str, object]]:
    txt = str(raw or "all").strip()
    if not txt or txt.lower() in {"all", "*"}:
        return list(COMPONENT_SETS)

    requested = [x.strip() for x in txt.split(",") if x.strip()]
    selected = []
    unknown = []
    for key in requested:
        canonical = COMPONENT_SET_ALIASES.get(key, key)
        if canonical in COMPONENT_SET_BY_KEY:
            selected.append(COMPONENT_SET_BY_KEY[canonical])
        else:
            unknown.append(key)
    if unknown:
        print(f"[warning] Unknown component set(s): {', '.join(unknown)}", file=sys.stderr)
        print("[warning] Available sets: " + ", ".join(s["key"] for s in COMPONENT_SETS), file=sys.stderr)
    return selected


def print_available_sets() -> None:
    for s in COMPONENT_SETS:
        print(f"{s['key']}	{s['title']}")


def group_rows(rows: List[Dict[str, str]]) -> Dict[Tuple[str, str], List[Dict[str, str]]]:
    groups: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (row.get("job_label", ""), row.get("centrality", ""))
        groups[key].append(row)
    for key in groups:
        groups[key].sort(key=lambda r: parse_float(r.get("bin", "nan")))
    return dict(groups)


def active_columns(rows: List[Dict[str, str]], component_cols: List[str], mode: str, min_abs: float) -> List[str]:
    active = []
    for col in component_cols:
        vals = []
        for r in rows:
            v = parse_float(r.get(col, ""))
            if finite(v):
                if mode == "abs":
                    v = abs(v)
                vals.append(v)
        if not vals:
            continue
        if max(abs(v) for v in vals) <= min_abs:
            continue
        active.append(col)
    return active


def x_values_and_labels(rows: List[Dict[str, str]], x_mode: str) -> Tuple[List[float], List[str], str]:
    bins = [parse_float(r.get("bin", "nan")) for r in rows]
    lows = [parse_float(r.get("bin_low", "nan")) for r in rows]
    highs = [parse_float(r.get("bin_high", "nan")) for r in rows]
    labels = []
    for lo, hi, b in zip(lows, highs, bins):
        if finite(lo) and finite(hi):
            labels.append(f"{lo:g}-{hi:g}")
        else:
            labels.append(f"{int(b)}" if finite(b) else "")
    if x_mode == "center" and all(finite(lo) and finite(hi) for lo, hi in zip(lows, highs)):
        return [0.5 * (lo + hi) for lo, hi in zip(lows, highs)], labels, "bin center"
    return bins, labels, "bin"


def make_title(rows: List[Dict[str, str]]) -> str:
    r0 = rows[0]
    job = r0.get("job_label", "")
    obs = r0.get("observable_pretty", "") or r0.get("observable", "")
    cent = r0.get("centrality_label", "") or r0.get("centrality", "")
    return f"{job}: {obs}, {cent}"


def plot_one(
    rows: List[Dict[str, str]],
    cols: List[str],
    out_png: Path,
    pdf: PdfPages,
    mode: str,
    x_mode: str,
    ymax: float,
    legend_cols: int,
    dpi: int,
    title_suffix: str = "",
) -> None:
    x, xlabels, xlabel = x_values_and_labels(rows, x_mode)
    title = make_title(rows)
    if title_suffix:
        title = f"{title}\n{title_suffix}"

    ncols = max(1, int(legend_cols))
    fig_width = 13.5
    fig_height = 7.5 if len(cols) <= 18 else 8.8
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    colors = distinct_colors(len(cols))

    for i, col in enumerate(cols):
        y = []
        for r in rows:
            v = parse_float(r.get(col, ""))
            if finite(v) and mode == "abs":
                v = abs(v)
            y.append(v if finite(v) else float("nan"))

        group = component_group(col)
        ax.plot(
            x, y,
            color=colors[i],
            linestyle=LINESTYLE_BY_GROUP.get(group, "-"),
            marker=MARKERS[i % len(MARKERS)],
            linewidth=1.45,
            markersize=3.8,
            alpha=0.95,
            label=pretty_label(col),
        )

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("|contribution| [%]" if mode == "abs" else "contribution [%]")
    ax.grid(True, alpha=0.3)
    if mode == "signed":
        ax.axhline(0.0, linewidth=1.0)
    else:
        ax.set_ylim(bottom=0.0)
    if ymax > 0:
        if mode == "signed":
            ax.set_ylim(-ymax, ymax)
        else:
            ax.set_ylim(0.0, ymax)

    if len(x) <= 14:
        ax.set_xticks(x)
        ax.set_xticklabels(xlabels, rotation=35, ha="right")

    if cols:
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0,
            fontsize=7.5,
            ncol=ncols,
            handlelength=2.8,
            handletextpad=0.6,
            columnspacing=1.0,
        )
    fig.tight_layout(rect=[0.0, 0.0, 0.76 if ncols == 1 else 0.70, 1.0])

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi)
    pdf.savefig(fig)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot individual systematic contributions from wide TSV table.")
    p.add_argument("--input", default="Systematics/systematics_components_wide_v3.tsv", help="Wide TSV table.")
    p.add_argument("--out-dir", default="Systematics/Plots/SystematicsComponents", help="Output directory.")
    p.add_argument("--mode", choices=["abs", "signed"], default="abs", help="Plot absolute or signed percentages.")
    p.add_argument("--x-mode", choices=["center", "bin"], default="center", help="Use bin center or bin index on x-axis.")
    p.add_argument("--groups", default="iteration,binning,priorshape,sweight,jetsreco,d0meson", help="Comma-separated component groups to plot. Empty = all.")
    p.add_argument("--regex", default="", help="Only plot component columns matching this regex.")
    p.add_argument("--include-stat", action="store_true", help="Also plot nominal_stat_pct as a component.")
    p.add_argument("--min-abs", type=float, default=0.0, help="Skip components whose max absolute value is <= this threshold.")
    p.add_argument("--ymax", type=float, default=0.0, help="Manual y-axis max. 0 = automatic.")
    p.add_argument("--legend-cols", type=int, default=1, help="Number of legend columns.")
    p.add_argument("--dpi", type=int, default=180, help="PNG DPI.")
    p.add_argument("--no-png", action="store_true", help="Only write the multipage PDF.")
    p.add_argument("--no-sets", action="store_true", help="Do not write the additional per-systematic-set plots.")
    p.add_argument("--sets", default="all", help="Comma-separated category-set keys to plot: unfolding,splot,d0meson,jetsreco. Use 'all' for all four.")
    p.add_argument("--list-sets", action="store_true", help="Print available component-set keys and exit.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.list_sets:
        print_available_sets()
        return 0

    input_path = Path(args.input)
    out_dir = Path(args.out_dir)
    png_dir = out_dir / "png"
    mode_tag = "abs" if args.mode == "abs" else "signed"
    pdf_path = out_dir / f"all_systematics_components_{mode_tag}.pdf"

    if not input_path.exists():
        print(f"[error] Missing input file: {input_path}", file=sys.stderr)
        return 1

    header, rows = read_wide_table(input_path)
    if not rows:
        print(f"[error] No rows in {input_path}", file=sys.stderr)
        return 1

    groups = [g.strip().lower() for g in args.groups.split(",") if g.strip()] if args.groups.strip() else []
    component_cols = select_component_columns(header, include_stat=args.include_stat, groups=groups, regex=args.regex)
    grouped = group_rows(rows)

    out_dir.mkdir(parents=True, exist_ok=True)
    if not args.no_png:
        png_dir.mkdir(parents=True, exist_ok=True)

    n_plots = 0
    skipped = 0
    with PdfPages(pdf_path) as pdf:
        for (job, cent), gr in sorted(grouped.items(), key=lambda kv: (kv[0][0], parse_float(kv[0][1]))):
            cols = active_columns(gr, component_cols, args.mode, args.min_abs)
            if not cols:
                skipped += 1
                continue
            cent_label = gr[0].get("centrality_label", cent)
            base = sanitize_filename(f"{job}_cent{cent}_{cent_label}_{mode_tag}")
            out_png = png_dir / f"{base}.png" if not args.no_png else out_dir / "_dummy.png"
            plot_one(
                rows=gr,
                cols=cols,
                out_png=out_png,
                pdf=pdf,
                mode=args.mode,
                x_mode=args.x_mode,
                ymax=args.ymax,
                legend_cols=args.legend_cols,
                dpi=args.dpi,
                title_suffix="All individual components",
            )
            n_plots += 1

    n_set_plots = 0
    n_set_pages = 0
    if not args.no_sets:
        set_root = out_dir / "by_set"
        set_root.mkdir(parents=True, exist_ok=True)
        selected_sets = selected_component_sets(args.sets)
        for set_def in selected_sets:
            set_key = str(set_def["key"])
            set_title = str(set_def["title"])
            set_cols_all = component_set_columns(component_cols, set_def)
            if not set_cols_all:
                print(f"[warning] Set '{set_key}' has no matching columns in the input table.", file=sys.stderr)
                continue

            set_pdf_path = set_root / f"{set_key}_{mode_tag}.pdf"
            set_png_dir = set_root / "png" / set_key
            if not args.no_png:
                set_png_dir.mkdir(parents=True, exist_ok=True)

            set_pages = 0
            set_skipped = 0
            with PdfPages(set_pdf_path) as pdf:
                for (job, cent), gr in sorted(grouped.items(), key=lambda kv: (kv[0][0], parse_float(kv[0][1]))):
                    cols = active_columns(gr, set_cols_all, args.mode, args.min_abs)
                    if not cols:
                        set_skipped += 1
                        continue
                    cent_label = gr[0].get("centrality_label", cent)
                    base = sanitize_filename(f"{set_key}_{job}_cent{cent}_{cent_label}_{mode_tag}")
                    out_png = set_png_dir / f"{base}.png" if not args.no_png else set_root / "_dummy.png"
                    plot_one(
                        rows=gr,
                        cols=cols,
                        out_png=out_png,
                        pdf=pdf,
                        mode=args.mode,
                        x_mode=args.x_mode,
                        ymax=args.ymax,
                        legend_cols=args.legend_cols,
                        dpi=args.dpi,
                        title_suffix=set_title,
                    )
                    set_pages += 1

            n_set_plots += 1
            n_set_pages += set_pages
            print(f"[info] Set PDF: {set_pdf_path} pages={set_pages}, skipped={set_skipped}")

    print(f"[info] Input rows: {len(rows)}")
    print(f"[info] Component columns available: {len(component_cols)}")
    print(f"[info] All-component plots written: {n_plots}")
    if not args.no_sets:
        print(f"[info] Set PDFs written: {n_set_plots}")
        print(f"[info] Set pages written: {n_set_pages}")
    print(f"[info] Skipped empty groups: {skipped}")
    print(f"[info] PDF: {pdf_path}")
    if not args.no_png:
        print(f"[info] PNG dir: {png_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
