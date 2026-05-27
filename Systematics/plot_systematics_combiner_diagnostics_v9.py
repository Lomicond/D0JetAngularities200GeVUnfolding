#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_systematics_combiner_diagnostics_v7.py

Diagnostic plots for the final systematics combiner output.

Version v7: adds presentation-style colors for main sets, unfolding breakdown, and JetsReco Barlow diagnostics. Based on v6; fixes JetsReco Barlow diagnostics. Main-set components (sPlot, D0 meson, Jets reco) can be resolved
either from one aggregate column (e.g. sweight_sym_pct) or from several
sub-component columns (e.g. sweight_bkg_model_sym_pct, ...), which are then
combined in quadrature for the diagnostic plot. Barlow diagnostics now include split GUI group names (e.g. jetsreco_jet_dca), compute the threshold directly from stat errors when possible, and set the y-axis after plotting so percent values are not clipped to 0--1.

Inputs:
  1) final_systematics_gui_v*.tsv or final_systematics_v*.tsv
     produced by systematics_final_combiner_v*.py / GUI config workflow

  2) optional final_systematics_gui_barlow_details_v*.tsv
     or final_systematics_barlow_details_v*.tsv
     produced by systematics_final_combiner_v*.py when write_details=true

Outputs:
  * set subtotal / breakdown plots:
      Unfolding, sPlot, D0 meson, Jets reco, and their quadrature sum

  * optional unfolding-internal breakdown plots:
      iteration, binning 1st/2nd, prior shape 1st/2nd, and their quadrature sum

  * JetsReco Barlow diagnostics, if the details TSV is available:
      before Barlow, Barlow threshold, after Barlow

Example:
  python3 Systematics/plot_systematics_combiner_diagnostics_v9.py \
    --final Systematics/final_systematics_v9.tsv \
    --details Systematics/final_systematics_barlow_details_v9.tsv

If --details is omitted, the script tries to infer the details file from the
version number in --final.  If it is not found, only subtotal/breakdown plots
are produced.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# -----------------------------------------------------------------------------
# User-facing labels and grouping definitions
# -----------------------------------------------------------------------------

# Main sets from the summary slide.
# The last entry is computed here as a quadrature sum of the four sets.
SET_COMPONENTS = [
    ("unfolding", "Unfolding"),
    ("sweight", "sPlot"),
    ("d0meson", r"D$^{0}$ meson"),
    ("jetsreco", "Jets reco (Barlow)"),
    ("set_quadrature_sum", r"Final systematic"),
]

# Internal decomposition of the Unfolding set. This is kept as a separate,
# optional diagnostic so it is not confused with the four main sets.
UNFOLDING_BREAKDOWN_COMPONENTS = [
    ("binning_first_variable", "Binning - 1st var"),
    ("binning_second_variable", "Binning - 2nd var"),
    ("iteration", "Iterations"),
    ("prior_shape_first_variable", "Prior shape - 1st var"),
    ("prior_shape_second_variable", "Prior shape - 2nd var"),
    ("unfolding", r"Unfolding quadrature sum"),
]

# Presentation-style colors.  Change these dictionaries if you want to tune
# the diagnostic plots to match your slides.
MAIN_SET_COLORS = {
    # keys used internally
    "unfolding": "#d62728",              # red
    "sweight": "#0000ff",               # blue
    "d0meson": "#008f3a",               # green
    "jetsreco": "#f5a623",              # orange
    "set_quadrature_sum": "black",
    # labels shown in the legend, kept for convenience
    "Unfolding": "#d62728",
    "sPlot": "#0000ff",
    r"D$^{0}$ meson": "#008f3a",
    "D0 meson": "#008f3a",
    "Jets reco": "#f5a623",
    r"Final Systematic sum": "black",
    "Quadrature sum": "black",
}

UNFOLDING_COLORS = {
    "binning_first_variable": "#d62728",
    "binning_second_variable": "#d62728",
    "iteration": "#d62728",
    "prior_shape_first_variable": "#d62728",
    "prior_shape_second_variable": "#d62728",
    "unfolding": "black",
    "Binning - 1st var": "#d62728",
    "Binning - 2nd var": "#d62728",
    "Iterations": "#d62728",
    "Prior shape - 1st var": "#d62728",
    "Prior shape - 2nd var": "#d62728",
    "Unfolding quadrature sum": "black",
}

JETSRECO_COLORS = {
    "Tracking eff.": "#1f77b4",
    "nHitsFit": "#9467bd",
    "nHitsFit 13": "#9467bd",
    "nHitsFit 17": "#9467bd",
    "kT dropped": "#2ca02c",
    "kT-Dropped": "#2ca02c",
    "DCA": "#d62728",
    "DCA 2.8": "#d62728",
    "DCA 3.2": "#d62728",
    "Hadr. correction": "#ff7f0e",
}

UNFOLDING_PARTS = [
    "binning_first_variable",
    "binning_second_variable",
    "iteration",
    "prior_shape_first_variable",
    "prior_shape_second_variable",
]

# Column-name aliases.  GUI configs can rename groups, so the output columns
# do not always have the same prefix as the built-in template.  The plotter
# resolves each logical component to the first matching *_<side>_pct column.
COMPONENT_ALIASES = {
    "sweight": ["sweight", "splot", "s_plot"],
    "d0meson": ["d0meson", "d0_meson", "d0", "d0meson_sys"],
    "jetsreco": ["jetsreco", "jets_reco", "jet_reco", "jet_reconstruction", "jets_reconstruction"],
    "iteration": ["iteration", "iterations"],
    "binning_first_variable": ["binning_first_variable", "binning_1st_variable", "binning_1st_var", "binning_first_var"],
    "binning_second_variable": ["binning_second_variable", "binning_2nd_variable", "binning_2nd_var", "binning_second_var"],
    "prior_shape_first_variable": ["prior_shape_first_variable", "prior_shape_1st_variable", "prior_shape_1st_var", "prior_shape_first_var"],
    "prior_shape_second_variable": ["prior_shape_second_variable", "prior_shape_2nd_variable", "prior_shape_2nd_var", "prior_shape_second_var"],
}

# These logical main sets may be represented either as one aggregate column
# (e.g. jetsreco_sym_pct) or as several GUI-created sub-components
# (e.g. jetsreco_jet_dca_sym_pct, jetsreco_jet_nhitsfit_sym_pct, ...).
# When sub-components are present, the diagnostic plot combines them in
# quadrature and does not double-count the aggregate column.
AGGREGATE_SET_COMPONENTS = {"sweight", "d0meson", "jetsreco"}

PRETTY_VARIATION = {
    "jet_rec_efficiency": "Tracking eff.",
    "jet_recoefficiency": "Tracking eff.",
    "jet_tracking_eff": "Tracking eff.",
    "jet_nhitsfit13": "nHitsFit 13",
    "jet_nhitsfit_13": "nHitsFit 13",
    "jet_nhitsfit17": "nHitsFit 17",
    "jet_nhitsfit_17": "nHitsFit 17",
    "jet_nhitsfit": "nHitsFit",
    "jet_ktdrop": "kT dropped",
    "jet_kt_dropped": "kT dropped",
    "jet_dca2_8": "DCA 2.8",
    "jet_dca_2_8": "DCA 2.8",
    "jet_dca3_2": "DCA 3.2",
    "jet_dca_3_2": "DCA 3.2",
    "jet_dca": "DCA",
    "jet_hadroniccorr": "Hadr. correction",
    "jet_hadronic_corr": "Hadr. correction",
}

# JetsReco variation-code map used by several scan/summary formats.  This is
# only used as a fallback when the name itself is too generic, e.g. "jet_dca".
JETSRECO_CODE_LABELS = {
    "21": "Tracking eff.",
    "22": "nHitsFit 13",
    "23": "nHitsFit 17",
    "24": "kT dropped",
    "25": "DCA 2.8",
    "26": "DCA 3.2",
    "27": "Hadr. correction",
}

MARKERS = ("o", "s", "^", "v", "D", "P", "X", "h", "<", ">", "p", "*")


def color_for_component(key: str, label: str) -> Optional[str]:
    """Color for subtotal and unfolding-breakdown plots."""
    if key in MAIN_SET_COLORS:
        return MAIN_SET_COLORS[key]
    if label in MAIN_SET_COLORS:
        return MAIN_SET_COLORS[label]
    if key in UNFOLDING_COLORS:
        return UNFOLDING_COLORS[key]
    if label in UNFOLDING_COLORS:
        return UNFOLDING_COLORS[label]
    return None


def color_for_jetsreco_variation(label: str) -> Optional[str]:
    """Use the same color for paired JetsReco variations."""
    if label in JETSRECO_COLORS:
        return JETSRECO_COLORS[label]
    txt = str(label).lower()
    if "nhits" in txt:
        return JETSRECO_COLORS.get("nHitsFit")
    if "dca" in txt:
        return JETSRECO_COLORS.get("DCA")
    if "tracking" in txt or "eff" in txt:
        return JETSRECO_COLORS.get("Tracking eff.")
    if "kt" in txt and "drop" in txt:
        return JETSRECO_COLORS.get("kT dropped")
    if "hadr" in txt or "hadronic" in txt:
        return JETSRECO_COLORS.get("Hadr. correction")
    return None


# -----------------------------------------------------------------------------
# Basic helpers
# -----------------------------------------------------------------------------

def warn(msg: str) -> None:
    print(f"[warning] {msg}", file=sys.stderr)


def info(msg: str) -> None:
    print(f"[info] {msg}")


def parse_float(x: Any, default: float = math.nan) -> float:
    try:
        if x is None:
            return default
        s = str(x).strip()
        if not s or s.lower() in {"nan", "none", "<na>"}:
            return default
        return float(s)
    except Exception:
        return default


def finite(x: Any) -> bool:
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def safe_name(text: Any) -> str:
    s = str(text or "").strip()
    s = s.replace("%", "pct")
    s = re.sub(r"[^A-Za-z0-9_.+-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "plot"


def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        header = reader.fieldnames or []
    return header, rows


def print_component_column_diagnostics(header: Sequence[str], side: str) -> None:
    suffix = f"_{side}_pct"
    cols = [c for c in header if c.endswith(suffix)]
    info(f"Available *{suffix} columns: " + (", ".join(cols) if cols else "none"))
    if cols:
        dummy = {c: "" for c in header}
        for comp, _ in SET_COMPONENTS:
            if comp == "set_quadrature_sum":
                continue
            if comp == "unfolding":
                parts = [resolve_component_col(dummy, part, side) for part in UNFOLDING_PARTS]
                parts = [p for p in parts if p]
                info("Resolved unfolding: computed from " + (", ".join(parts) if parts else "NO PARTS FOUND"))
                continue
            resolved_cols = resolve_component_cols(dummy, comp, side)
            if not resolved_cols:
                info(f"Resolved {comp}: NOT FOUND")
            elif len(resolved_cols) == 1:
                info(f"Resolved {comp}: {resolved_cols[0]}")
            else:
                info(f"Resolved {comp}: quadrature of " + ", ".join(resolved_cols))

def group_by_job_cent(rows: Sequence[Dict[str, str]]) -> Dict[Tuple[str, str], List[Dict[str, str]]]:
    groups: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for r in rows:
        groups[(str(r.get("job_label", "")), str(r.get("centrality", "")))].append(r)
    for key in groups:
        groups[key].sort(key=row_sort_key)
    return dict(groups)


def row_sort_key(row: Dict[str, str]) -> Tuple[float, float]:
    b = parse_float(row.get("bin"), math.nan)
    lo = parse_float(row.get("bin_low"), math.nan)
    return (b if finite(b) else 1e9, lo if finite(lo) else 0.0)


def x_values_and_labels(rows: Sequence[Dict[str, str]], x_mode: str) -> Tuple[np.ndarray, List[str], str]:
    xs: List[float] = []
    labels: List[str] = []
    for i, r in enumerate(rows, start=1):
        b = parse_float(r.get("bin"), math.nan)
        lo = parse_float(r.get("bin_low"), math.nan)
        hi = parse_float(r.get("bin_high"), math.nan)
        if x_mode == "center" and finite(lo) and finite(hi):
            xs.append(0.5 * (lo + hi))
        elif finite(b):
            xs.append(b)
        else:
            xs.append(float(i))

        if finite(lo) and finite(hi):
            labels.append(f"{lo:g}-{hi:g}")
        elif finite(b):
            labels.append(f"{int(b)}")
        else:
            labels.append(str(i))
    xlabel = "bin center" if x_mode == "center" else "bin"
    return np.asarray(xs, dtype=float), labels, xlabel


def make_title(rows: Sequence[Dict[str, str]], suffix: str = "") -> str:
    r0 = rows[0]
    job = str(r0.get("job_label", ""))
    obs = str(r0.get("observable_pretty") or r0.get("observable") or job)
    cent = str(r0.get("centrality_label") or r0.get("centrality") or "")
    title = f"{job}: {obs}, {cent}"
    if suffix:
        title += f"\n{suffix}"
    return title


def autoset_xticks(ax: Any, x: np.ndarray, labels: List[str]) -> None:
    if len(x) <= 14:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right")


def component_col(component: str, side: str) -> str:
    return f"{component}_{side}_pct"


def possible_component_cols(component: str, side: str) -> List[str]:
    prefixes = COMPONENT_ALIASES.get(component, [component])
    out = []
    for prefix in prefixes:
        out.append(f"{prefix}_{side}_pct")
    # Also try the literal component name last, in case the alias list was incomplete.
    literal = component_col(component, side)
    if literal not in out:
        out.append(literal)
    return out


def resolve_component_col(row: Dict[str, str], component: str, side: str) -> Optional[str]:
    """Resolve a component represented by exactly one column."""
    for col in possible_component_cols(component, side):
        if col in row:
            return col
    return None


def _alias_prefixes(component: str) -> List[str]:
    prefixes = COMPONENT_ALIASES.get(component, [component])
    literal = component
    if literal not in prefixes:
        prefixes.append(literal)
    return prefixes


def resolve_component_cols(row: Dict[str, str], component: str, side: str) -> List[str]:
    """Resolve columns used for one logical component.

    For normal components this returns either [component_side_pct] or [].
    For the main sets sweight/d0meson/jetsreco it also supports GUI split
    columns such as sweight_bkg_model_sym_pct. If split columns are present,
    they are returned and the aggregate column is not included, preventing
    double counting. If no split columns exist, the aggregate column is used.
    """
    if component not in AGGREGATE_SET_COMPONENTS:
        col = resolve_component_col(row, component, side)
        return [col] if col else []

    suffix = f"_{side}_pct"
    exact = [c for c in possible_component_cols(component, side) if c in row]

    split: List[str] = []
    for col in row.keys():
        if not col.endswith(suffix):
            continue
        if col in exact:
            continue
        for prefix in _alias_prefixes(component):
            if col.startswith(prefix + "_"):
                split.append(col)
                break

    # Deterministic order, no duplicates.
    split = sorted(dict.fromkeys(split))
    exact = sorted(dict.fromkeys(exact))
    return split if split else exact


def component_row_set_value(row: Dict[str, str], component: str, side: str) -> float:
    cols = resolve_component_cols(row, component, side)
    if not cols:
        return math.nan
    if len(cols) == 1:
        return row_value(row, cols[0])
    return quadrature(row_value(row, c) for c in cols)

def row_value(row: Dict[str, str], col: str) -> float:
    return parse_float(row.get(col, ""), math.nan)


def component_row_value(row: Dict[str, str], component: str, side: str) -> float:
    return component_row_set_value(row, component, side)


def available_component(row: Dict[str, str], component: str, side: str) -> bool:
    return bool(resolve_component_cols(row, component, side))


def quadrature(values: Iterable[float]) -> float:
    vals = [float(v) for v in values if finite(v)]
    if not vals:
        return math.nan
    return math.sqrt(sum(v * v for v in vals))


def component_values(rows: Sequence[Dict[str, str]], component: str, side: str) -> np.ndarray:
    """Return one component as an array in percent.

    component="unfolding" is constructed from the five unfolding
    subcomponents.  component="set_quadrature_sum" is the quadrature sum of
    the four main summary-slide sets: Unfolding, sPlot, D0 meson, Jets reco.
    """
    vals: List[float] = []
    if component == "unfolding":
        for r in rows:
            vals.append(quadrature(component_row_value(r, part, side) for part in UNFOLDING_PARTS))
        return np.asarray(vals, dtype=float)

    if component == "set_quadrature_sum":
        for r in rows:
            unfolding = quadrature(component_row_value(r, part, side) for part in UNFOLDING_PARTS)
            vals.append(quadrature([
                unfolding,
                component_row_value(r, "sweight", side),
                component_row_value(r, "d0meson", side),
                component_row_value(r, "jetsreco", side),
            ]))
        return np.asarray(vals, dtype=float)

    vals = [component_row_set_value(r, component, side) for r in rows]
    return np.asarray(vals, dtype=float)


def component_exists(rows: Sequence[Dict[str, str]], component: str, side: str) -> bool:
    if not rows:
        return False
    if component == "unfolding":
        return any(resolve_component_col(rows[0], part, side) is not None for part in UNFOLDING_PARTS)
    if component == "set_quadrature_sum":
        return (component_exists(rows, "unfolding", side)
                and component_exists(rows, "sweight", side)
                and component_exists(rows, "d0meson", side)
                and component_exists(rows, "jetsreco", side))
    return bool(resolve_component_cols(rows[0], component, side))


# -----------------------------------------------------------------------------
# Subtotal and breakdown plotting
# -----------------------------------------------------------------------------

def plot_component_collection(
    rows: Sequence[Dict[str, str]],
    components: Sequence[Tuple[str, str]],
    out_png: Optional[Path],
    pdf: PdfPages,
    side: str,
    x_mode: str,
    ylabel: str,
    title_suffix: str,
    ymax: float,
    dpi: int,
) -> bool:
    if not rows:
        return False

    x, xlabels, xlabel = x_values_and_labels(rows, x_mode)
    fig, ax = plt.subplots(figsize=(10.5, 6.2))

    n_drawn = 0
    for i, (key, label) in enumerate(components):
        if not component_exists(rows, key, side):
            continue
        y = component_values(rows, key, side)
        if not np.any(np.isfinite(y)):
            continue
        ax.plot(
            x,
            y,
            marker=MARKERS[i % len(MARKERS)],
            linewidth=1.6,
            markersize=4.2,
            label=label,
            color=color_for_component(key, label),
        )
        n_drawn += 1

    if n_drawn == 0:
        plt.close(fig)
        return False

    ax.set_title(make_title(rows, title_suffix))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.28)
    ax.set_ylim(bottom=0.0)
    if ymax > 0:
        ax.set_ylim(0.0, ymax)
    autoset_xticks(ax, x, xlabels)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0, fontsize=8)
    fig.tight_layout(rect=[0.0, 0.0, 0.78, 1.0])

    if out_png is not None:
        out_png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_png, dpi=dpi)
    pdf.savefig(fig)
    plt.close(fig)
    return True


def plot_subtotals(
    rows: List[Dict[str, str]],
    out_dir: Path,
    side: str,
    x_mode: str,
    ymax: float,
    dpi: int,
    write_png: bool,
) -> int:
    grouped = group_by_job_cent(rows)
    pdf_path = out_dir / f"subtotals_{side}.pdf"
    png_dir = out_dir / "png" / f"subtotals_{side}"
    n_pages = 0
    with PdfPages(pdf_path) as pdf:
        for (job, cent), gr in sorted(grouped.items(), key=lambda kv: (kv[0][0], parse_float(kv[0][1], 9999))):
            cent_label = gr[0].get("centrality_label", cent)
            out_png = png_dir / f"{safe_name(job)}_cent{safe_name(cent)}_{safe_name(cent_label)}.png" if write_png else None
            ok = plot_component_collection(
                rows=gr,
                components=SET_COMPONENTS,
                out_png=out_png,
                pdf=pdf,
                side=side,
                x_mode=x_mode,
                ylabel=f"systematic uncertainty [{side}, %]",
                title_suffix="Main systematic sets and final systematic",
                ymax=ymax,
                dpi=dpi,
            )
            if ok:
                n_pages += 1
    info(f"Wrote subtotal PDF: {pdf_path} pages={n_pages}")
    return n_pages


def plot_unfolding_breakdown(
    rows: List[Dict[str, str]],
    out_dir: Path,
    side: str,
    x_mode: str,
    ymax: float,
    dpi: int,
    write_png: bool,
) -> int:
    """Plot only the internal decomposition of the Unfolding set."""
    grouped = group_by_job_cent(rows)
    pdf_path = out_dir / f"unfolding_breakdown_{side}.pdf"
    png_dir = out_dir / "png" / f"unfolding_breakdown_{side}"
    n_pages = 0
    with PdfPages(pdf_path) as pdf:
        for (job, cent), gr in sorted(grouped.items(), key=lambda kv: (kv[0][0], parse_float(kv[0][1], 9999))):
            cent_label = gr[0].get("centrality_label", cent)
            out_png = png_dir / f"{safe_name(job)}_cent{safe_name(cent)}_{safe_name(cent_label)}.png" if write_png else None
            ok = plot_component_collection(
                rows=gr,
                components=UNFOLDING_BREAKDOWN_COMPONENTS,
                out_png=out_png,
                pdf=pdf,
                side=side,
                x_mode=x_mode,
                ylabel=f"systematic uncertainty [{side}, %]",
                title_suffix="Unfolding internal breakdown",
                ymax=ymax,
                dpi=dpi,
            )
            if ok:
                n_pages += 1
    info(f"Wrote unfolding breakdown PDF: {pdf_path} pages={n_pages}")
    return n_pages


# -----------------------------------------------------------------------------
# Barlow diagnostics
# -----------------------------------------------------------------------------

def latest_versioned_file(patterns: Sequence[str]) -> Optional[Path]:
    """Return newest file matching patterns with _v<number>, or newest by mtime."""
    candidates: List[Path] = []
    for pattern in patterns:
        candidates.extend(Path().glob(pattern))

    candidates = [p for p in candidates if p.exists() and p.is_file() and p.stat().st_size > 0]
    if not candidates:
        return None

    def key(p: Path) -> Tuple[int, float, str]:
        m = re.search(r"_v(\d+)\.", p.name)
        version = int(m.group(1)) if m else -1
        return (version, p.stat().st_mtime, str(p))

    return sorted(candidates, key=key)[-1]


def find_default_final_path() -> Optional[Path]:
    """Find the newest final table, preferring GUI-named outputs."""
    search_sets = [
        ["Systematics/final_systematics_gui_v*.tsv", "final_systematics_gui_v*.tsv"],
        ["Systematics/final_systematics_v*.tsv", "final_systematics_v*.tsv"],
    ]
    for patterns in search_sets:
        p = latest_versioned_file(patterns)
        if p is not None:
            return p
    return None


def infer_details_path(final_path: Path) -> Optional[Path]:
    """Infer details TSV from final TSV path and version.

    Supports both naming schemes:
      final_systematics_gui_v2.tsv -> final_systematics_gui_barlow_details_v2.tsv
      final_systematics_v9.tsv     -> final_systematics_barlow_details_v9.tsv
    """
    parent = final_path.parent
    name = final_path.name
    candidates: List[Path] = []

    m_gui = re.search(r"final_systematics_gui_v(\d+)\.tsv$", name)
    if m_gui:
        v = m_gui.group(1)
        candidates.append(parent / f"final_systematics_gui_barlow_details_v{v}.tsv")
        candidates.append(parent / f"final_systematics_barlow_details_v{v}.tsv")

    m_plain = re.search(r"final_systematics_v(\d+)\.tsv$", name)
    if m_plain:
        v = m_plain.group(1)
        candidates.append(parent / f"final_systematics_barlow_details_v{v}.tsv")
        candidates.append(parent / f"final_systematics_gui_barlow_details_v{v}.tsv")

    candidates.append(parent / "final_systematics_gui_barlow_details.tsv")
    candidates.append(parent / "final_systematics_barlow_details.tsv")
    candidates.extend(sorted(parent.glob("final_systematics_gui_barlow_details_v*.tsv")))
    candidates.extend(sorted(parent.glob("final_systematics_barlow_details_v*.tsv")))

    seen = set()
    unique: List[Path] = []
    for p in candidates:
        if p in seen:
            continue
        seen.add(p)
        unique.append(p)

    for p in unique:
        if p.exists() and p.stat().st_size > 0:
            return p
    return None


def _pretty_variation_name(text: str) -> str:
    """Human-readable variation label with case/spacing normalization."""
    raw = str(text or "").strip()
    if not raw:
        return ""
    key = raw.strip().lower()
    key = key.replace("-", "_").replace(" ", "_")
    key = re.sub(r"_+", "_", key)
    if key in PRETTY_VARIATION:
        return PRETTY_VARIATION[key]
    # A few forgiving pattern matches for names generated by different scripts.
    if "nhits" in key and "13" in key:
        return "nHitsFit 13"
    if "nhits" in key and "17" in key:
        return "nHitsFit 17"
    if "dca" in key and ("2_8" in key or "2.8" in key):
        return "DCA 2.8"
    if "dca" in key and ("3_2" in key or "3.2" in key):
        return "DCA 3.2"
    if "hadronic" in key or "hadr" in key:
        return "Hadr. correction"
    if "kt" in key and "drop" in key:
        return "kT dropped"
    if "rec" in key and "eff" in key:
        return "Tracking eff."
    return raw


def variation_label(row: Dict[str, str]) -> str:
    """Label one JetsReco variation.

    Prefer a specific human-readable label.  Some GUI details rows for paired
    variations only contain a generic source name such as ``jet_dca`` or
    ``jet_nhitsfit``; in that case the numeric JetsReco variation code is used
    to recover DCA2.8/DCA3.2 or nHitsFit13/nHitsFit17.
    """
    code = str(row.get("variation_code", "")).strip()
    code_label = JETSRECO_CODE_LABELS.get(code, "")

    candidates = [
        str(row.get("variation_name", "")).strip(),
        str(row.get("source_name", "")).strip(),
        str(row.get("group_rule", "")).strip(),
        code,
    ]

    best = ""
    for val in candidates:
        if not val:
            continue
        pretty = _pretty_variation_name(val)
        if not best:
            best = pretty
        # If the text already identifies the pair member, keep it immediately.
        if pretty not in {"DCA", "nHitsFit", "variation", ""}:
            return pretty

    # If the best text was only the pair family, recover the member from the code.
    if best in {"DCA", "nHitsFit"} and code_label:
        return code_label
    if code_label:
        return code_label
    return best or "variation"

def parse_stat_diff_from_note(note: str) -> float:
    m = re.search(r"(?:^|;)stat_diff=([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)", str(note or ""))
    if not m:
        return math.nan
    return parse_float(m.group(1), math.nan)


def threshold_pct(row: Dict[str, str]) -> float:
    """Return Barlow threshold in percent.

    Priority:
      1) compute directly from nominal/variation statistical errors,
      2) parse stat_diff from barlow_note,
      3) use raw_abs_pct / barlow_value.

    This makes the threshold panel robust even if the details table was produced
    by an older GUI script where barlow_note was incomplete.
    """
    nominal = abs(parse_float(row.get("nominal_value"), math.nan))

    nom_stat = parse_float(row.get("nominal_stat_abs"), math.nan)
    var_stat = parse_float(row.get("variation_stat_abs"), math.nan)
    if finite(nominal) and nominal > 0 and finite(nom_stat) and finite(var_stat):
        note = str(row.get("barlow_note", "")).lower()
        if "independent" in note:
            stat_diff2 = nom_stat * nom_stat + var_stat * var_stat
        else:
            # Default used by the combiner/presentation-style Barlow.
            stat_diff2 = abs(var_stat * var_stat - nom_stat * nom_stat)
        return 100.0 * math.sqrt(max(0.0, stat_diff2)) / nominal

    stat_diff = parse_stat_diff_from_note(row.get("barlow_note", ""))
    if finite(stat_diff) and finite(nominal) and nominal > 0:
        return 100.0 * stat_diff / nominal

    raw = parse_float(row.get("raw_abs_pct"), math.nan)
    if not finite(raw):
        raw = abs(parse_float(row.get("raw_diff_pct"), math.nan))
    b = parse_float(row.get("barlow_value"), math.nan)
    if finite(raw) and finite(b) and b > 0:
        return raw / b
    return math.nan


def detail_raw_abs(row: Dict[str, str]) -> float:
    v = parse_float(row.get("raw_abs_pct"), math.nan)
    if finite(v):
        return abs(v)
    return abs(parse_float(row.get("raw_diff_pct"), math.nan))


def detail_used_abs(row: Dict[str, str]) -> float:
    v = parse_float(row.get("used_abs_pct"), math.nan)
    if finite(v):
        return abs(v)
    return abs(parse_float(row.get("used_diff_pct"), math.nan))


def _text_has_jetsreco_source_pattern(text: str) -> bool:
    """True for JetsReco source names, but not for PriorShape jet-pT tilts."""
    t = str(text or "").strip().lower()
    if not t:
        return False
    t = t.replace("-", "_").replace(" ", "_")
    t = re.sub(r"_+", "_", t)

    # Explicitly reject prior-shape jet-pT tilts that caused false positives.
    if re.search(r"jet_?pt.*(plus|minus|\+|-)\s*20", t):
        return False
    if re.search(r"\b(pt|jetpt|jet_pt)_(plus|minus)?20\b", t):
        return False

    patterns = [
        r"nhits\s*fit|nhitsfit|n_hits_fit",
        r"dca\s*(2[._]?8|3[._]?2)?\b|jet_dca",
        r"kt\s*drop|ktdrop|k_t_drop|ktdropped|dropped",
        r"hadronic|hadr",
        r"rec_?eff|tracking_?eff|track_?eff|efficiency",
    ]
    return any(re.search(pat, t) for pat in patterns)


def _is_prior_shape_like(row: Dict[str, str]) -> bool:
    txt = " ".join(str(row.get(k, "")) for k in [
        "group_rule", "source_group", "source_name", "variation_name", "variation_code"
    ]).lower()
    return ("prior" in txt or "priorshape" in txt or "prior_shape" in txt or
            "jet_pt_plus20" in txt or "jet_pt_minus20" in txt or
            "jetpt_plus20" in txt or "jetpt_minus20" in txt)


def filter_jetsreco_details(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    """Select JetsReco detail rows without pulling in prior-shape jet-pT tilts.

    Different GUI/config versions store paired sources differently.  Some use
    group names like ``jetsreco_jet_dca``; others may only expose ``jet_dca`` or
    ``jet_nhitsfit`` in source/variation names.  Therefore we accept genuine
    JetsReco patterns (nHitsFit, DCA, kT dropped, tracking efficiency,
    hadronic correction), while explicitly rejecting prior-shape jet-pT tilts
    such as ``jet_pt_plus20`` and ``jet_pt_minus20``.
    """
    out = []
    for r in rows:
        group_rule = str(r.get("group_rule", "")).strip().lower()
        source_group = str(r.get("source_group", "")).strip().lower()
        source_name = str(r.get("source_name", "")).strip().lower()
        variation_name = str(r.get("variation_name", "")).strip().lower()
        variation_code = str(r.get("variation_code", "")).strip()

        is_jetsreco_group_rule = (
            group_rule == "jetsreco"
            or group_rule.startswith("jetsreco_")
            or group_rule.startswith("jets_reco")
            or group_rule.startswith("jets_reconstruction")
            or _text_has_jetsreco_source_pattern(group_rule)
        )
        is_jetsreco_source_group = source_group in {
            "jetsreco",
            "jets reco",
            "jets reconstruction",
            "jetreco",
            "jet reco",
        }
        has_jetsreco_source_pattern = (
            _text_has_jetsreco_source_pattern(source_name)
            or _text_has_jetsreco_source_pattern(variation_name)
            or variation_code in JETSRECO_CODE_LABELS
        )

        # If it is clearly a prior-shape row, keep it out even when names start
        # with jet_.  However, do not reject rows explicitly selected by a
        # JetsReco group/source-group rule.
        if _is_prior_shape_like(r) and not (is_jetsreco_group_rule or is_jetsreco_source_group):
            continue

        if is_jetsreco_group_rule or is_jetsreco_source_group or has_jetsreco_source_pattern:
            out.append(r)
    return out

def group_details(rows: Sequence[Dict[str, str]]) -> Dict[Tuple[str, str], List[Dict[str, str]]]:
    groups: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for r in rows:
        groups[(str(r.get("job_label", "")), str(r.get("centrality", "")))].append(r)
    return dict(groups)


def plot_barlow_one(
    rows: Sequence[Dict[str, str]],
    out_png: Optional[Path],
    pdf: PdfPages,
    x_mode: str,
    ymax: float,
    dpi: int,
) -> bool:
    if not rows:
        return False

    # Use representative bin rows for x labels.
    by_bin: Dict[str, Dict[str, str]] = {}
    for r in rows:
        by_bin.setdefault(str(r.get("bin", "")), r)
    bin_rows = sorted(by_bin.values(), key=row_sort_key)
    x_ref, xlabels, xlabel = x_values_and_labels(bin_rows, x_mode)
    x_by_bin = {str(r.get("bin", "")): x for r, x in zip(bin_rows, x_ref)}

    # Group values by variation.
    var_rows: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for r in rows:
        var_rows[variation_label(r)].append(r)

    if not var_rows:
        return False

    fig, axes = plt.subplots(3, 1, figsize=(11.0, 10.0), sharex=True)
    titles = ["Before Barlow: |variation - nominal|", "Barlow threshold", "After Barlow"]
    ylabels = ["raw [%]", r"$\sigma_{\mathrm{stat,diff}}$ [%]", "used [%]"]

    for ax, title, ylabel in zip(axes, titles, ylabels):
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.28)

    for i, (label, sub) in enumerate(sorted(var_rows.items(), key=lambda kv: kv[0])):
        sub = sorted(sub, key=row_sort_key)
        xs: List[float] = []
        raw: List[float] = []
        thr: List[float] = []
        used: List[float] = []
        for r in sub:
            bkey = str(r.get("bin", ""))
            if bkey not in x_by_bin:
                continue
            xs.append(x_by_bin[bkey])
            raw.append(detail_raw_abs(r))
            thr.append(threshold_pct(r))
            used.append(detail_used_abs(r))

        if not xs:
            continue
        marker = MARKERS[i % len(MARKERS)]
        color = color_for_jetsreco_variation(label)
        axes[0].plot(xs, raw, marker=marker, linewidth=1.45, markersize=4.0, label=label, color=color)
        axes[1].plot(xs, thr, marker=marker, linewidth=1.45, markersize=4.0, label=label, color=color)
        axes[2].plot(xs, used, marker=marker, linewidth=1.45, markersize=4.0, label=label, color=color)

    axes[-1].set_xlabel(xlabel)
    autoset_xticks(axes[-1], x_ref, xlabels)
    axes[0].set_title(make_title(rows, "JetsReco Barlow diagnostics") + "\n" + titles[0], fontsize=10)

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(
            handles,
            labels,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0,
            fontsize=7.5,
            frameon=False,
        )
    fig.tight_layout(rect=[0.0, 0.0, 0.78, 1.0])

    if out_png is not None:
        out_png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_png, dpi=dpi)
    pdf.savefig(fig)
    plt.close(fig)
    return True


def plot_barlow_diagnostics(
    details_rows: List[Dict[str, str]],
    out_dir: Path,
    x_mode: str,
    ymax: float,
    dpi: int,
    write_png: bool,
) -> int:
    jr = filter_jetsreco_details(details_rows)
    if not jr:
        warn("No JetsReco rows found in Barlow details table; Barlow diagnostics skipped.")
        return 0

    grouped = group_details(jr)
    pdf_path = out_dir / "jetsreco_barlow_diagnostics.pdf"
    png_dir = out_dir / "png" / "jetsreco_barlow"
    n_pages = 0
    with PdfPages(pdf_path) as pdf:
        for (job, cent), gr in sorted(grouped.items(), key=lambda kv: (kv[0][0], parse_float(kv[0][1], 9999))):
            cent_label = gr[0].get("centrality_label", cent)
            out_png = png_dir / f"{safe_name(job)}_cent{safe_name(cent)}_{safe_name(cent_label)}.png" if write_png else None
            ok = plot_barlow_one(gr, out_png, pdf, x_mode=x_mode, ymax=ymax, dpi=dpi)
            if ok:
                n_pages += 1
    info(f"Wrote JetsReco Barlow PDF: {pdf_path} pages={n_pages}")
    return n_pages



def print_barlow_summary(details_rows: List[Dict[str, str]]) -> None:
    jr = filter_jetsreco_details(details_rows)
    info(f"JetsReco detail rows selected: {len(jr)} / {len(details_rows)}")
    counts: Dict[Tuple[str, str, str, str], int] = defaultdict(int)
    for r in jr:
        key = (
            str(r.get("group_rule", "")),
            str(r.get("source_group", "")),
            str(r.get("source_name", "")),
            str(r.get("variation_name", "")),
        )
        counts[key] += 1
    for (group_rule, source_group, source_name, variation_name), n in sorted(counts.items(), key=lambda kv: (kv[0][0], kv[0][2], kv[0][3])):
        info(f"  n={n:5d}  group_rule={group_rule!r}  source_group={source_group!r}  source_name={source_name!r}  variation_name={variation_name!r}")

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot final-combiner diagnostic systematics plots.")
    p.add_argument("--final", default="", help="Final systematics TSV. If empty, newest final_systematics_gui_v*.tsv is used first, then final_systematics_v*.tsv.")
    p.add_argument("--details", default="", help="Optional Barlow details TSV. If empty, inferred from --final.")
    p.add_argument("--out-dir", default="", help="Output directory. If empty, use <final-table-dir>/Plots/SystematicsDiagnostics.")
    p.add_argument("--side", choices=["sym", "down", "up"], default="sym", help="Which uncertainty side to plot for subtotals/breakdown.")
    p.add_argument("--x-mode", choices=["center", "bin"], default="center", help="Use bin centers or bin index on x-axis.")
    p.add_argument("--ymax", type=float, default=0.0, help="Manual y-axis max for subtotal/breakdown plots. 0 = automatic.")
    p.add_argument("--barlow-ymax", type=float, default=0.0, help="Manual y-axis max for Barlow diagnostics. 0 = automatic.")
    p.add_argument("--dpi", type=int, default=180, help="PNG DPI.")
    p.add_argument("--no-png", action="store_true", help="Only write PDFs.")
    p.add_argument("--skip-subtotals", action="store_true", help="Do not plot intermediate sums.")
    p.add_argument("--skip-unfolding-breakdown", action="store_true", help="Do not plot the internal Unfolding-set breakdown.")
    p.add_argument("--print-columns", action="store_true", help="Print resolved *_<side>_pct columns and exit.")
    p.add_argument("--skip-barlow", action="store_true", help="Do not plot JetsReco Barlow diagnostics.")
    p.add_argument("--print-barlow-summary", action="store_true", help="Print selected JetsReco Barlow detail groups and exit.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.final:
        final_path = Path(args.final)
    else:
        found = find_default_final_path()
        if found is None:
            print("[error] Could not find final_systematics_gui_v*.tsv or final_systematics_v*.tsv.", file=sys.stderr)
            return 1
        final_path = found
        info(f"Auto-selected final table: {final_path}")

    if args.out_dir:
        out_dir = Path(args.out_dir)
    else:
        out_dir = final_path.parent / "Plots" / "SystematicsDiagnostics"

    write_png = not bool(args.no_png)

    if not final_path.exists():
        print(f"[error] Missing final systematics table: {final_path}", file=sys.stderr)
        return 1

    header, final_rows = read_tsv(final_path)
    if not final_rows:
        print(f"[error] No rows in final systematics table: {final_path}", file=sys.stderr)
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)
    info(f"Read final rows: {len(final_rows)} from {final_path}")
    if args.print_columns:
        print_component_column_diagnostics(header, args.side)
        return 0

    if not args.skip_subtotals:
        plot_subtotals(
            rows=final_rows,
            out_dir=out_dir,
            side=args.side,
            x_mode=args.x_mode,
            ymax=args.ymax,
            dpi=args.dpi,
            write_png=write_png,
        )

    if not args.skip_unfolding_breakdown:
        plot_unfolding_breakdown(
            rows=final_rows,
            out_dir=out_dir,
            side=args.side,
            x_mode=args.x_mode,
            ymax=args.ymax,
            dpi=args.dpi,
            write_png=write_png,
        )

    if not args.skip_barlow:
        details_path: Optional[Path]
        if args.details:
            details_path = Path(args.details)
        else:
            details_path = infer_details_path(final_path)

        if details_path is None or not details_path.exists():
            warn("Barlow details TSV not found. Run the combiner with write_details=true or pass --details explicitly.")
        else:
            _, details_rows = read_tsv(details_path)
            info(f"Read detail rows: {len(details_rows)} from {details_path}")
            if args.print_barlow_summary:
                print_barlow_summary(details_rows)
                return 0
            plot_barlow_diagnostics(
                details_rows=details_rows,
                out_dir=out_dir,
                x_mode=args.x_mode,
                ymax=args.barlow_ymax,
                dpi=args.dpi,
                write_png=write_png,
            )

    info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
