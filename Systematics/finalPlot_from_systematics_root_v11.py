#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
finalPlot_from_systematics_root.py

Plot final preliminary-style figures from the ROOT file produced by
systematics_final_combiner_v7.py.

Expected ROOT structure, for example:
  l11/g_value_stat_0_10pct
  l11/g_value_syst_0_10pct
  rcp_l11_5_20/g_value_stat_0_10pct_40_80pct
  rcp_l11_5_20/g_value_syst_0_10pct_40_80pct

The plotting style follows the older finalPlot.py script, but the uncertainties
are now read from the final systematic ROOT file instead of recomputing the
iteration envelope inside this script.

Optional LIDO/model comparison:
  python3 finalPlot_from_systematics_root_v9.py --show-lido

This reads Systematics/lido_model_rcp.root by default and draws it only in
results_RCP_v9.pdf.
"""

from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

try:
    import uproot  # type: ignore
except Exception:
    uproot = None

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerBase
from matplotlib.patches import Polygon

# ---------------- default user inputs ----------------
DEFAULT_ROOT_FILE = "Systematics/final_systematics_results_v10.root"
DEFAULT_OUTPUT_DIR = "FinalPlot"
DEFAULT_ALL_GRAPHS_SUBDIR = "AllGraphs"
VERSION_TAG = "_v9"

# Optional LIDO/model comparison. It is drawn only in results_RCP*_v9.pdf,
# i.e. in plot_rcp_pages(), not in the overlay or spectra plots.
SHOW_LIDO_MODEL_DEFAULT = False
DEFAULT_LIDO_ROOT_FILE = "Systematics/lido_model_rcp.root"
LIDO_MODEL_LABEL = "LIDO model"

# These are only used in plot annotations now. The input ROOT file already
# contains the final selected iteration and final systematic uncertainties.
Method = "ICS"
jet_pt_range_label = r"$5 < p_{\mathrm{T,jet}}$ [GeV/$c$] $< 20$, $|\eta_{\mathrm{jet}}| < 0.6$"
d0_pt_range_label = r"$1 < p_{\mathrm{T,D}^{0}}$ [GeV/$c$] $< 10$, $|y_{\mathrm{D}^{0}}| < 0.6$"

# Keep the same choices as in the old plotting script.
lambdas_to_plot = [5, 1, 2, 3, 4]
#lambdas_overlay = [5, 1, 2, 3, 4]
lambdas_overlay = [5, 1, 4]
rcp_lambdas_to_plot = [5, 1, 2, 3, 4]
rcp_overlay_lambdas = [5, 1, 2, 3, 4]

# --- LIDO drawing control ---
CLIP_LIDO_TO_DATA_RANGE = True

# Optional manual x-ranges for LIDO.
# Keys are iLamb:
#   1 = lambda_1^1
#   2 = lambda_1.5^1
#   3 = lambda_2^1
#   4 = lambda_3^1
#   5 = lambda_0.5^1
#
# If value is None, the range is taken from measured data.
LIDO_MANUAL_X_RANGE = {
    # 5: (0.0, 1.0),
    # 1: (0.0, 0.5),
    # 2: (0.0, 0.4),
    # 3: (0.0, 0.4),
    # 4: (0.0, 0.2),
}

# Branching ratio uncertainty.
# The spectra are divided by BrRatioD0, so the relative uncertainty is
# BR_RATIO_D0_ERR / BR_RATIO_D0.  By default this is added in quadrature
# to the final spectrum systematic uncertainty.  It is not applied to R_CP,
# where the common branching-ratio normalization cancels between numerator
# and denominator.
APPLY_BR_UNCERTAINTY_TO_SPECTRA = True
APPLY_BR_UNCERTAINTY_TO_RCP = False
BR_RATIO_D0 = 0.0395
BR_RATIO_D0_ERR = 0.0003

# Horizontal legend settings for overlay plots with all angularities.
# The legend is kept inside the axes; only its entries are arranged horizontally.
OVERLAY_LEGEND_NCOL = 3
OVERLAY_LEGEND_FONTSIZE = 8
OVERLAY_LEGEND_LOC = "upper right"
OVERLAY_LEGEND_BBOX = (0.985, 0.8)

# Separate placement for the R_CP angularity-overlay plot.
# This affects only results_RCP_overlay*.pdf.
RCP_OVERLAY_LEGEND_NCOL = 3
RCP_OVERLAY_LEGEND_FONTSIZE = 8
RCP_OVERLAY_LEGEND_LOC = "lower right"
RCP_OVERLAY_LEGEND_BBOX = (0.985, 0.02)
RCP_OVERLAY_HEADER_X = 0.97
RCP_OVERLAY_HEADER_Y = 0.18
RCP_OVERLAY_HEADER_FONTSIZE = 8.5

# Manual bin-visibility switch for each plotted variable.
# Set first/last to False if the corresponding edge bin should be hidden
# in spectra and R_CP plots.  The values are applied after reading the ROOT
# graphs, so the input ROOT file is not modified.
SHOW_EDGE_BINS = {
    0: {"first": True, "last": True},  # z
    1: {"first": True, "last": False},  # lambda_1^1
    2: {"first": True, "last": False},  # lambda_1.5^1
    3: {"first": True, "last": False},  # lambda_2^1
    4: {"first": True, "last": False},  # lambda_3^1
    5: {"first": True, "last": False},  # lambda_0.5^1
    6: {"first": True, "last": True},  # p_T^D
}

# centrality labels for iCent = 0,1,2
cent_labels = ["0–10%", "10–40%", "40–80%"]
cent_tags = ["0_10pct", "10_40pct", "40_80pct"]

# R_CP labels/tags for centrality-like index in the final ROOT file
rcp_types = ["0/1", "0/2", "1/2"]
rcp_titles = {
    "0/1": "0–10% / 10–40%",
    "0/2": "0–10% / 40–80%",
    "1/2": "10–40% / 40–80%",
}
rcp_tags = {
    "0/1": "0_10pct_10_40pct",
    "0/2": "0_10pct_40_80pct",
    "1/2": "10_40pct_40_80pct",
}

# N_coll values used for R_CP.
# The uncertainty is drawn separately as a global normalization box.
# It is not included in the bin-by-bin systematic boxes from the ROOT file.
ncoll_values = np.array([951.99764, 397.03073, 58.07356], dtype=float)
ncoll_errors = np.array([28.20515, 31.64791, 15.19217], dtype=float)

rcp_cent_indices = {
    "0/1": (0, 1),
    "0/2": (0, 2),
    "1/2": (1, 2),
}

lambda_labels = {
    0: r"$z$",
    1: r"$\lambda_{1}^{1}$",
    2: r"$\lambda_{1.5}^{1}$",
    3: r"$\lambda_{2}^{1}$",
    4: r"$\lambda_{3}^{1}$",
    5: r"$\lambda_{0.5}^{1}$",
    6: r"$p_{\mathrm{T}}^{\mathrm{D}}$",
}

lambda_titles = {
    0: "z",
    1: r"$\lambda_{1}^{1}$",
    2: r"$\lambda_{1.5}^{1}$",
    3: r"$\lambda_{2}^{1}$",
    4: r"$\lambda_{3}^{1}$",
    5: r"$\lambda_{0.5}^{1}$",
    6: r"$p_{\mathrm{T}}^{\mathrm{D}}$",
}

lambda_colors = {
    1: "tab:blue",
    2: "tab:orange",
    3: "tab:green",
    4: "tab:red",
    5: "tab:purple",
    6: "tab:brown",
    0: "tab:gray",
}

# Mapping from the old lambda index to the job directory in the final ROOT file.
# The final ROOT file uses safe_name(job_label), so pTD becomes ptd.
spectrum_job_dir = {
    0: "z",
    1: "l11",
    2: "l15",
    3: "l21",
    4: "l31",
    5: "l051",
    6: "ptd",
}

rcp_job_dir = {
    0: "rcp_z_5_20",
    1: "rcp_l11_5_20",
    2: "rcp_l15_5_20",
    3: "rcp_l21_5_20",
    4: "rcp_l31_5_20",
    5: "rcp_l051_5_20",
    6: "rcp_ptd_5_20",
}

# Names used by plot_lido_rcp_v3.py in Systematics/lido_model_rcp.root.
# Example path:
#   RCP_0_10_over_40_80/gRCP_0_10_over_40_80_lambda_1d0
lido_rcp_dirs = {
    "0/1": "RCP_0_10_over_10_40",
    "0/2": "RCP_0_10_over_40_80",
    "1/2": "RCP_10_40_over_40_80",
}

lido_observable_keys = {
    0: "z",
    1: "lambda_1d0",
    2: "lambda_1d5",
    3: "lambda_2d0",
    4: "lambda_3d0",
    5: "lambda_0d5",
    6: "pT",
}

ANGULARITY_OVERLAY_XLABEL = r"$\lambda_{\alpha}^{1}$"

SHOW_PRELIMINARY_LABEL = True
PRELIMINARY_LABEL = "STAR Preliminary"
PRELIMINARY_LABEL_COLOR = "red"
PRELIMINARY_LABEL_FONTSIZE = 15
PRELIMINARY_LABEL_X = 0.02
PRELIMINARY_LABEL_Y = 0.97

cent_colors = {
    0: "tab:blue",
    1: "tab:orange",
    2: "tab:green",
}

cent_markers = {
    0: "o",  # circle
    1: "s",  # square
    2: "D",  # triangle
}

rcp_colors = {
    "0/1": "#6A3D9A",  # purple
    "0/2": "tab:red",  # red
    "1/2": "#B15928",  # brown
}

rcp_markers = {
    "0/1": "o",
    "0/2": "s",
    "1/2": "D",
}

rcp_axis_labels = {
    "0/1": r"$R_{\mathrm{CP}}$ (0–10% / 10–40%)",
    "0/2": r"$R_{\mathrm{CP}}$ (0–10% / 40–80%)",
    "1/2": r"$R_{\mathrm{CP}}$ (10–40% / 40–80%)",
}

MAIN_LEGEND_HEADER_X = 0.45
MAIN_LEGEND_CENT_X   = 0.438
MAIN_LEGEND_SYST_X   = 0.438   # posun Syst. unc. doleva

MAIN_LEGEND_HEADER_Y = 0.97
MAIN_LEGEND_CENT_Y   = 0.69
MAIN_LEGEND_SYST_Y   = 0.64
MAIN_LEGEND_FONTSIZE = 11.5

AXIS_LABEL_FONTSIZE = 14
TICK_LABEL_FONTSIZE = 14

RCP_LEGEND_HEADER_X = 0.04
RCP_LEGEND_HEADER_Y = 0.965

RCP_LEGEND_UNC_X = 0.027
RCP_LEGEND_UNC_Y = 0.74

RCP_LEGEND_FONTSIZE = 12

RCP_AXIS_LABEL_FONTSIZE = 16
RCP_TICK_LABEL_FONTSIZE = 16
RCP_YLABEL_PAD = 4

RCP_OVERLAY_HEADER_X = 0.56
RCP_OVERLAY_HEADER_Y = 0.345
RCP_OVERLAY_HEADER_FONTSIZE = 12

RCP_OVERLAY_ANG_X = 0.55
RCP_OVERLAY_ANG_Y = 0.125
RCP_OVERLAY_ANG_FONTSIZE = 12

RCP_OVERLAY_UNC_X = 0.55
RCP_OVERLAY_UNC_Y = 0.07
RCP_OVERLAY_UNC_FONTSIZE = 12

# Separate legend placement for results_by_centrality*.pdf.
# It is intentionally independent of the R_CP overlay settings.
SPECTRA_OVERLAY_HEADER_X = 0.59
SPECTRA_OVERLAY_HEADER_Y = 0.97
SPECTRA_OVERLAY_HEADER_FONTSIZE = 12

SPECTRA_OVERLAY_ANG_X = 0.58
SPECTRA_OVERLAY_ANG_Y = 0.76
SPECTRA_OVERLAY_ANG_FONTSIZE = 12

SPECTRA_OVERLAY_UNC_X = 0.58
SPECTRA_OVERLAY_UNC_Y = 0.71
SPECTRA_OVERLAY_UNC_FONTSIZE = 12

ANGULARITY_OVERLAY_LABEL_FONTSIZE = 13
ANGULARITY_OVERLAY_TICK_LABEL_FONTSIZE = 12

SPECTRA_OVERLAY_CENT_X = 0.96
SPECTRA_OVERLAY_CENT_Y = 0.09
SPECTRA_OVERLAY_CENT_FONTSIZE = 19
# -----------------------------------------------------
def add_spectra_overlay_centrality_legend(ax: Any, cent_label: str) -> None:
    """Add centrality label aligned by its right edge."""
    ax.text(
        SPECTRA_OVERLAY_CENT_X,
        SPECTRA_OVERLAY_CENT_Y,
        cent_label,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=SPECTRA_OVERLAY_CENT_FONTSIZE
#        fontweight="bold",
    )

def add_rcp_overlay_legend(ax: Any, handles: list, labels: list) -> None:
    # Header
    ax.text(
        RCP_OVERLAY_HEADER_X,
        RCP_OVERLAY_HEADER_Y,
        "\n".join(legend_header_labels()),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=RCP_OVERLAY_HEADER_FONTSIZE,
        linespacing=1.30,
    )

    # Angularity legend
    leg_ang = ax.legend(
        handles,
        labels,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(RCP_OVERLAY_ANG_X, RCP_OVERLAY_ANG_Y),
        bbox_transform=ax.transAxes,
        fontsize=RCP_OVERLAY_ANG_FONTSIZE,
        ncol=5,
        columnspacing=0.7,
        handlelength=1.1,
        handletextpad=0.30,
        borderaxespad=0.0,
    )
    ax.add_artist(leg_ang)

    # Uncertainty legend
    sys_patch = MultiColorSysBox(
        [lambda_colors[i] for i in rcp_overlay_lambdas],
        alpha=0.55,
        skew=0.18,
        yshift=-0.15,
    )

    unc_handles = [
        sys_patch,
        ncoll_legend_patch(),
    ]

    unc_labels = [
        r"$\mathrm{Syst.\ unc.}$",
        r"$\mathrm{Rel.\ }N_{\mathrm{coll}}\mathrm{\ unc.}$",
    ]

    ax.legend(
        unc_handles,
        unc_labels,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(RCP_OVERLAY_UNC_X, RCP_OVERLAY_UNC_Y),
        bbox_transform=ax.transAxes,
        fontsize=RCP_OVERLAY_UNC_FONTSIZE,
        ncol=2,
        columnspacing=0.9,
        handlelength=2.3,
        handleheight=0.9,
        handletextpad=0.35,
        borderaxespad=0.0,
        handler_map={MultiColorSysBox: HandlerMultiColorSysBox()},
    )
    
def lido_legend_handle() -> LidoBandLegend:
    return LidoBandLegend(
        facecolor="tab:blue",
        alpha=0.25,
        linecolor="tab:blue",
        linewidth=1.5,
    )


def add_rcp_page_legend(ax: Any, color: str, marker: str, show_lido_model: bool = False) -> None:
    ax.text(
        RCP_LEGEND_HEADER_X,
        RCP_LEGEND_HEADER_Y,
        "\n".join(legend_header_labels()),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=RCP_LEGEND_FONTSIZE,
        linespacing=1.30,
    )

    unc_handles = [
        Patch(facecolor=color, edgecolor="none", alpha=0.18),
        ncoll_legend_patch(),
    ]

    unc_labels = [
        r"$\mathrm{Syst.\ unc.}_{\ }$",
        r"$\mathrm{Rel.\ }N_{\mathrm{coll}}\mathrm{\ unc.}$",
    ]

    ncol = 2
    if show_lido_model:
        unc_handles = [
            unc_handles[0],        # Syst. unc.
            lido_legend_handle(),  # LIDO -> visually second row
            unc_handles[1],        # N_coll -> visually first row, second column
        ]

        unc_labels = [
            unc_labels[0],
            LIDO_MODEL_LABEL,
            unc_labels[1],
        ]

    ax.legend(
        unc_handles,
        unc_labels,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(RCP_LEGEND_UNC_X, RCP_LEGEND_UNC_Y),
        bbox_transform=ax.transAxes,
        fontsize=RCP_LEGEND_FONTSIZE,
        ncol=ncol,
        columnspacing=0.8,
        handlelength=1.6,
        handletextpad=0.35,
        borderaxespad=0.0,
        handler_map={
            LidoBandLegend: HandlerLidoBandLegend(),
        },
    )

def add_rcp_ratio_label(ax: Any, rcp_type: str, color: str) -> None:
    ax.text(
        0.06,
        0.92,
        rcp_titles[rcp_type],
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11,
        color=color,
        fontweight="bold",
    )

def add_spectrum_main_legend(ax: Any, cent_handles: list, cent_labels: list, sys_handle: Any) -> None:
    ax.text(
        MAIN_LEGEND_HEADER_X,
        MAIN_LEGEND_HEADER_Y,
        "\n".join(legend_header_labels()),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=MAIN_LEGEND_FONTSIZE,
        linespacing=1.35,
    )

    leg_cent = ax.legend(
        cent_handles,
        cent_labels,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(MAIN_LEGEND_CENT_X, MAIN_LEGEND_CENT_Y),
        bbox_transform=ax.transAxes,
        fontsize=MAIN_LEGEND_FONTSIZE,
        ncol=3,
        columnspacing=0.5,
        handlelength=1.1,
        handletextpad=0.3,
        borderaxespad=0.0,
    )
    ax.add_artist(leg_cent)

    ax.legend(
        [sys_handle],
        [r"Syst. unc."],
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(MAIN_LEGEND_SYST_X, MAIN_LEGEND_SYST_Y),
        bbox_transform=ax.transAxes,
        fontsize=MAIN_LEGEND_FONTSIZE,
        handlelength=2.8,
        handleheight=0.9,
        handletextpad=0.3,
        borderaxespad=0.0,
        handler_map={MultiColorSysBox: HandlerMultiColorSysBox()},
    )


def add_spectra_overlay_legend(
    ax: Any,
    handles: list,
    labels: list,
    plotted_lambdas: list,
) -> None:
    """Split the results_by_centrality legend into angularities and systematics.

    The systematic-uncertainty legend entry is a multi-color box built directly
    from the currently plotted lambda indices.  Therefore it automatically
    follows the user choice in lambdas_overlay.
    """
    ax.text(
        SPECTRA_OVERLAY_HEADER_X,
        SPECTRA_OVERLAY_HEADER_Y,
        "\n".join(legend_header_labels()),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=SPECTRA_OVERLAY_HEADER_FONTSIZE,
        linespacing=1.25,
    )

    n_ang = max(1, len(handles))
    # Keep it compact for any reasonable choice of lambdas_overlay.
    ncol_ang = min(n_ang, 5)

    leg_ang = ax.legend(
        handles,
        labels,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(SPECTRA_OVERLAY_ANG_X, SPECTRA_OVERLAY_ANG_Y),
        bbox_transform=ax.transAxes,
        fontsize=SPECTRA_OVERLAY_ANG_FONTSIZE,
        ncol=ncol_ang,
        columnspacing=0.7,
        handlelength=1.1,
        handletextpad=0.30,
        borderaxespad=0.0,
    )
    ax.add_artist(leg_ang)

    sys_colors = [lambda_colors[i] for i in plotted_lambdas if i in lambda_colors]
    if not sys_colors:
        sys_colors = ["0.5"]

    sys_patch = MultiColorSysBox(
        sys_colors,
        alpha=0.55,
        skew=0.18,
        yshift=-0.15,
    )

    ax.legend(
        [sys_patch],
        [r"$\mathrm{Syst.\ unc.}$"],
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(SPECTRA_OVERLAY_UNC_X, SPECTRA_OVERLAY_UNC_Y),
        bbox_transform=ax.transAxes,
        fontsize=SPECTRA_OVERLAY_UNC_FONTSIZE,
        ncol=1,
        columnspacing=0.9,
        handlelength=2.3,
        handleheight=0.9,
        handletextpad=0.35,
        borderaxespad=0.0,
        handler_map={MultiColorSysBox: HandlerMultiColorSysBox()},
    )
    
class MultiColorSysBox:
    def __init__(self, colors, alpha=0.55, skew=0.18, yshift=-0.05):
        self.colors = colors
        self.alpha = alpha
        self.skew = skew
        self.yshift = yshift


class HandlerMultiColorSysBox(HandlerBase):
    def create_artists(
        self, legend, orig_handle,
        xdescent, ydescent, width, height, fontsize, trans
    ):
        colors = orig_handle.colors
        n = len(colors)

        x0 = xdescent
        w = width

        box_h = 0.70 * height
        y0 = ydescent + 0.5 * (height - box_h) + orig_handle.yshift * height
        h = box_h

        # total horizontal tilt of the internal boundaries
        skew = orig_handle.skew * w

        # boundary positions at bottom and top
        bottom = [x0]
        top = [x0]

        for j in range(1, n):
            xc = x0 + w * j / n
            bottom.append(xc - 0.5 * skew)
            top.append(xc + 0.5 * skew)

        bottom.append(x0 + w)
        top.append(x0 + w)

        artists = []
        for i, col in enumerate(colors):
            pts = [
                (bottom[i], y0),
                (bottom[i + 1], y0),
                (top[i + 1], y0 + h),
                (top[i], y0 + h),
            ]

            p = Polygon(
                pts,
                closed=True,
                facecolor=col,
                edgecolor="none",
                alpha=orig_handle.alpha,
                transform=trans,
            )
            artists.append(p)

        return artists

class LidoBandLegend:
    def __init__(self, facecolor="tab:blue", alpha=0.25, linecolor="tab:blue", linewidth=1.5):
        self.facecolor = facecolor
        self.alpha = alpha
        self.linecolor = linecolor
        self.linewidth = linewidth        

class HandlerLidoBandLegend(HandlerBase):
    def create_artists(
        self, legend, orig_handle,
        xdescent, ydescent, width, height, fontsize, trans
    ):
        # band
        rect = Rectangle(
            (xdescent, ydescent + 0.15 * height),
            width,
            0.70 * height,
            facecolor=orig_handle.facecolor,
            edgecolor="none",
            alpha=orig_handle.alpha,
            transform=trans,
        )

        # central line
        ymid = ydescent + 0.5 * height
        line = Line2D(
            [xdescent, xdescent + width],
            [ymid, ymid],
            color=orig_handle.linecolor,
            linewidth=orig_handle.linewidth,
            transform=trans,
        )

        return [rect, line]

def _to_str(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, bytes):
        return x.decode("utf-8", errors="ignore")
    return str(x)


def make_ylabel_from_xlabel(xlabel: str) -> str:
    """Return the final spectra y-axis label.

    This follows the normalization used in NormalizeFinalSpectra():
      * divided by the jet-acceptance width Delta eta_jet = 1.2,
      * divided by the plotted-variable bin width,
      * divided by 2 for D0 + anti-D0,
      * divided by the D0 -> K pi branching ratio.

    The old preliminary-style factor 1/(2*pi*x) is deliberately not used.
    """
    xclean = str(xlabel).strip()
    if xclean.startswith("$") and xclean.endswith("$"):
        xclean = xclean[1:-1]
    return (
        rf"$"
        rf"\frac{{1}}{{N_{{\mathrm{{evt}}}}}}"
        rf"\frac{{1}}{{\mathrm{{BR}}}}"
        rf"\frac{{\mathrm{{d}}^2N}}{{\mathrm{{d}}{xclean}\,\mathrm{{d}}\eta_{{\mathrm{{jet}}}}}}$"
    )


def safe_path_exists(f: Any, path: str) -> bool:
    try:
        f[path]
        return True
    except Exception:
        try:
            f[path + ";1"]
            return True
        except Exception:
            return False


def get_root_object(f: Any, path: str) -> Optional[Any]:
    try:
        return f[path]
    except Exception:
        try:
            return f[path + ";1"]
        except Exception:
            return None


def _member_array(obj: Any, name: str, n: Optional[int] = None, default: float = 0.0) -> np.ndarray:
    try:
        arr = np.asarray(obj.member(name), dtype=float)
        if n is not None:
            if len(arr) < n:
                arr = np.pad(arr, (0, n - len(arr)), constant_values=default)
            elif len(arr) > n:
                arr = arr[:n]
        return arr
    except Exception:
        if n is None:
            n = 0
        return np.full(n, default, dtype=float)


def read_tgraph_asymm(obj: Any) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return x, y, exlow, exhigh, eylow, eyhigh from TGraphAsymmErrors."""
    x = _member_array(obj, "fX")
    y = _member_array(obj, "fY")
    n = min(len(x), len(y))
    x = x[:n]
    y = y[:n]

    exl = _member_array(obj, "fEXlow", n=n, default=0.0)
    exh = _member_array(obj, "fEXhigh", n=n, default=0.0)
    eyl = _member_array(obj, "fEYlow", n=n, default=0.0)
    eyh = _member_array(obj, "fEYhigh", n=n, default=0.0)
    return x, y, exl, exh, eyl, eyh


def read_hist_values(obj: Any) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return x, y, symmetric xerr, yerr from TH1-like object."""
    y, edges = obj.to_numpy(flow=False)
    x = 0.5 * (edges[1:] + edges[:-1])
    xerr = 0.5 * (edges[1:] - edges[:-1])
    var = obj.variances(flow=False)
    if var is None:
        yerr = np.zeros_like(y, dtype=float)
    else:
        yerr = np.sqrt(np.clip(np.asarray(var, dtype=float), 0.0, None))
    return np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(xerr, dtype=float), np.asarray(yerr, dtype=float)


def get_final_result(f: Any, job_dir: str, tag: str, xlabel: str) -> Optional[Dict[str, np.ndarray | str]]:
    """Read central values, stat. errors and final systematics from final ROOT file."""
    stat_obj = get_root_object(f, f"{job_dir}/g_value_stat_{tag}")
    syst_obj = get_root_object(f, f"{job_dir}/g_value_syst_{tag}")

    if stat_obj is not None:
        x, y, exl, exh, stat_low, stat_high = read_tgraph_asymm(stat_obj)
    else:
        h_stat = get_root_object(f, f"{job_dir}/h_value_stat_{tag}")
        if h_stat is None:
            return None
        x, y, xerr, yerr = read_hist_values(h_stat)
        exl = xerr
        exh = xerr
        stat_low = yerr
        stat_high = yerr

    if syst_obj is not None:
        xs, ys, exls, exhs, syst_low, syst_high = read_tgraph_asymm(syst_obj)
        # Prefer x-binning from stat graph, but use syst errors if the stat graph had none.
        n = min(len(x), len(xs), len(y), len(ys))
        x, y, exl, exh = x[:n], y[:n], exl[:n], exh[:n]
        stat_low, stat_high = stat_low[:n], stat_high[:n]
        syst_low, syst_high = syst_low[:n], syst_high[:n]
        if np.all(exl == 0) and np.any(exls[:n] > 0):
            exl = exls[:n]
        if np.all(exh == 0) and np.any(exhs[:n] > 0):
            exh = exhs[:n]
    else:
        h_down = get_root_object(f, f"{job_dir}/h_syst_down_abs_{tag}")
        h_up = get_root_object(f, f"{job_dir}/h_syst_up_abs_{tag}")
        if h_down is None or h_up is None:
            syst_low = np.zeros_like(y, dtype=float)
            syst_high = np.zeros_like(y, dtype=float)
        else:
            _, syst_low, _, _ = read_hist_values(h_down)
            _, syst_high, _, _ = read_hist_values(h_up)
            n = min(len(x), len(y), len(syst_low), len(syst_high))
            x, y, exl, exh = x[:n], y[:n], exl[:n], exh[:n]
            stat_low, stat_high = stat_low[:n], stat_high[:n]
            syst_low, syst_high = syst_low[:n], syst_high[:n]

    return {
        "x": x,
        "y": y,
        "exl": exl,
        "exh": exh,
        "stat_low": stat_low,
        "stat_high": stat_high,
        "syst_low": syst_low,
        "syst_high": syst_high,
        "xlabel": xlabel,
        "ylabel": make_ylabel_from_xlabel(xlabel),
    }


def get_lido_model_result(
    lido_f: Optional[Any],
    rcp_type: str,
    iLamb: int,
) -> Optional[Dict[str, np.ndarray | str]]:
    """Read one LIDO/model R_CP graph from Systematics/lido_model_rcp.root.

    The model is intentionally kept separate from the measured data. It is used
    only as an optional visual comparison in plot_rcp_pages().
    """
    if lido_f is None:
        return None

    lido_dir = lido_rcp_dirs.get(rcp_type)
    obs_key = lido_observable_keys.get(iLamb)
    if not lido_dir or not obs_key:
        return None

    graph_path = f"{lido_dir}/g{lido_dir}_{obs_key}"
    obj = get_root_object(lido_f, graph_path)
    if obj is not None:
        x, y, exl, exh, eyl, eyh = read_tgraph_asymm(obj)
        return {
            "x": x,
            "y": y,
            "exl": exl,
            "exh": exh,
            "stat_low": eyl,
            "stat_high": eyh,
        }

    # Fallback for explicit-binning observables, where plot_lido_rcp_v3.py also
    # writes a TH1D. This is mainly useful if the graph was removed by hand.
    hist_path = f"{lido_dir}/h{lido_dir}_{obs_key}"
    hobj = get_root_object(lido_f, hist_path)
    if hobj is not None:
        x, y, xerr, yerr = read_hist_values(hobj)
        return {
            "x": x,
            "y": y,
            "exl": xerr,
            "exh": xerr,
            "stat_low": yerr,
            "stat_high": yerr,
        }

    print(f"[WARN] Missing LIDO model graph: {graph_path}")
    return None


def draw_lido_model(
    ax: Any,
    lido_res: Optional[Dict[str, np.ndarray]],
    keep_xlim=None,
    x_range: Optional[Tuple[float, float]] = None,
) -> bool:
    """Draw LIDO/model R_CP as one continuous blue statistical band.

    If x_range is given, the LIDO curve and band are clipped to that range.
    Boundary points are added by interpolation, so the filled band ends exactly
    at the requested limits instead of at the nearest LIDO point.
    """
    if lido_res is None:
        return False

    x = np.asarray(lido_res.get("x", []), dtype=float)
    y = np.asarray(lido_res.get("y", []), dtype=float)

    stat_low = np.asarray(
        lido_res.get("stat_low", lido_res.get("eyl", np.zeros_like(y))),
        dtype=float,
    )
    stat_high = np.asarray(
        lido_res.get("stat_high", lido_res.get("eyh", np.zeros_like(y))),
        dtype=float,
    )

    n = min(len(x), len(y), len(stat_low), len(stat_high))
    if n == 0:
        return False

    x = x[:n]
    y = y[:n]
    stat_low = stat_low[:n]
    stat_high = stat_high[:n]

    # First clean and sort the model points.  This must happen before clipping,
    # because interpolation expects a monotonic x array and y_low/y_high must
    # already exist.
    good = (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(stat_low)
        & np.isfinite(stat_high)
    )
    if not np.any(good):
        return False

    x = x[good]
    y = y[good]
    stat_low = stat_low[good]
    stat_high = stat_high[good]

    order = np.argsort(x)
    x = x[order]
    y = y[order]
    stat_low = stat_low[order]
    stat_high = stat_high[order]

    y_low = y - stat_low
    y_high = y + stat_high

    if x_range is not None:
        xmin, xmax = map(float, x_range)
        if xmax <= xmin:
            return False

        # No overlap between LIDO range and requested drawing range.
        if len(x) == 0 or xmax < x[0] or xmin > x[-1]:
            return False

        xmin_clip = max(xmin, float(x[0]))
        xmax_clip = min(xmax, float(x[-1]))
        if xmax_clip <= xmin_clip:
            return False

        mask_range = (x >= xmin_clip) & (x <= xmax_clip)

        x_new = []
        y_new = []
        ylow_new = []
        yhigh_new = []

        # Add left boundary by interpolation when it is not already present.
        if len(x) >= 2 and xmin_clip > x[0] and xmin_clip < x[-1]:
            x_new.append(xmin_clip)
            y_new.append(np.interp(xmin_clip, x, y))
            ylow_new.append(np.interp(xmin_clip, x, y_low))
            yhigh_new.append(np.interp(xmin_clip, x, y_high))

        x_new.extend(x[mask_range])
        y_new.extend(y[mask_range])
        ylow_new.extend(y_low[mask_range])
        yhigh_new.extend(y_high[mask_range])

        # Add right boundary by interpolation when it is not already present.
        if len(x) >= 2 and xmax_clip > x[0] and xmax_clip < x[-1]:
            x_new.append(xmax_clip)
            y_new.append(np.interp(xmax_clip, x, y))
            ylow_new.append(np.interp(xmax_clip, x, y_low))
            yhigh_new.append(np.interp(xmax_clip, x, y_high))

        x = np.asarray(x_new, dtype=float)
        y = np.asarray(y_new, dtype=float)
        y_low = np.asarray(ylow_new, dtype=float)
        y_high = np.asarray(yhigh_new, dtype=float)

        # Remove possible duplicate boundary points and sort again.
        if len(x) < 2:
            return False
        order = np.argsort(x)
        x = x[order]
        y = y[order]
        y_low = y_low[order]
        y_high = y_high[order]
        _, unique_idx = np.unique(x, return_index=True)
        unique_idx = np.sort(unique_idx)
        x = x[unique_idx]
        y = y[unique_idx]
        y_low = y_low[unique_idx]
        y_high = y_high[unique_idx]
        if len(x) < 2:
            return False

    ax.fill_between(
        x,
        y_low,
        y_high,
        color="tab:blue",
        alpha=0.25,
        linewidth=0.0,
        zorder=0.6,
        label="LIDO stat. unc.",
    )

    ax.plot(
        x,
        y,
        color="tab:blue",
        linewidth=1.5,
        linestyle="-",
        zorder=0.8,
    )

    if keep_xlim is not None:
        ax.set_xlim(keep_xlim)

    return True

def apply_edge_bin_visibility(res: Dict[str, np.ndarray | str], iLamb: int) -> Dict[str, np.ndarray | str]:
    """Hide the first and/or last bin for a selected variable.

    Controlled by SHOW_EDGE_BINS at the top of the file.  The function returns
    a shallow copy with all array-like fields sliced consistently.
    """
    cfg = SHOW_EDGE_BINS.get(iLamb, {"first": True, "last": True})
    show_first = bool(cfg.get("first", True))
    show_last = bool(cfg.get("last", True))

    x = np.asarray(res.get("x", []), dtype=float)
    n = len(x)
    if n == 0:
        return res

    mask = np.ones(n, dtype=bool)
    if not show_first:
        mask[0] = False
    if not show_last and n > 1:
        mask[-1] = False

    out: Dict[str, np.ndarray | str] = dict(res)
    for key in ["x", "y", "exl", "exh", "stat_low", "stat_high", "syst_low", "syst_high"]:
        val = out.get(key)
        if isinstance(val, np.ndarray) and len(val) == n:
            out[key] = val[mask]
    return out



def br_relative_uncertainty() -> float:
    """Relative uncertainty from Br(D0 -> K pi)."""
    if BR_RATIO_D0 <= 0 or BR_RATIO_D0_ERR < 0:
        return 0.0
    return float(BR_RATIO_D0_ERR / BR_RATIO_D0)


def add_relative_uncertainty_to_result(
    res: Dict[str, np.ndarray | str],
    rel_unc: float,
) -> Dict[str, np.ndarray | str]:
    """Add a symmetric multiplicative uncertainty to syst_low/high.

    For a quantity y divided by a normalization factor B, the relative
    uncertainty from B is sigma_B/B.  This function adds (rel_unc*y) in
    quadrature to both systematic sides.
    """
    if not math.isfinite(float(rel_unc)) or rel_unc <= 0:
        return res

    out: Dict[str, np.ndarray | str] = dict(res)
    y = np.asarray(out.get("y", []), dtype=float)
    syst_low = np.asarray(out.get("syst_low", np.zeros_like(y)), dtype=float)
    syst_high = np.asarray(out.get("syst_high", np.zeros_like(y)), dtype=float)

    add = np.abs(y) * float(rel_unc)
    n = min(len(y), len(syst_low), len(syst_high), len(add))
    syst_low = syst_low[:n]
    syst_high = syst_high[:n]
    add = add[:n]
    out["syst_low"] = np.sqrt(np.clip(syst_low * syst_low + add * add, 0.0, None))
    out["syst_high"] = np.sqrt(np.clip(syst_high * syst_high + add * add, 0.0, None))

    # Keep all other arrays consistent if some pathological input had a
    # different length.
    for key in ["x", "y", "exl", "exh", "stat_low", "stat_high"]:
        val = out.get(key)
        if isinstance(val, np.ndarray) and len(val) > n:
            out[key] = val[:n]
    return out


def apply_spectrum_global_uncertainties(res: Dict[str, np.ndarray | str]) -> Dict[str, np.ndarray | str]:
    """Apply global normalization uncertainties relevant for spectra."""
    if APPLY_BR_UNCERTAINTY_TO_SPECTRA:
        res = add_relative_uncertainty_to_result(res, br_relative_uncertainty())
    return res


def apply_rcp_global_uncertainties(res: Dict[str, np.ndarray | str]) -> Dict[str, np.ndarray | str]:
    """Optional global uncertainties for R_CP.

    The branching-ratio uncertainty is off by default because it cancels in
    R_CP.  N_coll is drawn separately as a normalization box and is not folded
    into the bin-by-bin systematic boxes here.
    """
    if APPLY_BR_UNCERTAINTY_TO_RCP:
        res = add_relative_uncertainty_to_result(res, br_relative_uncertainty())
    return res


def add_sys_boxes(
    ax: Any,
    x: np.ndarray,
    y: np.ndarray,
    exl: np.ndarray,
    exh: np.ndarray,
    sys_low: np.ndarray,
    sys_high: np.ndarray,
    alpha: float = 0.22,
    zorder: int = 1,
    color: str = "0.5",
    skip_nonpositive: bool = False,
) -> None:
    """Draw systematic uncertainty as rectangles with asymmetric y-errors."""
    for xi, yi, xl, xh, yl, yh in zip(x, y, exl, exh, sys_low, sys_high):
        if not all(math.isfinite(float(v)) for v in [xi, yi, xl, xh, yl, yh]):
            continue
        if yl <= 0 and yh <= 0:
            continue
        y0 = yi - yl
        height = yl + yh
        if skip_nonpositive and y0 <= 0:
            continue
        rect = Rectangle(
            (xi - xl, y0),
            xl + xh,
            height,
            facecolor=color,
            edgecolor="none",
            alpha=alpha,
            zorder=zorder,
        )
        ax.add_patch(rect)




def rcp_ncoll_relative_uncertainty(rcp_type: str) -> float:
    """Relative global N_coll uncertainty for a given R_CP centrality pair."""
    if rcp_type not in rcp_cent_indices:
        return 0.0
    i_num, i_den = rcp_cent_indices[rcp_type]
    rel_num = ncoll_errors[i_num] / ncoll_values[i_num]
    rel_den = ncoll_errors[i_den] / ncoll_values[i_den]
    return float(np.sqrt(rel_num * rel_num + rel_den * rel_den))


def add_ncoll_box(
    ax: Any,
    rcp_type: str,
    *,
    y_center: float = 1.0,
    x_frac: float = 0.965,
    width_frac: float = 0.035,
    color: str = "teal",
    alpha: float = 0.45,
    zorder: int = 2,
) -> Optional[Patch]:
    """Draw the global N_coll scale uncertainty as a small box at R_CP=1.

    The x-position is in axes coordinates, so the box remains at the left side
    of the plot independent of the x range.  The y-size is in data units.
    """
    rel = rcp_ncoll_relative_uncertainty(rcp_type)
    if not math.isfinite(rel) or rel <= 0:
        return None
    rect = Rectangle(
        (x_frac, y_center - rel),
        width_frac,
        2.0 * rel,
        transform=ax.get_yaxis_transform(),
        facecolor=color,
        edgecolor="none",
        alpha=alpha,
        zorder=zorder,
        clip_on=False,
    )
    ax.add_patch(rect)
    return rect


def ncoll_legend_patch() -> Patch:
    return Patch(facecolor="teal", edgecolor="none", alpha=0.45)


def legend_header_labels() -> list:
    return [
        r"Au+Au $\sqrt{s_{\mathrm{NN}}}=200$ GeV",
        "Inclusive D$^{0}$-meson-tagged jets",
        r"Anti-$k_{\mathrm{T}}$ algorithm, R = 0.4",
        jet_pt_range_label,
        d0_pt_range_label,
    ]


def add_header_text(ax: Any, fontsize: float = 8.5) -> None:
    """Place the common preliminary header directly in the plot.

    This keeps the angularity overlay legend compact and horizontal.
    """
    ax.text(
        0.72,
        0.97,
        "\n".join(legend_header_labels()),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=fontsize,
        linespacing=1.25,
    )
    
def add_preliminary_label(
    ax: Any,
    x: Optional[float] = None,
    y: Optional[float] = None,
    ha: Optional[str] = None,
    va: Optional[str] = None,
) -> None:
    if not SHOW_PRELIMINARY_LABEL:
        return

    if x is None:
        x = PRELIMINARY_LABEL_X
    if y is None:
        y = PRELIMINARY_LABEL_Y
    if ha is None:
        ha = "left"
    if va is None:
        va = "top"

    ax.text(
        x,
        y,
        PRELIMINARY_LABEL,
        transform=ax.transAxes,
        ha=ha,
        va=va,
        fontsize=PRELIMINARY_LABEL_FONTSIZE,
        color=PRELIMINARY_LABEL_COLOR,
        fontstyle="italic",
        fontweight="bold",
    )

def add_rcp_overlay_header_text(ax: Any) -> None:
    """Place the common header in the lower-right corner for R_CP overlays only."""
    ax.text(
        RCP_OVERLAY_HEADER_X,
        RCP_OVERLAY_HEADER_Y,
        "\n".join(legend_header_labels()),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=RCP_OVERLAY_HEADER_FONTSIZE,
        linespacing=1.25,
    )


def set_horizontal_legend(ax: Any, handles: list, labels: list, ncol: Optional[int] = None) -> None:
    """Place a compact horizontal legend inside the plotting area.

    Previous versions moved the angularity legend below the axes.  This version
    keeps the legend inside the plot and only changes the internal layout of the
    legend entries using ncol.
    """
    if ncol is None:
        ncol = OVERLAY_LEGEND_NCOL
    ax.legend(
        handles,
        labels,
        frameon=False,
        loc=OVERLAY_LEGEND_LOC,
        bbox_to_anchor=OVERLAY_LEGEND_BBOX,
        bbox_transform=ax.transAxes,
        fontsize=OVERLAY_LEGEND_FONTSIZE,
        ncol=max(1, int(ncol)),
        columnspacing=1.0,
        handlelength=2.0,
        handletextpad=0.45,
        borderaxespad=0.0,
    )


def set_rcp_overlay_legend(ax: Any, handles: list, labels: list, ncol: Optional[int] = None) -> None:
    """Place the R_CP overlay legend in the lower-right corner only."""
    if ncol is None:
        ncol = RCP_OVERLAY_LEGEND_NCOL
    ax.legend(
        handles,
        labels,
        frameon=False,
        loc=RCP_OVERLAY_LEGEND_LOC,
        bbox_to_anchor=RCP_OVERLAY_LEGEND_BBOX,
        bbox_transform=ax.transAxes,
        fontsize=RCP_OVERLAY_LEGEND_FONTSIZE,
        ncol=max(1, int(ncol)),
        columnspacing=1.0,
        handlelength=2.0,
        handletextpad=0.45,
        borderaxespad=0.0,
    )


def legend_header_handles_and_labels() -> Tuple[list, list]:
    labels = legend_header_labels()
    handles = [Line2D([], [], linestyle="none") for _ in labels]
    return handles, labels


def apply_preliminary_layout(fig: Any) -> None:
    try:
        fig.tight_layout()
    except Exception as e:
        print(f"[WARN] tight_layout failed: {e}")
        fig.subplots_adjust(left=0.16, right=0.98, bottom=0.16, top=0.90)


def set_logy_with_margin(ax: Any, y_values: Iterable[np.ndarray], low_errors: Iterable[np.ndarray], factor: float = 10.0) -> None:
    ymin_candidates = []
    ymax_candidates = []
    for y, lo in zip(y_values, low_errors):
        y = np.asarray(y, dtype=float)
        lo = np.asarray(lo, dtype=float)
        good = np.isfinite(y) & (y > 0)
        if np.any(good):
            ymin_candidates.extend(y[good].tolist())
            ymax_candidates.extend(y[good].tolist())
        good_low = np.isfinite(y - lo) & ((y - lo) > 0)
        if np.any(good_low):
            ymin_candidates.extend((y[good_low] - lo[good_low]).tolist())
    ax.set_yscale("log")
    if ymin_candidates and ymax_candidates:
        ymin = min(ymin_candidates)
        ymax = max(ymax_candidates)
        if ymin > 0 and ymax > ymin:
            ax.set_ylim(ymin * 0.5, ymax * factor)


def safe_filename_part(text: Any) -> str:
    """Return a filesystem-safe short tag for individual plot names."""
    out = _to_str(text)
    out = out.replace("–", "-").replace("/", "_over_")
    out = re.sub(r"[^A-Za-z0-9_.+\-]+", "_", out)
    out = re.sub(r"_+", "_", out).strip("_")
    return out or "plot"


def lambda_file_tag(iLamb: int) -> str:
    """Stable short tag for observable names in output filenames."""
    return safe_filename_part(spectrum_job_dir.get(iLamb, f"lambda{iLamb}"))


def cent_file_tag(tag: str) -> str:
    return safe_filename_part(tag.replace("pct", ""))


def rcp_file_tag(rcp_type: str) -> str:
    return safe_filename_part(rcp_tags.get(rcp_type, rcp_type).replace("pct", ""))


def save_individual_figure(fig: Any, all_graphs_dir: Optional[Path], stem: str) -> None:
    """Save one already prepared matplotlib figure as both PDF and PNG.

    This is deliberately additive: the original multi-page PdfPages output is
    written elsewhere exactly as before.  This helper only mirrors the current
    figure into FinalPlot/AllGraphs (or a user-selected directory).
    """
    if all_graphs_dir is None:
        return

    all_graphs_dir.mkdir(parents=True, exist_ok=True)
    clean_stem = safe_filename_part(stem)
    pdf_path = all_graphs_dir / f"{clean_stem}.pdf"
    png_path = all_graphs_dir / f"{clean_stem}.png"

    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")



def plot_spectra_by_observable(f: Any, output_pdf: Path, all_graphs_dir: Optional[Path] = None) -> None:
    with PdfPages(output_pdf) as pdf:
        for iLamb in lambdas_to_plot:
            fig, ax = plt.subplots()
            found_any = False
            xlabel = lambda_labels.get(iLamb, r"$x$")
            ylabel = make_ylabel_from_xlabel(xlabel)
            y_for_limits = []
            lo_for_limits = []

            for iCent, (cent_label, tag) in enumerate(zip(cent_labels, cent_tags)):
                job_dir = spectrum_job_dir[iLamb]
                res = get_final_result(f, job_dir, tag, xlabel=xlabel)
                if res is None:
                    print(f"[WARN] Missing spectrum: {job_dir}/g_value_stat_{tag}")
                    continue
                res = apply_edge_bin_visibility(res, iLamb)
                res = apply_spectrum_global_uncertainties(res)

                x = res["x"]
                y = res["y"]
                exl = res["exl"]
                exh = res["exh"]
                stat_low = res["stat_low"]
                stat_high = res["stat_high"]
                syst_low = res["syst_low"]
                syst_high = res["syst_high"]

                add_sys_boxes(ax, x, y, exl, exh, syst_low, syst_high, alpha=0.22, zorder=1, skip_nonpositive=True)

                color = cent_colors.get(iCent, None)

                add_sys_boxes(
                    ax,
                    x,
                    y,
                    exl,
                    exh,
                    syst_low,
                    syst_high,
                    alpha=0.18,
                    zorder=1,
                    color=color if color is not None else "0.5",
                    skip_nonpositive=True,
                )

                ax.errorbar(
                    x,
                    y,
                    xerr=[exl, exh],
                    yerr=[stat_low, stat_high],
                    fmt=cent_markers.get(iCent, "o"),
                    capsize=3,
                    elinewidth=1,
                    capthick=1,
                    markersize=4,
                    linestyle="none",
                    color=color,
                    label=cent_label,
                    zorder=3,
                )

                found_any = True
                y_for_limits.append(y)
                lo_for_limits.append(syst_low)
                ax.set_xlabel(xlabel, fontsize=AXIS_LABEL_FONTSIZE)
                ax.set_ylabel(ylabel, fontsize=AXIS_LABEL_FONTSIZE + 5, labelpad=0)

                ax.tick_params(axis="both", which="major", labelsize=TICK_LABEL_FONTSIZE)
                ax.tick_params(axis="both", which="minor", labelsize=TICK_LABEL_FONTSIZE)
                xlabel = str(res["xlabel"])
                ylabel = str(res["ylabel"])

            if not found_any:
                plt.close(fig)
                continue

            #ax.set_title(f"{lambda_titles.get(iLamb, f'Lambda{iLamb}')}   ({Method}; final syst.)")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            set_logy_with_margin(ax, y_for_limits, lo_for_limits, factor=200.0)
            
            handles, labels = ax.get_legend_handles_labels()

            sys_patch = MultiColorSysBox(
                [cent_colors[0], cent_colors[1], cent_colors[2]],
                alpha=0.55,
                skew=0.18,
                yshift=-0.20,
            )

            add_spectrum_main_legend(ax, handles, labels, sys_patch)
#            handles, labels = ax.get_legend_handles_labels()
#            header_handles, header_labels = legend_header_handles_and_labels()
#            sys_patch = MultiColorSysBox(
#                [cent_colors[0], cent_colors[1], cent_colors[2]],
#                alpha=0.55,
#                skew=0.18,
#                yshift=-0.15,
#            )
#           legend_handles = header_handles + handles + [sys_patch]
#           legend_labels = header_labels + labels + [r"Syst. unc. "]

#            ax.legend(
#                legend_handles,
#                legend_labels,
#                frameon=False,
#                loc="best",
#                fontsize=9,
#                handlelength=2.8,
#                handleheight=0.9,
#                handler_map={MultiColorSysBox: HandlerMultiColorSysBox()},
#            )
            apply_preliminary_layout(fig)
            add_preliminary_label(ax)
            pdf.savefig(fig)
            save_individual_figure(fig, all_graphs_dir, f"spectra_{lambda_file_tag(iLamb)}")
            plt.close(fig)

    print(f"[OK] Saved: {output_pdf}")


def plot_spectra_overlay_by_centrality(f: Any, output_pdf: Path, all_graphs_dir: Optional[Path] = None) -> None:
    with PdfPages(output_pdf) as pdf:
        for iCent, (cent_label, tag) in enumerate(zip(cent_labels, cent_tags)):
            fig, ax = plt.subplots(figsize=(8.5, 6.0))
            found_any = False
            xlabel = r"$\lambda_{\alpha}^{1}$"
            ylabel = "Counts"
            y_for_limits = []
            lo_for_limits = []
            plotted_lambdas = []

            for iLamb in lambdas_overlay:
                job_dir = spectrum_job_dir[iLamb]
                xl = lambda_labels.get(iLamb, r"$x$")
                res = get_final_result(f, job_dir, tag, xlabel=xl)
                if res is None:
                    print(f"[WARN] Missing overlay spectrum: {job_dir}/g_value_stat_{tag}")
                    continue
                res = apply_edge_bin_visibility(res, iLamb)
                res = apply_spectrum_global_uncertainties(res)

                x = res["x"]
                y = res["y"]
                exl = res["exl"]
                exh = res["exh"]
                stat_low = res["stat_low"]
                stat_high = res["stat_high"]
                syst_low = res["syst_low"]
                syst_high = res["syst_high"]
                color = lambda_colors[iLamb]

                add_sys_boxes(ax, x, y, exl, exh, syst_low, syst_high, alpha=0.18, zorder=1, color=color, skip_nonpositive=True)

                ax.errorbar(
                    x,
                    y,
                    xerr=[exl, exh],
                    yerr=[stat_low, stat_high],
                    fmt="o",
                    capsize=3,
                    elinewidth=1.6,
                    capthick=1.6,
                    markersize=5,
                    linestyle="none",
                    linewidth=0.0,
                    color=color,
                    label=lambda_labels[iLamb],
                    zorder=3,
                )

                if not found_any:
                    xlabel = ANGULARITY_OVERLAY_XLABEL
                    ylabel = str(res["ylabel"])
                found_any = True
                plotted_lambdas.append(iLamb)
                y_for_limits.append(y)
                lo_for_limits.append(syst_low)

            if not found_any:
                plt.close(fig)
                continue

            #ax.set_title(f"{cent_label}   ({Method}; final syst.)")
            ax.set_xlabel(xlabel, fontsize=ANGULARITY_OVERLAY_LABEL_FONTSIZE)
            ax.set_ylabel(ylabel, fontsize=ANGULARITY_OVERLAY_LABEL_FONTSIZE+3)
            ax.tick_params(axis="both", which="major", labelsize=ANGULARITY_OVERLAY_TICK_LABEL_FONTSIZE)
            ax.tick_params(axis="both", which="minor", labelsize=ANGULARITY_OVERLAY_TICK_LABEL_FONTSIZE)
            set_logy_with_margin(ax, y_for_limits, lo_for_limits, factor=2.0)

            handles, labels = ax.get_legend_handles_labels()
            add_spectra_overlay_legend(ax, handles, labels, plotted_lambdas)
            add_spectra_overlay_centrality_legend(ax, cent_label)
            add_preliminary_label(ax)
            try:
                fig.tight_layout()
            except Exception as e:
                print(f"[WARN] tight_layout failed: {e}")
                fig.subplots_adjust(left=0.14, right=0.98, bottom=0.14, top=0.92)
            pdf.savefig(fig)
            save_individual_figure(fig, all_graphs_dir, f"spectra_overlay_{cent_file_tag(tag)}")
            plt.close(fig)

    print(f"[OK] Saved: {output_pdf}")


def plot_rcp_pages(f: Any, output_pdf: Path, lido_f: Optional[Any] = None, show_lido_model: bool = False) -> None:
    with PdfPages(output_pdf) as pdf:
        for iLamb in rcp_lambdas_to_plot:
            fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=False)
            found_any = False
            xlabel = lambda_labels.get(iLamb, r"$x$")

            for ax, rcp_type in zip(axes, rcp_types):
                tag = rcp_tags[rcp_type]
                job_dir = rcp_job_dir[iLamb]
                res = get_final_result(f, job_dir, tag, xlabel=xlabel)
                color = rcp_colors.get(rcp_type, "black")
                marker = rcp_markers.get(rcp_type, "o")
                if res is None:
                    print(f"[WARN] Missing RCP: {job_dir}/g_value_stat_{tag}")
                    #ax.set_title(rcp_titles[rcp_type])
                    ax.text(0.5, 0.5, "Missing graph", ha="center", va="center", transform=ax.transAxes)
                    ax.set_ylim(0, 2.25)
                    continue
                res = apply_edge_bin_visibility(res, iLamb)
                res = apply_rcp_global_uncertainties(res)

                x = res["x"]
                y = res["y"]
                exl = res["exl"]
                exh = res["exh"]
                stat_low = res["stat_low"]
                stat_high = res["stat_high"]
                syst_low = res["syst_low"]
                syst_high = res["syst_high"]

                add_sys_boxes(
                    ax,
                    x,
                    y,
                    exl,
                    exh,
                    syst_low,
                    syst_high,
                    alpha=0.18,
                    zorder=1,
                    color=color,
                    skip_nonpositive=False,
                )
                add_ncoll_box(ax, rcp_type)

                ax.errorbar(
                    x,
                    y,
                    xerr=[exl, exh],
                    yerr=[stat_low, stat_high],
                    fmt=marker,
                    capsize=2.5,
                    elinewidth=1.2,
                    capthick=1.2,
                    markersize=5,
                    linestyle="none",
                    color=color,
                    zorder=3,
                    label=rcp_titles[rcp_type],
                )

                # Store the data-defined x-range before adding the optional
                # model. If LIDO has a longer x-range, it is clipped instead of
                # expanding the measured-data plot.
                data_xlim = ax.get_xlim()

                lido_x_range = None
                manual_range = LIDO_MANUAL_X_RANGE.get(iLamb, None)

                if manual_range is not None:
                    lido_x_range = manual_range
                elif CLIP_LIDO_TO_DATA_RANGE:
                    # Use the actually visible measured-data bin range, not the matplotlib margin.
                    data_xmin = float(np.nanmin(x))
                    data_xmax = float(np.nanmax(x))   # center of the last visible data bin
                    lido_x_range = (data_xmin, data_xmax)

                lido_drawn = False
                if show_lido_model and lido_f is not None:
                    lido_res = get_lido_model_result(lido_f, rcp_type, iLamb)
                    lido_drawn = draw_lido_model(
                        ax,
                        lido_res,
                        keep_xlim=data_xlim,
                        x_range=lido_x_range,
                    )

                ax.set_ylim(0, 4)
                if iLamb not in (0, 6):   # 0 = z, 6 = p_T^D
                    ax.set_ylim(0, 2.4)
                ax.axhline(1.0, linestyle="--", linewidth=1.2, color="black", zorder=0)
                #add_rcp_ratio_label(ax, rcp_type, color)
                #ax.set_title(rcp_titles[rcp_type])
                ax.set_xlabel(str(res["xlabel"]), fontsize=RCP_AXIS_LABEL_FONTSIZE)
                
                ax.set_ylabel(
                rcp_axis_labels.get(rcp_type, r"$R_{\mathrm{CP}}$"),
                fontsize=RCP_AXIS_LABEL_FONTSIZE,
                labelpad=RCP_YLABEL_PAD,
                )
                
                ax.tick_params(axis="both", which="major", labelsize=RCP_TICK_LABEL_FONTSIZE)
                ax.tick_params(axis="both", which="minor", labelsize=RCP_TICK_LABEL_FONTSIZE)
                #if ax is axes[0]:
                #ax.set_ylabel(r"$R_{\mathrm{CP}}$")
                ax.set_ylabel(rcp_axis_labels.get(rcp_type, r"$R_{\mathrm{CP}}$"))

                found_any = True

#                handles = [
#                    Line2D([], [], linestyle="none"),
#                    Line2D([], [], linestyle="none"),
#                    Line2D([], [], linestyle="none"),
#                    Line2D([], [], linestyle="none"),
#                    Line2D([], [], linestyle="none"),
#                    Line2D([], [], marker=marker, linestyle="none", color=color),
#                    Patch(facecolor=color, edgecolor="none", alpha=0.18),
#                    ncoll_legend_patch(),
#                    Line2D([], [], linestyle="--", color="black"),
#                ]
#                labels = [
#                    r"STAR Au+Au, $\sqrt{s_{\mathrm{NN}}}=200$ GeV",
#                    "Inclusive D$^{0}$-meson-tagged jets",
#                    r"Anti-$k_{\mathrm{T}}$ algorithm, R = 0.4",
#                    jet_pt_range_label,
#                    d0_pt_range_label,
#                    "Stat. unc.",
#                    r"Syst. unc.",
#                    r"$N_{\mathrm{coll}}$ unc."
                 #   r"$R_{\mathrm{CP}} = 1$",
#                ]
#                ax.legend(handles, labels, frameon=False, loc="best", fontsize=8)
                add_rcp_page_legend(ax, color, marker, show_lido_model=lido_drawn)
                add_preliminary_label(ax, x=0.96, y=0.04, ha="right", va="bottom")

            if not found_any:
                plt.close(fig)
                continue

            #fig.suptitle(f"{lambda_titles.get(iLamb, f'Lambda{iLamb}')}   ({Method}; final syst.)", y=0.98)
            try:
                fig.tight_layout(rect=[0, 0, 1, 0.95])
            except Exception as e:
                print(f"[WARN] tight_layout failed for RCP page: {e}")
                fig.subplots_adjust(left=0.06, right=0.98, bottom=0.12, top=0.88, wspace=0.28)
            pdf.savefig(fig)
            plt.close(fig)


    print(f"[OK] Saved: {output_pdf}")


def plot_rcp_individual_graphs(
    f: Any,
    all_graphs_dir: Optional[Path],
    lido_f: Optional[Any] = None,
    show_lido_model: bool = False,
) -> None:
    """Save each R_CP panel from results_RCP as a standalone PDF and PNG.

    The regular results_RCP*.pdf stays unchanged with three panels per page.
    This function only creates additional one-panel exports in AllGraphs.
    """
    if all_graphs_dir is None:
        return

    saved = 0
    for iLamb in rcp_lambdas_to_plot:
        xlabel = lambda_labels.get(iLamb, r"$x$")
        job_dir = rcp_job_dir[iLamb]

        for rcp_type in rcp_types:
            tag = rcp_tags[rcp_type]
            res = get_final_result(f, job_dir, tag, xlabel=xlabel)
            color = rcp_colors.get(rcp_type, "black")
            marker = rcp_markers.get(rcp_type, "o")

            if res is None:
                print(f"[WARN] Missing individual RCP: {job_dir}/g_value_stat_{tag}")
                continue

            res = apply_edge_bin_visibility(res, iLamb)
            res = apply_rcp_global_uncertainties(res)

            x = res["x"]
            y = res["y"]
            exl = res["exl"]
            exh = res["exh"]
            stat_low = res["stat_low"]
            stat_high = res["stat_high"]
            syst_low = res["syst_low"]
            syst_high = res["syst_high"]

            fig, ax = plt.subplots(figsize=(8.5, 6.0))

            add_sys_boxes(
                ax,
                x,
                y,
                exl,
                exh,
                syst_low,
                syst_high,
                alpha=0.18,
                zorder=1,
                color=color,
                skip_nonpositive=False,
            )
            add_ncoll_box(ax, rcp_type)

            ax.errorbar(
                x,
                y,
                xerr=[exl, exh],
                yerr=[stat_low, stat_high],
                fmt=marker,
                capsize=2.5,
                elinewidth=1.2,
                capthick=1.2,
                markersize=5,
                linestyle="none",
                color=color,
                zorder=3,
                label=rcp_titles[rcp_type],
            )

            data_xlim = ax.get_xlim()

            lido_x_range = None
            manual_range = LIDO_MANUAL_X_RANGE.get(iLamb, None)
            if manual_range is not None:
                lido_x_range = manual_range
            elif CLIP_LIDO_TO_DATA_RANGE and len(x) > 0:
                data_xmin = float(np.nanmin(x))
                data_xmax = float(np.nanmax(x))
                lido_x_range = (data_xmin, data_xmax)

            lido_drawn = False
            if show_lido_model and lido_f is not None:
                lido_res = get_lido_model_result(lido_f, rcp_type, iLamb)
                lido_drawn = draw_lido_model(
                    ax,
                    lido_res,
                    keep_xlim=data_xlim,
                    x_range=lido_x_range,
                )

            ax.set_ylim(0, 4)
            if iLamb not in (0, 6):
                ax.set_ylim(0, 2.4)

            ax.axhline(1.0, linestyle="--", linewidth=1.2, color="black", zorder=0)
            ax.set_xlabel(str(res["xlabel"]), fontsize=RCP_AXIS_LABEL_FONTSIZE)
            ax.set_ylabel(
                rcp_axis_labels.get(rcp_type, r"$R_{\mathrm{CP}}$"),
                fontsize=RCP_AXIS_LABEL_FONTSIZE,
                labelpad=RCP_YLABEL_PAD,
            )
            ax.tick_params(axis="both", which="major", labelsize=RCP_TICK_LABEL_FONTSIZE)
            ax.tick_params(axis="both", which="minor", labelsize=RCP_TICK_LABEL_FONTSIZE)

            add_rcp_page_legend(ax, color, marker, show_lido_model=lido_drawn)
            add_preliminary_label(ax, x=0.96, y=0.04, ha="right", va="bottom")

            try:
                fig.tight_layout()
            except Exception as e:
                print(f"[WARN] tight_layout failed for individual RCP: {e}")
                fig.subplots_adjust(left=0.14, right=0.98, bottom=0.14, top=0.92)

            save_individual_figure(
                fig,
                all_graphs_dir,
                f"rcp_{lambda_file_tag(iLamb)}_{rcp_file_tag(rcp_type)}",
            )
            saved += 1
            plt.close(fig)

    if saved:
        print(f"[OK] Saved individual R_CP graphs to: {all_graphs_dir}")



def plot_rcp_overlay(f: Any, output_pdf: Path, all_graphs_dir: Optional[Path] = None) -> None:
    with PdfPages(output_pdf) as pdf:
        for rcp_type in rcp_types:
            fig, ax = plt.subplots(figsize=(8.5, 6.0))
            found_any = False
            xlabel = r"$x$"
            tag = rcp_tags[rcp_type]

            for iLamb in rcp_overlay_lambdas:
                job_dir = rcp_job_dir[iLamb]
                xl = lambda_labels.get(iLamb, r"$x$")
                res = get_final_result(f, job_dir, tag, xlabel=xl)
                if res is None:
                    print(f"[WARN] Missing RCP overlay: {job_dir}/g_value_stat_{tag}")
                    continue
                res = apply_edge_bin_visibility(res, iLamb)
                res = apply_rcp_global_uncertainties(res)

                x = res["x"]
                y = res["y"]
                exl = res["exl"]
                exh = res["exh"]
                stat_low = res["stat_low"]
                stat_high = res["stat_high"]
                syst_low = res["syst_low"]
                syst_high = res["syst_high"]
                color = lambda_colors[iLamb]

                add_sys_boxes(
                    ax,
                    x,
                    y,
                    exl,
                    exh,
                    syst_low,
                    syst_high,
                    alpha=0.18,
                    zorder=1,
                    color=color,
                    skip_nonpositive=False,
                )

                ax.errorbar(
                    x,
                    y,
                    xerr=[exl, exh],
                    yerr=[stat_low, stat_high],
                    fmt="o",
                    capsize=3,
                    elinewidth=1.6,
                    capthick=1.6,
                    markersize=5,
                    linestyle="none",
                    linewidth=1.6,
                    color=color,
                    label=lambda_labels[iLamb],
                    zorder=3,
                )

                if not found_any:
                    xlabel = ANGULARITY_OVERLAY_XLABEL
                found_any = True

            if not found_any:
                plt.close(fig)
                continue

            ax.axhline(1.0, linestyle="--", linewidth=1.2, color="black", zorder=0)
            add_ncoll_box(ax, rcp_type)
            ax.set_ylim(0, 1.8)
            #ax.set_title(f"{rcp_titles[rcp_type]}   ({Method})")
            ax.set_xlabel(xlabel, fontsize=RCP_AXIS_LABEL_FONTSIZE+2)
            ax.set_ylabel(
                rcp_axis_labels.get(rcp_type, r"$R_{\mathrm{CP}}$"),
                fontsize=RCP_AXIS_LABEL_FONTSIZE+2,
                labelpad=RCP_YLABEL_PAD,
            )
                
            ax.tick_params(axis="both", which="major", labelsize=RCP_TICK_LABEL_FONTSIZE)
            ax.tick_params(axis="both", which="minor", labelsize=RCP_TICK_LABEL_FONTSIZE)
                #if ax is axes[0]:
                #ax.set_ylabel(r"$R_{\mathrm{CP}}$")
            #ax.set_ylabel(rcp_axis_labels.get(rcp_type, r"$R_{\mathrm{CP}}$"))

            handles, labels = ax.get_legend_handles_labels()

            add_rcp_overlay_legend(ax, handles, labels)

            add_preliminary_label(
                ax,
                x=0.02,
                y=0.97,
                ha="left",
                va="top",
            )

            try:
                fig.tight_layout()
            except Exception as e:
                print(f"[WARN] tight_layout failed: {e}")
                fig.subplots_adjust(left=0.14, right=0.98, bottom=0.14, top=0.92)
            pdf.savefig(fig)
            save_individual_figure(fig, all_graphs_dir, f"rcp_overlay_{rcp_file_tag(rcp_type)}")
            plt.close(fig)

    print(f"[OK] Saved: {output_pdf}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot final results from the fresh ROOT file produced by systematics_final_combiner_v7.py")
    p.add_argument("--root-file", default=DEFAULT_ROOT_FILE, help="Final ROOT file from systematics_final_combiner_v7.py")
    p.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR, help="Directory for output PDFs")
    p.add_argument(
        "--all-graphs-dir",
        default=None,
        help=(
            "Directory for individual PDF/PNG graphs. "
            "Default: <output-dir>/AllGraphs"
        ),
    )
    p.add_argument(
        "--no-all-graphs",
        action="store_true",
        help="Disable saving individual PDF/PNG graphs.",
    )
    p.add_argument("--method-label", default=Method, help="Text label used in plot titles")
    p.add_argument(
        "--show-lido",
        action="store_true",
        default=SHOW_LIDO_MODEL_DEFAULT,
        help="Draw LIDO/model comparison in results_RCP_v9.pdf only",
    )
    p.add_argument(
        "--lido-root-file",
        default=DEFAULT_LIDO_ROOT_FILE,
        help="ROOT file with LIDO/model R_CP graphs",
    )
    return p.parse_args()


def main() -> int:
    global Method
    args = parse_args()
    Method = args.method_label

    root_file = Path(args.root_file)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.no_all_graphs:
        all_graphs_dir = None
    elif args.all_graphs_dir is not None:
        all_graphs_dir = Path(args.all_graphs_dir)
    else:
        all_graphs_dir = out_dir / DEFAULT_ALL_GRAPHS_SUBDIR

    if all_graphs_dir is not None:
        all_graphs_dir.mkdir(parents=True, exist_ok=True)

    if not root_file.exists():
        raise FileNotFoundError(f"Missing input ROOT file: {root_file}")
    if uproot is None:
        raise RuntimeError("This plotter needs the Python package 'uproot' to read the final ROOT file.")

    output_pdf = out_dir / f"results{VERSION_TAG}.pdf"
    output_pdf_overlay = out_dir / f"results_by_centrality{VERSION_TAG}.pdf"
    output_pdf_rcp = out_dir / f"results_RCP{VERSION_TAG}.pdf"
    output_pdf_rcp_overlay = out_dir / f"results_RCP_overlay{VERSION_TAG}.pdf"

    lido_root_file = Path(args.lido_root_file)
    use_lido = bool(args.show_lido)
    if use_lido and not lido_root_file.exists():
        print(f"[WARN] --show-lido was requested, but the LIDO ROOT file is missing: {lido_root_file}")
        print("[WARN] Continuing without LIDO/model comparison.")
        use_lido = False

    with uproot.open(str(root_file)) as f:
        if use_lido:
            with uproot.open(str(lido_root_file)) as lido_f:
                plot_spectra_by_observable(f, output_pdf, all_graphs_dir=all_graphs_dir)
                plot_spectra_overlay_by_centrality(f, output_pdf_overlay, all_graphs_dir=all_graphs_dir)
                plot_rcp_pages(f, output_pdf_rcp, lido_f=lido_f, show_lido_model=True)
                plot_rcp_individual_graphs(f, all_graphs_dir, lido_f=lido_f, show_lido_model=True)
                plot_rcp_overlay(f, output_pdf_rcp_overlay, all_graphs_dir=all_graphs_dir)
        else:
            plot_spectra_by_observable(f, output_pdf, all_graphs_dir=all_graphs_dir)
            plot_spectra_overlay_by_centrality(f, output_pdf_overlay, all_graphs_dir=all_graphs_dir)
            plot_rcp_pages(f, output_pdf_rcp)
            plot_rcp_individual_graphs(f, all_graphs_dir)
            plot_rcp_overlay(f, output_pdf_rcp_overlay, all_graphs_dir=all_graphs_dir)

    if all_graphs_dir is not None:
        print(f"[OK] Individual graphs directory: {all_graphs_dir}")
    print("[OK] Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
