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
  python3 finalPlot_from_systematics_root_v12_lido_spectra.py --show-lido

This reads Systematics/lido_model_rcp.root by default.  In addition to the
legacy outputs, the script writes one configurable stacked panel with spectra
above one separate row per selected R_CP and a shared x-axis in each
selected-variable column.
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
lambdas_overlay = [5, 1]
rcp_lambdas_to_plot = [5, 1, 2, 3, 4]
rcp_overlay_lambdas = [5, 2, 4]

# --- Combined spectra + R_CP panel ---
# One selected variable = one column.  The upper row contains the three
# centrality spectra.  Every selected R_CP ratio gets its own row below it.
#
# Variable indices:
#   0 = z, 1 = lambda_1^1, 2 = lambda_1.5^1, 3 = lambda_2^1,
#   4 = lambda_3^1, 5 = lambda_0.5^1, 6 = p_T^D
COMBINED_PANEL_VARIABLES = [5, 1, 2, 3, 4]

# R_CP selection indices:
#   0 = 0-10% / 10-40%
#   1 = 0-10% / 40-80%
#   2 = 10-40% / 40-80%
# Entries can also be given directly as "0/1", "0/2", and "1/2".
COMBINED_PANEL_RCP_SELECTION = [0, 1, 2]

# Share the y-scale independently within the spectra row and the complete R_CP
# block.  The x-axis is always shared by all axes in each variable column.
COMBINED_PANEL_SHARE_SPECTRA_Y = True
COMBINED_PANEL_SHARE_RCP_Y = True
COMBINED_PANEL_LOG_SPECTRA_Y = True

# Optional fixed y-ranges.  Use None for an automatic range.
COMBINED_PANEL_SPECTRA_YLIM = None
COMBINED_PANEL_RCP_YLIM = (0.0, 2.0)
COMBINED_PANEL_SPECTRA_YMAX_FACTOR = 5.0

# Ignore model points outside the final visible data range when determining the
# automatic spectrum y-limits.  They may still exist in the input graph, but
# they must not pull the displayed log scale down if the x-axis clips them.
COMBINED_PANEL_SPECTRA_Y_LIMIT_VISIBLE_X_ONLY = True

# Column titles repeat the observable already shown on the common x-axis.
COMBINED_PANEL_SHOW_COLUMN_TITLES = False

# X labels are retained at both sides of every column.  Their separation is
# provided by the wider columns and automatic x padding below.
COMBINED_PANEL_HIDE_TOUCHING_X_TICK_LABELS = False
COMBINED_PANEL_HIDE_TOUCHING_Y_TICK_LABELS = True

# Keep the lowest spectrum y label.  At the spectra/R_CP boundary, the upper
# label of the first R_CP row is hidden instead.  R_CP-to-R_CP boundaries keep
# the previous behavior.
COMBINED_PANEL_KEEP_SPECTRA_BOTTOM_Y_TICK_LABEL = True

# Optional manual x-ranges keyed by the variable index.  Unlisted variables
# use the visible measured bins.  Example: 5: (0.0, 0.45)
COMBINED_PANEL_X_RANGES = {
}

# Extra space around automatically determined x-ranges keeps the first visible
# tick away from a touching column boundary.  Manual ranges above are used
# exactly as written.
COMBINED_PANEL_X_PADDING_FRACTION = 0.05

# Geometry is calculated from the numbers of selected variables and R_CP rows.
COMBINED_PANEL_COLUMN_WIDTH = 5.2
COMBINED_PANEL_HEIGHT_PER_RATIO_UNIT = 2.30
COMBINED_PANEL_SPECTRA_HEIGHT_RATIO = 1.35
COMBINED_PANEL_RCP_HEIGHT_RATIO = 1.0

# Zero gives a continuous panel with touching axes.
COMBINED_PANEL_WSPACE = 0.0
COMBINED_PANEL_HSPACE = 0.0

# Draw tick marks on the bottom, left and right sides of every cell.  The top
# side is intentionally excluded; numeric labels remain controlled separately.
COMBINED_PANEL_TICKS_ON_ALL_PANELS = True
COMBINED_PANEL_TICK_DIRECTION = "in"

# LIDO is drawn in every available upper/lower panel when --show-lido is used.
# The upper-row drawing additionally respects --no-lido-spectra.
COMBINED_PANEL_SHOW_LIDO_SPECTRA = True
COMBINED_PANEL_SHOW_LIDO_RCP = True

# Legend placement.  COLUMN is zero-based; negative values count from the
# right.  Set it to None to hide the corresponding legend.
COMBINED_PANEL_SPECTRA_LEGEND_COLUMN = 0
COMBINED_PANEL_SPECTRA_LEGEND_LOC = "upper right"
COMBINED_PANEL_SPECTRA_LEGEND_BBOX = (0.98, 0.98)
COMBINED_PANEL_SPECTRA_LEGEND_NCOL = 2
COMBINED_PANEL_SPECTRA_LEGEND_FONTSIZE = 11.0

COMBINED_PANEL_RCP_LEGEND_COLUMN = 0
# The centrality-pair entry is drawn in every selected R_CP row.  ROW chooses
# Every row additionally contains its own matching-color Syst. and optional
# LIDO entries.  ROW chooses where the common relative N_coll explanation is
# added.  Indices are zero-based; negative values count from the bottom.  Set
# COLUMN to None to hide all R_CP legends, or ROW to None to omit only the
# N_coll explanation.
COMBINED_PANEL_RCP_LEGEND_ROW = 0
COMBINED_PANEL_RCP_LEGEND_LOC = "upper right"
COMBINED_PANEL_RCP_LEGEND_BBOX = (0.98, 0.98)
COMBINED_PANEL_RCP_LEGEND_NCOL = 2
COMBINED_PANEL_RCP_LEGEND_FONTSIZE = 11.0

# Common experiment header and preliminary label.  Their coordinates are in
# the selected axes coordinates, so they can be moved independently of the
# legends above.
COMBINED_PANEL_HEADER_COLUMN = -1
COMBINED_PANEL_HEADER_X = 0.97
COMBINED_PANEL_HEADER_Y = 0.97
COMBINED_PANEL_HEADER_HA = "right"
COMBINED_PANEL_HEADER_VA = "top"
COMBINED_PANEL_HEADER_FONTSIZE = 11.0

COMBINED_PANEL_PRELIMINARY_COLUMN = -1
# Zero-based index within the selected R_CP rows; -1 is the bottom row.  If no
# R_CP is selected, the label is placed in the spectra row instead.
COMBINED_PANEL_PRELIMINARY_RCP_ROW = -1
COMBINED_PANEL_PRELIMINARY_X = 0.97
COMBINED_PANEL_PRELIMINARY_Y = 0.04
COMBINED_PANEL_PRELIMINARY_HA = "right"
COMBINED_PANEL_PRELIMINARY_VA = "bottom"

COMBINED_PANEL_AXIS_LABEL_FONTSIZE = 13
COMBINED_PANEL_TICK_LABEL_FONTSIZE = 11
COMBINED_PANEL_TITLE_FONTSIZE = 14

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

# LIDO spectra were added to lido_model_rcp.root by make_lido_rcp_from_spectra_v7.py.
# Example path:
#   Spectra_0_10/gSpectra_0_10_lambda_0d5
lido_spectra_dirs = {
    0: "Spectra_0_10",
    1: "Spectra_10_40",
    2: "Spectra_40_80",
}

# Draw LIDO spectra in results_by_centrality*.pdf when --show-lido is used.
# The default scale is intentionally 1.0 for a first direct check of the model output.
# Later, if needed, this can absorb factors such as 2, 2*pi*x, BR, or Ncoll.
LIDO_SPECTRA_Y_SCALE = 1.0
LIDO_SPECTRA_MULTIPLY_BY_NCOLL = False
LIDO_SPECTRA_LINESTYLE = "-"
LIDO_SPECTRA_LINEWIDTH = 1.8
LIDO_SPECTRA_BAND_ALPHA = 0.16
LIDO_SPECTRA_LINE_ALPHA = 0.95

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

RCP_OVERLAY_MARKERS = {
    5: "o",  # lambda_0.5^1
    1: "s",  # lambda_1^1
    2: "X",  # lambda_1.5^1
    3: "P",  # lambda_2^1
    4: "D",  # lambda_3^1
}

# Separate legend placement for results_by_centrality*.pdf.
# It is intentionally independent of the R_CP overlay settings.
SPECTRA_OVERLAY_HEADER_X = 0.45
SPECTRA_OVERLAY_HEADER_Y = 0.97
SPECTRA_OVERLAY_HEADER_FONTSIZE = 16

SPECTRA_OVERLAY_ANG_X = 0.44
SPECTRA_OVERLAY_ANG_Y = 0.70
SPECTRA_OVERLAY_ANG_FONTSIZE = 16

SPECTRA_OVERLAY_UNC_X = 0.44
SPECTRA_OVERLAY_UNC_Y = 0.64
SPECTRA_OVERLAY_UNC_FONTSIZE = 16

ANGULARITY_OVERLAY_LABEL_FONTSIZE = 14
ANGULARITY_OVERLAY_TICK_LABEL_FONTSIZE = 14

SPECTRA_OVERLAY_CENT_X = 0.96
SPECTRA_OVERLAY_CENT_Y = 0.09
SPECTRA_OVERLAY_CENT_FONTSIZE = 19


# --- LIDO spectra normalization switches ---
LIDO_REMOVE_X_FACTOR = False

# 0.0 = no Ncoll scaling, as in Lin's direct ratio
# 1.0 = binary-collision scaling for absolute Au+Au per-event spectra
# values between 0 and 1 can be used as a diagnostic
LIDO_NCOLL_POWER = 1

# Global centrality-independent scale.
# Use this only for convention factors such as 1/(2*pi), 1/2, BR, etc.
#LIDO_GLOBAL_NORM = 1/43.82#/1.2
LIDO_GLOBAL_NORM = 1/42#/1.2
# -----------------------------------------------------
def lido_spectrum_total_factor(x: np.ndarray, iCent: int) -> np.ndarray:
    """Centrality-dependent conversion factor for LIDO spectra.

    LIDO spectra are given as:
        1/N_evt d^2N / (X dX d eta)

    This function converts them to the plotting convention and optionally
    applies binary-collision scaling.
    """
    x = np.asarray(x, dtype=float)

    if LIDO_REMOVE_X_FACTOR:
        x_factor = x
    else:
        x_factor = np.ones_like(x)

    if 0 <= iCent < len(ncoll_values):
        ncoll_factor = float(ncoll_values[iCent]) ** float(LIDO_NCOLL_POWER)
    else:
        ncoll_factor = 1.0

    return x_factor * ncoll_factor * float(LIDO_GLOBAL_NORM)


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
    show_lido_model: bool = False,
) -> None:
    """Split the results_by_centrality legend into angularities and uncertainties.

    If LIDO spectra are drawn, add one generic LIDO entry. The LIDO curves are
    colored by angularity, so the legend entry describes the style rather than
    repeating every color.
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

    unc_handles = [sys_patch]
    unc_labels = [r"$\mathrm{Syst.\ unc.}$"]

    if show_lido_model:
        unc_handles.append(
            LidoBandLegend(
                facecolor="0.5",
                alpha=LIDO_SPECTRA_BAND_ALPHA,
                linecolor="0.2",
                linewidth=LIDO_SPECTRA_LINEWIDTH,
            )
        )
        unc_labels.append(LIDO_MODEL_LABEL)

    ax.legend(
        unc_handles,
        unc_labels,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(SPECTRA_OVERLAY_UNC_X, SPECTRA_OVERLAY_UNC_Y),
        bbox_transform=ax.transAxes,
        fontsize=SPECTRA_OVERLAY_UNC_FONTSIZE,
        ncol=len(unc_handles),
        columnspacing=0.9,
        handlelength=2.3,
        handleheight=0.9,
        handletextpad=0.35,
        borderaxespad=0.0,
        handler_map={
            MultiColorSysBox: HandlerMultiColorSysBox(),
            LidoBandLegend: HandlerLidoBandLegend(),
        },
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


class MultiColorLineLegend:
    def __init__(self, colors, linewidth=1.5):
        self.colors = list(colors)
        self.linewidth = linewidth


class HandlerMultiColorLineLegend(HandlerBase):
    def create_artists(
        self, legend, orig_handle,
        xdescent, ydescent, width, height, fontsize, trans
    ):
        colors = orig_handle.colors or ["0.25"]
        n = len(colors)
        ymid = ydescent + 0.5 * height
        artists = []
        for index, color in enumerate(colors):
            x0 = xdescent + width * index / n
            x1 = xdescent + width * (index + 1) / n
            artists.append(
                Line2D(
                    [x0, x1],
                    [ymid, ymid],
                    color=color,
                    linewidth=orig_handle.linewidth,
                    solid_capstyle="butt",
                    transform=trans,
                )
            )
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
    """Return x, y, exlow, exhigh, eylow, eyhigh from TGraphAsymmErrors or TGraphErrors."""
    x = _member_array(obj, "fX")
    y = _member_array(obj, "fY")
    n = min(len(x), len(y))
    x = x[:n]
    y = y[:n]

    # TGraphAsymmErrors
    exl = _member_array(obj, "fEXlow", n=n, default=0.0)
    exh = _member_array(obj, "fEXhigh", n=n, default=0.0)
    eyl = _member_array(obj, "fEYlow", n=n, default=0.0)
    eyh = _member_array(obj, "fEYhigh", n=n, default=0.0)

    # TGraphErrors fallback: symmetric errors are stored as fEX/fEY.
    ex = _member_array(obj, "fEX", n=n, default=0.0)
    ey = _member_array(obj, "fEY", n=n, default=0.0)

    if np.all(exl == 0.0) and np.all(exh == 0.0) and np.any(ex != 0.0):
        exl = ex.copy()
        exh = ex.copy()

    if np.all(eyl == 0.0) and np.all(eyh == 0.0) and np.any(ey != 0.0):
        eyl = ey.copy()
        eyh = ey.copy()

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


def lido_spectrum_scale(iCent: int) -> float:
    """Return optional scale factor for LIDO spectra.

    Kept as 1 by default for a direct visual check. Toggle
    LIDO_SPECTRA_MULTIPLY_BY_NCOLL or change LIDO_SPECTRA_Y_SCALE at the top
    if we later decide that a normalization conversion is needed.
    """
    scale = float(LIDO_SPECTRA_Y_SCALE)
    if LIDO_SPECTRA_MULTIPLY_BY_NCOLL and 0 <= iCent < len(ncoll_values):
        scale *= float(ncoll_values[iCent])
    return scale


def get_lido_spectrum_result(
    lido_f: Optional[Any],
    iCent: int,
    iLamb: int,
) -> Optional[Dict[str, np.ndarray | str]]:
    """Read one LIDO spectrum from lido_model_rcp.root.

    Expected paths from make_lido_rcp_from_spectra_v7.py:
      Spectra_0_10/gSpectra_0_10_lambda_0d5
      Spectra_10_40/gSpectra_10_40_lambda_0d5
      Spectra_40_80/gSpectra_40_80_lambda_0d5
    """
    if lido_f is None:
        return None

    lido_dir = lido_spectra_dirs.get(iCent)
    obs_key = lido_observable_keys.get(iLamb)
    if not lido_dir or not obs_key:
        return None

    scale = lido_spectrum_scale(iCent)

    graph_path = f"{lido_dir}/g{lido_dir}_{obs_key}"
    obj = get_root_object(lido_f, graph_path)
    if obj is not None:
        x, y, exl, exh, eyl, eyh = read_tgraph_asymm(obj)
        lambda_factor = x
        #ncoll_factor = ncoll_values[iCent]
        ncoll_ref = ncoll_values[2]  # 40-80%
        ncoll_factor = ncoll_values[iCent] / ncoll_ref
        br_factor = 0.0395 / 2
        total_factor = lido_spectrum_total_factor(x, iCent)

        return {
            "x": x,
            "y": y * scale * total_factor,
            "exl": exl,
            "exh": exh,
            "stat_low": eyl * abs(scale) * np.abs(total_factor),
            "stat_high": eyh * abs(scale) * np.abs(total_factor),
        }

    hist_path = f"{lido_dir}/h{lido_dir}_{obs_key}"
    hobj = get_root_object(lido_f, hist_path)
    if hobj is not None:
        x, y, xerr, yerr = read_hist_values(hobj)
        lambda_factor = x
        #ncoll_factor = ncoll_values[iCent]
        ncoll_ref = ncoll_values[2]  # 40-80%
        ncoll_factor = ncoll_values[iCent] / ncoll_ref        
        br_factor = 0.0395 / 2
        total_factor = lido_spectrum_total_factor(x, iCent)
        return {
            "x": x,
            "y": y * scale * total_factor,
            "exl": xerr,
            "exh": xerr,
            "stat_low": yerr * abs(scale) * np.abs(total_factor),
            "stat_high": yerr * abs(scale) * np.abs(total_factor),
        }

    print(f"[WARN] Missing LIDO spectrum graph: {graph_path}")
    return None


def draw_lido_spectrum_model(
    ax: Any,
    lido_res: Optional[Dict[str, np.ndarray]],
    color: str,
    x_range: Optional[Tuple[float, float]] = None,
) -> bool:
    """Draw one LIDO spectrum as a colored line plus statistical band.

    If x_range is given, interpolate both band edges at its boundaries so the
    model ends exactly like the R_CP model curves.
    """
    if lido_res is None:
        return False

    x = np.asarray(lido_res.get("x", []), dtype=float)
    y = np.asarray(lido_res.get("y", []), dtype=float)
    stat_low = np.asarray(lido_res.get("stat_low", np.zeros_like(y)), dtype=float)
    stat_high = np.asarray(lido_res.get("stat_high", np.zeros_like(y)), dtype=float)

    n = min(len(x), len(y), len(stat_low), len(stat_high))
    if n == 0:
        return False

    x = x[:n]
    y = y[:n]
    stat_low = stat_low[:n]
    stat_high = stat_high[:n]

    good = (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(stat_low)
        & np.isfinite(stat_high)
        & (y > 0.0)
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
    # results_by_centrality is drawn on log-y scale. Avoid negative lower band edges.
    y_low = np.maximum(y_low, y * 1.0e-6)
    y_high = np.maximum(y_high, y_low)

    if x_range is not None:
        xmin, xmax = map(float, x_range)
        if xmax <= xmin:
            return False
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

        if len(x) >= 2 and xmin_clip > x[0] and xmin_clip < x[-1]:
            x_new.append(xmin_clip)
            y_new.append(np.interp(xmin_clip, x, y))
            ylow_new.append(np.interp(xmin_clip, x, y_low))
            yhigh_new.append(np.interp(xmin_clip, x, y_high))

        x_new.extend(x[mask_range])
        y_new.extend(y[mask_range])
        ylow_new.extend(y_low[mask_range])
        yhigh_new.extend(y_high[mask_range])

        if len(x) >= 2 and xmax_clip > x[0] and xmax_clip < x[-1]:
            x_new.append(xmax_clip)
            y_new.append(np.interp(xmax_clip, x, y))
            ylow_new.append(np.interp(xmax_clip, x, y_low))
            yhigh_new.append(np.interp(xmax_clip, x, y_high))

        x = np.asarray(x_new, dtype=float)
        y = np.asarray(y_new, dtype=float)
        y_low = np.asarray(ylow_new, dtype=float)
        y_high = np.asarray(yhigh_new, dtype=float)
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
        color=color,
        alpha=LIDO_SPECTRA_BAND_ALPHA,
        linewidth=0.0,
        zorder=1.8,
    )
    ax.plot(
        x,
        y,
        color=color,
        alpha=LIDO_SPECTRA_LINE_ALPHA,
        linewidth=LIDO_SPECTRA_LINEWIDTH,
        linestyle=LIDO_SPECTRA_LINESTYLE,
        zorder=2.0,
    )

    return True


def draw_lido_model(
    ax: Any,
    lido_res: Optional[Dict[str, np.ndarray]],
    keep_xlim=None,
    x_range: Optional[Tuple[float, float]] = None,
    color: str = "tab:blue",
    band_alpha: float = 0.25,
    line_alpha: float = 1.0,
    linestyle: str = "-",
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
        color=color,
        alpha=band_alpha,
        linewidth=0.0,
        zorder=0.6,
        label="LIDO stat. unc.",
    )

    ax.plot(
        x,
        y,
        color=color,
        alpha=line_alpha,
        linewidth=1.5,
        linestyle=linestyle,
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
            ax.set_ylim(ymin * 0.5, ymax * factor) #0.5


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


def combined_panel_rcp_types() -> list:
    """Resolve the user-facing 0/1/2 R_CP selection to ROOT-file keys."""
    selected = []
    for item in COMBINED_PANEL_RCP_SELECTION:
        resolved = None
        if isinstance(item, (int, np.integer)):
            idx = int(item)
            if 0 <= idx < len(rcp_types):
                resolved = rcp_types[idx]
        else:
            key = str(item).strip()
            if key in rcp_types:
                resolved = key

        if resolved is None:
            print(f"[WARN] Ignoring unknown combined-panel R_CP selection: {item!r}")
            continue
        if resolved not in selected:
            selected.append(resolved)
    return selected


def combined_panel_axis(axes: list, column: Optional[int], what: str) -> Optional[Any]:
    """Return a configured panel axis, allowing negative column indices."""
    idx = combined_panel_index(column, len(axes), what)
    if idx is None:
        return None
    return axes[idx]


def combined_panel_index(value: Optional[int], length: int, what: str) -> Optional[int]:
    """Resolve a possibly negative configured row/column index."""
    if value is None or length <= 0:
        return None
    try:
        idx = int(value)
    except Exception:
        print(f"[WARN] Invalid {what}: {value!r}")
        return None
    if idx < 0:
        idx += length
    if not 0 <= idx < length:
        print(f"[WARN] {what} {value!r} is outside 0..{length - 1}")
        return None
    return idx


def combined_panel_cell(
    axes_by_row: list,
    row: Optional[int],
    column: Optional[int],
    what: str,
) -> Optional[Any]:
    """Return one configured cell from a row-major list of panel axes."""
    row_idx = combined_panel_index(row, len(axes_by_row), f"{what} row")
    if row_idx is None:
        return None
    return combined_panel_axis(axes_by_row[row_idx], column, f"{what} column")


def hide_edge_major_tick_label(ax: Any, axis: str, edge: str) -> None:
    """Hide the first/last visible major tick label inside the current limits."""
    if axis == "x":
        axis_object = ax.xaxis
        limits = ax.get_xlim()
    elif axis == "y":
        axis_object = ax.yaxis
        limits = ax.get_ylim()
    else:
        raise ValueError(f"Unknown axis {axis!r}")

    lower, upper = sorted(map(float, limits))
    tolerance = 1.0e-9 * max(1.0, abs(lower), abs(upper), abs(upper - lower))
    visible_ticks = []
    for tick in axis_object.get_major_ticks():
        location = float(tick.get_loc())
        label = tick.label1
        if not math.isfinite(location):
            continue
        if location < lower - tolerance or location > upper + tolerance:
            continue
        if not label.get_visible() or not label.get_text():
            continue
        visible_ticks.append((location, label))

    if not visible_ticks:
        return
    target = min(visible_ticks, key=lambda item: item[0]) if edge == "first" else max(
        visible_ticks,
        key=lambda item: item[0],
    )
    target[1].set_visible(False)


def clean_touching_panel_tick_labels(fig: Any, axes_by_row: list) -> None:
    """Remove only the tick labels that collide where zero-gap axes meet."""
    if not axes_by_row:
        return

    # Tick text and locations are finalized lazily by Matplotlib.
    fig.canvas.draw()

    if COMBINED_PANEL_HIDE_TOUCHING_X_TICK_LABELS:
        bottom_axes = axes_by_row[-1]
        for ax in bottom_axes[:-1]:
            hide_edge_major_tick_label(ax, "x", "last")

    if COMBINED_PANEL_HIDE_TOUCHING_Y_TICK_LABELS:
        if COMBINED_PANEL_KEEP_SPECTRA_BOTTOM_Y_TICK_LABEL and len(axes_by_row) > 1:
            # Preserve the lowest spectrum label and resolve only this boundary
            # from the R_CP side by hiding the first R_CP row's upper label.
            for ax in axes_by_row[1]:
                hide_edge_major_tick_label(ax, "y", "last")
            rows_with_hidden_lower_label = axes_by_row[1:-1]
        else:
            rows_with_hidden_lower_label = axes_by_row[:-1]

        # Inside the R_CP stack, retain the upper tick of the lower panel and
        # hide the lower tick of the panel immediately above it.
        for row_axes in rows_with_hidden_lower_label:
            for ax in row_axes:
                hide_edge_major_tick_label(ax, "y", "first")


def style_combined_panel_ticks(axes_by_row: list) -> None:
    """Show tick marks on every panel side except the top side."""
    if not COMBINED_PANEL_TICKS_ON_ALL_PANELS:
        return
    for row_axes in axes_by_row:
        for ax in row_axes:
            ax.tick_params(
                axis="both",
                which="both",
                bottom=True,
                top=False,
                left=True,
                right=True,
                labeltop=False,
                labelright=False,
                direction=COMBINED_PANEL_TICK_DIRECTION,
            )


def result_x_bounds(res: Dict[str, np.ndarray | str]) -> Optional[Tuple[float, float]]:
    """Return the full visible bin range, including asymmetric x-errors."""
    x = np.asarray(res.get("x", []), dtype=float)
    exl = np.asarray(res.get("exl", np.zeros_like(x)), dtype=float)
    exh = np.asarray(res.get("exh", np.zeros_like(x)), dtype=float)
    n = min(len(x), len(exl), len(exh))
    if n == 0:
        return None
    low = x[:n] - exl[:n]
    high = x[:n] + exh[:n]
    good_low = np.isfinite(low)
    good_high = np.isfinite(high)
    if not np.any(good_low) or not np.any(good_high):
        return None
    return float(np.min(low[good_low])), float(np.max(high[good_high]))


def merge_x_bounds(bounds: list) -> Optional[Tuple[float, float]]:
    valid = [b for b in bounds if b is not None and np.isfinite(b[0]) and np.isfinite(b[1])]
    if not valid:
        return None
    xmin = min(b[0] for b in valid)
    xmax = max(b[1] for b in valid)
    if xmax <= xmin:
        pad = max(abs(xmin) * 0.05, 0.05)
        return xmin - pad, xmax + pad
    pad_fraction = max(0.0, float(COMBINED_PANEL_X_PADDING_FRACTION))
    pad = pad_fraction * (xmax - xmin)
    return xmin - pad, xmax + pad


def y_arrays_in_x_range(
    x: np.ndarray,
    y: np.ndarray,
    low: np.ndarray,
    high: np.ndarray,
    x_range: Optional[Tuple[float, float]],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return aligned y/error arrays restricted to a visible x interval."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    low = np.asarray(low, dtype=float)
    high = np.asarray(high, dtype=float)
    n = min(len(x), len(y), len(low), len(high))
    x, y, low, high = x[:n], y[:n], low[:n], high[:n]
    if x_range is None or n == 0:
        return y, low, high

    xmin, xmax = sorted(map(float, x_range))
    tolerance = 1.0e-12 * max(1.0, abs(xmin), abs(xmax), abs(xmax - xmin))
    mask = np.isfinite(x) & (x >= xmin - tolerance) & (x <= xmax + tolerance)
    return y[mask], low[mask], high[mask]


def set_linear_y_with_margin(
    ax: Any,
    y_values: Iterable[np.ndarray],
    low_errors: Iterable[np.ndarray],
    high_errors: Iterable[np.ndarray],
    include_one: bool = False,
) -> None:
    """Set an automatic linear y-range from central values and errors."""
    lows = []
    highs = []
    for y, lo, hi in zip(y_values, low_errors, high_errors):
        y = np.asarray(y, dtype=float)
        lo = np.asarray(lo, dtype=float)
        hi = np.asarray(hi, dtype=float)
        n = min(len(y), len(lo), len(hi))
        if n == 0:
            continue
        low = y[:n] - lo[:n]
        high = y[:n] + hi[:n]
        lows.extend(low[np.isfinite(low)].tolist())
        highs.extend(high[np.isfinite(high)].tolist())
    if include_one:
        lows.append(1.0)
        highs.append(1.0)
    if not lows or not highs:
        return
    ymin = min(lows)
    ymax = max(highs)
    span = max(ymax - ymin, abs(ymax) * 0.10, 0.1)
    lower = 0.0 if ymin >= 0.0 else ymin - 0.08 * span
    ax.set_ylim(lower, ymax + 0.12 * span)


def combined_spectrum_ylabel() -> str:
    return (
        r"$\frac{1}{N_{\mathrm{evt}}}\frac{1}{\mathrm{BR}}"
        r"\frac{\mathrm{d}^{2}N}{\mathrm{d}X\,\mathrm{d}\eta_{\mathrm{jet}}}$"
    )


def add_combined_spectra_legend(ax: Any, show_lido_model: bool) -> None:
    handles = [
        Line2D(
            [], [],
            marker=cent_markers.get(iCent, "o"),
            linestyle="none",
            color=cent_colors.get(iCent, "black"),
            markersize=5,
        )
        for iCent in range(len(cent_labels))
    ]
    labels = list(cent_labels)
    handles.append(
        MultiColorSysBox(
            [cent_colors[i] for i in range(len(cent_labels))],
            alpha=0.55,
            skew=0.18,
            yshift=-0.15,
        )
    )
    labels.append(r"$\mathrm{Syst.\ unc.}$")
    if show_lido_model:
        handles.append(
            MultiColorLineLegend(
                [cent_colors[i] for i in range(len(cent_labels))],
                linewidth=LIDO_SPECTRA_LINEWIDTH,
            )
        )
        labels.append(LIDO_MODEL_LABEL)

    ax.legend(
        handles,
        labels,
        frameon=False,
        loc=COMBINED_PANEL_SPECTRA_LEGEND_LOC,
        bbox_to_anchor=COMBINED_PANEL_SPECTRA_LEGEND_BBOX,
        bbox_transform=ax.transAxes,
        fontsize=COMBINED_PANEL_SPECTRA_LEGEND_FONTSIZE,
        ncol=max(1, int(COMBINED_PANEL_SPECTRA_LEGEND_NCOL)),
        columnspacing=0.8,
        handlelength=1.8,
        handletextpad=0.35,
        borderaxespad=0.0,
        handler_map={
            MultiColorSysBox: HandlerMultiColorSysBox(),
            MultiColorLineLegend: HandlerMultiColorLineLegend(),
        },
    )


def add_combined_rcp_legend(
    ax: Any,
    selected_rcp_types: list,
    show_lido_model: bool,
    include_uncertainties: bool = True,
    include_ncoll: bool = True,
) -> None:
    handles = [
        Line2D(
            [], [],
            marker=rcp_markers.get(rcp_type, "o"),
            linestyle="none",
            color=rcp_colors.get(rcp_type, "black"),
            markersize=5,
        )
        for rcp_type in selected_rcp_types
    ]
    labels = [rcp_titles[rcp_type] for rcp_type in selected_rcp_types]

    if include_uncertainties:
        sys_colors = [rcp_colors.get(rcp_type, "0.5") for rcp_type in selected_rcp_types]
        if not sys_colors:
            sys_colors = ["0.5"]
        handles.append(MultiColorSysBox(sys_colors, alpha=0.55, skew=0.18, yshift=-0.15))
        labels.append(r"$\mathrm{Syst.\ unc.}$")

        if include_ncoll:
            handles.append(ncoll_legend_patch())
            labels.append(r"$\mathrm{Rel.\ }N_{\mathrm{coll}}\mathrm{\ unc.}$")

        if show_lido_model:
            lido_color = sys_colors[0]
            handles.append(
                LidoBandLegend(
                    facecolor=lido_color,
                    alpha=0.20,
                    linecolor=lido_color,
                    linewidth=1.5,
                )
            )
            labels.append(LIDO_MODEL_LABEL)

    ax.legend(
        handles,
        labels,
        frameon=False,
        loc=COMBINED_PANEL_RCP_LEGEND_LOC,
        bbox_to_anchor=COMBINED_PANEL_RCP_LEGEND_BBOX,
        bbox_transform=ax.transAxes,
        fontsize=COMBINED_PANEL_RCP_LEGEND_FONTSIZE,
        ncol=max(1, int(COMBINED_PANEL_RCP_LEGEND_NCOL)),
        columnspacing=0.8,
        handlelength=1.8,
        handletextpad=0.35,
        borderaxespad=0.0,
        handler_map={
            MultiColorSysBox: HandlerMultiColorSysBox(),
            LidoBandLegend: HandlerLidoBandLegend(),
        },
    )


def plot_combined_spectra_rcp_panel(
    f: Any,
    output_pdf: Path,
    all_graphs_dir: Optional[Path] = None,
    lido_f: Optional[Any] = None,
    show_lido_spectra: bool = False,
    show_lido_rcp: bool = False,
) -> None:
    """Draw spectra above one separate row for every selected R_CP ratio."""
    variables = []
    for value in COMBINED_PANEL_VARIABLES:
        try:
            iLamb = int(value)
        except Exception:
            print(f"[WARN] Ignoring unknown combined-panel variable: {value!r}")
            continue
        if iLamb not in spectrum_job_dir or iLamb not in rcp_job_dir:
            print(f"[WARN] Ignoring unknown combined-panel variable: {value!r}")
            continue
        if iLamb not in variables:
            variables.append(iLamb)

    if not variables:
        print("[WARN] Combined panel skipped: COMBINED_PANEL_VARIABLES is empty.")
        return

    selected_rcp_types = combined_panel_rcp_types()
    if not selected_rcp_types:
        print("[INFO] Combined panel has no selected R_CP rows; drawing spectra only.")

    ncols = len(variables)
    n_rcp_rows = len(selected_rcp_types)
    height_ratios = [COMBINED_PANEL_SPECTRA_HEIGHT_RATIO]
    height_ratios.extend([COMBINED_PANEL_RCP_HEIGHT_RATIO] * n_rcp_rows)
    figure_height = max(
        1.0,
        COMBINED_PANEL_HEIGHT_PER_RATIO_UNIT * float(sum(height_ratios)),
    )
    fig = plt.figure(
        figsize=(max(1.0, COMBINED_PANEL_COLUMN_WIDTH * ncols), figure_height)
    )
    grid = fig.add_gridspec(
        1 + n_rcp_rows,
        ncols,
        height_ratios=height_ratios,
        wspace=COMBINED_PANEL_WSPACE,
        hspace=COMBINED_PANEL_HSPACE,
    )

    spectra_axes = []
    for icol in range(ncols):
        share_spectrum_ax = spectra_axes[0] if COMBINED_PANEL_SHARE_SPECTRA_Y and spectra_axes else None
        ax_spectrum = fig.add_subplot(grid[0, icol], sharey=share_spectrum_ax)
        spectra_axes.append(ax_spectrum)

    # Row-major layout: rcp_axes[ircp][icol].  When requested, every R_CP cell
    # shares one common y-scale; x is shared only inside the same variable
    # column because different columns represent different observables.
    rcp_axes = []
    shared_rcp_y_axis = None
    for ircp in range(n_rcp_rows):
        row_axes = []
        for icol in range(ncols):
            share_rcp_ax = shared_rcp_y_axis if COMBINED_PANEL_SHARE_RCP_Y else None
            ax_rcp = fig.add_subplot(
                grid[1 + ircp, icol],
                sharex=spectra_axes[icol],
                sharey=share_rcp_ax,
            )
            if shared_rcp_y_axis is None:
                shared_rcp_y_axis = ax_rcp
            row_axes.append(ax_rcp)
        rcp_axes.append(row_axes)

    spectra_y = [[] for _ in variables]
    spectra_low = [[] for _ in variables]
    spectra_high = [[] for _ in variables]
    rcp_y = [[[] for _ in variables] for _ in selected_rcp_types]
    rcp_low = [[[] for _ in variables] for _ in selected_rcp_types]
    rcp_high = [[[] for _ in variables] for _ in selected_rcp_types]
    lido_spectra_drawn_any = False
    lido_rcp_drawn_rows = [False for _ in selected_rcp_types]
    found_any = False

    for icol, (iLamb, ax_spectrum) in enumerate(zip(variables, spectra_axes)):
        xlabel = lambda_labels.get(iLamb, r"$x$")
        x_bounds = []
        spectrum_found = False
        rcp_found = [False for _ in selected_rcp_types]

        for iCent, (cent_label, tag) in enumerate(zip(cent_labels, cent_tags)):
            res = get_final_result(f, spectrum_job_dir[iLamb], tag, xlabel=xlabel)
            if res is None:
                print(f"[WARN] Missing combined-panel spectrum: {spectrum_job_dir[iLamb]}/g_value_stat_{tag}")
                continue
            res = apply_edge_bin_visibility(res, iLamb)
            res = apply_spectrum_global_uncertainties(res)

            x = np.asarray(res["x"], dtype=float)
            y = np.asarray(res["y"], dtype=float)
            exl = np.asarray(res["exl"], dtype=float)
            exh = np.asarray(res["exh"], dtype=float)
            stat_low = np.asarray(res["stat_low"], dtype=float)
            stat_high = np.asarray(res["stat_high"], dtype=float)
            syst_low = np.asarray(res["syst_low"], dtype=float)
            syst_high = np.asarray(res["syst_high"], dtype=float)
            color = cent_colors.get(iCent, "black")

            add_sys_boxes(
                ax_spectrum,
                x,
                y,
                exl,
                exh,
                syst_low,
                syst_high,
                alpha=0.18,
                zorder=1,
                color=color,
                skip_nonpositive=COMBINED_PANEL_LOG_SPECTRA_Y,
            )
            ax_spectrum.errorbar(
                x,
                y,
                xerr=[exl, exh],
                yerr=[stat_low, stat_high],
                fmt=cent_markers.get(iCent, "o"),
                capsize=2.5,
                elinewidth=1.0,
                capthick=1.0,
                markersize=4.5,
                linestyle="none",
                color=color,
                zorder=3,
            )
            x_bounds.append(result_x_bounds(res))
            spectra_y[icol].append(y)
            spectra_low[icol].append(np.maximum(stat_low, syst_low))
            spectra_high[icol].append(np.maximum(stat_high, syst_high))
            spectrum_found = True
            found_any = True

            if show_lido_spectra and lido_f is not None:
                lido_spectrum_x_range = COMBINED_PANEL_X_RANGES.get(iLamb)
                if lido_spectrum_x_range is None:
                    manual_range = LIDO_MANUAL_X_RANGE.get(iLamb)
                    if manual_range is not None:
                        lido_spectrum_x_range = manual_range
                    elif CLIP_LIDO_TO_DATA_RANGE and len(x) > 0:
                        lido_spectrum_x_range = (float(np.nanmin(x)), float(np.nanmax(x)))
                lido_res = get_lido_spectrum_result(lido_f, iCent, iLamb)
                if draw_lido_spectrum_model(
                    ax_spectrum,
                    lido_res,
                    color=color,
                    x_range=lido_spectrum_x_range,
                ):
                    lido_spectra_drawn_any = True
                    lx = np.asarray(lido_res.get("x", []), dtype=float)
                    ly = np.asarray(lido_res.get("y", []), dtype=float)
                    llo = np.asarray(lido_res.get("stat_low", np.zeros_like(ly)), dtype=float)
                    lhi = np.asarray(lido_res.get("stat_high", np.zeros_like(ly)), dtype=float)
                    lido_y_limit_range = None
                    if COMBINED_PANEL_SPECTRA_Y_LIMIT_VISIBLE_X_ONLY:
                        lido_y_limit_range = lido_spectrum_x_range
                        if lido_y_limit_range is None:
                            lido_y_limit_range = result_x_bounds(res)
                    ly_visible, llo_visible, lhi_visible = y_arrays_in_x_range(
                        lx,
                        ly,
                        llo,
                        lhi,
                        lido_y_limit_range,
                    )
                    if len(ly_visible) > 0:
                        spectra_y[icol].append(ly_visible)
                        spectra_low[icol].append(llo_visible)
                        spectra_high[icol].append(lhi_visible)

        for ircp, rcp_type in enumerate(selected_rcp_types):
            ax_rcp = rcp_axes[ircp][icol]
            tag = rcp_tags[rcp_type]
            res = get_final_result(f, rcp_job_dir[iLamb], tag, xlabel=xlabel)
            if res is None:
                print(f"[WARN] Missing combined-panel RCP: {rcp_job_dir[iLamb]}/g_value_stat_{tag}")
                continue
            res = apply_edge_bin_visibility(res, iLamb)
            res = apply_rcp_global_uncertainties(res)

            x = np.asarray(res["x"], dtype=float)
            y = np.asarray(res["y"], dtype=float)
            exl = np.asarray(res["exl"], dtype=float)
            exh = np.asarray(res["exh"], dtype=float)
            stat_low = np.asarray(res["stat_low"], dtype=float)
            stat_high = np.asarray(res["stat_high"], dtype=float)
            syst_low = np.asarray(res["syst_low"], dtype=float)
            syst_high = np.asarray(res["syst_high"], dtype=float)
            color = rcp_colors.get(rcp_type, "black")

            add_sys_boxes(
                ax_rcp,
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
            ax_rcp.errorbar(
                x,
                y,
                xerr=[exl, exh],
                yerr=[stat_low, stat_high],
                fmt=rcp_markers.get(rcp_type, "o"),
                capsize=2.5,
                elinewidth=1.0,
                capthick=1.0,
                markersize=4.5,
                linestyle="none",
                color=color,
                zorder=3,
            )

            add_ncoll_box(ax_rcp, rcp_type, color="teal", alpha=0.38)

            x_bounds.append(result_x_bounds(res))
            rcp_y[ircp][icol].append(y)
            rcp_low[ircp][icol].append(np.maximum(stat_low, syst_low))
            rcp_high[ircp][icol].append(np.maximum(stat_high, syst_high))
            rcp_found[ircp] = True
            found_any = True

            if show_lido_rcp and lido_f is not None:
                lido_x_range = COMBINED_PANEL_X_RANGES.get(iLamb)
                if lido_x_range is None:
                    manual_range = LIDO_MANUAL_X_RANGE.get(iLamb)
                    if manual_range is not None:
                        lido_x_range = manual_range
                    elif CLIP_LIDO_TO_DATA_RANGE and len(x) > 0:
                        lido_x_range = (float(np.nanmin(x)), float(np.nanmax(x)))
                lido_res = get_lido_model_result(lido_f, rcp_type, iLamb)
                if draw_lido_model(
                    ax_rcp,
                    lido_res,
                    x_range=lido_x_range,
                    color=color,
                    band_alpha=0.16,
                    line_alpha=0.90,
                ):
                    lido_rcp_drawn_rows[ircp] = True
                    ly = np.asarray(lido_res.get("y", []), dtype=float)
                    llo = np.asarray(lido_res.get("stat_low", np.zeros_like(ly)), dtype=float)
                    lhi = np.asarray(lido_res.get("stat_high", np.zeros_like(ly)), dtype=float)
                    rcp_y[ircp][icol].append(ly)
                    rcp_low[ircp][icol].append(llo)
                    rcp_high[ircp][icol].append(lhi)

        if not spectrum_found:
            ax_spectrum.text(0.5, 0.5, "Missing spectra", ha="center", va="center", transform=ax_spectrum.transAxes)
        for ircp, was_found in enumerate(rcp_found):
            if not was_found:
                ax_rcp = rcp_axes[ircp][icol]
                ax_rcp.text(
                    0.5,
                    0.5,
                    "Missing R$_{CP}$",
                    ha="center",
                    va="center",
                    transform=ax_rcp.transAxes,
                )

        if COMBINED_PANEL_SHOW_COLUMN_TITLES:
            ax_spectrum.set_title(
                lambda_titles.get(iLamb, f"Variable {iLamb}"),
                fontsize=COMBINED_PANEL_TITLE_FONTSIZE,
            )
        ax_spectrum.tick_params(
            axis="both",
            which="both",
            labelsize=COMBINED_PANEL_TICK_LABEL_FONTSIZE,
            labelbottom=(n_rcp_rows == 0),
        )
        if n_rcp_rows == 0:
            ax_spectrum.set_xlabel(xlabel, fontsize=COMBINED_PANEL_AXIS_LABEL_FONTSIZE)

        for ircp in range(n_rcp_rows):
            ax_rcp = rcp_axes[ircp][icol]
            is_bottom_row = ircp == n_rcp_rows - 1
            ax_rcp.tick_params(
                axis="both",
                which="both",
                labelsize=COMBINED_PANEL_TICK_LABEL_FONTSIZE,
                labelbottom=is_bottom_row,
            )
            if is_bottom_row:
                ax_rcp.set_xlabel(xlabel, fontsize=COMBINED_PANEL_AXIS_LABEL_FONTSIZE)
            ax_rcp.axhline(1.0, linestyle="--", linewidth=1.0, color="black", zorder=0)

        manual_xlim = COMBINED_PANEL_X_RANGES.get(iLamb)
        xlim = manual_xlim if manual_xlim is not None else merge_x_bounds(x_bounds)
        if xlim is not None:
            bottom_axis = rcp_axes[-1][icol] if rcp_axes else ax_spectrum
            bottom_axis.set_xlim(xlim)

    if not found_any:
        plt.close(fig)
        print("[WARN] Combined panel skipped: no requested ROOT graphs were found.")
        return

    if COMBINED_PANEL_LOG_SPECTRA_Y:
        if COMBINED_PANEL_SPECTRA_YLIM is not None:
            ymin, ymax = COMBINED_PANEL_SPECTRA_YLIM
            if ymin <= 0 or ymax <= ymin:
                raise ValueError("COMBINED_PANEL_SPECTRA_YLIM must satisfy 0 < ymin < ymax on a log scale")
            for ax in spectra_axes:
                ax.set_yscale("log")
                ax.set_ylim(ymin, ymax)
        elif COMBINED_PANEL_SHARE_SPECTRA_Y:
            set_logy_with_margin(
                spectra_axes[0],
                [arr for col in spectra_y for arr in col],
                [arr for col in spectra_low for arr in col],
                factor=COMBINED_PANEL_SPECTRA_YMAX_FACTOR,
            )
        else:
            for ax, ys, lows in zip(spectra_axes, spectra_y, spectra_low):
                set_logy_with_margin(
                    ax,
                    ys,
                    lows,
                    factor=COMBINED_PANEL_SPECTRA_YMAX_FACTOR,
                )
    else:
        if COMBINED_PANEL_SPECTRA_YLIM is not None:
            for ax in spectra_axes:
                ax.set_ylim(COMBINED_PANEL_SPECTRA_YLIM)
        elif COMBINED_PANEL_SHARE_SPECTRA_Y:
            set_linear_y_with_margin(
                spectra_axes[0],
                [arr for col in spectra_y for arr in col],
                [arr for col in spectra_low for arr in col],
                [arr for col in spectra_high for arr in col],
            )
        else:
            for ax, ys, lows, highs in zip(spectra_axes, spectra_y, spectra_low, spectra_high):
                set_linear_y_with_margin(ax, ys, lows, highs)

    flat_rcp_axes = [ax for row_axes in rcp_axes for ax in row_axes]
    if flat_rcp_axes:
        if COMBINED_PANEL_RCP_YLIM is not None:
            for ax in flat_rcp_axes:
                ax.set_ylim(COMBINED_PANEL_RCP_YLIM)
        elif COMBINED_PANEL_SHARE_RCP_Y:
            set_linear_y_with_margin(
                flat_rcp_axes[0],
                [arr for row in rcp_y for col in row for arr in col],
                [arr for row in rcp_low for col in row for arr in col],
                [arr for row in rcp_high for col in row for arr in col],
                include_one=True,
            )
        else:
            for ircp, row_axes in enumerate(rcp_axes):
                for icol, ax in enumerate(row_axes):
                    set_linear_y_with_margin(
                        ax,
                        rcp_y[ircp][icol],
                        rcp_low[ircp][icol],
                        rcp_high[ircp][icol],
                        include_one=True,
                    )

    if COMBINED_PANEL_SHARE_SPECTRA_Y:
        spectra_axes[0].set_ylabel(combined_spectrum_ylabel(), fontsize=COMBINED_PANEL_AXIS_LABEL_FONTSIZE)
        for ax in spectra_axes[1:]:
            ax.tick_params(axis="y", labelleft=False)
    else:
        for ax, iLamb in zip(spectra_axes, variables):
            ax.set_ylabel(make_ylabel_from_xlabel(lambda_labels.get(iLamb, r"$x$")), fontsize=COMBINED_PANEL_AXIS_LABEL_FONTSIZE)

    for ircp, (rcp_type, row_axes) in enumerate(zip(selected_rcp_types, rcp_axes)):
        row_axes[0].set_ylabel(r"$R_{\mathrm{CP}}$", fontsize=COMBINED_PANEL_AXIS_LABEL_FONTSIZE)
        if COMBINED_PANEL_SHARE_RCP_Y:
            for ax in row_axes[1:]:
                ax.tick_params(axis="y", labelleft=False)

    spectrum_legend_ax = combined_panel_axis(
        spectra_axes,
        COMBINED_PANEL_SPECTRA_LEGEND_COLUMN,
        "combined-panel spectra legend",
    )
    if spectrum_legend_ax is not None:
        add_combined_spectra_legend(spectrum_legend_ax, show_lido_model=lido_spectra_drawn_any)

    full_rcp_legend_row = combined_panel_index(
        COMBINED_PANEL_RCP_LEGEND_ROW,
        n_rcp_rows,
        "combined-panel R_CP legend row",
    )
    if COMBINED_PANEL_RCP_LEGEND_COLUMN is not None:
        for ircp, rcp_type in enumerate(selected_rcp_types):
            rcp_legend_ax = combined_panel_axis(
                rcp_axes[ircp],
                COMBINED_PANEL_RCP_LEGEND_COLUMN,
                f"combined-panel R_CP legend column for row {ircp}",
            )
            if rcp_legend_ax is None:
                continue
            add_combined_rcp_legend(
                rcp_legend_ax,
                [rcp_type],
                show_lido_model=lido_rcp_drawn_rows[ircp],
                include_uncertainties=True,
                include_ncoll=(ircp == full_rcp_legend_row),
            )

    header_ax = combined_panel_axis(spectra_axes, COMBINED_PANEL_HEADER_COLUMN, "combined-panel header")
    if header_ax is not None:
        header_ax.text(
            COMBINED_PANEL_HEADER_X,
            COMBINED_PANEL_HEADER_Y,
            "\n".join(legend_header_labels()),
            transform=header_ax.transAxes,
            ha=COMBINED_PANEL_HEADER_HA,
            va=COMBINED_PANEL_HEADER_VA,
            fontsize=COMBINED_PANEL_HEADER_FONTSIZE,
            linespacing=1.20,
        )

    if rcp_axes:
        preliminary_ax = combined_panel_cell(
            rcp_axes,
            COMBINED_PANEL_PRELIMINARY_RCP_ROW,
            COMBINED_PANEL_PRELIMINARY_COLUMN,
            "combined-panel preliminary label",
        )
    else:
        preliminary_ax = combined_panel_axis(
            spectra_axes,
            COMBINED_PANEL_PRELIMINARY_COLUMN,
            "combined-panel preliminary label column",
        )
    if preliminary_ax is not None:
        add_preliminary_label(
            preliminary_ax,
            x=COMBINED_PANEL_PRELIMINARY_X,
            y=COMBINED_PANEL_PRELIMINARY_Y,
            ha=COMBINED_PANEL_PRELIMINARY_HA,
            va=COMBINED_PANEL_PRELIMINARY_VA,
        )

    fig.subplots_adjust(
        left=0.075,
        right=0.995,
        bottom=0.075,
        top=0.995,
        wspace=COMBINED_PANEL_WSPACE,
        hspace=COMBINED_PANEL_HSPACE,
    )
    all_panel_axes = [spectra_axes] + rcp_axes
    style_combined_panel_ticks(all_panel_axes)
    clean_touching_panel_tick_labels(fig, all_panel_axes)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_pdf, bbox_inches="tight")
    output_png = output_pdf.with_suffix(".png")
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    save_individual_figure(fig, all_graphs_dir, "spectra_rcp_combined_panel")
    plt.close(fig)
    print(f"[OK] Saved: {output_pdf}")
    print(f"[OK] Saved: {output_png}")



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


def plot_spectra_overlay_by_centrality(
    f: Any,
    output_pdf: Path,
    all_graphs_dir: Optional[Path] = None,
    lido_f: Optional[Any] = None,
    show_lido_model: bool = False,
) -> None:
    with PdfPages(output_pdf) as pdf:
        for iCent, (cent_label, tag) in enumerate(zip(cent_labels, cent_tags)):
            fig, ax = plt.subplots(figsize=(8.5, 6.0))
            found_any = False
            xlabel = r"$\lambda_{\alpha}^{1}$"
            ylabel = "Counts"
            y_for_limits = []
            lo_for_limits = []
            plotted_lambdas = []
            lido_drawn_any = False

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

                if show_lido_model and lido_f is not None:
                    lido_res = get_lido_spectrum_result(lido_f, iCent, iLamb)
                    if lido_res is not None:
                        lido_res = apply_edge_bin_visibility(lido_res, iLamb)
                        if draw_lido_spectrum_model(ax, lido_res, color=color):
                            lido_drawn_any = True
                            y_lido = np.asarray(lido_res.get("y", []), dtype=float)
                            stat_lido = np.asarray(lido_res.get("stat_low", np.zeros_like(y_lido)), dtype=float)
                            y_for_limits.append(y_lido)
                            lo_for_limits.append(stat_lido)

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
            set_logy_with_margin(ax, y_for_limits, lo_for_limits, factor=10.0)

            handles, labels = ax.get_legend_handles_labels()
            add_spectra_overlay_legend(ax, handles, labels, plotted_lambdas, show_lido_model=lido_drawn_any)
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

                marker = RCP_OVERLAY_MARKERS.get(iLamb, "o")

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
                    fmt=marker,
                    capsize=3,
                    elinewidth=1.6,
                    capthick=1.6,
                    markersize=7,
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
        help="Draw LIDO/model comparison in R_CP, spectra overlays, and the combined spectra/R_CP panel.",
    )
    p.add_argument(
        "--lido-root-file",
        default=DEFAULT_LIDO_ROOT_FILE,
        help="ROOT file with LIDO/model R_CP graphs and optional Spectra_* directories",
    )
    p.add_argument(
        "--no-lido-spectra",
        action="store_true",
        help="With --show-lido, do not draw LIDO spectra in spectra overlays or the combined panel.",
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
    output_pdf_combined_panel = out_dir / f"results_spectra_RCP_panel{VERSION_TAG}.pdf"

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
                plot_spectra_overlay_by_centrality(
                    f,
                    output_pdf_overlay,
                    all_graphs_dir=all_graphs_dir,
                    lido_f=lido_f,
                    show_lido_model=(not args.no_lido_spectra),
                )
                plot_rcp_pages(f, output_pdf_rcp, lido_f=lido_f, show_lido_model=True)
                plot_rcp_individual_graphs(f, all_graphs_dir, lido_f=lido_f, show_lido_model=True)
                plot_rcp_overlay(f, output_pdf_rcp_overlay, all_graphs_dir=all_graphs_dir)
                plot_combined_spectra_rcp_panel(
                    f,
                    output_pdf_combined_panel,
                    all_graphs_dir=all_graphs_dir,
                    lido_f=lido_f,
                    show_lido_spectra=(
                        COMBINED_PANEL_SHOW_LIDO_SPECTRA and not args.no_lido_spectra
                    ),
                    show_lido_rcp=COMBINED_PANEL_SHOW_LIDO_RCP,
                )
        else:
            plot_spectra_by_observable(f, output_pdf, all_graphs_dir=all_graphs_dir)
            plot_spectra_overlay_by_centrality(f, output_pdf_overlay, all_graphs_dir=all_graphs_dir)
            plot_rcp_pages(f, output_pdf_rcp)
            plot_rcp_individual_graphs(f, all_graphs_dir)
            plot_rcp_overlay(f, output_pdf_rcp_overlay, all_graphs_dir=all_graphs_dir)
            plot_combined_spectra_rcp_panel(
                f,
                output_pdf_combined_panel,
                all_graphs_dir=all_graphs_dir,
            )

    if all_graphs_dir is not None:
        print(f"[OK] Individual graphs directory: {all_graphs_dir}")
    print("[OK] Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
