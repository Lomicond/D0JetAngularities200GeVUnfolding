#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import uproot
except Exception:
    uproot = None

try:
    import ROOT
except Exception:
    ROOT = None


METRIC_COLS = [
    "worst_bin_pct",
    "mean_abs_drift_pct",
    "rms_drift_pct",
    "weighted_drift_pct",
    "unfolded_to_mc_pct",
]

METRIC_LABELS = {
    "worst_bin_pct": "Worst bin [%]",
    "mean_abs_drift_pct": "Mean abs drift [%]",
    "rms_drift_pct": "RMS drift [%]",
    "weighted_drift_pct": "Weighted drift [%]",
    "unfolded_to_mc_pct": "Unfolded/MC TVD [%]",
}

SUMMARY_CANONICAL_COLS = [
    "run_id",
    "reco_pt_min",
    "reco_pt_max",
    "reco_n_bins",
    "true_pt_min",
    "true_pt_max",
    "true_n_bins",
    "min_width",
    "max_width",
    "step",
    "trend",
    "reco_edges",
    "true_edges",
]

RESULT_COLS = [
    "run_id",
    "cent",
    "dim",
    "observable",
    "axis",
    "reco_pt_min",
    "reco_pt_max",
    "reco_n_bins",
    "true_pt_min",
    "true_pt_max",
    "true_n_bins",
    "worst_bin_pct",
    "mean_abs_drift_pct",
    "rms_drift_pct",
    "weighted_drift_pct",
    "unfolded_to_mc_pct",
]

FILTER_COLS = [
    "cent",
    "dim",
    "observable",
    "axis",
    "reco_pt_min",
    "reco_pt_max",
    "reco_n_bins",
    "true_pt_min",
    "true_pt_max",
    "true_n_bins",
    "trend",
]

BASE_FILTER_COLS = ["cent", "dim", "observable", "axis"]

FILTER_LABELS = {
    "cent": "Centrality",
    "dim": "Dim",
    "observable": "Observable",
    "axis": "Axis",
    "reco_pt_min": "Reco min",
    "reco_pt_max": "Reco max",
    "reco_n_bins": "Reco N bins",
    "true_pt_min": "True min",
    "true_pt_max": "True max",
    "true_n_bins": "True N bins",
    "trend": "Trend",
}

DISPLAY_LABELS = {
    "run_id": "Run",
    "cent": "Cent",
    "dim": "Dim",
    "observable": "Observable",
    "axis": "Axis",
    "reco_pt_min": "Reco min",
    "reco_pt_max": "Reco max",
    "reco_n_bins": "Reco N",
    "true_pt_min": "True min",
    "true_pt_max": "True max",
    "true_n_bins": "True N",
    "worst_bin_pct": "Worst [%]",
    "mean_abs_drift_pct": "Mean abs [%]",
    "rms_drift_pct": "RMS [%]",
    "weighted_drift_pct": "Weighted [%]",
    "unfolded_to_mc_pct": "TVD [%]",
    "count": "Count",
}

CROSS_DISPLAY_MODES = [
    "All rows of matching (run, cent)",
    "Only rows matching Filter B",
]

RESULT_HIST_KIND_MAP = {
    "PT (1D)": "d0pt",
    "PT, z": "d0ptLambda0",
    "PT, #lambda^{1}_{1}": "d0ptLambda1",
    "PT, #lambda^{1}_{1.5}": "d0ptLambda2",
    "PT, #lambda^{1}_{2}": "d0ptLambda3",
    "PT, #lambda^{1}_{3}": "d0ptLambda4",
    "PT, #lambda^{1}_{0.5}": "d0ptLambda5",
    "PT, P_{T}^{D}": "d0ptLambda6",
    "z": "Lambda0",
    "#lambda^{1}_{1}": "Lambda1",
    "#lambda^{1}_{1.5}": "Lambda2",
    "#lambda^{1}_{2}": "Lambda3",
    "#lambda^{1}_{3}": "Lambda4",
    "#lambda^{1}_{0.5}": "Lambda5",
    "P_{T}^{D}": "Lambda6",
    "R_CP (5-20), z": "RCP_5_20_Lambda0",
    "R_CP (5-20), #lambda^{1}_{1}": "RCP_5_20_Lambda1",
    "R_CP (5-20), #lambda^{1}_{1.5}": "RCP_5_20_Lambda2",
    "R_CP (5-20), #lambda^{1}_{2}": "RCP_5_20_Lambda3",
    "R_CP (5-20), #lambda^{1}_{3}": "RCP_5_20_Lambda4",
    "R_CP (5-20), #lambda^{1}_{0.5}": "RCP_5_20_Lambda5",
    "R_CP (5-20), P_{T}^{D}": "RCP_5_20_Lambda6",
}

RESULT_HIST_KIND_LABELS = list(RESULT_HIST_KIND_MAP.keys())

HEATMAP_PARAM_CHOICES = [
    "reco_pt_min",
    "reco_pt_max",
    "reco_n_bins",
    "true_pt_min",
    "true_pt_max",
    "true_n_bins",
]


def maybe_float(x):
    try:
        return float(x)
    except Exception:
        return None


def format_value(x):
    if pd.isna(x):
        return ""
    if isinstance(x, (float, np.floating)):
        if math.isfinite(float(x)):
            # Use significant digits so tiny spectrum values are not displayed as 0.
            return f"{float(x):.8g}"
        return str(x)
    return str(x)


def load_stability(path: str) -> pd.DataFrame:
    base_cols = [
        "run_id",
        "cent",
        "dim",
        "observable",
        "axis",
        "worst_bin_pct",
        "mean_abs_drift_pct",
        "rms_drift_pct",
        "weighted_drift_pct",
    ]
    optional_cols = ["unfolded_to_mc_pct"]
    all_cols = base_cols + optional_cols

    raw = pd.read_csv(path, sep="\t", header=None)
    ncols = raw.shape[1]

    if ncols < len(base_cols):
        raise ValueError(
            f"stability.tsv has only {ncols} columns, but I expect at least {len(base_cols)}."
        )
    if ncols > len(all_cols):
        raise ValueError(
            f"stability.tsv has {ncols} columns, but this script supports at most {len(all_cols)}."
        )

    raw.columns = all_cols[:ncols]
    df = raw.copy()

    for col in all_cols:
        if col not in df.columns:
            df[col] = pd.NA

    df["run_id"] = df["run_id"].astype(str)
    df["cent"] = pd.to_numeric(df["cent"], errors="coerce").astype("Int64")
    df["dim"] = df["dim"].fillna("").astype(str)
    df["observable"] = df["observable"].fillna("").astype(str)
    df["axis"] = df["axis"].fillna("").astype(str)

    for c in METRIC_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["cent"]).copy()
    df["cent"] = df["cent"].astype(int)
    return df


def load_summary(path: str) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", header=0)
    cols = list(raw.columns)

    if "run_id" not in cols:
        raise ValueError(
            "summary.tsv must contain the 'run_id' column.\n"
            f"Found columns:\n{cols}"
        )

    # Convert older supported formats to unified names
    if cols == SUMMARY_CANONICAL_COLS:
        df = raw.copy()
    elif cols == ["run_id", "pt_min", "pt_max", "n_bins", "min_width", "step", "trend", "edges"]:
        df = pd.DataFrame({
            "run_id": raw["run_id"],
            "reco_pt_min": raw["pt_min"],
            "reco_pt_max": raw["pt_max"],
            "reco_n_bins": raw["n_bins"],
            "true_pt_min": pd.NA,
            "true_pt_max": pd.NA,
            "true_n_bins": pd.NA,
            "min_width": raw["min_width"],
            "max_width": pd.NA,
            "step": raw["step"],
            "trend": raw["trend"],
            "reco_edges": raw["edges"],
            "true_edges": pd.NA,
        })
    elif cols == ["run_id", "pt_min", "pt_max", "n_bins", "min_width", "max_width", "step", "trend", "edges"]:
        df = pd.DataFrame({
            "run_id": raw["run_id"],
            "reco_pt_min": raw["pt_min"],
            "reco_pt_max": raw["pt_max"],
            "reco_n_bins": raw["n_bins"],
            "true_pt_min": pd.NA,
            "true_pt_max": pd.NA,
            "true_n_bins": pd.NA,
            "min_width": raw["min_width"],
            "max_width": raw["max_width"],
            "step": raw["step"],
            "trend": raw["trend"],
            "reco_edges": raw["edges"],
            "true_edges": pd.NA,
        })
    else:
        # General format: keep all columns exactly as they are in summary
        df = raw.copy()

    df["run_id"] = df["run_id"].astype(str)

    # Keep 'edges' columns as text; convert the others to numbers only when
    # they are truly fully numeric
    for c in df.columns:
        if c == "run_id":
            continue
        if "edges" in c.lower():
            df[c] = df[c].fillna("").astype(str)
            continue

        converted = pd.to_numeric(df[c], errors="coerce")
        non_empty = df[c].notna().sum()
        numeric_ok = converted.notna().sum()

        if non_empty > 0 and numeric_ok == non_empty:
            df[c] = converted
        else:
            df[c] = df[c].fillna("").astype(str)

    return df.copy()


def merge_data(summary_path: str, stability_path: str) -> pd.DataFrame:
    sdf = load_summary(summary_path)
    tdf = load_stability(stability_path)
    return tdf.merge(sdf, on="run_id", how="left", validate="many_to_one")



def run_id_to_int(run_id: str) -> int:
    m = re.search(r"(\d+)$", str(run_id))
    return int(m.group(1)) if m else -1


def normalize_run_id(run_text: str) -> str:
    txt = str(run_text).strip()
    if not txt:
        return ""
    if txt.startswith("r") and txt[1:].isdigit():
        return txt
    if txt.isdigit():
        return f"r{int(txt):06d}"
    m = re.search(r"(\d+)$", txt)
    if m:
        return f"r{int(m.group(1)):06d}"
    return txt


def build_result_hist_name(kind_label: str, cent: int, iteration_display: int, method: str) -> str:
    """
    The GUI uses human-friendly iteration numbering: entering 4 means the 4th
    iteration. Histogram names in the ROOT files are zero-based, so this is
    loaded as it3.
    """
    if kind_label not in RESULT_HIST_KIND_MAP:
        raise KeyError(f"Unknown ROOT result: {kind_label}")
    if iteration_display < 1:
        raise ValueError("Iteration is entered as 1-based: for example, 4 loads the ROOT histogram *_it3_*" )
    kind = RESULT_HIST_KIND_MAP[kind_label]
    iteration_root = iteration_display - 1
    return f"{kind}_{cent}_it{iteration_root}_{method}"


def centers_from_edges(edges: np.ndarray) -> np.ndarray:
    return 0.5 * (edges[:-1] + edges[1:])


def step_plot(ax, edges, values, **kwargs):
    ax.stairs(values, edges, **kwargs)


def fill_between_bins(ax, edges, y_low, y_high, **kwargs):
    """Draw a per-bin band using exact bin edges instead of only bin centers."""
    x = np.repeat(np.asarray(edges, dtype=float), 2)[1:-1]
    low = np.repeat(np.asarray(y_low, dtype=float), 2)
    high = np.repeat(np.asarray(y_high, dtype=float), 2)
    ax.fill_between(x, low, high, **kwargs)


def rebin_density_to_edges(old_edges: np.ndarray, old_values: np.ndarray, new_edges: np.ndarray) -> np.ndarray:
    out = np.full(len(new_edges) - 1, np.nan, dtype=float)
    for j in range(len(new_edges) - 1):
        lo = float(new_edges[j])
        hi = float(new_edges[j + 1])
        width = hi - lo
        if width <= 0:
            continue
        acc = 0.0
        overlap_total = 0.0
        for i in range(len(old_values)):
            olo = float(old_edges[i])
            ohi = float(old_edges[i + 1])
            overlap = min(hi, ohi) - max(lo, olo)
            if overlap > 0:
                acc += float(old_values[i]) * overlap
                overlap_total += overlap
        if overlap_total > 0:
            out[j] = acc / width
    return out


def find_root_object_name(file_obj, hist_name: str):
    if hist_name in file_obj:
        return hist_name
    prefix = hist_name + ";"
    for key in file_obj.keys():
        if str(key).startswith(prefix):
            return key
    return None


def read_hist_arrays(root_path: str, hist_name: str):
    if uproot is not None:
        with uproot.open(root_path) as f:
            key = find_root_object_name(f, hist_name)
            if key is None:
                raise KeyError(f"Histogram {hist_name} was not found in {root_path}")
            h = f[key]
            values = np.asarray(h.values(flow=False), dtype=float)
            edges = np.asarray(h.axis().edges(), dtype=float)
            return edges, values

    if ROOT is not None:
        f = ROOT.TFile.Open(root_path, "READ")
        if not f or f.IsZombie():
            raise OSError(f"Cannot open ROOT file: {root_path}")
        try:
            obj = f.Get(hist_name)
            if obj is None:
                key = None
                for k in f.GetListOfKeys():
                    name = k.GetName()
                    if name == hist_name:
                        key = name
                        break
                if key is None:
                    raise KeyError(f"Histogram {hist_name} was not found in {root_path}")
                obj = f.Get(key)
            nb = obj.GetNbinsX()
            edges = np.array([obj.GetBinLowEdge(i + 1) for i in range(nb)] + [obj.GetBinLowEdge(nb + 1)], dtype=float)
            values = np.array([obj.GetBinContent(i + 1) for i in range(nb)], dtype=float)
            return edges, values
        finally:
            f.Close()

    raise ImportError("Neither uproot nor ROOT (PyROOT) is available.")


class StabilityBrowser(tk.Tk):
    def __init__(self, summary_path="", stability_path=""):
        super().__init__()
        self.title("Stability browser v12")
        self.geometry("1800x1080")

        self.summary_path_var = tk.StringVar(value=summary_path)
        self.stability_path_var = tk.StringVar(value=stability_path)

        self.df_all = pd.DataFrame()
        self.df_filtered = pd.DataFrame()

        self.filter_cols = BASE_FILTER_COLS.copy()
        self.result_cols = ["run_id"] + self.filter_cols + METRIC_COLS
        self.summary_cols = []
        self.summary_filter_cols = []
        self.summary_edge_cols = []
        self.heatmap_param_choices = []

        self.filter_vars = {}
        self.filter_boxes = {}
        self.metric_min_vars = {m: tk.StringVar(value="") for m in METRIC_COLS}
        self.metric_max_vars = {m: tk.StringVar(value="") for m in METRIC_COLS}

        self.cross_vars = {}
        for tag in ("A", "B"):
            self.cross_vars[tag] = {
                "observable": tk.StringVar(value="(any)"),
                "axis": tk.StringVar(value="(any)"),
                "metric": tk.StringVar(value="(none)"),
                "min": tk.StringVar(value=""),
                "max": tk.StringVar(value=""),
            }

        self.cross_display_mode_var = tk.StringVar(value=CROSS_DISPLAY_MODES[0])

        self.group_metric_var = tk.StringVar(value=METRIC_COLS[0])
        self.group_agg_var = tk.StringVar(value="mean")

        self.results_root_pattern_var = tk.StringVar(value="./pTCheck/OutputSpectra{run_id}.root")
        self.results_kind_var = tk.StringVar(value="PT (1D)")
        self.results_cent_var = tk.StringVar(value="0")
        self.results_iter_var = tk.StringVar(value="4")
        self.results_method_var = tk.StringVar(value="ICS")
        self.results_ref_run_var = tk.StringVar(value="")
        self.results_heatmap_x_var = tk.StringVar(value="reco_pt_min")
        self.results_heatmap_y_var = tk.StringVar(value="reco_pt_max")
        self.results_heatmap_bin_var = tk.StringVar(value="1")
        self.results_status_var = tk.StringVar(value="")

        self._build_ui()

        if summary_path and stability_path:
            self.load_files()

    def _build_ui(self):
        self.columnconfigure(0, weight=1)
        self.rowconfigure(2, weight=1)
        self.rowconfigure(3, weight=1)

        top = ttk.LabelFrame(self, text="Input files")
        top.grid(row=0, column=0, sticky="ew", padx=8, pady=8)
        top.columnconfigure(1, weight=1)
        top.columnconfigure(4, weight=1)

        ttk.Label(top, text="summary.tsv").grid(row=0, column=0, sticky="w", padx=4, pady=4)
        ttk.Entry(top, textvariable=self.summary_path_var).grid(row=0, column=1, sticky="ew", padx=4, pady=4)
        ttk.Button(top, text="Browse", command=self.browse_summary).grid(row=0, column=2, padx=4, pady=4)

        ttk.Label(top, text="stability.tsv").grid(row=0, column=3, sticky="w", padx=12, pady=4)
        ttk.Entry(top, textvariable=self.stability_path_var).grid(row=0, column=4, sticky="ew", padx=4, pady=4)
        ttk.Button(top, text="Browse", command=self.browse_stability).grid(row=0, column=5, padx=4, pady=4)

        ttk.Button(top, text="Load", command=self.load_files).grid(row=0, column=6, padx=10, pady=4)

        filters = ttk.LabelFrame(self, text="Filters")
        filters.grid(row=1, column=0, sticky="ew", padx=8, pady=4)
        for i in range(6):
            filters.columnconfigure(i, weight=1)

        self.cat_frame = ttk.Frame(filters)
        self.cat_frame.grid(row=0, column=0, columnspan=6, sticky="ew", padx=4, pady=4)

        metric_frame = ttk.LabelFrame(filters, text="Metric ranges by row")
        metric_frame.grid(row=1, column=0, columnspan=6, sticky="ew", padx=4, pady=6)

        for i, m in enumerate(METRIC_COLS):
            metric_frame.columnconfigure(i, weight=1)
            f = ttk.Frame(metric_frame)
            f.grid(row=0, column=i, sticky="ew", padx=4, pady=4)
            ttk.Label(f, text=METRIC_LABELS[m]).grid(row=0, column=0, columnspan=2, sticky="w")
            ttk.Label(f, text="min").grid(row=1, column=0, sticky="w")
            ttk.Entry(f, textvariable=self.metric_min_vars[m], width=8).grid(row=1, column=1, sticky="ew", padx=2)
            ttk.Label(f, text="max").grid(row=2, column=0, sticky="w")
            ttk.Entry(f, textvariable=self.metric_max_vars[m], width=8).grid(row=2, column=1, sticky="ew", padx=2)

        cross_frame = ttk.LabelFrame(filters, text="Cross-filters on the same (run_id, cent)")
        cross_frame.grid(row=2, column=0, columnspan=6, sticky="ew", padx=4, pady=6)
        for i in range(2):
            cross_frame.columnconfigure(i, weight=1)

        self.cross_boxes = {}
        for idx, tag in enumerate(("A", "B")):
            lf = ttk.LabelFrame(cross_frame, text=f"Filter {tag}")
            lf.grid(row=0, column=idx, sticky="nsew", padx=4, pady=4)
            lf.columnconfigure(1, weight=1)

            ttk.Label(lf, text="Observable").grid(row=0, column=0, sticky="w", padx=4, pady=2)
            cb_obs = ttk.Combobox(lf, textvariable=self.cross_vars[tag]["observable"], state="readonly")
            cb_obs["values"] = ["(any)"]
            cb_obs.grid(row=0, column=1, sticky="ew", padx=4, pady=2)

            ttk.Label(lf, text="Axis").grid(row=1, column=0, sticky="w", padx=4, pady=2)
            cb_axis = ttk.Combobox(lf, textvariable=self.cross_vars[tag]["axis"], state="readonly")
            cb_axis["values"] = ["(any)"]
            cb_axis.grid(row=1, column=1, sticky="ew", padx=4, pady=2)

            ttk.Label(lf, text="Metric").grid(row=2, column=0, sticky="w", padx=4, pady=2)
            cb_metric = ttk.Combobox(lf, textvariable=self.cross_vars[tag]["metric"], state="readonly")
            cb_metric["values"] = ["(none)"] + METRIC_COLS
            cb_metric.grid(row=2, column=1, sticky="ew", padx=4, pady=2)

            ttk.Label(lf, text="min").grid(row=3, column=0, sticky="w", padx=4, pady=2)
            ttk.Entry(lf, textvariable=self.cross_vars[tag]["min"]).grid(row=3, column=1, sticky="ew", padx=4, pady=2)

            ttk.Label(lf, text="max").grid(row=4, column=0, sticky="w", padx=4, pady=2)
            ttk.Entry(lf, textvariable=self.cross_vars[tag]["max"]).grid(row=4, column=1, sticky="ew", padx=4, pady=2)

            self.cross_boxes[tag] = {
                "observable": cb_obs,
                "axis": cb_axis,
                "metric": cb_metric,
            }

        mode_frame = ttk.Frame(cross_frame)
        mode_frame.grid(row=1, column=0, columnspan=2, sticky="ew", padx=4, pady=4)
        ttk.Label(mode_frame, text="After cross-filter, display").pack(side="left", padx=4)
        ttk.Combobox(
            mode_frame,
            textvariable=self.cross_display_mode_var,
            values=CROSS_DISPLAY_MODES,
            state="readonly",
            width=34,
        ).pack(side="left", padx=4)

        btns = ttk.Frame(filters)
        btns.grid(row=3, column=0, columnspan=6, sticky="ew", padx=4, pady=4)
        ttk.Button(btns, text="Apply filters", command=self.apply_filters).pack(side="left", padx=4)
        ttk.Button(btns, text="Reset filters", command=self.reset_filters).pack(side="left", padx=4)
        ttk.Button(btns, text="Export filtered CSV", command=self.export_filtered_csv).pack(side="left", padx=4)

        self.status_var = tk.StringVar(value="No data loaded.")
        ttk.Label(btns, textvariable=self.status_var).pack(side="right", padx=4)

        result_frame = ttk.LabelFrame(self, text="Matching rows")
        result_frame.grid(row=2, column=0, sticky="nsew", padx=8, pady=4)
        result_frame.columnconfigure(0, weight=1)
        result_frame.rowconfigure(0, weight=1)

        self.result_tree = ttk.Treeview(
            result_frame,
            columns=(),
            show="headings",
            selectmode="browse",
        )
        self.result_tree.grid(row=0, column=0, sticky="nsew")
        self.result_tree.bind("<<TreeviewSelect>>", self.on_result_select)

        ysb = ttk.Scrollbar(result_frame, orient="vertical", command=self.result_tree.yview)
        xsb = ttk.Scrollbar(result_frame, orient="horizontal", command=self.result_tree.xview)
        self.result_tree.configure(yscrollcommand=ysb.set, xscrollcommand=xsb.set)
        ysb.grid(row=0, column=1, sticky="ns")
        xsb.grid(row=1, column=0, sticky="ew")

        bottom = ttk.Notebook(self)
        bottom.grid(row=3, column=0, sticky="nsew", padx=8, pady=4)

        detail_tab = ttk.Frame(bottom)
        detail_tab.columnconfigure(0, weight=1)
        detail_tab.rowconfigure(0, weight=1)
        self.detail_text = tk.Text(detail_tab, wrap="word", height=12)
        self.detail_text.grid(row=0, column=0, sticky="nsew")
        detail_scroll = ttk.Scrollbar(detail_tab, orient="vertical", command=self.detail_text.yview)
        self.detail_text.configure(yscrollcommand=detail_scroll.set)
        detail_scroll.grid(row=0, column=1, sticky="ns")
        bottom.add(detail_tab, text="Selected row detail")

        agg_tab = ttk.Frame(bottom)
        agg_tab.columnconfigure(1, weight=1)
        agg_tab.rowconfigure(1, weight=1)
        bottom.add(agg_tab, text="Aggregation / impact mapping")

        left_agg = ttk.LabelFrame(agg_tab, text="Aggregation settings")
        left_agg.grid(row=0, column=0, sticky="nsw", padx=4, pady=4)

        ttk.Label(left_agg, text="Group by (Ctrl/Shift for multiple)").grid(row=0, column=0, sticky="w", padx=4, pady=4)
        self.groupby_listbox = tk.Listbox(left_agg, selectmode="extended", height=10, exportselection=False)
        self.groupby_listbox.grid(row=1, column=0, sticky="nsew", padx=4, pady=4)

        ttk.Label(left_agg, text="Metric").grid(row=2, column=0, sticky="w", padx=4, pady=2)
        ttk.Combobox(left_agg, textvariable=self.group_metric_var, values=METRIC_COLS, state="readonly").grid(row=3, column=0, sticky="ew", padx=4, pady=2)

        ttk.Label(left_agg, text="Aggregation").grid(row=4, column=0, sticky="w", padx=4, pady=2)
        ttk.Combobox(left_agg, textvariable=self.group_agg_var, values=["mean", "median", "min", "max"], state="readonly").grid(row=5, column=0, sticky="ew", padx=4, pady=2)

        ttk.Button(left_agg, text="Compute aggregation", command=self.compute_group_summary).grid(row=6, column=0, sticky="ew", padx=4, pady=8)
        ttk.Button(left_agg, text="Export aggregation CSV", command=self.export_group_csv).grid(row=7, column=0, sticky="ew", padx=4, pady=2)

        right_agg = ttk.LabelFrame(agg_tab, text="Aggregation result")
        right_agg.grid(row=0, column=1, rowspan=2, sticky="nsew", padx=4, pady=4)
        right_agg.columnconfigure(0, weight=1)
        right_agg.rowconfigure(0, weight=1)

        self.group_tree = ttk.Treeview(right_agg, show="headings")
        self.group_tree.grid(row=0, column=0, sticky="nsew")

        gy = ttk.Scrollbar(right_agg, orient="vertical", command=self.group_tree.yview)
        gx = ttk.Scrollbar(right_agg, orient="horizontal", command=self.group_tree.xview)
        self.group_tree.configure(yscrollcommand=gy.set, xscrollcommand=gx.set)
        gy.grid(row=0, column=1, sticky="ns")
        gx.grid(row=1, column=0, sticky="ew")

        results_tab = ttk.Frame(bottom)
        results_tab.columnconfigure(0, weight=1)
        bottom.add(results_tab, text="Results")

        res_top = ttk.LabelFrame(results_tab, text="ROOT results")
        res_top.grid(row=0, column=0, sticky="ew", padx=4, pady=4)
        for i in range(6):
            res_top.columnconfigure(i, weight=1)

        ttk.Label(res_top, text="ROOT file pattern").grid(row=0, column=0, sticky="w", padx=4, pady=2)
        ttk.Entry(res_top, textvariable=self.results_root_pattern_var).grid(row=0, column=1, columnspan=5, sticky="ew", padx=4, pady=2)

        ttk.Label(res_top, text="ROOT result").grid(row=1, column=0, sticky="w", padx=4, pady=2)
        ttk.Combobox(res_top, textvariable=self.results_kind_var, values=RESULT_HIST_KIND_LABELS, state="readonly").grid(row=1, column=1, sticky="ew", padx=4, pady=2)

        ttk.Label(res_top, text="Centrality / R_CP index").grid(row=1, column=2, sticky="w", padx=4, pady=2)
        ttk.Combobox(res_top, textvariable=self.results_cent_var, values=["0", "1", "2"], state="readonly").grid(row=1, column=3, sticky="ew", padx=4, pady=2)

        ttk.Label(res_top, text="Iteration (enter 4 -> it3)").grid(row=1, column=4, sticky="w", padx=4, pady=2)
        ttk.Entry(res_top, textvariable=self.results_iter_var).grid(row=1, column=5, sticky="ew", padx=4, pady=2)

        ttk.Label(res_top, text="Method").grid(row=2, column=0, sticky="w", padx=4, pady=2)
        ttk.Entry(res_top, textvariable=self.results_method_var).grid(row=2, column=1, sticky="ew", padx=4, pady=2)

        ttk.Label(res_top, text="Reference run").grid(row=2, column=2, sticky="w", padx=4, pady=2)
        ttk.Entry(res_top, textvariable=self.results_ref_run_var).grid(row=2, column=3, sticky="ew", padx=4, pady=2)

        ttk.Label(res_top, text="Heatmap X").grid(row=2, column=4, sticky="w", padx=4, pady=2)
        self.results_heatmap_x_cb = ttk.Combobox(res_top, textvariable=self.results_heatmap_x_var, values=[], state="readonly")
        self.results_heatmap_x_cb.grid(row=2, column=5, sticky="ew", padx=4, pady=2)

        ttk.Label(res_top, text="Heatmap Y").grid(row=3, column=0, sticky="w", padx=4, pady=2)
        self.results_heatmap_y_cb = ttk.Combobox(res_top, textvariable=self.results_heatmap_y_var, values=[], state="readonly")
        self.results_heatmap_y_cb.grid(row=3, column=1, sticky="ew", padx=4, pady=2)

        ttk.Label(res_top, text="Heatmap bin (1-based)").grid(row=3, column=2, sticky="w", padx=4, pady=2)
        ttk.Entry(res_top, textvariable=self.results_heatmap_bin_var).grid(row=3, column=3, sticky="ew", padx=4, pady=2)

        btn_row = ttk.Frame(res_top)
        btn_row.grid(row=4, column=0, columnspan=6, sticky="ew", padx=4, pady=4)
        ttk.Button(btn_row, text="1) Overlay + ratio", command=self.plot_overlay_and_ratio).pack(side="left", padx=4)
        ttk.Button(btn_row, text="2) Mean +/- std", command=self.plot_mean_and_std).pack(side="left", padx=4)
        ttk.Button(btn_row, text="3) Bin evolution", command=self.plot_bin_evolution).pack(side="left", padx=4)
        ttk.Button(btn_row, text="4) Heatmap", command=self.plot_heatmap).pack(side="left", padx=4)
        ttk.Button(btn_row, text="5) Ref syst. by bin", command=self.plot_reference_systematics_by_bin).pack(side="left", padx=4)
        ttk.Button(btn_row, text="6) Iter +/-1", command=self.plot_iteration_stability).pack(side="left", padx=4)
        ttk.Button(btn_row, text="All 6", command=self.plot_all_results).pack(side="left", padx=10)

        ttk.Label(results_tab, textvariable=self.results_status_var).grid(row=1, column=0, sticky="ew", padx=6, pady=4)

        self.group_df = pd.DataFrame()
        self.results_syst_df = pd.DataFrame()
        self.results_iter_syst_df = pd.DataFrame()

    def configure_dynamic_columns(self):
        merged_cols = list(self.df_all.columns)
        self.summary_cols = [
            c for c in merged_cols
            if c not in (["run_id"] + BASE_FILTER_COLS + METRIC_COLS)
        ]
        self.summary_edge_cols = [c for c in self.summary_cols if "edges" in c.lower()]
        self.summary_filter_cols = [c for c in self.summary_cols if c not in self.summary_edge_cols]

        self.filter_cols = BASE_FILTER_COLS + [c for c in self.summary_filter_cols if c not in BASE_FILTER_COLS]
        self.result_cols = ["run_id"] + [c for c in self.filter_cols if c != "run_id"] + METRIC_COLS
        self.heatmap_param_choices = [
            c for c in self.summary_filter_cols
            if pd.api.types.is_numeric_dtype(self.df_all[c])
        ]

    def rebuild_filter_widgets(self):
        for child in self.cat_frame.winfo_children():
            child.destroy()

        self.filter_vars = {}
        self.filter_boxes = {}
        for idx, col in enumerate(self.filter_cols):
            self.filter_vars[col] = tk.StringVar(value="(all)")
            r = idx // 6
            c = idx % 6
            box_frame = ttk.Frame(self.cat_frame)
            box_frame.grid(row=r, column=c, sticky="ew", padx=4, pady=2)
            ttk.Label(box_frame, text=FILTER_LABELS.get(col, col)).pack(anchor="w")
            cb = ttk.Combobox(box_frame, textvariable=self.filter_vars[col], state="readonly", width=18)
            cb["values"] = ["(all)"]
            cb.pack(fill="x")
            self.filter_boxes[col] = cb

    def configure_result_tree(self):
        self.result_tree["columns"] = self.result_cols
        for col in self.result_cols:
            self.result_tree.heading(
                col,
                text=DISPLAY_LABELS.get(col, col),
                command=lambda c=col: self.sort_tree(self.result_tree, c, False),
            )
            width = 95 if "pct" not in col else 110
            if col == "observable":
                width = 140
            if col == "run_id":
                width = 90
            if "edges" in col.lower():
                width = 220
            self.result_tree.column(col, width=width, anchor="center")

    def populate_groupby_listbox(self):
        self.groupby_listbox.delete(0, "end")
        for col in self.filter_cols:
            self.groupby_listbox.insert("end", col)

    def update_results_controls(self):
        self.results_heatmap_x_cb["values"] = self.heatmap_param_choices
        self.results_heatmap_y_cb["values"] = self.heatmap_param_choices
        if self.heatmap_param_choices:
            if self.results_heatmap_x_var.get() not in self.heatmap_param_choices:
                self.results_heatmap_x_var.set(self.heatmap_param_choices[0])
            if self.results_heatmap_y_var.get() not in self.heatmap_param_choices:
                self.results_heatmap_y_var.set(self.heatmap_param_choices[min(1, len(self.heatmap_param_choices)-1)])
        else:
            self.results_heatmap_x_var.set("")
            self.results_heatmap_y_var.set("")

    def browse_summary(self):
        path = filedialog.askopenfilename(
            title="Select summary.tsv",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")]
        )
        if path:
            self.summary_path_var.set(path)

    def browse_stability(self):
        path = filedialog.askopenfilename(
            title="Select stability.tsv",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")]
        )
        if path:
            self.stability_path_var.set(path)

    def load_files(self):
        summary_path = self.summary_path_var.get().strip()
        stability_path = self.stability_path_var.get().strip()

        if not summary_path or not stability_path:
            messagebox.showerror("Error", "Select both summary.tsv and stability.tsv.")
            return

        try:
            self.df_all = merge_data(summary_path, stability_path)
            self.configure_dynamic_columns()
            self.rebuild_filter_widgets()
            self.configure_result_tree()
            self.populate_groupby_listbox()
            self.update_results_controls()
            self.populate_filter_boxes()
            self.apply_filters()
            self.init_results_defaults()
        except Exception as e:
            messagebox.showerror("Loading error", str(e))

    def populate_filter_boxes(self):
        if self.df_all.empty:
            return

        for col, cb in self.filter_boxes.items():
            vals = self.df_all[col].dropna().unique().tolist()
            vals = [str(v) for v in vals if str(v) != ""]
            vals = sorted(vals, key=self._sort_value_key)
            cb["values"] = ["(all)"] + vals
            self.filter_vars[col].set("(all)")

        obs_vals = self.df_all["observable"].dropna().astype(str).unique().tolist()
        obs_vals = [v for v in obs_vals if v != ""]
        obs_vals = sorted(obs_vals, key=self._sort_value_key)

        axis_vals = self.df_all["axis"].dropna().astype(str).unique().tolist()
        axis_vals = [v for v in axis_vals if v != ""]
        axis_vals = sorted(axis_vals, key=self._sort_value_key)

        for tag in ("A", "B"):
            self.cross_boxes[tag]["observable"]["values"] = ["(any)"] + obs_vals
            self.cross_boxes[tag]["axis"]["values"] = ["(any)"] + axis_vals
            self.cross_boxes[tag]["metric"]["values"] = ["(none)"] + METRIC_COLS

            self.cross_vars[tag]["observable"].set("(any)")
            self.cross_vars[tag]["axis"].set("(any)")
            self.cross_vars[tag]["metric"].set("(none)")
            self.cross_vars[tag]["min"].set("")
            self.cross_vars[tag]["max"].set("")

        for m in METRIC_COLS:
            self.metric_min_vars[m].set("")
            self.metric_max_vars[m].set("")

    def _sort_value_key(self, x):
        xf = maybe_float(x)
        return (0, xf) if xf is not None else (1, str(x))

    def _parse_metric_bounds(self, min_text: str, max_text: str, label: str):
        out_min = None
        out_max = None

        if min_text.strip():
            try:
                out_min = float(min_text)
            except ValueError:
                raise ValueError(f"Invalid minimum for {label}: {min_text}")

        if max_text.strip():
            try:
                out_max = float(max_text)
            except ValueError:
                raise ValueError(f"Invalid maximum for {label}: {max_text}")

        return out_min, out_max

    def _row_mask_for_cross_filter(self, df: pd.DataFrame, tag: str) -> pd.Series:
        cfg = self.cross_vars[tag]
        metric = cfg["metric"].get().strip()

        mask = pd.Series(True, index=df.index)

        obs = cfg["observable"].get().strip()
        axis = cfg["axis"].get().strip()

        if obs and obs != "(any)":
            mask &= (df["observable"].astype(str) == obs)

        if axis and axis != "(any)":
            mask &= (df["axis"].astype(str) == axis)

        if metric and metric != "(none)":
            vmin, vmax = self._parse_metric_bounds(cfg["min"].get(), cfg["max"].get(), f"Filter {tag} / {metric}")
            if vmin is not None:
                mask &= (df[metric] >= vmin)
            if vmax is not None:
                mask &= (df[metric] <= vmax)

        return mask

    def _cross_filter_active(self, tag: str) -> bool:
        return self.cross_vars[tag]["metric"].get().strip() not in ("", "(none)")

    def apply_filters(self):
        if self.df_all.empty:
            self.status_var.set("No data loaded.")
            return

        df = self.df_all.copy()

        try:
            for col, var in self.filter_vars.items():
                val = var.get().strip()
                if val and val != "(all)":
                    df = df[df[col].astype(str) == val]

            for metric in METRIC_COLS:
                vmin_txt = self.metric_min_vars[metric].get().strip()
                vmax_txt = self.metric_max_vars[metric].get().strip()
                vmin, vmax = self._parse_metric_bounds(vmin_txt, vmax_txt, metric)

                if vmin is not None:
                    df = df[df[metric] >= vmin]
                if vmax is not None:
                    df = df[df[metric] <= vmax]

            active_a = self._cross_filter_active("A")
            active_b = self._cross_filter_active("B")
            groups_kept = None

            if active_a or active_b:
                if df.empty:
                    self.df_filtered = df.reset_index(drop=True)
                    self.fill_result_tree()
                    self.status_var.set(f"Loaded {len(self.df_all)} rows, filter selected 0 rows.")
                    return

                mask_a = self._row_mask_for_cross_filter(df, "A") if active_a else pd.Series(True, index=df.index)
                mask_b = self._row_mask_for_cross_filter(df, "B") if active_b else pd.Series(True, index=df.index)

                groups = df[["run_id", "cent"]].copy()
                groups["a"] = mask_a.values if active_a else True
                groups["b"] = mask_b.values if active_b else True

                group_summary = groups.groupby(["run_id", "cent"], dropna=False).agg({"a": "any", "b": "any"}).reset_index()
                group_summary["keep"] = group_summary["a"] & group_summary["b"]

                groups_kept = set(zip(group_summary.loc[group_summary["keep"], "run_id"], group_summary.loc[group_summary["keep"], "cent"]))
                group_mask = pd.Series([(run_id, cent) in groups_kept for run_id, cent in zip(df["run_id"], df["cent"])], index=df.index)

                if self.cross_display_mode_var.get().strip() == CROSS_DISPLAY_MODES[1] and active_b:
                    df = df[group_mask & mask_b]
                else:
                    df = df[group_mask]

            self.df_filtered = df.reset_index(drop=True)
            self.fill_result_tree()

            msg = f"Loaded {len(self.df_all)} rows, filter selected {len(self.df_filtered)} rows."
            if groups_kept is not None:
                msg += f" Cross-filter kept {len(groups_kept)} groups (run_id, cent)."
            self.status_var.set(msg)

        except ValueError as e:
            messagebox.showerror("Error", str(e))
            return

    def reset_filters(self):
        for col in self.filter_vars:
            self.filter_vars[col].set("(all)")
        for m in METRIC_COLS:
            self.metric_min_vars[m].set("")
            self.metric_max_vars[m].set("")

        for tag in ("A", "B"):
            self.cross_vars[tag]["observable"].set("(any)")
            self.cross_vars[tag]["axis"].set("(any)")
            self.cross_vars[tag]["metric"].set("(none)")
            self.cross_vars[tag]["min"].set("")
            self.cross_vars[tag]["max"].set("")

        self.cross_display_mode_var.set(CROSS_DISPLAY_MODES[0])
        self.apply_filters()

    def fill_result_tree(self):
        self.result_tree.delete(*self.result_tree.get_children())

        for idx, row in self.df_filtered.iterrows():
            vals = [format_value(row.get(c, "")) for c in self.result_cols]
            self.result_tree.insert("", "end", iid=str(idx), values=vals)

        self.detail_text.delete("1.0", "end")
        self.group_tree.delete(*self.group_tree.get_children())
        self.group_df = pd.DataFrame()

    def on_result_select(self, _event=None):
        sel = self.result_tree.selection()
        if not sel:
            return

        idx = int(sel[0])
        if idx < 0 or idx >= len(self.df_filtered):
            return

        row = self.df_filtered.iloc[idx]

        lines = []
        for key in ["run_id", "cent", "dim", "observable", "axis"]:
            if key in row.index:
                lines.append(f"{key}: {format_value(row.get(key))}")

        if self.summary_cols:
            lines.append("")
            lines.append("Summary columns:")
            for c in self.summary_cols:
                lines.append(f"  {c}: {format_value(row.get(c))}")

        lines.append("")
        lines.append("Metrics:")
        for m in METRIC_COLS:
            if m in row.index:
                lines.append(f"  {METRIC_LABELS[m]}: {format_value(row.get(m))}")

        self.detail_text.delete("1.0", "end")
        self.detail_text.insert("1.0", "\n".join(lines))

    def compute_group_summary(self):
        if self.df_filtered.empty:
            messagebox.showinfo("Info", "There are no filtered data.")
            return

        selected = [self.groupby_listbox.get(i) for i in self.groupby_listbox.curselection()]
        metric = self.group_metric_var.get()
        agg = self.group_agg_var.get()

        if metric not in METRIC_COLS:
            messagebox.showerror("Error", "Select a valid metric.")
            return

        df = self.df_filtered.copy()

        if selected:
            grp = (
                df.groupby(selected, dropna=False)[metric]
                .agg([("count", "size"), (f"{agg}_{metric}", agg)])
                .reset_index()
            )
        else:
            value = getattr(df[metric], agg)()
            grp = pd.DataFrame({
                "count": [len(df)],
                f"{agg}_{metric}": [value],
            })

        self.group_df = grp
        self.fill_group_tree(grp)

    def fill_group_tree(self, df: pd.DataFrame):
        self.group_tree.delete(*self.group_tree.get_children())

        cols = list(df.columns)
        self.group_tree["columns"] = cols

        for col in cols:
            self.group_tree.heading(
                col,
                text=DISPLAY_LABELS.get(col, col),
                command=lambda c=col: self.sort_tree(self.group_tree, c, False),
            )
            self.group_tree.column(col, width=120, anchor="center")

        for idx, (_, row) in enumerate(df.iterrows()):
            vals = [format_value(row[c]) for c in cols]
            self.group_tree.insert("", "end", iid=str(idx), values=vals)

    def export_filtered_csv(self):
        if self.df_filtered.empty:
            messagebox.showinfo("Info", "Nothing to export.")
            return

        path = filedialog.asksaveasfilename(
            title="Save filtered data",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not path:
            return

        self.df_filtered.to_csv(path, index=False)
        messagebox.showinfo("Done", f"Filtered data saved to:\n{path}")

    def export_group_csv(self):
        if self.group_df.empty:
            messagebox.showinfo("Info", "Compute the aggregation first.")
            return

        path = filedialog.asksaveasfilename(
            title="Save aggregation",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not path:
            return

        self.group_df.to_csv(path, index=False)
        messagebox.showinfo("Done", f"Aggregation saved to:\n{path}")

    def sort_tree(self, tree: ttk.Treeview, col: str, reverse: bool):
        data = [(tree.set(k, col), k) for k in tree.get_children("")]

        def keyfunc(item):
            value = item[0]
            fv = maybe_float(value)
            if fv is not None:
                return (0, fv)
            return (1, value)

        data.sort(key=keyfunc, reverse=reverse)

        for idx, (_, k) in enumerate(data):
            tree.move(k, "", idx)

        tree.heading(col, command=lambda: self.sort_tree(tree, col, not reverse))


    def init_results_defaults(self):
        run_ids = self.get_filtered_run_ids()
        if run_ids and not self.results_ref_run_var.get().strip():
            self.results_ref_run_var.set(run_ids[0])

    def get_filtered_run_ids(self):
        if self.df_filtered.empty or "run_id" not in self.df_filtered.columns:
            return []
        run_ids = sorted(self.df_filtered["run_id"].dropna().astype(str).unique().tolist(), key=run_id_to_int)
        return run_ids

    def get_results_hist_name(self):
        try:
            cent = int(self.results_cent_var.get().strip())
            iteration = int(self.results_iter_var.get().strip())
        except ValueError:
            raise ValueError("Centrality and iteration must be integers. Iteration is entered as 1-based: for example, 4 loads *_it3_*.")
        kind_label = self.results_kind_var.get().strip()
        method = self.results_method_var.get().strip()
        if not kind_label or not method:
            raise ValueError("ROOT result and Method must be filled in.")
        if kind_label not in RESULT_HIST_KIND_MAP:
            raise ValueError(f"Unknown ROOT result: {kind_label}")
        return build_result_hist_name(kind_label, cent, iteration, method)

    def collect_histograms_for_filtered_runs(self):
        run_ids = self.get_filtered_run_ids()
        if not run_ids:
            raise ValueError("No runs are available after the current filters.")

        hist_name = self.get_results_hist_name()
        pattern = self.results_root_pattern_var.get().strip()
        if not pattern:
            raise ValueError("ROOT file pattern must not be empty.")

        ref_run = normalize_run_id(self.results_ref_run_var.get())
        if not ref_run:
            ref_run = run_ids[0]
            self.results_ref_run_var.set(ref_run)

        data = {}
        skipped = []

        self.results_status_var.set(f"Loading {len(run_ids)} ROOT files for {hist_name}...")
        self.update_idletasks()

        for idx, run_id in enumerate(run_ids, start=1):
            root_path = pattern.format(run_id=run_id)
            try:
                edges, values = read_hist_arrays(root_path, hist_name)
                data[run_id] = {
                    "edges": np.asarray(edges, dtype=float),
                    "values": np.asarray(values, dtype=float),
                    "path": root_path,
                }
            except Exception as e:
                skipped.append(f"{run_id}: {e}")

            if idx == 1 or idx % 10 == 0 or idx == len(run_ids):
                self.results_status_var.set(f"Loaded {idx}/{len(run_ids)} files, successful {len(data)}, skipped {len(skipped)}")
                self.update_idletasks()

        if ref_run not in data:
            raise ValueError(
                f"Reference run {ref_run} could not be loaded.\n"
                + ("\n".join(skipped[:8]) if skipped else "")
            )

        ref_edges = data[ref_run]["edges"]
        ref_values = data[ref_run]["values"]

        mapped = {}
        for run_id, item in data.items():
            mapped[run_id] = rebin_density_to_edges(item["edges"], item["values"], ref_edges)

        meta_cols = ["run_id"] + [c for c in self.summary_filter_cols if c in self.df_filtered.columns]
        meta = self.df_filtered.drop_duplicates(subset=["run_id"])[meta_cols].copy()
        meta["run_num"] = meta["run_id"].map(run_id_to_int)
        meta = meta.sort_values("run_num")

        self.results_status_var.set(
            f"Histogram {hist_name}: loaded {len(data)} runs, skipped {len(skipped)}. "
            f"Reference = {ref_run}."
        )

        return {
            "hist_name": hist_name,
            "ref_run": ref_run,
            "ref_edges": ref_edges,
            "ref_values": ref_values,
            "data": data,
            "mapped": mapped,
            "meta": meta,
            "skipped": skipped,
        }

    def get_results_settings(self):
        try:
            cent = int(self.results_cent_var.get().strip())
            iteration = int(self.results_iter_var.get().strip())
        except ValueError:
            raise ValueError("Centrality and iteration must be integers. Iteration is entered as 1-based: for example, 4 loads *_it3_*.")

        kind_label = self.results_kind_var.get().strip()
        method = self.results_method_var.get().strip()
        if not kind_label or not method:
            raise ValueError("ROOT result and Method must be filled in.")
        if kind_label not in RESULT_HIST_KIND_MAP:
            raise ValueError(f"Unknown ROOT result: {kind_label}")
        if iteration < 1:
            raise ValueError("Iteration is entered as 1-based: for example, 4 loads the ROOT histogram *_it3_*")

        return cent, iteration, kind_label, method

    def collect_iteration_variations_for_reference_run(self):
        """Load the selected iteration and its +/-1 neighbours for the selected reference run."""
        run_ids = self.get_filtered_run_ids()
        pattern = self.results_root_pattern_var.get().strip()
        if not pattern:
            raise ValueError("ROOT file pattern must not be empty.")

        ref_run = normalize_run_id(self.results_ref_run_var.get())
        if not ref_run:
            if not run_ids:
                raise ValueError("No reference run is selected and no runs are available after filtering.")
            ref_run = run_ids[0]
            self.results_ref_run_var.set(ref_run)

        cent, iteration_nominal, kind_label, method = self.get_results_settings()
        iter_candidates = [iteration_nominal - 1, iteration_nominal, iteration_nominal + 1]
        iter_candidates = [it for it in iter_candidates if it >= 1]

        root_path = pattern.format(run_id=ref_run)
        data = {}
        skipped = []

        self.results_status_var.set(
            f"Loading iterations {iter_candidates} for {ref_run}..."
        )
        self.update_idletasks()

        for it in iter_candidates:
            hist_name = build_result_hist_name(kind_label, cent, it, method)
            try:
                edges, values = read_hist_arrays(root_path, hist_name)
                data[it] = {
                    "hist_name": hist_name,
                    "edges": np.asarray(edges, dtype=float),
                    "values": np.asarray(values, dtype=float),
                    "path": root_path,
                }
            except Exception as e:
                skipped.append(f"iteration {it} ({hist_name}): {e}")

        if iteration_nominal not in data:
            raise ValueError(
                f"Nominal iteration {iteration_nominal} could not be loaded.\n"
                + ("\n".join(skipped[:8]) if skipped else "")
            )

        variation_iters = [it for it in sorted(data.keys()) if it != iteration_nominal]
        if not variation_iters:
            raise ValueError(
                "Could not load any neighboring iteration.\n"
                + ("\n".join(skipped[:8]) if skipped else "")
            )

        ref_edges = data[iteration_nominal]["edges"]
        ref_values = data[iteration_nominal]["values"]

        mapped = {}
        for it, item in data.items():
            mapped[it] = rebin_density_to_edges(item["edges"], item["values"], ref_edges)

        self.results_status_var.set(
            f"Iteration stability: run {ref_run}, nominal iter {iteration_nominal}, "
            f"loaded {sorted(data.keys())}, skipped {len(skipped)}."
        )

        return {
            "ref_run": ref_run,
            "root_path": root_path,
            "nominal_iteration": iteration_nominal,
            "kind_label": kind_label,
            "method": method,
            "cent": cent,
            "nominal_hist_name": data[iteration_nominal]["hist_name"],
            "ref_edges": ref_edges,
            "ref_values": ref_values,
            "data": data,
            "mapped": mapped,
            "variation_iters": variation_iters,
            "skipped": skipped,
        }

    def compute_iteration_stability_by_bin(self, pack):
        ref_edges = np.asarray(pack["ref_edges"], dtype=float)
        ref_values = np.asarray(pack["ref_values"], dtype=float)
        mapped = pack["mapped"]
        nominal_iteration = pack["nominal_iteration"]
        centers = centers_from_edges(ref_edges)

        minus_it = nominal_iteration - 1
        plus_it = nominal_iteration + 1
        minus_vals = mapped.get(minus_it)
        plus_vals = mapped.get(plus_it)

        rows = []
        for ib in range(len(ref_values)):
            ref_val = ref_values[ib]

            def get_var(vals):
                if vals is None or ib >= len(vals):
                    return np.nan, np.nan
                val = vals[ib]
                if not np.isfinite(val) or not np.isfinite(ref_val) or ref_val == 0:
                    return val, np.nan
                return val, 100.0 * (val / ref_val - 1.0)

            val_minus, diff_minus = get_var(minus_vals)
            val_plus, diff_plus = get_var(plus_vals)

            diffs = np.asarray([d for d in (diff_minus, diff_plus) if np.isfinite(d)], dtype=float)
            if diffs.size == 0:
                envelope_down_pct = np.nan
                envelope_up_pct = np.nan
                max_abs_pct = np.nan
            else:
                env_min = float(np.min(diffs))
                env_max = float(np.max(diffs))
                envelope_down_pct = max(0.0, -env_min)
                envelope_up_pct = max(0.0, env_max)
                max_abs_pct = float(np.max(np.abs(diffs)))

            rows.append({
                "bin": ib + 1,
                "x_low": ref_edges[ib],
                "x_high": ref_edges[ib + 1],
                "x_center": centers[ib],
                "nominal_iter": nominal_iteration,
                "nominal_value": ref_val,
                "iter_minus": minus_it if minus_it in mapped else pd.NA,
                "value_minus": val_minus,
                "diff_minus_pct": diff_minus,
                "iter_plus": plus_it if plus_it in mapped else pd.NA,
                "value_plus": val_plus,
                "diff_plus_pct": diff_plus,
                "n_var": int(diffs.size),
                "envelope_down_pct": envelope_down_pct,
                "envelope_up_pct": envelope_up_pct,
                "max_abs_pct": max_abs_pct,
                "envelope_down_abs": abs(ref_val) * envelope_down_pct / 100.0 if np.isfinite(envelope_down_pct) else np.nan,
                "envelope_up_abs": abs(ref_val) * envelope_up_pct / 100.0 if np.isfinite(envelope_up_pct) else np.nan,
            })

        return pd.DataFrame(rows)

    def show_iteration_stability_table(self, iter_df: pd.DataFrame, pack):
        self.results_iter_syst_df = iter_df.copy()

        win = tk.Toplevel(self)
        win.title(
            f"Iteration stability: {pack['ref_run']} / iter {pack['nominal_iteration']}"
        )
        win.geometry("1550x520")
        win.columnconfigure(0, weight=1)
        win.rowconfigure(1, weight=1)

        loaded = sorted(pack["data"].keys())
        info = (
            f"Run = {pack['ref_run']} | Histogram = {pack['nominal_hist_name']} | "
            f"nominal iter = {pack['nominal_iteration']} | loaded iterations = {loaded} | "
            "envelope columns are absolute magnitudes"
        )
        ttk.Label(win, text=info).grid(row=0, column=0, sticky="ew", padx=8, pady=6)

        cols = list(iter_df.columns)
        tree = ttk.Treeview(win, columns=cols, show="headings")
        tree.grid(row=1, column=0, sticky="nsew", padx=8, pady=4)

        ysb = ttk.Scrollbar(win, orient="vertical", command=tree.yview)
        xsb = ttk.Scrollbar(win, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=ysb.set, xscrollcommand=xsb.set)
        ysb.grid(row=1, column=1, sticky="ns")
        xsb.grid(row=2, column=0, sticky="ew", padx=8)

        for col in cols:
            tree.heading(col, text=col, command=lambda c=col: self.sort_tree(tree, c, False))
            width = 95
            if col in ("bin", "n_var", "iter_minus", "iter_plus", "nominal_iter"):
                width = 75
            elif col in ("nominal_value", "value_minus", "value_plus", "envelope_down_abs", "envelope_up_abs"):
                width = 130
            elif "pct" in col:
                width = 130
            tree.column(col, width=width, anchor="center")

        for idx, (_, row) in enumerate(iter_df.iterrows()):
            vals = [format_value(row[c]) for c in cols]
            tree.insert("", "end", iid=str(idx), values=vals)

        btns = ttk.Frame(win)
        btns.grid(row=3, column=0, sticky="ew", padx=8, pady=8)
        ttk.Button(btns, text="Copy TSV", command=lambda: self.copy_dataframe_to_clipboard(iter_df)).pack(side="left", padx=4)
        ttk.Button(btns, text="Export CSV", command=lambda: self.export_dataframe_csv(iter_df, "iteration_stability_systematics.csv")).pack(side="left", padx=4)

    def plot_iteration_stability(self):
        try:
            pack = self.collect_iteration_variations_for_reference_run()
            iter_df = self.compute_iteration_stability_by_bin(pack)
        except Exception as e:
            messagebox.showerror("Error", str(e))
            return

        ref_edges = pack["ref_edges"]
        ref_values = np.asarray(pack["ref_values"], dtype=float)
        centers = centers_from_edges(ref_edges)
        nominal_iteration = pack["nominal_iteration"]
        mapped = pack["mapped"]
        data = pack["data"]

        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(11, 8), sharex=True,
            gridspec_kw={"height_ratios": [2.2, 1.0]}
        )

        for it in sorted(data.keys()):
            item = data[it]
            if it == nominal_iteration:
                step_plot(ax1, item["edges"], item["values"], label=f"iter {it} (nominal)", linewidth=2.2, color="black")
            else:
                step_plot(ax1, item["edges"], item["values"], label=f"iter {it}", linewidth=1.4, alpha=0.8)

        ax1.set_title(
            f"Iteration stability: {pack['nominal_hist_name']} / {pack['ref_run']}"
        )
        ax1.set_ylabel("bin content")
        ax1.set_yscale("log")
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        for it in sorted(data.keys()):
            if it == nominal_iteration:
                continue
            vals = mapped[it]
            with np.errstate(divide='ignore', invalid='ignore'):
                diff_pct = np.where(
                    np.isfinite(vals) & np.isfinite(ref_values) & (ref_values != 0),
                    100.0 * (vals / ref_values - 1.0),
                    np.nan,
                )
            ax2.plot(centers, diff_pct, marker="o", linewidth=1.0, label=f"iter {it} / iter {nominal_iteration} - 1")

        env_down = iter_df["envelope_down_pct"].to_numpy(dtype=float)
        env_up = iter_df["envelope_up_pct"].to_numpy(dtype=float)
        ax2.plot(centers, env_up, marker="o", linewidth=1.0, alpha=0.55, label="envelope up [%]")
        ax2.plot(centers, -env_down, marker="o", linewidth=1.0, alpha=0.55, label="envelope down [%]")
        ax2.axhline(0.0, linestyle="--", linewidth=1.0)
        ax2.set_xlim(ref_edges[0], ref_edges[-1])
        ax2.set_xlabel(f"Nominal iteration {nominal_iteration} bins")
        ax2.set_ylabel("relative variation [%]")
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=8, ncol=2)

        fig.tight_layout()
        plt.show(block=False)

        self.show_iteration_stability_table(iter_df, pack)


    def plot_overlay_and_ratio(self):
        try:
            pack = self.collect_histograms_for_filtered_runs()
        except Exception as e:
            messagebox.showerror("Error", str(e))
            return

        ref_run = pack["ref_run"]
        ref_edges = pack["ref_edges"]
        ref_values = pack["ref_values"]
        data = pack["data"]
        mapped = pack["mapped"]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True, gridspec_kw={"height_ratios": [2.2, 1.0]})

        sorted_runs = sorted(data.keys(), key=run_id_to_int)
        for run_id in sorted_runs:
            item = data[run_id]
            if run_id == ref_run:
                step_plot(ax1, item["edges"], item["values"], label=f"{run_id} (ref)", linewidth=2.2, color="black")
            else:
                step_plot(ax1, item["edges"], item["values"], label=run_id, linewidth=1.0, alpha=0.35)

        ax1.set_title(f"Overlay spectra: {pack['hist_name']}")
        ax1.set_ylabel("bin content")
        ax1.set_yscale("log")
        ax1.grid(True, alpha=0.3)
        if len(sorted_runs) <= 18:
            ax1.legend(fontsize=8, ncol=2)

        ref_centers = centers_from_edges(ref_edges)
        for run_id in sorted_runs:
            y = mapped[run_id]
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = np.where(np.isfinite(y) & (ref_values != 0), y / ref_values, np.nan)
            if run_id == ref_run:
                ax2.plot(ref_centers, ratio, marker='o', markersize=3, linewidth=1.8, color='black')
            else:
                ax2.plot(ref_centers, ratio, marker='o', markersize=2.5, linewidth=1.0, alpha=0.45)

        ax2.axhline(1.0, color='black', linestyle='--', linewidth=1.0)
        ax2.set_xlabel(f"Reference bins of {ref_run}")
        ax2.set_ylabel("run/ref")
        ax2.grid(True, alpha=0.3)

        fig.tight_layout()
        plt.show()

    def plot_mean_and_std(self):
        try:
            pack = self.collect_histograms_for_filtered_runs()
        except Exception as e:
            messagebox.showerror("Error", str(e))
            return

        ref_edges = pack["ref_edges"]
        mapped = pack["mapped"]
        arr = np.array([mapped[r] for r in sorted(mapped.keys(), key=run_id_to_int)], dtype=float)
        mean = np.nanmean(arr, axis=0)
        std = np.nanstd(arr, axis=0)
        centers = centers_from_edges(ref_edges)

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True, gridspec_kw={"height_ratios": [2.2, 1.0]})
        step_plot(ax1, ref_edges, mean, label='mean', linewidth=1.8)

        y_low = np.clip(mean - std, 1e-300, None)
        y_high = np.clip(mean + std, 1e-300, None)
        ax1.fill_between(centers, y_low, y_high, alpha=0.3, label='mean ± std')
        ax1.set_title(f"Mean ± std on reference binning: {pack['hist_name']} ({pack['ref_run']})")
        ax1.set_ylabel("mapped content")
        ax1.set_yscale("log")
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        rel = np.where(np.isfinite(mean) & (mean != 0), 100.0 * std / np.abs(mean), np.nan)
        ax2.plot(centers, rel, marker='o', linewidth=1.0)
        ax1.set_xlim(ref_edges[0], ref_edges[-1])
        ax2.set_xlim(ref_edges[0], ref_edges[-1])
        ax2.set_xlabel(f"Reference bins of {pack['ref_run']}")
        ax2.set_ylabel("std/|mean| [%]")
        ax2.grid(True, alpha=0.3)

        fig.tight_layout()
        plt.show()

    def plot_bin_evolution(self):
        try:
            pack = self.collect_histograms_for_filtered_runs()
        except Exception as e:
            messagebox.showerror("Error", str(e))
            return

        meta = pack["meta"]
        ref_edges = pack["ref_edges"]
        run_ids = [r for r in meta["run_id"].tolist() if r in pack["mapped"]]
        if not run_ids:
            messagebox.showerror("Error", "No loaded data were found for bin evolution.")
            return

        x = [run_id_to_int(r) for r in run_ids]
        arr = np.array([pack["mapped"][r] for r in run_ids], dtype=float)

        fig, ax = plt.subplots(figsize=(12, 7))
        for ib in range(arr.shape[1]):
            label = f"bin {ib+1}: [{ref_edges[ib]}, {ref_edges[ib+1]})"
            ax.plot(x, arr[:, ib], marker='o', markersize=2.5, linewidth=1.0, label=label)

        ax.set_title(f"Bin evolution vs run number: {pack['hist_name']} on ref bins")
        ax.set_xlabel("run number")
        ax.set_yscale("log")
        ax.set_ylabel("mapped content")
        ax.grid(True, alpha=0.3)
        if arr.shape[1] <= 14:
            ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        plt.show()

    def compute_reference_systematics_by_bin(self, pack):
        ref_run = pack["ref_run"]
        ref_edges = np.asarray(pack["ref_edges"], dtype=float)
        ref_values = np.asarray(pack["ref_values"], dtype=float)
        mapped = pack["mapped"]

        alt_runs = [r for r in sorted(mapped.keys(), key=run_id_to_int) if r != ref_run]
        if not alt_runs:
            raise ValueError("At least one alternative run outside the reference is needed to compute the systematic uncertainty.")

        centers = centers_from_edges(ref_edges)
        rows = []

        for ib in range(len(ref_values)):
            ref_val = ref_values[ib]
            diffs = []

            for run_id in alt_runs:
                vals = mapped.get(run_id)
                if vals is None or ib >= len(vals):
                    continue

                val = vals[ib]
                if not np.isfinite(val) or not np.isfinite(ref_val) or ref_val == 0:
                    continue

                diffs.append(100.0 * (val / ref_val - 1.0))

            diffs = np.asarray(diffs, dtype=float)
            if diffs.size == 0:
                sigma_pct = np.nan
                sigma_down_pct = np.nan
                sigma_up_pct = np.nan
                envelope_down_pct = np.nan
                envelope_up_pct = np.nan
                max_abs_pct = np.nan
                n_down = 0
                n_up = 0
            else:
                # Symmetric RMS deviation with respect to the selected reference binning,
                # not the standard deviation with respect to the mean.
                sigma_pct = float(np.sqrt(np.mean(diffs * diffs)))

                # Asymmetric RMS variant: negative and positive deviations separately.
                # Store values as positive magnitudes of the systematic uncertainty.
                diffs_down = diffs[diffs < 0.0]
                diffs_up = diffs[diffs > 0.0]
                n_down = int(diffs_down.size)
                n_up = int(diffs_up.size)
                sigma_down_pct = float(np.sqrt(np.mean(diffs_down * diffs_down))) if n_down > 0 else 0.0
                sigma_up_pct = float(np.sqrt(np.mean(diffs_up * diffs_up))) if n_up > 0 else 0.0

                # Store the envelope as positive magnitudes as well, not as signed values.
                env_min = float(np.min(diffs))
                env_max = float(np.max(diffs))
                envelope_down_pct = max(0.0, -env_min)
                envelope_up_pct = max(0.0, env_max)
                max_abs_pct = float(np.max(np.abs(diffs)))

            rows.append({
                "bin": ib + 1,
                "x_low": ref_edges[ib],
                "x_high": ref_edges[ib + 1],
                "x_center": centers[ib],
                "ref_value": ref_val,
                "n_alt": int(diffs.size),
                "n_down": n_down,
                "n_up": n_up,
                "sigma_pct": sigma_pct,
                "sigma_down_pct": sigma_down_pct,
                "sigma_up_pct": sigma_up_pct,
                "envelope_down_pct": envelope_down_pct,
                "envelope_up_pct": envelope_up_pct,
                "max_abs_pct": max_abs_pct,
                "sigma_abs": abs(ref_val) * sigma_pct / 100.0 if np.isfinite(sigma_pct) else np.nan,
                "sigma_down_abs": abs(ref_val) * sigma_down_pct / 100.0 if np.isfinite(sigma_down_pct) else np.nan,
                "sigma_up_abs": abs(ref_val) * sigma_up_pct / 100.0 if np.isfinite(sigma_up_pct) else np.nan,
                "envelope_down_abs": abs(ref_val) * envelope_down_pct / 100.0 if np.isfinite(envelope_down_pct) else np.nan,
                "envelope_up_abs": abs(ref_val) * envelope_up_pct / 100.0 if np.isfinite(envelope_up_pct) else np.nan,
            })

        return pd.DataFrame(rows)

    def copy_dataframe_to_clipboard(self, df: pd.DataFrame):
        if df is None or df.empty:
            messagebox.showinfo("Info", "Nothing to copy.")
            return
        txt = df.to_csv(sep="\t", index=False, float_format="%.8g")
        self.clipboard_clear()
        self.clipboard_append(txt)
        self.update()
        messagebox.showinfo("Done", "The table was copied to the clipboard as TSV.")

    def export_dataframe_csv(self, df: pd.DataFrame, default_name="table.csv"):
        if df is None or df.empty:
            messagebox.showinfo("Info", "Nothing to export.")
            return

        path = filedialog.asksaveasfilename(
            title="Save table",
            initialfile=default_name,
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not path:
            return

        df.to_csv(path, index=False, float_format="%.10g")
        messagebox.showinfo("Done", f"Table saved to:\n{path}")

    def show_reference_systematics_table(self, syst_df: pd.DataFrame, pack):
        self.results_syst_df = syst_df.copy()

        win = tk.Toplevel(self)
        win.title(f"Reference systematics by bin: {pack['hist_name']} / {pack['ref_run']}")
        win.geometry("1500x520")
        win.columnconfigure(0, weight=1)
        win.rowconfigure(1, weight=1)

        info = (
            f"Reference = {pack['ref_run']} | Histogram = {pack['hist_name']} | "
            "reference run is excluded from RMS/envelope | envelope columns are absolute magnitudes"
        )
        ttk.Label(win, text=info).grid(row=0, column=0, sticky="ew", padx=8, pady=6)

        cols = list(syst_df.columns)
        tree = ttk.Treeview(win, columns=cols, show="headings")
        tree.grid(row=1, column=0, sticky="nsew", padx=8, pady=4)

        ysb = ttk.Scrollbar(win, orient="vertical", command=tree.yview)
        xsb = ttk.Scrollbar(win, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=ysb.set, xscrollcommand=xsb.set)
        ysb.grid(row=1, column=1, sticky="ns")
        xsb.grid(row=2, column=0, sticky="ew", padx=8)

        for col in cols:
            tree.heading(col, text=col, command=lambda c=col: self.sort_tree(tree, c, False))
            width = 95
            if col in ("bin", "n_alt"):
                width = 60
            elif col in ("ref_value", "sigma_abs", "sys_down_abs", "sys_up_abs"):
                width = 130
            elif "pct" in col:
                width = 130
            tree.column(col, width=width, anchor="center")

        for idx, (_, row) in enumerate(syst_df.iterrows()):
            vals = [format_value(row[c]) for c in cols]
            tree.insert("", "end", iid=str(idx), values=vals)

        btns = ttk.Frame(win)
        btns.grid(row=3, column=0, sticky="ew", padx=8, pady=8)
        ttk.Button(btns, text="Copy TSV", command=lambda: self.copy_dataframe_to_clipboard(syst_df)).pack(side="left", padx=4)
        ttk.Button(btns, text="Export CSV", command=lambda: self.export_dataframe_csv(syst_df, "reference_binning_systematics.csv")).pack(side="left", padx=4)

    def plot_reference_systematics_by_bin(self):
        try:
            pack = self.collect_histograms_for_filtered_runs()
            syst_df = self.compute_reference_systematics_by_bin(pack)
        except Exception as e:
            messagebox.showerror("Error", str(e))
            return

        ref_edges = pack["ref_edges"]
        ref_values = np.asarray(pack["ref_values"], dtype=float)
        centers = centers_from_edges(ref_edges)

        sigma_pct = syst_df["sigma_pct"].to_numpy(dtype=float)
        sigma_down_pct = syst_df["sigma_down_pct"].to_numpy(dtype=float)
        sigma_up_pct = syst_df["sigma_up_pct"].to_numpy(dtype=float)
        env_down_pct = syst_df["envelope_down_pct"].to_numpy(dtype=float)
        env_up_pct = syst_df["envelope_up_pct"].to_numpy(dtype=float)

        sigma_low = ref_values * (1.0 - sigma_down_pct / 100.0)
        sigma_high = ref_values * (1.0 + sigma_up_pct / 100.0)
        env_low = ref_values * (1.0 - env_down_pct / 100.0)
        env_high = ref_values * (1.0 + env_up_pct / 100.0)

        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(11, 8), sharex=True,
            gridspec_kw={"height_ratios": [2.2, 1.0]}
        )

        step_plot(ax1, ref_edges, ref_values, label=f"reference {pack['ref_run']}", linewidth=1.8)
        fill_between_bins(ax1, ref_edges, np.clip(env_low, 1e-300, None), np.clip(env_high, 1e-300, None), alpha=0.25, label="envelope")
        fill_between_bins(ax1, ref_edges, np.clip(sigma_low, 1e-300, None), np.clip(sigma_high, 1e-300, None), alpha=0.25, label="asym. RMS to ref")
        ax1.set_title(f"Binning uncertainty vs reference: {pack['hist_name']} ({pack['ref_run']})")
        ax1.set_ylabel("reference content")
        ax1.set_yscale("log")
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        ax2.plot(centers, sigma_pct, marker="o", linewidth=1.0, label="symmetric RMS [%]")
        ax2.plot(centers, sigma_up_pct, marker="o", linewidth=1.0, label="asym. RMS up [%]")
        ax2.plot(centers, -sigma_down_pct, marker="o", linewidth=1.0, label="asym. RMS down [%]")
        ax2.plot(centers, env_up_pct, marker="o", linewidth=1.0, alpha=0.55, label="envelope up [%]")
        ax2.plot(centers, -env_down_pct, marker="o", linewidth=1.0, alpha=0.55, label="envelope down [%]")
        ax2.axhline(0.0, linestyle="--", linewidth=1.0)
        ax2.set_xlim(ref_edges[0], ref_edges[-1])
        ax2.set_xlabel(f"Reference bins of {pack['ref_run']}")
        ax2.set_ylabel("relative variation [%]")
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=8, ncol=2)

        fig.tight_layout()
        plt.show(block=False)

        self.show_reference_systematics_table(syst_df, pack)

    def plot_heatmap(self):
        try:
            pack = self.collect_histograms_for_filtered_runs()
        except Exception as e:
            messagebox.showerror("Error", str(e))
            return

        try:
            bin_index = int(self.results_heatmap_bin_var.get().strip())
        except ValueError:
            messagebox.showerror("Error", "Heatmap bin must be an integer (1-based).")
            return

        ref_edges = pack["ref_edges"]
        n_bins = len(ref_edges) - 1
        if bin_index < 1 or bin_index > n_bins:
            messagebox.showerror("Error", f"Heatmap bin must be between 1 a {n_bins}.")
            return

        x_col = self.results_heatmap_x_var.get().strip()
        y_col = self.results_heatmap_y_var.get().strip()
        if x_col not in self.heatmap_param_choices or y_col not in self.heatmap_param_choices:
            messagebox.showerror("Error", "Select valid heatmap parameters.")
            return
        if x_col == y_col:
            messagebox.showerror("Error", "Heatmap X and Y must be different parameters.")
            return

        ref_val = pack["ref_values"][bin_index - 1]
        rows = []
        meta = pack["meta"].set_index("run_id")
        for run_id, vals in pack["mapped"].items():
            if run_id not in meta.index:
                continue
            row = meta.loc[run_id]
            val = vals[bin_index - 1]
            rel = np.nan
            if np.isfinite(val) and ref_val != 0:
                rel = 100.0 * (val / ref_val - 1.0)
            rows.append({
                "run_id": run_id,
                x_col: row[x_col],
                y_col: row[y_col],
                "rel_diff_pct": rel,
            })

        hdf = pd.DataFrame(rows).dropna(subset=[x_col, y_col])
        if hdf.empty:
            messagebox.showerror("Error", "No data are available for the heatmap.")
            return

        pivot = hdf.pivot_table(index=y_col, columns=x_col, values="rel_diff_pct", aggfunc="mean")
        if pivot.empty:
            messagebox.showerror("Error", "Heatmap pivot is empty.")
            return

        fig, ax = plt.subplots(figsize=(9, 7))
        im = ax.imshow(pivot.values, aspect='auto', origin='lower')
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels([format_value(v) for v in pivot.columns])
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels([format_value(v) for v in pivot.index])
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        ax.set_title(
            f"Heatmap: bin {bin_index} [{ref_edges[bin_index-1]}, {ref_edges[bin_index]}), "
            f"100*(run/ref - 1) [%]"
        )
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("relative diff to reference [%]")
        fig.tight_layout()
        plt.show()

    def plot_all_results(self):
        self.plot_overlay_and_ratio()
        self.plot_mean_and_std()
        self.plot_bin_evolution()
        self.plot_heatmap()
        self.plot_reference_systematics_by_bin()
        self.plot_iteration_stability()


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Interactive browser for summary.tsv and stability.tsv. "
            "Allows filtering by reco/true range, by metric limits "
            "and additionally via a pair of A/B cross-filters on the same (run_id, cent)."
        )
    )
    p.add_argument("--summary", default="", help="Cesta k summary.tsv")
    p.add_argument("--stability", default="", help="Cesta k stability.tsv")
    return p.parse_args()


def main():
    args = parse_args()
    app = StabilityBrowser(summary_path=args.summary, stability_path=args.stability)
    app.mainloop()


if __name__ == "__main__":
    main()
