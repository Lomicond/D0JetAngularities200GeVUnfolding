#!/usr/bin/env python3
"""
machine_canvas_browser_v1.py

Small GUI for browsing TCanvas pages saved by Machine.C and composing selected
ROOT pads into a clean Matplotlib multi-panel figure.

Designed for ROOT files containing keys such as:
    page_0001  TCanvas
    page_0002  TCanvas
    ...

Main ideas
----------
* PyROOT is used only to read the stored TCanvas/TPad/TH1/TH2/TGraph objects.
* The original ROOT page/pad can be previewed as rendered by ROOT.
* Selected LEAF pads are re-drawn from their numerical ROOT objects with
  Matplotlib, so the final layout is editable and publication-friendly.
* Figure selections can be saved to / loaded from JSON.

Supported for Matplotlib re-drawing in v1:
    TH1 / TProfile-like 1D histograms
    TH2
    TGraph / TGraphErrors / TGraphAsymmErrors
    TMultiGraph
    THStack
    TLine

ROOT decorations such as TLatex/TLegend/TGaxis are not copied visually.
TLegend labels are, when possible, used to label Matplotlib series.

Usage
-----
    python3 machine_canvas_browser_v1.py OutputPdf/..._canvases.root

or simply:
    python3 machine_canvas_browser_v1.py

Optional:
    python3 machine_canvas_browser_v1.py --config my_figure.json

Requirements
------------
    PyROOT, numpy, matplotlib
Optional:
    Pillow  (used only for nicer resizing of ROOT preview PNGs)
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import tempfile
import traceback
from typing import Any, Dict, List, Optional, Sequence, Tuple

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

try:
    import numpy as np
except Exception as exc:
    raise SystemExit(f"numpy is required: {exc}")

try:
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
except Exception as exc:
    raise SystemExit(f"matplotlib with Tk support is required: {exc}")

try:
    import ROOT
except Exception as exc:
    raise SystemExit(
        "PyROOT could not be imported.\n"
        "Run this script in an environment where ROOT/PyROOT is available.\n"
        f"Original error: {exc}"
    )

ROOT.gROOT.SetBatch(True)

try:
    from PIL import Image, ImageTk
    HAVE_PIL = True
except Exception:
    HAVE_PIL = False


# ---------------------------------------------------------------------------
# Small ROOT helpers
# ---------------------------------------------------------------------------

def is_root_class(obj: Any, class_name: str) -> bool:
    try:
        return bool(obj and obj.InheritsFrom(class_name))
    except Exception:
        return False


def root_name(obj: Any) -> str:
    try:
        return str(obj.GetName())
    except Exception:
        return obj.__class__.__name__


def root_title(obj: Any) -> str:
    try:
        return str(obj.GetTitle())
    except Exception:
        return ""


def safe_label(obj: Any) -> str:
    title = root_title(obj).strip()
    name = root_name(obj).strip()
    if title and title not in (name, "Graph", "Histogram"):
        return title
    return name


def page_sort_key(name: str):
    m = re.fullmatch(r"page_(\d+)", name)
    if m:
        return (0, int(m.group(1)))
    return (1, name)


def root_to_mpl_text(text: str) -> str:
    """
    Lightweight conversion of common ROOT TLatex/axis syntax to Matplotlib mathtext.
    It intentionally stays conservative; panel x/y/title overrides are available
    in the GUI for final polishing.
    """
    if not text:
        return ""

    s = str(text)
    repl = {
        "#lambda": r"\lambda",
        "#eta": r"\eta",
        "#phi": r"\phi",
        "#Delta": r"\Delta",
        "#delta": r"\delta",
        "#sigma": r"\sigma",
        "#Sigma": r"\Sigma",
        "#pi": r"\pi",
        "#rho": r"\rho",
        "#alpha": r"\alpha",
        "#beta": r"\beta",
        "#gamma": r"\gamma",
        "#sqrt": r"\sqrt",
        "#pm": r"\pm",
        "#times": r"\times",
    }
    for a, b in repl.items():
        s = s.replace(a, b)

    s = s.replace("#it{", r"\mathit{")
    s = s.replace("#bf{", r"\mathbf{")
    s = s.replace("#frac{", r"\frac{")

    # ROOT sometimes uses plain braces/subscripts even outside TLatex symbols.
    if ("\\" in s or "_" in s or "^" in s) and "$" not in s:
        return f"${s}$"
    return s


def iter_primitives(pad: Any):
    try:
        prims = pad.GetListOfPrimitives()
        if not prims:
            return
        for obj in prims:
            yield obj
    except Exception:
        return


def direct_child_pads(pad: Any) -> List[Any]:
    return [obj for obj in iter_primitives(pad) if is_root_class(obj, "TPad")]


def all_nonpad_objects(pad: Any) -> List[Any]:
    return [obj for obj in iter_primitives(pad) if not is_root_class(obj, "TPad")]


def get_axis_edges(axis: Any) -> np.ndarray:
    n = int(axis.GetNbins())
    return np.asarray(
        [axis.GetBinLowEdge(i) for i in range(1, n + 1)]
        + [axis.GetBinUpEdge(n)],
        dtype=float,
    )


def graph_xy(graph: Any) -> Tuple[np.ndarray, np.ndarray]:
    n = int(graph.GetN())

    # ROOT 6 has GetPointX/GetPointY.
    if hasattr(graph, "GetPointX") and hasattr(graph, "GetPointY"):
        x = np.asarray([float(graph.GetPointX(i)) for i in range(n)], dtype=float)
        y = np.asarray([float(graph.GetPointY(i)) for i in range(n)], dtype=float)
        return x, y

    # ROOT 5 fallback.
    xs, ys = [], []
    for i in range(n):
        xv = ROOT.Double(0.0)
        yv = ROOT.Double(0.0)
        graph.GetPoint(i, xv, yv)
        xs.append(float(xv))
        ys.append(float(yv))
    return np.asarray(xs), np.asarray(ys)


def graph_errors(graph: Any):
    n = int(graph.GetN())

    if is_root_class(graph, "TGraphAsymmErrors"):
        exl = np.asarray([float(graph.GetErrorXlow(i)) for i in range(n)])
        exh = np.asarray([float(graph.GetErrorXhigh(i)) for i in range(n)])
        eyl = np.asarray([float(graph.GetErrorYlow(i)) for i in range(n)])
        eyh = np.asarray([float(graph.GetErrorYhigh(i)) for i in range(n)])
        return np.vstack([exl, exh]), np.vstack([eyl, eyh])

    if is_root_class(graph, "TGraphErrors"):
        ex = np.asarray([float(graph.GetErrorX(i)) for i in range(n)])
        ey = np.asarray([float(graph.GetErrorY(i)) for i in range(n)])
        return ex, ey

    return None, None


def legend_label_map(pad: Any) -> Dict[str, str]:
    """
    Map ROOT object names to labels found in a TLegend in the selected pad.
    """
    out: Dict[str, str] = {}
    for obj in iter_primitives(pad):
        if not is_root_class(obj, "TLegend"):
            continue
        try:
            entries = obj.GetListOfPrimitives()
            if not entries:
                continue
            for entry in entries:
                target = entry.GetObject()
                label = str(entry.GetLabel())
                if target and label:
                    out[root_name(target)] = label
        except Exception:
            pass
    return out


def root_logy(pad: Any) -> bool:
    try:
        return bool(pad.GetLogy())
    except Exception:
        return False


def root_logx(pad: Any) -> bool:
    try:
        return bool(pad.GetLogx())
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Matplotlib rendering of one selected ROOT pad
# ---------------------------------------------------------------------------

def marker_from_root(style: int) -> str:
    # Small pragmatic mapping; unknown ROOT styles fall back to "o".
    mapping = {
        1: "",
        2: "+",
        3: "*",
        4: "o",
        5: "x",
        20: "o",
        21: "s",
        22: "^",
        23: "v",
        24: "o",
        25: "s",
        26: "^",
        27: "D",
        28: "*",
        29: "*",
        30: "*",
        33: "D",
        34: "P",
    }
    return mapping.get(int(style), "o")


def linestyle_from_root(style: int) -> str:
    mapping = {
        1: "-",
        2: "--",
        3: ":",
        4: "-.",
        5: "--",
        6: "-.",
        7: ":",
        8: "--",
        9: "-.",
        10: ":",
    }
    return mapping.get(int(style), "-")


def object_axis_titles(obj: Any) -> Tuple[str, str]:
    try:
        return (
            root_to_mpl_text(str(obj.GetXaxis().GetTitle())),
            root_to_mpl_text(str(obj.GetYaxis().GetTitle())),
        )
    except Exception:
        return "", ""


def render_th1(ax, h, label: str):
    n = int(h.GetNbinsX())
    edges = get_axis_edges(h.GetXaxis())
    centers = 0.5 * (edges[:-1] + edges[1:])
    values = np.asarray([float(h.GetBinContent(i)) for i in range(1, n + 1)])
    errors = np.asarray([float(h.GetBinError(i)) for i in range(1, n + 1)])

    try:
        marker = marker_from_root(int(h.GetMarkerStyle()))
    except Exception:
        marker = "o"

    try:
        ls = linestyle_from_root(int(h.GetLineStyle()))
    except Exception:
        ls = "-"

    # ROOT histograms used here are mostly marker/error plots. If there are no
    # useful errors, a step curve is cleaner.
    finite_err = np.isfinite(errors)
    if np.any(finite_err & (errors > 0)):
        ax.errorbar(
            centers,
            values,
            yerr=errors,
            xerr=None,
            fmt=marker if marker else "none",
            linestyle=ls,
            capsize=0,
            label=label,
        )
    else:
        # "post" needs one extra y value.
        ystep = np.r_[values, values[-1] if len(values) else np.nan]
        ax.step(edges, ystep, where="post", label=label)

    return object_axis_titles(h)


def render_th2(fig, ax, h):
    nx = int(h.GetNbinsX())
    ny = int(h.GetNbinsY())
    xedges = get_axis_edges(h.GetXaxis())
    yedges = get_axis_edges(h.GetYaxis())

    z = np.empty((ny, nx), dtype=float)
    for iy in range(1, ny + 1):
        for ix in range(1, nx + 1):
            z[iy - 1, ix - 1] = float(h.GetBinContent(ix, iy))

    mesh = ax.pcolormesh(xedges, yedges, z, shading="auto")
    fig.colorbar(mesh, ax=ax, pad=0.02)
    return object_axis_titles(h)


def render_graph(ax, graph, label: str):
    x, y = graph_xy(graph)
    xerr, yerr = graph_errors(graph)

    try:
        marker = marker_from_root(int(graph.GetMarkerStyle()))
    except Exception:
        marker = "o"

    try:
        ls = linestyle_from_root(int(graph.GetLineStyle()))
    except Exception:
        ls = "-"

    if xerr is not None or yerr is not None:
        ax.errorbar(
            x,
            y,
            xerr=xerr,
            yerr=yerr,
            fmt=marker if marker else "none",
            linestyle=ls,
            capsize=0,
            label=label,
        )
    else:
        ax.plot(x, y, marker=marker if marker else None, linestyle=ls, label=label)

    # TGraph axes may exist only after drawing; guard access.
    try:
        hist = graph.GetHistogram()
        if hist:
            return object_axis_titles(hist)
    except Exception:
        pass
    return "", ""


def render_line(ax, line):
    try:
        ax.plot(
            [float(line.GetX1()), float(line.GetX2())],
            [float(line.GetY1()), float(line.GetY2())],
            linestyle=linestyle_from_root(int(line.GetLineStyle())),
        )
    except Exception:
        pass


def flatten_drawables(pad: Any) -> List[Any]:
    """
    Return direct drawable objects in a leaf pad, expanding THStack/TMultiGraph.
    """
    out: List[Any] = []
    for obj in iter_primitives(pad):
        if is_root_class(obj, "TPad"):
            continue

        if is_root_class(obj, "THStack"):
            try:
                hs = obj.GetHists()
                if hs:
                    out.extend(list(hs))
            except Exception:
                pass
            continue

        if is_root_class(obj, "TMultiGraph"):
            try:
                gs = obj.GetListOfGraphs()
                if gs:
                    out.extend(list(gs))
            except Exception:
                pass
            continue

        out.append(obj)
    return out


def render_pad_to_axis(
    fig,
    ax,
    pad: Any,
    panel_cfg: Dict[str, Any],
):
    objects = flatten_drawables(pad)
    labels = legend_label_map(pad)

    default_xlabel = ""
    default_ylabel = ""
    have_data = False
    have_labeled_artist = False

    # TH2 first: a typical ROOT pad should have one main TH2.
    th2s = [o for o in objects if is_root_class(o, "TH2")]
    if th2s:
        xlabel, ylabel = render_th2(fig, ax, th2s[0])
        default_xlabel = default_xlabel or xlabel
        default_ylabel = default_ylabel or ylabel
        have_data = True

    # TH1 but not TH2.
    for obj in objects:
        if is_root_class(obj, "TH2"):
            continue
        if is_root_class(obj, "TH1"):
            label = labels.get(root_name(obj), safe_label(obj))
            xlabel, ylabel = render_th1(ax, obj, label)
            default_xlabel = default_xlabel or xlabel
            default_ylabel = default_ylabel or ylabel
            have_data = True
            if label and not label.startswith("h"):
                have_labeled_artist = True

    for obj in objects:
        if is_root_class(obj, "TGraph"):
            label = labels.get(root_name(obj), safe_label(obj))
            xlabel, ylabel = render_graph(ax, obj, label)
            default_xlabel = default_xlabel or xlabel
            default_ylabel = default_ylabel or ylabel
            have_data = True
            if label and not label.startswith("Graph"):
                have_labeled_artist = True
        elif is_root_class(obj, "TLine"):
            render_line(ax, obj)

    if not have_data:
        ax.text(
            0.5, 0.5,
            "No supported TH1/TH2/TGraph objects\nfound in this pad.",
            ha="center", va="center", transform=ax.transAxes,
        )

    title = panel_cfg.get("title", "").strip()
    xlabel = panel_cfg.get("xlabel", "").strip() or default_xlabel
    ylabel = panel_cfg.get("ylabel", "").strip() or default_ylabel

    if title:
        ax.set_title(root_to_mpl_text(title))
    if xlabel:
        ax.set_xlabel(root_to_mpl_text(xlabel))
    if ylabel:
        ax.set_ylabel(root_to_mpl_text(ylabel))

    logy_mode = panel_cfg.get("logy", "root")
    if logy_mode == "yes" or (logy_mode == "root" and root_logy(pad)):
        ax.set_yscale("log")

    logx_mode = panel_cfg.get("logx", "root")
    if logx_mode == "yes" or (logx_mode == "root" and root_logx(pad)):
        ax.set_xscale("log")

    if panel_cfg.get("legend", True) and have_labeled_artist:
        ax.legend(frameon=False)

    label = panel_cfg.get("panel_label", "").strip()
    if label:
        ax.text(
            0.02, 0.98, label,
            ha="left", va="top",
            transform=ax.transAxes,
            fontweight="bold",
        )

    ax.tick_params(direction="in", top=True, right=True)


# ---------------------------------------------------------------------------
# ROOT file model
# ---------------------------------------------------------------------------

class RootCanvasFile:
    def __init__(self):
        self.path: Optional[Path] = None
        self.file = None
        self.canvases: Dict[str, Any] = {}

    def close(self):
        self.canvases.clear()
        if self.file:
            try:
                self.file.Close()
            except Exception:
                pass
        self.file = None
        self.path = None

    def open(self, path: str):
        self.close()
        p = Path(path).expanduser().resolve()
        f = ROOT.TFile.Open(str(p), "READ")
        if not f or f.IsZombie():
            raise RuntimeError(f"Cannot open ROOT file: {p}")

        self.path = p
        self.file = f

        keys = []
        for key in f.GetListOfKeys():
            try:
                class_name = str(key.GetClassName())
                if "TCanvas" not in class_name:
                    continue
                keys.append(str(key.GetName()))
            except Exception:
                continue

        for name in sorted(keys, key=page_sort_key):
            obj = f.Get(name)
            if obj and is_root_class(obj, "TCanvas"):
                self.canvases[name] = obj

        if not self.canvases:
            raise RuntimeError(
                f"No TCanvas keys found in {p}.\n"
                "Expected objects such as page_0001, page_0002, ..."
            )

    def pages(self) -> List[str]:
        return sorted(self.canvases.keys(), key=page_sort_key)

    def get_page(self, page: str):
        return self.canvases[page]

    def get_pad(self, page: str, pad_path: Sequence[str]):
        obj = self.get_page(page)
        for name in pad_path:
            found = None
            for child in direct_child_pads(obj):
                if root_name(child) == name:
                    found = child
                    break
            if found is None:
                raise KeyError(f"Pad path not found: {page}/{'/'.join(pad_path)}")
            obj = found
        return obj


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------

class CanvasBrowserApp(tk.Tk):
    def __init__(self, initial_root: Optional[str] = None, initial_config: Optional[str] = None):
        super().__init__()
        self.title("Machine canvas browser / figure builder")
        self.geometry("1600x950")
        self.minsize(1150, 700)

        self.model = RootCanvasFile()
        self.tree_meta: Dict[str, Dict[str, Any]] = {}
        self.panels: List[Dict[str, Any]] = []

        self.preview_photo = None
        self.preview_temp_dir = Path(tempfile.mkdtemp(prefix="machine_canvas_browser_"))
        self.mpl_canvas = None
        self.current_figure = None

        self._build_ui()
        self.protocol("WM_DELETE_WINDOW", self._on_close)

        if initial_config:
            self.after(100, lambda: self.load_config(initial_config))
        elif initial_root:
            self.after(100, lambda: self.open_root(initial_root))

    # ---------------- UI ----------------

    def _build_ui(self):
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)

        top = ttk.Frame(self, padding=6)
        top.grid(row=0, column=0, sticky="ew")
        top.columnconfigure(2, weight=1)

        ttk.Button(top, text="Open ROOT", command=self.choose_root).grid(row=0, column=0, padx=(0, 6))
        ttk.Button(top, text="Load figure config", command=self.choose_load_config).grid(row=0, column=1, padx=(0, 10))

        self.file_var = tk.StringVar(value="No ROOT file loaded")
        ttk.Label(top, textvariable=self.file_var).grid(row=0, column=2, sticky="w")

        ttk.Label(top, text="Go to PDF page:").grid(row=0, column=3, padx=(10, 4))
        self.goto_page_var = tk.StringVar()
        goto = ttk.Entry(top, textvariable=self.goto_page_var, width=8)
        goto.grid(row=0, column=4)
        goto.bind("<Return>", lambda _e: self.go_to_page())
        ttk.Button(top, text="Go", command=self.go_to_page).grid(row=0, column=5, padx=(4, 0))

        main = ttk.Panedwindow(self, orient=tk.HORIZONTAL)
        main.grid(row=1, column=0, sticky="nsew", padx=6, pady=(0, 6))

        # Left browser
        left = ttk.Frame(main, padding=4)
        main.add(left, weight=1)
        left.rowconfigure(1, weight=1)
        left.columnconfigure(0, weight=1)

        ttk.Label(left, text="ROOT pages / pads").grid(row=0, column=0, sticky="w")

        self.tree = ttk.Treeview(left, columns=("kind",), show="tree headings", selectmode="browse")
        self.tree.heading("#0", text="Page / pad")
        self.tree.heading("kind", text="Type")
        self.tree.column("#0", width=240, stretch=True)
        self.tree.column("kind", width=90, stretch=False)
        self.tree.grid(row=1, column=0, sticky="nsew")
        sb = ttk.Scrollbar(left, orient="vertical", command=self.tree.yview)
        sb.grid(row=1, column=1, sticky="ns")
        self.tree.configure(yscrollcommand=sb.set)
        self.tree.bind("<<TreeviewSelect>>", self.on_tree_select)

        ttk.Button(left, text="Add selected leaf pad → figure", command=self.add_selected_pad).grid(
            row=2, column=0, columnspan=2, sticky="ew", pady=(6, 4)
        )

        ttk.Label(left, text="Objects in selected pad").grid(row=3, column=0, sticky="w", pady=(8, 2))
        self.object_list = tk.Listbox(left, height=10)
        self.object_list.grid(row=4, column=0, columnspan=2, sticky="nsew")
        left.rowconfigure(4, weight=0)

        # Middle preview
        middle = ttk.Frame(main, padding=4)
        main.add(middle, weight=2)
        middle.rowconfigure(1, weight=1)
        middle.columnconfigure(0, weight=1)

        self.preview_title_var = tk.StringVar(value="ROOT preview")
        ttk.Label(middle, textvariable=self.preview_title_var).grid(row=0, column=0, sticky="w")

        self.preview_frame = ttk.Frame(middle, relief="sunken")
        self.preview_frame.grid(row=1, column=0, sticky="nsew")
        self.preview_frame.rowconfigure(0, weight=1)
        self.preview_frame.columnconfigure(0, weight=1)

        self.preview_label = ttk.Label(self.preview_frame, text="Open a ROOT file.")
        self.preview_label.grid(row=0, column=0, sticky="nsew")

        # Right figure builder
        right = ttk.Frame(main, padding=4)
        main.add(right, weight=2)
        right.columnconfigure(0, weight=1)
        right.rowconfigure(1, weight=1)

        ttk.Label(right, text="Selected figure panels").grid(row=0, column=0, sticky="w")

        self.panel_list = tk.Listbox(right, height=10, exportselection=False)
        self.panel_list.grid(row=1, column=0, sticky="nsew")
        self.panel_list.bind("<<ListboxSelect>>", self.on_panel_select)

        buttons = ttk.Frame(right)
        buttons.grid(row=2, column=0, sticky="ew", pady=4)
        ttk.Button(buttons, text="↑", width=4, command=lambda: self.move_panel(-1)).pack(side="left")
        ttk.Button(buttons, text="↓", width=4, command=lambda: self.move_panel(+1)).pack(side="left", padx=2)
        ttk.Button(buttons, text="Remove", command=self.remove_panel).pack(side="left", padx=(8, 0))
        ttk.Button(buttons, text="Clear", command=self.clear_panels).pack(side="left", padx=2)

        # Panel properties
        props = ttk.LabelFrame(right, text="Selected panel settings", padding=6)
        props.grid(row=3, column=0, sticky="ew", pady=(4, 6))
        props.columnconfigure(1, weight=1)

        self.panel_label_var = tk.StringVar()
        self.panel_title_var = tk.StringVar()
        self.panel_xlabel_var = tk.StringVar()
        self.panel_ylabel_var = tk.StringVar()
        self.panel_logy_var = tk.StringVar(value="root")
        self.panel_logx_var = tk.StringVar(value="root")
        self.panel_legend_var = tk.BooleanVar(value=True)

        fields = [
            ("Panel label", self.panel_label_var),
            ("Title override", self.panel_title_var),
            ("X label override", self.panel_xlabel_var),
            ("Y label override", self.panel_ylabel_var),
        ]
        for r, (label, var) in enumerate(fields):
            ttk.Label(props, text=label + ":").grid(row=r, column=0, sticky="w", padx=(0, 4))
            ttk.Entry(props, textvariable=var).grid(row=r, column=1, sticky="ew")

        ttk.Label(props, text="Log Y:").grid(row=4, column=0, sticky="w")
        ttk.Combobox(
            props, textvariable=self.panel_logy_var,
            values=("root", "yes", "no"), state="readonly", width=8
        ).grid(row=4, column=1, sticky="w")

        ttk.Label(props, text="Log X:").grid(row=5, column=0, sticky="w")
        ttk.Combobox(
            props, textvariable=self.panel_logx_var,
            values=("root", "yes", "no"), state="readonly", width=8
        ).grid(row=5, column=1, sticky="w")

        ttk.Checkbutton(props, text="Legend", variable=self.panel_legend_var).grid(
            row=6, column=1, sticky="w"
        )
        ttk.Button(props, text="Apply panel settings", command=self.apply_panel_settings).grid(
            row=7, column=0, columnspan=2, sticky="ew", pady=(5, 0)
        )

        # Layout settings
        layout = ttk.LabelFrame(right, text="Figure layout", padding=6)
        layout.grid(row=4, column=0, sticky="ew")
        for c in range(4):
            layout.columnconfigure(c, weight=1)

        self.ncols_var = tk.IntVar(value=2)
        self.figw_var = tk.DoubleVar(value=9.0)
        self.figh_var = tk.DoubleVar(value=7.0)
        self.sharex_var = tk.BooleanVar(value=False)
        self.sharey_var = tk.BooleanVar(value=False)

        ttk.Label(layout, text="Columns").grid(row=0, column=0)
        ttk.Spinbox(layout, from_=1, to=8, textvariable=self.ncols_var, width=5).grid(row=0, column=1)
        ttk.Label(layout, text="Width").grid(row=0, column=2)
        ttk.Entry(layout, textvariable=self.figw_var, width=7).grid(row=0, column=3)

        ttk.Label(layout, text="Height").grid(row=1, column=0)
        ttk.Entry(layout, textvariable=self.figh_var, width=7).grid(row=1, column=1)
        ttk.Checkbutton(layout, text="Shared X", variable=self.sharex_var).grid(row=1, column=2)
        ttk.Checkbutton(layout, text="Shared Y", variable=self.sharey_var).grid(row=1, column=3)

        outbuttons = ttk.Frame(right)
        outbuttons.grid(row=5, column=0, sticky="ew", pady=(6, 0))
        ttk.Button(outbuttons, text="Preview figure", command=self.preview_figure).pack(side="left")
        ttk.Button(outbuttons, text="Save PDF", command=lambda: self.save_figure("pdf")).pack(side="left", padx=3)
        ttk.Button(outbuttons, text="Save PNG", command=lambda: self.save_figure("png")).pack(side="left")
        ttk.Button(outbuttons, text="Save config", command=self.choose_save_config).pack(side="right")

        # Bottom Matplotlib preview spanning window
        mplbox = ttk.LabelFrame(self, text="Matplotlib figure preview", padding=4)
        mplbox.grid(row=2, column=0, sticky="nsew", padx=6, pady=(0, 6))
        mplbox.columnconfigure(0, weight=1)
        mplbox.rowconfigure(0, weight=1)
        self.mpl_frame = mplbox

        self.rowconfigure(2, weight=1)

    # ---------------- ROOT loading ----------------

    def choose_root(self):
        path = filedialog.askopenfilename(
            title="Open Machine canvas ROOT file",
            filetypes=[("ROOT files", "*.root"), ("All files", "*.*")],
        )
        if path:
            self.open_root(path)

    def open_root(self, path: str):
        try:
            self.model.open(path)
            self.file_var.set(str(self.model.path))
            self.populate_tree()
            pages = self.model.pages()
            if pages:
                iid = self._iid_for_page(pages[0])
                if self.tree.exists(iid):
                    self.tree.selection_set(iid)
                    self.tree.focus(iid)
                    self.tree.see(iid)
                    self.on_tree_select()
        except Exception as exc:
            self._show_exception("Cannot open ROOT file", exc)

    def _iid_for_page(self, page: str) -> str:
        return f"page::{page}"

    def populate_tree(self):
        self.tree.delete(*self.tree.get_children())
        self.tree_meta.clear()

        for page in self.model.pages():
            iid = self._iid_for_page(page)
            self.tree.insert("", "end", iid=iid, text=page, values=("TCanvas",), open=False)
            self.tree_meta[iid] = {"page": page, "pad_path": []}
            self._insert_child_pads(iid, self.model.get_page(page), page, [])

    def _insert_child_pads(self, parent_iid: str, parent_pad: Any, page: str, path: List[str]):
        counts: Dict[str, int] = {}
        for child in direct_child_pads(parent_pad):
            name = root_name(child)
            counts[name] = counts.get(name, 0) + 1
            display = name if counts[name] == 1 else f"{name} [{counts[name]}]"

            child_path = path + [name]
            iid = f"{parent_iid}::{len(self.tree.get_children(parent_iid))}::{name}"
            self.tree.insert(parent_iid, "end", iid=iid, text=display, values=(child.ClassName(),))
            self.tree_meta[iid] = {"page": page, "pad_path": child_path}
            self._insert_child_pads(iid, child, page, child_path)

    def go_to_page(self):
        raw = self.goto_page_var.get().strip()
        if not raw:
            return
        try:
            num = int(raw)
        except ValueError:
            messagebox.showwarning("Page", "Enter a numeric PDF page, e.g. 63.")
            return

        page = f"page_{num:04d}"
        iid = self._iid_for_page(page)
        if not self.tree.exists(iid):
            messagebox.showwarning("Page", f"{page} was not found in the ROOT file.")
            return

        self.tree.selection_set(iid)
        self.tree.focus(iid)
        self.tree.see(iid)
        self.on_tree_select()

    # ---------------- Selection / ROOT preview ----------------

    def selected_tree_meta(self) -> Optional[Dict[str, Any]]:
        sel = self.tree.selection()
        if not sel:
            return None
        return self.tree_meta.get(sel[0])

    def on_tree_select(self, *_args):
        meta = self.selected_tree_meta()
        if not meta:
            return

        try:
            pad = self.model.get_pad(meta["page"], meta["pad_path"])
            label = meta["page"]
            if meta["pad_path"]:
                label += " / " + " / ".join(meta["pad_path"])
            self.preview_title_var.set(label)
            self.populate_object_list(pad)
            self.render_root_preview(pad, label)
        except Exception as exc:
            self._show_exception("Preview failed", exc)

    def populate_object_list(self, pad: Any):
        self.object_list.delete(0, tk.END)
        for obj in iter_primitives(pad):
            self.object_list.insert(
                tk.END,
                f"{obj.ClassName():<22}  {root_name(obj)}"
            )

    def render_root_preview(self, pad: Any, label: str):
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", label)
        png = self.preview_temp_dir / f"{safe}.png"

        try:
            pad.Modified()
            pad.Update()
        except Exception:
            pass

        # TCanvas and TPad both inherit from TVirtualPad and can Print().
        pad.Print(str(png))

        self.update_idletasks()
        max_w = max(300, self.preview_frame.winfo_width() - 10)
        max_h = max(250, self.preview_frame.winfo_height() - 10)

        if HAVE_PIL:
            img = Image.open(png)
            img.thumbnail((max_w, max_h))
            self.preview_photo = ImageTk.PhotoImage(img)
        else:
            self.preview_photo = tk.PhotoImage(file=str(png))

        self.preview_label.configure(image=self.preview_photo, text="")

    # ---------------- Panel builder ----------------

    def add_selected_pad(self):
        meta = self.selected_tree_meta()
        if not meta:
            return

        try:
            pad = self.model.get_pad(meta["page"], meta["pad_path"])
        except Exception as exc:
            self._show_exception("Cannot add pad", exc)
            return

        children = direct_child_pads(pad)
        if children:
            messagebox.showinfo(
                "Select a leaf pad",
                "This selection contains nested TPad objects.\n\n"
                "For an editable Matplotlib panel, select one of its leaf pads "
                "in the tree instead of the parent canvas/pad."
            )
            return

        panel = {
            "page": meta["page"],
            "pad_path": list(meta["pad_path"]),
            "panel_label": self._next_panel_label(),
            "title": "",
            "xlabel": "",
            "ylabel": "",
            "logy": "root",
            "logx": "root",
            "legend": True,
        }
        self.panels.append(panel)
        self.refresh_panel_list(select=len(self.panels) - 1)

    def _next_panel_label(self) -> str:
        i = len(self.panels)
        if i < 26:
            return f"({chr(ord('a') + i)})"
        return f"({i + 1})"

    def panel_display_name(self, p: Dict[str, Any]) -> str:
        path = " / ".join(p.get("pad_path", []))
        return p["page"] + (f" / {path}" if path else "")

    def refresh_panel_list(self, select: Optional[int] = None):
        self.panel_list.delete(0, tk.END)
        for i, p in enumerate(self.panels, start=1):
            self.panel_list.insert(tk.END, f"{i:02d}. {self.panel_display_name(p)}")

        if select is not None and 0 <= select < len(self.panels):
            self.panel_list.selection_set(select)
            self.panel_list.activate(select)
            self.on_panel_select()

    def selected_panel_index(self) -> Optional[int]:
        sel = self.panel_list.curselection()
        if not sel:
            return None
        return int(sel[0])

    def on_panel_select(self, *_args):
        idx = self.selected_panel_index()
        if idx is None:
            return
        p = self.panels[idx]
        self.panel_label_var.set(p.get("panel_label", ""))
        self.panel_title_var.set(p.get("title", ""))
        self.panel_xlabel_var.set(p.get("xlabel", ""))
        self.panel_ylabel_var.set(p.get("ylabel", ""))
        self.panel_logy_var.set(p.get("logy", "root"))
        self.panel_logx_var.set(p.get("logx", "root"))
        self.panel_legend_var.set(bool(p.get("legend", True)))

    def apply_panel_settings(self):
        idx = self.selected_panel_index()
        if idx is None:
            return
        p = self.panels[idx]
        p["panel_label"] = self.panel_label_var.get()
        p["title"] = self.panel_title_var.get()
        p["xlabel"] = self.panel_xlabel_var.get()
        p["ylabel"] = self.panel_ylabel_var.get()
        p["logy"] = self.panel_logy_var.get()
        p["logx"] = self.panel_logx_var.get()
        p["legend"] = bool(self.panel_legend_var.get())

    def move_panel(self, delta: int):
        idx = self.selected_panel_index()
        if idx is None:
            return
        j = idx + delta
        if j < 0 or j >= len(self.panels):
            return
        self.panels[idx], self.panels[j] = self.panels[j], self.panels[idx]
        self.refresh_panel_list(select=j)

    def remove_panel(self):
        idx = self.selected_panel_index()
        if idx is None:
            return
        self.panels.pop(idx)
        self.refresh_panel_list(select=min(idx, len(self.panels) - 1) if self.panels else None)

    def clear_panels(self):
        self.panels.clear()
        self.refresh_panel_list()

    # ---------------- Matplotlib figure ----------------

    def build_figure(self):
        if not self.panels:
            raise RuntimeError("No panels selected.")

        self.apply_panel_settings()

        n = len(self.panels)
        ncols = max(1, int(self.ncols_var.get()))
        ncols = min(ncols, n)
        nrows = int(math.ceil(n / ncols))

        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(float(self.figw_var.get()), float(self.figh_var.get())),
            squeeze=False,
            sharex=bool(self.sharex_var.get()),
            sharey=bool(self.sharey_var.get()),
        )

        flat_axes = axes.ravel()

        for i, p in enumerate(self.panels):
            pad = self.model.get_pad(p["page"], p.get("pad_path", []))
            render_pad_to_axis(fig, flat_axes[i], pad, p)

        for ax in flat_axes[n:]:
            ax.set_visible(False)

        fig.tight_layout()
        return fig

    def preview_figure(self):
        try:
            fig = self.build_figure()
        except Exception as exc:
            self._show_exception("Cannot build figure", exc)
            return

        if self.mpl_canvas:
            try:
                self.mpl_canvas.get_tk_widget().destroy()
            except Exception:
                pass

        if self.current_figure is not None:
            try:
                plt.close(self.current_figure)
            except Exception:
                pass

        self.current_figure = fig
        self.mpl_canvas = FigureCanvasTkAgg(fig, master=self.mpl_frame)
        self.mpl_canvas.draw()
        self.mpl_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

    def save_figure(self, ext: str):
        if not self.panels:
            messagebox.showwarning("Figure", "No panels selected.")
            return

        path = filedialog.asksaveasfilename(
            title=f"Save figure as {ext.upper()}",
            defaultextension=f".{ext}",
            filetypes=[(ext.upper(), f"*.{ext}"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            fig = self.build_figure()
            fig.savefig(path, bbox_inches="tight")
            plt.close(fig)
            messagebox.showinfo("Saved", f"Saved:\n{path}")
        except Exception as exc:
            self._show_exception("Cannot save figure", exc)

    # ---------------- JSON config ----------------

    def current_config(self) -> Dict[str, Any]:
        self.apply_panel_settings()
        return {
            "version": 1,
            "input_root": str(self.model.path) if self.model.path else "",
            "layout": {
                "ncols": int(self.ncols_var.get()),
                "fig_width": float(self.figw_var.get()),
                "fig_height": float(self.figh_var.get()),
                "sharex": bool(self.sharex_var.get()),
                "sharey": bool(self.sharey_var.get()),
            },
            "panels": self.panels,
        }

    def choose_save_config(self):
        if not self.model.path:
            messagebox.showwarning("Config", "Open a ROOT file first.")
            return
        path = filedialog.asksaveasfilename(
            title="Save figure config",
            defaultextension=".json",
            filetypes=[("JSON", "*.json"), ("All files", "*.*")],
        )
        if path:
            self.save_config(path)

    def save_config(self, path: str):
        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(self.current_config(), f, indent=2)
                f.write("\n")
            messagebox.showinfo("Saved", f"Saved figure config:\n{path}")
        except Exception as exc:
            self._show_exception("Cannot save config", exc)

    def choose_load_config(self):
        path = filedialog.askopenfilename(
            title="Load figure config",
            filetypes=[("JSON", "*.json"), ("All files", "*.*")],
        )
        if path:
            self.load_config(path)

    def load_config(self, path: str):
        try:
            with open(path, "r", encoding="utf-8") as f:
                cfg = json.load(f)

            root_path = cfg.get("input_root", "")
            if not root_path:
                raise RuntimeError("Config has no input_root.")

            # If a relative ROOT path was saved, resolve it relative to the JSON.
            rp = Path(root_path)
            if not rp.is_absolute():
                rp = (Path(path).resolve().parent / rp).resolve()

            self.open_root(str(rp))

            layout = cfg.get("layout", {})
            self.ncols_var.set(int(layout.get("ncols", 2)))
            self.figw_var.set(float(layout.get("fig_width", 9.0)))
            self.figh_var.set(float(layout.get("fig_height", 7.0)))
            self.sharex_var.set(bool(layout.get("sharex", False)))
            self.sharey_var.set(bool(layout.get("sharey", False)))

            panels = cfg.get("panels", [])
            validated = []
            for p in panels:
                # Validate locator now, so stale configs fail clearly.
                self.model.get_pad(p["page"], p.get("pad_path", []))
                q = {
                    "page": p["page"],
                    "pad_path": list(p.get("pad_path", [])),
                    "panel_label": p.get("panel_label", ""),
                    "title": p.get("title", ""),
                    "xlabel": p.get("xlabel", ""),
                    "ylabel": p.get("ylabel", ""),
                    "logy": p.get("logy", "root"),
                    "logx": p.get("logx", "root"),
                    "legend": bool(p.get("legend", True)),
                }
                validated.append(q)

            self.panels = validated
            self.refresh_panel_list(select=0 if self.panels else None)

        except Exception as exc:
            self._show_exception("Cannot load config", exc)

    # ---------------- misc ----------------

    def _show_exception(self, title: str, exc: Exception):
        traceback.print_exc()
        messagebox.showerror(title, f"{exc}")

    def _on_close(self):
        try:
            self.model.close()
        except Exception:
            pass
        try:
            if self.current_figure is not None:
                plt.close(self.current_figure)
        except Exception:
            pass
        self.destroy()


def parse_args():
    ap = argparse.ArgumentParser(
        description="Browse Machine.C TCanvas pages and build Matplotlib figures."
    )
    ap.add_argument(
        "root_file", nargs="?",
        help="ROOT file containing page_XXXX TCanvas objects."
    )
    ap.add_argument(
        "--config",
        help="Load a previously saved figure-builder JSON config."
    )
    return ap.parse_args()


def main():
    args = parse_args()
    app = CanvasBrowserApp(
        initial_root=args.root_file,
        initial_config=args.config,
    )
    app.mainloop()


if __name__ == "__main__":
    main()
