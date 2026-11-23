#!/usr/bin/env python3
"""
machine_canvas_browser_v4.py

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
* Per-panel axis limits and series labels are editable in the GUI.
* Points can be connected/disconnected per panel.
* Free text can be placed in any of four panel corners.
* Multiple selected panels can receive shared settings/legend labels at once.
* Panels can be assigned to multiple output PDF pages.

Supported for Matplotlib re-drawing in v4:
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
    python3 machine_canvas_browser_v4.py OutputPdf/..._canvases.root

or simply:
    python3 machine_canvas_browser_v4.py

Optional:
    python3 machine_canvas_browser_v4.py --config my_figure.json

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
import subprocess
import sys
import tempfile
import traceback
from typing import Any, Dict, List, Optional, Sequence, Tuple

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

try:
    import numpy as np
except Exception as exc:
    raise SystemExit(f"numpy is required: {exc}")

HELPER_MODE = "--_render-preview" in sys.argv

try:
    import matplotlib
    matplotlib.use("Agg" if HELPER_MODE else "TkAgg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.backends.backend_pdf import PdfPages
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


def legend_label_map(pad: Any, recursive: bool = False) -> Dict[str, str]:
    """
    Map ROOT object names to labels found in TLegend objects.

    In Machine.C the legend is often drawn on the parent canvas while the
    histogram itself lives in a nested TPad.  v1/v2 only searched the leaf pad,
    which is why labels such as "Unfold dummy1" could appear.  When recursive
    is True we search the complete canvas/pad hierarchy.
    """
    out: Dict[str, str] = {}

    def walk(p):
        for obj in iter_primitives(p):
            if is_root_class(obj, "TLegend"):
                try:
                    entries = obj.GetListOfPrimitives()
                    if entries:
                        for entry in entries:
                            target = entry.GetObject()
                            label = str(entry.GetLabel()).strip()
                            if target and label:
                                out[root_name(target)] = label
                except Exception:
                    pass

            if recursive and is_root_class(obj, "TPad"):
                walk(obj)

    walk(pad)
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


def pad_abs_ndc_bbox(pad: Any) -> Optional[Tuple[float, float, float, float]]:
    """
    Absolute pad rectangle in canvas NDC coordinates:
        (x_low, y_low, width, height)

    Reading these geometry values does not repaint the ROOT canvas.
    """
    try:
        x = float(pad.GetAbsXlowNDC())
        y = float(pad.GetAbsYlowNDC())
        w = float(pad.GetAbsWNDC())
        h = float(pad.GetAbsHNDC())
        vals = (x, y, w, h)
        if not all(np.isfinite(v) for v in vals):
            return None
        if w <= 0 or h <= 0:
            return None
        return vals
    except Exception:
        return None


def render_page_preview_in_child(root_path: str, page: str, output_png: str) -> int:
    """
    Internal helper mode.

    ROOT's painting of a deserialized TCanvas can occasionally segfault inside
    C++.  Therefore the GUI never calls TCanvas::Print itself.  This function
    is executed in a short-lived subprocess; if ROOT crashes, only that helper
    process dies and the GUI remains alive.
    """
    f = ROOT.TFile.Open(str(root_path), "READ")
    if not f or f.IsZombie():
        print(f"[preview-helper] Cannot open ROOT file: {root_path}", file=sys.stderr, flush=True)
        return 2

    can = f.Get(page)
    if not can or not is_root_class(can, "TCanvas"):
        print(f"[preview-helper] TCanvas not found: {page}", file=sys.stderr, flush=True)
        try:
            f.Close()
        except Exception:
            pass
        return 3

    try:
        can.Modified()
        can.Update()
        # Flush a marker before ROOT painting: useful if ROOT dies afterwards.
        print(f"[preview-helper] printing {page}", flush=True)
        can.Print(str(output_png))
        print(f"[preview-helper] done {page}", flush=True)
    finally:
        try:
            f.Close()
        except Exception:
            pass

    return 0


def optional_float(value: Any) -> Optional[float]:
    """Blank/None -> None, otherwise finite float."""
    if value is None:
        return None
    s = str(value).strip()
    if not s:
        return None
    v = float(s)
    if not np.isfinite(v):
        raise ValueError(f"Axis limit is not finite: {value}")
    return v


def drawable_series_objects(pad: Any) -> List[Any]:
    """
    Ordered TH1/TGraph objects that become legend-able 1D series.
    TH2 objects are excluded. THStack/TMultiGraph are expanded.
    """
    out: List[Any] = []
    for obj in flatten_drawables(pad):
        if is_root_class(obj, "TH2"):
            continue
        if is_root_class(obj, "TH1") or is_root_class(obj, "TGraph"):
            out.append(obj)
    return out


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


def render_th1(ax, h, label: str, connect_points: bool = True):
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
            linestyle=ls if connect_points else "none",
            capsize=0,
            label=label,
        )
    else:
        if connect_points:
            # "post" needs one extra y value.
            ystep = np.r_[values, values[-1] if len(values) else np.nan]
            ax.step(edges, ystep, where="post", label=label)
        else:
            ax.plot(
                centers,
                values,
                marker=marker if marker else "o",
                linestyle="none",
                label=label,
            )

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


def render_graph(ax, graph, label: str, connect_points: bool = True):
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
            linestyle=ls if connect_points else "none",
            capsize=0,
            label=label,
        )
    else:
        ax.plot(
            x,
            y,
            marker=marker if marker else None,
            linestyle=ls if connect_points else "none",
            label=label,
        )

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
    page_legend_labels: Optional[Dict[str, str]] = None,
):
    objects = flatten_drawables(pad)

    # First use legends found anywhere on the page, then let a legend drawn
    # directly in the leaf pad override them, and finally let GUI edits win.
    labels: Dict[str, str] = {}
    if page_legend_labels:
        labels.update(page_legend_labels)
    labels.update(legend_label_map(pad, recursive=False))
    custom_labels = panel_cfg.get("series_labels", {}) or {}

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
            obj_name = root_name(obj)
            label = custom_labels.get(obj_name, labels.get(obj_name, safe_label(obj)))
            xlabel, ylabel = render_th1(
                ax, obj, label,
                connect_points=bool(panel_cfg.get("connect_points", True))
            )
            default_xlabel = default_xlabel or xlabel
            default_ylabel = default_ylabel or ylabel
            have_data = True
            if label and not label.startswith("h"):
                have_labeled_artist = True

    for obj in objects:
        if is_root_class(obj, "TGraph"):
            obj_name = root_name(obj)
            label = custom_labels.get(obj_name, labels.get(obj_name, safe_label(obj)))
            xlabel, ylabel = render_graph(
                ax, obj, label,
                connect_points=bool(panel_cfg.get("connect_points", True))
            )
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

    # Explicit per-panel axis limits. Blank means Matplotlib/ROOT-derived auto.
    xmin = optional_float(panel_cfg.get("xmin", ""))
    xmax = optional_float(panel_cfg.get("xmax", ""))
    ymin = optional_float(panel_cfg.get("ymin", ""))
    ymax = optional_float(panel_cfg.get("ymax", ""))

    if xmin is not None or xmax is not None:
        current = ax.get_xlim()
        ax.set_xlim(
            left=xmin if xmin is not None else current[0],
            right=xmax if xmax is not None else current[1],
        )

    if ymin is not None or ymax is not None:
        current = ax.get_ylim()
        ax.set_ylim(
            bottom=ymin if ymin is not None else current[0],
            top=ymax if ymax is not None else current[1],
        )

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

    annotation = panel_cfg.get("annotation_text", "").strip()
    if annotation:
        pos = panel_cfg.get("annotation_position", "upper left")
        positions = {
            "upper left":  (0.04, 0.90, "left",  "top"),
            "upper right": (0.96, 0.90, "right", "top"),
            "lower left":  (0.04, 0.06, "left",  "bottom"),
            "lower right": (0.96, 0.06, "right", "bottom"),
        }
        x, y, ha, va = positions.get(pos, positions["upper left"])
        ax.text(
            x, y,
            root_to_mpl_text(annotation),
            ha=ha, va=va,
            transform=ax.transAxes,
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

        ttk.Label(
            right,
            text="Selected figure panels   (Ctrl/Shift = multi-select)"
        ).grid(row=0, column=0, sticky="w")

        self.panel_list = tk.Listbox(
            right,
            height=10,
            exportselection=False,
            selectmode=tk.EXTENDED,
        )
        self.panel_list.grid(row=1, column=0, sticky="nsew")
        self.panel_list.bind("<<ListboxSelect>>", self.on_panel_select)

        buttons = ttk.Frame(right)
        buttons.grid(row=2, column=0, sticky="ew", pady=4)
        ttk.Button(buttons, text="↑", width=4, command=lambda: self.move_panel(-1)).pack(side="left")
        ttk.Button(buttons, text="↓", width=4, command=lambda: self.move_panel(+1)).pack(side="left", padx=2)
        ttk.Button(buttons, text="Remove", command=self.remove_panel).pack(side="left", padx=(8, 0))
        ttk.Button(buttons, text="Select all", command=self.select_all_panels).pack(side="left", padx=2)
        ttk.Button(buttons, text="Clear", command=self.clear_panels).pack(side="left", padx=2)

        # Panel properties
        props = ttk.LabelFrame(right, text="Selected panel settings", padding=6)
        props.grid(row=3, column=0, sticky="ew", pady=(4, 6))
        props.columnconfigure(1, weight=1)

        self.panel_page_var = tk.IntVar(value=1)
        self.panel_label_var = tk.StringVar()
        self.panel_title_var = tk.StringVar()
        self.panel_xlabel_var = tk.StringVar()
        self.panel_ylabel_var = tk.StringVar()
        self.panel_xmin_var = tk.StringVar()
        self.panel_xmax_var = tk.StringVar()
        self.panel_ymin_var = tk.StringVar()
        self.panel_ymax_var = tk.StringVar()
        self.panel_logy_var = tk.StringVar(value="root")
        self.panel_logx_var = tk.StringVar(value="root")
        self.panel_legend_var = tk.BooleanVar(value=True)
        self.panel_connect_var = tk.BooleanVar(value=True)
        self.panel_annotation_text_var = tk.StringVar()
        self.panel_annotation_position_var = tk.StringVar(value="upper left")

        ttk.Label(props, text="PDF page:").grid(row=0, column=0, sticky="w", padx=(0, 4))
        ttk.Spinbox(props, from_=1, to=999, textvariable=self.panel_page_var, width=7).grid(
            row=0, column=1, sticky="w"
        )

        fields = [
            ("Panel label", self.panel_label_var),
            ("Title override", self.panel_title_var),
            ("X label override", self.panel_xlabel_var),
            ("Y label override", self.panel_ylabel_var),
        ]
        for r, (label, var) in enumerate(fields, start=1):
            ttk.Label(props, text=label + ":").grid(row=r, column=0, sticky="w", padx=(0, 4))
            ttk.Entry(props, textvariable=var).grid(row=r, column=1, sticky="ew")

        # Axis limits in a compact 2x2 block.
        limits = ttk.Frame(props)
        limits.grid(row=5, column=0, columnspan=2, sticky="ew", pady=(4, 2))
        for c in range(4):
            limits.columnconfigure(c, weight=1)

        ttk.Label(limits, text="X min").grid(row=0, column=0)
        ttk.Entry(limits, textvariable=self.panel_xmin_var, width=8).grid(row=0, column=1, sticky="ew")
        ttk.Label(limits, text="X max").grid(row=0, column=2)
        ttk.Entry(limits, textvariable=self.panel_xmax_var, width=8).grid(row=0, column=3, sticky="ew")
        ttk.Label(limits, text="Y min").grid(row=1, column=0)
        ttk.Entry(limits, textvariable=self.panel_ymin_var, width=8).grid(row=1, column=1, sticky="ew")
        ttk.Label(limits, text="Y max").grid(row=1, column=2)
        ttk.Entry(limits, textvariable=self.panel_ymax_var, width=8).grid(row=1, column=3, sticky="ew")

        modes = ttk.Frame(props)
        modes.grid(row=6, column=0, columnspan=2, sticky="ew")
        ttk.Label(modes, text="Log Y:").pack(side="left")
        ttk.Combobox(
            modes, textvariable=self.panel_logy_var,
            values=("root", "yes", "no"), state="readonly", width=6
        ).pack(side="left", padx=(2, 8))
        ttk.Label(modes, text="Log X:").pack(side="left")
        ttk.Combobox(
            modes, textvariable=self.panel_logx_var,
            values=("root", "yes", "no"), state="readonly", width=6
        ).pack(side="left", padx=2)
        ttk.Checkbutton(modes, text="Legend", variable=self.panel_legend_var).pack(side="right")
        ttk.Checkbutton(
            modes, text="Connect points", variable=self.panel_connect_var
        ).pack(side="right", padx=(0, 8))

        annotation = ttk.Frame(props)
        annotation.grid(row=7, column=0, columnspan=2, sticky="ew", pady=(4, 1))
        annotation.columnconfigure(1, weight=1)
        ttk.Label(annotation, text="Text:").grid(row=0, column=0, sticky="w")
        ttk.Entry(
            annotation, textvariable=self.panel_annotation_text_var
        ).grid(row=0, column=1, sticky="ew", padx=(4, 4))
        ttk.Combobox(
            annotation,
            textvariable=self.panel_annotation_position_var,
            values=("upper left", "upper right", "lower left", "lower right"),
            state="readonly",
            width=12,
        ).grid(row=0, column=2, sticky="e")

        # Editable legend/series labels.
        series_box = ttk.LabelFrame(props, text="Series / legend labels", padding=3)
        series_box.grid(row=8, column=0, columnspan=2, sticky="ew", pady=(5, 2))
        series_box.columnconfigure(0, weight=1)

        self.series_tree = ttk.Treeview(
            series_box, columns=("label",), show="tree headings", height=5, selectmode="browse"
        )
        self.series_tree.heading("#0", text="ROOT object")
        self.series_tree.heading("label", text="Legend label")
        self.series_tree.column("#0", width=185, stretch=True)
        self.series_tree.column("label", width=150, stretch=True)
        self.series_tree.grid(row=0, column=0, sticky="ew")
        self.series_tree.bind("<<TreeviewSelect>>", self.on_series_select)

        series_edit = ttk.Frame(series_box)
        series_edit.grid(row=1, column=0, sticky="ew", pady=(3, 0))
        series_edit.columnconfigure(0, weight=1)
        self.series_label_var = tk.StringVar()
        ttk.Entry(series_edit, textvariable=self.series_label_var).grid(row=0, column=0, sticky="ew")
        ttk.Button(
            series_edit,
            text="Apply label → selected panels",
            command=self.apply_series_label,
        ).grid(row=0, column=1, padx=(4, 0))
        ttk.Button(
            series_edit,
            text="Copy all labels → selected panels",
            command=self.copy_all_series_labels_to_selected,
        ).grid(row=1, column=0, columnspan=2, sticky="ew", pady=(3, 0))

        ttk.Button(
            props,
            text="Apply settings → selected panels",
            command=lambda: self.apply_panel_settings(batch=True),
        ).grid(row=9, column=0, columnspan=2, sticky="ew", pady=(5, 0))

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

        self.preview_output_page_var = tk.IntVar(value=1)
        ttk.Label(outbuttons, text="Preview page").pack(side="left")
        ttk.Spinbox(
            outbuttons, from_=1, to=999,
            textvariable=self.preview_output_page_var, width=5
        ).pack(side="left", padx=(3, 6))

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

            # Cached PNGs are page-name based.  Remove them when a different
            # ROOT file is opened so page_0001 from two files cannot mix.
            for old_png in self.preview_temp_dir.glob("*.png"):
                try:
                    old_png.unlink()
                except Exception:
                    pass

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
            self.render_root_preview(meta, pad, label)
        except Exception as exc:
            self._show_exception("Preview failed", exc)

    def populate_object_list(self, pad: Any):
        self.object_list.delete(0, tk.END)
        for obj in iter_primitives(pad):
            self.object_list.insert(
                tk.END,
                f"{obj.ClassName():<22}  {root_name(obj)}"
            )

    def render_root_preview(self, meta: Dict[str, Any], pad: Any, label: str):
        """
        Safely preview a ROOT page/pad.

        Important: the GUI process never calls TCanvas::Print/TPad::Print.
        The complete page is rendered in a subprocess.  A selected sub-pad is
        then obtained by cropping that page PNG using the saved ROOT pad NDC
        geometry.  Thus even a ROOT painter segmentation fault cannot kill the
        GUI.
        """
        page = meta["page"]
        page_png = self.preview_temp_dir / f"{page}.png"
        renderer_note = ""

        if not page_png.exists() or page_png.stat().st_size == 0:
            cmd = [
                sys.executable,
                str(Path(__file__).resolve()),
                str(self.model.path),
                "--_render-preview",
                "--preview-page", page,
                "--preview-output", str(page_png),
            ]

            try:
                proc = subprocess.run(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=60,
                )
            except subprocess.TimeoutExpired:
                proc = None
                renderer_note = "ROOT preview renderer timed out"

            if proc is not None and proc.returncode != 0:
                renderer_note = f"ROOT preview renderer exited with code {proc.returncode}"

            # ROOT sometimes writes a complete PNG and then crashes while
            # cleaning up the deserialized canvas.  If the PNG exists and
            # Pillow can read it, it is still perfectly useful as a preview.
            if not page_png.exists() or page_png.stat().st_size == 0:
                details = ""
                if proc is not None:
                    details = (proc.stderr or proc.stdout or "").strip()
                    if len(details) > 1200:
                        details = details[-1200:]
                raise RuntimeError(
                    f"Could not render {page} in the isolated ROOT helper."
                    + (f"\n\n{details}" if details else "")
                )

        self.update_idletasks()
        max_w = max(300, self.preview_frame.winfo_width() - 10)
        max_h = max(250, self.preview_frame.winfo_height() - 10)

        if HAVE_PIL:
            try:
                img = Image.open(page_png)
                img.load()
            except Exception:
                # Do not keep a half-written PNG in the cache.
                try:
                    page_png.unlink()
                except Exception:
                    pass
                raise

            # For a sub-pad, crop from the already-rendered FULL canvas.
            # No ROOT repaint is needed.
            if meta.get("pad_path"):
                bbox = pad_abs_ndc_bbox(pad)
                if bbox is not None:
                    x, y, w, h = bbox
                    iw, ih = img.size

                    left   = max(0, min(iw, int(round(x * iw))))
                    right  = max(0, min(iw, int(round((x + w) * iw))))
                    top    = max(0, min(ih, int(round((1.0 - (y + h)) * ih))))
                    bottom = max(0, min(ih, int(round((1.0 - y) * ih))))

                    if right > left and bottom > top:
                        img = img.crop((left, top, right, bottom))

            img.thumbnail((max_w, max_h))
            self.preview_photo = ImageTk.PhotoImage(img)
        else:
            # Without Pillow we cannot crop a sub-pad safely, so show the
            # complete page.  This still avoids calling ROOT painter in GUI.
            self.preview_photo = tk.PhotoImage(file=str(page_png))

        title = label
        if renderer_note:
            title += f"   [{renderer_note}; PNG recovered]"
        self.preview_title_var.set(title)
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

        # Try to recover legend labels from anywhere on the full ROOT page.
        page_labels = legend_label_map(self.model.get_page(meta["page"]), recursive=True)
        series_labels = {}
        for obj in drawable_series_objects(pad):
            name = root_name(obj)
            if name in page_labels:
                series_labels[name] = page_labels[name]

        panel = {
            "page": meta["page"],
            "pad_path": list(meta["pad_path"]),
            "figure_page": max([int(p.get("figure_page", 1)) for p in self.panels], default=1),
            "panel_label": self._next_panel_label(),
            "title": "",
            "xlabel": "",
            "ylabel": "",
            "xmin": "",
            "xmax": "",
            "ymin": "",
            "ymax": "",
            "logy": "root",
            "logx": "root",
            "legend": True,
            "connect_points": True,
            "annotation_text": "",
            "annotation_position": "upper left",
            "series_labels": series_labels,
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
        where = p["page"] + (f" / {path}" if path else "")
        return f"[P{int(p.get('figure_page', 1))}] {where}"

    def refresh_panel_list(self, select: Optional[int] = None):
        self.panel_list.delete(0, tk.END)
        for i, p in enumerate(self.panels, start=1):
            self.panel_list.insert(tk.END, f"{i:02d}. {self.panel_display_name(p)}")

        if select is not None and 0 <= select < len(self.panels):
            self.panel_list.selection_set(select)
            self.panel_list.activate(select)
            self.on_panel_select()

    def selected_panel_indices(self) -> List[int]:
        return [int(i) for i in self.panel_list.curselection()]

    def selected_panel_index(self) -> Optional[int]:
        sel = self.selected_panel_indices()
        if not sel:
            return None

        # Prefer the active row for editing if it belongs to the selection.
        try:
            active = int(self.panel_list.index(tk.ACTIVE))
            if active in sel:
                return active
        except Exception:
            pass

        return sel[0]

    def select_all_panels(self):
        if not self.panels:
            return
        self.panel_list.selection_set(0, tk.END)
        self.panel_list.activate(0)
        self.on_panel_select()

    def on_panel_select(self, *_args):
        idx = self.selected_panel_index()
        if idx is None:
            return
        p = self.panels[idx]
        self.panel_page_var.set(int(p.get("figure_page", 1)))
        self.panel_label_var.set(p.get("panel_label", ""))
        self.panel_title_var.set(p.get("title", ""))
        self.panel_xlabel_var.set(p.get("xlabel", ""))
        self.panel_ylabel_var.set(p.get("ylabel", ""))
        self.panel_xmin_var.set(str(p.get("xmin", "")))
        self.panel_xmax_var.set(str(p.get("xmax", "")))
        self.panel_ymin_var.set(str(p.get("ymin", "")))
        self.panel_ymax_var.set(str(p.get("ymax", "")))
        self.panel_logy_var.set(p.get("logy", "root"))
        self.panel_logx_var.set(p.get("logx", "root"))
        self.panel_legend_var.set(bool(p.get("legend", True)))
        self.panel_connect_var.set(bool(p.get("connect_points", True)))
        self.panel_annotation_text_var.set(p.get("annotation_text", ""))
        self.panel_annotation_position_var.set(
            p.get("annotation_position", "upper left")
        )
        self.refresh_series_tree()

    def apply_panel_settings(self, batch: bool = False):
        active_idx = self.selected_panel_index()
        if active_idx is None:
            return

        # Validate numeric limits here, so typos are caught in the GUI instead
        # of later during PDF export.
        for label, var in (
            ("X min", self.panel_xmin_var),
            ("X max", self.panel_xmax_var),
            ("Y min", self.panel_ymin_var),
            ("Y max", self.panel_ymax_var),
        ):
            raw = var.get().strip()
            if raw:
                try:
                    optional_float(raw)
                except Exception:
                    messagebox.showerror(
                        "Axis range",
                        f"{label} is not a valid finite number: {raw}"
                    )
                    return

        if batch:
            indices = self.selected_panel_indices()
            if not indices:
                indices = [active_idx]
        else:
            indices = [active_idx]

        page_no = max(1, int(self.panel_page_var.get()))

        for idx in indices:
            p = self.panels[idx]

            # Panel labels are intentionally kept unique during a multi-panel
            # batch operation.  With a single target they remain editable.
            if len(indices) == 1:
                p["panel_label"] = self.panel_label_var.get()

            p["figure_page"] = page_no
            p["title"] = self.panel_title_var.get()
            p["xlabel"] = self.panel_xlabel_var.get()
            p["ylabel"] = self.panel_ylabel_var.get()
            p["xmin"] = self.panel_xmin_var.get().strip()
            p["xmax"] = self.panel_xmax_var.get().strip()
            p["ymin"] = self.panel_ymin_var.get().strip()
            p["ymax"] = self.panel_ymax_var.get().strip()
            p["logy"] = self.panel_logy_var.get()
            p["logx"] = self.panel_logx_var.get()
            p["legend"] = bool(self.panel_legend_var.get())
            p["connect_points"] = bool(self.panel_connect_var.get())
            p["annotation_text"] = self.panel_annotation_text_var.get()
            p["annotation_position"] = self.panel_annotation_position_var.get()

        # Preserve the multi-selection after refreshing the displayed names.
        selected = list(indices)
        self.refresh_panel_list()
        for idx in selected:
            if 0 <= idx < len(self.panels):
                self.panel_list.selection_set(idx)
        if active_idx < len(self.panels):
            self.panel_list.activate(active_idx)
            self.panel_list.see(active_idx)

    def refresh_series_tree(self):
        self.series_tree.delete(*self.series_tree.get_children())
        idx = self.selected_panel_index()
        if idx is None:
            return

        p = self.panels[idx]
        try:
            pad = self.model.get_pad(p["page"], p.get("pad_path", []))
            page_labels = legend_label_map(self.model.get_page(p["page"]), recursive=True)
        except Exception:
            return

        custom = p.setdefault("series_labels", {})

        for i, obj in enumerate(drawable_series_objects(pad)):
            name = root_name(obj)
            label = custom.get(name, page_labels.get(name, safe_label(obj)))
            iid = f"series::{i}"
            self.series_tree.insert("", "end", iid=iid, text=name, values=(label,))

    def on_series_select(self, *_args):
        sel = self.series_tree.selection()
        if not sel:
            return
        vals = self.series_tree.item(sel[0], "values")
        self.series_label_var.set(vals[0] if vals else "")

    def apply_series_label(self):
        active_idx = self.selected_panel_index()
        sel = self.series_tree.selection()
        if active_idx is None or not sel:
            return

        iid = sel[0]
        try:
            series_index = int(iid.split("::", 1)[1])
        except Exception:
            return

        label = self.series_label_var.get().strip()
        target_indices = self.selected_panel_indices() or [active_idx]

        for idx in target_indices:
            p = self.panels[idx]
            try:
                pad = self.model.get_pad(p["page"], p.get("pad_path", []))
                objs = drawable_series_objects(pad)
                if series_index >= len(objs):
                    continue
                obj_name = root_name(objs[series_index])
                p.setdefault("series_labels", {})[obj_name] = label
            except Exception:
                continue

        self.refresh_series_tree()

        # Re-select the same series row.
        target_iid = f"series::{series_index}"
        if self.series_tree.exists(target_iid):
            self.series_tree.selection_set(target_iid)
            self.series_tree.focus(target_iid)

    def copy_all_series_labels_to_selected(self):
        active_idx = self.selected_panel_index()
        if active_idx is None:
            return

        source_panel = self.panels[active_idx]
        try:
            source_pad = self.model.get_pad(
                source_panel["page"], source_panel.get("pad_path", [])
            )
            source_objs = drawable_series_objects(source_pad)
            page_labels = legend_label_map(
                self.model.get_page(source_panel["page"]),
                recursive=True,
            )
        except Exception as exc:
            self._show_exception("Cannot read source legend labels", exc)
            return

        source_custom = source_panel.setdefault("series_labels", {})
        effective_labels: List[str] = []
        for obj in source_objs:
            name = root_name(obj)
            effective_labels.append(
                source_custom.get(name, page_labels.get(name, safe_label(obj)))
            )

        target_indices = self.selected_panel_indices() or [active_idx]

        for idx in target_indices:
            p = self.panels[idx]
            try:
                pad = self.model.get_pad(p["page"], p.get("pad_path", []))
                objs = drawable_series_objects(pad)
            except Exception:
                continue

            custom = p.setdefault("series_labels", {})
            for i, label in enumerate(effective_labels):
                if i >= len(objs):
                    break
                custom[root_name(objs[i])] = label

        self.refresh_series_tree()

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
        indices = self.selected_panel_indices()
        if not indices:
            return
        first = min(indices)
        for idx in sorted(indices, reverse=True):
            self.panels.pop(idx)
        self.refresh_panel_list(
            select=min(first, len(self.panels) - 1) if self.panels else None
        )

    def clear_panels(self):
        self.panels.clear()
        self.refresh_panel_list()

    # ---------------- Matplotlib figure ----------------

    def figure_page_numbers(self) -> List[int]:
        return sorted({max(1, int(p.get("figure_page", 1))) for p in self.panels})

    def panels_for_figure_page(self, page_number: int) -> List[Dict[str, Any]]:
        return [
            p for p in self.panels
            if max(1, int(p.get("figure_page", 1))) == int(page_number)
        ]

    def build_figure(self, page_number: Optional[int] = None):
        if not self.panels:
            raise RuntimeError("No panels selected.")

        self.apply_panel_settings()

        page_numbers = self.figure_page_numbers()
        if page_number is None:
            page_number = page_numbers[0]

        page_panels = self.panels_for_figure_page(page_number)
        if not page_panels:
            raise RuntimeError(
                f"No panels are assigned to output PDF page {page_number}. "
                f"Available pages: {page_numbers}"
            )

        n = len(page_panels)
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

        # Cache one legend-name map per stored ROOT page.
        page_legend_cache: Dict[str, Dict[str, str]] = {}

        for i, p in enumerate(page_panels):
            stored_page = p["page"]
            pad = self.model.get_pad(stored_page, p.get("pad_path", []))

            if stored_page not in page_legend_cache:
                page_legend_cache[stored_page] = legend_label_map(
                    self.model.get_page(stored_page),
                    recursive=True,
                )

            render_pad_to_axis(
                fig,
                flat_axes[i],
                pad,
                p,
                page_legend_labels=page_legend_cache[stored_page],
            )

        for ax in flat_axes[n:]:
            ax.set_visible(False)

        fig.tight_layout()
        return fig

    def preview_figure(self):
        try:
            page_number = int(self.preview_output_page_var.get())
            fig = self.build_figure(page_number=page_number)
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
            pages = self.figure_page_numbers()

            if ext.lower() == "pdf":
                # A true multi-page PDF. One Matplotlib figure is generated for
                # every distinct panel "PDF page" number.
                with PdfPages(path) as pdf:
                    for page_number in pages:
                        fig = self.build_figure(page_number=page_number)
                        pdf.savefig(fig, bbox_inches="tight")
                        plt.close(fig)

            elif ext.lower() == "png":
                # PNG itself cannot be multi-page. For >1 page, save a numbered
                # family: name_page01.png, name_page02.png, ...
                base = Path(path)
                if len(pages) == 1:
                    fig = self.build_figure(page_number=pages[0])
                    fig.savefig(str(base), bbox_inches="tight")
                    plt.close(fig)
                else:
                    for page_number in pages:
                        out = base.with_name(
                            f"{base.stem}_page{page_number:02d}{base.suffix}"
                        )
                        fig = self.build_figure(page_number=page_number)
                        fig.savefig(str(out), bbox_inches="tight")
                        plt.close(fig)
            else:
                raise RuntimeError(f"Unsupported output format: {ext}")

            messagebox.showinfo(
                "Saved",
                f"Saved {len(pages)} figure page(s):\n{path}"
            )

        except Exception as exc:
            self._show_exception("Cannot save figure", exc)

    # ---------------- JSON config ----------------

    def current_config(self) -> Dict[str, Any]:
        self.apply_panel_settings()
        return {
            "version": 4,
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
                    "figure_page": max(1, int(p.get("figure_page", 1))),
                    "panel_label": p.get("panel_label", ""),
                    "title": p.get("title", ""),
                    "xlabel": p.get("xlabel", ""),
                    "ylabel": p.get("ylabel", ""),
                    "xmin": str(p.get("xmin", "")),
                    "xmax": str(p.get("xmax", "")),
                    "ymin": str(p.get("ymin", "")),
                    "ymax": str(p.get("ymax", "")),
                    "logy": p.get("logy", "root"),
                    "logx": p.get("logx", "root"),
                    "legend": bool(p.get("legend", True)),
                    "connect_points": bool(p.get("connect_points", True)),
                    "annotation_text": p.get("annotation_text", ""),
                    "annotation_position": p.get(
                        "annotation_position", "upper left"
                    ),
                    "series_labels": dict(p.get("series_labels", {}) or {}),
                }
                validated.append(q)

            self.panels = validated
            self.refresh_panel_list(select=0 if self.panels else None)
            if self.panels:
                self.preview_output_page_var.set(self.figure_page_numbers()[0])

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

    # Internal crash-isolated ROOT preview mode.  These options are hidden from
    # ordinary --help output because users should never need to call them.
    ap.add_argument("--_render-preview", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--preview-page", help=argparse.SUPPRESS)
    ap.add_argument("--preview-output", help=argparse.SUPPRESS)

    return ap.parse_args()


def main():
    args = parse_args()

    if args._render_preview:
        if not args.root_file or not args.preview_page or not args.preview_output:
            print("[preview-helper] Missing helper arguments.", file=sys.stderr)
            raise SystemExit(2)

        code = render_page_preview_in_child(
            args.root_file,
            args.preview_page,
            args.preview_output,
        )
        raise SystemExit(code)

    app = CanvasBrowserApp(
        initial_root=args.root_file,
        initial_config=args.config,
    )
    app.mainloop()


if __name__ == "__main__":
    main()
