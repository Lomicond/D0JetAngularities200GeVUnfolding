#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_root_hist_bins_gui.py

Small GUI/CLI helper to inspect selected TH1 bins in a ROOT file.

Default checks:
  Lambda3_2_it3_ICS                 bin 7
  RCP_5_20_Lambda_4_0_it3_ICS       bin 1
  RCP_5_20_Lambda_5_0_it3_ICS       bin 7

The bin number is ROOT-style, i.e. 1 = first visible bin.
The script prints bin content, absolute statistical error, and relative
statistical error in percent: 100 * bin_error / abs(bin_content).

Dependencies:
  Preferred: PyROOT
  Fallback: uproot
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


DEFAULT_CHECKS = [
    ("Lambda3_2_it3_ICS", 7),
    ("RCP_5_20_Lambda_4_0_it3_ICS", 1),
    ("RCP_5_20_Lambda_5_2_it3_ICS", 7),
]


@dataclass
class BinResult:
    requested_name: str
    resolved_name: str
    bin_index: int
    value: float
    stat_abs: float
    stat_rel_pct: float
    note: str = ""


def finite(x: float) -> bool:
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def rel_stat_pct(value: float, err: float) -> float:
    if not finite(value) or not finite(err) or value == 0.0:
        return math.nan
    return 100.0 * abs(err) / abs(value)


def candidate_names(name: str) -> List[str]:
    """Try exact name first, then common Lambda4 vs Lambda_4 variants."""
    out: List[str] = []

    def add(x: str) -> None:
        if x and x not in out:
            out.append(x)

    s = str(name).strip()
    add(s)

    # ROOT sometimes stores objects as name;1.  TFile.Get usually does not need
    # this, but uproot can use it.
    if ";" not in s:
        add(s + ";1")

    # User-friendly fallback:
    #   RCP_5_20_Lambda_4_0_it3_ICS -> RCP_5_20_Lambda4_0_it3_ICS
    compact = re.sub(r"Lambda_(\d+)", r"Lambda\1", s)
    add(compact)
    if ";" not in compact:
        add(compact + ";1")

    # Reverse fallback:
    #   Lambda4_0_it3_ICS -> Lambda_4_0_it3_ICS
    expanded = re.sub(r"Lambda(\d+)", r"Lambda_\1", s)
    add(expanded)
    if ";" not in expanded:
        add(expanded + ";1")

    return out


class RootBackend:
    def __init__(self, path: str):
        self.path = str(path)
        self.backend = ""
        self.root_file = None
        self.uproot_file = None

        try:
            import ROOT  # type: ignore

            ROOT.gROOT.SetBatch(True)
            f = ROOT.TFile.Open(self.path, "READ")
            if not f or f.IsZombie():
                raise OSError(f"Cannot open ROOT file with PyROOT: {self.path}")
            self.backend = "PyROOT"
            self.root_file = f
            self.ROOT = ROOT
            return
        except Exception as e_root:
            self.root_error = e_root

        try:
            import uproot  # type: ignore

            self.backend = "uproot"
            self.uproot_file = uproot.open(self.path)
            return
        except Exception as e_uproot:
            raise RuntimeError(
                "Could not open ROOT file. Neither PyROOT nor uproot worked.\n"
                f"PyROOT error: {self.root_error}\n"
                f"uproot error: {e_uproot}"
            )

    def close(self) -> None:
        try:
            if self.root_file is not None:
                self.root_file.Close()
        except Exception:
            pass
        try:
            if self.uproot_file is not None:
                self.uproot_file.close()
        except Exception:
            pass

    def _get_pyroot_hist(self, name: str):
        assert self.root_file is not None
        for cand in candidate_names(name):
            obj = self.root_file.Get(cand)
            if obj:
                return cand.split(";", 1)[0], obj
        return "", None

    def _get_uproot_hist(self, name: str):
        assert self.uproot_file is not None
        for cand in candidate_names(name):
            try:
                obj = self.uproot_file[cand]
                return cand, obj
            except Exception:
                continue
        return "", None

    def read_bin(self, hist_name: str, bin_index: int) -> BinResult:
        if int(bin_index) < 1:
            raise ValueError(f"Bin index must be ROOT-style 1-based, got {bin_index}")
        ibin = int(bin_index)

        if self.backend == "PyROOT":
            resolved, h = self._get_pyroot_hist(hist_name)
            if h is None:
                raise KeyError(f"Histogram not found: {hist_name}")
            nb = int(h.GetNbinsX())
            if ibin > nb:
                raise IndexError(f"Requested bin {ibin}, but histogram {resolved} has only {nb} bins")
            value = float(h.GetBinContent(ibin))
            err = float(h.GetBinError(ibin))
            return BinResult(hist_name, resolved, ibin, value, err, rel_stat_pct(value, err))

        resolved, h = self._get_uproot_hist(hist_name)
        if h is None:
            raise KeyError(f"Histogram not found: {hist_name}")
        try:
            values = h.values(flow=False)
            variances = h.variances(flow=False)
        except Exception as e:
            raise TypeError(f"Object {resolved} does not look like a TH1 histogram: {e}")
        nb = len(values)
        if ibin > nb:
            raise IndexError(f"Requested bin {ibin}, but histogram {resolved} has only {nb} bins")
        value = float(values[ibin - 1])
        if variances is None:
            err = math.nan
            note = "missing variances/errors"
        else:
            err = math.sqrt(max(0.0, float(variances[ibin - 1])))
            note = ""
        return BinResult(hist_name, resolved.split(";", 1)[0], ibin, value, err, rel_stat_pct(value, err), note)


def format_results(path: str, results: List[BinResult], backend: str) -> str:
    lines = []
    lines.append(f"ROOT file: {path}")
    lines.append(f"Backend  : {backend}")
    lines.append("")
    header = f"{'requested hist':42s} {'resolved hist':42s} {'bin':>4s} {'value':>16s} {'stat abs':>16s} {'stat rel [%]':>14s}"
    lines.append(header)
    lines.append("-" * len(header))
    for r in results:
        rel = "nan" if not finite(r.stat_rel_pct) else f"{r.stat_rel_pct:.6g}"
        val = "nan" if not finite(r.value) else f"{r.value:.10g}"
        err = "nan" if not finite(r.stat_abs) else f"{r.stat_abs:.10g}"
        lines.append(
            f"{r.requested_name:42s} {r.resolved_name:42s} {r.bin_index:4d} {val:>16s} {err:>16s} {rel:>14s}"
        )
        if r.note:
            lines.append(f"  note for {r.resolved_name}: {r.note}")
    return "\n".join(lines)


def run_checks(path: str, checks: List[Tuple[str, int]]) -> str:
    rb = RootBackend(path)
    try:
        results = [rb.read_bin(h, b) for h, b in checks]
        return format_results(path, results, rb.backend)
    finally:
        rb.close()


def parse_check(text: str) -> Tuple[str, int]:
    """Parse 'hist:bin' or 'hist,bin'."""
    if ":" in text:
        hist, b = text.rsplit(":", 1)
    elif "," in text:
        hist, b = text.rsplit(",", 1)
    else:
        raise ValueError(f"Check must be hist:bin, got {text!r}")
    hist = hist.strip()
    if not hist:
        raise ValueError("Empty histogram name")
    return hist, int(str(b).strip())


def launch_gui(initial_file: str = "") -> None:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    root = tk.Tk()
    root.title("ROOT histogram bin checker")
    root.geometry("1180x610")

    file_var = tk.StringVar(value=initial_file)
    hist_vars: List[tk.StringVar] = []
    bin_vars: List[tk.StringVar] = []

    top = ttk.Frame(root, padding=8)
    top.pack(fill="x")
    top.columnconfigure(1, weight=1)

    ttk.Label(top, text="ROOT file:").grid(row=0, column=0, sticky="w", padx=(0, 6))
    ttk.Entry(top, textvariable=file_var).grid(row=0, column=1, sticky="ew", padx=(0, 6))

    def browse() -> None:
        chosen = filedialog.askopenfilename(
            title="Select ROOT file",
            filetypes=[("ROOT files", "*.root"), ("All files", "*")],
        )
        if chosen:
            file_var.set(chosen)

    ttk.Button(top, text="Browse", command=browse).grid(row=0, column=2, sticky="ew")

    table = ttk.LabelFrame(root, text="Checks: histogram name + ROOT-style bin index", padding=8)
    table.pack(fill="x", padx=8, pady=4)
    table.columnconfigure(0, weight=1)

    ttk.Label(table, text="Histogram name").grid(row=0, column=0, sticky="w")
    ttk.Label(table, text="Bin").grid(row=0, column=1, sticky="w", padx=(8, 0))

    rows_frame = ttk.Frame(table)
    rows_frame.grid(row=1, column=0, columnspan=3, sticky="ew")
    rows_frame.columnconfigure(0, weight=1)

    def redraw_rows() -> None:
        for child in rows_frame.winfo_children():
            child.destroy()
        for i, (hv, bv) in enumerate(zip(hist_vars, bin_vars)):
            ttk.Entry(rows_frame, textvariable=hv).grid(row=i, column=0, sticky="ew", pady=2)
            ttk.Entry(rows_frame, textvariable=bv, width=8).grid(row=i, column=1, sticky="w", padx=(8, 0), pady=2)
            ttk.Button(rows_frame, text="×", width=3, command=lambda idx=i: remove_row(idx)).grid(row=i, column=2, sticky="w", padx=(6, 0), pady=2)
        rows_frame.columnconfigure(0, weight=1)

    def add_row(hist: str = "", bin_index: int = 1) -> None:
        hist_vars.append(tk.StringVar(value=hist))
        bin_vars.append(tk.StringVar(value=str(bin_index)))
        redraw_rows()

    def remove_row(idx: int) -> None:
        if 0 <= idx < len(hist_vars):
            hist_vars.pop(idx)
            bin_vars.pop(idx)
            redraw_rows()

    for h, b in DEFAULT_CHECKS:
        add_row(h, b)

    btns = ttk.Frame(root, padding=(8, 2))
    btns.pack(fill="x")

    output = tk.Text(root, wrap="none", height=23)
    output.pack(fill="both", expand=True, padx=8, pady=(4, 8))

    xscroll = ttk.Scrollbar(output, orient="horizontal", command=output.xview)
    output.configure(xscrollcommand=xscroll.set)

    def evaluate() -> None:
        path = file_var.get().strip()
        if not path:
            messagebox.showerror("Missing file", "Select a ROOT file first.")
            return
        if not Path(path).exists():
            messagebox.showerror("Missing file", f"File does not exist:\n{path}")
            return

        checks: List[Tuple[str, int]] = []
        try:
            for hv, bv in zip(hist_vars, bin_vars):
                h = hv.get().strip()
                if not h:
                    continue
                checks.append((h, int(bv.get().strip())))
            if not checks:
                raise ValueError("No checks defined.")
            text = run_checks(path, checks)
        except Exception as e:
            text = f"ERROR:\n{e}"

        output.delete("1.0", "end")
        output.insert("1.0", text)

    ttk.Button(btns, text="Add row", command=lambda: add_row("", 1)).pack(side="left")
    ttk.Button(btns, text="Evaluate", command=evaluate).pack(side="left", padx=8)
    ttk.Button(btns, text="Quit", command=root.destroy).pack(side="right")

    root.mainloop()


def main() -> int:
    p = argparse.ArgumentParser(description="Inspect selected TH1 bin values and relative statistical errors.")
    p.add_argument("root_file", nargs="?", default="", help="ROOT file. If omitted, GUI is opened.")
    p.add_argument(
        "--check",
        action="append",
        default=[],
        help="Check in format hist:bin. Can be given multiple times. Default uses the three built-in checks.",
    )
    p.add_argument("--gui", action="store_true", help="Force GUI mode even if root_file is given.")
    p.add_argument("--no-gui", action="store_true", help="Force command-line mode.")
    args = p.parse_args()

    if args.gui or (not args.no_gui and not args.root_file):
        launch_gui(args.root_file)
        return 0

    if not args.root_file:
        p.error("root_file is required in --no-gui mode")

    checks = [parse_check(x) for x in args.check] if args.check else list(DEFAULT_CHECKS)
    print(run_checks(args.root_file, checks))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
