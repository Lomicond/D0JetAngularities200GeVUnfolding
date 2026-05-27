#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
systematics_final_combiner_gui_v7.py

Small Tkinter GUI for editing JSON rules for the final systematic combiner.

Main workflow:
  1) Load Systematics/systematics_components_v15.tsv.
  2) Auto-discover systematic sources and create editable rules.
  3) Choose include/symmetrize/correlation/Barlow/combine mode per source.
  4) Save a JSON config.
  5) Optionally run systematics_final_combiner_v9.py.

This GUI writes a config compatible with systematics_final_combiner_v9.py.
It also remains mostly compatible with v2, but v2 ignores the new fields
"symmetrize", "total_mode" and "correlation".
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import threading
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Any, Dict, List, Tuple

import tkinter as tk
from tkinter import ttk, filedialog, messagebox


COMBINE_MODES = [
    "envelope",
    "max_abs",
    "take_rms",
    "rms",
    "quadrature",
    "linear",
    "pair_envelope_then_quadrature",
]

TOTAL_MODES = [
    "quadrature",   # uncorrelated
    "linear",       # correlated / conservative linear sum
    "envelope",     # take maximum bucket
    "exclude",
]

TOTAL_MODE_LABELS = {
    "quadrature": "uncorrelated/quadrature",
    "linear": "correlated/linear",
    "envelope": "envelope",
    "exclude": "exclude from total",
}

SYM_MODES = ["none", "max_abs", "average_abs"]
BARLOW_MODES = ["correlated", "independent"]
BARLOW_ACTIONS = ["flag_only", "zero_if_not_significant", "subtract_stat_diff"]

BASE_KEY_COLS = [
    "job_label", "observable", "observable_pretty", "hist_name",
    "centrality", "centrality_label", "method", "iteration_display",
    "iteration_root", "bin", "bin_low", "bin_high",
]


def safe_name(x: Any) -> str:
    s = str(x or "").strip().lower()
    s = s.replace("#", "")
    s = re.sub(r"[^a-zA-Z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "source"


def truthy(x: Any, default: bool = False) -> bool:
    if x is None:
        return default
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in ("1", "true", "yes", "y", "on")


def resolve_path(text: str, base: Path | None = None) -> Path:
    p = Path(os.path.expandvars(os.path.expanduser(str(text))))
    if not p.is_absolute() and base is not None:
        p = base / p
    return p


@dataclass
class Rule:
    include: bool = True
    name: str = ""
    source_group: str = ""
    source_name: str = ""
    variation_name: str = ""
    combine_variations: str = "envelope"
    total_mode: str = "quadrature"
    symmetrize: str = "none"
    apply_barlow: bool = True
    pair_by: str = ""
    notes: str = ""

    def to_group_config(self) -> Dict[str, Any]:
        select: Dict[str, Any] = {}
        if self.source_group:
            select["source_group"] = self.source_group
        if self.source_name:
            select["source_name"] = self.source_name
        if self.variation_name:
            select["variation_name"] = self.variation_name

        out: Dict[str, Any] = {
            "name": self.name,
            "select": select,
            "combine_variations": self.combine_variations,
            "total_mode": self.total_mode,
            "correlation": "correlated" if self.total_mode == "linear" else ("uncorrelated" if self.total_mode == "quadrature" else self.total_mode),
            "symmetrize": self.symmetrize,
            "apply_barlow": bool(self.apply_barlow),
        }
        if self.pair_by:
            out["pair_by"] = self.pair_by
        if self.notes:
            out["notes"] = self.notes
        return out

    @staticmethod
    def from_group_config(g: Dict[str, Any], include: bool = True) -> "Rule":
        sel = dict(g.get("select", {}))
        mode = str(g.get("total_mode", "") or "").strip()
        if not mode:
            corr = str(g.get("correlation", "") or "").strip().lower()
            if corr in ("correlated", "corr", "linear"):
                mode = "linear"
            elif corr in ("uncorrelated", "uncorr", "quadrature", "quad"):
                mode = "quadrature"
            elif corr in ("envelope", "max"):
                mode = "envelope"
            elif corr in ("exclude", "none", "off"):
                mode = "exclude"
            else:
                mode = "quadrature"
        return Rule(
            include=include,
            name=str(g.get("name", "")),
            source_group=str(sel.get("source_group", "")),
            source_name=str(sel.get("source_name", "")),
            variation_name=str(sel.get("variation_name", "")),
            combine_variations=str(g.get("combine_variations", "envelope")),
            total_mode=mode if mode in TOTAL_MODES else "quadrature",
            symmetrize=str(g.get("symmetrize", "none")),
            apply_barlow=truthy(g.get("apply_barlow", True), True),
            pair_by=str(g.get("pair_by", "")),
            notes=str(g.get("notes", "")),
        )


def read_unique_sources(path: Path) -> Tuple[List[Dict[str, str]], Dict[str, int]]:
    if not path.exists():
        raise FileNotFoundError(path)
    counts: Dict[Tuple[str, str, str], int] = {}
    rows_sample: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    total = 0
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            total += 1
            key = (row.get("source_group", ""), row.get("source_name", ""), row.get("variation_name", ""))
            counts[key] = counts.get(key, 0) + 1
            rows_sample.setdefault(key, row)
    samples = []
    for key, cnt in sorted(counts.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2])):
        row = dict(rows_sample[key])
        row["_count"] = str(cnt)
        samples.append(row)
    return samples, {"input_rows": total, "unique_components": len(samples)}


def discover_rules(path: Path, detailed: bool = True) -> List[Rule]:
    samples, _ = read_unique_sources(path)
    sg_to_names: Dict[str, set] = {}
    sg_sn_to_vars: Dict[Tuple[str, str], set] = {}
    for r in samples:
        sg = r.get("source_group", "")
        sn = r.get("source_name", "")
        vn = r.get("variation_name", "")
        if not sg:
            continue
        sg_to_names.setdefault(sg, set()).add(sn)
        sg_sn_to_vars.setdefault((sg, sn), set()).add(vn)

    rules: List[Rule] = []

    def add(rule: Rule):
        if rule.name and not any(r.name == rule.name for r in rules):
            rules.append(rule)

    # Nominal iteration variations.
    if ("Nominal", "iteration") in sg_sn_to_vars:
        add(Rule(
            name="iteration",
            source_group="Nominal",
            source_name="iteration",
            combine_variations="envelope",
            total_mode="quadrature",
            symmetrize="none",
            apply_barlow=False,
            notes="iteration +/-1 envelope",
        ))

    # Binning RMS sources.
    for sn in sorted(sg_to_names.get("Binning", set())):
        add(Rule(
            name=safe_name(sn),
            source_group="Binning",
            source_name=sn,
            combine_variations="take_rms",
            total_mode="quadrature",
            symmetrize="max_abs",
            apply_barlow=False,
            notes="RMS from binning scan",
        ))

    # PriorShape: usually first/second variable; keep as source-level rules.
    for sn in sorted(sg_to_names.get("PriorShape", set())):
        nice = safe_name(sn).replace("prior_shape_", "prior_shape_")
        add(Rule(
            name=nice,
            source_group="PriorShape",
            source_name=sn,
            combine_variations="envelope",
            total_mode="quadrature",
            symmetrize="none",
            apply_barlow=True,
            notes="prior-shape plus/minus envelope",
        ))

    # sWeight: detailed rules per source_name are useful for GUI control.
    for sn in sorted(sg_to_names.get("sWeight", set())):
        add(Rule(
            name="sweight_" + safe_name(sn),
            source_group="sWeight",
            source_name=sn,
            combine_variations="max_abs",
            total_mode="quadrature",
            symmetrize="max_abs",
            apply_barlow=True,
            notes="single sWeight variation source",
        ))

    # JetsReco: one source per reconstruction variation.
    for sn in sorted(sg_to_names.get("JetsReco", set())):
        add(Rule(
            name="jetsreco_" + safe_name(sn),
            source_group="JetsReco",
            source_name=sn,
            combine_variations="max_abs",
            total_mode="quadrature",
            symmetrize="max_abs",
            apply_barlow=True,
            notes="single jet-reconstruction variation source",
        ))

    # D0Meson: source_name groups the up/down pair, envelope inside each source.
    for sn in sorted(sg_to_names.get("D0Meson", set())):
        add(Rule(
            name="d0meson_" + safe_name(sn),
            source_group="D0Meson",
            source_name=sn,
            combine_variations="envelope",
            total_mode="quadrature",
            symmetrize="none",
            apply_barlow=True,
            notes="D0-meson up/down source envelope",
        ))

    if not detailed:
        # Compact preset compatible with the old v2 template.
        compact = [
            Rule(True, "iteration", "Nominal", "iteration", "", "envelope", "quadrature", "none", False),
            Rule(True, "binning_first_variable", "Binning", "binning_first_variable", "", "take_rms", "quadrature", "max_abs", False),
            Rule(True, "binning_second_variable", "Binning", "binning_second_variable", "", "take_rms", "quadrature", "max_abs", False),
            Rule(True, "prior_shape_first_variable", "PriorShape", "prior_shape_first_variable", "", "envelope", "quadrature", "none", True),
            Rule(True, "prior_shape_second_variable", "PriorShape", "prior_shape_second_variable", "", "envelope", "quadrature", "none", True),
            Rule(True, "sweight", "sWeight", "", "", "max_abs", "quadrature", "max_abs", True),
            Rule(True, "jetsreco", "JetsReco", "", "", "quadrature", "quadrature", "none", True),
            Rule(True, "d0meson", "D0Meson", "", "", "pair_envelope_then_quadrature", "quadrature", "none", True, "source_name"),
        ]
        return [r for r in compact if r.source_group in sg_to_names]

    return rules


class CombinerGui(tk.Tk):
    def __init__(self, initial_args: argparse.Namespace | None = None):
        super().__init__()
        self.title("Final systematics combiner GUI")
        self.geometry("1750x950")

        self.rules: List[Rule] = []
        self.loaded_config: Dict[str, Any] = {}
        self.stats: Dict[str, int] = {}

        self._build_vars()
        self._build_ui()
        self.apply_initial_args(initial_args)

    def _sv(self, value=""):
        return tk.StringVar(value=str(value))

    def _build_vars(self):
        self.input_var = self._sv("Systematics/systematics_components_v15.tsv")
        self.config_var = self._sv("Systematics/systematics_final_combiner_gui_config_v7.json")
        self.combiner_var = self._sv("Systematics/systematics_final_combiner_v9.py")
        self.output_tsv_var = self._sv("Systematics/final_systematics_gui_v7.tsv")
        self.output_csv_var = self._sv("Systematics/final_systematics_gui_v7.csv")
        self.details_tsv_var = self._sv("Systematics/final_systematics_gui_barlow_details_v7.tsv")
        self.output_root_var = self._sv("Systematics/final_systematics_results_v9.root")

        self.global_barlow_enabled = tk.BooleanVar(value=True)
        self.global_barlow_mode = self._sv("correlated")
        self.global_barlow_threshold = self._sv("1.0")
        self.global_barlow_action = self._sv("subtract_stat_diff")
        self.total_mode_var = self._sv("mixed")
        self.include_stat_var = tk.BooleanVar(value=True)
        self.default_group_mode_var = self._sv("quadrature")
        self.plots_enabled_var = tk.BooleanVar(value=True)

        # Editor variables.
        self.e_include = tk.BooleanVar(value=True)
        self.e_name = self._sv("")
        self.e_source_group = self._sv("")
        self.e_source_name = self._sv("")
        self.e_variation_name = self._sv("")
        self.e_combine = self._sv("envelope")
        self.e_total_mode = self._sv("quadrature")
        self.e_sym = self._sv("none")
        self.e_barlow = tk.BooleanVar(value=True)
        self.e_pair_by = self._sv("")
        self.e_notes = self._sv("")

        self.status_var = self._sv("Load components or config to start.")

    def _build_ui(self):
        self.columnconfigure(0, weight=1)
        self.rowconfigure(2, weight=1)
        self.rowconfigure(3, weight=0)

        top = ttk.LabelFrame(self, text="Files")
        top.grid(row=0, column=0, sticky="ew", padx=8, pady=6)
        for c in range(8):
            top.columnconfigure(c, weight=1)

        self._path_row(top, 0, "Input components TSV", self.input_var, self.browse_input)
        self._path_row(top, 1, "Config JSON", self.config_var, self.browse_config)
        self._path_row(top, 2, "Combiner script", self.combiner_var, self.browse_combiner)
        self._path_row(top, 3, "Output TSV", self.output_tsv_var, lambda: self.browse_save(self.output_tsv_var, ".tsv"))
        self._path_row(top, 4, "Output CSV", self.output_csv_var, lambda: self.browse_save(self.output_csv_var, ".csv"))
        self._path_row(top, 5, "Barlow details TSV", self.details_tsv_var, lambda: self.browse_save(self.details_tsv_var, ".tsv"))
        self._path_row(top, 6, "Output ROOT", self.output_root_var, lambda: self.browse_save(self.output_root_var, ".root"))

        btns = ttk.Frame(top)
        btns.grid(row=7, column=0, columnspan=8, sticky="ew", padx=4, pady=4)
        ttk.Button(btns, text="Discover detailed rules", command=lambda: self.discover_rules(True)).pack(side="left", padx=4)
        ttk.Button(btns, text="Discover compact preset", command=lambda: self.discover_rules(False)).pack(side="left", padx=4)
        ttk.Button(btns, text="Load config", command=self.load_config).pack(side="left", padx=12)
        ttk.Button(btns, text="Save config", command=self.save_config).pack(side="left", padx=4)
        ttk.Button(btns, text="Run combiner", command=self.run_combiner).pack(side="left", padx=12)
        ttk.Button(btns, text="Copy run command", command=self.copy_run_command).pack(side="left", padx=4)
        ttk.Label(btns, textvariable=self.status_var).pack(side="left", padx=16)

        opts = ttk.LabelFrame(self, text="Global settings")
        opts.grid(row=1, column=0, sticky="ew", padx=8, pady=4)
        for c in range(12):
            opts.columnconfigure(c, weight=1)
        ttk.Checkbutton(opts, text="Barlow enabled", variable=self.global_barlow_enabled).grid(row=0, column=0, sticky="w", padx=4, pady=2)
        ttk.Label(opts, text="Barlow mode").grid(row=0, column=1, sticky="e")
        ttk.Combobox(opts, textvariable=self.global_barlow_mode, values=BARLOW_MODES, state="readonly", width=14).grid(row=0, column=2, sticky="ew", padx=4)
        ttk.Label(opts, text="threshold").grid(row=0, column=3, sticky="e")
        ttk.Entry(opts, textvariable=self.global_barlow_threshold, width=8).grid(row=0, column=4, sticky="ew", padx=4)
        ttk.Label(opts, text="action").grid(row=0, column=5, sticky="e")
        ttk.Combobox(opts, textvariable=self.global_barlow_action, values=BARLOW_ACTIONS, state="readonly", width=22).grid(row=0, column=6, sticky="ew", padx=4)
        ttk.Label(opts, text="Total mode").grid(row=0, column=7, sticky="e")
        ttk.Combobox(opts, textvariable=self.total_mode_var, values=["mixed", "quadrature", "linear", "envelope"], state="readonly", width=12).grid(row=0, column=8, sticky="ew", padx=4)
        ttk.Label(opts, text="default group mode").grid(row=0, column=9, sticky="e")
        ttk.Combobox(opts, textvariable=self.default_group_mode_var, values=TOTAL_MODES[:-1], state="readonly", width=12).grid(row=0, column=10, sticky="ew", padx=4)
        ttk.Checkbutton(opts, text="include stat", variable=self.include_stat_var).grid(row=0, column=11, sticky="w", padx=4)
        ttk.Checkbutton(opts, text="make plots", variable=self.plots_enabled_var).grid(row=1, column=0, sticky="w", padx=4, pady=2)
        ttk.Label(opts, text="Mixed total convention: total = linear bucket + sqrt(quadrature bucket² + envelope bucket²). Use total_mode per row.").grid(row=1, column=1, columnspan=10, sticky="w", padx=4)

        main = ttk.PanedWindow(self, orient="horizontal")
        main.grid(row=2, column=0, sticky="nsew", padx=8, pady=4)

        left = ttk.LabelFrame(main, text="Rules / components")
        left.columnconfigure(0, weight=1)
        left.rowconfigure(0, weight=1)
        main.add(left, weight=4)

        cols = ["include", "name", "source_group", "source_name", "variation", "combine", "total", "sym", "barlow", "pair_by"]
        self.tree = ttk.Treeview(left, columns=cols, show="headings", selectmode="browse")
        self.tree.grid(row=0, column=0, sticky="nsew")
        widths = {"include": 70, "name": 220, "source_group": 100, "source_name": 170, "variation": 140, "combine": 170, "total": 100, "sym": 90, "barlow": 70, "pair_by": 90}
        for col in cols:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=widths.get(col, 120), anchor="center")
        self.tree.bind("<<TreeviewSelect>>", self.load_selected_to_editor)
        self.tree.bind("<Double-1>", self.double_click_toggle)
        y = ttk.Scrollbar(left, orient="vertical", command=self.tree.yview)
        x = ttk.Scrollbar(left, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=y.set, xscrollcommand=x.set)
        y.grid(row=0, column=1, sticky="ns")
        x.grid(row=1, column=0, sticky="ew")

        row_btns = ttk.Frame(left)
        row_btns.grid(row=2, column=0, columnspan=2, sticky="ew", pady=4)
        ttk.Button(row_btns, text="Toggle include", command=self.toggle_include).pack(side="left", padx=4)
        ttk.Button(row_btns, text="Toggle Barlow", command=self.toggle_barlow).pack(side="left", padx=4)
        ttk.Button(row_btns, text="Add", command=self.add_rule).pack(side="left", padx=12)
        ttk.Button(row_btns, text="Delete", command=self.delete_rule).pack(side="left", padx=4)
        ttk.Button(row_btns, text="Move up", command=lambda: self.move_rule(-1)).pack(side="left", padx=12)
        ttk.Button(row_btns, text="Move down", command=lambda: self.move_rule(+1)).pack(side="left", padx=4)

        right = ttk.LabelFrame(main, text="Selected rule editor")
        for c in range(4):
            right.columnconfigure(c, weight=1)
        main.add(right, weight=2)

        self._editor_row(right, 0, "Include", ttk.Checkbutton(right, variable=self.e_include))
        self._editor_entry(right, 1, "Rule name", self.e_name)
        self._editor_entry(right, 2, "source_group", self.e_source_group)
        self._editor_entry(right, 3, "source_name", self.e_source_name)
        self._editor_entry(right, 4, "variation_name", self.e_variation_name)
        self._editor_combo(right, 5, "combine variations", self.e_combine, COMBINE_MODES)
        self._editor_combo(right, 6, "total mode / correlation", self.e_total_mode, TOTAL_MODES)
        self._editor_combo(right, 7, "symmetrize", self.e_sym, SYM_MODES)
        self._editor_row(right, 8, "Apply Barlow", ttk.Checkbutton(right, variable=self.e_barlow))
        self._editor_entry(right, 9, "pair_by", self.e_pair_by)
        self._editor_entry(right, 10, "notes", self.e_notes)
        ttk.Button(right, text="Apply to selected", command=self.apply_editor_to_selected).grid(row=11, column=0, columnspan=2, sticky="ew", padx=4, pady=8)
        ttk.Button(right, text="Duplicate selected", command=self.duplicate_selected).grid(row=11, column=2, columnspan=2, sticky="ew", padx=4, pady=8)

        help_text = (
            "Hints:\n"
            "- include controls whether the rule enters total.include_groups.\n"
            "- total mode: quadrature = uncorrelated, linear = correlated, envelope = max bucket.\n"
            "- symmetrize=max_abs forces down/up to the larger side after the group is combined.\n"
            "- Barlow action subtract_stat_diff applies the presentation-style sqrt(Delta^2 - sigma_stat,diff^2) correction.\n"
            "- Barlow applies only where the combiner has variation and stat errors.\n"
            "- Double-click a row to toggle include."
        )
        ttk.Label(right, text=help_text, justify="left").grid(row=12, column=0, columnspan=4, sticky="ew", padx=4, pady=8)

        log_frame = ttk.LabelFrame(self, text="Run log")
        log_frame.grid(row=3, column=0, sticky="ew", padx=8, pady=4)
        self.log_text = tk.Text(log_frame, height=7, wrap="none")
        self.log_text.pack(fill="both", expand=True, padx=4, pady=4)

    def _path_row(self, parent, row, label, var, browse_cmd):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=4, pady=2)
        ttk.Entry(parent, textvariable=var).grid(row=row, column=1, columnspan=6, sticky="ew", padx=4, pady=2)
        ttk.Button(parent, text="...", command=browse_cmd).grid(row=row, column=7, sticky="ew", padx=4, pady=2)

    def _editor_row(self, parent, row, label, widget):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=4, pady=3)
        widget.grid(row=row, column=1, columnspan=3, sticky="w", padx=4, pady=3)

    def _editor_entry(self, parent, row, label, var):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=4, pady=3)
        ttk.Entry(parent, textvariable=var).grid(row=row, column=1, columnspan=3, sticky="ew", padx=4, pady=3)

    def _editor_combo(self, parent, row, label, var, values):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=4, pady=3)
        ttk.Combobox(parent, textvariable=var, values=values, state="readonly").grid(row=row, column=1, columnspan=3, sticky="ew", padx=4, pady=3)

    def browse_input(self):
        p = filedialog.askopenfilename(filetypes=[("TSV", "*.tsv"), ("All", "*.*")])
        if p:
            self.input_var.set(p)

    def browse_config(self):
        p = filedialog.askopenfilename(filetypes=[("JSON", "*.json"), ("All", "*.*")])
        if p:
            self.config_var.set(p)

    def browse_combiner(self):
        p = filedialog.askopenfilename(filetypes=[("Python", "*.py"), ("All", "*.*")])
        if p:
            self.combiner_var.set(p)

    def browse_save(self, var, ext):
        p = filedialog.asksaveasfilename(defaultextension=ext, filetypes=[("Output", f"*{ext}"), ("All", "*.*")])
        if p:
            var.set(p)

    def apply_initial_args(self, args: argparse.Namespace | None):
        if args is None:
            return
        if getattr(args, "config", ""):
            self.config_var.set(args.config)
        if getattr(args, "input", ""):
            self.input_var.set(args.input)
        if getattr(args, "combiner", ""):
            self.combiner_var.set(args.combiner)
        if getattr(args, "output_root", ""):
            self.output_root_var.set(args.output_root)

        cfg_path = Path(self.config_var.get().strip())
        if cfg_path.exists() and not getattr(args, "no_auto_load", False):
            self.load_config()
            # Command-line values should have final priority after loading.
            if getattr(args, "input", ""):
                self.input_var.set(args.input)
            if getattr(args, "combiner", ""):
                self.combiner_var.set(args.combiner)
            if getattr(args, "output_root", ""):
                self.output_root_var.set(args.output_root)

    def log(self, msg: str):
        self.log_text.insert("end", str(msg) + "\n")
        self.log_text.see("end")
        self.update_idletasks()

    def discover_rules(self, detailed: bool):
        try:
            path = Path(self.input_var.get().strip())
            self.rules = discover_rules(path, detailed=detailed)
            _, self.stats = read_unique_sources(path)
            self.refresh_tree()
            self.status_var.set(f"Discovered {len(self.rules)} rules from {self.stats.get('input_rows', 0)} rows.")
        except Exception as e:
            messagebox.showerror("Discover failed", str(e))

    def refresh_tree(self):
        self.tree.delete(*self.tree.get_children())
        for i, r in enumerate(self.rules):
            vals = [
                "yes" if r.include else "no",
                r.name,
                r.source_group,
                r.source_name,
                r.variation_name,
                r.combine_variations,
                r.total_mode,
                r.symmetrize,
                "yes" if r.apply_barlow else "no",
                r.pair_by,
            ]
            self.tree.insert("", "end", iid=str(i), values=vals)

    def selected_index(self):
        sel = self.tree.selection()
        if not sel:
            return None
        return int(sel[0])

    def load_selected_to_editor(self, event=None):
        idx = self.selected_index()
        if idx is None or idx >= len(self.rules):
            return
        r = self.rules[idx]
        self.e_include.set(r.include)
        self.e_name.set(r.name)
        self.e_source_group.set(r.source_group)
        self.e_source_name.set(r.source_name)
        self.e_variation_name.set(r.variation_name)
        self.e_combine.set(r.combine_variations)
        self.e_total_mode.set(r.total_mode)
        self.e_sym.set(r.symmetrize)
        self.e_barlow.set(r.apply_barlow)
        self.e_pair_by.set(r.pair_by)
        self.e_notes.set(r.notes)

    def editor_to_rule(self) -> Rule:
        return Rule(
            include=bool(self.e_include.get()),
            name=self.e_name.get().strip(),
            source_group=self.e_source_group.get().strip(),
            source_name=self.e_source_name.get().strip(),
            variation_name=self.e_variation_name.get().strip(),
            combine_variations=self.e_combine.get().strip() or "envelope",
            total_mode=self.e_total_mode.get().strip() or "quadrature",
            symmetrize=self.e_sym.get().strip() or "none",
            apply_barlow=bool(self.e_barlow.get()),
            pair_by=self.e_pair_by.get().strip(),
            notes=self.e_notes.get().strip(),
        )

    def apply_editor_to_selected(self):
        idx = self.selected_index()
        if idx is None:
            messagebox.showinfo("Info", "Select a rule first.")
            return
        r = self.editor_to_rule()
        if not r.name:
            messagebox.showerror("Invalid rule", "Rule name cannot be empty.")
            return
        self.rules[idx] = r
        self.refresh_tree()
        self.tree.selection_set(str(idx))

    def double_click_toggle(self, event=None):
        self.toggle_include()

    def toggle_include(self):
        idx = self.selected_index()
        if idx is None:
            return
        self.rules[idx].include = not self.rules[idx].include
        self.refresh_tree()
        self.tree.selection_set(str(idx))

    def toggle_barlow(self):
        idx = self.selected_index()
        if idx is None:
            return
        self.rules[idx].apply_barlow = not self.rules[idx].apply_barlow
        self.refresh_tree()
        self.tree.selection_set(str(idx))

    def add_rule(self):
        base = self.editor_to_rule()
        if not base.name:
            base.name = f"rule_{len(self.rules)+1}"
        self.rules.append(base)
        self.refresh_tree()
        self.tree.selection_set(str(len(self.rules)-1))

    def duplicate_selected(self):
        idx = self.selected_index()
        if idx is None:
            return
        r = Rule(**asdict(self.rules[idx]))
        r.name = r.name + "_copy"
        self.rules.insert(idx + 1, r)
        self.refresh_tree()
        self.tree.selection_set(str(idx + 1))

    def delete_rule(self):
        idx = self.selected_index()
        if idx is None:
            return
        del self.rules[idx]
        self.refresh_tree()

    def move_rule(self, delta: int):
        idx = self.selected_index()
        if idx is None:
            return
        j = idx + delta
        if j < 0 or j >= len(self.rules):
            return
        self.rules[idx], self.rules[j] = self.rules[j], self.rules[idx]
        self.refresh_tree()
        self.tree.selection_set(str(j))

    def load_config(self):
        try:
            path = Path(self.config_var.get().strip())
            with path.open("r", encoding="utf-8") as f:
                cfg = json.load(f)
            self.loaded_config = cfg
            self.input_var.set(str(cfg.get("input", self.input_var.get())))
            self.output_tsv_var.set(str(cfg.get("output_tsv", self.output_tsv_var.get())))
            self.output_csv_var.set(str(cfg.get("output_csv", self.output_csv_var.get())))
            self.details_tsv_var.set(str(cfg.get("details_tsv", self.details_tsv_var.get())))
            self.output_root_var.set(str(cfg.get("output_root", self.output_root_var.get())))

            b = dict(cfg.get("barlow", {}))
            self.global_barlow_enabled.set(bool(b.get("enabled", True)))
            self.global_barlow_mode.set(str(b.get("mode", "correlated")))
            self.global_barlow_threshold.set(str(b.get("threshold", "1.0")))
            self.global_barlow_action.set(str(b.get("action", "subtract_stat_diff")))

            total = dict(cfg.get("total", {}))
            self.total_mode_var.set(str(total.get("combine_groups", "mixed")))
            self.default_group_mode_var.set(str(total.get("default_group_mode", "quadrature")))
            self.include_stat_var.set(bool(total.get("include_stat", True)))

            plots = dict(cfg.get("plots", {}))
            self.plots_enabled_var.set(bool(plots.get("enabled", True)))

            include = set(str(x) for x in total.get("include_groups", []))
            self.rules = []
            for g in cfg.get("groups", []):
                name = str(g.get("name", ""))
                inc = (name in include) if include else True
                self.rules.append(Rule.from_group_config(g, inc))
            self.refresh_tree()
            self.status_var.set(f"Loaded {len(self.rules)} rules from {path}")
        except Exception as e:
            messagebox.showerror("Load config failed", str(e))

    def build_config(self) -> Dict[str, Any]:
        cfg = dict(self.loaded_config) if self.loaded_config else {}
        cfg["project_dir"] = str(cfg.get("project_dir", "."))
        cfg["input"] = self.input_var.get().strip()
        cfg["output_tsv"] = self.output_tsv_var.get().strip()
        cfg["output_csv"] = self.output_csv_var.get().strip()
        cfg["details_tsv"] = self.details_tsv_var.get().strip()
        cfg["output_root"] = self.output_root_var.get().strip() or "Systematics/final_systematics_results_v9.root"
        cfg["write_root"] = True
        cfg.setdefault("root_command", "root")
        cfg["write_details"] = True

        cfg["barlow"] = {
            "enabled": bool(self.global_barlow_enabled.get()),
            "mode": self.global_barlow_mode.get().strip() or "correlated",
            "threshold": float(self.global_barlow_threshold.get().strip() or "1.0"),
            "action": self.global_barlow_action.get().strip() or "subtract_stat_diff",
            "missing_stat_action": "keep",
            "denominator_floor": 1e-30,
            "apply_to_component_types": ["variation"],
            "default_apply_to_groups": True,
        }

        cfg["groups"] = [r.to_group_config() for r in self.rules]
        cfg["total"] = {
            "combine_groups": self.total_mode_var.get().strip() or "mixed",
            "default_group_mode": self.default_group_mode_var.get().strip() or "quadrature",
            "include_stat": bool(self.include_stat_var.get()),
            "include_groups": [r.name for r in self.rules if r.include],
        }

        plots = dict(cfg.get("plots", {}))
        plots.setdefault("output_dir", "Systematics/FinalPlots")
        plots.setdefault("pdf", "final_systematics_gui_values_v4.pdf")
        plots.setdefault("write_png", True)
        plots.setdefault("figsize", [7.0, 5.0])
        plots.setdefault("dpi", 150)
        plots.setdefault("draw_syst_boxes", True)
        plots.setdefault("draw_stat_errors", True)
        plots.setdefault("draw_points", True)
        plots.setdefault("draw_unity_for_rcp", True)
        plots.setdefault("grid", True)
        plots.setdefault("legend", True)
        plots.setdefault("title_template", "{job_label}, {centrality_label}")
        plots.setdefault("stat_label", "stat.")
        plots.setdefault("syst_label", "syst.")

        # Plot-axis options for systematics_final_combiner_v9.py.
        # They are intentionally not exposed in the GUI yet; edit them directly
        # in the saved JSON config if needed.
        plots["logy"] = False
        plots["logy_spectra"] = True
        plots["force_default_logy_spectra"] = True
        plots["logy_rcp"] = False
        plots.setdefault("ylim", [])
        plots.setdefault("spectra_ylim", [])
        plots["rcp_ylim"] = [0.0, 2.0]
        plots["force_default_rcp_ylim"] = True
        plots.setdefault("ylim_by_job", {})
        plots.setdefault("ylim_by_observable", {})
        plots.setdefault("ylim_by_plot", {})
        plots.setdefault("rcp_ylim_by_job", {})

        plots["enabled"] = bool(self.plots_enabled_var.get())
        cfg["plots"] = plots
        return cfg

    def save_config(self) -> bool:
        try:
            cfg = self.build_config()
            path = Path(self.config_var.get().strip())
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("w", encoding="utf-8") as f:
                json.dump(cfg, f, indent=2, ensure_ascii=False)
                f.write("\n")
            self.status_var.set(f"Saved config: {path}")
            return True
        except Exception as e:
            messagebox.showerror("Save config failed", str(e))
            return False

    def run_command(self) -> List[str]:
        return [
            sys.executable,
            self.combiner_var.get().strip(),
            "--config",
            self.config_var.get().strip(),
            "--force-root",
        ]

    def copy_run_command(self):
        cmd = " ".join([repr(x) if " " in x else x for x in self.run_command()])
        self.clipboard_clear()
        self.clipboard_append(cmd)
        self.status_var.set("Run command copied to clipboard.")

    def run_combiner(self):
        if not self.save_config():
            return
        cmd = self.run_command()
        self.log("$ " + " ".join(cmd))
        def worker():
            try:
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
                assert proc.stdout is not None
                for line in proc.stdout:
                    self.after(0, self.log, line.rstrip("\n"))
                rc = proc.wait()
                self.after(0, self.status_var.set, f"Combiner finished with exit code {rc}")
            except Exception as e:
                self.after(0, self.log, f"[error] {e}")
                self.after(0, self.status_var.set, "Combiner failed.")
        threading.Thread(target=worker, daemon=True).start()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="GUI editor for final systematic-combiner config.")
    p.add_argument("--config", default="", help="Config JSON to load/save.")
    p.add_argument("--input", default="", help="Input long components TSV to place into the config.")
    p.add_argument("--combiner", default="", help="Non-GUI combiner script used by the Run combiner button.")
    p.add_argument("--output-root", default="", help="Output ROOT file to place into the config.")
    p.add_argument("--no-auto-load", action="store_true", help="Do not auto-load --config on startup.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    app = CombinerGui(args)
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
