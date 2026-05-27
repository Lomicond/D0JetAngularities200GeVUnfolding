#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
systematics_components_builder_v15.py

First command-line draft for collecting individual systematic components into one
long TSV table.

Main idea:
  * read one nominal ROOT file as the central reference,
  * compare named systematic variations to the nominal file,
  * keep every source as a separate row/component,
  * do NOT decide here how the sources are combined.

The output is intentionally "long format": one row = one observable, centrality,
bin and source/variation.  A later script can decide whether to combine sources
quadratically, by envelope, by correlation groups, etc.

Typical use from the main project directory:

  python3 Systematics/systematics_components_builder_v11.py \
      --write-template Systematics/systematics_components_config_v1.json

  # edit the JSON paths/labels if needed, then
  python3 Systematics/systematics_components_builder_v11.py \
      --config Systematics/systematics_components_config_v1.json

Dependencies:
  * numpy, pandas
  * uproot preferred; PyROOT fallback is supported for TH1 reading

Notes:
  * GUI iteration convention from older scripts is preserved:
      iteration = 4  -> ROOT histogram contains _it3_
  * ROOT histogram names are built as:
      {hist_prefix}_{centrality}_it{iteration-1}_{method}
    e.g. Lambda0_0_it3_ICS.
  * R_CP names are supported via prefixes like RCP_5_20_Lambda0.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import re
import sys
import traceback
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import uproot  # type: ignore
except Exception:
    uproot = None

try:
    import ROOT  # type: ignore
except Exception:
    ROOT = None


# -----------------------------------------------------------------------------
# Histogram naming
# -----------------------------------------------------------------------------

RESULT_HIST_KIND_MAP: Dict[str, str] = {
    # 1D / projected spectra
    "PT (1D)": "d0pt",
    "z": "Lambda0",
    "#lambda^{1}_{1}": "Lambda1",
    "#lambda^{1}_{1.5}": "Lambda2",
    "#lambda^{1}_{2}": "Lambda3",
    "#lambda^{1}_{3}": "Lambda4",
    "#lambda^{1}_{0.5}": "Lambda5",
    "P_{T}^{D}": "Lambda6",
    # 2D-result projections, if you want to read them directly
    "PT, z": "d0ptLambda0",
    "PT, #lambda^{1}_{1}": "d0ptLambda1",
    "PT, #lambda^{1}_{1.5}": "d0ptLambda2",
    "PT, #lambda^{1}_{2}": "d0ptLambda3",
    "PT, #lambda^{1}_{3}": "d0ptLambda4",
    "PT, #lambda^{1}_{0.5}": "d0ptLambda5",
    "PT, P_{T}^{D}": "d0ptLambda6",
    # R_CP for jet pT range 5--20 GeV/c; extend in the JSON if needed
    "RCP 5-20 z": "RCP_5_20_Lambda0",
    "RCP 5-20 #lambda^{1}_{1}": "RCP_5_20_Lambda1",
    "RCP 5-20 #lambda^{1}_{1.5}": "RCP_5_20_Lambda2",
    "RCP 5-20 #lambda^{1}_{2}": "RCP_5_20_Lambda3",
    "RCP 5-20 #lambda^{1}_{3}": "RCP_5_20_Lambda4",
    "RCP 5-20 #lambda^{1}_{0.5}": "RCP_5_20_Lambda5",
    "RCP 5-20 P_{T}^{D}": "RCP_5_20_Lambda6",
}

CENT_LABELS = {
    0: "0-10%",
    1: "10-40%",
    2: "40-80%",
}

RCP_LABELS = {
    0: "0-10%/10-40%",
    1: "0-10%/40-80%",
    2: "10-40%/40-80%",
}


def centrality_label_for(observable: str, cent: int) -> str:
    """Return a label for the centrality-like index.

    For ordinary spectra, the index is a centrality bin.
    For R_CP histograms, the same index labels the numerator/denominator pair:
      0 -> 0-10 / 10-40
      1 -> 0-10 / 40-80
      2 -> 10-40 / 40-80
    """
    obs = str(observable or "")
    if obs.startswith("RCP") or obs.startswith("R_{CP}") or "RCP" in obs:
        return RCP_LABELS.get(int(cent), str(cent))
    return CENT_LABELS.get(int(cent), str(cent))

DEFAULT_OUTPUT_COLUMNS = [
    "job_label",
    "observable",
    "observable_pretty",
    "hist_name",
    "centrality",
    "centrality_label",
    "method",
    "iteration_display",
    "iteration_root",
    "bin",
    "bin_low",
    "bin_high",
    "nominal_value",
    "nominal_stat_abs",
    "nominal_stat_pct",
    "source_group",
    "source_name",
    "variation_name",
    "variation_code",
    "component_type",
    "direction",
    "variation_value",
    "variation_stat_abs",
    "diff",
    "diff_pct",
    "abs_pct",
    "n_variations",
    "reference_file",
    "variation_file",
    "run_id",
    "notes",
]

METRIC_COLS = [
    "worst_bin_pct",
    "mean_abs_drift_pct",
    "rms_drift_pct",
    "weighted_drift_pct",
    "unfolded_to_mc_pct",
]

STABILITY_BASE_COLS = [
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


# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------

def warn(msg: str) -> None:
    print(f"[warning] {msg}", file=sys.stderr)


def info(msg: str) -> None:
    print(f"[info] {msg}")


def maybe_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        if isinstance(x, str) and not x.strip():
            return None
        return float(x)
    except Exception:
        return None


def is_finite_number(x: Any) -> bool:
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def normalize_run_id(run_text: Any) -> str:
    txt = str(run_text or "").strip()
    if not txt:
        return ""
    if txt.startswith("r") and txt[1:].isdigit():
        return f"r{int(txt[1:]):06d}"
    if txt.isdigit():
        return f"r{int(txt):06d}"
    m = re.search(r"(\d+)$", txt)
    if m:
        return f"r{int(m.group(1)):06d}"
    return txt


def normalize_code_text(value: Any) -> str:
    """Return stable text for codes loaded from TSV/JSON, e.g. 0, 0.0 -> '0'."""
    txt = str(value if value is not None else "").strip()
    if txt == "" or txt.lower() in ("nan", "none", "<na>"):
        return ""
    try:
        f = float(txt)
        if math.isfinite(f) and abs(f - round(f)) < 1e-12:
            return str(int(round(f)))
    except Exception:
        pass
    return txt




PRIOR_SHAPE_CATEGORY_LABELS = {
    0: "jet_pt",
    1: "z",
    2: "lambda11",
    3: "lambda15",
    4: "lambda21",
    5: "lambda31",
    6: "lambda051",
    7: "pTD",
}

PRIOR_SHAPE_SECOND_CATEGORY_BY_RESULT = {
    "z": 1,
    "RCP 5-20 z": 1,
    "#lambda^{1}_{1}": 2,
    "RCP 5-20 #lambda^{1}_{1}": 2,
    "#lambda^{1}_{1.5}": 3,
    "RCP 5-20 #lambda^{1}_{1.5}": 3,
    "#lambda^{1}_{2}": 4,
    "RCP 5-20 #lambda^{1}_{2}": 4,
    "#lambda^{1}_{3}": 5,
    "RCP 5-20 #lambda^{1}_{3}": 5,
    "#lambda^{1}_{0.5}": 6,
    "RCP 5-20 #lambda^{1}_{0.5}": 6,
    "P_{T}^{D}": 7,
    "RCP 5-20 P_{T}^{D}": 7,
}


def decode_prior_shape_code(code_text: str, scheme: str = "zero_based") -> Tuple[Optional[int], Optional[int]]:
    """
    Decode usePriorShapeWeighting code into (variable_category, sign_digit).

    Default zero_based scheme follows the convention used in our current macros:
      code = 10 * variable_index + sign_digit
      variable_index: 0=jet pT, 1=z, 2=lambda11, ..., 7=pTD

    The alternative one_based scheme is kept as a config switch in case an older
    scan used 1=jet pT, 2=z, ... as the tens category.
    """
    txt = normalize_code_text(code_text)
    if not txt:
        return None, None
    try:
        code = int(txt)
    except Exception:
        return None, None

    sign_digit = abs(code) % 10
    tens = abs(code) // 10
    scheme = str(scheme or "zero_based").strip().lower()
    if scheme in ("one_based", "1_based", "one-based"):
        variable_category = tens - 1
    else:
        variable_category = tens
    return variable_category, sign_digit


def prior_shape_label_from_code(
    code_text: str,
    result_kind: str,
    scheme: str = "zero_based",
    sign_map: Optional[Dict[str, str]] = None,
    filter_relevant: bool = True,
) -> Optional[Tuple[str, str, str]]:
    """
    Return (source_name, variation_name, notes) for prior-shape variations.

    source_name is intentionally reduced to either prior_shape_first_variable or
    prior_shape_second_variable so the later combination script can group the two
    +/- tilts into the source requested by the analysis note.

    Supported code schemes:
      * zero_based:       10*variable_index + sign, 0=jet pT, 1=z, ...
      * one_based:        10*(variable_index+1) + sign
      * first_second:     10/11 = jet pT +/-; 20/21 = second variable +/-
                          This matches scanPriorShape summary.tsv with labels
                          nominal, jetPt_plus20, jetPt_minus20,
                          secondVar_plus20, secondVar_minus20.
    """
    scheme_norm = str(scheme or "zero_based").strip().lower()
    second_category = PRIOR_SHAPE_SECOND_CATEGORY_BY_RESULT.get(str(result_kind))

    if scheme_norm in (
        "first_second",
        "first-second",
        "two_slot",
        "two-slot",
        "jetpt_secondvar",
        "jetpt-secondvar",
        "scan_prior_shape",
        "scanpriorshape",
    ):
        txt = normalize_code_text(code_text)
        if not txt:
            return None
        try:
            code = int(txt)
        except Exception:
            return None
        sign_digit = abs(code) % 10
        slot = abs(code) // 10
        if slot == 1:
            variable_category = 0
        elif slot == 2:
            if second_category is None:
                return None
            variable_category = second_category
        else:
            return None
    else:
        variable_category, sign_digit = decode_prior_shape_code(code_text, scheme_norm)
        if variable_category is None or sign_digit is None:
            return None

    if filter_relevant and variable_category not in (0, second_category):
        return None

    if sign_map is None:
        sign_map = {
            "0": "plus20",
            "1": "minus20",
            "2": "minus20",
            "+": "plus20",
            "-": "minus20",
        }
    sign_label = sign_map.get(str(sign_digit), f"sign{sign_digit}")
    var_label = PRIOR_SHAPE_CATEGORY_LABELS.get(variable_category, f"var{variable_category}")

    if variable_category == 0:
        source_name = "prior_shape_first_variable"
    elif second_category is not None and variable_category == second_category:
        source_name = "prior_shape_second_variable"
    else:
        source_name = f"prior_shape_{var_label}"

    variation_name = f"{var_label}_{sign_label}"
    notes = f"auto-decoded prior shape code {code_text}: variable={var_label}, sign={sign_label}, scheme={scheme_norm}"
    return source_name, variation_name, notes


def d0meson_label_from_paper_label(label_text: str) -> Optional[Tuple[str, str, str]]:
    """Return (source_name, variation_name, notes) for scanPaperSys paper_label values."""
    raw = str(label_text or "").strip()
    if not raw or raw.lower() == "nominal":
        return None

    m = re.match(r"^paper_(.+)_(up|down)$", raw)
    if not m:
        return None

    core, sign = m.group(1), m.group(2)
    source_map = {
        "tpc_track": "tpc_track",
        "pid": "pid",
        "single_track_pt": "single_track",
        "topo_eff": "topo_eff",
        "double_counting": "double_counting",
        "vertex_corr": "vertex_correction",
        "secondary_track": "secondary_track",
    }
    variation_core_map = {
        "tpc_track": "tpc_track",
        "pid": "pid",
        "single_track_pt": "single_track_pt",
        "topo_eff": "topo_eff",
        "double_counting": "double_counting",
        "vertex_corr": "vertex_correction",
        "secondary_track": "secondary_track",
    }

    source_name = source_map.get(core, core)
    variation_core = variation_core_map.get(core, source_name)
    if core == "topo_eff":
        variation_name = "topo_eff_150pct" if sign == "up" else "topo_eff_50pct"
    else:
        variation_name = f"{variation_core}_{sign}"

    notes = f"auto-decoded D0Meson paper_label {raw}: source={source_name}, variation={variation_name}"
    return source_name, variation_name, notes


def run_id_to_int(run_id: Any) -> int:
    m = re.search(r"(\d+)$", str(run_id))
    return int(m.group(1)) if m else -1


def as_list(x: Any) -> List[Any]:
    if x is None:
        return []
    if isinstance(x, list):
        return x
    return [x]


def read_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: Path, data: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")


def resolve_path(path: str, project_dir: Path) -> Path:
    p = Path(os.path.expandvars(os.path.expanduser(str(path))))
    if not p.is_absolute():
        p = project_dir / p
    return p


def expand_template(text: str, context: Dict[str, Any]) -> str:
    """
    Replace only named placeholders present in *context*.

    This deliberately does NOT use str.format/format_map, because ROOT/LaTeX
    labels such as ``#lambda^{1}_{1}`` contain braces with positional fields
    (``{1}``) and would otherwise raise:
        ValueError: Format string contains positional fields

    Supported placeholders are simple named fields such as {scan_dir}, {cent},
    {job_label}, {base_job_label}, {run_id}, ... . Unknown placeholders are
    left untouched.
    """
    if text is None:
        return ""

    def repl(match: re.Match) -> str:
        key = match.group(1)
        if key in context:
            val = context.get(key)
            return "" if val is None else str(val)
        return match.group(0)

    return re.sub(r"\{([A-Za-z_][A-Za-z0-9_]*)\}", repl, str(text))


def first_existing_glob(pattern: str, project_dir: Path) -> Optional[Path]:
    ptn = str(resolve_path(pattern, project_dir))
    matches = sorted(glob.glob(ptn))
    if not matches:
        return None
    # Prefer final spectra-like names over auxiliary ROOT files.
    preferred = []
    for m in matches:
        name = Path(m).name
        score = 0
        if name.startswith("OutputSpectra"):
            score -= 20
        if "OutputSpectraOutput" in name:
            score -= 10
        if "root" in name:
            score -= 1
        preferred.append((score, m))
    preferred.sort(key=lambda x: (x[0], x[1]))
    return Path(preferred[0][1])


def find_root_file_from_pattern(
    pattern: str,
    project_dir: Path,
    context: Dict[str, Any],
    required: bool = False,
) -> Optional[Path]:
    """
    Resolve either an exact path or a glob/pattern path.
    Supports {run_id}, {scan_dir}, {source_name}, etc. from context.
    """
    expanded = expand_template(pattern, context)
    p = resolve_path(expanded, project_dir)
    if p.exists():
        return p
    if any(ch in str(p) for ch in "*?["):
        found = first_existing_glob(str(p), Path("/")) if p.is_absolute() else first_existing_glob(str(p), project_dir)
        if found is not None:
            return found
    if required:
        raise FileNotFoundError(f"ROOT file/pattern not found: {expanded}")
    return None


def build_hist_name(result_kind: str, cent: int, iteration_display: int, method: str, extra_map: Optional[Dict[str, str]] = None) -> str:
    mapping = dict(RESULT_HIST_KIND_MAP)
    if extra_map:
        mapping.update(extra_map)
    if result_kind in mapping:
        prefix = mapping[result_kind]
    else:
        # Treat unknown result_kind as direct ROOT prefix.
        prefix = result_kind
    if int(iteration_display) < 1:
        raise ValueError("iteration_display is 1-based: e.g. 4 loads *_it3_*.")
    iteration_root = int(iteration_display) - 1
    return f"{prefix}_{int(cent)}_it{iteration_root}_{method}"


def centers_from_edges(edges: np.ndarray) -> np.ndarray:
    return 0.5 * (edges[:-1] + edges[1:])


def rebin_density_to_edges(old_edges: np.ndarray, old_values: np.ndarray, new_edges: np.ndarray) -> np.ndarray:
    """
    Map density-like histogram bin contents to another binning by overlap.
    This is the same convention as the older binning-systematics tool.
    """
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


def map_to_reference_edges(old_edges: np.ndarray, old_values: np.ndarray, ref_edges: np.ndarray) -> np.ndarray:
    if len(old_edges) == len(ref_edges) and np.allclose(old_edges, ref_edges, rtol=0, atol=1e-12):
        return np.asarray(old_values, dtype=float)
    return rebin_density_to_edges(old_edges, old_values, ref_edges)


# -----------------------------------------------------------------------------
# ROOT histogram reading
# -----------------------------------------------------------------------------

def _uproot_find_key(file_obj: Any, hist_name: str) -> Optional[str]:
    if hist_name in file_obj:
        return hist_name
    prefix = hist_name + ";"
    for key in file_obj.keys():
        if str(key).startswith(prefix):
            return str(key)
    return None


def read_hist_arrays(root_path: Path, hist_name: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return edges, values, errors for a TH1-like object."""
    root_path = Path(root_path)
    if not root_path.exists():
        raise FileNotFoundError(f"Missing ROOT file: {root_path}")

    if uproot is not None:
        with uproot.open(str(root_path)) as f:
            key = _uproot_find_key(f, hist_name)
            if key is None:
                raise KeyError(f"Histogram {hist_name!r} not found in {root_path}")
            h = f[key]
            values = np.asarray(h.values(flow=False), dtype=float)
            edges = np.asarray(h.axis().edges(), dtype=float)
            try:
                variances = h.variances(flow=False)
                if variances is None:
                    errors = np.full_like(values, np.nan, dtype=float)
                else:
                    errors = np.sqrt(np.clip(np.asarray(variances, dtype=float), 0, None))
            except Exception:
                errors = np.full_like(values, np.nan, dtype=float)
            return edges, values, errors

    if ROOT is not None:
        f = ROOT.TFile.Open(str(root_path), "READ")
        if not f or f.IsZombie():
            raise OSError(f"Cannot open ROOT file: {root_path}")
        try:
            obj = f.Get(hist_name)
            if obj is None:
                raise KeyError(f"Histogram {hist_name!r} not found in {root_path}")
            nb = obj.GetNbinsX()
            edges = np.array([obj.GetBinLowEdge(i + 1) for i in range(nb)] + [obj.GetBinLowEdge(nb + 1)], dtype=float)
            values = np.array([obj.GetBinContent(i + 1) for i in range(nb)], dtype=float)
            errors = np.array([obj.GetBinError(i + 1) for i in range(nb)], dtype=float)
            return edges, values, errors
        finally:
            f.Close()

    raise ImportError("Neither uproot nor PyROOT is available. Install uproot or run in a ROOT/PyROOT environment.")


@dataclass
class HistData:
    path: Path
    hist_name: str
    edges: np.ndarray
    values: np.ndarray
    errors: np.ndarray


def load_hist(path: Path, hist_name: str) -> HistData:
    edges, values, errors = read_hist_arrays(path, hist_name)
    return HistData(path=path, hist_name=hist_name, edges=edges, values=values, errors=errors)


# -----------------------------------------------------------------------------
# Output row construction
# -----------------------------------------------------------------------------

def pct(numer: float, denom: float) -> float:
    if not is_finite_number(numer) or not is_finite_number(denom) or float(denom) == 0.0:
        return np.nan
    return 100.0 * float(numer) / float(denom)


def make_row(
    *,
    job_label: str,
    observable: str,
    observable_pretty: str,
    hist_name: str,
    cent: int,
    method: str,
    iteration_display: int,
    bin_idx: int,
    bin_low: float,
    bin_high: float,
    nominal_value: float,
    nominal_stat_abs: float,
    source_group: str,
    source_name: str,
    variation_name: str,
    component_type: str,
    direction: str,
    variation_value: float = np.nan,
    variation_stat_abs: float = np.nan,
    variation_code: Any = "",
    diff: float = np.nan,
    diff_pct: float = np.nan,
    abs_pct: float = np.nan,
    n_variations: Any = "",
    reference_file: Any = "",
    variation_file: Any = "",
    run_id: Any = "",
    notes: str = "",
) -> Dict[str, Any]:
    stat_pct = pct(nominal_stat_abs, nominal_value)
    if not is_finite_number(abs_pct) and is_finite_number(diff_pct):
        abs_pct = abs(float(diff_pct))
    return {
        "job_label": job_label,
        "observable": observable,
        "observable_pretty": observable_pretty,
        "hist_name": hist_name,
        "centrality": int(cent),
        "centrality_label": centrality_label_for(observable, int(cent)),
        "method": method,
        "iteration_display": int(iteration_display),
        "iteration_root": int(iteration_display) - 1,
        "bin": int(bin_idx),
        "bin_low": float(bin_low),
        "bin_high": float(bin_high),
        "nominal_value": float(nominal_value) if is_finite_number(nominal_value) else np.nan,
        "nominal_stat_abs": float(nominal_stat_abs) if is_finite_number(nominal_stat_abs) else np.nan,
        "nominal_stat_pct": float(stat_pct) if is_finite_number(stat_pct) else np.nan,
        "source_group": source_group,
        "source_name": source_name,
        "variation_name": variation_name,
        "variation_code": variation_code,
        "component_type": component_type,
        "direction": direction,
        "variation_value": float(variation_value) if is_finite_number(variation_value) else np.nan,
        "variation_stat_abs": float(variation_stat_abs) if is_finite_number(variation_stat_abs) else np.nan,
        "diff": float(diff) if is_finite_number(diff) else np.nan,
        "diff_pct": float(diff_pct) if is_finite_number(diff_pct) else np.nan,
        "abs_pct": float(abs_pct) if is_finite_number(abs_pct) else np.nan,
        "n_variations": n_variations,
        "reference_file": str(reference_file),
        "variation_file": str(variation_file),
        "run_id": run_id,
        "notes": notes,
    }


# -----------------------------------------------------------------------------
# Config model helpers
# -----------------------------------------------------------------------------

@dataclass
class Job:
    label: str
    result_kind: str
    centralities: List[int] = field(default_factory=lambda: [0, 1, 2])
    method: str = "ICS"
    iteration: int = 4
    observable_pretty: str = ""
    enabled: bool = True

    @staticmethod
    def from_dict(d: Dict[str, Any], global_method: str, global_iteration: int) -> "Job":
        return Job(
            label=str(d.get("label") or d.get("result_kind") or "job"),
            result_kind=str(d.get("result_kind") or d.get("observable") or d.get("label")),
            centralities=[int(c) for c in as_list(d.get("centralities", [0, 1, 2]))],
            method=str(d.get("method", global_method)),
            iteration=int(d.get("iteration", global_iteration)),
            observable_pretty=str(d.get("observable_pretty", d.get("result_kind", d.get("label", "")))),
            enabled=bool(d.get("enabled", True)),
        )


@dataclass
class NominalReference:
    file: Path
    hist: HistData


# -----------------------------------------------------------------------------
# Summary/stability support for simple groups and binning
# -----------------------------------------------------------------------------

def load_summary_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing summary table: {path}")
    df = pd.read_csv(path, sep="\t", header=0)
    if "run_id" in df.columns:
        df["run_id"] = df["run_id"].astype(str).map(normalize_run_id)
    return df


def column_candidates(group_cfg: Dict[str, Any], singular_key: str, plural_key: str, default: str = "") -> List[str]:
    """Return ordered candidate column names from either a string or a list in JSON."""
    vals: List[str] = []
    if plural_key in group_cfg:
        vals.extend(str(x).strip() for x in as_list(group_cfg.get(plural_key)) if str(x).strip())
    if singular_key in group_cfg:
        vals.extend(str(x).strip() for x in as_list(group_cfg.get(singular_key)) if str(x).strip())
    if default:
        vals.append(default)
    # Preserve order but remove duplicates.
    out: List[str] = []
    seen = set()
    for v in vals:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def first_existing_column(row_dict: Dict[str, Any], candidates: List[str]) -> str:
    for col in candidates:
        if col in row_dict:
            return col
    return ""


def first_row_value(row_dict: Dict[str, Any], candidates: List[str], default: Any = "") -> Any:
    col = first_existing_column(row_dict, candidates)
    if not col:
        return default
    val = row_dict.get(col, default)
    if val is None:
        return default
    try:
        if pd.isna(val):
            return default
    except Exception:
        pass
    return val


def available_column_note(prefix: str, candidates: List[str], row_dict: Dict[str, Any]) -> str:
    col = first_existing_column(row_dict, candidates)
    if col:
        return f"{prefix}={col}"
    return f"{prefix}=missing({','.join(candidates)})"


def load_stability_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing stability table: {path}")
    raw = pd.read_csv(path, sep="\t", header=None)
    optional_cols = ["unfolded_to_mc_pct"]
    cols = STABILITY_BASE_COLS + optional_cols
    if raw.shape[1] < len(STABILITY_BASE_COLS):
        raise ValueError(f"{path} has {raw.shape[1]} columns; expected at least {len(STABILITY_BASE_COLS)}.")
    if raw.shape[1] > len(cols):
        raise ValueError(f"{path} has {raw.shape[1]} columns; this draft supports at most {len(cols)}.")
    raw.columns = cols[: raw.shape[1]]
    for c in cols:
        if c not in raw.columns:
            raw[c] = pd.NA
    raw["run_id"] = raw["run_id"].astype(str).map(normalize_run_id)
    raw["cent"] = pd.to_numeric(raw["cent"], errors="coerce").astype("Int64")
    for c in METRIC_COLS:
        raw[c] = pd.to_numeric(raw[c], errors="coerce")
    raw = raw.dropna(subset=["cent"]).copy()
    raw["cent"] = raw["cent"].astype(int)
    return raw


def apply_exact_filter(df: pd.DataFrame, col: str, value: Any) -> pd.DataFrame:
    if value is None or str(value).strip() == "":
        return df
    if col not in df.columns:
        raise KeyError(f"Filter column {col!r} not found. Available: {list(df.columns)}")
    value_s = str(value).strip()
    value_f = maybe_float(value_s)
    col_num = pd.to_numeric(df[col], errors="coerce")
    if value_f is not None and col_num.notna().any():
        return df[np.isclose(col_num.to_numpy(dtype=float), value_f, rtol=0.0, atol=1e-12)]
    return df[df[col].astype(str) == value_s]


def filter_stability_for_source(df: pd.DataFrame, selection: Dict[str, Any], context: Dict[str, Any]) -> pd.DataFrame:
    out = df.copy()
    # Common exact filters.
    for col in ["cent", "dim", "observable", "axis"]:
        if col in selection:
            val = expand_template(str(selection.get(col, "")), context)
            out = apply_exact_filter(out, col, val)

    # Additional exact filters.
    for col, raw_val in selection.get("exact", {}).items():
        val = expand_template(str(raw_val), context)
        out = apply_exact_filter(out, col, val)

    # Optional metric lower/upper cuts.
    min_col = str(selection.get("metric_min_col", "") or "").strip()
    min_val = str(selection.get("metric_min_value", "") or "").strip()
    max_col = str(selection.get("metric_max_col", "") or "").strip()
    max_val = str(selection.get("metric_max_value", "") or "").strip()

    if min_col and min_val:
        if min_col not in METRIC_COLS:
            raise ValueError(f"Unknown metric_min_col {min_col!r}")
        out = out[pd.to_numeric(out[min_col], errors="coerce") >= float(min_val)]
    if max_col and max_val:
        if max_col not in METRIC_COLS:
            raise ValueError(f"Unknown metric_max_col {max_col!r}")
        out = out[pd.to_numeric(out[max_col], errors="coerce") <= float(max_val)]

    return out.reset_index(drop=True)


def merge_nested_dict(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge dictionaries without modifying the inputs."""
    out = dict(base or {})
    for key, value in dict(override or {}).items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = merge_nested_dict(out[key], value)
        else:
            out[key] = value
    return out


def value_matches_pattern(value: Any, pattern: Any, context: Dict[str, Any]) -> bool:
    """Match helper for binning override rules. Supports scalars, lists and '*'."""
    if pattern is None:
        return True
    if isinstance(pattern, list):
        return any(value_matches_pattern(value, item, context) for item in pattern)
    pat = expand_template(str(pattern), context).strip()
    if pat in ("", "*"):
        return True
    val = str(value).strip()

    pat_f = maybe_float(pat)
    val_f = maybe_float(val)
    if pat_f is not None and val_f is not None:
        return abs(pat_f - val_f) < 1e-12
    return val == pat


def binning_rule_matches(rule: Dict[str, Any], context: Dict[str, Any]) -> bool:
    """
    Decide whether a binning metric-cut override applies.

    Supported match keys:
      job_label / label, result_kind / observable, centrality / cent, source_name.
    Missing match block means the rule applies to everything.
    """
    match = dict(rule.get("match", {}))
    if not match:
        # Also allow short form: {"job_label": ..., "centrality": ..., "metric_cuts": ...}
        for key in ["job_label", "label", "result_kind", "observable", "centrality", "cent", "source_name"]:
            if key in rule:
                match[key] = rule[key]
    for key, pattern in match.items():
        if key in ("job_label", "label"):
            value = context.get("job_label", "")
        elif key in ("result_kind", "observable"):
            value = context.get("result_kind", "")
        elif key in ("centrality", "cent"):
            value = context.get("cent", "")
        elif key == "source_name":
            value = context.get("source_name", "")
        else:
            # Unknown match keys are deliberately strict: typo should not silently match.
            return False
        if not value_matches_pattern(value, pattern, context):
            return False
    return True


def normalize_metric_cut_value(value: Any, context: Dict[str, Any]) -> str:
    if value is None:
        return ""
    return expand_template(str(value), context).strip()


def metric_cuts_from_legacy_selection(selection: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """Convert old metric_min_col/value + metric_max_col/value syntax to metric_cuts."""
    cuts: Dict[str, Dict[str, Any]] = {}
    min_col = str(selection.get("metric_min_col", "") or "").strip()
    min_val = str(selection.get("metric_min_value", "") or "").strip()
    max_col = str(selection.get("metric_max_col", "") or "").strip()
    max_val = str(selection.get("metric_max_value", "") or "").strip()
    if min_col and min_val:
        cuts.setdefault(min_col, {})["min"] = min_val
    if max_col and max_val:
        cuts.setdefault(max_col, {})["max"] = max_val
    return cuts


def resolve_binning_metric_cuts(bcfg: Dict[str, Any], source: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """
    Build final metric cuts for this (source, job, centrality).

    Precedence, from lowest to highest:
      1) binning.default_metric_cuts
      2) source.default_metric_cuts / source.metric_cuts
      3) legacy metric cuts inside source.selection
      4) binning.metric_cut_overrides matching this job/centrality/source
      5) source.metric_cut_overrides matching this job/centrality/source

    Missing min/max values are inherited from lower-precedence defaults.
    Empty string disables that side of the cut.
    """
    cuts: Dict[str, Dict[str, Any]] = {}
    cuts = merge_nested_dict(cuts, dict(bcfg.get("default_metric_cuts", {})))
    cuts = merge_nested_dict(cuts, dict(source.get("default_metric_cuts", {})))
    cuts = merge_nested_dict(cuts, dict(source.get("metric_cuts", {})))
    cuts = merge_nested_dict(cuts, metric_cuts_from_legacy_selection(dict(source.get("selection", {}))))

    for rule in list(bcfg.get("metric_cut_overrides", [])) + list(source.get("metric_cut_overrides", [])):
        if binning_rule_matches(dict(rule), context):
            cuts = merge_nested_dict(cuts, dict(rule.get("metric_cuts", {})))

    # Normalize to {metric: {min: ..., max: ...}} and drop fully empty cuts.
    out: Dict[str, Dict[str, Any]] = {}
    for metric, limits in cuts.items():
        metric = str(metric).strip()
        if not metric:
            continue
        if metric not in METRIC_COLS:
            raise ValueError(f"Unknown binning metric {metric!r}. Allowed: {METRIC_COLS}")
        if isinstance(limits, dict):
            vmin = normalize_metric_cut_value(limits.get("min", ""), context)
            vmax = normalize_metric_cut_value(limits.get("max", ""), context)
        else:
            # Short form: "metric": 5.0 means max <= 5.0.
            vmin = ""
            vmax = normalize_metric_cut_value(limits, context)
        if vmin or vmax:
            out[metric] = {"min": vmin, "max": vmax}
    return out


def apply_binning_metric_cuts(df: pd.DataFrame, cuts: Dict[str, Dict[str, Any]]) -> pd.DataFrame:
    out = df.copy()
    for metric, limits in cuts.items():
        vals = pd.to_numeric(out[metric], errors="coerce")
        vmin = str(limits.get("min", "") or "").strip()
        vmax = str(limits.get("max", "") or "").strip()
        if vmin:
            out = out[vals >= float(vmin)]
            vals = pd.to_numeric(out[metric], errors="coerce")
        if vmax:
            out = out[vals <= float(vmax)]
    return out.reset_index(drop=True)


def metric_cuts_to_note(cuts: Dict[str, Dict[str, Any]]) -> str:
    parts: List[str] = []
    for metric in sorted(cuts.keys()):
        limits = cuts[metric]
        vmin = str(limits.get("min", "") or "").strip()
        vmax = str(limits.get("max", "") or "").strip()
        if vmin:
            parts.append(f"{metric}>={vmin}")
        if vmax:
            parts.append(f"{metric}<={vmax}")
    return ", ".join(parts) if parts else "no metric cuts"


def extract_first_last_number(text: Any) -> Tuple[Optional[float], Optional[float]]:
    """Extract first and last numeric value from an edge-list-like string."""
    txt = str(text or "")
    nums = re.findall(r"[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?", txt)
    if len(nums) < 2:
        return None, None
    try:
        return float(nums[0]), float(nums[-1])
    except Exception:
        return None, None


def series_has_nonempty_values(s: pd.Series) -> bool:
    try:
        vals = s.dropna().astype(str).str.strip()
        return bool(vals.ne("").any())
    except Exception:
        return False


def first_nonempty_col(df: pd.DataFrame, candidates: List[str]) -> str:
    for c in candidates:
        if c in df.columns and series_has_nonempty_values(df[c]):
            return c
    return ""



def _source_prefers_second_axis(source_name: str) -> bool:
    txt = str(source_name or "").lower()
    return any(k in txt for k in ("second", "sec", "var")) and "first" not in txt


def ensure_binning_endpoint_columns(df: pd.DataFrame, level: str, source_name: str = "") -> pd.DataFrame:
    """
    Add canonical endpoint columns for true/reco binning filters.

    The resulting columns are:
      _true_edge_min, _true_edge_max
      _reco_edge_min, _reco_edge_max

    The helper tries several summary.tsv conventions.  If the source is the
    second variable, sec/var columns are preferred; if it is the first variable,
    main/pt columns are preferred.  This avoids accidentally applying the PT
    range to the second variable when both are present in a scan summary.
    """
    level = str(level).strip().lower()
    if level not in ("true", "reco"):
        raise ValueError(f"Unknown endpoint level {level!r}; expected 'true' or 'reco'.")

    out = df.copy()
    c_out_min = f"_{level}_edge_min"
    c_out_max = f"_{level}_edge_max"
    if c_out_min not in out.columns:
        out[c_out_min] = np.nan
    if c_out_max not in out.columns:
        out[c_out_max] = np.nan

    cols = list(out.columns)
    lower = {c: c.lower() for c in cols}
    second = _source_prefers_second_axis(source_name)

    if level == "true":
        main_min = ["main_true_min", "true_pt_min", "truth_pt_min", "mc_pt_min", "pt_true_min", "true_x_min", "x_true_min"]
        main_max = ["main_true_max", "true_pt_max", "truth_pt_max", "mc_pt_max", "pt_true_max", "true_x_max", "x_true_max"]
        var_min = ["sec_true_min", "second_true_min", "true_var_min", "var_true_min", "true_y_min", "y_true_min"]
        var_max = ["sec_true_max", "second_true_max", "true_var_max", "var_true_max", "true_y_max", "y_true_max"]
        generic_min = ["true_min", "truth_min", "mc_min"]
        generic_max = ["true_max", "truth_max", "mc_max"]
        level_like_words = ("true", "truth", "mc")
        main_edges = ["main_true_edges", "true_pt_edges", "pt_true_edges"]
        var_edges = ["sec_true_edges", "second_true_edges", "true_var_edges", "var_true_edges", "true_y_edges", "y_true_edges"]
        generic_edges = ["true_edges", "truth_edges", "mc_edges", "true_bins", "truth_bins", "mc_bins"]
    else:
        main_min = ["main_reco_min", "reco_pt_min", "rec_pt_min", "meas_pt_min", "measured_pt_min", "det_pt_min", "pt_reco_min", "reco_x_min", "x_reco_min"]
        main_max = ["main_reco_max", "reco_pt_max", "rec_pt_max", "meas_pt_max", "measured_pt_max", "det_pt_max", "pt_reco_max", "reco_x_max", "x_reco_max"]
        var_min = ["sec_reco_min", "second_reco_min", "reco_var_min", "var_reco_min", "reco_y_min", "y_reco_min"]
        var_max = ["sec_reco_max", "second_reco_max", "reco_var_max", "var_reco_max", "reco_y_max", "y_reco_max"]
        generic_min = ["reco_min", "rec_min", "meas_min", "measured_min", "det_min"]
        generic_max = ["reco_max", "rec_max", "meas_max", "measured_max", "det_max"]
        level_like_words = ("reco", "rec", "meas", "measured", "det")
        main_edges = ["main_reco_edges", "reco_pt_edges", "pt_reco_edges", "rec_pt_edges", "meas_pt_edges", "measured_pt_edges", "det_pt_edges"]
        var_edges = ["sec_reco_edges", "second_reco_edges", "reco_var_edges", "var_reco_edges", "reco_y_edges", "y_reco_edges"]
        generic_edges = ["reco_edges", "rec_edges", "meas_edges", "measured_edges", "det_edges", "reco_bins", "rec_bins", "meas_bins", "det_bins"]

    min_candidates = (var_min + main_min if second else main_min + var_min) + generic_min
    max_candidates = (var_max + main_max if second else main_max + var_max) + generic_max

    # Add any additional obvious columns from the summary.
    for c, lc in lower.items():
        level_like = any(w in lc for w in level_like_words)
        if level_like and ("min" in lc or "start" in lc) and c not in min_candidates:
            # Source-specific columns first, generic later.
            if second and ("var" in lc or "sec" in lc or "second" in lc or "_y" in lc):
                min_candidates.insert(0, c)
            elif (not second) and ("pt" in lc or "main" in lc or "_x" in lc):
                min_candidates.insert(0, c)
            else:
                min_candidates.append(c)
        if level_like and "max" in lc and c not in max_candidates:
            if second and ("var" in lc or "sec" in lc or "second" in lc or "_y" in lc):
                max_candidates.insert(0, c)
            elif (not second) and ("pt" in lc or "main" in lc or "_x" in lc):
                max_candidates.insert(0, c)
            else:
                max_candidates.append(c)

    cmin = first_nonempty_col(out, min_candidates)
    cmax = first_nonempty_col(out, max_candidates)
    if cmin:
        out[c_out_min] = pd.to_numeric(out[cmin], errors="coerce")
    if cmax:
        out[c_out_max] = pd.to_numeric(out[cmax], errors="coerce")

    need_edges = (not pd.to_numeric(out[c_out_min], errors="coerce").notna().any()) or (not pd.to_numeric(out[c_out_max], errors="coerce").notna().any())
    if need_edges:
        edge_candidates = (var_edges + main_edges if second else main_edges + var_edges) + generic_edges
        for c, lc in lower.items():
            level_like = any(w in lc for w in level_like_words)
            if level_like and ("edge" in lc or "bin" in lc) and c not in edge_candidates:
                if second and ("var" in lc or "sec" in lc or "second" in lc or "_y" in lc):
                    edge_candidates.insert(0, c)
                elif (not second) and ("pt" in lc or "main" in lc or "_x" in lc):
                    edge_candidates.insert(0, c)
                else:
                    edge_candidates.append(c)
        # Last-resort fallback for older one-axis summaries.
        for c in ["edges", "bins"]:
            if c in out.columns and c not in edge_candidates:
                edge_candidates.append(c)

        for edge_col in edge_candidates:
            if edge_col not in out.columns or not series_has_nonempty_values(out[edge_col]):
                continue
            mins: List[Any] = []
            maxs: List[Any] = []
            ok_any = False
            for val in out[edge_col]:
                lo, hi = extract_first_last_number(val)
                if lo is None or hi is None:
                    mins.append(np.nan)
                    maxs.append(np.nan)
                else:
                    mins.append(lo)
                    maxs.append(hi)
                    ok_any = True
            if ok_any:
                if not pd.to_numeric(out[c_out_min], errors="coerce").notna().any():
                    out[c_out_min] = mins
                if not pd.to_numeric(out[c_out_max], errors="coerce").notna().any():
                    out[c_out_max] = maxs
                break

    return out


def ensure_binning_true_range_columns(df: pd.DataFrame) -> pd.DataFrame:
    return ensure_binning_endpoint_columns(df, "true", "")


def _range_values(raw: Any) -> List[float]:
    if raw is None:
        return []
    if isinstance(raw, (list, tuple, set)):
        vals = list(raw)
    else:
        txt = str(raw).strip()
        if not txt:
            return []
        vals = [x.strip() for x in re.split(r"[,;|]+", txt) if x.strip()]
    out: List[float] = []
    for v in vals:
        try:
            f = float(str(v).strip())
            if math.isfinite(f):
                out.append(f)
        except Exception:
            pass
    return out


def endpoint_filter_to_note(endpoint_filter: Dict[str, Any], level: str) -> str:
    level = str(level).strip().lower()
    mode = str(endpoint_filter.get("mode", "ignore") or "ignore").strip().lower()
    enabled = bool(endpoint_filter.get("enabled", False))
    if not enabled or mode in ("", "ignore", "none", "all"):
        return f"{level} range ignored"

    min_key = f"{level}_min"
    max_key = f"{level}_max"
    mins = _range_values(endpoint_filter.get(min_key, endpoint_filter.get("min", "")))
    maxs = _range_values(endpoint_filter.get(max_key, endpoint_filter.get("max", "")))
    if not mins and not maxs:
        return f"{level} range ignored"

    def fmt(vals: List[float]) -> str:
        return "/".join(f"{v:g}" for v in vals)

    if mode == "range":
        parts = []
        if mins:
            parts.append(f"{level}_min>={mins[0]:g}")
        if maxs:
            parts.append(f"{level}_max<={maxs[0]:g}")
        return ", ".join(parts)

    if mode in ("min_exact_max_le", "start_exact_max_le", "min_allowed_max_le"):
        parts = []
        if mins:
            parts.append(f"{level}_min=={fmt(mins)}")
        if maxs:
            parts.append(f"{level}_max<={maxs[0]:g}")
        return ", ".join(parts)

    parts = []
    if mins:
        op = "in" if len(mins) > 1 or mode in ("allowed", "any", "exact_any") else "=="
        parts.append(f"{level}_min{op}{fmt(mins)}")
    if maxs:
        op = "in" if len(maxs) > 1 or mode in ("allowed", "any", "exact_any") else "=="
        parts.append(f"{level}_max{op}{fmt(maxs)}")
    return ", ".join(parts)


def true_range_filter_to_note(true_filter: Dict[str, Any]) -> str:
    return endpoint_filter_to_note(true_filter, "true")


def reco_range_filter_to_note(reco_filter: Dict[str, Any]) -> str:
    return endpoint_filter_to_note(reco_filter, "reco")


def resolve_binning_endpoint_filter(bcfg: Dict[str, Any], source: Dict[str, Any], context: Dict[str, Any], level: str) -> Dict[str, Any]:
    """
    Resolve optional true/reco endpoint filter.

    Supported modes:
      ignore/all/none: no endpoint cut
      exact: min and/or max must match one value exactly
      allowed/exact_any/any: min and/or max may match any value from a list
      range: keep min >= first min value and max <= first max value
      min_exact_max_le/start_exact_max_le: keep min exactly equal to allowed min(s), max <= first max value

    Precedence:
      binning.default_<level>_range_filter -> source.<level>_range_filter ->
      binning.<level>_range_overrides/source.<level>_range_overrides.
    """
    level = str(level).strip().lower()
    filt: Dict[str, Any] = {}
    filt = merge_nested_dict(filt, dict(bcfg.get(f"default_{level}_range_filter", {})))
    filt = merge_nested_dict(filt, dict(source.get(f"{level}_range_filter", {})))

    for rule in list(bcfg.get(f"{level}_range_overrides", [])) + list(source.get(f"{level}_range_overrides", [])):
        if binning_rule_matches(dict(rule), context):
            filt = merge_nested_dict(filt, dict(rule.get(f"{level}_range_filter", {})))

    # Backwards compatibility: true_range_overrides from v12.
    if level == "true":
        for rule in list(bcfg.get("true_range_overrides", [])) + list(source.get("true_range_overrides", [])):
            if binning_rule_matches(dict(rule), context):
                filt = merge_nested_dict(filt, dict(rule.get("true_range_filter", {})))

    out: Dict[str, Any] = dict(filt)
    min_key = f"{level}_min"
    max_key = f"{level}_max"
    if "min" in out and min_key not in out:
        out[min_key] = out.get("min")
    if "max" in out and max_key not in out:
        out[max_key] = out.get("max")

    # Normalize templates inside strings/lists.
    def norm_value(v: Any) -> Any:
        if isinstance(v, list):
            return [normalize_metric_cut_value(x, context) for x in v]
        return normalize_metric_cut_value(v, context)

    out[min_key] = norm_value(out.get(min_key, ""))
    out[max_key] = norm_value(out.get(max_key, ""))
    out["mode"] = str(out.get("mode", "exact") or "exact").strip().lower()
    if "enabled" not in out:
        out["enabled"] = bool(_range_values(out.get(min_key, "")) or _range_values(out.get(max_key, "")))
    return out


def resolve_binning_true_range_filter(bcfg: Dict[str, Any], source: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    return resolve_binning_endpoint_filter(bcfg, source, context, "true")


def resolve_binning_reco_range_filter(bcfg: Dict[str, Any], source: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    return resolve_binning_endpoint_filter(bcfg, source, context, "reco")


def apply_binning_endpoint_filter(df: pd.DataFrame, endpoint_filter: Dict[str, Any], level: str, source_name: str = "") -> pd.DataFrame:
    level = str(level).strip().lower()
    enabled = bool(endpoint_filter.get("enabled", False))
    mode = str(endpoint_filter.get("mode", "ignore") or "ignore").strip().lower()
    if not enabled or mode in ("", "ignore", "none", "all"):
        return df.reset_index(drop=True)

    min_key = f"{level}_min"
    max_key = f"{level}_max"
    mins = _range_values(endpoint_filter.get(min_key, endpoint_filter.get("min", "")))
    maxs = _range_values(endpoint_filter.get(max_key, endpoint_filter.get("max", "")))
    if not mins and not maxs:
        return df.reset_index(drop=True)

    work = ensure_binning_endpoint_columns(df, level, source_name)
    cmin = f"_{level}_edge_min"
    cmax = f"_{level}_edge_max"
    if not pd.to_numeric(work[cmin], errors="coerce").notna().any() and not pd.to_numeric(work[cmax], errors="coerce").notna().any():
        policy = str(endpoint_filter.get("on_missing_columns", "error") or "error").strip().lower()
        if policy in ("ignore", "no_filter", "pass"):
            return df.reset_index(drop=True)
        raise ValueError(f"{level.capitalize()}-level range filter requested, but no {level}_min/{level}_max/{level}_edges-like columns were found in summary.tsv.")

    out = work

    def apply_allowed(current: pd.DataFrame, col: str, allowed: List[float]) -> pd.DataFrame:
        if not allowed:
            return current
        vals = pd.to_numeric(current[col], errors="coerce").to_numpy(dtype=float)
        mask = np.zeros(len(current), dtype=bool)
        for a in allowed:
            mask |= np.isclose(vals, float(a), rtol=0.0, atol=1e-12)
        return current[mask]

    if mode == "range":
        if mins:
            vals = pd.to_numeric(out[cmin], errors="coerce")
            out = out[vals >= float(mins[0])]
        if maxs:
            vals = pd.to_numeric(out[cmax], errors="coerce")
            out = out[vals <= float(maxs[0])]
    elif mode in ("min_exact_max_le", "start_exact_max_le", "min_allowed_max_le"):
        # Useful for variable scans where all starts must remain fixed
        # but the upper reco endpoint may be any value inside the nominal range.
        out = apply_allowed(out, cmin, mins)
        if maxs:
            vals = pd.to_numeric(out[cmax], errors="coerce")
            out = out[vals <= float(maxs[0]) + 1e-12]
    else:
        # exact and allowed are both implemented as membership in allowed lists.
        out = apply_allowed(out, cmin, mins)
        out = apply_allowed(out, cmax, maxs)

    return out.reset_index(drop=True)


def apply_binning_true_range_filter(df: pd.DataFrame, true_filter: Dict[str, Any], source_name: str = "") -> pd.DataFrame:
    return apply_binning_endpoint_filter(df, true_filter, "true", source_name)


def apply_binning_reco_range_filter(df: pd.DataFrame, reco_filter: Dict[str, Any], source_name: str = "") -> pd.DataFrame:
    return apply_binning_endpoint_filter(df, reco_filter, "reco", source_name)

def binning_base_job_label(job_label: str) -> str:
    """Return the underlying observable label for RCP jobs, otherwise job_label."""
    label = str(job_label or "").strip()
    m = re.match(r"^RCP_(.+?)_\d+_\d+$", label)
    if m:
        return m.group(1)
    return label


def resolve_binning_scan_dir(source: Dict[str, Any], job: Job, cent: int, context0: Dict[str, Any]) -> str:
    """Resolve scan_dir, optionally using source.scan_dir_map keyed by job/base label/result_kind."""
    scan_map = dict(source.get("scan_dir_map", {}) or {})
    base_label = binning_base_job_label(job.label)
    keys = [
        job.label,
        base_label,
        job.result_kind,
        str(cent),
        "default",
    ]
    raw = ""
    for key in keys:
        if key in scan_map:
            raw = str(scan_map[key])
            break
    if not raw:
        raw = str(source.get("scan_dir", ""))
    return expand_template(raw, {**context0, "base_job_label": base_label})


def default_binning_observable_candidates(job: Job) -> List[str]:
    """Candidate values for the stability.tsv observable column."""
    base = binning_base_job_label(job.label)
    by_base = {
        "z": ["z", "Lambda0", "PT, z", "d0ptLambda0"],
        "l11": ["#lambda^{1}_{1}", "Lambda1", "PT, #lambda^{1}_{1}", "d0ptLambda1", "l11"],
        "l15": ["#lambda^{1}_{1.5}", "Lambda2", "PT, #lambda^{1}_{1.5}", "d0ptLambda2", "l15", "l11p5"],
        "l21": ["#lambda^{1}_{2}", "Lambda3", "PT, #lambda^{1}_{2}", "d0ptLambda3", "l21", "l12"],
        "l31": ["#lambda^{1}_{3}", "Lambda4", "PT, #lambda^{1}_{3}", "d0ptLambda4", "l31", "l13"],
        "l051": ["#lambda^{1}_{0.5}", "Lambda5", "PT, #lambda^{1}_{0.5}", "d0ptLambda5", "l051", "l10p5"],
        "pTD": ["p_{T}^{D}", "P_{T}^{D}", "Lambda6", "PT, p_{T}^{D}", "PT, P_{T}^{D}", "d0ptLambda6", "pTD", "ptd"],
    }
    out: List[str] = []
    for item in by_base.get(base, []):
        if item not in out:
            out.append(item)
    for item in [job.result_kind, job.observable_pretty if hasattr(job, "observable_pretty") else "", job.label, base]:
        if item and item not in out:
            out.append(str(item))
    return out


def configured_binning_observable_candidates(source: Dict[str, Any], job: Job, context: Dict[str, Any]) -> List[str]:
    """Resolve optional source.observable_map and append robust defaults."""
    base = binning_base_job_label(job.label)
    obs_map = dict(source.get("observable_map", {}) or {})
    raw: Any = None
    for key in [job.label, base, job.result_kind, "default"]:
        if key in obs_map:
            raw = obs_map[key]
            break
    vals: List[str] = []
    if raw is not None:
        for x in as_list(raw):
            y = expand_template(str(x), context).strip()
            if y and y not in vals:
                vals.append(y)
    for x in default_binning_observable_candidates(job):
        y = expand_template(str(x), context).strip()
        if y and y not in vals:
            vals.append(y)
    return vals


def apply_auto_observable_filter(df: pd.DataFrame, source: Dict[str, Any], job: Job, context: Dict[str, Any]) -> Tuple[pd.DataFrame, str]:
    """Auto-select the first observable candidate present in stability.tsv.

    If no candidate matches, the behaviour is configurable through
    source["on_no_observable_match"]:
      * "empty"     -> select zero rows (strict/old behaviour),
      * "no_filter" -> keep all rows that survived previous filters.

    The no-filter fallback is useful for dedicated scans where the directory
    itself already identifies the varied observable, or for quick first-pass
    RMS estimates before exact metric cuts are tuned.
    """
    if "observable" not in df.columns:
        return df, "observable column absent"
    candidates = configured_binning_observable_candidates(source, job, context)
    available = set(df["observable"].dropna().astype(str).unique().tolist())
    for cand in candidates:
        if cand in available:
            return df[df["observable"].astype(str) == cand].reset_index(drop=True), cand
    if len(available) == 1:
        only = sorted(available)[0]
        return df[df["observable"].astype(str) == only].reset_index(drop=True), f"{only} (single available fallback)"

    mode = str(source.get("on_no_observable_match", "empty") or "empty").strip().lower()
    if mode in ("no_filter", "keep", "all", "ignore"):
        sample = ", ".join(sorted(available)[:8])
        return df.reset_index(drop=True), "no match -> no observable cut; tried: " + ", ".join(candidates[:8]) + (f"; available sample: {sample}" if sample else "")

    return df.iloc[0:0].copy(), "no match; tried: " + ", ".join(candidates[:12])


def default_binning_axis_candidates(source_name: str) -> List[str]:
    """Candidate values for the stability.tsv axis column."""
    sname = str(source_name or "").lower()
    if "first" in sname or "jet" in sname or "main" in sname:
        return ["PT", "pt", "pT", "JetPt", "jetPt", "jet_pt", "X", "x", "main", "first"]
    if "second" in sname or "sec" in sname:
        return ["var", "Y", "y", "secondary", "second", "observable"]
    return ["X", "Y", "x", "y"]


def configured_binning_axis_candidates(source: Dict[str, Any], context: Dict[str, Any]) -> List[str]:
    """Resolve optional source.axis_candidates and append robust defaults."""
    vals: List[str] = []
    raw = source.get("axis_candidates", None)
    if raw is not None:
        for x in as_list(raw):
            y = expand_template(str(x), context).strip()
            if y and y not in vals:
                vals.append(y)
    for x in default_binning_axis_candidates(str(context.get("source_name", ""))):
        y = expand_template(str(x), context).strip()
        if y and y not in vals:
            vals.append(y)
    return vals


def apply_auto_axis_filter(df: pd.DataFrame, source: Dict[str, Any], context: Dict[str, Any]) -> Tuple[pd.DataFrame, str]:
    """Auto-select the first axis candidate present in stability.tsv.

    If no candidate matches, source["on_no_axis_match"] controls the fallback:
      * "empty"     -> select zero rows,
      * "no_filter" -> keep all rows that survived previous filters.

    For the current binning scans, the axis labels are not guaranteed to be
    stable across older summary/stability versions, so the template uses
    no_filter and records the fallback in the notes.
    """
    if "axis" not in df.columns:
        return df, "axis column absent"
    candidates = configured_binning_axis_candidates(source, context)
    available = set(df["axis"].dropna().astype(str).unique().tolist())
    for cand in candidates:
        if cand in available:
            return df[df["axis"].astype(str) == cand].reset_index(drop=True), cand
    if len(available) == 1:
        only = sorted(available)[0]
        return df[df["axis"].astype(str) == only].reset_index(drop=True), f"{only} (single available fallback)"

    mode = str(source.get("on_no_axis_match", "empty") or "empty").strip().lower()
    if mode in ("no_filter", "keep", "all", "ignore"):
        sample = ", ".join(sorted(available)[:8])
        return df.reset_index(drop=True), "no match -> no axis cut; tried: " + ", ".join(candidates[:8]) + (f"; available sample: {sample}" if sample else "")

    return df.iloc[0:0].copy(), "no match; tried: " + ", ".join(candidates[:12])


# -----------------------------------------------------------------------------
# Component collectors
# -----------------------------------------------------------------------------

def get_nominal_file(config: Dict[str, Any], project_dir: Path) -> Path:
    nominal_cfg = config.get("nominal", {})
    if not nominal_cfg.get("enabled", True):
        raise ValueError("nominal.enabled is false; a nominal reference is required in this draft.")

    if nominal_cfg.get("file"):
        p = resolve_path(str(nominal_cfg["file"]), project_dir)
        if p.exists():
            return p
        raise FileNotFoundError(f"nominal.file does not exist: {p}")

    if nominal_cfg.get("file_glob"):
        found = first_existing_glob(str(nominal_cfg["file_glob"]), project_dir)
        if found is not None:
            return found
        raise FileNotFoundError(f"nominal.file_glob matched nothing: {nominal_cfg['file_glob']}")

    # Conservative default.
    default_globs = [
        "Systematics/Nominal/scanNominal/Output/OutputSpectra*.root",
        "Systematics/Nominal/scanNominal/Output/*.root",
    ]
    for ptn in default_globs:
        found = first_existing_glob(ptn, project_dir)
        if found is not None:
            return found
    raise FileNotFoundError("Could not find nominal ROOT file. Set nominal.file or nominal.file_glob in the JSON config.")


def load_nominal_reference(
    config: Dict[str, Any],
    project_dir: Path,
    job: Job,
    cent: int,
    extra_hist_map: Dict[str, str],
) -> NominalReference:
    nominal_file = get_nominal_file(config, project_dir)
    hist_name = build_hist_name(job.result_kind, cent, job.iteration, job.method, extra_hist_map)
    return NominalReference(file=nominal_file, hist=load_hist(nominal_file, hist_name))


def collect_nominal_rows(
    config: Dict[str, Any],
    project_dir: Path,
    job: Job,
    cent: int,
    nominal: NominalReference,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    h = nominal.hist
    for ib, nominal_value in enumerate(h.values):
        nominal_err = h.errors[ib] if ib < len(h.errors) else np.nan
        rows.append(make_row(
            job_label=job.label,
            observable=job.result_kind,
            observable_pretty=job.observable_pretty or job.result_kind,
            hist_name=h.hist_name,
            cent=cent,
            method=job.method,
            iteration_display=job.iteration,
            bin_idx=ib + 1,
            bin_low=h.edges[ib],
            bin_high=h.edges[ib + 1],
            nominal_value=nominal_value,
            nominal_stat_abs=nominal_err,
            source_group="Nominal",
            source_name="nominal_value",
            variation_name="nominal",
            component_type="nominal",
            direction="none",
            variation_value=nominal_value,
            variation_stat_abs=nominal_err,
            diff=0.0,
            diff_pct=0.0,
            abs_pct=0.0,
            reference_file=nominal.file,
            variation_file=nominal.file,
            notes="central reference value",
        ))

        # Statistical uncertainty as its own component, to make later combination easy.
        rows.append(make_row(
            job_label=job.label,
            observable=job.result_kind,
            observable_pretty=job.observable_pretty or job.result_kind,
            hist_name=h.hist_name,
            cent=cent,
            method=job.method,
            iteration_display=job.iteration,
            bin_idx=ib + 1,
            bin_low=h.edges[ib],
            bin_high=h.edges[ib + 1],
            nominal_value=nominal_value,
            nominal_stat_abs=nominal_err,
            source_group="Nominal",
            source_name="statistical_uncertainty",
            variation_name="statistical",
            component_type="uncertainty",
            direction="symmetric",
            variation_value=np.nan,
            diff=nominal_err,
            diff_pct=pct(nominal_err, nominal_value),
            abs_pct=abs(pct(nominal_err, nominal_value)) if is_finite_number(pct(nominal_err, nominal_value)) else np.nan,
            reference_file=nominal.file,
            variation_file=nominal.file,
            notes="TH1 bin error from nominal histogram",
        ))
    return rows


def collect_iteration_rows(
    config: Dict[str, Any],
    project_dir: Path,
    job: Job,
    cent: int,
    nominal: NominalReference,
    extra_hist_map: Dict[str, str],
    strict: bool,
) -> List[Dict[str, Any]]:
    nominal_cfg = config.get("nominal", {})
    if not nominal_cfg.get("include_iteration_variations", True):
        return []

    rows: List[Dict[str, Any]] = []
    nominal_file = nominal.file
    hnom = nominal.hist
    for delta, name in [(-1, "iteration_minus_1"), (+1, "iteration_plus_1")]:
        it = int(job.iteration) + delta
        if it < 1:
            continue
        hist_name = build_hist_name(job.result_kind, cent, it, job.method, extra_hist_map)
        try:
            hvar = load_hist(nominal_file, hist_name)
        except Exception as e:
            msg = f"Cannot load {name} histogram {hist_name} from {nominal_file}: {e}"
            if strict:
                raise
            warn(msg)
            continue

        var_values = map_to_reference_edges(hvar.edges, hvar.values, hnom.edges)
        var_errors = map_to_reference_edges(hvar.edges, hvar.errors, hnom.edges)
        for ib, nominal_value in enumerate(hnom.values):
            var_value = var_values[ib] if ib < len(var_values) else np.nan
            var_err = var_errors[ib] if ib < len(var_errors) else np.nan
            diff = var_value - nominal_value if is_finite_number(var_value) and is_finite_number(nominal_value) else np.nan
            diff_pct = pct(diff, nominal_value)
            direction = "up" if is_finite_number(diff) and diff > 0 else ("down" if is_finite_number(diff) and diff < 0 else "zero")
            rows.append(make_row(
                job_label=job.label,
                observable=job.result_kind,
                observable_pretty=job.observable_pretty or job.result_kind,
                hist_name=hnom.hist_name,
                cent=cent,
                method=job.method,
                iteration_display=job.iteration,
                bin_idx=ib + 1,
                bin_low=hnom.edges[ib],
                bin_high=hnom.edges[ib + 1],
                nominal_value=nominal_value,
                nominal_stat_abs=hnom.errors[ib],
                source_group="Nominal",
                source_name="iteration",
                variation_name=name,
                variation_code=it,
                component_type="variation",
                direction=direction,
                variation_value=var_value,
                variation_stat_abs=var_err,
                diff=diff,
                diff_pct=diff_pct,
                abs_pct=abs(diff_pct) if is_finite_number(diff_pct) else np.nan,
                reference_file=nominal_file,
                variation_file=nominal_file,
                notes=f"display iteration {it} -> ROOT _it{it-1}_",
            ))
    return rows


def rows_from_variation_hist(
    *,
    job: Job,
    cent: int,
    nominal: NominalReference,
    hvar: HistData,
    source_group: str,
    source_name: str,
    variation_name: str,
    variation_code: Any,
    component_type: str,
    run_id: str = "",
    notes: str = "",
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    hnom = nominal.hist
    var_values = map_to_reference_edges(hvar.edges, hvar.values, hnom.edges)
    var_errors = map_to_reference_edges(hvar.edges, hvar.errors, hnom.edges)
    for ib, nominal_value in enumerate(hnom.values):
        var_value = var_values[ib] if ib < len(var_values) else np.nan
        var_err = var_errors[ib] if ib < len(var_errors) else np.nan
        diff = var_value - nominal_value if is_finite_number(var_value) and is_finite_number(nominal_value) else np.nan
        diff_pct = pct(diff, nominal_value)
        direction = "up" if is_finite_number(diff) and diff > 0 else ("down" if is_finite_number(diff) and diff < 0 else "zero")
        rows.append(make_row(
            job_label=job.label,
            observable=job.result_kind,
            observable_pretty=job.observable_pretty or job.result_kind,
            hist_name=hnom.hist_name,
            cent=cent,
            method=job.method,
            iteration_display=job.iteration,
            bin_idx=ib + 1,
            bin_low=hnom.edges[ib],
            bin_high=hnom.edges[ib + 1],
            nominal_value=nominal_value,
            nominal_stat_abs=hnom.errors[ib],
            source_group=source_group,
            source_name=source_name,
            variation_name=variation_name,
            variation_code=variation_code,
            component_type=component_type,
            direction=direction,
            variation_value=var_value,
            variation_stat_abs=var_err,
            diff=diff,
            diff_pct=diff_pct,
            abs_pct=abs(diff_pct) if is_finite_number(diff_pct) else np.nan,
            reference_file=nominal.file,
            variation_file=hvar.path,
            run_id=run_id,
            notes=notes,
        ))
    return rows


def collect_simple_group_rows(
    config: Dict[str, Any],
    project_dir: Path,
    group_cfg: Dict[str, Any],
    job: Job,
    cent: int,
    nominal: NominalReference,
    extra_hist_map: Dict[str, str],
    strict: bool,
) -> List[Dict[str, Any]]:
    if not group_cfg.get("enabled", True):
        return []

    source_group = str(group_cfg.get("source_group", group_cfg.get("name", "UnknownGroup")))
    rows: List[Dict[str, Any]] = []

    scan_dir = str(group_cfg.get("scan_dir", ""))
    context_base = {
        "project_dir": str(project_dir),
        "scan_dir": scan_dir,
        "job_label": job.label,
        "result_kind": job.result_kind,
        "cent": cent,
        "method": job.method,
        "iteration": job.iteration,
        "iteration_root": job.iteration - 1,
    }

    variations: List[Dict[str, Any]] = []

    # Option A: explicit variation list in JSON.
    for v in group_cfg.get("variations", []):
        vv = dict(v)
        vv.setdefault("source_group", source_group)
        variations.append(vv)

    # Option B: read variation list from summary.tsv.
    summary_path_raw = group_cfg.get("summary", "")
    if summary_path_raw:
        summary_path = resolve_path(expand_template(str(summary_path_raw), context_base), project_dir)
        try:
            sdf = load_summary_table(summary_path)
        except Exception:
            if strict:
                raise
            warn(f"Skipping {source_group}: cannot read summary {summary_path}")
            sdf = pd.DataFrame()

        if not sdf.empty:
            skip_names = set(str(x).strip().lower() for x in group_cfg.get("skip_source_names", ["nominal"]))
            skip_codes = set(normalize_code_text(x) for x in group_cfg.get("skip_variation_codes", []))

            source_name_cols = column_candidates(group_cfg, "source_name_column", "source_name_columns", "source_name")
            variation_name_cols = column_candidates(group_cfg, "variation_name_column", "variation_name_columns", "")
            if not variation_name_cols:
                variation_name_cols = list(source_name_cols)
            variation_code_cols = column_candidates(group_cfg, "variation_code_column", "variation_code_columns", "")

            # Optional mappings are useful when summary.tsv only contains run_id + numeric code.
            # Example: "code_label_map": {"1": "cheby2_bkg", "2": "double_gauss_sig"}.
            code_label_map = {normalize_code_text(k): str(v) for k, v in dict(group_cfg.get("code_label_map", {})).items()}
            source_name_map = {str(k): str(v) for k, v in dict(group_cfg.get("source_name_map", {})).items()}
            variation_name_map = {str(k): str(v) for k, v in dict(group_cfg.get("variation_name_map", {})).items()}

            prior_shape_auto_labels = bool(group_cfg.get("prior_shape_auto_labels", False))
            prior_shape_scheme = str(group_cfg.get("prior_shape_code_scheme", "zero_based"))
            prior_shape_filter_relevant = bool(group_cfg.get("prior_shape_filter_relevant_to_job", True))
            prior_shape_sign_map = {str(k): str(v) for k, v in dict(group_cfg.get("prior_shape_sign_map", {})).items()} or None
            d0meson_auto_labels = bool(group_cfg.get("d0meson_auto_labels", False))

            enabled_filter = group_cfg.get("summary_filter", {})
            if enabled_filter:
                for col, raw_val in enabled_filter.items():
                    if col not in sdf.columns:
                        raise KeyError(f"{source_group}: summary_filter column {col!r} not found in {summary_path}")
                    sdf = sdf[sdf[col].astype(str) == str(raw_val)]

            for _, row in sdf.iterrows():
                row_dict = row.to_dict()
                run_id = normalize_run_id(row_dict.get("run_id", ""))
                vcode = first_row_value(row_dict, variation_code_cols, "")
                vcode_txt = normalize_code_text(vcode)
                if vcode_txt in skip_codes:
                    continue

                prior_notes = ""
                prior_auto = None
                if prior_shape_auto_labels:
                    prior_auto = prior_shape_label_from_code(
                        vcode_txt,
                        job.result_kind,
                        scheme=prior_shape_scheme,
                        sign_map=prior_shape_sign_map,
                        filter_relevant=prior_shape_filter_relevant,
                    )
                    if prior_auto is None:
                        continue

                source_name_raw = str(first_row_value(row_dict, source_name_cols, ""))
                variation_name_raw = str(first_row_value(row_dict, variation_name_cols, source_name_raw))

                d0_notes = ""
                d0_auto = None
                if d0meson_auto_labels:
                    d0_auto = d0meson_label_from_paper_label(variation_name_raw or source_name_raw)

                source_name = source_name_map.get(source_name_raw, source_name_raw).strip()
                variation_name = variation_name_map.get(variation_name_raw, variation_name_raw).strip()

                if prior_auto is not None:
                    source_name, variation_name, prior_notes = prior_auto
                elif d0_auto is not None:
                    source_name, variation_name, d0_notes = d0_auto
                else:
                    if (not source_name or source_name == run_id) and vcode_txt in code_label_map:
                        source_name = code_label_map[vcode_txt]
                    if (not variation_name or variation_name == run_id) and vcode_txt in code_label_map:
                        variation_name = code_label_map[vcode_txt]

                source_name = source_name or variation_name or run_id
                variation_name = variation_name or source_name or run_id

                if source_name.strip().lower() in skip_names or variation_name.strip().lower() in skip_names:
                    continue

                column_notes = "; ".join([
                    available_column_note("source_col", source_name_cols, row_dict),
                    available_column_note("variation_col", variation_name_cols, row_dict),
                    available_column_note("code_col", variation_code_cols, row_dict) if variation_code_cols else "code_col=disabled",
                ])
                notes_joined = "; ".join([x for x in [prior_notes, d0_notes, column_notes] if x])

                variations.append({
                    "run_id": run_id,
                    "source_name": source_name,
                    "variation_name": variation_name,
                    "variation_code": vcode_txt,
                    "summary_row": row_dict,
                    "notes": notes_joined,
                })

    if not variations:
        return []

    default_file_pattern = str(group_cfg.get("file_pattern", ""))
    default_file_glob = str(group_cfg.get("file_glob", ""))
    default_hist_name = group_cfg.get("hist_name", "")
    default_iteration = int(group_cfg.get("iteration", job.iteration))
    default_method = str(group_cfg.get("method", job.method))

    for variation in variations:
        vctx = dict(context_base)
        vctx.update({k: v for k, v in variation.items() if not isinstance(v, dict)})
        if "summary_row" in variation:
            for k, v in variation["summary_row"].items():
                vctx[str(k)] = v
        run_id = normalize_run_id(vctx.get("run_id", ""))
        vctx["run_id"] = run_id

        source_name = str(variation.get("source_name", variation.get("name", run_id or "variation")))
        variation_name = str(variation.get("variation_name", source_name))
        variation_code = variation.get("variation_code", variation.get("code", ""))
        component_type = str(variation.get("component_type", group_cfg.get("component_type", "variation")))

        # Resolve variation ROOT file.
        var_file: Optional[Path] = None
        if variation.get("file"):
            var_file = find_root_file_from_pattern(str(variation["file"]), project_dir, vctx, required=False)
        if var_file is None and variation.get("file_glob"):
            var_file = find_root_file_from_pattern(str(variation["file_glob"]), project_dir, vctx, required=False)
        if var_file is None and default_file_pattern:
            var_file = find_root_file_from_pattern(default_file_pattern, project_dir, vctx, required=False)
        if var_file is None and default_file_glob:
            var_file = find_root_file_from_pattern(default_file_glob, project_dir, vctx, required=False)

        # Fallbacks for common scan layout.
        if var_file is None and scan_dir and run_id:
            fallback_patterns = [
                f"{scan_dir}/Output/OutputSpectra{run_id}.root",
                f"{scan_dir}/Output/*{run_id}*.root",
                f"{scan_dir}/runs/{run_id}/OutputSpectra*.root",
                f"{scan_dir}/runs/{run_id}/*.root",
            ]
            for ptn in fallback_patterns:
                var_file = find_root_file_from_pattern(ptn, project_dir, vctx, required=False)
                if var_file is not None:
                    break

        if var_file is None:
            msg = f"{source_group}/{source_name}: could not resolve variation ROOT file; run_id={run_id}"
            if strict:
                raise FileNotFoundError(msg)
            warn(msg)
            continue

        v_iteration = int(variation.get("iteration", default_iteration))
        v_method = str(variation.get("method", default_method))
        if default_hist_name:
            hist_name = expand_template(str(default_hist_name), {**vctx, "iteration": v_iteration, "iteration_root": v_iteration - 1})
        else:
            hist_name = build_hist_name(job.result_kind, cent, v_iteration, v_method, extra_hist_map)

        try:
            hvar = load_hist(var_file, hist_name)
        except Exception as e:
            msg = f"{source_group}/{source_name}: cannot read {hist_name} from {var_file}: {e}"
            if strict:
                raise
            warn(msg)
            continue

        rows.extend(rows_from_variation_hist(
            job=job,
            cent=cent,
            nominal=nominal,
            hvar=hvar,
            source_group=source_group,
            source_name=source_name,
            variation_name=variation_name,
            variation_code=variation_code,
            component_type=component_type,
            run_id=run_id,
            notes="; ".join([x for x in [str(group_cfg.get("notes", "")), str(variation.get("notes", ""))] if x]),
        ))

    return rows


def collect_binning_rows(
    config: Dict[str, Any],
    project_dir: Path,
    job: Job,
    cent: int,
    nominal: NominalReference,
    extra_hist_map: Dict[str, str],
    strict: bool,
) -> List[Dict[str, Any]]:
    """
    Binning-systematics reader.

    It selects run_ids from stability.tsv using exact filters plus metric cuts.
    Metric cuts can be defined globally, per source, and overridden per
    job_label/result_kind/centrality/source_name in the JSON config.

    The default output is one RMS row per bin and source, where RMS is computed
    from percent differences of selected binning variations relative to the
    independent nominal reference histogram:

        diff_pct(run, bin) = 100 * (variation / nominal - 1)
        rms(bin)           = sqrt(mean(diff_pct^2))

    Envelope/asymmetric rows can be re-enabled through aggregate_modes, but the
    recommended first-pass table uses only RMS.
    """
    bcfg = config.get("binning", {})
    if not bcfg.get("enabled", False):
        return []

    rows: List[Dict[str, Any]] = []
    hnom = nominal.hist
    global_aggregate_modes = [str(x).strip() for x in as_list(bcfg.get("aggregate_modes", ["rms"])) if str(x).strip()]

    for source in bcfg.get("sources", []):
        if not source.get("enabled", True):
            continue
        source_group = str(source.get("source_group", "Binning"))
        source_name = str(source.get("source_name", source.get("name", "binning_source")))
        context0 = {
            "project_dir": str(project_dir),
            "job_label": job.label,
            "base_job_label": binning_base_job_label(job.label),
            "result_kind": job.result_kind,
            "cent": cent,
            "centrality": cent,
            "source_name": source_name,
            "method": job.method,
            "iteration": job.iteration,
            "iteration_root": job.iteration - 1,
        }
        scan_dir = resolve_binning_scan_dir(source, job, cent, context0)
        context = dict(context0)
        context["scan_dir"] = scan_dir

        summary_path = resolve_path(expand_template(str(source.get("summary", "")), context), project_dir)
        stability_path = resolve_path(expand_template(str(source.get("stability", "")), context), project_dir)
        if not summary_path.exists() or not stability_path.exists():
            msg = f"{source_group}/{source_name}: missing summary/stability ({summary_path}, {stability_path})"
            if strict:
                raise FileNotFoundError(msg)
            warn(msg)
            continue

        try:
            stability = load_stability_table(stability_path)
            summary = load_summary_table(summary_path)
            if "run_id" in summary.columns:
                df = stability.merge(summary, on="run_id", how="left", validate="many_to_one")
            else:
                df = stability

            selection = dict(source.get("selection", {}))
            selection.setdefault("cent", "{cent}")
            # Do exact filters first.  Metric cuts are resolved separately below
            # so per-job/per-centrality overrides can inherit defaults cleanly.
            exact_selection = dict(selection)
            for k in ["metric_min_col", "metric_min_value", "metric_max_col", "metric_max_value"]:
                exact_selection.pop(k, None)

            auto_notes: List[str] = []

            obs_value = str(exact_selection.get("observable", "") or "").strip()
            use_auto_obs = obs_value in ("{auto_observable}", "auto", "AUTO", "*")
            if use_auto_obs:
                exact_selection.pop("observable", None)

            axis_value = str(exact_selection.get("axis", "") or "").strip()
            use_auto_axis = axis_value in ("{auto_axis}", "auto", "AUTO", "*")
            if use_auto_axis:
                exact_selection.pop("axis", None)

            selected = filter_stability_for_source(df, exact_selection, context)

            if use_auto_obs:
                selected, note = apply_auto_observable_filter(selected, source, job, context)
                auto_notes.append(f"observable={note}")
            if use_auto_axis:
                selected, note = apply_auto_axis_filter(selected, source, context)
                auto_notes.append(f"axis={note}")

            n_before_true_range = len(selected)
            true_range_filter = resolve_binning_true_range_filter(bcfg, source, context)
            selected = apply_binning_true_range_filter(selected, true_range_filter, source_name)

            n_before_reco_range = len(selected)
            reco_range_filter = resolve_binning_reco_range_filter(bcfg, source, context)
            selected = apply_binning_reco_range_filter(selected, reco_range_filter, source_name)

            n_before_metric = len(selected)
            metric_cuts = resolve_binning_metric_cuts(bcfg, source, context)
            selected = apply_binning_metric_cuts(selected, metric_cuts)
        except Exception as e:
            if strict:
                raise
            warn(f"{source_group}/{source_name}: selection failed: {e}")
            continue

        run_ids = sorted(selected["run_id"].dropna().astype(str).unique().tolist(), key=run_id_to_int)
        if not run_ids:
            note = "; ".join(auto_notes) if 'auto_notes' in locals() and auto_notes else ""
            warn(
                f"{source_group}/{source_name}: selected zero binning variations for {job.label}, cent={cent}; "
                f"before_true_range={locals().get('n_before_true_range', 'n/a')}; "
                f"before_reco_range={locals().get('n_before_reco_range', 'n/a')}; "
                f"before_metric={locals().get('n_before_metric', 'n/a')}; "
                f"true_range: {true_range_filter_to_note(locals().get('true_range_filter', {}))}; "
                f"reco_range: {reco_range_filter_to_note(locals().get('reco_range_filter', {}))}; "
                f"cuts: {metric_cuts_to_note(metric_cuts)}"
                + (f"; {note}" if note else "")
            )
            continue

        file_pattern = str(source.get("file_pattern", ""))
        file_glob = str(source.get("file_glob", ""))
        hist_name_template = str(source.get("hist_name", ""))
        diffs_by_bin: List[List[float]] = [[] for _ in range(len(hnom.values))]
        loaded = 0
        skipped: List[str] = []

        for run_id in run_ids:
            vctx = dict(context)
            vctx["run_id"] = run_id
            var_file: Optional[Path] = None
            if file_pattern:
                var_file = find_root_file_from_pattern(file_pattern, project_dir, vctx, required=False)
            if var_file is None and file_glob:
                var_file = find_root_file_from_pattern(file_glob, project_dir, vctx, required=False)
            if var_file is None and scan_dir:
                for ptn in [
                    f"{scan_dir}/Output/OutputSpectra{run_id}.root",
                    f"{scan_dir}/Output/*{run_id}*.root",
                    f"{scan_dir}/runs/{run_id}/OutputSpectra*.root",
                    f"{scan_dir}/runs/{run_id}/*.root",
                ]:
                    var_file = find_root_file_from_pattern(ptn, project_dir, vctx, required=False)
                    if var_file is not None:
                        break
            if var_file is None:
                skipped.append(f"{run_id}: missing ROOT file")
                continue

            if hist_name_template:
                hist_name = expand_template(hist_name_template, vctx)
            else:
                hist_name = build_hist_name(job.result_kind, cent, job.iteration, job.method, extra_hist_map)
            try:
                hvar = load_hist(var_file, hist_name)
                var_values = map_to_reference_edges(hvar.edges, hvar.values, hnom.edges)
            except Exception as e:
                skipped.append(f"{run_id}: {e}")
                continue

            loaded += 1
            for ib, nominal_value in enumerate(hnom.values):
                var_value = var_values[ib] if ib < len(var_values) else np.nan
                d_pct = pct(var_value - nominal_value, nominal_value)
                if is_finite_number(d_pct):
                    diffs_by_bin[ib].append(float(d_pct))

        if loaded == 0:
            msg = f"{source_group}/{source_name}: no selected ROOT files could be loaded. First skipped: {skipped[:3]}"
            if strict:
                raise RuntimeError(msg)
            warn(msg)
            continue

        aggregate_modes = [str(x).strip() for x in as_list(source.get("aggregate_modes", global_aggregate_modes)) if str(x).strip()]
        if not aggregate_modes:
            aggregate_modes = ["rms"]
        allowed_modes = {"rms", "envelope", "envelope_down", "envelope_up"}
        unknown_modes = [m for m in aggregate_modes if m not in allowed_modes]
        if unknown_modes:
            msg = f"{source_group}/{source_name}: unknown aggregate_modes {unknown_modes}; allowed={sorted(allowed_modes)}"
            if strict:
                raise ValueError(msg)
            warn(msg)
            aggregate_modes = [m for m in aggregate_modes if m in allowed_modes] or ["rms"]

        for ib, diffs in enumerate(diffs_by_bin):
            nominal_value = hnom.values[ib]
            nominal_err = hnom.errors[ib]
            arr = np.asarray(diffs, dtype=float)
            if arr.size:
                env_down = max(0.0, -float(np.min(arr)))
                env_up = max(0.0, float(np.max(arr)))
                rms = float(np.sqrt(np.mean(arr * arr)))
            else:
                env_down = env_up = rms = np.nan

            common = dict(
                job_label=job.label,
                observable=job.result_kind,
                observable_pretty=job.observable_pretty or job.result_kind,
                hist_name=hnom.hist_name,
                cent=cent,
                method=job.method,
                iteration_display=job.iteration,
                bin_idx=ib + 1,
                bin_low=hnom.edges[ib],
                bin_high=hnom.edges[ib + 1],
                nominal_value=nominal_value,
                nominal_stat_abs=nominal_err,
                source_group=source_group,
                source_name=source_name,
                variation_code="",
                component_type="aggregate",
                variation_value=np.nan,
                variation_stat_abs=np.nan,
                n_variations=arr.size,
                reference_file=nominal.file,
                variation_file=scan_dir,
                run_id="",
                notes=(
                    f"RMS/envelope from binning variations relative to independent nominal; "
                    f"selected={len(run_ids)}, loaded={loaded}, skipped={len(skipped)}; "
                    f"auto_filters={'; '.join(auto_notes) if auto_notes else 'not used'}; "
                    f"true_range: {true_range_filter_to_note(true_range_filter)}; "
                    f"reco_range: {reco_range_filter_to_note(reco_range_filter)}; "
                    f"cuts: {metric_cuts_to_note(metric_cuts)}"
                ),
            )
            if "rms" in aggregate_modes:
                rows.append(make_row(**common, variation_name="rms", direction="symmetric", diff=np.nan, diff_pct=rms, abs_pct=rms))
            if "envelope" in aggregate_modes or "envelope_down" in aggregate_modes:
                rows.append(make_row(**common, variation_name="envelope_down", direction="down", diff=np.nan, diff_pct=-env_down, abs_pct=env_down))
            if "envelope" in aggregate_modes or "envelope_up" in aggregate_modes:
                rows.append(make_row(**common, variation_name="envelope_up", direction="up", diff=np.nan, diff_pct=env_up, abs_pct=env_up))

    return rows


# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------

def default_template_config() -> Dict[str, Any]:
    return {
        "project_dir": ".",
        "output_tsv": "Systematics/systematics_components_v14.tsv",
        "output_csv": "Systematics/systematics_components_v14.csv",
        "strict": False,
        "method": "ICS",
        "iteration": 4,
        "extra_hist_map": {},
        "jobs": [
            {"label": "z", "result_kind": "z", "observable_pretty": "z", "centralities": [0, 1, 2]},
            {"label": "l11", "result_kind": "#lambda^{1}_{1}", "observable_pretty": "#lambda^{1}_{1}", "centralities": [0, 1, 2]},
            {"label": "l15", "result_kind": "#lambda^{1}_{1.5}", "observable_pretty": "#lambda^{1}_{1.5}", "centralities": [0, 1, 2]},
            {"label": "l21", "result_kind": "#lambda^{1}_{2}", "observable_pretty": "#lambda^{1}_{2}", "centralities": [0, 1, 2]},
            {"label": "l31", "result_kind": "#lambda^{1}_{3}", "observable_pretty": "#lambda^{1}_{3}", "centralities": [0, 1, 2]},
            {"label": "l051", "result_kind": "#lambda^{1}_{0.5}", "observable_pretty": "#lambda^{1}_{0.5}", "centralities": [0, 1, 2]},
            {"label": "pTD", "result_kind": "P_{T}^{D}", "observable_pretty": "P_{T}^{D}", "centralities": [0, 1, 2]},
            {"label": "RCP_z_5_20", "result_kind": "RCP 5-20 z", "observable_pretty": "R_{CP}(z), 5<p_{T,Jet}<20", "centralities": [0]},
            {"label": "RCP_l11_5_20", "result_kind": "RCP 5-20 #lambda^{1}_{1}", "observable_pretty": "R_{CP}(#lambda^{1}_{1}), 5<p_{T,Jet}<20", "centralities": [0]},
            {"label": "RCP_l15_5_20", "result_kind": "RCP 5-20 #lambda^{1}_{1.5}", "observable_pretty": "R_{CP}(#lambda^{1}_{1.5}), 5<p_{T,Jet}<20", "centralities": [0]},
            {"label": "RCP_l21_5_20", "result_kind": "RCP 5-20 #lambda^{1}_{2}", "observable_pretty": "R_{CP}(#lambda^{1}_{2}), 5<p_{T,Jet}<20", "centralities": [0]},
            {"label": "RCP_l31_5_20", "result_kind": "RCP 5-20 #lambda^{1}_{3}", "observable_pretty": "R_{CP}(#lambda^{1}_{3}), 5<p_{T,Jet}<20", "centralities": [0]},
            {"label": "RCP_l051_5_20", "result_kind": "RCP 5-20 #lambda^{1}_{0.5}", "observable_pretty": "R_{CP}(#lambda^{1}_{0.5}), 5<p_{T,Jet}<20", "centralities": [0]},
            {"label": "RCP_pTD_5_20", "result_kind": "RCP 5-20 P_{T}^{D}", "observable_pretty": "R_{CP}(P_{T}^{D}), 5<p_{T,Jet}<20", "centralities": [0]},
        ],
        "nominal": {
            "enabled": True,
            "file_glob": "Systematics/Nominal/scanNominal/Output/OutputSpectra*.root",
            "include_iteration_variations": True,
        },
        "groups": [
            {
                "enabled": True,
                "source_group": "JetsReco",
                "scan_dir": "Systematics/JetsReco/scanJetsReco",
                "summary": "{scan_dir}/summary.tsv",
                "file_pattern": "{scan_dir}/Output/OutputSpectra{run_id}.root",
                "source_name_column": "jets_reco_label",
                "variation_name_column": "jets_reco_label",
                "variation_code_column": "systematicSPlot",
                "skip_source_names": ["nominal"],
                "component_type": "variation",
                "notes": "auto-loaded from JetsReco summary.tsv",
            },
            {
                "enabled": True,
                "source_group": "sWeight",
                "scan_dir": "Systematics/sWeight/scanSWeight",
                "summary": "{scan_dir}/summary.tsv",
                "file_pattern": "{scan_dir}/Output/OutputSpectra{run_id}.root",
                "source_name_column": "sweight_label",
                "variation_name_column": "sweight_label",
                "variation_code_column": "systematicSPlot",
                "skip_source_names": ["nominal"],
                "skip_variation_codes": ["0"],
                "code_label_map": {
                    "1": "cheby2_bkg",
                    "2": "double_gauss_sig",
                    "3": "student_t_sig",
                    "4": "narrow_fit_range",
                    "5": "wide_fit_range",
                    "6": "zeroed_negative_bins"
                },
                "component_type": "variation",
                "notes": "auto-labels numeric systematicSPlot codes when sweight_label is absent",
            },
            {
                "enabled": True,
                "source_group": "PriorShape",
                "scan_dir": "Systematics/PriorShape/scanPriorShape",
                "summary": "{scan_dir}/summary.tsv",
                "file_pattern": "{scan_dir}/Output/OutputSpectra{run_id}.root",
                "source_name_column": "prior_label",
                "variation_name_column": "prior_label",
                "variation_code_column": "usePriorShapeWeighting",
                "skip_source_names": ["nominal"],
                "skip_variation_codes": ["0"],
                "prior_shape_auto_labels": True,
                "prior_shape_code_scheme": "first_second",
                "prior_shape_filter_relevant_to_job": True,
                "prior_shape_sign_map": {
                    "0": "plus20",
                    "1": "minus20",
                    "2": "minus20"
                },
                "component_type": "variation",
                "notes": "scanPriorShape scheme: 10/11=jet pT +/-; 20/21=second variable +/-",
            },
            {
                "enabled": True,
                "source_group": "D0Meson",
                "scan_dir": "Systematics/D0Meson/scanPaperSys",
                "summary": "{scan_dir}/summary.tsv",
                "file_pattern": "{scan_dir}/Output/OutputSpectra{run_id}.root",
                "source_name_columns": [
                    "paper_label",
                    "d0meson_label",
                    "d0_meson_label",
                    "paperSys_label",
                    "paper_sys_label",
                    "d0_label",
                    "systematic_label",
                    "label"
                ],
                "variation_name_columns": [
                    "paper_label",
                    "d0meson_label",
                    "d0_meson_label",
                    "paperSys_label",
                    "paper_sys_label",
                    "d0_label",
                    "systematic_label",
                    "label"
                ],
                "variation_code_columns": [
                    "systematicSPlot",
                    "paperSys",
                    "paper_sys",
                    "systematicD0Meson",
                    "systematicD0",
                    "systematicPaper",
                    "systematicSPlot",
                    "systematic"
                ],
                "skip_source_names": ["nominal"],
                "skip_variation_codes": ["0"],
                "d0meson_auto_labels": True,
                "component_type": "variation",
                "notes": "scanPaperSys: paper_label is parsed into paired D0-meson systematic sources",
            },
        ],
        "binning": {
            "enabled": True,
            "aggregate_modes": ["rms"],
            "default_metric_cuts": {},
            "default_true_range_filter": {
                "enabled": False,
                "mode": "ignore",
                "true_min": "",
                "true_max": ""
            },
            "true_range_overrides": [
                {
                    "match": {
                        "job_label": "l11",
                        "centrality": 0,
                        "source_name": "binning_first_variable"
                    },
                    "true_range_filter": {
                        "enabled": False,
                        "mode": "exact",
                        "true_min": "",
                        "true_max": ""
                    },
                    "notes": "example true-level endpoint filter. enabled=false or mode=ignore means use all true ranges."
                }
            ],
            "metric_cut_overrides": [
                {
                    "match": {"job_label": "l11", "centrality": 0, "source_name": "binning_first_variable"},
                    "metric_cuts": {
                        "weighted_drift_pct": {"min": "", "max": ""},
                        "rms_drift_pct": {"min": "", "max": ""},
                        "worst_bin_pct": {"min": "", "max": ""}
                    },
                    "notes": "example rule; empty min/max means disabled. Add concrete values per observable/centrality/source."
                }
            ],
            "sources": [
                {
                    "enabled": True,
                    "source_group": "Binning",
                    "source_name": "binning_first_variable",
                    "scan_dir": "Systematics/Binning/scanJetPt",
                    "summary": "{scan_dir}/summary.tsv",
                    "stability": "{scan_dir}/stability.tsv",
                    "file_pattern": "{scan_dir}/Output/OutputSpectra{run_id}.root",
                    "aggregate_modes": ["rms"],
                    "selection": {
                        "cent": "{cent}",
                        "dim": "2D",
                        "observable": "{auto_observable}",
                        "axis": "{auto_axis}",
                        "exact": {}
                    },
                    "metric_cuts": {},
                    "observable_map": {
                        "default": "{auto_observable}"
                    },
                    "axis_candidates": ["PT", "pt", "pT", "JetPt", "jetPt", "jet_pt", "X", "x", "main", "first"],
                    "on_no_axis_match": "no_filter",
                    "on_no_observable_match": "no_filter"
                },
                {
                    "enabled": True,
                    "source_group": "Binning",
                    "source_name": "binning_second_variable",
                    "scan_dir_map": {
                        "z": "Systematics/Binning/scanZ",
                        "l11": "Systematics/Binning/scanL11",
                        "l15": "Systematics/Binning/scanL11p5",
                        "l21": "Systematics/Binning/scanL12",
                        "l31": "Systematics/Binning/scanL13",
                        "l051": "Systematics/Binning/scanL10p5",
                        "pTD": "Systematics/Binning/scanPTD",
                        "default": "Systematics/Binning/scanZ"
                    },
                    "summary": "{scan_dir}/summary.tsv",
                    "stability": "{scan_dir}/stability.tsv",
                    "file_pattern": "{scan_dir}/Output/OutputSpectra{run_id}.root",
                    "aggregate_modes": ["rms"],
                    "selection": {
                        "cent": "{cent}",
                        "dim": "2D",
                        "observable": "{auto_observable}",
                        "axis": "{auto_axis}",
                        "exact": {}
                    },
                    "metric_cuts": {},
                    "observable_map": {
                        "default": "{auto_observable}"
                    },
                    "axis_candidates": ["var", "Y", "y", "secondary", "second", "observable"],
                    "on_no_axis_match": "no_filter",
                    "on_no_observable_match": "no_filter"
                }
            ]
        }
    }


def build_components(config: Dict[str, Any]) -> pd.DataFrame:
    project_dir = resolve_path(str(config.get("project_dir", ".")), Path.cwd()).resolve()
    strict = bool(config.get("strict", False))
    global_method = str(config.get("method", "ICS"))
    global_iteration = int(config.get("iteration", 4))
    extra_hist_map = {str(k): str(v) for k, v in dict(config.get("extra_hist_map", {})).items()}

    jobs = [Job.from_dict(j, global_method, global_iteration) for j in config.get("jobs", [])]
    jobs = [j for j in jobs if j.enabled]
    if not jobs:
        raise ValueError("No enabled jobs in config['jobs'].")

    all_rows: List[Dict[str, Any]] = []
    errors: List[str] = []

    total_tasks = sum(len(j.centralities) for j in jobs)
    task_idx = 0
    for job in jobs:
        for cent in job.centralities:
            task_idx += 1
            info(f"[{task_idx}/{total_tasks}] {job.label}, cent={cent}, result={job.result_kind}")
            try:
                nominal = load_nominal_reference(config, project_dir, job, cent, extra_hist_map)
            except Exception as e:
                msg = f"Nominal load failed for {job.label}, cent={cent}: {e}"
                if strict:
                    raise
                warn(msg)
                errors.append(msg)
                continue

            try:
                all_rows.extend(collect_nominal_rows(config, project_dir, job, cent, nominal))
                all_rows.extend(collect_iteration_rows(config, project_dir, job, cent, nominal, extra_hist_map, strict))
            except Exception as e:
                msg = f"Nominal/iteration rows failed for {job.label}, cent={cent}: {e}"
                if strict:
                    raise
                warn(msg)
                errors.append(msg)

            for group_cfg in config.get("groups", []):
                try:
                    all_rows.extend(collect_simple_group_rows(
                        config=config,
                        project_dir=project_dir,
                        group_cfg=group_cfg,
                        job=job,
                        cent=cent,
                        nominal=nominal,
                        extra_hist_map=extra_hist_map,
                        strict=strict,
                    ))
                except Exception as e:
                    msg = f"Group {group_cfg.get('source_group', group_cfg.get('name', 'unknown'))} failed for {job.label}, cent={cent}: {e}"
                    if strict:
                        raise
                    warn(msg)
                    errors.append(msg)

            try:
                all_rows.extend(collect_binning_rows(config, project_dir, job, cent, nominal, extra_hist_map, strict))
            except Exception as e:
                msg = f"Binning failed for {job.label}, cent={cent}: {e}"
                if strict:
                    raise
                warn(msg)
                errors.append(msg)

    df = pd.DataFrame(all_rows, columns=DEFAULT_OUTPUT_COLUMNS)
    if errors:
        df.attrs["errors"] = errors
    return df


def write_outputs(df: pd.DataFrame, config: Dict[str, Any]) -> None:
    project_dir = resolve_path(str(config.get("project_dir", ".")), Path.cwd()).resolve()
    output_tsv = config.get("output_tsv", "")
    output_csv = config.get("output_csv", "")

    if output_tsv:
        out = resolve_path(str(output_tsv), project_dir)
        out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, sep="\t", index=False, float_format="%.10g")
        info(f"Wrote TSV: {out}")

    if output_csv:
        out = resolve_path(str(output_csv), project_dir)
        out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, index=False, float_format="%.10g")
        info(f"Wrote CSV: {out}")

    errors = df.attrs.get("errors", [])
    if errors:
        err_path = None
        if output_tsv:
            err_path = resolve_path(str(output_tsv), project_dir).with_suffix(".errors.txt")
        elif output_csv:
            err_path = resolve_path(str(output_csv), project_dir).with_suffix(".errors.txt")
        if err_path:
            with err_path.open("w", encoding="utf-8") as f:
                for e in errors:
                    f.write(str(e) + "\n")
            warn(f"Some inputs were skipped. Details: {err_path}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build a long table of individual systematic components.")
    p.add_argument("--config", default="", help="JSON config file.")
    p.add_argument("--write-template", default="", help="Write a template JSON config and exit.")
    p.add_argument("--output", default="", help="Override output_tsv from config.")
    p.add_argument("--strict", action="store_true", help="Fail immediately on missing files/histograms.")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if args.write_template:
        out = Path(args.write_template)
        write_json(out, default_template_config())
        info(f"Template config written to: {out}")
        return 0

    if not args.config:
        print("[error] Use --config CONFIG.json or --write-template CONFIG.json", file=sys.stderr)
        return 2

    config_path = Path(args.config)
    config = read_json(config_path)
    if args.output:
        config["output_tsv"] = args.output
    if args.strict:
        config["strict"] = True

    try:
        df = build_components(config)
        write_outputs(df, config)
        info(f"Rows: {len(df)}")
        return 0
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
