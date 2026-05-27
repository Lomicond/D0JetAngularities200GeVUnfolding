#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
systematics_final_combiner_v9.py

Phase-3 combiner + final plotting for D0-jet generalized-angularity systematics.

Input:
  Long component table, e.g. Systematics/systematics_components_v15.tsv

Output:
  One row per histogram/bin with final grouped systematic uncertainties.
  Optional final plots: nominal value with statistical error bars and systematic uncertainty boxes.

Main ideas:
  * each group is selected from the long table via JSON rules,
  * variations inside each group are combined by a configurable method,
  * optional Barlow test can be applied before group combination,
  * final groups are combined into total systematic uncertainty.

Typical use:
  python3 Systematics/systematics_final_combiner_v9.py \
    --write-template Systematics/systematics_final_combiner_config_v1.json

  # edit JSON if needed, then
  python3 Systematics/systematics_final_combiner_v9.py \
    --config Systematics/systematics_final_combiner_config_v1.json

Notes:
  * Percent uncertainties are relative to nominal_value.
  * Up/down are stored as positive magnitudes.
  * Symmetric value is max(up, down), unless the method itself is symmetric.
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


# -----------------------------------------------------------------------------
# Basic helpers
# -----------------------------------------------------------------------------

def info(msg: str) -> None:
    print(f"[info] {msg}")


def warn(msg: str) -> None:
    print(f"[warning] {msg}", file=sys.stderr)


def read_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: Path, data: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")


def resolve_path(path: str, base: Path) -> Path:
    p = Path(os.path.expandvars(os.path.expanduser(str(path))))
    if not p.is_absolute():
        p = base / p
    return p


def to_float(x: Any, default: float = math.nan) -> float:
    try:
        if x is None:
            return default
        s = str(x).strip()
        if not s or s.lower() in ("nan", "none", "<na>"):
            return default
        return float(s)
    except Exception:
        return default


def isfinite(x: Any) -> bool:
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def fmt(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, float):
        if not math.isfinite(x):
            return ""
        return f"{x:.10g}"
    try:
        f = float(x)
        if math.isfinite(f):
            return f"{f:.10g}"
    except Exception:
        pass
    return str(x)


def pct_to_abs(pct_value: float, nominal_value: float) -> float:
    if not isfinite(pct_value) or not isfinite(nominal_value):
        return math.nan
    return abs(float(nominal_value)) * float(pct_value) / 100.0


def safe_name(text: Any) -> str:
    s = str(text or "").strip().lower()
    s = s.replace("%", "pct")
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "unnamed"


# -----------------------------------------------------------------------------
# Input format
# -----------------------------------------------------------------------------

KEY_COLS = [
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
]

BASE_NUM_COLS = ["nominal_value", "nominal_stat_abs", "nominal_stat_pct"]


def make_key(row: Dict[str, str]) -> Tuple[str, ...]:
    return tuple(str(row.get(c, "")) for c in KEY_COLS)


def read_components(path: Path) -> Tuple[List[str], Dict[Tuple[str, ...], List[Dict[str, str]]]]:
    bins: Dict[Tuple[str, ...], List[Dict[str, str]]] = defaultdict(list)
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Input table {path} has no header.")
        fieldnames = list(reader.fieldnames)
        for row in reader:
            bins[make_key(row)].append(row)
    return fieldnames, bins


def base_info_from_rows(key: Tuple[str, ...], rows: List[Dict[str, str]]) -> Dict[str, Any]:
    out = {c: key[i] for i, c in enumerate(KEY_COLS)}

    nominal_row = None
    stat_row = None
    for r in rows:
        if r.get("source_group") == "Nominal" and r.get("source_name") == "nominal_value":
            nominal_row = r
            break
    for r in rows:
        if r.get("source_group") == "Nominal" and r.get("source_name") == "statistical_uncertainty":
            stat_row = r
            break

    source = nominal_row or rows[0]
    for c in BASE_NUM_COLS:
        out[c] = source.get(c, "")

    if stat_row is not None:
        out["stat_abs"] = stat_row.get("diff", stat_row.get("nominal_stat_abs", ""))
        out["stat_pct"] = stat_row.get("diff_pct", stat_row.get("nominal_stat_pct", ""))
    else:
        out["stat_abs"] = source.get("nominal_stat_abs", "")
        out["stat_pct"] = source.get("nominal_stat_pct", "")

    return out


# -----------------------------------------------------------------------------
# Selection rules
# -----------------------------------------------------------------------------

def _parse_list_like_expected(expected: Any) -> Any:
    """Accept real JSON lists and also GUI-produced stringified lists.

    The GUI editor stores free-text fields as strings.  If a user writes
    ["a", "b"] or ['a', 'b'] in source_name, older combiner versions
    treated the whole text as one literal source name.  This helper converts
    such strings back to a Python list before matching.
    """
    if not isinstance(expected, str):
        return expected
    s = expected.strip()
    if not (s.startswith("[") and s.endswith("]")):
        return expected
    try:
        parsed = ast.literal_eval(s)
    except Exception:
        return expected
    if isinstance(parsed, (list, tuple, set)):
        return [str(x) for x in parsed]
    return expected


def value_matches(actual: str, expected: Any) -> bool:
    actual_s = str(actual or "")
    expected = _parse_list_like_expected(expected)
    if isinstance(expected, list):
        return actual_s in [str(x) for x in expected]
    if isinstance(expected, dict):
        if "regex" in expected:
            return re.search(str(expected["regex"]), actual_s) is not None
        if "not_regex" in expected:
            return re.search(str(expected["not_regex"]), actual_s) is None
        if "in" in expected:
            return actual_s in [str(x) for x in expected["in"]]
        if "not_in" in expected:
            return actual_s not in [str(x) for x in expected["not_in"]]
    return actual_s == str(expected)


def row_matches(row: Dict[str, str], select: Dict[str, Any]) -> bool:
    for k, expected in select.items():
        if k == "regex":
            for col, pat in expected.items():
                if re.search(str(pat), str(row.get(col, ""))) is None:
                    return False
            continue
        if k == "not_regex":
            for col, pat in expected.items():
                if re.search(str(pat), str(row.get(col, ""))) is not None:
                    return False
            continue
        if not value_matches(row.get(k, ""), expected):
            return False
    return True


def select_rows(rows: List[Dict[str, str]], select: Dict[str, Any]) -> List[Dict[str, str]]:
    return [r for r in rows if row_matches(r, select)]


# -----------------------------------------------------------------------------
# Barlow test and entries
# -----------------------------------------------------------------------------

class Entry:
    __slots__ = (
        "row", "raw_signed", "raw_abs", "used_signed", "used_abs", "direction",
        "barlow", "barlow_pass", "barlow_note", "zeroed_by_barlow"
    )

    def __init__(self, row: Dict[str, str]):
        self.row = row
        self.raw_signed = to_float(row.get("diff_pct", ""))
        self.raw_abs = to_float(row.get("abs_pct", ""))
        if not isfinite(self.raw_abs) and isfinite(self.raw_signed):
            self.raw_abs = abs(self.raw_signed)
        self.used_signed = self.raw_signed
        self.used_abs = self.raw_abs
        self.direction = str(row.get("direction", ""))
        self.barlow = math.nan
        self.barlow_pass = ""
        self.barlow_note = ""
        self.zeroed_by_barlow = False


def compute_barlow(row: Dict[str, str], cfg: Dict[str, Any]) -> Tuple[float, Optional[bool], str]:
    """
    Return (barlow_value, pass_or_None, note).
    None means Barlow could not be evaluated and the row should be kept unless configured otherwise.
    """
    diff = to_float(row.get("diff", ""))
    nom_stat = to_float(row.get("nominal_stat_abs", ""))
    var_stat = to_float(row.get("variation_stat_abs", ""))
    if not isfinite(diff):
        return math.nan, None, "missing_diff"
    if not isfinite(nom_stat) or not isfinite(var_stat):
        return math.nan, None, "missing_stat"

    mode = str(cfg.get("mode", "correlated")).lower()
    floor = float(cfg.get("denominator_floor", 1e-30))

    if mode == "independent":
        denom2 = nom_stat * nom_stat + var_stat * var_stat
    elif mode == "correlated":
        denom2 = abs(var_stat * var_stat - nom_stat * nom_stat)
    else:
        raise ValueError(f"Unknown Barlow mode: {mode}")

    if denom2 <= floor * floor:
        if abs(diff) <= floor:
            return 0.0, True, "zero_denom_zero_diff"
        return math.inf, True, "zero_denom_nonzero_diff"

    b = abs(diff) / math.sqrt(denom2)
    threshold = float(cfg.get("threshold", 1.0))
    return b, bool(b >= threshold), "ok"


def apply_barlow_subtract_stat_diff(e: Entry, row: Dict[str, str], cfg: Dict[str, Any]) -> Tuple[bool, str]:
    """
    Apply the presentation-style Barlow correction in absolute units and convert
    the corrected result back to percent for the combiner.

      sigma_sys = 0                                 if |Delta| <= sigma_stat,diff
      sigma_sys = sqrt(Delta^2 - sigma_stat,diff^2) if |Delta| >  sigma_stat,diff

    In the default "correlated" mode this uses the convention shown on the slide,

      sigma_stat,diff = sqrt(abs(var_stat * var_stat - nom_stat * nom_stat)).

    The older "independent" option is also supported and uses
      sigma_stat,diff^2 = sigma_stat,var^2 + sigma_stat,nom^2.
    """
    diff = to_float(row.get("diff", ""))
    nominal = abs(to_float(row.get("nominal_value", "")))
    nom_stat = to_float(row.get("nominal_stat_abs", ""))
    var_stat = to_float(row.get("variation_stat_abs", ""))
    floor = float(cfg.get("denominator_floor", 1e-30))

    if not isfinite(diff):
        return False, "barlow_subtract_missing_diff"
    if not isfinite(nominal) or nominal <= 0.0:
        return False, "barlow_subtract_missing_nominal"
    if not isfinite(nom_stat) or not isfinite(var_stat):
        return False, "barlow_subtract_missing_stat"

    mode = str(cfg.get("mode", "correlated")).lower()
    if mode == "independent":
        stat_diff2 = nom_stat * nom_stat + var_stat * var_stat
        stat_note = "independent"
    elif mode == "correlated":
        stat_diff2 = abs(var_stat * var_stat - nom_stat * nom_stat)
        stat_note = "correlated_slide_max0"
    else:
        raise ValueError(f"Unknown Barlow mode: {mode}")

    abs_diff = abs(diff)
    if stat_diff2 <= floor * floor:
        stat_diff = 0.0
        corrected_abs = abs_diff
        passed = True
        bvalue = math.inf if abs_diff > floor else 0.0
        note = f"barlow_subtract_zero_statdiff:{stat_note}"
    else:
        stat_diff = math.sqrt(stat_diff2)
        bvalue = abs_diff / stat_diff
        passed = bool(abs_diff > stat_diff)
        if passed:
            corrected_abs = math.sqrt(max(0.0, diff * diff - stat_diff2))
            note = f"barlow_subtract_applied:{stat_note}"
        else:
            corrected_abs = 0.0
            note = f"barlow_subtract_zeroed:{stat_note}"

    corrected_pct = 100.0 * corrected_abs / nominal
    sign = 1.0 if diff >= 0.0 else -1.0

    e.used_abs = corrected_pct
    e.used_signed = sign * corrected_pct
    e.barlow = bvalue
    e.barlow_pass = "1" if passed else "0"
    e.zeroed_by_barlow = (corrected_abs == 0.0 and abs_diff > 0.0)

    # Keep the absolute quantities in the details note so the output TSV can be
    # audited without re-reading the input ROOT histograms.
    return True, f"{note};abs_delta={abs_diff:.10g};stat_diff={stat_diff:.10g};used_abs={corrected_abs:.10g}"


def barlow_applies(row: Dict[str, str], global_cfg: Dict[str, Any], group_cfg: Dict[str, Any]) -> bool:
    if not bool(global_cfg.get("enabled", False)):
        return False
    if "apply_barlow" in group_cfg:
        if not bool(group_cfg.get("apply_barlow")):
            return False
    elif not bool(global_cfg.get("default_apply_to_groups", True)):
        return False

    ctype = str(row.get("component_type", ""))
    allowed_types = [str(x) for x in global_cfg.get("apply_to_component_types", ["variation"])]
    return ctype in allowed_types


def make_entries(rows: List[Dict[str, str]], global_barlow: Dict[str, Any], group_cfg: Dict[str, Any]) -> List[Entry]:
    entries: List[Entry] = []
    action = str(global_barlow.get("action", "flag_only")).lower()
    missing_action = str(global_barlow.get("missing_stat_action", "keep")).lower()

    for row in rows:
        e = Entry(row)
        if barlow_applies(row, global_barlow, group_cfg):
            b, passed, note = compute_barlow(row, global_barlow)
            e.barlow = b
            e.barlow_pass = "" if passed is None else ("1" if passed else "0")
            e.barlow_note = note
            if action in ("subtract_stat_diff", "barlow_subtract", "presentation"):
                ok, extra_note = apply_barlow_subtract_stat_diff(e, row, global_barlow)
                if extra_note:
                    e.barlow_note = f"{e.barlow_note};{extra_note}" if e.barlow_note else extra_note
                if not ok and missing_action == "zero":
                    e.used_signed = 0.0
                    e.used_abs = 0.0
                    e.zeroed_by_barlow = True
            elif passed is False and action == "zero_if_not_significant":
                e.used_signed = 0.0
                e.used_abs = 0.0
                e.zeroed_by_barlow = True
            elif passed is None and missing_action == "zero":
                e.used_signed = 0.0
                e.used_abs = 0.0
                e.zeroed_by_barlow = True
        entries.append(e)
    return entries


# -----------------------------------------------------------------------------
# Combination methods
# -----------------------------------------------------------------------------

def entry_up_down(e: Entry) -> Tuple[float, float]:
    """Return positive magnitudes (down, up)."""
    direction = str(e.direction or "").lower()
    abs_v = e.used_abs if isfinite(e.used_abs) else (abs(e.used_signed) if isfinite(e.used_signed) else math.nan)
    signed = e.used_signed

    if direction in ("symmetric", "both"):
        v = abs_v if isfinite(abs_v) else 0.0
        return v, v
    if direction == "up":
        v = abs_v if isfinite(abs_v) else (abs(signed) if isfinite(signed) else 0.0)
        return 0.0, v
    if direction == "down":
        v = abs_v if isfinite(abs_v) else (abs(signed) if isfinite(signed) else 0.0)
        return v, 0.0

    if isfinite(signed):
        if signed > 0:
            return 0.0, abs(signed)
        if signed < 0:
            return abs(signed), 0.0
        return 0.0, 0.0

    if isfinite(abs_v):
        return abs_v, abs_v
    return 0.0, 0.0


def combine_envelope(entries: List[Entry]) -> Tuple[float, float, str]:
    downs, ups = [], []
    for e in entries:
        d, u = entry_up_down(e)
        downs.append(d)
        ups.append(u)
    return (max(downs) if downs else math.nan, max(ups) if ups else math.nan, "envelope")


def combine_max_abs(entries: List[Entry]) -> Tuple[float, float, str]:
    vals = [e.used_abs for e in entries if isfinite(e.used_abs)]
    if not vals:
        vals = [abs(e.used_signed) for e in entries if isfinite(e.used_signed)]
    m = max(vals) if vals else math.nan
    return m, m, "max_abs"


def combine_take_rms(entries: List[Entry]) -> Tuple[float, float, str]:
    # Binning rows should already contain one RMS value per bin/source.
    vals = [e.used_abs for e in entries if isfinite(e.used_abs)]
    if not vals:
        vals = [abs(e.used_signed) for e in entries if isfinite(e.used_signed)]
    if not vals:
        return math.nan, math.nan, "take_rms_no_values"
    # Use maximum to be safe if duplicate RMS rows exist for the same bin.
    m = max(vals)
    return m, m, "take_rms"


def combine_rms(entries: List[Entry]) -> Tuple[float, float, str]:
    vals = []
    for e in entries:
        if isfinite(e.used_abs):
            vals.append(e.used_abs)
        elif isfinite(e.used_signed):
            vals.append(abs(e.used_signed))
    if not vals:
        return math.nan, math.nan, "rms_no_values"
    r = math.sqrt(sum(v * v for v in vals) / len(vals))
    return r, r, "rms"


def combine_quadrature(entries: List[Entry]) -> Tuple[float, float, str]:
    d2 = 0.0
    u2 = 0.0
    any_v = False
    for e in entries:
        d, u = entry_up_down(e)
        d2 += d * d
        u2 += u * u
        any_v = True
    if not any_v:
        return math.nan, math.nan, "quadrature_no_values"
    return math.sqrt(d2), math.sqrt(u2), "quadrature"


def combine_linear(entries: List[Entry]) -> Tuple[float, float, str]:
    dsum = 0.0
    usum = 0.0
    any_v = False
    for e in entries:
        d, u = entry_up_down(e)
        dsum += d
        usum += u
        any_v = True
    if not any_v:
        return math.nan, math.nan, "linear_no_values"
    return dsum, usum, "linear"


def combine_pair_envelope_then_quadrature(entries: List[Entry], pair_by: str) -> Tuple[float, float, str]:
    grouped: Dict[str, List[Entry]] = defaultdict(list)
    for e in entries:
        grouped[str(e.row.get(pair_by, ""))].append(e)
    d2 = 0.0
    u2 = 0.0
    n_pairs = 0
    for _, sub in grouped.items():
        d, u, _ = combine_envelope(sub)
        if isfinite(d) or isfinite(u):
            d = d if isfinite(d) else 0.0
            u = u if isfinite(u) else 0.0
            d2 += d * d
            u2 += u * u
            n_pairs += 1
    if n_pairs == 0:
        return math.nan, math.nan, "pair_envelope_then_quadrature_no_pairs"
    return math.sqrt(d2), math.sqrt(u2), f"pair_envelope_then_quadrature:n_pairs={n_pairs}"


def combine_entries(entries: List[Entry], group_cfg: Dict[str, Any]) -> Tuple[float, float, float, str]:
    method = str(group_cfg.get("combine_variations", "envelope")).lower()
    if not entries:
        return math.nan, math.nan, math.nan, "no_entries"

    if method == "envelope":
        down, up, note = combine_envelope(entries)
    elif method == "max_abs":
        down, up, note = combine_max_abs(entries)
    elif method == "take_rms":
        down, up, note = combine_take_rms(entries)
    elif method == "rms":
        down, up, note = combine_rms(entries)
    elif method == "quadrature":
        down, up, note = combine_quadrature(entries)
    elif method in ("linear", "linear_sum", "correlated_linear"):
        down, up, note = combine_linear(entries)
    elif method == "pair_envelope_then_quadrature":
        pair_by = str(group_cfg.get("pair_by", "source_name"))
        down, up, note = combine_pair_envelope_then_quadrature(entries, pair_by)
    else:
        raise ValueError(f"Unknown combine_variations method: {method}")

    sym = max([v for v in [down, up] if isfinite(v)], default=math.nan)

    # Optional per-group symmetrization.  This is useful for GUI-driven studies
    # where the user wants a given source to enter the final table as a symmetric
    # uncertainty even if its raw envelope is asymmetric.
    sym_mode = str(group_cfg.get("symmetrize", "none")).lower()
    if sym_mode in ("1", "true", "yes", "y", "max", "max_abs", "symmetric"):
        down = sym
        up = sym
        note = f"{note};symmetrize=max_abs"
    elif sym_mode in ("average", "avg", "average_abs", "mean_abs"):
        vals = [v for v in (down, up) if isfinite(v)]
        avg = sum(vals) / len(vals) if vals else math.nan
        down = avg
        up = avg
        sym = avg
        note = f"{note};symmetrize=average_abs"

    return down, up, sym, note


def combine_group_totals(group_values: Dict[str, Tuple[float, float, float]], include_groups: List[str], method: str) -> Tuple[float, float, float]:
    method = str(method or "quadrature").lower()
    vals = []
    for name in include_groups:
        if name not in group_values:
            continue
        vals.append(group_values[name])

    if not vals:
        return math.nan, math.nan, math.nan

    if method == "quadrature":
        d2 = sum((d if isfinite(d) else 0.0) ** 2 for d, _, _ in vals)
        u2 = sum((u if isfinite(u) else 0.0) ** 2 for _, u, _ in vals)
        down = math.sqrt(d2)
        up = math.sqrt(u2)
    elif method in ("linear", "linear_sum", "correlated_linear"):
        down = sum((d if isfinite(d) else 0.0) for d, _, _ in vals)
        up = sum((u if isfinite(u) else 0.0) for _, u, _ in vals)
    elif method == "envelope":
        down = max((d if isfinite(d) else 0.0) for d, _, _ in vals)
        up = max((u if isfinite(u) else 0.0) for _, u, _ in vals)
    else:
        raise ValueError(f"Unknown total.combine_groups method: {method}")

    sym = max(down, up)
    return down, up, sym


def _rule_total_mode(rule: Dict[str, Any], default_mode: str = "quadrature") -> str:
    mode = str(rule.get("total_mode", "") or "").strip().lower()
    if not mode:
        corr = str(rule.get("correlation", "") or "").strip().lower()
        if corr in ("correlated", "corr", "linear"):
            mode = "linear"
        elif corr in ("uncorrelated", "uncorr", "quadrature", "quad"):
            mode = "quadrature"
        elif corr in ("envelope", "max"):
            mode = "envelope"
        elif corr in ("exclude", "none", "off"):
            mode = "exclude"
        else:
            mode = default_mode
    aliases = {
        "uncorrelated": "quadrature",
        "uncorr": "quadrature",
        "quad": "quadrature",
        "correlated": "linear",
        "corr": "linear",
        "linear_sum": "linear",
        "off": "exclude",
        "none": "exclude",
    }
    return aliases.get(mode, mode)


def combine_group_totals_mixed(
    group_values: Dict[str, Tuple[float, float, float]],
    group_rules: List[Dict[str, Any]],
    include_groups: List[str],
    total_cfg: Dict[str, Any],
) -> Tuple[float, float, float]:
    """
    Mixed total combination used by the GUI.

    Per group, total_mode/correlation controls how that already-combined group
    enters the final total:
      * quadrature / uncorrelated: sum in quadrature
      * linear / correlated: add linearly
      * envelope: take the maximum envelope bucket
      * exclude: ignore, even if present in include_groups

    The final convention is:
      total = linear_bucket + sqrt(quadrature_bucket^2 + envelope_bucket^2)
    separately for down and up. This keeps correlated sources conservative while
    still treating the remaining buckets as independent.
    """
    include_set = set(include_groups)
    default_mode = str(total_cfg.get("default_group_mode", "quadrature"))
    rule_by_name = {str(g.get("name", "")): g for g in group_rules}

    qd2 = qu2 = 0.0
    ld = lu = 0.0
    ed = eu = 0.0
    any_v = False

    for name in include_groups:
        if name not in group_values:
            continue
        d, u, _ = group_values[name]
        d = d if isfinite(d) else 0.0
        u = u if isfinite(u) else 0.0
        rule = rule_by_name.get(name, {})
        mode = _rule_total_mode(rule, default_mode)
        if mode == "exclude":
            continue
        if mode == "quadrature":
            qd2 += d * d
            qu2 += u * u
            any_v = True
        elif mode == "linear":
            ld += d
            lu += u
            any_v = True
        elif mode == "envelope":
            ed = max(ed, d)
            eu = max(eu, u)
            any_v = True
        else:
            raise ValueError(f"Unknown group total_mode/correlation for {name}: {mode}")

    if not any_v:
        return math.nan, math.nan, math.nan

    down = ld + math.sqrt(qd2 + ed * ed)
    up = lu + math.sqrt(qu2 + eu * eu)
    return down, up, max(down, up)


# -----------------------------------------------------------------------------
# Details output
# -----------------------------------------------------------------------------

DETAIL_COLUMNS = [
    *KEY_COLS,
    "group_rule",
    "source_group",
    "source_name",
    "variation_name",
    "variation_code",
    "component_type",
    "direction",
    "nominal_value",
    "nominal_stat_abs",
    "variation_value",
    "variation_stat_abs",
    "raw_diff_pct",
    "raw_abs_pct",
    "used_diff_pct",
    "used_abs_pct",
    "barlow_value",
    "barlow_pass",
    "barlow_note",
    "zeroed_by_barlow",
]


def make_detail_row(key: Tuple[str, ...], group_name: str, e: Entry) -> Dict[str, Any]:
    row = {c: key[i] for i, c in enumerate(KEY_COLS)}
    r = e.row
    row.update({
        "group_rule": group_name,
        "source_group": r.get("source_group", ""),
        "source_name": r.get("source_name", ""),
        "variation_name": r.get("variation_name", ""),
        "variation_code": r.get("variation_code", ""),
        "component_type": r.get("component_type", ""),
        "direction": r.get("direction", ""),
        "nominal_value": r.get("nominal_value", ""),
        "nominal_stat_abs": r.get("nominal_stat_abs", ""),
        "variation_value": r.get("variation_value", ""),
        "variation_stat_abs": r.get("variation_stat_abs", ""),
        "raw_diff_pct": e.raw_signed,
        "raw_abs_pct": e.raw_abs,
        "used_diff_pct": e.used_signed,
        "used_abs_pct": e.used_abs,
        "barlow_value": e.barlow,
        "barlow_pass": e.barlow_pass,
        "barlow_note": e.barlow_note,
        "zeroed_by_barlow": "1" if e.zeroed_by_barlow else "0",
    })
    return row


# -----------------------------------------------------------------------------
# Main combiner
# -----------------------------------------------------------------------------

def build_final_table(config: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[str]]:
    project_dir = resolve_path(str(config.get("project_dir", ".")), Path.cwd()).resolve()
    input_path = resolve_path(str(config.get("input", "Systematics/systematics_components_v15.tsv")), project_dir)

    if not input_path.exists():
        raise FileNotFoundError(f"Missing input table: {input_path}")

    _, bins = read_components(input_path)
    info(f"Read bins: {len(bins)} from {input_path}")

    group_rules = list(config.get("groups", []))
    if not group_rules:
        raise ValueError("No group rules found in config['groups'].")

    global_barlow = dict(config.get("barlow", {}))
    total_cfg = dict(config.get("total", {}))
    include_groups = [str(x) for x in total_cfg.get("include_groups", [g.get("name") for g in group_rules])]
    total_method = str(total_cfg.get("combine_groups", "quadrature"))
    include_stat = bool(total_cfg.get("include_stat", True))

    output_rows: List[Dict[str, Any]] = []
    detail_rows: List[Dict[str, Any]] = []
    group_names: List[str] = []

    for g in group_rules:
        name = str(g.get("name", "")).strip()
        if not name:
            raise ValueError("Each group rule needs a non-empty 'name'.")
        if name not in group_names:
            group_names.append(name)

    for key in sorted(bins.keys(), key=lambda k: (k[0], int(float(k[4])) if str(k[4]).replace('.', '', 1).isdigit() else k[4], int(float(k[9])) if str(k[9]).replace('.', '', 1).isdigit() else k[9])):
        rows = bins[key]
        out = base_info_from_rows(key, rows)
        nominal_value = to_float(out.get("nominal_value", ""))
        stat_pct = to_float(out.get("stat_pct", out.get("nominal_stat_pct", "")))

        group_values: Dict[str, Tuple[float, float, float]] = {}
        notes: List[str] = []

        for g in group_rules:
            name = str(g.get("name", ""))
            selected = select_rows(rows, dict(g.get("select", {})))
            entries = make_entries(selected, global_barlow, g)
            down, up, sym, note = combine_entries(entries, g)
            group_values[name] = (down, up, sym)

            out[f"{safe_name(name)}_down_pct"] = down
            out[f"{safe_name(name)}_up_pct"] = up
            out[f"{safe_name(name)}_sym_pct"] = sym
            out[f"{safe_name(name)}_n"] = len(entries)
            out[f"{safe_name(name)}_barlow_zeroed_n"] = sum(1 for e in entries if e.zeroed_by_barlow)
            out[f"{safe_name(name)}_barlow_failed_n"] = sum(1 for e in entries if e.barlow_pass == "0")

            if note and note != "no_entries":
                notes.append(f"{name}:{note}")
            if not entries:
                notes.append(f"{name}:no_entries")

            if bool(config.get("write_details", True)):
                for e in entries:
                    detail_rows.append(make_detail_row(key, name, e))

        if total_method.lower() == "mixed":
            syst_down, syst_up, syst_sym = combine_group_totals_mixed(group_values, group_rules, include_groups, total_cfg)
        else:
            syst_down, syst_up, syst_sym = combine_group_totals(group_values, include_groups, total_method)
        out["total_syst_down_pct"] = syst_down
        out["total_syst_up_pct"] = syst_up
        out["total_syst_sym_pct"] = syst_sym
        out["total_syst_down_abs"] = pct_to_abs(syst_down, nominal_value)
        out["total_syst_up_abs"] = pct_to_abs(syst_up, nominal_value)
        out["total_syst_sym_abs"] = pct_to_abs(syst_sym, nominal_value)

        if include_stat and isfinite(stat_pct):
            out["total_with_stat_down_pct"] = math.sqrt((syst_down if isfinite(syst_down) else 0.0) ** 2 + stat_pct ** 2)
            out["total_with_stat_up_pct"] = math.sqrt((syst_up if isfinite(syst_up) else 0.0) ** 2 + stat_pct ** 2)
            out["total_with_stat_sym_pct"] = math.sqrt((syst_sym if isfinite(syst_sym) else 0.0) ** 2 + stat_pct ** 2)
        else:
            out["total_with_stat_down_pct"] = syst_down
            out["total_with_stat_up_pct"] = syst_up
            out["total_with_stat_sym_pct"] = syst_sym

        out["total_with_stat_down_abs"] = pct_to_abs(to_float(out["total_with_stat_down_pct"]), nominal_value)
        out["total_with_stat_up_abs"] = pct_to_abs(to_float(out["total_with_stat_up_pct"]), nominal_value)
        out["total_with_stat_sym_abs"] = pct_to_abs(to_float(out["total_with_stat_sym_pct"]), nominal_value)
        out["combination_notes"] = "; ".join(notes)
        output_rows.append(out)

    return output_rows, detail_rows, group_names


def build_output_columns(rows: List[Dict[str, Any]], group_names: List[str]) -> List[str]:
    base = [*KEY_COLS, *BASE_NUM_COLS, "stat_abs", "stat_pct"]
    group_cols: List[str] = []
    for name in group_names:
        s = safe_name(name)
        group_cols.extend([
            f"{s}_down_pct",
            f"{s}_up_pct",
            f"{s}_sym_pct",
            f"{s}_n",
            f"{s}_barlow_zeroed_n",
            f"{s}_barlow_failed_n",
        ])
    total_cols = [
        "total_syst_down_pct",
        "total_syst_up_pct",
        "total_syst_sym_pct",
        "total_syst_down_abs",
        "total_syst_up_abs",
        "total_syst_sym_abs",
        "total_with_stat_down_pct",
        "total_with_stat_up_pct",
        "total_with_stat_sym_pct",
        "total_with_stat_down_abs",
        "total_with_stat_up_abs",
        "total_with_stat_sym_abs",
        "combination_notes",
    ]
    cols = [c for c in [*base, *group_cols, *total_cols] if c]
    # Preserve any surprise columns at the end.
    seen = set(cols)
    for r in rows:
        for k in r:
            if k not in seen:
                cols.append(k)
                seen.add(k)
    return cols


def write_table(path: Path, rows: List[Dict[str, Any]], columns: List[str], delimiter: str = "\t") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow({c: fmt(r.get(c, "")) for c in columns})



# -----------------------------------------------------------------------------
# Final plotting
# -----------------------------------------------------------------------------

def plot_safe_filename(*parts: Any) -> str:
    name = "_".join(safe_name(p) for p in parts if str(p or "").strip())
    return name or "plot"


def default_x_title(row: Dict[str, Any]) -> str:
    obs = str(row.get("observable_pretty") or row.get("observable") or row.get("job_label") or "x")
    if str(row.get("job_label", "")).startswith("RCP_"):
        # Remove leading R_CP(...) wording only if the user later wants a custom map.
        return obs
    return obs


def default_y_title(row: Dict[str, Any]) -> str:
    job = str(row.get("job_label", ""))
    if job.startswith("RCP_"):
        return "R_CP"
    return "Value"


def get_axis_title(row: Dict[str, Any], plots_cfg: Dict[str, Any], which: str) -> str:
    job = str(row.get("job_label", ""))
    obs = str(row.get("observable", ""))
    key = "x_titles" if which == "x" else "y_titles"
    mapping = dict(plots_cfg.get(key, {}))
    if job in mapping:
        return str(mapping[job])
    if obs in mapping:
        return str(mapping[obs])
    return default_x_title(row) if which == "x" else default_y_title(row)


def cfg_bool(value: Any, default: bool = False) -> bool:
    """Boolean helper accepting JSON booleans and common string values."""
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def parse_ylim(value: Any) -> Tuple[Optional[float], Optional[float]]:
    """Return (ymin, ymax). Empty string / null means automatic side."""
    if value is None or value == "":
        return (None, None)
    if isinstance(value, str):
        parts = [x.strip() for x in re.split(r"[,;: ]+", value) if x.strip()]
    elif isinstance(value, (list, tuple)):
        parts = list(value)
    else:
        return (None, None)
    if len(parts) < 2:
        return (None, None)

    out: List[Optional[float]] = []
    for item in parts[:2]:
        if item is None or str(item).strip() == "":
            out.append(None)
            continue
        try:
            f = float(item)
            out.append(f if math.isfinite(f) else None)
        except Exception:
            out.append(None)
    return (out[0], out[1])


def get_plot_ylim(row: Dict[str, Any], plots_cfg: Dict[str, Any]) -> Tuple[Optional[float], Optional[float]]:
    """Resolve y-range for the current final plot.

    v5 default behaviour is intentionally conservative and hard to miss:
      * all R_CP plots use y in [0, 2] unless a more specific config value is given,
      * ordinary spectra are left automatic, unless spectra_ylim/ylim is set.

    Precedence:
      1) ylim_by_plot["job_label|centrality_label"] or ["job_label|centrality"]
      2) ylim_by_job[job_label]
      3) ylim_by_observable[observable]
      4) rcp_ylim_by_job[job_label] for RCP plots
      5) rcp_ylim for all RCP plots
      6) v5 hard default [0, 2] for RCP plots
      7) spectra_ylim for all non-RCP spectra
      8) ylim for all plots
    """
    job = str(row.get("job_label", ""))
    obs = str(row.get("observable", ""))
    cent = str(row.get("centrality", ""))
    cent_label = str(row.get("centrality_label", ""))
    is_rcp = job.startswith("RCP_")

    def nonempty_ylim(value: Any) -> Tuple[Optional[float], Optional[float]]:
        ymin, ymax = parse_ylim(value)
        if ymin is not None or ymax is not None:
            return ymin, ymax
        return (None, None)

    by_plot = dict(plots_cfg.get("ylim_by_plot", {}))
    for key in (f"{job}|{cent_label}", f"{job}|{cent}"):
        if key in by_plot:
            ymin, ymax = nonempty_ylim(by_plot[key])
            if ymin is not None or ymax is not None:
                return ymin, ymax

    by_job = dict(plots_cfg.get("ylim_by_job", {}))
    if job in by_job:
        ymin, ymax = nonempty_ylim(by_job[job])
        if ymin is not None or ymax is not None:
            return ymin, ymax

    by_obs = dict(plots_cfg.get("ylim_by_observable", {}))
    if obs in by_obs:
        ymin, ymax = nonempty_ylim(by_obs[obs])
        if ymin is not None or ymax is not None:
            return ymin, ymax

    if is_rcp:
        rcp_by_job = dict(plots_cfg.get("rcp_ylim_by_job", {}))
        if job in rcp_by_job:
            ymin, ymax = nonempty_ylim(rcp_by_job[job])
            if ymin is not None or ymax is not None:
                return ymin, ymax
        ymin, ymax = nonempty_ylim(plots_cfg.get("rcp_ylim", None))
        if ymin is not None or ymax is not None:
            return ymin, ymax
        # Hard v5 default requested for final R_CP plots.
        if cfg_bool(plots_cfg.get("force_default_rcp_ylim", True), True):
            return (0.0, 2.0)
        return (None, None)

    ymin, ymax = nonempty_ylim(plots_cfg.get("spectra_ylim", None))
    if ymin is not None or ymax is not None:
        return ymin, ymax
    return parse_ylim(plots_cfg.get("ylim", None))


def get_plot_logy(row: Dict[str, Any], plots_cfg: Dict[str, Any]) -> bool:
    """Resolve log-y setting for spectra/RCP final plots.

    v5 default behaviour:
      * spectra: log-y ON by default,
      * R_CP: log-y OFF by default.

    The optional force_default_logy_spectra flag keeps spectra logarithmic even
    if an older GUI config still contains "logy_spectra": false.
    """
    job = str(row.get("job_label", ""))
    is_rcp = job.startswith("RCP_")

    # In v5, R_CP is linear by default even if an old config contains
    # the legacy global "logy": true. Use logy_rcp=true only if you really
    # want logarithmic R_CP plots.
    if is_rcp:
        return cfg_bool(plots_cfg.get("logy_rcp", False), False)

    if cfg_bool(plots_cfg.get("force_default_logy_spectra", True), True):
        return True
    if cfg_bool(plots_cfg.get("logy_spectra", True), True):
        return True
    return cfg_bool(plots_cfg.get("logy", False), False)

def final_plot_groups(rows: List[Dict[str, Any]]) -> Dict[Tuple[str, str], List[Dict[str, Any]]]:
    groups: Dict[Tuple[str, str], List[Dict[str, Any]]] = defaultdict(list)
    for r in rows:
        groups[(str(r.get("job_label", "")), str(r.get("centrality", "")))].append(r)
    return groups


def row_sort_bin(row: Dict[str, Any]) -> Tuple[int, float]:
    b = to_float(row.get("bin", ""), math.nan)
    x = to_float(row.get("bin_low", ""), math.nan)
    return (int(b) if isfinite(b) else 10**9, x if isfinite(x) else 0.0)


def plot_one_final_result(ax: Any, subrows: List[Dict[str, Any]], plots_cfg: Dict[str, Any]) -> None:
    import numpy as _np
    from matplotlib.patches import Rectangle

    subrows = sorted(subrows, key=row_sort_bin)
    first = subrows[0]

    x_centers: List[float] = []
    x_lows: List[float] = []
    x_highs: List[float] = []
    y_vals: List[float] = []
    stat_vals: List[float] = []
    syst_down_vals: List[float] = []
    syst_up_vals: List[float] = []

    for i, r in enumerate(subrows, start=1):
        lo = to_float(r.get("bin_low", ""), math.nan)
        hi = to_float(r.get("bin_high", ""), math.nan)
        if not isfinite(lo) or not isfinite(hi) or hi <= lo:
            lo = float(i) - 0.4
            hi = float(i) + 0.4
        x_lows.append(lo)
        x_highs.append(hi)
        x_centers.append(0.5 * (lo + hi))
        y_vals.append(to_float(r.get("nominal_value", ""), math.nan))
        stat_vals.append(abs(to_float(r.get("stat_abs", r.get("nominal_stat_abs", "")), math.nan)))
        syst_down_vals.append(abs(to_float(r.get("total_syst_down_abs", ""), math.nan)))
        syst_up_vals.append(abs(to_float(r.get("total_syst_up_abs", ""), math.nan)))

    draw_syst = bool(plots_cfg.get("draw_syst_boxes", True))
    draw_stat = bool(plots_cfg.get("draw_stat_errors", True))
    draw_points = bool(plots_cfg.get("draw_points", True))
    stat_label = str(plots_cfg.get("stat_label", "stat."))
    syst_label = str(plots_cfg.get("syst_label", "syst."))

    # Draw systematic uncertainty boxes first, so points/stat bars remain visible.
    syst_handle = None
    if draw_syst:
        for lo, hi, y, yd, yu in zip(x_lows, x_highs, y_vals, syst_down_vals, syst_up_vals):
            if not (isfinite(y) and isfinite(yd) and isfinite(yu)):
                continue
            height = yd + yu
            if height < 0:
                continue
            rect = Rectangle(
                (lo, y - yd),
                hi - lo,
                height,
                facecolor="0.75",
                edgecolor="0.35",
                alpha=float(plots_cfg.get("syst_box_alpha", 0.35)),
                linewidth=float(plots_cfg.get("syst_box_linewidth", 0.8)),
            )
            ax.add_patch(rect)
            if syst_handle is None:
                syst_handle = rect

    # Statistical errors and central values.
    valid = [i for i, (x, y) in enumerate(zip(x_centers, y_vals)) if isfinite(x) and isfinite(y)]
    if valid:
        x = _np.asarray([x_centers[i] for i in valid], dtype=float)
        y = _np.asarray([y_vals[i] for i in valid], dtype=float)
        stat = _np.asarray([stat_vals[i] if isfinite(stat_vals[i]) else 0.0 for i in valid], dtype=float)
        if draw_stat:
            stat_handle = ax.errorbar(
                x,
                y,
                yerr=stat,
                fmt="o" if draw_points else "none",
                markersize=float(plots_cfg.get("marker_size", 4.0)),
                linewidth=float(plots_cfg.get("stat_linewidth", 1.0)),
                capsize=float(plots_cfg.get("stat_capsize", 2.0)),
                label=stat_label,
                color="black",
            )
        elif draw_points:
            stat_handle = ax.plot(x, y, "o", markersize=float(plots_cfg.get("marker_size", 4.0)), color="black", label="value")
        else:
            stat_handle = None
    else:
        stat_handle = None

    if bool(plots_cfg.get("draw_unity_for_rcp", True)) and str(first.get("job_label", "")).startswith("RCP_"):
        ax.axhline(1.0, linestyle="--", linewidth=1.0, color="0.4")

    title_template = str(plots_cfg.get("title_template", "{job_label}, {centrality_label}"))
    title = title_template.format(
        job_label=first.get("job_label", ""),
        observable=first.get("observable", ""),
        observable_pretty=first.get("observable_pretty", ""),
        centrality=first.get("centrality", ""),
        centrality_label=first.get("centrality_label", ""),
        method=first.get("method", ""),
    )
    ax.set_title(title)
    ax.set_xlabel(get_axis_title(first, plots_cfg, "x"))
    ax.set_ylabel(get_axis_title(first, plots_cfg, "y"))

    use_logy = get_plot_logy(first, plots_cfg)
    if use_logy:
        positive = [y for y in y_vals if isfinite(y) and y > 0]
        if positive:
            ax.set_yscale("log")

    ymin, ymax = get_plot_ylim(first, plots_cfg)
    if ymin is not None or ymax is not None:
        # For log-y plots, matplotlib requires positive lower bounds.
        if use_logy and ymin is not None and ymin <= 0.0:
            positive_candidates = []
            for y, yd in zip(y_vals, syst_down_vals):
                if isfinite(y) and y > 0:
                    low = y - yd if isfinite(yd) else y
                    if low > 0:
                        positive_candidates.append(low)
                    positive_candidates.append(y)
            ymin = 0.5 * min(positive_candidates) if positive_candidates else None
        ax.set_ylim(bottom=ymin, top=ymax)

    ax.grid(bool(plots_cfg.get("grid", True)), alpha=0.25)
    ax.margins(x=0.03)

    handles = []
    labels = []
    if syst_handle is not None:
        handles.append(syst_handle)
        labels.append(syst_label)
    if stat_handle is not None:
        # ErrorbarContainer is accepted by legend.
        handles.append(stat_handle)
        labels.append(stat_label if draw_stat else "value")
    if handles and bool(plots_cfg.get("legend", True)):
        ax.legend(handles, labels, frameon=False, fontsize=str(plots_cfg.get("legend_fontsize", "small")))


def write_final_plots(config: Dict[str, Any], rows: List[Dict[str, Any]]) -> None:
    plots_cfg = dict(config.get("plots", {}))
    if not bool(plots_cfg.get("enabled", False)):
        return

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except Exception as e:
        warn(f"Plots requested, but matplotlib could not be imported: {e}")
        return

    project_dir = resolve_path(str(config.get("project_dir", ".")), Path.cwd()).resolve()
    out_dir = resolve_path(str(plots_cfg.get("output_dir", "Systematics/FinalPlots")), project_dir)
    png_dir = out_dir / "png"
    out_dir.mkdir(parents=True, exist_ok=True)
    if bool(plots_cfg.get("write_png", True)):
        png_dir.mkdir(parents=True, exist_ok=True)

    pdf_name = str(plots_cfg.get("pdf", "final_systematics_values.pdf"))
    pdf_path = out_dir / pdf_name
    figsize = plots_cfg.get("figsize", [7.0, 5.0])
    try:
        figsize_tuple = (float(figsize[0]), float(figsize[1]))
    except Exception:
        figsize_tuple = (7.0, 5.0)
    dpi = int(plots_cfg.get("dpi", 150))

    groups = final_plot_groups(rows)
    written_png = 0
    with PdfPages(pdf_path) as pdf:
        for (job_label, cent), subrows in sorted(groups.items(), key=lambda kv: (kv[0][0], to_float(kv[0][1], 9999))):
            if not subrows:
                continue
            fig, ax = plt.subplots(figsize=figsize_tuple)
            plot_one_final_result(ax, subrows, plots_cfg)
            fig.tight_layout()
            pdf.savefig(fig)
            if bool(plots_cfg.get("write_png", True)):
                fname = plot_safe_filename(job_label, "cent", cent) + ".png"
                fig.savefig(png_dir / fname, dpi=dpi)
                written_png += 1
            plt.close(fig)

    info(f"Wrote final plots PDF: {pdf_path}")
    if bool(plots_cfg.get("write_png", True)):
        info(f"Wrote final plot PNGs: {png_dir} ({written_png} files)")


# -----------------------------------------------------------------------------
# Final ROOT output
# -----------------------------------------------------------------------------

def _root_import():
    try:
        import ROOT  # type: ignore
        return ROOT
    except Exception as e:
        warn(f"PyROOT could not be imported: {e}. Trying ROOT command-line fallback.")
        return None


def _rows_to_edges(subrows: List[Dict[str, Any]]) -> List[float]:
    """Return monotonically increasing variable-bin edges for a group of rows."""
    subrows = sorted(subrows, key=row_sort_bin)
    edges: List[float] = []
    for i, r in enumerate(subrows):
        lo = to_float(r.get("bin_low", ""), math.nan)
        hi = to_float(r.get("bin_high", ""), math.nan)
        if not isfinite(lo) or not isfinite(hi) or hi <= lo:
            # Fallback: one artificial unit-width bin per output row.
            lo = float(i)
            hi = float(i + 1)
        if i == 0:
            edges.append(float(lo))
        edges.append(float(hi))

    # TH1D requires strictly increasing edges. If an input is malformed, fall back
    # to a simple sequential axis, but keep the graph x-errors below faithful when possible.
    if len(edges) < 2 or any(edges[i + 1] <= edges[i] for i in range(len(edges) - 1)):
        return [float(i) for i in range(len(subrows) + 1)]
    return edges


def _make_th1(ROOT: Any, name: str, title: str, edges: List[float]) -> Any:
    from array import array
    arr = array('d', [float(x) for x in edges])
    h = ROOT.TH1D(name, title, len(edges) - 1, arr)
    h.SetDirectory(0)
    return h


def _fill_hist(h: Any, values: List[float], errors: Optional[List[float]] = None) -> None:
    for i, v in enumerate(values, start=1):
        h.SetBinContent(i, float(v) if isfinite(v) else 0.0)
        if errors is not None:
            e = errors[i - 1] if i - 1 < len(errors) else 0.0
            h.SetBinError(i, float(e) if isfinite(e) else 0.0)


def _make_graph(ROOT: Any, name: str, title: str, subrows: List[Dict[str, Any]], ey_low_col: str, ey_high_col: str) -> Any:
    g = ROOT.TGraphAsymmErrors(len(subrows))
    g.SetName(name)
    g.SetTitle(title)
    for i, r in enumerate(sorted(subrows, key=row_sort_bin)):
        lo = to_float(r.get("bin_low", ""), math.nan)
        hi = to_float(r.get("bin_high", ""), math.nan)
        if not isfinite(lo) or not isfinite(hi) or hi <= lo:
            lo = float(i)
            hi = float(i + 1)
        x = 0.5 * (lo + hi)
        exl = x - lo
        exh = hi - x
        y = to_float(r.get("nominal_value", ""), 0.0)
        eyl = abs(to_float(r.get(ey_low_col, ""), 0.0))
        eyh = abs(to_float(r.get(ey_high_col, ""), 0.0))
        g.SetPoint(i, x, y)
        g.SetPointError(i, exl, exh, eyl, eyh)
    return g


def _write_final_root_file_pyroot(config: Dict[str, Any], rows: List[Dict[str, Any]], root_path: Path) -> bool:
    """Write a ROOT file with final values and uncertainties using PyROOT.

    The file is intended as the single source for preliminary-result plotting.
    For each (job_label, centrality/R_CP index) it writes:
      * h_value_stat_*              central values with statistical TH1 errors,
      * h_syst_{down,up,sym}_{abs,pct}_*,
      * h_total_with_stat_{down,up,sym}_{abs,pct}_*,
      * g_value_stat_*              TGraphAsymmErrors with stat errors,
      * g_value_syst_*              TGraphAsymmErrors with syst errors,
      * g_value_total_with_stat_*   TGraphAsymmErrors with stat+syst errors.
    """
    ROOT = _root_import()
    if ROOT is None:
        return False

    root_path.parent.mkdir(parents=True, exist_ok=True)

    f = ROOT.TFile.Open(str(root_path), "RECREATE")
    if not f or f.IsZombie():
        warn(f"Could not create ROOT output file with PyROOT: {root_path}")
        return False

    try:
        try:
            ROOT.TH1.AddDirectory(False)
        except Exception:
            pass

        groups = final_plot_groups(rows)
        for (job_label, cent), subrows in sorted(groups.items(), key=lambda kv: (kv[0][0], to_float(kv[0][1], 9999))):
            if not subrows:
                continue
            subrows = sorted(subrows, key=row_sort_bin)
            first = subrows[0]
            job_dir_name = safe_name(job_label)
            job_dir = f.GetDirectory(job_dir_name)
            if not job_dir:
                job_dir = f.mkdir(job_dir_name)
            job_dir.cd()

            cent_label = str(first.get("centrality_label", cent))
            tag = safe_name(cent_label)
            title_base = f"{job_label}, {cent_label}"
            edges = _rows_to_edges(subrows)

            y = [to_float(r.get("nominal_value", ""), 0.0) for r in subrows]
            stat = [abs(to_float(r.get("stat_abs", r.get("nominal_stat_abs", "")), 0.0)) for r in subrows]
            syst_down_abs = [abs(to_float(r.get("total_syst_down_abs", ""), 0.0)) for r in subrows]
            syst_up_abs = [abs(to_float(r.get("total_syst_up_abs", ""), 0.0)) for r in subrows]
            syst_sym_abs = [abs(to_float(r.get("total_syst_sym_abs", ""), 0.0)) for r in subrows]
            syst_down_pct = [abs(to_float(r.get("total_syst_down_pct", ""), 0.0)) for r in subrows]
            syst_up_pct = [abs(to_float(r.get("total_syst_up_pct", ""), 0.0)) for r in subrows]
            syst_sym_pct = [abs(to_float(r.get("total_syst_sym_pct", ""), 0.0)) for r in subrows]
            total_down_abs = [abs(to_float(r.get("total_with_stat_down_abs", ""), 0.0)) for r in subrows]
            total_up_abs = [abs(to_float(r.get("total_with_stat_up_abs", ""), 0.0)) for r in subrows]
            total_sym_abs = [abs(to_float(r.get("total_with_stat_sym_abs", ""), 0.0)) for r in subrows]
            total_down_pct = [abs(to_float(r.get("total_with_stat_down_pct", ""), 0.0)) for r in subrows]
            total_up_pct = [abs(to_float(r.get("total_with_stat_up_pct", ""), 0.0)) for r in subrows]
            total_sym_pct = [abs(to_float(r.get("total_with_stat_sym_pct", ""), 0.0)) for r in subrows]

            hist_specs = [
                (f"h_value_stat_{tag}", y, stat, "central values; stat. as TH1 errors"),
                (f"h_syst_down_abs_{tag}", syst_down_abs, None, "absolute systematic down uncertainty"),
                (f"h_syst_up_abs_{tag}", syst_up_abs, None, "absolute systematic up uncertainty"),
                (f"h_syst_sym_abs_{tag}", syst_sym_abs, None, "absolute symmetric systematic uncertainty"),
                (f"h_syst_down_pct_{tag}", syst_down_pct, None, "percent systematic down uncertainty"),
                (f"h_syst_up_pct_{tag}", syst_up_pct, None, "percent systematic up uncertainty"),
                (f"h_syst_sym_pct_{tag}", syst_sym_pct, None, "percent symmetric systematic uncertainty"),
                (f"h_total_with_stat_down_abs_{tag}", total_down_abs, None, "absolute stat+syst down uncertainty"),
                (f"h_total_with_stat_up_abs_{tag}", total_up_abs, None, "absolute stat+syst up uncertainty"),
                (f"h_total_with_stat_sym_abs_{tag}", total_sym_abs, None, "absolute symmetric stat+syst uncertainty"),
                (f"h_total_with_stat_down_pct_{tag}", total_down_pct, None, "percent stat+syst down uncertainty"),
                (f"h_total_with_stat_up_pct_{tag}", total_up_pct, None, "percent stat+syst up uncertainty"),
                (f"h_total_with_stat_sym_pct_{tag}", total_sym_pct, None, "percent symmetric stat+syst uncertainty"),
            ]
            for hname, values, errors, desc in hist_specs:
                h = _make_th1(ROOT, hname, f"{title_base};bin;{desc}", edges)
                _fill_hist(h, values, errors)
                h.Write()

            g_stat = _make_graph(ROOT, f"g_value_stat_{tag}", f"{title_base};x;value", subrows, "stat_abs", "stat_abs")
            g_syst = _make_graph(ROOT, f"g_value_syst_{tag}", f"{title_base};x;value", subrows, "total_syst_down_abs", "total_syst_up_abs")
            g_total = _make_graph(ROOT, f"g_value_total_with_stat_{tag}", f"{title_base};x;value", subrows, "total_with_stat_down_abs", "total_with_stat_up_abs")
            g_stat.Write()
            g_syst.Write()
            g_total.Write()

            # Also write a small note per directory to document the object tag.
            try:
                note = ROOT.TNamed(f"info_{tag}", f"job_label={job_label}; centrality_label={cent_label}; tag={tag}")
                note.Write()
            except Exception:
                pass

            f.cd()

        try:
            note = ROOT.TNamed(
                "README",
                "Final systematic results. Use g_value_stat_*, g_value_syst_* and g_value_total_with_stat_* for preliminary plotting."
            )
            note.Write()
        except Exception:
            pass
    finally:
        f.Close()

    info(f"Wrote ROOT results with PyROOT: {root_path}")
    return True



def _cpp_quote(text: Any) -> str:
    s = str(text if text is not None else "")
    s = s.replace('\\', '\\\\').replace('"', '\\"').replace('\n', ' ')
    return '"' + s + '"'


def _cpp_float(x: Any) -> str:
    v = to_float(x, 0.0)
    if not isfinite(v):
        v = 0.0
    return f"{float(v):.17g}"


def _cpp_array(name: str, values: List[Any]) -> str:
    vals = ", ".join(_cpp_float(v) for v in values)
    return f"double {name}[] = {{{vals}}};"


def _root_cli_available(root_command: str) -> bool:
    import shutil
    return shutil.which(root_command) is not None or Path(root_command).exists()


def _write_final_root_file_root_cli(config: Dict[str, Any], rows: List[Dict[str, Any]], root_path: Path) -> bool:
    # Fallback for environments where the ROOT executable exists, but PyROOT is
    # not importable in the Python used to run this script.
    import subprocess

    root_command = str(config.get("root_command", "root"))
    if not _root_cli_available(root_command):
        warn(f"ROOT command-line fallback requested, but command is not available: {root_command}")
        return False

    project_dir = resolve_path(str(config.get("project_dir", ".")), Path.cwd()).resolve()
    macro_path = resolve_path(str(config.get("root_macro", "Systematics/generated_final_systematics_root_writer.C")), project_dir)
    macro_path.parent.mkdir(parents=True, exist_ok=True)
    root_path.parent.mkdir(parents=True, exist_ok=True)

    macro_func = re.sub(r"[^A-Za-z0-9_]", "_", macro_path.stem)
    if not macro_func or not re.match(r"^[A-Za-z_]", macro_func):
        macro_func = "generated_final_systematics_root_writer"

    helper = r'''
#include <vector>
#include <string>
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TNamed.h"
#include "TROOT.h"

void writeHist(const char* name, const char* title, double* edges, double* vals, double* errs, int n)
{
    TH1D* h = new TH1D(name, title, n, edges);
    h->SetDirectory(0);
    for (int i = 0; i < n; ++i) {
        h->SetBinContent(i + 1, vals[i]);
        if (errs) h->SetBinError(i + 1, errs[i]);
    }
    h->Write();
    delete h;
}

void writeGraph(const char* name, const char* title,
                double* x, double* exl, double* exh,
                double* y, double* eyl, double* eyh, int n)
{
    TGraphAsymmErrors* g = new TGraphAsymmErrors(n);
    g->SetName(name);
    g->SetTitle(title);
    for (int i = 0; i < n; ++i) {
        g->SetPoint(i, x[i], y[i]);
        g->SetPointError(i, exl[i], exh[i], eyl[i], eyh[i]);
    }
    g->Write();
    delete g;
}

void __MACRO_FUNC__()
{
    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);
'''.replace('__MACRO_FUNC__', macro_func)
    lines: List[str] = [helper]
    lines.append(f"    TFile* f = TFile::Open({_cpp_quote(str(root_path))}, \"RECREATE\");")
    lines.append('    if (!f || f->IsZombie()) { printf("Could not create ROOT file.\\n"); return; }')
    lines.append('    TNamed readme("README", "Final systematic results. Use g_value_stat_*, g_value_syst_* and g_value_total_with_stat_* for preliminary plotting.");')
    lines.append('    readme.Write();')

    groups = final_plot_groups(rows)
    block_id = 0
    for (job_label, cent), subrows0 in sorted(groups.items(), key=lambda kv: (kv[0][0], to_float(kv[0][1], 9999))):
        if not subrows0:
            continue
        subrows = sorted(subrows0, key=row_sort_bin)
        first = subrows[0]
        job_dir_name = safe_name(job_label)
        cent_label = str(first.get("centrality_label", cent))
        tag = safe_name(cent_label)
        title_base = f"{job_label}, {cent_label}"
        edges = _rows_to_edges(subrows)

        x_vals: List[float] = []
        exl_vals: List[float] = []
        exh_vals: List[float] = []
        for i, r in enumerate(subrows):
            lo = to_float(r.get("bin_low", ""), math.nan)
            hi = to_float(r.get("bin_high", ""), math.nan)
            if not isfinite(lo) or not isfinite(hi) or hi <= lo:
                lo = float(i)
                hi = float(i + 1)
            x = 0.5 * (lo + hi)
            x_vals.append(x)
            exl_vals.append(x - lo)
            exh_vals.append(hi - x)

        y = [to_float(r.get("nominal_value", ""), 0.0) for r in subrows]
        stat = [abs(to_float(r.get("stat_abs", r.get("nominal_stat_abs", "")), 0.0)) for r in subrows]
        syst_down_abs = [abs(to_float(r.get("total_syst_down_abs", ""), 0.0)) for r in subrows]
        syst_up_abs = [abs(to_float(r.get("total_syst_up_abs", ""), 0.0)) for r in subrows]
        syst_sym_abs = [abs(to_float(r.get("total_syst_sym_abs", ""), 0.0)) for r in subrows]
        syst_down_pct = [abs(to_float(r.get("total_syst_down_pct", ""), 0.0)) for r in subrows]
        syst_up_pct = [abs(to_float(r.get("total_syst_up_pct", ""), 0.0)) for r in subrows]
        syst_sym_pct = [abs(to_float(r.get("total_syst_sym_pct", ""), 0.0)) for r in subrows]
        total_down_abs = [abs(to_float(r.get("total_with_stat_down_abs", ""), 0.0)) for r in subrows]
        total_up_abs = [abs(to_float(r.get("total_with_stat_up_abs", ""), 0.0)) for r in subrows]
        total_sym_abs = [abs(to_float(r.get("total_with_stat_sym_abs", ""), 0.0)) for r in subrows]
        total_down_pct = [abs(to_float(r.get("total_with_stat_down_pct", ""), 0.0)) for r in subrows]
        total_up_pct = [abs(to_float(r.get("total_with_stat_up_pct", ""), 0.0)) for r in subrows]
        total_sym_pct = [abs(to_float(r.get("total_with_stat_sym_pct", ""), 0.0)) for r in subrows]

        arrays = {
            "edges": edges,
            "x": x_vals,
            "exl": exl_vals,
            "exh": exh_vals,
            "y": y,
            "stat": stat,
            "syst_down_abs": syst_down_abs,
            "syst_up_abs": syst_up_abs,
            "syst_sym_abs": syst_sym_abs,
            "syst_down_pct": syst_down_pct,
            "syst_up_pct": syst_up_pct,
            "syst_sym_pct": syst_sym_pct,
            "total_down_abs": total_down_abs,
            "total_up_abs": total_up_abs,
            "total_sym_abs": total_sym_abs,
            "total_down_pct": total_down_pct,
            "total_up_pct": total_up_pct,
            "total_sym_pct": total_sym_pct,
        }

        lines.append("    {")
        lines.append(f"        TDirectory* dir = f->GetDirectory({_cpp_quote(job_dir_name)});")
        lines.append(f"        if (!dir) dir = f->mkdir({_cpp_quote(job_dir_name)});")
        lines.append("        dir->cd();")
        n = len(subrows)
        for aname, vals in arrays.items():
            lines.append("        " + _cpp_array(f"{aname}_{block_id}", vals))
        lines.append(f"        int n_{block_id} = {n};")

        hist_specs = [
            (f"h_value_stat_{tag}", "y", "stat", "central values; stat. as TH1 errors"),
            (f"h_syst_down_abs_{tag}", "syst_down_abs", "0", "absolute systematic down uncertainty"),
            (f"h_syst_up_abs_{tag}", "syst_up_abs", "0", "absolute systematic up uncertainty"),
            (f"h_syst_sym_abs_{tag}", "syst_sym_abs", "0", "absolute symmetric systematic uncertainty"),
            (f"h_syst_down_pct_{tag}", "syst_down_pct", "0", "percent systematic down uncertainty"),
            (f"h_syst_up_pct_{tag}", "syst_up_pct", "0", "percent systematic up uncertainty"),
            (f"h_syst_sym_pct_{tag}", "syst_sym_pct", "0", "percent symmetric systematic uncertainty"),
            (f"h_total_with_stat_down_abs_{tag}", "total_down_abs", "0", "absolute stat+syst down uncertainty"),
            (f"h_total_with_stat_up_abs_{tag}", "total_up_abs", "0", "absolute stat+syst up uncertainty"),
            (f"h_total_with_stat_sym_abs_{tag}", "total_sym_abs", "0", "absolute symmetric stat+syst uncertainty"),
            (f"h_total_with_stat_down_pct_{tag}", "total_down_pct", "0", "percent stat+syst down uncertainty"),
            (f"h_total_with_stat_up_pct_{tag}", "total_up_pct", "0", "percent stat+syst up uncertainty"),
            (f"h_total_with_stat_sym_pct_{tag}", "total_sym_pct", "0", "percent symmetric stat+syst uncertainty"),
        ]
        for hname, vals_name, errs_name, desc in hist_specs:
            errs_expr = "0" if errs_name == "0" else f"{errs_name}_{block_id}"
            lines.append(
                f"        writeHist({_cpp_quote(hname)}, {_cpp_quote(title_base + ';bin;' + desc)}, "
                f"edges_{block_id}, {vals_name}_{block_id}, {errs_expr}, n_{block_id});"
            )

        lines.append(
            f"        writeGraph({_cpp_quote('g_value_stat_' + tag)}, {_cpp_quote(title_base + ';x;value')}, "
            f"x_{block_id}, exl_{block_id}, exh_{block_id}, y_{block_id}, stat_{block_id}, stat_{block_id}, n_{block_id});"
        )
        lines.append(
            f"        writeGraph({_cpp_quote('g_value_syst_' + tag)}, {_cpp_quote(title_base + ';x;value')}, "
            f"x_{block_id}, exl_{block_id}, exh_{block_id}, y_{block_id}, syst_down_abs_{block_id}, syst_up_abs_{block_id}, n_{block_id});"
        )
        lines.append(
            f"        writeGraph({_cpp_quote('g_value_total_with_stat_' + tag)}, {_cpp_quote(title_base + ';x;value')}, "
            f"x_{block_id}, exl_{block_id}, exh_{block_id}, y_{block_id}, total_down_abs_{block_id}, total_up_abs_{block_id}, n_{block_id});"
        )
        lines.append(f"        TNamed info_{block_id}({_cpp_quote('info_' + tag)}, {_cpp_quote(f'job_label={job_label}; centrality_label={cent_label}; tag={tag}')});")
        lines.append(f"        info_{block_id}.Write();")
        lines.append("        f->cd();")
        lines.append("    }")
        block_id += 1

    lines.append('    f->Close();')
    lines.append('    delete f;')
    lines.append('}')
    # The generated function name matches the macro filename, so ROOT executes
    # it with the standard form: root -l -b -q generated_macro.C.
    # Do not append a free-standing function call here; Cling treats that as
    # invalid C++ in this macro context.
    lines.append('')

    macro_path.write_text("\n".join(lines), encoding="utf-8")
    cmd = [root_command, "-l", "-b", "-q", str(macro_path)]
    info("Running ROOT fallback: " + " ".join(cmd))
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        warn("ROOT command-line fallback failed. Last output lines:\n" + "\n".join(proc.stdout.splitlines()[-30:]))
        return False
    info(f"Wrote ROOT results with ROOT command-line fallback: {root_path}")
    info(f"ROOT macro kept for inspection: {macro_path}")
    return True


def _file_nonempty(path: Path) -> bool:
    try:
        return path.exists() and path.is_file() and path.stat().st_size > 0
    except Exception:
        return False


def write_final_root_file(config: Dict[str, Any], rows: List[Dict[str, Any]]) -> Optional[Path]:
    """Write the final ROOT file and fail loudly if this is impossible.

    v7 deliberately avoids using an old, already-existing ROOT file as a
    successful output.  The ROOT writers first create a temporary file and only
    after a successful non-empty write is it atomically moved to output_root.
    This prevents the downstream preliminary plotter from silently reading stale
    results when PyROOT/ROOT fails.
    """
    enabled = cfg_bool(config.get("write_root", True), True)
    if not enabled:
        info("ROOT output disabled by config write_root=false.")
        return None

    if not rows:
        raise RuntimeError("Cannot write ROOT output because the final table has zero rows.")

    project_dir = resolve_path(str(config.get("project_dir", ".")), Path.cwd()).resolve()
    root_path = resolve_path(str(config.get("output_root", "Systematics/final_systematics_results_v9.root")), project_dir)
    root_path.parent.mkdir(parents=True, exist_ok=True)

    tmp_root = root_path.with_name(f".{root_path.name}.tmp.{os.getpid()}")
    if tmp_root.exists():
        tmp_root.unlink()

    ok = False
    try:
        ok = _write_final_root_file_pyroot(config, rows, tmp_root)
        if not ok:
            ok = _write_final_root_file_root_cli(config, rows, tmp_root)

        if not ok or not _file_nonempty(tmp_root):
            raise RuntimeError(
                "Could not create a fresh final ROOT output. "
                "PyROOT failed and/or the ROOT command-line fallback failed. "
                "Run inside a ROOT/PyROOT environment or set root_command in the config."
            )

        os.replace(str(tmp_root), str(root_path))
        if not _file_nonempty(root_path):
            raise RuntimeError(f"Final ROOT file was moved into place but is empty: {root_path}")

        info(f"Fresh ROOT output is ready: {root_path}")
        return root_path
    finally:
        if tmp_root.exists():
            try:
                tmp_root.unlink()
            except Exception:
                pass

def write_outputs(config: Dict[str, Any], rows: List[Dict[str, Any]], detail_rows: List[Dict[str, Any]], group_names: List[str]) -> None:
    project_dir = resolve_path(str(config.get("project_dir", ".")), Path.cwd()).resolve()
    output_tsv = resolve_path(str(config.get("output_tsv", "Systematics/final_systematics_v9.tsv")), project_dir)
    output_csv = resolve_path(str(config.get("output_csv", "Systematics/final_systematics_v9.csv")), project_dir)

    columns = build_output_columns(rows, group_names)
    write_table(output_tsv, rows, columns, delimiter="\t")
    info(f"Wrote TSV: {output_tsv}")
    write_table(output_csv, rows, columns, delimiter=",")
    info(f"Wrote CSV: {output_csv}")

    if bool(config.get("write_details", True)):
        details_path = resolve_path(str(config.get("details_tsv", "Systematics/final_systematics_barlow_details_v9.tsv")), project_dir)
        write_table(details_path, detail_rows, DETAIL_COLUMNS, delimiter="\t")
        info(f"Wrote details TSV: {details_path}")

    root_path = write_final_root_file(config, rows)
    if cfg_bool(config.get("write_root", True), True) and root_path is None:
        raise RuntimeError("ROOT output was requested but no ROOT path was produced.")
    write_final_plots(config, rows)


# -----------------------------------------------------------------------------
# Template config
# -----------------------------------------------------------------------------

def default_template_config() -> Dict[str, Any]:
    return {
        "project_dir": ".",
        "input": "Systematics/systematics_components_v15.tsv",
        "output_tsv": "Systematics/final_systematics_v9.tsv",
        "output_csv": "Systematics/final_systematics_v9.csv",
        "details_tsv": "Systematics/final_systematics_barlow_details_v9.tsv",
        "output_root": "Systematics/final_systematics_results_v9.root",
        "write_root": True,
        "write_details": True,
        "barlow": {
            "enabled": True,
            "mode": "correlated",
            "threshold": 1.0,
            "action": "subtract_stat_diff",
            "missing_stat_action": "keep",
            "denominator_floor": 1e-30,
            "apply_to_component_types": ["variation"],
            "default_apply_to_groups": True,
            "comment": "Use subtract_stat_diff for the presentation-style Barlow correction: 0 if |Delta| <= sigma_stat,diff, else sqrt(Delta^2 - sigma_stat,diff^2). Use flag_only for diagnostics only."
        },
        "groups": [
            {
                "name": "iteration",
                "select": {"source_group": "Nominal", "source_name": "iteration"},
                "combine_variations": "envelope",
                "apply_barlow": False
            },
            {
                "name": "binning_first_variable",
                "select": {"source_group": "Binning", "source_name": "binning_first_variable"},
                "combine_variations": "take_rms",
                "apply_barlow": False
            },
            {
                "name": "binning_second_variable",
                "select": {"source_group": "Binning", "source_name": "binning_second_variable"},
                "combine_variations": "take_rms",
                "apply_barlow": False
            },
            {
                "name": "prior_shape_first_variable",
                "select": {"source_group": "PriorShape", "source_name": "prior_shape_first_variable"},
                "combine_variations": "envelope"
            },
            {
                "name": "prior_shape_second_variable",
                "select": {"source_group": "PriorShape", "source_name": "prior_shape_second_variable"},
                "combine_variations": "envelope"
            },
            {
                "name": "sweight",
                "select": {"source_group": "sWeight"},
                "combine_variations": "max_abs"
            },
            {
                "name": "jetsreco",
                "select": {"source_group": "JetsReco"},
                "combine_variations": "quadrature"
            },
            {
                "name": "d0meson",
                "select": {"source_group": "D0Meson"},
                "combine_variations": "pair_envelope_then_quadrature",
                "pair_by": "source_name"
            }
        ],
        "total": {
            "combine_groups": "mixed",
            "default_group_mode": "quadrature",
            "include_stat": True,
            "include_groups": [
                "iteration",
                "binning_first_variable",
                "binning_second_variable",
                "prior_shape_first_variable",
                "prior_shape_second_variable",
                "sweight",
                "jetsreco",
                "d0meson"
            ]
        },
        "plots": {
            "enabled": True,
            "output_dir": "Systematics/FinalPlots",
            "pdf": "final_systematics_values_v9.pdf",
            "write_png": True,
            "figsize": [7.0, 5.0],
            "dpi": 150,
            "draw_syst_boxes": True,
            "draw_stat_errors": True,
            "draw_points": True,
            "draw_unity_for_rcp": True,
            "grid": True,
            "legend": True,
            "title_template": "{job_label}, {centrality_label}",
            "stat_label": "stat.",
            "syst_label": "syst.",
            "logy": False,
            "logy_spectra": True,
            "force_default_logy_spectra": True,
            "logy_rcp": False,
            "ylim": [],
            "spectra_ylim": [],
            "rcp_ylim": [0.0, 2.0],
            "force_default_rcp_ylim": True,
            "ylim_by_job": {},
            "ylim_by_observable": {},
            "ylim_by_plot": {},
            "rcp_ylim_by_job": {},
            "x_titles": {},
            "y_titles": {}
        }
    }


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Combine long systematic components into final systematic table.")
    p.add_argument("--config", default="", help="JSON config file.")
    p.add_argument("--write-template", default="", help="Write template JSON config and exit.")
    p.add_argument("--input", default="", help="Override input TSV from config.")
    p.add_argument("--output", default="", help="Override output TSV from config.")
    p.add_argument(
        "--barlow-action",
        default="",
        choices=["", "flag_only", "zero_if_not_significant", "subtract_stat_diff", "barlow_subtract", "presentation"],
        help="Override barlow.action.",
    )
    p.add_argument("--barlow-mode", default="", choices=["", "correlated", "independent"], help="Override barlow.mode.")
    p.add_argument("--barlow-threshold", default="", help="Override barlow.threshold.")
    p.add_argument("--no-plots", action="store_true", help="Disable final plotting even if enabled in config.")
    p.add_argument("--output-root", default="", help="Override output ROOT file from config.")
    p.add_argument("--force-root", action="store_true", help="Force write_root=true even if the config disables ROOT output.")
    p.add_argument("--no-root", action="store_true", help="Disable ROOT output.")
    p.add_argument("--root-macro", default="", help="Override temporary ROOT CLI helper macro path.")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if args.write_template:
        out = Path(args.write_template)
        write_json(out, default_template_config())
        info(f"Template config written to: {out}")
        return 0

    if args.config:
        config = read_json(Path(args.config))
    else:
        config = default_template_config()
        warn("No --config specified; using built-in default config.")

    if args.input:
        config["input"] = args.input
    if args.output:
        config["output_tsv"] = args.output
    if args.barlow_action:
        config.setdefault("barlow", {})["action"] = args.barlow_action
    if args.barlow_mode:
        config.setdefault("barlow", {})["mode"] = args.barlow_mode
    if args.barlow_threshold:
        config.setdefault("barlow", {})["threshold"] = float(args.barlow_threshold)
    if args.no_plots:
        config.setdefault("plots", {})["enabled"] = False
    if args.output_root:
        config["output_root"] = args.output_root
    if args.force_root:
        config["write_root"] = True
    if args.no_root:
        config["write_root"] = False
    if args.root_macro:
        config["root_macro"] = args.root_macro

    rows, detail_rows, group_names = build_final_table(config)
    if not rows:
        raise RuntimeError("Final table has zero rows; refusing to write empty outputs.")
    write_outputs(config, rows, detail_rows, group_names)
    info(f"Final rows: {len(rows)}")
    info(f"Detail rows: {len(detail_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
