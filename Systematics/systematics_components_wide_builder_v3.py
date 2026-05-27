#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
systematics_components_wide_builder_v3.py

Memory-light converter with compact prior-shape column names from the long systematics component table to a wide
inspection table.

Input:
  Systematics/systematics_components_v15.tsv

Output:
  one row = one histogram / centrality(or R_CP index) / bin
  columns = nominal values, summary max-abs columns, individual component pct columns

This version intentionally avoids pandas pivoting.  It streams the TSV with the
standard csv module and stores only one dictionary per output bin, which is tiny
compared to the ROOT jobs/scans.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple


BASE_COLS = [
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
]

KEY_COLS = [
    "job_label",
    "hist_name",
    "centrality",
    "method",
    "iteration_display",
    "bin",
    "bin_low",
    "bin_high",
]

# Nice order for the most important group-level summaries.
SUMMARY_GROUP_ORDER = [
    "iteration_max_abs_pct",
    "binning_max_abs_pct",
    "priorshape_max_abs_pct",
    "sweight_max_abs_pct",
    "jetsreco_max_abs_pct",
    "d0meson_max_abs_pct",
]


def info(msg: str) -> None:
    print(f"[info] {msg}")


def warn(msg: str) -> None:
    print(f"[warning] {msg}")


def snake(text: Any) -> str:
    """Make a stable, readable column-name fragment."""
    s = str(text or "").strip()
    if not s:
        return "unknown"
    s = s.replace("+", "plus")
    s = s.replace("-", "minus")
    s = s.replace("%", "pct")
    s = s.replace("/", "_over_")
    s = s.replace("#", "")
    s = re.sub(r"\^\{?([^}\s]+)\}?", r"_pow_\1", s)
    s = re.sub(r"_\{?([^}\s]+)\}?", r"_\1", s)
    s = re.sub(r"[^0-9A-Za-z]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_").lower()
    return s or "unknown"


def group_prefix(source_group: str) -> str:
    sg = str(source_group or "").strip().lower()
    if sg == "priorshape":
        return "priorshape"
    if sg == "sweight":
        return "sweight"
    if sg == "jetsreco":
        return "jetsreco"
    if sg == "d0meson":
        return "d0meson"
    if sg == "binning":
        return "binning"
    return snake(source_group)


def compact_source_name(source_group: str, source_name: str) -> str:
    """Source-name fragment used for summary columns."""
    gp = group_prefix(source_group)
    s_sn = snake(source_name)
    if gp == "priorshape":
        if s_sn in {"prior_shape_first_variable", "first_variable"}:
            return "first_variable"
        if s_sn in {"prior_shape_second_variable", "second_variable"}:
            return "second_variable"
    return s_sn


def compact_priorshape_column(source_name: str, variation_name: str) -> str:
    """Collapse variable-specific prior-shape names to first/second-variable columns."""
    s_sn = snake(source_name)
    s_vn = snake(variation_name)

    if s_sn in {"prior_shape_first_variable", "first_variable"}:
        if "plus20" in s_vn:
            return "first_variable_plus20_pct"
        if "minus20" in s_vn:
            return "first_variable_minus20_pct"
        return f"first_variable_{s_vn}_pct"

    if s_sn in {"prior_shape_second_variable", "second_variable"}:
        if "plus20" in s_vn:
            return "second_variable_plus20_pct"
        if "minus20" in s_vn:
            return "second_variable_minus20_pct"
        return f"second_variable_{s_vn}_pct"

    return ""


def parse_float(value: Any) -> float:
    try:
        if value is None:
            return math.nan
        s = str(value).strip()
        if not s or s.lower() in {"nan", "none", "<na>"}:
            return math.nan
        return float(s)
    except Exception:
        return math.nan


def finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except Exception:
        return False


def fmt_float(value: float) -> str:
    if not finite(value):
        return ""
    return f"{float(value):.10g}"


def get_pct_value(row: Dict[str, str]) -> Tuple[str, float]:
    """Return printable signed pct value and numeric abs value for summaries."""
    diff_pct = parse_float(row.get("diff_pct", ""))
    abs_pct = parse_float(row.get("abs_pct", ""))

    # Prefer signed diff_pct for individual columns.  For aggregate RMS rows,
    # diff_pct is already positive/symmetric.  If missing, fall back to abs_pct.
    if finite(diff_pct):
        signed = diff_pct
    elif finite(abs_pct):
        signed = abs_pct
    else:
        signed = math.nan

    if finite(abs_pct):
        absval = abs_pct
    elif finite(diff_pct):
        absval = abs(diff_pct)
    else:
        absval = math.nan

    return fmt_float(signed), absval


def component_column(row: Dict[str, str]) -> str:
    """Column name for the individual component value."""
    sg = row.get("source_group", "")
    sn = row.get("source_name", "")
    vn = row.get("variation_name", "")
    gp = group_prefix(sg)
    s_sn = snake(sn)
    s_vn = snake(vn)

    if gp == "nominal" and s_sn == "iteration":
        return f"{s_vn}_pct"

    # Prior-shape second-variable variations are variable-specific in the long
    # table (z_plus20, lambda11_plus20, ...), but for each output row only the
    # relevant observable exists.  Collapse them to common inspection columns.
    if gp == "priorshape":
        compact = compact_priorshape_column(sn, vn)
        if compact:
            return compact

    # Avoid duplicated names like sweight_cheby2_bkg_cheby2_bkg_pct.
    if s_vn == s_sn or s_vn.startswith(s_sn + "_"):
        return f"{gp}_{s_vn}_pct"
    return f"{gp}_{s_sn}_{s_vn}_pct"


def component_n_column(pct_col: str) -> str:
    if pct_col.endswith("_pct"):
        return pct_col[:-4] + "_n"
    return pct_col + "_n"


def update_max(row_out: Dict[str, str], col: str, value: float) -> None:
    if not finite(value):
        return
    old = parse_float(row_out.get(col, ""))
    if not finite(old) or value > old:
        row_out[col] = fmt_float(value)


def should_skip_component(row: Dict[str, str]) -> bool:
    sg = str(row.get("source_group", "")).strip()
    sn = str(row.get("source_name", "")).strip()
    vn = str(row.get("variation_name", "")).strip()

    # Base nominal/statistical values are already represented by BASE_COLS.
    if sg == "Nominal" and sn in {"nominal_value", "statistical_uncertainty"}:
        return True
    # Empty/unknown source rows are not useful as wide components.
    if not sg or not sn or not vn:
        return True
    return False


def row_key(row: Dict[str, str]) -> Tuple[str, ...]:
    return tuple(str(row.get(c, "")) for c in KEY_COLS)


def ensure_base_row(rows: "OrderedDict[Tuple[str, ...], Dict[str, str]]", row: Dict[str, str]) -> Dict[str, str]:
    key = row_key(row)
    if key not in rows:
        out = OrderedDict()
        for c in BASE_COLS:
            out[c] = row.get(c, "")
        rows[key] = out
    else:
        out = rows[key]
        # Fill missing base metadata if a previous row had blanks.
        for c in BASE_COLS:
            if not out.get(c) and row.get(c):
                out[c] = row.get(c, "")
    return out


def build_wide(
    input_path: Path,
    include_abs_columns: bool = False,
    include_n_columns: bool = True,
) -> Tuple[List[Dict[str, str]], List[str]]:
    rows: "OrderedDict[Tuple[str, ...], Dict[str, str]]" = OrderedDict()
    component_cols: "OrderedDict[str, None]" = OrderedDict()
    abs_cols: "OrderedDict[str, None]" = OrderedDict()
    n_cols: "OrderedDict[str, None]" = OrderedDict()
    summary_cols: "OrderedDict[str, None]" = OrderedDict()

    n_read = 0
    n_components = 0

    with input_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = [c for c in BASE_COLS + ["source_group", "source_name", "variation_name", "diff_pct", "abs_pct"] if c not in (reader.fieldnames or [])]
        if missing:
            raise RuntimeError(f"Input is missing required columns: {missing}")

        for row in reader:
            n_read += 1
            out = ensure_base_row(rows, row)

            if should_skip_component(row):
                continue

            pct_col = component_column(row)
            pct_text, abs_val = get_pct_value(row)
            out[pct_col] = pct_text
            component_cols.setdefault(pct_col, None)
            n_components += 1

            if include_abs_columns:
                abs_col = pct_col[:-4] + "_abs_pct" if pct_col.endswith("_pct") else pct_col + "_abs_pct"
                out[abs_col] = fmt_float(abs_val)
                abs_cols.setdefault(abs_col, None)

            if include_n_columns:
                nvar = str(row.get("n_variations", "")).strip()
                if nvar:
                    n_col = component_n_column(pct_col)
                    out[n_col] = nvar
                    n_cols.setdefault(n_col, None)

            sg = row.get("source_group", "")
            sn = row.get("source_name", "")
            gp = group_prefix(sg)
            s_sn = compact_source_name(sg, sn)

            # Summary columns.
            if gp == "nominal" and s_sn == "iteration":
                scol = "iteration_max_abs_pct"
                update_max(out, scol, abs_val)
                summary_cols.setdefault(scol, None)
            elif gp != "nominal":
                gcol = f"{gp}_max_abs_pct"
                scol = f"{gp}_{s_sn}_max_abs_pct"
                update_max(out, gcol, abs_val)
                update_max(out, scol, abs_val)
                summary_cols.setdefault(gcol, None)
                summary_cols.setdefault(scol, None)

    info(f"Read rows: {n_read} from {input_path}")
    info(f"Output bins: {len(rows)}")
    info(f"Component entries used: {n_components}")
    info(f"Individual pct columns: {len(component_cols)}")
    info(f"Summary columns: {len(summary_cols)}")

    # Put the most common summary columns first, then source-specific summaries.
    summary_order: List[str] = []
    for c in SUMMARY_GROUP_ORDER:
        if c in summary_cols:
            summary_order.append(c)
    for c in summary_cols.keys():
        if c not in summary_order:
            summary_order.append(c)

    output_cols = BASE_COLS + summary_order + list(component_cols.keys()) + list(abs_cols.keys()) + list(n_cols.keys())

    return list(rows.values()), output_cols


def write_table(path: Path, rows: List[Dict[str, str]], columns: List[str], delimiter: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in columns})
    info(f"Wrote: {path}")


def default_outputs(input_path: Path) -> Tuple[Path, Path]:
    base_dir = input_path.parent
    return (
        base_dir / "systematics_components_wide_v3.tsv",
        base_dir / "systematics_components_wide_v3.csv",
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Memory-light wide table builder for systematics components.")
    p.add_argument("--input", required=True, help="Input long TSV, e.g. Systematics/systematics_components_v15.tsv")
    p.add_argument("--output-tsv", default="", help="Output wide TSV. Default: next to input.")
    p.add_argument("--output-csv", default="", help="Output wide CSV. Default: next to input.")
    p.add_argument("--include-abs-columns", action="store_true", help="Also write one *_abs_pct column per individual component.")
    p.add_argument("--no-n-columns", action="store_true", help="Do not write *_n columns for aggregate components such as binning RMS.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    default_tsv, default_csv = default_outputs(input_path)
    output_tsv = Path(args.output_tsv) if args.output_tsv else default_tsv
    output_csv = Path(args.output_csv) if args.output_csv else default_csv

    rows, columns = build_wide(
        input_path=input_path,
        include_abs_columns=bool(args.include_abs_columns),
        include_n_columns=not bool(args.no_n_columns),
    )

    write_table(output_tsv, rows, columns, delimiter="\t")
    write_table(output_csv, rows, columns, delimiter=",")
    info(f"Wide rows: {len(rows)}")
    info(f"Wide columns: {len(columns)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
