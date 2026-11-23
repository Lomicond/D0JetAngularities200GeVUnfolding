#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Effective radius from the SAME final ROOT file as
Systematics/finalPlot_from_systematics_root_v17_lido_spectra.py.

Stage 3: central values + statistical AND systematic uncertainties.

By default, both uncertainties use the DIAGONAL approximation. Systematics
are read from final g_value_syst graphs, i.e. the final combined per-bin
uncertainties. Their unknown bin/angularity correlations can materially affect
uncertainties on means and ratios. No correlations are invented.

Optional covariance NPZ inputs can be provided separately for stat and syst.
For systematic uncertainties, an even better eventual method is to propagate
individual coherent systematic variations through means and radius ratios.

By DEFAULT runs BOTH D0 pT ranges (1-10 and 5-10 GeV/c) from the
Systematics directory, writing a combined comparison table and plots.

Run from the unfolding project root, for example:
  python3 Systematics/effective_radius_from_final_root_v5.py
  python3 Systematics/effective_radius_from_final_root_v5.py --first-bin-zero
  python3 Systematics/effective_radius_from_final_root_v5.py --bin-policy plotter

Compare two separate final-systematics files:
  python3 Systematics/effective_radius_from_final_root_v5.py \
    --root /path/to/one_to_ten.root --label '1-10 GeV/c' \
    --root /path/to/five_to_ten.root --label '5-10 GeV/c' \
    --first-bin-zero

Input graph example: l051/g_value_stat_0_10pct
Output: three TSV tables and two PDF plots in Systematics/EffectiveRadius/.

The default statistical and systematic uncertainties use independent bins and
independent angularities (DIAGONAL approximations). These assumptions are NOT
exact. Optional --stat-cov-npz and --syst-cov-npz inputs can supply covariance
matrices for each spectrum and/or the five mean angularities. NPZ key format:
  sample0__0_10pct__lha   # n_bins x n_bins covariance in density units
  sample0__0_10pct__means # 5x5 mean covariance, order lha,girth,gti,thrust,jts
(sample0 = 1-10 and sample1 = 5-10 for the default ROOT files.)
If a mean covariance is provided, it replaces the corresponding per-angularity
errors for derived observables. These are separately supplied for stat and syst.

Physics:
  mean(lambda_a) = sum_i (bin_center_i * density_i * bin_width_i)
                   / sum_i (density_i * bin_width_i)
  R_eff(a,b) = R_jet * (mean(lambda_a) / mean(lambda_b)) ** (1/(a-b))

IMPORTANT: setting the FIRST bin to exactly zero sets ONLY its representative
x value to 0; the bin's yield still contributes to the normalization.

If spectra are not normalized over the same underlying jet sample because of
truncation, different bin ranges, or zero-bin approximations, cancellation of
single-constituent jets in ratios of means is NOT mathematically exact.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

try:
    import uproot
except ImportError:
    uproot = None

# ===================== USER SETTINGS =====================
HERE = Path(__file__).resolve().parent
PROJECT_DIR = HERE.parent
# Default: process both D0 pT selections automatically in one run.
# Pass --root (repeatable) to override this list for any custom inputs.
DEFAULT_SAMPLES = (
    (HERE / "final_systematics_results_v10.root", "1-10 GeV/c"),
    (HERE / "final_systematics_results_v10_5_10.root", "5-10 GeV/c"),
)
DEFAULT_OUTPUT = HERE / "EffectiveRadius"
R_JET_DEFAULT = 0.4
FIRST_BIN_IS_EXACT_ZERO_DEFAULT = False
BIN_POLICY_DEFAULT = "full"   # "full" or "plotter"
NEGATIVE_BINS_DEFAULT = "keep"  # "keep", "clip", "error"

CENTRALITY_TAGS = ("0_10pct", "10_40pct", "40_80pct")
CENTRALITY_LABELS = ("0-10%", "10-40%", "40-80%")
CENTRALITY_MIDPOINTS = (5.0, 25.0, 60.0)

# Same ROOT directory names as the existing finalPlot plotter.
ANGULARITIES = {
    "lha": (0.5, "l051"),
    "girth": (1.0, "l11"),
    "gti": (1.5, "l15"),
    "thrust": (2.0, "l21"),
    "jts": (3.0, "l31"),
}

ANGULARITY_LABELS = {
    "lha": "LHA",
    "girth": "Girth",
    "gti": "GTI",
    "thrust": "Thrust",
    "jts": "JTS",
}

# The existing plotter SHOW_EDGE_BINS switch (angularities have last=False).
PLOTTER_EDGE_POLICY = {
    "lha": {"first": True, "last": False},
    "girth": {"first": True, "last": False},
    "gti": {"first": True, "last": False},
    "thrust": {"first": True, "last": False},
    "jts": {"first": True, "last": False},
}

def build_radius_pairs() -> Tuple[Tuple[str, str, str], ...]:
    """
    Build all mathematically reasonable effective-radius combinations,
    i.e. all unique pairs with alpha_num > alpha_den.
    """
    items = sorted(
        ((name, cfg[0]) for name, cfg in ANGULARITIES.items()),
        key=lambda x: x[1]
    )

    pairs = []
    for i in range(len(items)):
        for j in range(i):
            num_name, num_alpha = items[i]
            den_name, den_alpha = items[j]
            key = f"{num_name}_over_{den_name}"
            pairs.append((key, num_name, den_name))
    return tuple(pairs)


RADIUS_PAIRS = build_radius_pairs()
# =========================================================


def resolve_path(raw: str) -> Path:
    """Allow paths relative to cwd, project root, or Systematics directory."""
    p = Path(raw).expanduser()
    if p.is_absolute():
        return p
    for q in (Path.cwd() / p, PROJECT_DIR / p, HERE / p):
        if q.exists():
            return q.resolve()
    return (Path.cwd() / p).resolve()


def member_array(obj: Any, name: str, n: Optional[int] = None) -> np.ndarray:
    try:
        data = np.asarray(obj.member(name), dtype=float)
    except Exception:
        data = np.zeros(0, dtype=float)
    if n is not None:
        if len(data) == 0:
            return np.zeros(n, dtype=float)
        if len(data) < n:
            raise ValueError(f"ROOT member {name}: {len(data)} points, expected {n}")
        data = data[:n]
    return data


def read_spectrum(root: Any, directory: str, tag: str) -> Dict[str, np.ndarray]:
    """Read final differential yields from stat TGraphAsymmErrors or TH1."""
    key = f"{directory}/g_value_stat_{tag}"
    try:
        graph = root[key]
    except KeyError:
        graph = None

    if graph is not None:
        x = member_array(graph, "fX")
        y = member_array(graph, "fY")
        n = min(len(x), len(y))
        if n == 0:
            raise ValueError(f"Empty graph {key}")
        x, y = x[:n], y[:n]
        exl = member_array(graph, "fEXlow", n)
        exh = member_array(graph, "fEXhigh", n)
        # Optional fallback if the graph is TGraphErrors (not asymmetric).
        if not (np.any(exl > 0) or np.any(exh > 0)):
            ex = member_array(graph, "fEX", n)
            exl, exh = ex.copy(), ex.copy()
        stat_low = member_array(graph, "fEYlow", n)
        stat_high = member_array(graph, "fEYhigh", n)
        if not (np.any(stat_low > 0) or np.any(stat_high > 0)):
            ey = member_array(graph, "fEY", n)
            stat_low, stat_high = ey.copy(), ey.copy()
        low, high = x - exl, x + exh
        if np.any(high <= low):
            raise ValueError(f"{key}: missing/invalid horizontal ROOT errors, cannot reconstruct bin widths")
    else:
        key = f"{directory}/h_value_stat_{tag}"
        try:
            h = root[key]
        except KeyError as exc:
            raise KeyError(f"Neither {directory}/g_value_stat_{tag} nor {key} found") from exc
        y, edges = h.to_numpy(flow=False)
        y = np.asarray(y, dtype=float)
        edges = np.asarray(edges, dtype=float)
        low, high = edges[:-1], edges[1:]
        x = (low + high) / 2.0
        try:
            var = h.variances(flow=False)
            stat_low = np.sqrt(np.maximum(var, 0)) if var is not None else np.zeros_like(y)
        except Exception:
            stat_low = np.zeros_like(y)
        stat_high = stat_low.copy()

    order = np.argsort(low)
    low, high, x, y = (np.asarray(v, dtype=float)[order] for v in (low, high, x, y))
    stat_low = np.asarray(stat_low, dtype=float)[order]
    stat_high = np.asarray(stat_high, dtype=float)[order]
    if not np.all(np.isfinite(low) & np.isfinite(high) & np.isfinite(x) & np.isfinite(y)):
        raise ValueError(f"{key}: nonfinite x or y")
    if np.any(high <= low) or np.any(low[1:] < high[:-1] - 1.0e-9):
        raise ValueError(f"{key}: invalid/overlapping bin boundaries")

    # Final combined per-bin systematic errors (same objects as preliminary plotting).
    syst_low = np.full(len(y), np.nan)
    syst_high = np.full(len(y), np.nan)
    try:
        gsys = root[f"{directory}/g_value_syst_{tag}"]
        sx = member_array(gsys, "fX")
        sy = member_array(gsys, "fY")
        sl = member_array(gsys, "fEYlow", len(sx))
        sh = member_array(gsys, "fEYhigh", len(sx))
        if len(sx) == len(x) and np.allclose(np.sort(sx), x, rtol=1e-6, atol=1e-8):
            sort_sys = np.argsort(sx)
            syst_low, syst_high = sl[sort_sys], sh[sort_sys]
        else:
            print(f"[WARN] {directory}/{tag}: systematic graph bins differ; uncertainties not copied", file=sys.stderr)
    except KeyError:
        pass

    # Consistent with finalPlot_from_systematics_root: TH1 fallback.
    # The histogram content is the ABSOLUTE uncertainty in density units.
    if not np.all(np.isfinite(syst_low) & np.isfinite(syst_high)):
        try:
            hd = root[f"{directory}/h_syst_down_abs_{tag}"]
            hu = root[f"{directory}/h_syst_up_abs_{tag}"]
            vd, ed = hd.to_numpy(flow=False)
            vu, eu = hu.to_numpy(flow=False)
            if (len(vd) == len(y) and len(vu) == len(y)
                    and np.allclose(ed[:-1], low, rtol=0.0, atol=1e-8)
                    and np.allclose(ed[1:], high, rtol=0.0, atol=1e-8)
                    and np.allclose(eu, ed, rtol=0.0, atol=1e-8)):
                syst_low = np.abs(np.asarray(vd, dtype=float))
                syst_high = np.abs(np.asarray(vu, dtype=float))
            else:
                print(f"[WARN] {directory}/{tag}: syst TH1 bins differ from nominal", file=sys.stderr)
        except KeyError:
            pass
    if not np.all(np.isfinite(syst_low) & np.isfinite(syst_high)):
        print(f"[WARN] {directory}/{tag}: systematic uncertainties unavailable", file=sys.stderr)

    return {
        "low": low, "high": high, "center": x, "density": y,
        "stat_low": stat_low, "stat_high": stat_high,
        "syst_low": syst_low, "syst_high": syst_high,
    }


def load_covariance(npz: Any, key: str, dim: int) -> Optional[np.ndarray]:
    """Optional covariance (statistical OR systematic) in density units (means for __means).

    NPZ keys for each sample (0=1-10, 1=5-10) and centrality:
      sample0__0_10pct__lha          # per-bin covariance, n_bins x n_bins
      sample0__0_10pct__means        # covariance of 5 means, 5x5
    Same for other samples, centralities, and observables. Order of the 5
    mean variables is ANGULARITIES insertion order: lha,girth,gti,thrust,jts.
    All covariance arrays must include exactly the same bins as the ROOT input.
    """
    if npz is None or key not in npz:
        return None
    cov = np.asarray(npz[key], dtype=float)
    if cov.shape != (dim, dim) or not np.all(np.isfinite(cov)):
        raise ValueError(f"{key}: expected finite {dim}x{dim} covariance, got {cov.shape}")
    if not np.allclose(cov, cov.T, atol=1.e-10, rtol=1.e-7):
        raise ValueError(f"{key}: covariance matrix is not symmetric")
    cov = (cov + cov.T) / 2.0
    eigen = np.linalg.eigvalsh(cov)
    if eigen[0] < -1.e-7 * max(1.e-20, float(np.max(np.abs(eigen)))):
        raise ValueError(f"{key}: covariance matrix has a significant negative eigenvalue")
    return cov


def signed_diagonal_error(grad: np.ndarray,
                          stat_low: np.ndarray,
                          stat_high: np.ndarray) -> Tuple[float, float]:
    """First-order downward/upward error for independent asymmetric errors."""
    # Increasing a bin with a negative derivative decreases the observable.
    err_up = np.where(grad >= 0.0, stat_high, stat_low)
    err_down = np.where(grad >= 0.0, stat_low, stat_high)
    if not np.all(np.isfinite(err_up) & np.isfinite(err_down)):
        raise ValueError("Nonfinite statistical bin errors")
    return float(np.linalg.norm(grad * err_down)), float(np.linalg.norm(grad * err_up))


def variance_from_cov(grad: np.ndarray, cov: np.ndarray) -> float:
    var = float(grad @ cov @ grad)
    if var < -1.0e-9 * max(1.0, float(np.max(np.abs(cov)))):
        raise ValueError(f"Negative propagated variance: {var}")
    return max(0.0, var)


def mean_from_spectrum(
    spectrum: Dict[str, np.ndarray], name: str, *,
    first_bin_zero: bool, bin_policy: str, negative_policy: str,
    covariance: Optional[np.ndarray] = None,
    syst_covariance: Optional[np.ndarray] = None,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    low = spectrum["low"]
    high = spectrum["high"]
    density = spectrum["density"].copy()
    widths = high - low
    centers = spectrum["center"].copy()
    n = len(centers)
    if not n:
        raise ValueError(f"{name}: empty spectrum")
    original_density = density.copy()
    if np.any(density < 0):
        count = int(np.count_nonzero(density < 0))
        if negative_policy == "error":
            raise ValueError(f"{name}: {count} negative input bins")
        if negative_policy == "clip":
            density = np.maximum(density, 0.0)
        else:
            print(f"[WARN] {name}: {count} negative bins retained; check physical interpretation", file=sys.stderr)

    include = np.ones(n, dtype=bool)
    if bin_policy == "plotter":
        policy = PLOTTER_EDGE_POLICY[name]
        if not policy["first"]:
            include[0] = False
        if not policy["last"]:
            include[-1] = False

    # The discrete approximation for the first bin changes ONLY the numerator.
    if first_bin_zero:
        centers[0] = 0.0

    bin_yield = density * widths
    weighted_yield = centers * bin_yield
    integral = float(np.sum(bin_yield[include]))
    moment = float(np.sum(weighted_yield[include]))
    if not (math.isfinite(integral) and integral > 0):
        raise ValueError(f"{name}: nonpositive/invalid total yield {integral:.6g}")
    mean = moment / integral
    if not math.isfinite(mean) or mean < 0:
        raise ValueError(f"{name}: invalid mean {mean}")

    # Analytic Jacobian of a ratio of two linear combinations.
    # d<lambda>/d(density_i) = width_i * (center_i - mean) / total_yield.
    grad = np.zeros(n, dtype=float)
    grad[include] = widths[include] * (centers[include] - mean) / integral
    if negative_policy == "clip":
        # Clipped bins are locally constant under small fluctuations.
        grad[original_density < 0.0] = 0.0

    if covariance is not None:
        if covariance.shape != (n, n):
            raise ValueError(f"{name}: covariance dimension does not match bin count")
        stat_low = stat_high = math.sqrt(variance_from_cov(grad, covariance))
        stat_source = "bin_covariance"
    else:
        stat_low, stat_high = signed_diagonal_error(
            grad, spectrum["stat_low"], spectrum["stat_high"]
        )
        stat_source = "diagonal_bins"
        if np.all(spectrum["stat_low"] == 0) and np.all(spectrum["stat_high"] == 0):
            print(f"[WARN] {name}: all input bin statistical errors are zero", file=sys.stderr)
    stat_sym = 0.5 * (stat_low + stat_high)

    # Same first-order Jacobian, but the systematic errors are not independent
    # in reality. Without a supplied systematic covariance, this is a useful
    # documented approximation, NOT a rigorous total systematic uncertainty.
    if syst_covariance is not None:
        if syst_covariance.shape != (n, n):
            raise ValueError(f"{name}: systematic covariance dimension does not match bin count")
        syst_low = syst_high = math.sqrt(variance_from_cov(grad, syst_covariance))
        syst_source = "bin_covariance"
    elif np.all(np.isfinite(spectrum["syst_low"]) & np.isfinite(spectrum["syst_high"])):
        syst_low, syst_high = signed_diagonal_error(
            grad, spectrum["syst_low"], spectrum["syst_high"]
        )
        syst_source = "diagonal_bins"
    else:
        syst_low = syst_high = math.nan
        syst_source = "missing"
    syst_sym = 0.5 * (syst_low + syst_high)

    details: List[Dict[str, Any]] = []
    for i in range(n):
        details.append({
            "observable": name,
            "bin_index": i + 1,
            "bin_low": float(low[i]),
            "bin_high": float(high[i]),
            "bin_width": float(widths[i]),
            "bin_center_original": float(spectrum["center"][i]),
            "bin_center_used": float(centers[i]),
            "density_original": float(original_density[i]),
            "density_used": float(density[i]),
            "bin_yield_used": float(bin_yield[i]),
            "bin_moment_used": float(weighted_yield[i]),
            "included": int(include[i]),
            "mean_stat_derivative_density": float(grad[i]),
            "mean_syst_derivative_density": float(grad[i]),
            "stat_low_density": float(spectrum["stat_low"][i]),
            "stat_high_density": float(spectrum["stat_high"][i]),
            "syst_low_density": float(spectrum["syst_low"][i]),
            "syst_high_density": float(spectrum["syst_high"][i]),
        })
    return {
        "alpha": ANGULARITIES[name][0],
        "mean": mean,
        "stat_err_low": stat_low,
        "stat_err_high": stat_high,
        "stat_err_sym": stat_sym,
        "stat_source": stat_source,
        "syst_err_low": syst_low,
        "syst_err_high": syst_high,
        "syst_err_sym": syst_sym,
        "syst_source": syst_source,
        "integral": integral,
        "moment": moment,
        "n_bins": int(np.count_nonzero(include)),
        "n_bins_total": n,
        "first_bin_zero": int(first_bin_zero),
        "bin_policy": bin_policy,
        "negative_policy": negative_policy,
    }, details


def derived_radii(means: Dict[str, float], radius: float) -> Dict[str, float]:
    ans: Dict[str, float] = {}
    for key, num, den in RADIUS_PAIRS:
        n, d = means[num], means[den]
        if d <= 0 or n <= 0:
            ans[key] = math.nan
        else:
            a = ANGULARITIES[num][0]
            b = ANGULARITIES[den][0]
            ans[key] = radius * float((n / d) ** (1.0 / (a - b)))

    if all(means[k] > 0 for k in ("lha", "girth", "gti")):
        ans["Q_girth_squared_over_lha_gti"] = (means["girth"] ** 2) / (means["lha"] * means["gti"])
    else:
        ans["Q_girth_squared_over_lha_gti"] = math.nan
    return ans


def propagated_derived_uncertainty(
    means: Dict[str, float],
    means_errors: Dict[str, Dict[str, float]],
    radius_values: Dict[str, float],
    mean_covariance: Optional[np.ndarray] = None,
    kind: str = "stat",
) -> Dict[str, Any]:
    """First-order propagation of the chosen uncertainty to all radii and Q.

    In diagonal mode, errors on different angularities are treated as
    independent. A supplied mean covariance (5x5) instead retains their
    correlations and produces symmetric errors.
    """
    names = tuple(ANGULARITIES)
    index = {name: i for i, name in enumerate(names)}
    if mean_covariance is not None:
        source = "full_mean_covariance"
    elif any(e.get("source") == "missing" for e in means_errors.values()):
        source = "missing_some_inputs"
    elif any(e.get("source") == "bin_covariance" for e in means_errors.values()):
        source = "independent_angularities_with_bin_covariance"
    else:
        source = "diagonal_bins_and_independent_angularities"
    out: Dict[str, Any] = {f"radius_{kind}_source": source}

    def errors(gradient: np.ndarray) -> Tuple[float, float]:
        if mean_covariance is not None:
            sigma = math.sqrt(variance_from_cov(gradient, mean_covariance))
            return sigma, sigma
        err_lo = np.array([means_errors[name]["low"] for name in names])
        err_hi = np.array([means_errors[name]["high"] for name in names])
        # Missing uncertainty on a variable with ZERO derivative should not
        # invalidate a radius that does not depend on that variable.
        active = np.abs(gradient) > 0.0
        if not np.all(np.isfinite(err_lo[active]) & np.isfinite(err_hi[active])):
            return math.nan, math.nan
        err_lo[~active] = 0.0
        err_hi[~active] = 0.0
        return signed_diagonal_error(gradient, err_lo, err_hi)

    for key, num, den in RADIUS_PAIRS:
        value = float(radius_values[key])
        grad = np.zeros(len(names), dtype=float)
        if math.isfinite(value) and means[num] > 0 and means[den] > 0:
            a, b = ANGULARITIES[num][0], ANGULARITIES[den][0]
            exponent = 1.0 / (a - b)
            grad[index[num]] = exponent * value / means[num]
            grad[index[den]] = -exponent * value / means[den]
            low, high = errors(grad)
        else:
            low = high = math.nan
        out[f"{key}__{kind}_low"] = low
        out[f"{key}__{kind}_high"] = high

    qkey = "Q_girth_squared_over_lha_gti"
    q = float(radius_values[qkey])
    if math.isfinite(q) and all(means[k] > 0 for k in ("lha", "girth", "gti")):
        grad = np.zeros(len(names), dtype=float)
        grad[index["lha"]] = -q / means["lha"]
        grad[index["girth"]] = 2 * q / means["girth"]
        grad[index["gti"]] = -q / means["gti"]
        low, high = errors(grad)
    else:
        low = high = math.nan
    out[f"{qkey}__{kind}_low"] = low
    out[f"{qkey}__{kind}_high"] = high
    return out


def propagated_derived_stat(means, means_stat, radius_values, mean_covariance=None):
    """Backward-compatible statistical wrapper."""
    return propagated_derived_uncertainty(means, means_stat, radius_values, mean_covariance, "stat")


def propagated_derived_syst(means, means_syst, radius_values, mean_covariance=None):
    return propagated_derived_uncertainty(means, means_syst, radius_values, mean_covariance, "syst")


def write_tsv(path: Path, rows: List[Dict[str, Any]], fields: List[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def make_plots(radii: List[Dict[str, Any]], means: List[Dict[str, Any]], outdir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import Rectangle

    sample_names = list(dict.fromkeys(row["sample"] for row in radii))
    markers = ("o", "s", "D", "^", "v", "x")

    # One page per radius pair; both D0 pT selections in the SAME axes.
    pdf_path = outdir / "effective_radii_all_pairs.pdf"
    with PdfPages(pdf_path) as pdf:
        for key, num, den in RADIUS_PAIRS:
            fig, ax = plt.subplots(figsize=(8.2, 5.8))
            for i, sample in enumerate(sample_names):
                rr = sorted(
                    (r for r in radii if r["sample"] == sample),
                    key=lambda r: CENTRALITY_TAGS.index(r["centrality_tag"]),
                )
                if not rr:
                    continue
                xs_nominal = np.asarray([CENTRALITY_MIDPOINTS[CENTRALITY_TAGS.index(r["centrality_tag"])] for r in rr])
                # Cosmetic horizontal offset separates overlapping samples.
                xs = xs_nominal + (i - 0.5 * (len(sample_names) - 1)) * 1.6
                ys = np.asarray([r.get(key, math.nan) for r in rr], dtype=float)
                lows = np.asarray([r[f"{key}__stat_low"] for r in rr], dtype=float)
                highs = np.asarray([r[f"{key}__stat_high"] for r in rr], dtype=float)
                syst_lows = np.asarray([r[f"{key}__syst_low"] for r in rr], dtype=float)
                syst_highs = np.asarray([r[f"{key}__syst_high"] for r in rr], dtype=float)
                # Stat.: narrow vertical bars. Syst.: semi-transparent boxes.
                line = ax.errorbar(
                    xs, ys, yerr=np.array([lows, highs]),
                    marker=markers[i % len(markers)], linewidth=1.8,
                    markersize=6, capsize=3, label=sample, zorder=3,
                )
                color = line.lines[0].get_color()
                for x, y, slo, shi in zip(xs, ys, syst_lows, syst_highs):
                    if np.all(np.isfinite([x, y, slo, shi])):
                        ax.add_patch(Rectangle(
                            (x - 1.2, y - slo), 2.4, slo + shi,
                            facecolor=color, edgecolor=color, linewidth=0.7,
                            alpha=0.20, zorder=2,
                        ))

            a = ANGULARITIES[num][0]
            b = ANGULARITIES[den][0]
            ax.set_title(
                rf"Effective radius from ratios: {ANGULARITY_LABELS[num]}/{ANGULARITY_LABELS[den]} "
                rf"$(\alpha={a},\;\beta={b})$"
            )
            ax.set_xticks(CENTRALITY_MIDPOINTS, CENTRALITY_LABELS)
            ax.set_xlabel("Centrality")
            ax.set_ylabel(r"$\Delta R_{\mathrm{eff}}$")
            ax.grid(alpha=0.2)
            ax.legend(title=r"$p_{\mathrm{T}}^{\mathrm{D}^{0}}$ interval")
            fig.text(
                0.5, 0.013,
                "Bars: stat.  |  Shaded boxes: syst.  |  Diagonal approximation unless covariance supplied",
                ha="center", fontsize=8.5,
            )
            fig.tight_layout(rect=[0, 0.04, 1, 1])
            pdf.savefig(fig)
            plt.close(fig)

    # Keep single-page overview of mean angularities; now with stat bars.
    fig, ax = plt.subplots(figsize=(9, 5.8))
    for si, sample in enumerate(sample_names):
        for obs in ANGULARITIES:
            vals = sorted(
                (r for r in means if r["sample"] == sample and r["observable"] == obs),
                key=lambda r: CENTRALITY_TAGS.index(r["centrality_tag"]),
            )
            xs = [CENTRALITY_MIDPOINTS[CENTRALITY_TAGS.index(r["centrality_tag"])] for r in vals]
            ys = [r["mean"] for r in vals]
            err = [[r["stat_err_low"] for r in vals], [r["stat_err_high"] for r in vals]]
            line = ax.errorbar(xs, ys, yerr=err, marker=markers[si % len(markers)],
                               linewidth=1.3, capsize=2, label=f"{sample}: {obs}", zorder=3)
            color = line.lines[0].get_color()
            for x, r in zip(xs, vals):
                lo, hi = r["syst_err_low"], r["syst_err_high"]
                if np.all(np.isfinite([x, r["mean"], lo, hi])):
                    ax.add_patch(Rectangle(
                        (x - 1.1, r["mean"] - lo), 2.2, lo + hi,
                        facecolor=color, edgecolor=color, linewidth=0.6,
                        alpha=0.16, zorder=2,
                    ))
    ax.set_xticks(CENTRALITY_MIDPOINTS, CENTRALITY_LABELS)
    ax.set_xlabel("Centrality")
    ax.set_ylabel(r"$\langle \lambda_{\alpha}^{1}\rangle$")
    ax.set_title("Mean generalized angularities (bars: stat., boxes: syst.)")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / "mean_angularities.pdf")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Mean angularities and effective radii from final systematic ROOT graphs")
    ap.add_argument("--root", action="append", help="Input ROOT file (repeat to compare two D0 pT selections)")
    ap.add_argument("--label", action="append", help="Name for each --root (e.g. '1-10 GeV/c')")
    ap.add_argument("--output-dir", default=str(DEFAULT_OUTPUT), help="Directory for TSV and PDF output")
    group = ap.add_mutually_exclusive_group()
    group.add_argument("--first-bin-zero", dest="first_bin_zero", action="store_true", help="Treat first angularity bin as a delta at lambda=0")
    group.add_argument("--first-bin-center", dest="first_bin_zero", action="store_false", help="Use ordinary bin center for first bin")
    ap.set_defaults(first_bin_zero=FIRST_BIN_IS_EXACT_ZERO_DEFAULT)
    ap.add_argument("--bin-policy", choices=("full", "plotter"), default=BIN_POLICY_DEFAULT,
                    help="full: all bins (recommended); plotter: hide edge bins like preliminary plots")
    ap.add_argument("--negative-bins", choices=("keep", "clip", "error"), default=NEGATIVE_BINS_DEFAULT)
    ap.add_argument("--jet-radius", type=float, default=R_JET_DEFAULT)
    ap.add_argument("--stat-cov-npz", default=None,
                    help="Optional NPZ with statistical covariance matrices (see script docstring)")
    ap.add_argument("--syst-cov-npz", default=None,
                    help="Optional NPZ with systematic covariance matrices (see script docstring)")
    ap.add_argument("--no-plots", action="store_true", help="Only write TSV files")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    if uproot is None:
        print("[ERROR] 'uproot' is required (same as the existing finalPlot script). Install: pip install uproot", file=sys.stderr)
        return 2
    if args.jet_radius <= 0:
        print("[ERROR] --jet-radius must be positive", file=sys.stderr)
        return 2
    if args.root:
        files = [resolve_path(r) for r in args.root]
        if args.label and len(args.label) != len(files):
            print("[ERROR] Number of --label arguments must equal number of --root arguments", file=sys.stderr)
            return 2
        labels = args.label or [p.stem for p in files]
    else:
        if args.label:
            print("[ERROR] --label requires corresponding --root arguments", file=sys.stderr)
            return 2
        files = [path for path, _ in DEFAULT_SAMPLES]
        labels = [label for _, label in DEFAULT_SAMPLES]
    output = resolve_path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)

    def open_cov_npz(raw: Optional[str], category: str):
        if not raw:
            return None
        cov_path = resolve_path(raw)
        if not cov_path.is_file():
            raise FileNotFoundError(f"{category} covariance NPZ not found: {cov_path}")
        print(f"[COV] {category}: {cov_path}")
        return np.load(cov_path, allow_pickle=False)

    try:
        cov_data = open_cov_npz(args.stat_cov_npz, "statistical")
        syst_cov_data = open_cov_npz(args.syst_cov_npz, "systematic")
    except (FileNotFoundError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2

    all_means: List[Dict[str, Any]] = []
    all_radii: List[Dict[str, Any]] = []
    all_bins: List[Dict[str, Any]] = []
    for sample_index, (path, sample) in enumerate(zip(files, labels)):
        if not path.is_file():
            print(f"[ERROR] ROOT file not found: {path}", file=sys.stderr)
            return 2
        print(f"\n[INPUT] {path}  (sample={sample})")
        with uproot.open(str(path)) as f:
            for tag, centrality in zip(CENTRALITY_TAGS, CENTRALITY_LABELS):
                selected: Dict[str, float] = {}
                selected_stat: Dict[str, Dict[str, float]] = {}
                selected_syst: Dict[str, Dict[str, float]] = {}
                local_mean_rows: Dict[str, Dict[str, Any]] = {}
                integrals: Dict[str, float] = {}
                for obs, (_alpha, directory) in ANGULARITIES.items():
                    try:
                        spectrum = read_spectrum(f, directory, tag)
                        stat_cov = load_covariance(
                            cov_data, f"sample{sample_index}__{tag}__{obs}", len(spectrum["center"])
                        )
                        syst_cov = load_covariance(
                            syst_cov_data, f"sample{sample_index}__{tag}__{obs}", len(spectrum["center"])
                        )
                        mean_info, details = mean_from_spectrum(
                            spectrum, obs,
                            first_bin_zero=args.first_bin_zero,
                            bin_policy=args.bin_policy,
                            negative_policy=args.negative_bins,
                            covariance=stat_cov,
                            syst_covariance=syst_cov,
                        )
                    except (KeyError, ValueError) as exc:
                        print(f"[ERROR] {sample}, {centrality}, {obs}: {exc}", file=sys.stderr)
                        return 2
                    mean_row = {"sample": sample, "input_root": str(path), "centrality": centrality,
                                "centrality_tag": tag, "observable": obs, **mean_info}
                    all_means.append(mean_row)
                    local_mean_rows[obs] = mean_row
                    selected_stat[obs] = {"low": mean_info["stat_err_low"],
                                          "high": mean_info["stat_err_high"],
                                          "source": mean_info["stat_source"]}
                    selected_syst[obs] = {"low": mean_info["syst_err_low"],
                                          "high": mean_info["syst_err_high"],
                                          "source": mean_info["syst_source"]}
                    for item in details:
                        all_bins.append({"sample": sample, "centrality": centrality,
                                         "centrality_tag": tag, "input_root": str(path), **item})
                    selected[obs] = mean_info["mean"]
                    integrals[obs] = mean_info["integral"]
                mean_cov = load_covariance(
                    cov_data, f"sample{sample_index}__{tag}__means", len(ANGULARITIES)
                )
                if mean_cov is not None:
                    for i, obs in enumerate(ANGULARITIES):
                        sigma = math.sqrt(max(float(mean_cov[i, i]), 0.0))
                        local_mean_rows[obs]["stat_err_low"] = sigma
                        local_mean_rows[obs]["stat_err_high"] = sigma
                        local_mean_rows[obs]["stat_err_sym"] = sigma
                        local_mean_rows[obs]["stat_source"] = "full_mean_covariance"
                        selected_stat[obs] = {"low": sigma, "high": sigma,
                                              "source": "full_mean_covariance"}
                syst_mean_cov = load_covariance(
                    syst_cov_data, f"sample{sample_index}__{tag}__means", len(ANGULARITIES)
                )
                if syst_mean_cov is not None:
                    for i, obs in enumerate(ANGULARITIES):
                        sigma = math.sqrt(max(float(syst_mean_cov[i, i]), 0.0))
                        local_mean_rows[obs]["syst_err_low"] = sigma
                        local_mean_rows[obs]["syst_err_high"] = sigma
                        local_mean_rows[obs]["syst_err_sym"] = sigma
                        local_mean_rows[obs]["syst_source"] = "full_mean_covariance"
                        selected_syst[obs] = {"low": sigma, "high": sigma,
                                              "source": "full_mean_covariance"}
                rad = derived_radii(selected, args.jet_radius)
                stat_rad = propagated_derived_stat(selected, selected_stat, rad, mean_cov)
                syst_rad = propagated_derived_syst(selected, selected_syst, rad, syst_mean_cov)
                all_radii.append({"sample": sample, "input_root": str(path), "centrality": centrality,
                                  "centrality_tag": tag, "jet_R": args.jet_radius,
                                  "first_bin_zero": int(args.first_bin_zero),
                                  "bin_policy": args.bin_policy, **rad, **stat_rad, **syst_rad})
                ratios = np.asarray(list(integrals.values()), dtype=float)
                relative_spread = float((max(ratios) - min(ratios)) / np.mean(ratios))
                if relative_spread > 0.05:
                    print(f"[WARN] {sample}, {centrality}: integrals differ across angularities by "
                          f"{relative_spread:.1%}; zero-jet cancellation is approximate", file=sys.stderr)
                radius_summary = ", ".join(
                    f"{key}={rad[key]:.6g}" for key, _, _ in RADIUS_PAIRS
                )
                print(
                    f"  {centrality:>6s}: "
                    f"LHA={selected['lha']:.6g}, girth={selected['girth']:.6g}, "
                    f"GTI={selected['gti']:.6g}, thrust={selected['thrust']:.6g}, jts={selected['jts']:.6g}, "
                    f"{radius_summary}, Q={rad['Q_girth_squared_over_lha_gti']:.5g}"
                )
                if rad["Q_girth_squared_over_lha_gti"] > 1.001:
                    print(f"[WARN] {sample}, {centrality}: Q > 1, check binning/truncation", file=sys.stderr)

    mean_cols = ["sample", "centrality", "centrality_tag", "observable", "alpha", "mean",
                 "stat_err_low", "stat_err_high", "stat_err_sym", "stat_source",
                 "syst_err_low", "syst_err_high", "syst_err_sym", "syst_source", "integral",
                 "moment", "n_bins", "n_bins_total", "first_bin_zero", "bin_policy", "negative_policy", "input_root"]
    rad_cols = ["sample", "centrality", "centrality_tag", "jet_R", "first_bin_zero", "bin_policy"] \
               + [key for key, _, _ in RADIUS_PAIRS] \
               + [f"{key}__stat_{side}" for key, _, _ in RADIUS_PAIRS for side in ("low", "high")] \
               + [f"{key}__syst_{side}" for key, _, _ in RADIUS_PAIRS for side in ("low", "high")] \
               + ["Q_girth_squared_over_lha_gti", "Q_girth_squared_over_lha_gti__stat_low",
                  "Q_girth_squared_over_lha_gti__stat_high", "Q_girth_squared_over_lha_gti__syst_low",
                  "Q_girth_squared_over_lha_gti__syst_high", "radius_stat_source", "radius_syst_source", "input_root"]
    bin_cols = ["sample", "centrality", "centrality_tag", "observable", "bin_index", "bin_low", "bin_high",
                "bin_width", "bin_center_original", "bin_center_used", "density_original", "density_used",
                "bin_yield_used", "bin_moment_used", "included", "mean_stat_derivative_density", "mean_syst_derivative_density", "stat_low_density", "stat_high_density",
                "syst_low_density", "syst_high_density", "input_root"]
    write_tsv(output / "angularity_means.tsv", all_means, mean_cols)
    write_tsv(output / "effective_radii.tsv", all_radii, rad_cols)
    write_tsv(output / "angularity_bins_QA.tsv", all_bins, bin_cols)
    if not args.no_plots:
        make_plots(all_radii, all_means, output)
    print(f"\n[OK] Wrote results to: {output}")
    print("[NOTE] Statistical errors propagated. Without full covariance, numbers are diagonal approximations.")
    print("[NOTE] Final combined per-bin systematics propagated with linearized Jacobians.")
    print("[NOTE] Without syst covariance, results assume independent bins and angularities.")
    print("[NOTE] Source-by-source coherent shifts are preferable for a final systematic result.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
