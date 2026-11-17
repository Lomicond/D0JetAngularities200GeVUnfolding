#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Effective radius from the SAME final ROOT file as
Systematics/finalPlot_from_systematics_root_v16_lido_spectra.py.

Stage 1: central values and bin-by-bin QA only. No propagated uncertainties yet.

By DEFAULT runs BOTH D0 pT ranges (1-10 and 5-10 GeV/c) from the
Systematics directory, writing a combined comparison table and plots.

Run from the unfolding project root, for example:
  python3 Systematics/effective_radius_from_final_root_v2.py
  python3 Systematics/effective_radius_from_final_root_v2.py --first-bin-zero
  python3 Systematics/effective_radius_from_final_root_v2.py --bin-policy plotter

Compare two separate final-systematics files:
  python3 Systematics/effective_radius_from_final_root_v2.py \
    --root /path/to/one_to_ten.root --label '1-10 GeV/c' \
    --root /path/to/five_to_ten.root --label '5-10 GeV/c' \
    --first-bin-zero

Input graph example: l051/g_value_stat_0_10pct
Output: three TSV tables and two PDF plots in Systematics/EffectiveRadius/.

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

# The existing plotter SHOW_EDGE_BINS switch (angularities have last=False).
PLOTTER_EDGE_POLICY = {
    "lha": {"first": True, "last": False},
    "girth": {"first": True, "last": False},
    "gti": {"first": True, "last": False},
    "thrust": {"first": True, "last": False},
    "jts": {"first": True, "last": False},
}

# (output column, numerator angularity, denominator angularity)
RADIUS_PAIRS = (
    ("girth_over_lha", "girth", "lha"),
    ("gti_over_lha", "gti", "lha"),
    ("gti_over_girth", "gti", "girth"),
)
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

    # Optional syst arrays are captured in QA now, not propagated in stage 1.
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

    return {
        "low": low, "high": high, "center": x, "density": y,
        "stat_low": stat_low, "stat_high": stat_high,
        "syst_low": syst_low, "syst_high": syst_high,
    }


def mean_from_spectrum(
    spectrum: Dict[str, np.ndarray], name: str, *,
    first_bin_zero: bool, bin_policy: str, negative_policy: str,
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
            "stat_low_density": float(spectrum["stat_low"][i]),
            "stat_high_density": float(spectrum["stat_high"][i]),
            "syst_low_density": float(spectrum["syst_low"][i]),
            "syst_high_density": float(spectrum["syst_high"][i]),
        })
    return {
        "alpha": ANGULARITIES[name][0],
        "mean": mean,
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


def write_tsv(path: Path, rows: List[Dict[str, Any]], fields: List[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def make_plots(radii: List[Dict[str, Any]], means: List[Dict[str, Any]], outdir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sample_names = list(dict.fromkeys(row["sample"] for row in radii))
    fig, ax = plt.subplots(figsize=(8, 5.2))
    pairs = (("girth_over_lha", r"$R_{\mathrm{eff}}^{(1,0.5)}$"),
             ("gti_over_lha", r"$R_{\mathrm{eff}}^{(1.5,0.5)}$"),
             ("gti_over_girth", r"$R_{\mathrm{eff}}^{(1.5,1)}$"))
    for sample in sample_names:
        rr = [r for r in radii if r["sample"] == sample]
        rr.sort(key=lambda r: CENTRALITY_TAGS.index(r["centrality_tag"]))
        xs = [CENTRALITY_MIDPOINTS[CENTRALITY_TAGS.index(r["centrality_tag"])] for r in rr]
        for key, caption in pairs:
            ys = [r[key] for r in rr]
            ax.plot(xs, ys, marker="o", label=f"{sample}: {caption}")
    ax.set_xticks(CENTRALITY_MIDPOINTS, CENTRALITY_LABELS)
    ax.set_xlabel("Centrality")
    ax.set_ylabel(r"$\Delta R_{\mathrm{eff}}$")
    ax.set_title("Effective radius from ratios of mean angularities")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "effective_radii.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5.2))
    for sample in sample_names:
        for obs in ("lha", "girth", "gti", "thrust", "jts"):
            vals = [r for r in means if r["sample"] == sample and r["observable"] == obs]
            vals.sort(key=lambda r: CENTRALITY_TAGS.index(r["centrality_tag"]))
            xs = [CENTRALITY_MIDPOINTS[CENTRALITY_TAGS.index(r["centrality_tag"])] for r in vals]
            ys = [r["mean"] for r in vals]
            ax.plot(xs, ys, marker="o", label=f"{sample}: {obs}")
    ax.set_xticks(CENTRALITY_MIDPOINTS, CENTRALITY_LABELS)
    ax.set_xlabel("Centrality")
    ax.set_ylabel(r"$\langle \lambda_{\alpha}^{1}\rangle$")
    ax.set_title("Mean generalized angularities")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
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

    all_means: List[Dict[str, Any]] = []
    all_radii: List[Dict[str, Any]] = []
    all_bins: List[Dict[str, Any]] = []
    for path, sample in zip(files, labels):
        if not path.is_file():
            print(f"[ERROR] ROOT file not found: {path}", file=sys.stderr)
            return 2
        print(f"\n[INPUT] {path}  (sample={sample})")
        with uproot.open(str(path)) as f:
            for tag, centrality in zip(CENTRALITY_TAGS, CENTRALITY_LABELS):
                selected: Dict[str, float] = {}
                integrals: Dict[str, float] = {}
                for obs, (_alpha, directory) in ANGULARITIES.items():
                    try:
                        spectrum = read_spectrum(f, directory, tag)
                        mean_info, details = mean_from_spectrum(
                            spectrum, obs,
                            first_bin_zero=args.first_bin_zero,
                            bin_policy=args.bin_policy,
                            negative_policy=args.negative_bins,
                        )
                    except (KeyError, ValueError) as exc:
                        print(f"[ERROR] {sample}, {centrality}, {obs}: {exc}", file=sys.stderr)
                        return 2
                    mean_row = {"sample": sample, "input_root": str(path), "centrality": centrality,
                                "centrality_tag": tag, "observable": obs, **mean_info}
                    all_means.append(mean_row)
                    for item in details:
                        all_bins.append({"sample": sample, "centrality": centrality,
                                         "centrality_tag": tag, "input_root": str(path), **item})
                    selected[obs] = mean_info["mean"]
                    integrals[obs] = mean_info["integral"]
                rad = derived_radii(selected, args.jet_radius)
                all_radii.append({"sample": sample, "input_root": str(path), "centrality": centrality,
                                  "centrality_tag": tag, "jet_R": args.jet_radius,
                                  "first_bin_zero": int(args.first_bin_zero),
                                  "bin_policy": args.bin_policy, **rad})
                ratios = np.asarray(list(integrals.values()), dtype=float)
                relative_spread = float((max(ratios) - min(ratios)) / np.mean(ratios))
                if relative_spread > 0.05:
                    print(f"[WARN] {sample}, {centrality}: integrals differ across angularities by "
                          f"{relative_spread:.1%}; zero-jet cancellation is approximate", file=sys.stderr)
                print(f"  {centrality:>6s}: LHA={selected['lha']:.6g}, girth={selected['girth']:.6g}, "
                      f"GTI={selected['gti']:.6g}, R(1/.5)={rad['girth_over_lha']:.6g}, "
                      f"R(1.5/.5)={rad['gti_over_lha']:.6g}, Q={rad['Q_girth_squared_over_lha_gti']:.5g}")
                if rad["Q_girth_squared_over_lha_gti"] > 1.001:
                    print(f"[WARN] {sample}, {centrality}: Q > 1, check binning/truncation", file=sys.stderr)

    mean_cols = ["sample", "centrality", "centrality_tag", "observable", "alpha", "mean", "integral",
                 "moment", "n_bins", "n_bins_total", "first_bin_zero", "bin_policy", "negative_policy", "input_root"]
    rad_cols = ["sample", "centrality", "centrality_tag", "jet_R", "first_bin_zero", "bin_policy",
                "girth_over_lha", "gti_over_lha", "gti_over_girth", "Q_girth_squared_over_lha_gti", "input_root"]
    bin_cols = ["sample", "centrality", "centrality_tag", "observable", "bin_index", "bin_low", "bin_high",
                "bin_width", "bin_center_original", "bin_center_used", "density_original", "density_used",
                "bin_yield_used", "bin_moment_used", "included", "stat_low_density", "stat_high_density",
                "syst_low_density", "syst_high_density", "input_root"]
    write_tsv(output / "angularity_means.tsv", all_means, mean_cols)
    write_tsv(output / "effective_radii.tsv", all_radii, rad_cols)
    write_tsv(output / "angularity_bins_QA.tsv", all_bins, bin_cols)
    if not args.no_plots:
        make_plots(all_radii, all_means, output)
    print(f"\n[OK] Wrote results to: {output}")
    print("[NOTE] Stage 1 only: statistical and systematic uncertainties are NOT propagated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
