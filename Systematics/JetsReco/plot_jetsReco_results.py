#!/usr/bin/env python3
import os
import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import uproot

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SCAN_DIR = os.path.join(SCRIPT_DIR, "scanJetsReco")
DEFAULT_SUMMARY = os.path.join(DEFAULT_SCAN_DIR, "summary.tsv")
DEFAULT_OUTPUT_DIR = os.path.join(DEFAULT_SCAN_DIR, "Output")
DEFAULT_OUT = os.path.join(DEFAULT_SCAN_DIR, "jetsReco_results.pdf")

# =========================================================
# Default settings
# =========================================================
CENT_LABELS = {
    0: "0-10%",
    1: "10-40%",
    2: "40-80%",
}

RCP_LABELS = {
    0: r"$R_{CP}$ 0-10/10-40%",
    1: r"$R_{CP}$ 0-10/40-80%",
    2: r"$R_{CP}$ 10-40/40-80%",
}

OBSERVABLES = {
    0: ("z",    r"$z$",                         r"$z$"),
    1: ("l11",  r"$\lambda_{1}^{1}$",          r"$\lambda_{1}^{1}$"),
    2: ("l1p5", r"$\lambda_{1.5}^{1}$",        r"$\lambda_{1.5}^{1}$"),
    3: ("l2",   r"$\lambda_{2}^{1}$",          r"$\lambda_{2}^{1}$"),
    4: ("l3",   r"$\lambda_{3}^{1}$",          r"$\lambda_{3}^{1}$"),
    5: ("l0p5", r"$\lambda_{0.5}^{1}$",        r"$\lambda_{0.5}^{1}$"),
    6: ("pTD",  r"$p_{T}^{D}$",                r"$p_{T}^{D}$"),
}

# summary label -> legend label
VARIATIONS = [
    ("nominal", "nominal"),

    ("jet_rec_efficiency", "jet reco efficiency"),
    ("jet_nHitsFit13",     r"jet $nHitsFit > 13$"),
    ("jet_nHitsFit17",     r"jet $nHitsFit > 17$"),
    ("jet_kTDrop",         r"jet $k_T$ drop"),
    ("jet_DCA2_8",         r"jet DCA $< 2.8$"),
    ("jet_DCA3_2",         r"jet DCA $< 3.2$"),
    ("jet_hadronicCorr",   "jet hadronic correction"),
]


def read_summary(summary_file):
    rows = []
    with open(summary_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def get_label(row):
    return (
        row.get("paper_label")
        or row.get("splot_label")
        or row.get("prior_label")
        or row.get("jets_reco_label")
        or row.get("label")
    )


def build_file_map(rows, output_dir):
    out = {}

    for row in rows:
        run_id = row.get("run_id", "").strip()
        label = get_label(row)

        if label is None:
            continue

        label = label.strip()

        if not run_id or not label:
            continue

        root_file = os.path.join(output_dir, f"OutputSpectra{run_id}.root")

        if os.path.isfile(root_file):
            out[label] = root_file
        else:
            print(f"[warning] Missing file for {label}: {root_file}")

    return out


def get_spectrum_hist_name(obs_idx, cent, iter_tag):
    return f"Lambda{obs_idx}_{cent}_{iter_tag}"


def get_rcp_hist_name(obs_idx, rcp_idx, iter_tag):
    return f"RCP_5_20_Lambda{obs_idx}_{rcp_idx}_{iter_tag}"


def load_hist(root_path, hist_name):
    with uproot.open(root_path) as f:
        if hist_name not in f:
            keys = list(f.keys())
            raise KeyError(
                f"Histogram '{hist_name}' not found in {root_path}\n"
                f"Available keys (first 30): {keys[:30]}"
            )

        h = f[hist_name]
        values = np.asarray(h.values(), dtype=float)

        variances = h.variances()
        if variances is None:
            errors = np.zeros_like(values)
        else:
            errors = np.sqrt(np.asarray(variances, dtype=float))

        edges = np.asarray(h.axis().edges(), dtype=float)
        centers = 0.5 * (edges[:-1] + edges[1:])
        widths = edges[1:] - edges[:-1]

        return centers, values, errors, widths, edges


def safe_ratio(num, den):
    out = np.full_like(num, np.nan, dtype=float)
    mask = den != 0
    out[mask] = num[mask] / den[mask]
    return out

def make_style_map(active_variations):
    """
    Assign unique colors/markers/linestyles to all active variations.
    """
    cmap = plt.get_cmap("tab20")

    markers = ["o", "s", "^", "v", "D", "P", "X", "*", "<", ">", "h", "H", "p", "8", "d"]
    linestyles = ["-", "--", "-.", ":"]

    style_map = {}

    for i, (label_key, _) in enumerate(active_variations):
        if label_key == "nominal":
            style_map[label_key] = {
                "color": "black",
                "marker": "o",
                "linestyle": "-",
                "linewidth": 1.8,
                "markersize": 4,
            }
        else:
            j = i - 1
            style_map[label_key] = {
                "color": cmap(j % 20),
                "marker": markers[j % len(markers)],
                "linestyle": linestyles[(j // len(markers)) % len(linestyles)],
                "linewidth": 1.1,
                "markersize": 3,
            }

    return style_map
    
def plot_variations_page(pdf, file_map, nominal_file, hist_name, xlab, title,
                         y_label, ratio_ymin, ratio_ymax, active_variations):
    x_nom, y_nom, e_nom, w_nom, edges_nom = load_hist(nominal_file, hist_name)

    fig = plt.figure(figsize=(8.5, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)

    style_map = make_style_map(active_variations)

    for i, (label_key, label_pretty) in enumerate(active_variations):
        root_file = file_map[label_key]
        x, y, e, w, edges = load_hist(root_file, hist_name)

        if len(edges) != len(edges_nom) or not np.allclose(edges, edges_nom):
            raise RuntimeError(f"Binning mismatch for {hist_name}, variation {label_key}")

        style = style_map[label_key]

        ax.errorbar(
            x, y, yerr=e,
            fmt=style["marker"] + style["linestyle"],
            ms=style["markersize"],
            lw=style["linewidth"],
            color=style["color"],
            label=label_pretty,
        )

        ratio = safe_ratio(y, y_nom)
        rax.plot(
            x, ratio,
            style["marker"] + style["linestyle"],
            ms=style["markersize"],
            lw=style["linewidth"],
            color=style["color"],
        )

    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, ncol=2)

    rax.axhline(1.0, color="black", lw=1)
    rax.set_xlabel(xlab)
    rax.set_ylabel("ratio")
    rax.set_ylim(ratio_ymin, ratio_ymax)
    rax.grid(True, alpha=0.3)

    plt.setp(ax.get_xticklabels(), visible=False)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Plot unfolded spectra and R_CP for jet-reconstruction systematic variations."
    )
    parser.add_argument("--summary", default=DEFAULT_SUMMARY)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--out", default=DEFAULT_OUT)
    parser.add_argument("--nominal-file", default="", help="Optional nominal OutputSpectra*.root file used as the ratio reference. Required if the summary does not contain a nominal row.")
    parser.add_argument("--iter-tag", default="it3_ICS", help="e.g. it3_ICS; note: 4iter is it3")
    parser.add_argument("--ratio-ymin", type=float, default=0.7)
    parser.add_argument("--ratio-ymax", type=float, default=1.3)
    parser.add_argument("--no-spectra", action="store_true", help="Skip unfolded spectra pages")
    parser.add_argument("--no-rcp", action="store_true", help="Skip R_CP pages")
    parser.add_argument("--strict", action="store_true", help="Require all listed variations")
    args = parser.parse_args()

    print(f"[info] Summary:    {args.summary}")
    print(f"[info] Output dir: {args.output_dir}")
    print(f"[info] Output PDF: {args.out}")

    if not os.path.isfile(args.summary):
        raise FileNotFoundError(f"Missing {args.summary}")

    rows = read_summary(args.summary)
    file_map = build_file_map(rows, args.output_dir)

    # JetsReco scans usually contain only the variations 21..27.
    # In that case, provide the nominal reference explicitly, e.g.
    #   --nominal-file ../sWeight/scanSWeight/Output/OutputSpectrar000001.root
    # or any other nominal OutputSpectra*.root produced with the same binning.
    if "nominal" not in file_map:
        if args.nominal_file:
            if not os.path.isfile(args.nominal_file):
                raise FileNotFoundError(f"Missing --nominal-file: {args.nominal_file}")
            file_map["nominal"] = args.nominal_file
        else:
            raise RuntimeError(
                "Missing nominal variation. JetsReco summary usually has only SYS=21..27, "
                "so pass --nominal-file /path/to/nominal/OutputSpectra*.root."
            )

    missing_labels = [lab for lab, _ in VARIATIONS if lab not in file_map]

    if missing_labels:
        msg = f"Missing variations: {missing_labels}"
        if args.strict:
            raise RuntimeError(msg)
        print(f"[warning] {msg}")
        print("[warning] Missing variations will be skipped.")

    active_variations = [
        (lab, pretty) for lab, pretty in VARIATIONS
        if lab in file_map
    ]

    nominal_file = file_map["nominal"]

    out_parent = os.path.dirname(os.path.abspath(args.out))
    if out_parent:
        os.makedirs(out_parent, exist_ok=True)

    with PdfPages(args.out) as pdf:
        if not args.no_spectra:
            for obs_idx, (obs_short, obs_pretty, xlab) in OBSERVABLES.items():
                for cent in [0, 1, 2]:
                    hist_name = get_spectrum_hist_name(obs_idx, cent, args.iter_tag)
                    title = f"{obs_pretty}, cent {CENT_LABELS[cent]}, {args.iter_tag}"

                    plot_variations_page(
                        pdf=pdf,
                        file_map=file_map,
                        nominal_file=nominal_file,
                        hist_name=hist_name,
                        xlab=xlab,
                        title=title,
                        y_label="Unfolded yield",
                        ratio_ymin=args.ratio_ymin,
                        ratio_ymax=args.ratio_ymax,
                        active_variations=active_variations,
                    )

                    print(f"[ok] {hist_name}")

        if not args.no_rcp:
            for obs_idx, (obs_short, obs_pretty, xlab) in OBSERVABLES.items():
                for rcp_idx in [0, 1, 2]:
                    hist_name = get_rcp_hist_name(obs_idx, rcp_idx, args.iter_tag)
                    title = (
                        f"$R_{{CP}}$, 5 < $p_{{T,Jet}}$ < 20 GeV/$c$, "
                        f"{obs_pretty}, {RCP_LABELS[rcp_idx]}, {args.iter_tag}"
                    )

                    plot_variations_page(
                        pdf=pdf,
                        file_map=file_map,
                        nominal_file=nominal_file,
                        hist_name=hist_name,
                        xlab=xlab,
                        title=title,
                        y_label=r"$R_{CP}$",
                        ratio_ymin=args.ratio_ymin,
                        ratio_ymax=args.ratio_ymax,
                        active_variations=active_variations,
                    )

                    print(f"[ok] {hist_name}")

    print(f"\nSaved: {args.out}")


if __name__ == "__main__":
    main()
