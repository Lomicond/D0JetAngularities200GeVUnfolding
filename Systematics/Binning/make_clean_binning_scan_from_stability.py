#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
from pathlib import Path
from collections import Counter, defaultdict

KEY_COLS = [
    "reco_pt_min",
    "reco_pt_max",
    "reco_n_bins",
    "true_pt_min",
    "true_pt_max",
    "true_n_bins",
    "min_width",
    "max_width",
    "step",
    "trend",
    "reco_edges",
    "true_edges",
]


def run_num(run_id):
    m = re.search(r"(\d+)$", str(run_id))
    return int(m.group(1)) if m else -1


def read_summary(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        return rows, reader.fieldnames


def config_key(row):
    return tuple(str(row.get(c, "")).strip() for c in KEY_COLS)


def read_stability(path):
    rows = []
    nf_counter = Counter()

    with open(path, "r", errors="replace") as f:
        for nr, line in enumerate(f, 1):
            line = line.rstrip("\n")
            fields = line.split("\t")
            nf_counter[len(fields)] += 1
            rows.append(fields)

    return rows, nf_counter


def write_dict_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in fieldnames})


def write_rows_tsv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scan_dir")
    ap.add_argument("--expected-total", type=int, default=9795)
    args = ap.parse_args()

    scan_dir = Path(args.scan_dir).resolve()
    summary_path = scan_dir / "summary.tsv"
    stability_path = scan_dir / "stability.tsv"

    out_dir = scan_dir / "clean_from_stability"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows, summary_cols = read_summary(summary_path)
    stability_rows, nf_counter = read_stability(stability_path)

    missing = [c for c in KEY_COLS if c not in summary_cols]
    if missing:
        raise SystemExit("[error] Missing columns in summary.tsv: " + ", ".join(missing))

    # Count stability rows per run.
    stability_count = Counter()
    for r in stability_rows:
        if r:
            stability_count[r[0]] += 1

    nonzero = list(stability_count.values())
    expected_rows_per_run = Counter(nonzero).most_common(1)[0][0]

    # Group runs by physical binning configuration.
    groups = defaultdict(list)
    row_by_run = {}

    for row in summary_rows:
        rid = row["run_id"]
        row_by_run[rid] = row
        groups[config_key(row)].append(rid)

    chosen_runs = []
    unresolved = []
    duplicates = []

    for key, runs in groups.items():
        complete = [r for r in runs if stability_count.get(r, 0) == expected_rows_per_run]
        complete = sorted(complete, key=run_num)

        if complete:
            chosen = complete[0]
            chosen_runs.append(chosen)

            for r in runs:
                if r != chosen:
                    duplicates.append({
                        "run_id": r,
                        "chosen_instead": chosen,
                        "stability_rows": stability_count.get(r, 0),
                        "expected_stability_rows": expected_rows_per_run,
                    })
        else:
            unresolved.append({
                "all_run_ids": ",".join(sorted(runs, key=run_num)),
                "max_stability_rows": max(stability_count.get(r, 0) for r in runs),
                "expected_stability_rows": expected_rows_per_run,
                **{c: v for c, v in zip(KEY_COLS, key)},
            })

    chosen_set = set(chosen_runs)

    # Clean summary: one chosen complete run per binning configuration.
    clean_summary = [row_by_run[r] for r in sorted(chosen_runs, key=run_num)]

    # Clean stability: only chosen runs.
    clean_stability = [r for r in stability_rows if r and r[0] in chosen_set]

    write_dict_tsv(out_dir / "summary_clean.tsv", clean_summary, summary_cols)
    write_rows_tsv(out_dir / "stability_clean.tsv", clean_stability)

    with open(out_dir / "kept_runs.list", "w") as f:
        for r in sorted(chosen_runs, key=run_num):
            f.write(r + "\n")

    write_dict_tsv(
        out_dir / "duplicates_ignored.tsv",
        duplicates,
        ["run_id", "chosen_instead", "stability_rows", "expected_stability_rows"],
    )

    unresolved_cols = [
        "all_run_ids",
        "max_stability_rows",
        "expected_stability_rows",
    ] + KEY_COLS

    write_dict_tsv(out_dir / "unresolved_configs.tsv", unresolved, unresolved_cols)

    print("")
    print("========== Clean from stability ==========")
    print(f"scan_dir:                    {scan_dir}")
    print(f"summary rows:                {len(summary_rows)}")
    print(f"unique bin configurations:   {len(groups)}")
    print(f"expected total configs:      {args.expected_total}")
    print(f"stability NF counts:         {dict(sorted(nf_counter.items()))}")
    print(f"expected stability rows/run: {expected_rows_per_run}")
    print(f"chosen complete configs:     {len(chosen_runs)}")
    print(f"unresolved configs:          {len(unresolved)}")
    print(f"ignored duplicate runs:      {len(duplicates)}")
    print("")
    print(f"Clean summary:               {out_dir / 'summary_clean.tsv'}")
    print(f"Clean stability:             {out_dir / 'stability_clean.tsv'}")
    print(f"Kept runs:                   {out_dir / 'kept_runs.list'}")
    print(f"Unresolved configs:          {out_dir / 'unresolved_configs.tsv'}")
    print(f"Ignored duplicates:          {out_dir / 'duplicates_ignored.tsv'}")

    if len(chosen_runs) == args.expected_total and not unresolved:
        print("")
        print("[OK] Čisté soubory by měly být použitelné v browseru.")
    else:
        print("")
        print("[warning] Něco nesedí. Podívej se hlavně do unresolved_configs.tsv.")

    print("==========================================")
    print("")


if __name__ == "__main__":
    main()
