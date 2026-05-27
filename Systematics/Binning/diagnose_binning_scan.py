#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
from pathlib import Path
from collections import Counter, defaultdict

SERIOUS_LOG_PATTERNS = [
    "no space left on device",
    "disk quota exceeded",
    "cannot allocate memory",
    "std::bad_alloc",
    "bad_alloc",
    "killed",
    "segmentation violation",
    "segmentation fault",
    "bus error",
    "aborted",
    "cgroup memory",
    "memory limit",
    "fatal",
]


def run_num(run_id):
    m = re.search(r"(\d+)$", str(run_id))
    return int(m.group(1)) if m else -1


def read_summary(summary_path):
    rows = []
    with open(summary_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader, start=2):
            row["_line"] = i
            rows.append(row)
    return rows


def detect_stability_columns(stability_path):
    cnt = Counter()
    bad = []
    with open(stability_path, "r", errors="replace") as f:
        for nr, line in enumerate(f, start=1):
            nf = len(line.rstrip("\n").split("\t"))
            cnt[nf] += 1
            if nf not in (9, 10):
                bad.append((nr, nf, line.rstrip("\n")[:250]))

    if cnt[10] >= cnt[9]:
        expected = 10
    else:
        expected = 9

    return expected, cnt, bad


def split_glued_stability_line(fields, expected_cols):
    """
    Typical crash/write collision:
      row1_col_last + row2_run_id are glued in fields[expected_cols-1]
    For expected_cols=10 this gives 19 fields instead of 20.
    """
    if len(fields) != 2 * expected_cols - 1:
        return None

    glued = fields[expected_cols - 1]
    m = re.match(r"^(.*?)(r[0-9]{6})$", glued)
    if not m:
        return None

    row1 = fields[:expected_cols - 1] + [m.group(1)]
    row2 = [m.group(2)] + fields[expected_cols:]

    if len(row1) != expected_cols or len(row2) != expected_cols:
        return None

    return row1, row2


def read_and_fix_stability(stability_path, out_fixed_path=None):
    expected_cols, nf_counter, bad_lines = detect_stability_columns(stability_path)

    rows = []
    fixed = []
    not_fixed = []

    with open(stability_path, "r", errors="replace") as f:
        for nr, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            fields = line.split("\t")

            if len(fields) == expected_cols:
                rows.append(fields)
                fixed.append(fields)
                continue

            if len(fields) in (9, 10):
                # Mixed 9/10-column file. Keep for diagnostics, but report later.
                rows.append(fields)
                fixed.append(fields)
                continue

            split = split_glued_stability_line(fields, expected_cols)
            if split is not None:
                r1, r2 = split
                rows.append(r1)
                rows.append(r2)
                fixed.append(r1)
                fixed.append(r2)
            else:
                not_fixed.append((nr, len(fields), line[:300]))

    if out_fixed_path is not None:
        with open(out_fixed_path, "w", newline="") as g:
            w = csv.writer(g, delimiter="\t", lineterminator="\n")
            for r in fixed:
                w.writerow(r)

    return rows, expected_cols, nf_counter, bad_lines, not_fixed


def config_key(row):
    key_cols = [
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
    return tuple(str(row.get(c, "")).strip() for c in key_cols)


def find_root_file(scan_dir, run_id):
    out_dir = scan_dir / "Output"

    candidates = [
        out_dir / f"OutputSpectra{run_id}.root",
        out_dir / f"OutputSpectra{run_id}.root.root",
        out_dir / f"{run_id}.root",
    ]

    for p in candidates:
        if p.exists() and p.stat().st_size > 0:
            return p

    matches = sorted(out_dir.glob(f"*{run_id}*.root"))
    matches = [p for p in matches if p.exists() and p.stat().st_size > 0]
    if matches:
        return matches[0]

    return None


def read_tail(path, max_bytes=200000):
    try:
        with open(path, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - max_bytes))
            return f.read().decode("utf-8", errors="replace")
    except Exception:
        return ""


def analyze_log(log_path):
    if not log_path or not Path(log_path).exists():
        return ["missing root.log"]

    txt = read_tail(log_path).lower()
    hits = []

    for pat in SERIOUS_LOG_PATTERNS:
        if pat in txt:
            hits.append(pat)

    return hits


def write_tsv(path, rows, columns):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=columns, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in columns})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scan_dir", help="Např. Systematics/Binning/scanJetPt")
    ap.add_argument("--no-root-check", action="store_true")
    args = ap.parse_args()

    scan_dir = Path(args.scan_dir).resolve()
    summary_path = scan_dir / "summary.tsv"
    stability_path = scan_dir / "stability.tsv"

    if not summary_path.exists():
        raise SystemExit(f"[error] Missing {summary_path}")
    if not stability_path.exists():
        raise SystemExit(f"[error] Missing {stability_path}")

    diag_dir = scan_dir / "diagnostics"
    diag_dir.mkdir(parents=True, exist_ok=True)

    fixed_stability_path = diag_dir / "stability_fixed.tsv"
    stability_rows, expected_cols, nf_counter, bad_lines, not_fixed = read_and_fix_stability(
        stability_path,
        out_fixed_path=fixed_stability_path,
    )

    stab_count = Counter()
    for r in stability_rows:
        if r:
            stab_count[r[0]] += 1

    nonzero_counts = [v for v in stab_count.values() if v > 0]
    if nonzero_counts:
        expected_stab_rows = Counter(nonzero_counts).most_common(1)[0][0]
    else:
        expected_stab_rows = 0

    summary_rows = read_summary(summary_path)
    summary_rows.sort(key=lambda r: run_num(r.get("run_id", "")))

    seen_config = {}
    duplicate_info = {}
    first_duplicate = None

    for idx, row in enumerate(summary_rows):
        rid = row.get("run_id", "")
        key = config_key(row)

        if key in seen_config:
            duplicate_info[rid] = seen_config[key]
            if first_duplicate is None:
                first_duplicate = idx
        else:
            seen_config[key] = rid

    status_rows = []
    first_problem = None

    for idx, row in enumerate(summary_rows):
        rid = row.get("run_id", "")
        issues = []

        nstab = stab_count.get(rid, 0)

        if expected_stab_rows > 0 and nstab != expected_stab_rows:
            issues.append(f"stability_rows={nstab}, expected={expected_stab_rows}")

        log_path = row.get("root_log", "")
        if not log_path:
            log_path = str(scan_dir / "runs" / rid / "root.log")

        log_hits = analyze_log(log_path)
        for h in log_hits:
            issues.append(f"log:{h}")

        root_file = ""
        if not args.no_root_check:
            rf = find_root_file(scan_dir, rid)
            if rf is None:
                issues.append("missing Output/OutputSpectra*.root")
            else:
                root_file = str(rf)

        dup_of = duplicate_info.get(rid, "")

        if issues:
            status = "BAD"
        else:
            status = "OK"

        if dup_of:
            status = "DUPLICATE_" + status

        if first_problem is None and issues:
            first_problem = idx

        status_rows.append({
            "summary_line": row.get("_line", ""),
            "run_id": rid,
            "status": status,
            "duplicate_of": dup_of,
            "stability_rows": nstab,
            "expected_stability_rows": expected_stab_rows,
            "issues": "; ".join(issues),
            "reco_pt_min": row.get("reco_pt_min", ""),
            "reco_pt_max": row.get("reco_pt_max", ""),
            "reco_n_bins": row.get("reco_n_bins", ""),
            "true_pt_min": row.get("true_pt_min", ""),
            "true_pt_max": row.get("true_pt_max", ""),
            "true_n_bins": row.get("true_n_bins", ""),
            "reco_edges": row.get("reco_edges", ""),
            "true_edges": row.get("true_edges", ""),
            "root_log": log_path,
            "root_file": root_file,
        })

    columns = [
        "summary_line",
        "run_id",
        "status",
        "duplicate_of",
        "stability_rows",
        "expected_stability_rows",
        "issues",
        "reco_pt_min",
        "reco_pt_max",
        "reco_n_bins",
        "true_pt_min",
        "true_pt_max",
        "true_n_bins",
        "reco_edges",
        "true_edges",
        "root_log",
        "root_file",
    ]

    write_tsv(diag_dir / "run_status.tsv", status_rows, columns)

    if first_problem is not None:
        write_tsv(diag_dir / "tail_from_first_problem.tsv", status_rows[first_problem:], columns)

    if first_duplicate is not None:
        write_tsv(diag_dir / "tail_from_first_duplicate.tsv", status_rows[first_duplicate:], columns)

    # Recommended keep/drop: keep first OK run for each bin configuration.
    by_key = defaultdict(list)
    row_by_run = {r["run_id"]: r for r in status_rows}

    for row in summary_rows:
        by_key[config_key(row)].append(row.get("run_id", ""))

    keep = set()
    drop = []

    for key, runs in by_key.items():
        ok_runs = [r for r in runs if row_by_run[r]["status"] in ("OK", "DUPLICATE_OK")]
        if ok_runs:
            chosen = sorted(ok_runs, key=run_num)[0]
        else:
            chosen = sorted(runs, key=lambda r: (row_by_run[r]["stability_rows"], run_num(r)), reverse=True)[0]

        keep.add(chosen)

        for r in runs:
            if r != chosen:
                drop.append({
                    "run_id": r,
                    "kept_instead": chosen,
                    "status": row_by_run[r]["status"],
                    "duplicate_of": row_by_run[r]["duplicate_of"],
                    "reco_edges": row_by_run[r]["reco_edges"],
                    "true_edges": row_by_run[r]["true_edges"],
                })

    write_tsv(
        diag_dir / "recommended_duplicate_drops.tsv",
        drop,
        ["run_id", "kept_instead", "status", "duplicate_of", "reco_edges", "true_edges"],
    )

    # Detect whether duplicate tail looks like restart from the beginning.
    restart_msg = "No duplicate tail detected."
    if first_duplicate is not None:
        dup_tail = status_rows[first_duplicate:]
        if dup_tail:
            first_dup = dup_tail[0]
            if first_dup["duplicate_of"] == status_rows[0]["run_id"]:
                n_seq = 0
                for j, r in enumerate(dup_tail):
                    if j < len(status_rows) and r["duplicate_of"] == status_rows[j]["run_id"]:
                        n_seq += 1
                    else:
                        break
                restart_msg = (
                    f"Duplicate tail starts at {first_dup['run_id']} and duplicates from "
                    f"{first_dup['duplicate_of']}. Sequential restart-like duplicates: {n_seq} runs."
                )
            else:
                restart_msg = (
                    f"First duplicate is {first_dup['run_id']} duplicating {first_dup['duplicate_of']}, "
                    "but it is not a clean restart from r000001."
                )

    print("")
    print("========== Scan diagnostics ==========")
    print(f"scan_dir: {scan_dir}")
    print(f"summary rows: {len(summary_rows)}")
    print(f"unique bin configurations: {len(by_key)}")
    print(f"stability NF counts: {dict(sorted(nf_counter.items()))}")
    print(f"expected stability columns: {expected_cols}")
    print(f"expected stability rows per successful run: {expected_stab_rows}")
    print(f"fixed stability written to: {fixed_stability_path}")

    if bad_lines:
        print("")
        print("[stability.tsv] malformed lines:")
        for nr, nf, txt in bad_lines[:10]:
            print(f"  line {nr}: NF={nf}: {txt[:160]}")
        if len(bad_lines) > 10:
            print(f"  ... plus {len(bad_lines)-10} more")

    if not_fixed:
        print("")
        print("[warning] Some malformed stability lines were NOT fixed:")
        for nr, nf, txt in not_fixed[:10]:
            print(f"  line {nr}: NF={nf}: {txt[:160]}")

    print("")
    if first_problem is None:
        print("First BAD run: none found by these checks.")
    else:
        r = status_rows[first_problem]
        print(f"First BAD run: {r['run_id']}  summary line {r['summary_line']}")
        print(f"  status: {r['status']}")
        print(f"  issues: {r['issues']}")
        print(f"  reco edges: {r['reco_edges']}")
        print(f"  true edges: {r['true_edges']}")
        print(f"  tail table: {diag_dir / 'tail_from_first_problem.tsv'}")

    if first_duplicate is None:
        print("First duplicate bin configuration: none")
    else:
        r = status_rows[first_duplicate]
        print(f"First duplicate bin configuration: {r['run_id']} duplicates {r['duplicate_of']}")
        print(f"  reco edges: {r['reco_edges']}")
        print(f"  true edges: {r['true_edges']}")
        print(f"  duplicate tail table: {diag_dir / 'tail_from_first_duplicate.tsv'}")

    print("")
    print(restart_msg)
    print("")
    print(f"Full status table: {diag_dir / 'run_status.tsv'}")
    print(f"Recommended duplicate drops: {diag_dir / 'recommended_duplicate_drops.tsv'}")
    print("======================================")
    print("")


if __name__ == "__main__":
    main()
