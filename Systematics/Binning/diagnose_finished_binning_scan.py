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

BAD_LOG_PATTERNS = [
    "no space left on device",
    "disk quota exceeded",
    "cannot allocate memory",
    "std::bad_alloc",
    "bad_alloc",
    "cgroup memory limit",
    "memory limit",
    "segmentation violation",
    "segmentation fault",
    "bus error",
    "aborted",
    "killed",
]


def run_num(run_id):
    m = re.search(r"(\d+)$", str(run_id))
    return int(m.group(1)) if m else -1


def read_summary(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = []
        for line, row in enumerate(reader, start=2):
            row["_summary_line"] = line
            rows.append(row)
        return rows, reader.fieldnames


def config_key(row, key_cols):
    return tuple(str(row.get(c, "")).strip() for c in key_cols)


def detect_stability_expected_cols(path):
    nf_counter = Counter()
    with open(path, "r", errors="replace") as f:
        for line in f:
            nf = len(line.rstrip("\n").split("\t"))
            nf_counter[nf] += 1

    # Browser podporuje 9 nebo 10 sloupců, pokud je někde 10, sjednotíme na 10.
    if nf_counter.get(10, 0) > 0:
        expected_cols = 10
    elif nf_counter.get(9, 0) > 0:
        expected_cols = 9
    else:
        expected_cols = nf_counter.most_common(1)[0][0] if nf_counter else 10

    return expected_cols, nf_counter


def split_glued_line(fields, expected_cols):
    """
    Oprava typického slepení:
      ... posledni_sloupec_runu1 + rXXXXXX ... run2 ...
    Pro 10 sloupců vznikne NF=19 místo 20.
    """
    if len(fields) != 2 * expected_cols - 1:
        return None

    glued = fields[expected_cols - 1]
    m = re.match(r"^(.*?)(r[0-9]{6})$", glued)
    if not m:
        return None

    row1 = fields[:expected_cols - 1] + [m.group(1)]
    row2 = [m.group(2)] + fields[expected_cols:]

    if len(row1) == expected_cols and len(row2) == expected_cols:
        return row1, row2

    return None


def read_stability_fixed(path, expected_cols):
    rows = []
    malformed = []
    fixed_glued = 0
    padded = 0

    with open(path, "r", errors="replace") as f:
        for line_nr, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            fields = line.split("\t")

            if len(fields) == expected_cols:
                rows.append(fields)
                continue

            # Smíšený 9/10 sloupcový soubor: doplníme prázdný poslední sloupec.
            if expected_cols == 10 and len(fields) == 9:
                rows.append(fields + [""])
                padded += 1
                continue

            split = split_glued_line(fields, expected_cols)
            if split is not None:
                r1, r2 = split
                rows.append(r1)
                rows.append(r2)
                fixed_glued += 1
                continue

            malformed.append({
                "line": line_nr,
                "nf": len(fields),
                "text": line[:300],
            })

    return rows, malformed, fixed_glued, padded


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
    p = Path(log_path)
    if not p.exists():
        return ["missing root.log"]

    txt = read_tail(p).lower()
    hits = []
    for pat in BAD_LOG_PATTERNS:
        if pat in txt:
            hits.append(pat)
    return hits


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

    if out_dir.exists():
        matches = sorted(out_dir.glob(f"*{run_id}*.root*"))
        matches = [p for p in matches if p.exists() and p.stat().st_size > 0]
        if matches:
            return matches[0]

    return None


def write_tsv_dict(path, rows, fieldnames):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in fieldnames})


def write_stability(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        for r in rows:
            w.writerow(r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scan_dir", help="Např. Systematics/Binning/scanJetPt")
    ap.add_argument("--expected-total", type=int, default=9795,
                    help="Očekávaný počet unikátních konfigurací. Default: 9795")
    ap.add_argument("--no-root-check", action="store_true",
                    help="Nekontrolovat existenci OutputSpectra*.root")
    args = ap.parse_args()

    scan_dir = Path(args.scan_dir).resolve()
    summary_path = scan_dir / "summary.tsv"
    stability_path = scan_dir / "stability.tsv"

    if not summary_path.exists():
        raise SystemExit(f"[error] Missing {summary_path}")
    if not stability_path.exists():
        raise SystemExit(f"[error] Missing {stability_path}")

    diag_dir = scan_dir / "diagnostics_finished"
    diag_dir.mkdir(parents=True, exist_ok=True)

    summary_rows, summary_fieldnames = read_summary(summary_path)

    missing_key_cols = [c for c in KEY_COLS if c not in summary_fieldnames]
    if missing_key_cols:
        raise SystemExit(
            "[error] summary.tsv nemá očekávané sloupce:\n"
            + "\n".join(missing_key_cols)
        )

    expected_cols, nf_counter = detect_stability_expected_cols(stability_path)
    stability_rows, malformed, fixed_glued, padded = read_stability_fixed(
        stability_path,
        expected_cols,
    )

    stability_fixed_path = diag_dir / "stability_fixed.tsv"
    write_stability(stability_fixed_path, stability_rows)

    stability_count = Counter()
    for r in stability_rows:
        if r:
            stability_count[r[0]] += 1

    nonzero_counts = [v for v in stability_count.values() if v > 0]
    if nonzero_counts:
        expected_stability_rows = Counter(nonzero_counts).most_common(1)[0][0]
    else:
        expected_stability_rows = 0

    # Nejdřív základní status bez rozhodnutí o duplicitách.
    base = {}
    rows_by_run = {}

    for row in summary_rows:
        rid = row["run_id"]
        rows_by_run[rid] = row

        issues = []

        nstab = stability_count.get(rid, 0)
        if expected_stability_rows > 0 and nstab != expected_stability_rows:
            issues.append(f"stability_rows={nstab}, expected={expected_stability_rows}")

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
                issues.append("missing OutputSpectra ROOT file")
            else:
                root_file = str(rf)

        base[rid] = {
            "run_id": rid,
            "summary_line": row.get("_summary_line", ""),
            "base_ok": len(issues) == 0,
            "issues": "; ".join(issues),
            "stability_rows": nstab,
            "expected_stability_rows": expected_stability_rows,
            "root_log": log_path,
            "root_file": root_file,
            "key": config_key(row, KEY_COLS),
        }

    # Skupiny konfigurací.
    groups = defaultdict(list)
    first_seen = {}
    for row in summary_rows:
        rid = row["run_id"]
        key = config_key(row, KEY_COLS)
        groups[key].append(rid)
        first_seen.setdefault(key, rid)

    chosen_by_key = {}
    unresolved_keys = []

    for key, runs in groups.items():
        ok_runs = [r for r in runs if base[r]["base_ok"]]
        if ok_runs:
            chosen_by_key[key] = sorted(ok_runs, key=run_num)[0]
        else:
            chosen_by_key[key] = ""
            unresolved_keys.append(key)

    keep_runs = set(r for r in chosen_by_key.values() if r)

    # Status tabulka.
    status_rows = []
    duplicate_rows = []
    unresolved_rows = []
    config_rows = []

    for key_index, (key, runs) in enumerate(groups.items(), start=1):
        chosen = chosen_by_key[key]
        ok_runs = [r for r in runs if base[r]["base_ok"]]
        bad_runs = [r for r in runs if not base[r]["base_ok"]]
        first_run = first_seen[key]

        cfg_row = {
            "config_index": key_index,
            "chosen_run_id": chosen,
            "first_run_id": first_run,
            "n_runs": len(runs),
            "n_ok": len(ok_runs),
            "n_bad": len(bad_runs),
            "all_run_ids": ",".join(runs),
        }
        for c, v in zip(KEY_COLS, key):
            cfg_row[c] = v
        config_rows.append(cfg_row)

        if not chosen:
            unresolved_rows.append(cfg_row)

    for row in summary_rows:
        rid = row["run_id"]
        key = config_key(row, KEY_COLS)
        chosen = chosen_by_key[key]
        first_run = first_seen[key]
        b = base[rid]

        if chosen == rid:
            status = "KEEP_OK"
        elif b["base_ok"] and chosen:
            status = "DUPLICATE_OK"
        elif (not b["base_ok"]) and chosen:
            status = "BAD_IGNORED_HAS_OK_DUPLICATE"
        else:
            status = "BAD_UNRESOLVED"

        out = {
            "summary_line": b["summary_line"],
            "run_id": rid,
            "status": status,
            "chosen_for_config": chosen,
            "first_run_for_config": first_run,
            "stability_rows": b["stability_rows"],
            "expected_stability_rows": b["expected_stability_rows"],
            "issues": b["issues"],
            "root_log": b["root_log"],
            "root_file": b["root_file"],
        }
        for c in KEY_COLS:
            out[c] = row.get(c, "")
        status_rows.append(out)

        if status in ("DUPLICATE_OK", "BAD_IGNORED_HAS_OK_DUPLICATE"):
            duplicate_rows.append(out)

    status_fieldnames = [
        "summary_line",
        "run_id",
        "status",
        "chosen_for_config",
        "first_run_for_config",
        "stability_rows",
        "expected_stability_rows",
        "issues",
    ] + KEY_COLS + [
        "root_log",
        "root_file",
    ]

    write_tsv_dict(diag_dir / "run_status.tsv", status_rows, status_fieldnames)
    write_tsv_dict(diag_dir / "duplicate_runs.tsv", duplicate_rows, status_fieldnames)

    config_fieldnames = [
        "config_index",
        "chosen_run_id",
        "first_run_id",
        "n_runs",
        "n_ok",
        "n_bad",
    ] + KEY_COLS + [
        "all_run_ids",
    ]

    write_tsv_dict(diag_dir / "config_summary.tsv", config_rows, config_fieldnames)
    write_tsv_dict(diag_dir / "unresolved_configs.tsv", unresolved_rows, config_fieldnames)

    # Čisté summary: jeden OK run pro každou konfiguraci.
    clean_summary_rows = []
    for cfg in config_rows:
        rid = cfg["chosen_run_id"]
        if rid:
            clean_summary_rows.append(rows_by_run[rid])

    clean_summary_rows.sort(key=lambda r: int(config_rows[
        list(chosen_by_key.values()).index(r["run_id"])
    ]["config_index"]) if r["run_id"] in chosen_by_key.values() else run_num(r["run_id"]))

    # Bez interního sloupce _summary_line.
    clean_fieldnames = [c for c in summary_fieldnames if c != "_summary_line"]
    clean_summary_out = []
    for r in clean_summary_rows:
        clean_summary_out.append({c: r.get(c, "") for c in clean_fieldnames})

    write_tsv_dict(diag_dir / "summary_clean.tsv", clean_summary_out, clean_fieldnames)

    # Čisté stability: jen keep runy.
    clean_stability_rows = [r for r in stability_rows if r and r[0] in keep_runs]
    write_stability(diag_dir / "stability_clean.tsv", clean_stability_rows)

    with open(diag_dir / "kept_runs.list", "w") as f:
        for r in sorted(keep_runs, key=run_num):
            f.write(r + "\n")

    status_counter = Counter(r["status"] for r in status_rows)

    print("")
    print("========== Finished scan diagnostics ==========")
    print(f"scan_dir:                     {scan_dir}")
    print(f"summary rows:                 {len(summary_rows)}")
    print(f"unique bin configurations:    {len(groups)}")
    print(f"expected total configurations:{args.expected_total}")
    print(f"chosen usable configurations: {len(keep_runs)}")
    print(f"unresolved configurations:    {len(unresolved_rows)}")
    print(f"stability NF counts:          {dict(sorted(nf_counter.items()))}")
    print(f"expected stability columns:   {expected_cols}")
    print(f"expected stability rows/run:  {expected_stability_rows}")
    print(f"fixed glued stability lines:  {fixed_glued}")
    print(f"padded 9->10 stability lines: {padded}")
    print(f"malformed stability lines:    {len(malformed)}")
    print("")
    print("Status counts:")
    for k, v in sorted(status_counter.items()):
        print(f"  {k:30s} {v}")

    print("")
    if len(groups) == args.expected_total and len(keep_runs) == args.expected_total and len(unresolved_rows) == 0:
        print("[OK] Máš kompletní unikátní scan. Pro browser použij clean soubory.")
    else:
        print("[warning] Počet konfigurací nesedí přesně na očekávání, koukni do diagnostics_finished/*.tsv.")

    if malformed:
        malformed_path = diag_dir / "malformed_stability_lines.tsv"
        write_tsv_dict(malformed_path, malformed, ["line", "nf", "text"])
        print(f"[warning] Malformed stability lines written to: {malformed_path}")

    print("")
    print(f"Clean summary:      {diag_dir / 'summary_clean.tsv'}")
    print(f"Clean stability:    {diag_dir / 'stability_clean.tsv'}")
    print(f"Full run status:    {diag_dir / 'run_status.tsv'}")
    print(f"Duplicate runs:     {diag_dir / 'duplicate_runs.tsv'}")
    print(f"Config summary:     {diag_dir / 'config_summary.tsv'}")
    print(f"Unresolved configs: {diag_dir / 'unresolved_configs.tsv'}")
    print("===============================================")
    print("")


if __name__ == "__main__":
    main()
