#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os
import re
import shutil
from pathlib import Path


def run_num(run_id):
    m = re.search(r"(\d+)$", str(run_id))
    return int(m.group(1)) if m else -1


def new_run_id(i):
    return f"r{i:06d}"


def copy_file(src, dst, hardlink=False):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if hardlink:
        os.link(src, dst)
    else:
        shutil.copy2(src, dst)


def copy_tree(src, dst, hardlink=False):
    if not src.exists():
        return False

    if dst.exists():
        raise RuntimeError(f"Target already exists: {dst}")

    copy_fun = os.link if hardlink else shutil.copy2
    shutil.copytree(src, dst, copy_function=copy_fun)
    return True


def read_summary(path):
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        cols = reader.fieldnames
    return rows, cols


def write_summary(path, rows, cols):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=cols, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def rewrite_stability(src_path, dst_path, run_map):
    n_in = 0
    n_out = 0
    unknown = set()

    with src_path.open("r", errors="replace") as f, dst_path.open("w", newline="") as g:
        writer = csv.writer(g, delimiter="\t", lineterminator="\n")

        for line in f:
            n_in += 1
            fields = line.rstrip("\n").split("\t")
            if not fields:
                continue

            old = fields[0]
            if old not in run_map:
                unknown.add(old)
                continue

            fields[0] = run_map[old]
            writer.writerow(fields)
            n_out += 1

    return n_in, n_out, unknown


def find_output_files(src_scan, old_run):
    found = []

    # ROOT outputs from Machine.C are normally here.
    for sub in ["Output", "OutputPdf"]:
        d = src_scan / sub
        if d.exists():
            found.extend(sorted(p for p in d.glob(f"*{old_run}*") if p.is_file()))

    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "scan_dir",
        help="Původní složka, např. Systematics/Binning/scanJetPt",
    )
    ap.add_argument(
        "--clean-input-dir",
        default="clean_from_stability",
        help="Podsložka s summary_clean.tsv a stability_clean.tsv",
    )
    ap.add_argument(
        "--out-dir",
        default="",
        help="Výstupní složka. Default: vedle původní složky jako scanJetPt_clean",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Smazat existující výstupní složku.",
    )
    ap.add_argument(
        "--hardlink",
        action="store_true",
        help="Použít hardlinky místo kopírování souborů. Šetří místo, funguje jen na stejném filesystému.",
    )
    ap.add_argument(
        "--strict",
        action="store_true",
        help="Skončit chybou, pokud chybí run dir, override nebo ROOT output.",
    )
    args = ap.parse_args()

    src_scan = Path(args.scan_dir).resolve()
    clean_in = src_scan / args.clean_input_dir

    src_summary = clean_in / "summary_clean.tsv"
    src_stability = clean_in / "stability_clean.tsv"

    if not src_summary.exists():
        raise SystemExit(f"[error] Missing {src_summary}")
    if not src_stability.exists():
        raise SystemExit(f"[error] Missing {src_stability}")

    if args.out_dir:
        dst_scan = Path(args.out_dir).resolve()
    else:
        dst_scan = src_scan.with_name(src_scan.name + "_clean")

    if dst_scan.exists():
        if args.force:
            shutil.rmtree(dst_scan)
        else:
            raise SystemExit(
                f"[error] Output directory already exists:\n  {dst_scan}\n"
                "Use --force if you want to recreate it."
            )

    dst_scan.mkdir(parents=True)

    dst_runs = dst_scan / "runs"
    dst_ovr = dst_scan / "overrides"
    dst_out = dst_scan / "Output"
    dst_pdf = dst_scan / "OutputPdf"

    dst_runs.mkdir()
    dst_ovr.mkdir()
    dst_out.mkdir()
    dst_pdf.mkdir()

    rows, cols = read_summary(src_summary)
    rows = sorted(rows, key=lambda r: run_num(r["run_id"]))

    run_map = {}
    map_rows = []

    for i, row in enumerate(rows, start=1):
        old = row["run_id"]
        new = new_run_id(i)
        run_map[old] = new
        map_rows.append({"old_run_id": old, "new_run_id": new})

    warnings = []

    new_summary_rows = []

    for row in rows:
        old = row["run_id"]
        new = run_map[old]

        new_row = dict(row)
        new_row["run_id"] = new

        # Rewrite paths in summary.
        if "override_file" in new_row:
            new_row["override_file"] = str(dst_ovr / f"override_{new}.C")

        if "root_log" in new_row:
            new_row["root_log"] = str(dst_runs / new / "root.log")

        new_summary_rows.append(new_row)

        # Copy run directory.
        old_run_dir = src_scan / "runs" / old
        new_run_dir = dst_runs / new

        if old_run_dir.exists():
            copy_tree(old_run_dir, new_run_dir, hardlink=args.hardlink)
        else:
            msg = f"missing run dir for {old}: {old_run_dir}"
            warnings.append(msg)
            if args.strict:
                raise SystemExit("[error] " + msg)

        # Copy override macro.
        old_override_candidates = []

        if row.get("override_file"):
            old_override_candidates.append(Path(row["override_file"]))

        old_override_candidates.append(src_scan / "overrides" / f"override_{old}.C")

        old_override = None
        for cand in old_override_candidates:
            if cand.exists():
                old_override = cand
                break

        if old_override is not None:
            copy_file(old_override, dst_ovr / f"override_{new}.C", hardlink=args.hardlink)
        else:
            msg = f"missing override for {old}"
            warnings.append(msg)
            if args.strict:
                raise SystemExit("[error] " + msg)

        # Copy ROOT/PDF outputs containing the run ID and rename them.
        out_files = find_output_files(src_scan, old)

        if not out_files:
            msg = f"no Output/OutputPdf files found for {old}"
            warnings.append(msg)
            if args.strict:
                raise SystemExit("[error] " + msg)

        for src_file in out_files:
            rel_subdir = src_file.parent.name
            dst_parent = dst_scan / rel_subdir
            new_name = src_file.name.replace(old, new)
            copy_file(src_file, dst_parent / new_name, hardlink=args.hardlink)

    # Write clean summary.
    write_summary(dst_scan / "summary.tsv", new_summary_rows, cols)

    # Rewrite stability with new run IDs.
    n_in, n_out, unknown = rewrite_stability(src_stability, dst_scan / "stability.tsv", run_map)

    if unknown:
        msg = f"stability rows with unknown run_id: {sorted(list(unknown))[:10]}"
        warnings.append(msg)
        if args.strict:
            raise SystemExit("[error] " + msg)

    # Write map old->new.
    with (dst_scan / "run_id_map.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=["old_run_id", "new_run_id"], lineterminator="\n")
        w.writeheader()
        w.writerows(map_rows)

    # Copy kept list in both old and new form.
    with (dst_scan / "kept_runs_old.list").open("w") as f:
        for r in rows:
            f.write(r["run_id"] + "\n")

    with (dst_scan / "kept_runs.list").open("w") as f:
        for i in range(1, len(rows) + 1):
            f.write(new_run_id(i) + "\n")

    # Small README.
    with (dst_scan / "README_clean_scan.txt").open("w") as f:
        f.write(
            "This folder was created from clean_from_stability/summary_clean.tsv "
            "and clean_from_stability/stability_clean.tsv.\n"
        )
        f.write("Run IDs were remapped to a clean continuous sequence.\n")
        f.write("See run_id_map.tsv for old -> new mapping.\n")

    print("")
    print("========== Clean scan folder created ==========")
    print(f"source scan:        {src_scan}")
    print(f"clean input:        {clean_in}")
    print(f"output scan:        {dst_scan}")
    print(f"runs mapped:        {len(run_map)}")
    print(f"stability rows in:  {n_in}")
    print(f"stability rows out: {n_out}")
    print(f"copy mode:          {'hardlink' if args.hardlink else 'copy'}")
    print("")
    print(f"New summary:        {dst_scan / 'summary.tsv'}")
    print(f"New stability:      {dst_scan / 'stability.tsv'}")
    print(f"Run map:            {dst_scan / 'run_id_map.tsv'}")

    if warnings:
        print("")
        print(f"[warning] {len(warnings)} warnings:")
        for w in warnings[:20]:
            print("  - " + w)
        if len(warnings) > 20:
            print(f"  ... plus {len(warnings) - 20} more")
    else:
        print("")
        print("[OK] No warnings.")

    print("==============================================")
    print("")


if __name__ == "__main__":
    main()
