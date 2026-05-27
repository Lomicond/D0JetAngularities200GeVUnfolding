#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import argparse
from pathlib import Path


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


def read_tsv_dict(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_lines(path, lines):
    with open(path, "w") as f:
        for line in lines:
            f.write(line.rstrip("\n") + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scan_dir", help="Např. Systematics/Binning/scanJetPt")
    ap.add_argument("original_script", help="Např. Systematics/Binning/BinningJetPt.sh")
    args = ap.parse_args()

    scan_dir = Path(args.scan_dir).resolve()
    original_script = Path(args.original_script).resolve()

    diag_dir = scan_dir / "diagnostics"
    run_status_path = diag_dir / "run_status.tsv"

    if not run_status_path.exists():
        raise SystemExit(f"[error] Missing {run_status_path}. Nejdřív spusť diagnose_binning_scan.py.")
    if not original_script.exists():
        raise SystemExit(f"[error] Missing {original_script}")

    rows = read_tsv_dict(run_status_path)

    done_keys = []
    kept_runs = []

    for r in rows:
        # Chceme jen první dobrý výskyt konfigurace.
        # DUPLICATE_OK nechceme dávat do hotových výstupů, protože už existuje originální OK run.
        if r.get("status") == "OK" and not r.get("duplicate_of"):
            key = "\t".join(str(r.get(c, "")).strip() for c in KEY_COLS)
            done_keys.append(key)
            kept_runs.append(r.get("run_id", ""))

    done_key_path = diag_dir / "done_config_keys.tsv"
    write_lines(done_key_path, done_keys)

    kept_runs_path = diag_dir / "kept_good_runs.list"
    write_lines(kept_runs_path, kept_runs)

    text = original_script.read_text()

    if "DONE_KEYS_FILE" in text or "CONFIG_KEY" in text:
        raise SystemExit(
            "[error] The script already seems to contain resume/skip logic. "
            "Make a clean copy of the original BinningJetPt.sh first."
        )

    marker1 = 'SUMMARY="${SCAN_DIR}/summary.tsv"\n'
    insert1 = r'''
# -------------------------
# Resume support
# -------------------------
# This file is created by prepare_resume_binning_scan.py.
# It contains only configurations that are already complete and should be skipped.
DONE_KEYS_FILE="${SCAN_DIR}/diagnostics/done_config_keys.tsv"
SKIP_DONE_CONFIGS=1

declare -A DONE_CONFIG_KEYS=()
if (( SKIP_DONE_CONFIGS )) && [[ -f "${DONE_KEYS_FILE}" ]]; then
  while IFS= read -r key; do
    [[ -n "${key}" ]] && DONE_CONFIG_KEYS["${key}"]=1
  done < "${DONE_KEYS_FILE}"
  echo "[resume] Loaded ${#DONE_CONFIG_KEYS[@]} completed bin configurations from ${DONE_KEYS_FILE}"
else
  echo "[resume] No completed-configuration list found, nothing will be skipped."
fi

'''
    if marker1 not in text:
        raise SystemExit("[error] Could not find insertion point after SUMMARY definition.")

    text = text.replace(marker1, marker1 + insert1, 1)

    marker2 = '''            # The true level must not have more bins than the reco level.
            if (( TRUE_N_BINS > RECO_N_BINS )); then
'''
    insert2 = r'''            # Resume: skip configurations that were already completed successfully.
            printf -v CONFIG_KEY "%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s" \
              "${RECO_PT_MIN}" "${RECO_PT_MAX}" "${RECO_N_BINS}" \
              "${TRUE_PT_MIN}" "${TRUE_PT_MAX}" "${TRUE_N_BINS}" \
              "${MIN_WIDTH}" "${MAX_WIDTH}" "${STEP}" "${WIDTH_TREND}" \
              "${RECO_EDGES}" "${TRUE_EDGES_STR}"

            if (( SKIP_DONE_CONFIGS )) && [[ -n "${DONE_CONFIG_KEYS[${CONFIG_KEY}]+x}" ]]; then
              echo "[resume skip] already complete: reco ${RECO_PT_MIN}-${RECO_PT_MAX}, N=${RECO_N_BINS} | true ${TRUE_PT_MIN}-${TRUE_PT_MAX}, N=${TRUE_N_BINS}"
              continue
            fi

'''
    if marker2 not in text:
        raise SystemExit("[error] Could not find insertion point before true/reco bin check.")

    text = text.replace(marker2, insert2 + marker2, 1)

    out_script = original_script.with_name(original_script.stem + "_resume.sh")
    out_script.write_text(text)
    out_script.chmod(0o755)

    print("")
    print("========== Resume preparation ==========")
    print(f"Good unique completed runs: {len(kept_runs)}")
    print(f"Done keys written to:       {done_key_path}")
    print(f"Good run list written to:   {kept_runs_path}")
    print(f"Resume script written to:   {out_script}")
    print("")
    print("Now run:")
    print(f"  bash {out_script}")
    print("========================================")
    print("")


if __name__ == "__main__":
    main()
