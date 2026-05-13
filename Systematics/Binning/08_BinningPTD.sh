#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Binning systematic scan: pTD
#
# Run from the main project directory, i.e. the directory containing:
#   config.h, config_hist.h, Unfolding/Machine.C, Data/...
#
# Recommended run command:
#   bash Systematics/Binning/08_BinningPTD.sh
# ============================================================

# -------------------------
# Working-directory checks
# -------------------------
PROJECT_DIR="$(pwd -P)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

[[ -f "${PROJECT_DIR}/config.h" ]]       || { echo "[error] Run this script from the main project directory: missing ./config.h"; exit 1; }
[[ -f "${PROJECT_DIR}/config_hist.h" ]]  || { echo "[error] Run this script from the main project directory: missing ./config_hist.h"; exit 1; }
[[ -f "${PROJECT_DIR}/Unfolding/Machine.C" ]] || { echo "[error] Missing ./Unfolding/Machine.C"; exit 1; }

command -v root    >/dev/null 2>&1 || { echo "[error] root is not in PATH"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "[error] python3 is not in PATH"; exit 1; }

# =========================
# User settings
# =========================
MACHINE_MACRO="${PROJECT_DIR}/Unfolding/Machine.C"
INPUT_FILE="${PROJECT_DIR}/Data/Output_real_final_01022026.root"

[[ -f "${INPUT_FILE}" ]] || { echo "[error] Missing input file: ${INPUT_FILE}"; exit 1; }

# -------------------------
# pTD scan only
# -------------------------
# Index 5 = p_T^D / pTD.
# The override below changes only angRecoBinsVec[*][5]
# and angMcBinsVecCustom[*][5].
# All other binnings remain exactly as defined in config.h/config_hist.h.
PTD_INDEX=5

PTD_TRUE_START_LIST=(0 0.3 0.5)
PTD_RECO_START_LIST=(0 0.3 0.5)
PTD_RECO_ALLOWED_WIDTHS=(0.1 0.15 0.2 0.3)
PTD_RECO_MIN_WIDTH=0.1
PTD_RIGHT_EDGE=1.01
PTD_LAST_BIN_LEFT=0.9   # The reco last bin is kept fixed as 0.9 -> 1.01.

# Edge nudging.
# According to the agreed setting for pTD, use a coherent shift of 0.1.
# The shift is applied only to middle internal reco edges.
# The first edge, the last edge, and the edge before the last edge are kept fixed.
# This keeps the last bin stable and avoids pathological tiny last bins.
PTD_EDGE_SHIFTS=(0 -0.01 0.01)

# -------------------------
# Machine() parameters
# Signature in the current Machine.C:
# Machine(fonllJet, CutOfNegative, minJetPtRecoCut, savedIter,
#         InputFile, OutputFile, minPtD0Cut, maxPtD0Cut,
#         OverrideMacro, ScanDir, usePriorShapeWeighting, systematicSPlot)
# -------------------------
FONLL_JET=1
CUT_NEG=1
MIN_JET_PT_RECO_CUT=-30
SAVED_ITER=4
MIN_PT_D0=1
MAX_PT_D0=10
USE_PRIOR_SHAPE_WEIGHTING=0
SYSTEMATIC_SPLOT=0

# Safety limit for the whole scan. 0 = no limit.
MAX_TOTAL_RUNS=0

# -------------------------
# Output location
# -------------------------
SCAN_DIR="${SCRIPT_DIR}/scanPTD"
OVR_DIR="${SCAN_DIR}/overrides"
RUN_DIR="${SCAN_DIR}/runs"
SUMMARY="${SCAN_DIR}/summary.tsv"
STABILITY="${SCAN_DIR}/stability.tsv"

mkdir -p "${OVR_DIR}" "${RUN_DIR}" "${PROJECT_DIR}/OutputPdf"

# =========================
# Helper functions
# =========================
run_in_file() {
  local run_id="$1"
  local file="$2"
  [[ -f "${file}" ]] || return 1
  grep -q "^${run_id}[[:space:]]" "${file}"
}

run_is_complete() {
  local run_id="$1"
  run_in_file "${run_id}" "${SUMMARY}" && run_in_file "${run_id}" "${STABILITY}"
}

remove_run_from_file() {
  local run_id="$1"
  local file="$2"
  [[ -f "${file}" ]] || return 0
  local tmp="${file}.tmp.$$"
  grep -v "^${run_id}[[:space:]]" "${file}" > "${tmp}" || true
  mv "${tmp}" "${file}"
}

format_true_edges_from_start() {
  local start="$1"
  python3 - "${start}" <<'PY'
import sys
start = float(sys.argv[1])
base = [0, 0.3, 0.5, 0.65, 0.75, 0.85, 1.01]
edges = [x for x in base if x >= start - 1e-9]

def fmt(x):
    s = f"{x:.10f}".rstrip('0').rstrip('.')
    return s if s else '0'

print(", ".join(fmt(x) for x in edges))
PY
}

count_edges_bins() {
  local edges_in="$1"
  awk -F',' '{print NF-1}' <<< "${edges_in}"
}

shift_edges() {
  local edges_in="$1"
  local shift="$2"
  python3 - "${edges_in}" "${shift}" <<'PY'
import sys

edges = [float(x.strip()) for x in sys.argv[1].split(',') if x.strip()]
shift = float(sys.argv[2])

if len(edges) < 2:
    raise SystemExit("Need at least two edges")

# Shift only middle internal edges.
# Keep fixed:
#   i = 0              first edge of the selected range
#   i = len(edges)-2   edge before the last edge, preserving the last bin
#   i = len(edges)-1   right edge
if abs(shift) > 0 and len(edges) > 3:
    for i in range(1, len(edges) - 2):
        edges[i] = round(edges[i] + shift, 10)

# The edges must remain strictly increasing.
for a, b in zip(edges, edges[1:]):
    if not (b > a):
        raise SystemExit(f"Invalid shifted edges: {edges}")

if edges[0] < -1e-9:
    raise SystemExit(f"Invalid shifted edges: {edges}")

if edges[-1] > 1.01 + 1e-9:
    raise SystemExit(f"Invalid shifted edges: {edges}")

def fmt(x):
    s = f"{x:.10f}".rstrip('0').rstrip('.')
    return s if s else '0'

print(", ".join(fmt(x) for x in edges))
PY
}

write_ptd_override() {
  local file="$1"
  local reco_edges="$2"
  local true_edges="$3"

  cat > "${file}" <<EOF_OVR
{
  cout << "[override] pTD reco edges: ${reco_edges}" << endl;
  cout << "[override] pTD true edges: ${true_edges}" << endl;

  for (int ic = 0; ic < nCentralityBins; ++ic) {
    angRecoBinsVec[ic][${PTD_INDEX}].clear();
  }
  double reco_edges[] = { ${reco_edges} };
  int nReco = (int)(sizeof(reco_edges)/sizeof(double));
  for (int ic = 0; ic < nCentralityBins; ++ic) {
    for (int i = 0; i < nReco; ++i) {
      angRecoBinsVec[ic][${PTD_INDEX}].push_back(reco_edges[i]);
    }
  }

  for (int ic = 0; ic < nCentralityBins; ++ic) {
    angMcBinsVecCustom[ic][${PTD_INDEX}].clear();
  }
  double true_edges[] = { ${true_edges} };
  int nTrue = (int)(sizeof(true_edges)/sizeof(double));
  for (int ic = 0; ic < nCentralityBins; ++ic) {
    for (int i = 0; i < nTrue; ++i) {
      angMcBinsVecCustom[ic][${PTD_INDEX}].push_back(true_edges[i]);
    }
  }
}
EOF_OVR
}

# =========================
# Summary header
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\tpTD_edge_shift\tpTD_true_start\tpTD_true_n_bins\tpTD_reco_start\tpTD_reco_n_bins\tpTD_reco_edges\tpTD_true_edges\toverride_file\troot_log\n" > "${SUMMARY}"
fi

# =========================
# Determine the last run ID
# =========================
run_counter=0
shopt -s nullglob
for d in "${RUN_DIR}"/r*; do
  [[ -d "${d}" ]] || continue
  bn=$(basename "${d}")
  if [[ "${bn}" =~ ^r([0-9]+)$ ]]; then
    num=$((10#${BASH_REMATCH[1]}))
    (( num > run_counter )) && run_counter=${num}
  fi
done
shopt -u nullglob

runs_started=0

# Use forced ACLiC rebuild only for the first real Machine() call.
# After that, use the already compiled library.
ACLIC_SUFFIX="++"

# =========================
# Main loop
# =========================
for PTD_TRUE_START in "${PTD_TRUE_START_LIST[@]}"; do
  PTD_TRUE_EDGES_BASE=$(format_true_edges_from_start "${PTD_TRUE_START}")

  for PTD_RECO_START in "${PTD_RECO_START_LIST[@]}"; do
    # Reco must be at least as wide as truth: reco_start <= true_start.
    if ! python3 - "${PTD_RECO_START}" "${PTD_TRUE_START}" <<'PY' >/dev/null
import sys
reco = float(sys.argv[1])
tru  = float(sys.argv[2])
raise SystemExit(0 if reco <= tru + 1e-9 else 1)
PY
    then
      continue
    fi

    while IFS= read -r PTD_RECO_EDGES_BASE; do
      [[ -n "${PTD_RECO_EDGES_BASE}" ]] || continue

      for PTD_EDGE_SHIFT in "${PTD_EDGE_SHIFTS[@]}"; do
        if ! PTD_RECO_EDGES=$(shift_edges "${PTD_RECO_EDGES_BASE}" "${PTD_EDGE_SHIFT}"); then
          echo "[skip] pTD shift=${PTD_EDGE_SHIFT}, true start=${PTD_TRUE_START}, reco start=${PTD_RECO_START} -> invalid reco edges"
          continue
        fi

	PTD_TRUE_EDGES="${PTD_TRUE_EDGES_BASE}"
        PTD_RECO_N_BINS=$(count_edges_bins "${PTD_RECO_EDGES}")
        PTD_TRUE_N_BINS=$(count_edges_bins "${PTD_TRUE_EDGES}")

        if (( PTD_RECO_N_BINS < PTD_TRUE_N_BINS )); then
          echo "[skip] pTD shift=${PTD_EDGE_SHIFT}, true start=${PTD_TRUE_START}, reco start=${PTD_RECO_START} -> Nreco=${PTD_RECO_N_BINS} < Ntrue=${PTD_TRUE_N_BINS}"
          continue
        fi

        ((run_counter += 1))
        printf -v RUN_ID "r%06d" "${run_counter}"

        if run_is_complete "${RUN_ID}"; then
          echo "[resume] ${RUN_ID} is already complete -> skip"
          continue
        fi

        if (( MAX_TOTAL_RUNS > 0 && runs_started >= MAX_TOTAL_RUNS )); then
          echo "Reached MAX_TOTAL_RUNS=${MAX_TOTAL_RUNS}, stopping."
          exit 0
        fi

        remove_run_from_file "${RUN_ID}" "${SUMMARY}"
        remove_run_from_file "${RUN_ID}" "${STABILITY}"

        ((runs_started += 1))

        OVR_FILE="${OVR_DIR}/override_${RUN_ID}.C"
        OUT_DIR="${RUN_DIR}/${RUN_ID}"
        ROOT_LOG="${OUT_DIR}/root.log"
        rm -rf "${OUT_DIR}"
        mkdir -p "${OUT_DIR}"

        write_ptd_override "${OVR_FILE}" "${PTD_RECO_EDGES}" "${PTD_TRUE_EDGES}"

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
          "${RUN_ID}" \
          "${PTD_EDGE_SHIFT}" \
          "${PTD_TRUE_START}" "${PTD_TRUE_N_BINS}" \
          "${PTD_RECO_START}" "${PTD_RECO_N_BINS}" \
          "${PTD_RECO_EDGES}" "${PTD_TRUE_EDGES}" \
          "${OVR_FILE}" "${ROOT_LOG}" \
          >> "${SUMMARY}"

        echo "[${RUN_ID}] pTD shift=${PTD_EDGE_SHIFT}, true start=${PTD_TRUE_START}, reco start=${PTD_RECO_START}, Nreco=${PTD_RECO_N_BINS}, Ntrue=${PTD_TRUE_N_BINS}"

        # OutputFile = RUN_ID because Machine.C uses runId = BaseName(OutputFile).
        # ScanDir    = SCAN_DIR so stability.tsv and Output/OutputSpectra*.root go there.
        if ! root -l -b -q \
          "${MACHINE_MACRO}${ACLIC_SUFFIX}( ${FONLL_JET}, ${CUT_NEG}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${RUN_ID}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\", \"${SCAN_DIR}\", ${USE_PRIOR_SHAPE_WEIGHTING}, ${SYSTEMATIC_SPLOT} )" \
          > "${ROOT_LOG}" 2>&1; then
          echo "[error] ROOT failed for ${RUN_ID}. Last 60 lines of the log:"
          tail -n 60 "${ROOT_LOG}" || true
          exit 1
        fi

        # Only the first run should force recompilation.
        ACLIC_SUFFIX="+"

      done
    done < <(
      python3 - "${PTD_RECO_START}" "${PTD_RIGHT_EDGE}" "${PTD_LAST_BIN_LEFT}" "${PTD_RECO_MIN_WIDTH}" "${PTD_RECO_ALLOWED_WIDTHS[@]}" <<'PY'
import sys

start = float(sys.argv[1])
right_edge = float(sys.argv[2])
last_left = float(sys.argv[3])
min_width = float(sys.argv[4])
allowed = [float(x) for x in sys.argv[5:]]
EPS = 1e-9

if start > last_left - EPS:
    raise SystemExit(0)

results = set()

def rec(prev_width_effective, acc_edges_desc):
    current_left = acc_edges_desc[-1]
    rem = current_left - start

    if rem >= prev_width_effective - EPS and rem >= min_width - EPS:
        edges_desc = acc_edges_desc + [start]
        edges = list(reversed(edges_desc))
        if abs(edges[-1] - right_edge) < 1e-8 and abs(edges[0] - start) < 1e-8:
            results.add(tuple(round(x, 10) for x in edges))

    for w in allowed:
        if w + EPS < prev_width_effective:
            continue
        new_left = current_left - w
        if new_left <= start + EPS:
            continue
        rec(w, acc_edges_desc + [new_left])

rec(min_width, [right_edge, last_left])

def fmt(x):
    s = f"{x:.10f}".rstrip('0').rstrip('.')
    return s if s else '0'

for edges in sorted(results):
    print(", ".join(fmt(x) for x in edges))
PY
    )

  done
done

echo "Done."
echo "Summary:   ${SUMMARY}"
echo "Stability: ${SCAN_DIR}/stability.tsv"
echo "Spectra:   ${SCAN_DIR}/Output/"
