#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Binning systematic scan: lambda_2^1 / l12
#
# Run from the main project directory, i.e. the directory containing:
#   config.h, config_hist.h, Unfolding/Machine.C, Data/...
#
# Recommended run command:
#   bash Systematics/Binning/05_BinningL12.sh
# ============================================================

# -------------------------
# Working-directory checks
# -------------------------
PROJECT_DIR="$(pwd -P)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

[[ -f "${PROJECT_DIR}/config.h" ]]            || { echo "[error] Run this script from the main project directory: missing ./config.h"; exit 1; }
[[ -f "${PROJECT_DIR}/config_hist.h" ]]       || { echo "[error] Run this script from the main project directory: missing ./config_hist.h"; exit 1; }
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
# lambda_2^1 / l12 scan only
# -------------------------
# Index 2 corresponds to lambda_2^1 / l12.
# The override macro changes only:
#   angRecoBinsVec[*][2]
#   angMcBinsVecCustom[*][2]
# All other binnings stay exactly as defined in config.h/config_hist.h.
L12_INDEX=2

L12_TRUE_MAX_LIST=(0.4 0.6 0.8 1.0)
L12_RECO_MAX_LIST=(0.4 0.6 0.8 1.0)
L12_DENSE_BIN_COUNTS=(4 5 6)
L12_MIN_WIDTH=0.025
L12_TAIL_WIDTHS=(0.1 0.2 0.3 0.4)
L12_TRUE_N_BINS=8

# Edge nudging.
# Use the same controlled edge shift as for the previous angularities: 0.005.
# The shift is applied coherently to internal reco and true edges.
# The first and last edges remain fixed.
L12_EDGE_SHIFTS=(0 -0.005 0.005)

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
SCAN_DIR="${SCRIPT_DIR}/scanL12"
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

format_edges_from_max() {
  local xmax="$1"
  python3 - "${xmax}" <<'PY'
import sys
xmax = float(sys.argv[1])
edges = [0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, xmax]

def fmt(x):
    s = f"{x:.10f}".rstrip('0').rstrip('.')
    return s if s else '0'

print(", ".join(fmt(x) for x in edges))
PY
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

# Shift internal edges only. Keep the first and last edges fixed.
if abs(shift) > 0:
    for i in range(1, len(edges) - 1):
        edges[i] = round(edges[i] + shift, 10)

# Safety check: edges must remain strictly increasing.
for a, b in zip(edges, edges[1:]):
    if not (b > a):
        raise SystemExit(f"Invalid shifted edges: {edges}")

def fmt(x):
    s = f"{x:.10f}".rstrip('0').rstrip('.')
    return s if s else '0'

print(", ".join(fmt(x) for x in edges))
PY
}

# =========================
# Summary header
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\tl12_edge_shift\tl12_reco_max\tl12_dense_bins\tl12_reco_n_bins\tl12_true_max\tl12_true_n_bins\tl12_reco_edges\tl12_true_edges\toverride_file\troot_log\n" > "${SUMMARY}"
fi

# =========================
# Resume counters
# =========================
run_counter=0
runs_started=0

# Use forced ACLiC rebuild only for the first real Machine() call.
# After that, use the already compiled library.
ACLIC_SUFFIX="++"

# =========================
# Main loop
# =========================
for L12_RECO_MAX in "${L12_RECO_MAX_LIST[@]}"; do
  for L12_DENSE_BINS in "${L12_DENSE_BIN_COUNTS[@]}"; do

    while IFS= read -r L12_RECO_EDGES_BASE; do
      [[ -n "${L12_RECO_EDGES_BASE}" ]] || continue

      for L12_EDGE_SHIFT in "${L12_EDGE_SHIFTS[@]}"; do
        L12_RECO_EDGES=$(shift_edges "${L12_RECO_EDGES_BASE}" "${L12_EDGE_SHIFT}")

        # Number of reco bins = number of edges - 1.
        L12_RECO_N_BINS=$(awk -F',' '{print NF-1}' <<< "${L12_RECO_EDGES}")

        for L12_TRUE_MAX in "${L12_TRUE_MAX_LIST[@]}"; do
          L12_TRUE_EDGES_BASE=$(format_edges_from_max "${L12_TRUE_MAX}")
          L12_TRUE_EDGES=$(shift_edges "${L12_TRUE_EDGES_BASE}" "${L12_EDGE_SHIFT}")
          L12_TRUE_N_BINS=${L12_TRUE_N_BINS}

          if (( L12_RECO_N_BINS < L12_TRUE_N_BINS )); then
            echo "[skip] l12 shift=${L12_EDGE_SHIFT}, reco max=${L12_RECO_MAX}, dense=${L12_DENSE_BINS} -> Nreco=${L12_RECO_N_BINS} < Ntrue=${L12_TRUE_N_BINS}"
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

          cat > "${OVR_FILE}" <<EOF_OVR
{
  cout << "[override] lambda_2^1 / l12 reco edges: ${L12_RECO_EDGES}" << endl;
  cout << "[override] lambda_2^1 / l12 true edges: ${L12_TRUE_EDGES}" << endl;

  // lambda_2^1 / l12 scan only, index ${L12_INDEX}.
  // All other binning variables stay exactly as defined in config.h/config_hist.h.

  double reco_edges[] = { ${L12_RECO_EDGES} };
  int nReco = (int)(sizeof(reco_edges)/sizeof(double));
  for (int ic = 0; ic < nCentralityBins; ++ic) {
    angRecoBinsVec[ic][${L12_INDEX}].clear();
    for (int i = 0; i < nReco; ++i) {
      angRecoBinsVec[ic][${L12_INDEX}].push_back(reco_edges[i]);
    }
  }

  double true_edges[] = { ${L12_TRUE_EDGES} };
  int nTrue = (int)(sizeof(true_edges)/sizeof(double));
  for (int ic = 0; ic < nCentralityBins; ++ic) {
    angMcBinsVecCustom[ic][${L12_INDEX}].clear();
    for (int i = 0; i < nTrue; ++i) {
      angMcBinsVecCustom[ic][${L12_INDEX}].push_back(true_edges[i]);
    }
  }
}
EOF_OVR

          echo "[${RUN_ID}] l12 shift=${L12_EDGE_SHIFT}, reco max=${L12_RECO_MAX}, dense=${L12_DENSE_BINS}, Nreco=${L12_RECO_N_BINS} | true max=${L12_TRUE_MAX}"

          if ! root -l -b -q \
            "${MACHINE_MACRO}${ACLIC_SUFFIX}( ${FONLL_JET}, ${CUT_NEG}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${RUN_ID}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\", \"${SCAN_DIR}\", ${USE_PRIOR_SHAPE_WEIGHTING}, ${SYSTEMATIC_SPLOT} )" \
            > "${ROOT_LOG}" 2>&1; then
            echo "[error] ROOT failed for ${RUN_ID}. Last 60 lines of the log:"
            tail -n 60 "${ROOT_LOG}" || true
            exit 1
          fi

          # Only the first run should force recompilation.
          ACLIC_SUFFIX="+"

          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${RUN_ID}" \
            "${L12_EDGE_SHIFT}" \
            "${L12_RECO_MAX}" "${L12_DENSE_BINS}" "${L12_RECO_N_BINS}" \
            "${L12_TRUE_MAX}" "${L12_TRUE_N_BINS}" \
            "${L12_RECO_EDGES}" "${L12_TRUE_EDGES}" \
            "${OVR_FILE}" "${ROOT_LOG}" \
            >> "${SUMMARY}"

        done
      done
    done < <(
      python3 - "${L12_RECO_MAX}" "${L12_DENSE_BINS}" "${L12_MIN_WIDTH}" "${L12_TRUE_N_BINS}" "${L12_TAIL_WIDTHS[@]}" <<'PY'
import sys

xmax = float(sys.argv[1])
dense_bins = int(sys.argv[2])
min_width = float(sys.argv[3])
min_n_bins = int(sys.argv[4])
allowed = [float(x) for x in sys.argv[5:]]
EPS = 1e-9

# Fixed dense prefix.
prefix = [0.0]
for _ in range(dense_bins):
    prefix.append(prefix[-1] + min_width)

fixed_sum = prefix[-1]
remaining = xmax - fixed_sum
if remaining < min(allowed) - EPS:
    raise SystemExit(0)

results = set()

def rec(prev, used, widths):
    # Tail bin count = len(widths)+1 because the final residual bin is still missing.
    tail_count_so_far = len(widths)
    total_bins_so_far = dense_bins + tail_count_so_far + 1  # +1 for final residual
    rem = remaining - used

    # Can we close now with one final residual bin?
    if total_bins_so_far >= min_n_bins:
        if rem + EPS >= max(prev, min(allowed)):
            tail = widths + [round(rem, 10)]
            edges = prefix[:]
            acc = fixed_sum
            for w in tail:
                acc += w
                edges.append(round(acc, 10))
            if abs(edges[-1] - xmax) < 1e-8:
                results.add(tuple(edges))

    # Or add another regular coarse bin from the allowed set.
    for w in allowed:
        if w + EPS < max(prev, min(allowed)):
            continue
        rem_after = rem - w
        if rem_after < max(w, min(allowed)) - EPS:
            continue
        rec(w, used + w, widths + [w])

rec(min(allowed), 0.0, [])

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
