#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Filtered jet-pT binning systematic scan for the D0-pT subrange
#
# Usage:
#   ./Systematics/Binning/01_BinningJetPtNew.sh 5 10
#
# Philosophy:
#   - keep the jet-pT binning style selected in the old large scan;
#   - use the three final reco-binning styles (one per centrality) as
#     reference styles, but—as in the old scan—apply one reco candidate
#     to all centralities in a given Machine.C run;
#   - for D0 pT >= 5 GeV/c, remove the unnecessary low-pT part of the
#     references.  The 10-40% reference intentionally keeps the 4.5 edge;
#   - keep the final true-level style 5,7,9,11,13,15,20;
#   - add coherent reco-edge "nudging" ("seveleni"): 0, -0.5, +0.5 GeV/c.
#     Only internal reco edges move; the first and last edges stay fixed.
#
# This gives 3 reco reference styles x 3 coherent shifts = 9 runs.
#
# Simple_splot.C is NOT run here.  This script assumes the nominal sPlot
# outputs already exist; Machine.C receives the requested D0-pT range.
# ============================================================

# -------------------------
# D0 pT range for Machine()
# -------------------------
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 MIN_PT_D0 MAX_PT_D0"
  echo "Example: $0 5 10"
  exit 1
fi

MIN_PT_D0="$1"
MAX_PT_D0="$2"

if ! [[ "${MIN_PT_D0}" =~ ^[0-9]+$ && "${MAX_PT_D0}" =~ ^[0-9]+$ ]]; then
  echo "[error] MIN_PT_D0 and MAX_PT_D0 must be integers."
  exit 1
fi

if (( MIN_PT_D0 >= MAX_PT_D0 )); then
  echo "[error] MIN_PT_D0 must be smaller than MAX_PT_D0."
  exit 1
fi

# The filtered candidate set below is intentionally constructed for
# the D0-pT lower cut of 5 GeV/c.
if (( MIN_PT_D0 != 5 )); then
  echo "[error] This filtered jet-pT candidate set is prepared for MIN_PT_D0=5."
  echo "        Requested MIN_PT_D0=${MIN_PT_D0}."
  exit 1
fi

# -------------------------
# Working-directory checks
# -------------------------
PROJECT_DIR="$(pwd -P)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

[[ -f "${PROJECT_DIR}/config.h" ]]              || { echo "[error] Run this script from the main project directory: missing ./config.h"; exit 1; }
[[ -f "${PROJECT_DIR}/config_hist.h" ]]         || { echo "[error] Run this script from the main project directory: missing ./config_hist.h"; exit 1; }
[[ -f "${PROJECT_DIR}/Unfolding/Machine.C" ]]   || { echo "[error] Missing ./Unfolding/Machine.C"; exit 1; }

command -v root    >/dev/null 2>&1 || { echo "[error] root is not in PATH"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "[error] python3 is not in PATH"; exit 1; }

# =========================
# User settings
# =========================
MACHINE_MACRO="${PROJECT_DIR}/Unfolding/Machine.C"
INPUT_FILE="${PROJECT_DIR}/Data/Output_real_final_01022026.root"

[[ -f "${INPUT_FILE}" ]] || { echo "[error] Missing input file: ${INPUT_FILE}"; exit 1; }

# -------------------------
# Filtered reco reference styles
# -------------------------
# These are the final reco styles from config.h, restricted to the useful
# low-pT region for D0 pT >= 5 GeV/c.
#
# Important:
#   * 0-10% and 40-80% naturally have an edge at 5.
#   * 10-40% keeps 4.5 as the first reco edge.  With D0 pT >= 5 GeV/c,
#     the 4.5-5 part is empty, and reco bins are not normalized by width
#     before unfolding.
#
# As in the old 01_BinningJetPt.sh, each candidate is applied to all
# centralities in one Machine.C run.  The labels only identify which
# final-centrality style the candidate came from.
RECO_REFERENCE_LABELS=(
  "ref_0_10"
  "ref_10_40"
  "ref_40_80"
)

RECO_REFERENCE_EDGES=(
  "5, 6, 7, 8, 9, 10, 13, 16, 25"
  "4.5, 7, 10, 13, 16, 19, 25"
  "5, 7, 9, 11, 13, 18, 25"
)

# Coherent reco-edge nudging, analogous to the angularity scans.
# Only internal edges are shifted.
RECO_EDGE_SHIFTS=(0 -0.5 0.5)

# -------------------------
# Filtered true-level style
# -------------------------
# This is the old final true binning with only the low-pT end removed.
# It is kept fixed during the reco "seveleni".
TRUE_EDGES="5, 7, 9, 11, 13, 15, 20"

# -------------------------
# Machine() parameters
# -------------------------
# Machine(fonllJet, CutOfNegative, minJetPtRecoCut, savedIter,
#         InputFile, OutputFile, minPtD0Cut, maxPtD0Cut,
#         OverrideMacro, ScanDir, usePriorShapeWeighting, systematicSPlot)
FONLL_JET=1
CUT_NEG=1
MIN_JET_PT_RECO_CUT=-30
SAVED_ITER=4
USE_PRIOR_SHAPE_WEIGHTING=0
SYSTEMATIC_SPLOT=0

# Safety limit for testing. 0 = no limit.
MAX_TOTAL_RUNS=0

# Fine jet-pT step used for the coherent edge shift.
EDGE_STEP=0.5

# -------------------------
# Output location
# -------------------------
SCAN_DIR="${SCRIPT_DIR}/scanJetPt_${MIN_PT_D0}_${MAX_PT_D0}"
OVR_DIR="${SCAN_DIR}/overrides"
RUN_DIR="${SCAN_DIR}/runs"
SUMMARY="${SCAN_DIR}/summary.tsv"
STABILITY="${SCAN_DIR}/stability.tsv"
PLAN="${SCAN_DIR}/planned_combinations.tsv"

mkdir -p \
  "${OVR_DIR}" \
  "${RUN_DIR}" \
  "${SCAN_DIR}/Output" \
  "${PROJECT_DIR}/OutputPdf"

# =========================
# Helper functions
# =========================
shift_edges() {
  local edges_in="$1"
  local shift="$2"

  python3 - "${edges_in}" "${shift}" <<'PY'
import sys

edges = [float(x.strip()) for x in sys.argv[1].split(",") if x.strip()]
shift = float(sys.argv[2])

if len(edges) < 2:
    raise SystemExit("Need at least two edges")

# Shift only internal reco edges; keep the first and last fixed.
if abs(shift) > 0:
    for i in range(1, len(edges) - 1):
        edges[i] = round(edges[i] + shift, 10)

for a, b in zip(edges, edges[1:]):
    if not b > a:
        raise SystemExit(f"Invalid shifted edges: {edges}")

def fmt(x):
    s = f"{x:.10f}".rstrip("0").rstrip(".")
    return s if s else "0"

print(", ".join(fmt(x) for x in edges))
PY
}

edge_stats() {
  local edges_in="$1"

  python3 - "${edges_in}" <<'PY'
import sys

edges = [float(x.strip()) for x in sys.argv[1].split(",") if x.strip()]
if len(edges) < 2:
    raise SystemExit("Need at least two edges")

widths = [b - a for a, b in zip(edges, edges[1:])]
if any(w <= 0 for w in widths):
    raise SystemExit(f"Non-increasing edges: {edges}")

def fmt(x):
    s = f"{x:.10f}".rstrip("0").rstrip(".")
    return s if s else "0"

print(fmt(edges[0]), fmt(edges[-1]), len(edges) - 1, fmt(min(widths)), fmt(max(widths)))
PY
}

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

# =========================
# Build a complete plan first
# =========================
printf "plan_index\treco_reference\tedge_shift\treco_edges_base\treco_edges\ttrue_edges\n" > "${PLAN}"

plan_index=0
for i in "${!RECO_REFERENCE_LABELS[@]}"; do
  RECO_REFERENCE="${RECO_REFERENCE_LABELS[$i]}"
  RECO_EDGES_BASE="${RECO_REFERENCE_EDGES[$i]}"

  for EDGE_SHIFT in "${RECO_EDGE_SHIFTS[@]}"; do
    RECO_EDGES="$(shift_edges "${RECO_EDGES_BASE}" "${EDGE_SHIFT}")"

    ((plan_index += 1))
    printf "%d\t%s\t%s\t%s\t%s\t%s\n" \
      "${plan_index}" \
      "${RECO_REFERENCE}" \
      "${EDGE_SHIFT}" \
      "${RECO_EDGES_BASE}" \
      "${RECO_EDGES}" \
      "${TRUE_EDGES}" \
      >> "${PLAN}"
  done
done

echo "========== Filtered jet-pT binning scan =========="
echo "D0 pT range:           ${MIN_PT_D0}-${MAX_PT_D0} GeV/c"
echo "Reco reference styles: ${#RECO_REFERENCE_LABELS[@]}"
echo "Reco edge shifts:      ${RECO_EDGE_SHIFTS[*]} GeV/c"
echo "True edges:            ${TRUE_EDGES}"
echo "Planned runs:          ${plan_index}"
echo "Plan:                  ${PLAN}"
echo "Scan directory:        ${SCAN_DIR}"
echo

# =========================
# Summary header
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\treco_pt_min\treco_pt_max\treco_n_bins\ttrue_pt_min\ttrue_pt_max\ttrue_n_bins\tmin_width\tmax_width\tstep\ttrend\treco_edges\ttrue_edges\toverride_file\troot_log\treco_reference\tedge_shift\treco_edges_base\n" > "${SUMMARY}"
fi

# =========================
# Main loop
# =========================
runs_started=0
ACLIC_SUFFIX="++"

while IFS=$'\t' read -r PLAN_INDEX RECO_REFERENCE EDGE_SHIFT RECO_EDGES_BASE RECO_EDGES TRUE_EDGES_THIS; do
  [[ "${PLAN_INDEX}" == "plan_index" ]] && continue

  printf -v RUN_ID "r%06d" "${PLAN_INDEX}"

  if run_is_complete "${RUN_ID}"; then
    echo "[resume] ${RUN_ID} ${RECO_REFERENCE}, shift=${EDGE_SHIFT} is already complete -> skip"
    continue
  fi

  if (( MAX_TOTAL_RUNS > 0 && runs_started >= MAX_TOTAL_RUNS )); then
    echo "Reached MAX_TOTAL_RUNS=${MAX_TOTAL_RUNS}, stopping."
    exit 0
  fi

  remove_run_from_file "${RUN_ID}" "${SUMMARY}"
  remove_run_from_file "${RUN_ID}" "${STABILITY}"

  read -r RECO_PT_MIN RECO_PT_MAX RECO_N_BINS RECO_MIN_WIDTH RECO_MAX_WIDTH \
    < <(edge_stats "${RECO_EDGES}")

  read -r TRUE_PT_MIN TRUE_PT_MAX TRUE_N_BINS TRUE_MIN_WIDTH TRUE_MAX_WIDTH \
    < <(edge_stats "${TRUE_EDGES_THIS}")

  # Keep the same safety condition as in the old jet-pT scan.
  if (( TRUE_N_BINS > RECO_N_BINS )); then
    echo "[skip] ${RUN_ID}: N_true=${TRUE_N_BINS} > N_reco=${RECO_N_BINS}"
    continue
  fi

  ((runs_started += 1))

  OVR_FILE="${OVR_DIR}/override_${RUN_ID}_${RECO_REFERENCE}_shift${EDGE_SHIFT}.C"
  OUT_DIR="${RUN_DIR}/${RUN_ID}"
  ROOT_LOG="${OUT_DIR}/root.log"

  rm -rf "${OUT_DIR}"
  mkdir -p "${OUT_DIR}"

  # -------------------------
  # Override macro for Machine.C
  # -------------------------
  cat > "${OVR_FILE}" <<EOF_OVR
{
  cout << "[override] filtered jet-pT binning scan" << endl;
  cout << "[override] source style: ${RECO_REFERENCE}" << endl;
  cout << "[override] coherent reco shift: ${EDGE_SHIFT} GeV/c" << endl;
  cout << "[override] jet-pT reco edges: ${RECO_EDGES}" << endl;
  cout << "[override] jet-pT true edges: ${TRUE_EDGES_THIS}" << endl;

  // Binning variations are built directly from the TTree response.
  UseCachedRM = true;
  FillStandardRM = false;
  FillCacheRM = false;

  useCustomPtRecoBins = true;

  // Match the old jet-pT scan: one reco candidate is applied to all
  // centralities in a given run.  Centrality-specific final styles can
  // therefore still be selected independently from the scan afterwards.
  for (int ic = 0; ic < nCentralityBins; ++ic) ptRecoBinsVec[ic].clear();

  double reco_edges[] = { ${RECO_EDGES} };
  int nReco = (int)(sizeof(reco_edges) / sizeof(double));

  for (int ic = 0; ic < nCentralityBins; ++ic) {
    for (int i = 0; i < nReco; ++i) {
      ptRecoBinsVec[ic].push_back(reco_edges[i]);
    }
  }

  for (int ic = 0; ic < nCentralityBins; ++ic) ptMcBinsVecCustom[ic].clear();

  double true_edges[] = { ${TRUE_EDGES_THIS} };
  int nTrue = (int)(sizeof(true_edges) / sizeof(double));

  for (int ic = 0; ic < nCentralityBins; ++ic) {
    for (int i = 0; i < nTrue; ++i) {
      ptMcBinsVecCustom[ic].push_back(true_edges[i]);
    }
  }
}
EOF_OVR

  TREND="filtered_${RECO_REFERENCE}_shift${EDGE_SHIFT}"

  printf "%s\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${RUN_ID}" \
    "${RECO_PT_MIN}" "${RECO_PT_MAX}" "${RECO_N_BINS}" \
    "${TRUE_PT_MIN}" "${TRUE_PT_MAX}" "${TRUE_N_BINS}" \
    "${RECO_MIN_WIDTH}" "${RECO_MAX_WIDTH}" "${EDGE_STEP}" "${TREND}" \
    "${RECO_EDGES}" "${TRUE_EDGES_THIS}" "${OVR_FILE}" "${ROOT_LOG}" \
    "${RECO_REFERENCE}" "${EDGE_SHIFT}" "${RECO_EDGES_BASE}" \
    >> "${SUMMARY}"

  ROOT_CMD="${MACHINE_MACRO}${ACLIC_SUFFIX}( ${FONLL_JET}, ${CUT_NEG}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${RUN_ID}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\", \"${SCAN_DIR}\", ${USE_PRIOR_SHAPE_WEIGHTING}, ${SYSTEMATIC_SPLOT} )"
  printf '%s\n' "root -l -b -q '${ROOT_CMD}'" > "${OUT_DIR}/root_command.txt"

  echo "[${RUN_ID}] ${RECO_REFERENCE}, shift=${EDGE_SHIFT} GeV/c"
  echo "           reco: ${RECO_EDGES}"
  echo "           true: ${TRUE_EDGES_THIS}"

  if ! root -l -b -q "${ROOT_CMD}" > "${ROOT_LOG}" 2>&1; then
    echo "[error] ROOT failed for ${RUN_ID}. Last 80 lines of the log:"
    tail -n 80 "${ROOT_LOG}" || true
    exit 1
  fi

  ACLIC_SUFFIX="+"

  echo "[done] ${RUN_ID} -> ${OUT_DIR}"
done < "${PLAN}"

echo
echo "Done."
echo "Planned combinations: ${PLAN}"
echo "Summary:              ${SUMMARY}"
echo "Stability:            ${STABILITY}"
echo "Machine runs:         ${RUN_DIR}"
echo "Spectra:              ${SCAN_DIR}/Output/"
