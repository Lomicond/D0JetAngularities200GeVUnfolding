#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Prior-shape systematic scan
#
# Run from the main project directory, i.e. the directory containing:
#   config.h, config_hist.h, Unfolding/Machine.C, Data/...
#
# Recommended run command:
#   bash Systematics/PriorShape/AllIn_priorShape_secondVar_compatV3.sh
# ============================================================

# -------------------------
# Working-directory checks
# -------------------------
PROJECT_DIR="$(pwd -P)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

[[ -f "${PROJECT_DIR}/config.h" ]]            || { echo "[error] Run this script from the main project directory: missing ./config.h"; exit 1; }
[[ -f "${PROJECT_DIR}/config_hist.h" ]]       || { echo "[error] Run this script from the main project directory: missing ./config_hist.h"; exit 1; }
[[ -f "${PROJECT_DIR}/Unfolding/Machine.C" ]] || { echo "[error] Missing ./Unfolding/Machine.C"; exit 1; }

command -v root >/dev/null 2>&1 || { echo "[error] root is not in PATH"; exit 1; }

# =========================
# User settings
# =========================
MACHINE_MACRO="${PROJECT_DIR}/Unfolding/Machine.C"
INPUT_FILE="${PROJECT_DIR}/Data/Output_real_final_01022026.root"

[[ -f "${INPUT_FILE}" ]] || { echo "[error] Missing input file: ${INPUT_FILE}"; exit 1; }

# -------------------------
# Prior-shape scan
# -------------------------
# Machine.C convention used here:
#   0  = no prior-shape weighting
#   1X = jet pT prior-shape variation
#   2X = second-observable prior-shape variation
#
# Sign convention:
#   X = 0 -> +20%
#   X = 1 -> -20%
#
# For this scan, only the nominal case, jet-pT prior variations,
# and the current second-observable prior variations are run.
PRIOR_CODES=(0 10 11 20 21)

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
SYSTEMATIC_SPLOT=0

# Safety limit for the whole scan. 0 = no limit.
MAX_TOTAL_RUNS=0

# -------------------------
# Output location
# -------------------------
SCAN_DIR="${SCRIPT_DIR}/scanPriorShape"
OVR_DIR="${SCAN_DIR}/overrides"
RUN_DIR="${SCAN_DIR}/runs"
SUMMARY="${SCAN_DIR}/summary.tsv"
STABILITY="${SCAN_DIR}/stability.tsv"

mkdir -p "${OVR_DIR}" "${RUN_DIR}" "${PROJECT_DIR}/OutputPdf"

# =========================
# Helper functions
# =========================
prior_label() {
  case "$1" in
    0)  printf 'nominal' ;;
    10) printf 'jetPt_plus20' ;;
    11) printf 'jetPt_minus20' ;;
    20) printf 'secondVar_plus20' ;;
    21) printf 'secondVar_minus20' ;;
    *)  printf 'unknown' ;;
  esac
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
# Summary header
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\tusePriorShapeWeighting\tprior_label\toverride_file\troot_log\n" > "${SUMMARY}"
fi

# =========================
# Main loop
# =========================
run_counter=0
runs_started=0

# Use forced ACLiC rebuild only for the first real Machine() call.
# After that, use the already compiled library.
ACLIC_SUFFIX="++"

for USE_PRIOR in "${PRIOR_CODES[@]}"; do
  ((run_counter += 1))
  printf -v RUN_ID "r%06d" "${run_counter}"

  LABEL="$(prior_label "${USE_PRIOR}")"

  if run_is_complete "${RUN_ID}"; then
    echo "[resume] ${RUN_ID} ${LABEL} is already complete -> skip"
    continue
  fi

  if (( MAX_TOTAL_RUNS > 0 && runs_started >= MAX_TOTAL_RUNS )); then
    echo "Reached MAX_TOTAL_RUNS=${MAX_TOTAL_RUNS}, stopping."
    exit 0
  fi

  remove_run_from_file "${RUN_ID}" "${SUMMARY}"
  remove_run_from_file "${RUN_ID}" "${STABILITY}"

  ((runs_started += 1))

  OVR_FILE="${OVR_DIR}/override_${RUN_ID}_${LABEL}.C"
  OUT_DIR="${RUN_DIR}/${RUN_ID}"
  ROOT_LOG="${OUT_DIR}/root.log"
  rm -rf "${OUT_DIR}"
  mkdir -p "${OUT_DIR}"

  cat > "${OVR_FILE}" <<EOF_OVR
{
  // Prior-shape weighting must be built directly from the TTree.
  // Machine.C exits if usePriorShapeWeighting != 0 while UseCachedRM is true.
  UseCachedRM = false;
  FillStandardRM = true;
  FillCacheRM = false;

  // No binning override is applied here.
}
EOF_OVR

  ROOT_CMD="${MACHINE_MACRO}${ACLIC_SUFFIX}( ${FONLL_JET}, ${CUT_NEG}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${RUN_ID}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\", \"${SCAN_DIR}\", ${USE_PRIOR}, ${SYSTEMATIC_SPLOT} )"
  printf '%s\n' "root -l -b -q '${ROOT_CMD}'" > "${OUT_DIR}/root_command.txt"

  printf "%s\t%s\t%s\t%s\t%s\n" \
    "${RUN_ID}" \
    "${USE_PRIOR}" \
    "${LABEL}" \
    "${OVR_FILE}" \
    "${ROOT_LOG}" \
    >> "${SUMMARY}"

  echo "[${RUN_ID}] usePriorShapeWeighting=${USE_PRIOR} (${LABEL})"

  if ! root -l -b -q "${ROOT_CMD}" > "${ROOT_LOG}" 2>&1; then
    echo "[error] ROOT failed for ${RUN_ID} (${LABEL}). Last 80 lines of the log:"
    tail -n 80 "${ROOT_LOG}" || true
    exit 1
  fi

  # Only the first run should force recompilation.
  ACLIC_SUFFIX="+"

  echo "[done] ${RUN_ID} -> ${OUT_DIR}"
done

echo "Done."
echo "Summary:   ${SUMMARY}"
echo "Stability: ${STABILITY}"
echo "Spectra:   ${SCAN_DIR}/Output/"
