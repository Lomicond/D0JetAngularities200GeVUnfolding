#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Nominal run
#
# Runs both phases with the nominal setup:
#   1) Simple_splot.C for _systematicSPlot = 0
#   2) Machine.C      for _systematicSPlot = 0
#
# The Machine phase is forced to build the response matrix from TTree,
# not from CacheRM.root, by writing UseCachedRM = false into the
# override file.
#
# Run from the main project directory, i.e. the directory containing:
#   config.h, config_hist.h, sPlot/Simple_splot.C, Unfolding/Machine.C, Data/...
#
# Recommended run command:
#   bash Systematics/Nominal/01_nominal.sh
# ============================================================

# -------------------------
# Working-directory checks
# -------------------------
PROJECT_DIR="$(pwd -P)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

[[ -f "${PROJECT_DIR}/config.h" ]]              || { echo "[error] Run this script from the main project directory: missing ./config.h"; exit 1; }
[[ -f "${PROJECT_DIR}/config_hist.h" ]]         || { echo "[error] Run this script from the main project directory: missing ./config_hist.h"; exit 1; }
[[ -f "${PROJECT_DIR}/sPlot/Simple_splot.C" ]]  || { echo "[error] Missing ./sPlot/Simple_splot.C"; exit 1; }
[[ -f "${PROJECT_DIR}/Unfolding/Machine.C" ]]   || { echo "[error] Missing ./Unfolding/Machine.C"; exit 1; }

command -v root >/dev/null 2>&1 || { echo "[error] root is not in PATH"; exit 1; }

# =========================
# User settings
# =========================
SPLOT_MACRO="${PROJECT_DIR}/sPlot/Simple_splot.C"
MACHINE_MACRO="${PROJECT_DIR}/Unfolding/Machine.C"
INPUT_FILE="${PROJECT_DIR}/Data/Output_real_final_01022026.root"

[[ -f "${INPUT_FILE}" ]] || { echo "[error] Missing input file: ${INPUT_FILE}"; exit 1; }

# Nominal setup
SYSTEMATIC_SPLOT=0
LABEL="nominal"
RUN_ID="r000001"

# -------------------------
# Simple_splot() and Machine() parameters
# -------------------------
MIN_PT_D0=1
MAX_PT_D0=10

FONLL_JET=1
CUT_NEG=1
MIN_JET_PT_RECO_CUT=-30
SAVED_ITER=4

# Nominal run does not use prior-shape weighting.
USE_PRIOR_SHAPE=0

# -------------------------
# Output location
# -------------------------
SCAN_DIR="${SCRIPT_DIR}/scanNominal"
OVR_DIR="${SCAN_DIR}/overrides"
RUN_DIR="${SCAN_DIR}/runs"
SPLOT_LOG_DIR="${SCAN_DIR}/splot_logs"
SUMMARY="${SCAN_DIR}/summary.tsv"
STABILITY="${SCAN_DIR}/stability.tsv"

mkdir -p \
  "${OVR_DIR}" \
  "${RUN_DIR}" \
  "${SPLOT_LOG_DIR}" \
  "${SCAN_DIR}/Output" \
  "${PROJECT_DIR}/Output" \
  "${PROJECT_DIR}/OutputPdf"

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

nominal_splot_outputs_exist() {
  [[ -f "${PROJECT_DIR}/Output/Output_ICS2.root" ]] && \
  [[ -f "${PROJECT_DIR}/Output/Output_AREA2.root" ]] && \
  [[ -f "${PROJECT_DIR}/Output/Output_D02.root" ]]
}

# =========================
# Summary header
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\tsystematicSPlot\tlabel\toverride_file\tsplot_log\tmachine_log\n" > "${SUMMARY}"
fi

# =========================
# PHASE 1: nominal sWeights
# =========================
echo "========== PHASE 1: Simple_splot.C nominal =========="

SPLOG="${SPLOT_LOG_DIR}/splot_sys${SYSTEMATIC_SPLOT}_${LABEL}.log"
DONE="${SPLOT_LOG_DIR}/splot_sys${SYSTEMATIC_SPLOT}_${LABEL}.done"

if [[ -f "${DONE}" ]] && nominal_splot_outputs_exist; then
  echo "[sPlot resume] nominal sPlot is already complete -> skip"
else
  rm -f "${DONE}"

  SPLOT_CMD="${SPLOT_MACRO}++( \"${INPUT_FILE}\", \"Output\", \"Output2\", ${MIN_PT_D0}, ${MAX_PT_D0}, ${SYSTEMATIC_SPLOT} )"
  printf '%s\n' "root -l -b -q '${SPLOT_CMD}'" > "${SPLOT_LOG_DIR}/splot_sys${SYSTEMATIC_SPLOT}_${LABEL}_command.txt"

  echo "[sPlot] systematicSPlot=${SYSTEMATIC_SPLOT} (${LABEL})"

  if ! root -l -b -q "${SPLOT_CMD}" > "${SPLOG}" 2>&1; then
    echo "[error] Simple_splot.C failed for nominal. Last 80 lines of the log:"
    tail -n 80 "${SPLOG}" || true
    exit 1
  fi

  if ! nominal_splot_outputs_exist; then
    echo "[error] Simple_splot.C finished, but expected nominal outputs are missing. See ${SPLOG}"
    echo "        Expected:"
    echo "          ${PROJECT_DIR}/Output/Output_ICS2.root"
    echo "          ${PROJECT_DIR}/Output/Output_AREA2.root"
    echo "          ${PROJECT_DIR}/Output/Output_D02.root"
    exit 1
  fi

  touch "${DONE}"
  echo "[sPlot done] nominal -> ${SPLOG}"
fi

# =========================
# PHASE 2: nominal unfolding
# =========================
echo "========== PHASE 2: Machine.C nominal =========="

if run_is_complete "${RUN_ID}"; then
  echo "[Machine resume] ${RUN_ID} nominal is already complete -> skip"
else
  if ! nominal_splot_outputs_exist; then
    echo "[error] Missing nominal sPlot outputs. Run/fix Simple_splot.C first."
    exit 1
  fi

  remove_run_from_file "${RUN_ID}" "${SUMMARY}"
  remove_run_from_file "${RUN_ID}" "${STABILITY}"

  OVR_FILE="${OVR_DIR}/override_${RUN_ID}_${LABEL}.C"
  OUT_DIR="${RUN_DIR}/${RUN_ID}"
  ROOT_LOG="${OUT_DIR}/root.log"

  rm -rf "${OUT_DIR}"
  mkdir -p "${OUT_DIR}"

  cat > "${OVR_FILE}" <<EOF_OVR
{
  // Nominal run, but force full TTree response-matrix loading.
  // This avoids reading the prepared CacheRM.root response histograms.
  UseCachedRM = false;
  FillStandardRM = true;
  FillCacheRM = false;
}
EOF_OVR

  ROOT_CMD="${MACHINE_MACRO}++( ${FONLL_JET}, ${CUT_NEG}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${RUN_ID}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\", \"${SCAN_DIR}\", ${USE_PRIOR_SHAPE}, ${SYSTEMATIC_SPLOT} )"
  printf '%s\n' "root -l -b -q '${ROOT_CMD}'" > "${OUT_DIR}/root_command.txt"

  echo "[Machine ${RUN_ID}] systematicSPlot=${SYSTEMATIC_SPLOT} (${LABEL})"
  echo "[Machine ${RUN_ID}] override: ${OVR_FILE}"
  echo "[Machine ${RUN_ID}] UseCachedRM forced to false"

  if ! root -l -b -q "${ROOT_CMD}" > "${ROOT_LOG}" 2>&1; then
    echo "[error] Machine.C failed for ${RUN_ID}, nominal. Last 80 lines of the log:"
    tail -n 80 "${ROOT_LOG}" || true
    exit 1
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${RUN_ID}" \
    "${SYSTEMATIC_SPLOT}" \
    "${LABEL}" \
    "${OVR_FILE}" \
    "${SPLOG}" \
    "${ROOT_LOG}" \
    >> "${SUMMARY}"

  echo "[Machine done] ${RUN_ID} -> ${OUT_DIR}"
fi

echo "Done."
echo "Summary:       ${SUMMARY}"
echo "Stability:     ${STABILITY}"
echo "sPlot logs:    ${SPLOT_LOG_DIR}"
echo "Machine runs:  ${RUN_DIR}"
echo "Spectra:       ${SCAN_DIR}/Output"
