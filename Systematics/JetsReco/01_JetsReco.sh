#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Jet-reconstruction systematic scan
#
# Runs both phases:
#   1) Simple_splot.C for _systematicSPlot = 21..27
#   2) Machine.C      for _systematicSPlot = 21..27
#
# Run from the main project directory, i.e. the directory containing:
#   config.h, config_hist.h, sPlot/Simple_splot.C, Unfolding/Machine.C, Data/...
#
# Recommended run command:
#   bash Systematics/JetsReco/01_JetsReco.sh
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

# -------------------------
# Jet-reconstruction systematics
# -------------------------
# These codes are handled inside Simple_splot.C and Machine.C.
# For SYS = 21..27, Machine.C switches from cached RM to TTree RM internally
# and selects the corresponding alternative MC jet file.
SYS_CODES=(0 21 22 23 24 25 26 27)

# Optional early checks for the MC files used by Machine.C for SYS=21..27.
# Set to 0 if you want ROOT/Machine.C to handle missing files instead.
CHECK_MC_INPUTS=1

# -------------------------
# Simple_splot() and Machine() parameters
# -------------------------
MIN_PT_D0=1
MAX_PT_D0=10

FONLL_JET=1
CUT_NEG=1
MIN_JET_PT_RECO_CUT=-30
SAVED_ITER=4

# Jet-reconstruction scan does not change prior-shape weighting.
USE_PRIOR_SHAPE=0

# Safety limit for the Machine phase. 0 = no limit.
MAX_TOTAL_RUNS=0

# -------------------------
# Output location
# -------------------------
SCAN_DIR="${SCRIPT_DIR}/scanJetsReco"
OVR_DIR="${SCAN_DIR}/overrides"
RUN_DIR="${SCAN_DIR}/runs"
SPLOT_LOG_DIR="${SCAN_DIR}/splot_logs"
SUMMARY="${SCAN_DIR}/summary.tsv"
STABILITY="${SCAN_DIR}/stability.tsv"

mkdir -p "${OVR_DIR}" "${RUN_DIR}" "${SPLOT_LOG_DIR}" "${SCAN_DIR}/Output" "${PROJECT_DIR}/Output" "${PROJECT_DIR}/OutputPdf"

# =========================
# Helper functions
# =========================
jets_reco_label() {
  case "$1" in
    0) printf 'nominal' ;;
    21) printf 'jet_rec_efficiency' ;;
    22) printf 'jet_nHitsFit13' ;;
    23) printf 'jet_nHitsFit17' ;;
    24) printf 'jet_kTDrop' ;;
    25) printf 'jet_DCA2_8' ;;
    26) printf 'jet_DCA3_2' ;;
    27) printf 'jet_hadronicCorr' ;;
    *)  printf 'unknown' ;;
  esac
}

jets_reco_suffix() {
  case "$1" in
    0)  printf '' ;;
    21) printf '_jetRecEfficiency' ;;
    22) printf '_jetnHitsFit13' ;;
    23) printf '_jetnHitsFit17' ;;
    24) printf '_jetKTDrop' ;;
    25) printf '_jetDCA2_8' ;;
    26) printf '_jetDCA3_2' ;;
    27) printf '_jetHadronicCorr' ;;
    *)  printf '' ;;
  esac
}

jets_reco_mc_file() {
  case "$1" in
    21) printf '%s/Data/Output_MC_MidLow_trackEff_05052026.root' "${PROJECT_DIR}" ;;
    22) printf '%s/Data/Output_MC_MidLow_nHitsFit13_06052026.root' "${PROJECT_DIR}" ;;
    23) printf '%s/Data/Output_MC_MidLow_nHitsFit17_06052026.root' "${PROJECT_DIR}" ;;
    24) printf '%s/Data/Output_MC_MidLow_kTDrop_07052026.root' "${PROJECT_DIR}" ;;
    25) printf '%s/Data/Output_MC_MidLow_DCA2_8_09052026.root' "${PROJECT_DIR}" ;;
    26) printf '%s/Data/Output_MC_MidLow_DCA3_2_10052026.root' "${PROJECT_DIR}" ;;
    27) printf '%s/Data/Output_MC_MidLow_hadrCorr_08052026.root' "${PROJECT_DIR}" ;;
    *)  printf '' ;;
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

jets_reco_outputs_exist() {
  local sys="$1"
  local suffix
  suffix="$(jets_reco_suffix "${sys}")"

  [[ -f "${PROJECT_DIR}/Output/Output_ICS${suffix}2.root" ]] && \
  [[ -f "${PROJECT_DIR}/Output/Output_AREA${suffix}2.root" ]] && \
  [[ -f "${PROJECT_DIR}/Output/Output_D0${suffix}2.root" ]]
}

# =========================
# Optional input checks
# =========================
if (( CHECK_MC_INPUTS )); then
  for SYS in "${SYS_CODES[@]}"; do
    MC_FILE="$(jets_reco_mc_file "${SYS}")"
    [[ -f "${MC_FILE}" ]] || { echo "[error] Missing MC input for SYS=${SYS}: ${MC_FILE}"; exit 1; }
  done
fi

# =========================
# Summary header
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\tsystematicSPlot\tjets_reco_label\tmc_file\toverride_file\tsplot_log\tmachine_log\n" > "${SUMMARY}"
fi

# =========================
# PHASE 1: prepare sWeights for all jet-reconstruction variations
# =========================
echo "========== PHASE 1: Simple_splot.C for jet-reconstruction systematics =========="

SPLOT_ACLIC_SUFFIX="++"

for SYS in "${SYS_CODES[@]}"; do
  LABEL="$(jets_reco_label "${SYS}")"
  SUFFIX="$(jets_reco_suffix "${SYS}")"
  SPLOG="${SPLOT_LOG_DIR}/splot_sys${SYS}_${LABEL}.log"
  DONE="${SPLOT_LOG_DIR}/splot_sys${SYS}_${LABEL}.done"

  if [[ -f "${DONE}" ]] && jets_reco_outputs_exist "${SYS}"; then
    echo "[sPlot resume] sys=${SYS} (${LABEL}) is already complete -> skip"
    SPLOT_ACLIC_SUFFIX="+"
    continue
  fi

  rm -f "${DONE}"

  SPLOT_CMD="${SPLOT_MACRO}${SPLOT_ACLIC_SUFFIX}( \"${INPUT_FILE}\", \"Output\", \"Output2\", ${MIN_PT_D0}, ${MAX_PT_D0}, ${SYS} )"
  printf '%s\n' "root -l -b -q '${SPLOT_CMD}'" > "${SPLOT_LOG_DIR}/splot_sys${SYS}_${LABEL}_command.txt"

  echo "[sPlot] sys=${SYS} (${LABEL}), suffix=${SUFFIX}"

  if ! root -l -b -q "${SPLOT_CMD}" > "${SPLOG}" 2>&1; then
    echo "[error] Simple_splot.C failed for sys=${SYS} (${LABEL}). Last 80 lines of the log:"
    tail -n 80 "${SPLOG}" || true
    exit 1
  fi

  SPLOT_ACLIC_SUFFIX="+"

  if ! jets_reco_outputs_exist "${SYS}"; then
    echo "[error] Simple_splot.C finished, but expected outputs are missing for sys=${SYS}. See ${SPLOG}"
    echo "        Expected for example: ${PROJECT_DIR}/Output/Output_ICS${SUFFIX}2.root"
    exit 1
  fi

  touch "${DONE}"
  echo "[sPlot done] sys=${SYS} (${LABEL}) -> ${SPLOG}"
done

# =========================
# PHASE 2: unfolding for all jet-reconstruction variations
# =========================
echo "========== PHASE 2: Machine.C for jet-reconstruction systematics =========="

run_counter=0
runs_started=0
MACHINE_ACLIC_SUFFIX="++"

for SYS in "${SYS_CODES[@]}"; do
  ((run_counter += 1))
  printf -v RUN_ID "r%06d" "${run_counter}"

  LABEL="$(jets_reco_label "${SYS}")"
  MC_FILE="$(jets_reco_mc_file "${SYS}")"
  SPLOG="${SPLOT_LOG_DIR}/splot_sys${SYS}_${LABEL}.log"

  if run_is_complete "${RUN_ID}"; then
    echo "[Machine resume] ${RUN_ID} sys=${SYS} (${LABEL}) is already complete -> skip"
    MACHINE_ACLIC_SUFFIX="+"
    continue
  fi

  if (( MAX_TOTAL_RUNS > 0 && runs_started >= MAX_TOTAL_RUNS )); then
    echo "Reached MAX_TOTAL_RUNS=${MAX_TOTAL_RUNS}, stopping."
    exit 0
  fi

  if ! jets_reco_outputs_exist "${SYS}"; then
    echo "[error] Missing sPlot outputs for sys=${SYS} (${LABEL}). Run/fix Simple_splot.C first."
    exit 1
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
  // Jet-reconstruction systematic scan.
  // Do not force UseCachedRM here.
  // For systematicSPlot = ${SYS}, Machine.C selects the corresponding MC file
  // and switches to TTree response-matrix loading internally.
}
EOF_OVR

  ROOT_CMD="${MACHINE_MACRO}${MACHINE_ACLIC_SUFFIX}( ${FONLL_JET}, ${CUT_NEG}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${RUN_ID}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\", \"${SCAN_DIR}\", ${USE_PRIOR_SHAPE}, ${SYS} )"
  printf '%s\n' "root -l -b -q '${ROOT_CMD}'" > "${OUT_DIR}/root_command.txt"

  echo "[Machine ${RUN_ID}] systematicSPlot=${SYS} (${LABEL})"
  echo "[Machine ${RUN_ID}] MC input: ${MC_FILE}"

  if ! root -l -b -q "${ROOT_CMD}" > "${ROOT_LOG}" 2>&1; then
    echo "[error] Machine.C failed for ${RUN_ID}, sys=${SYS} (${LABEL}). Last 80 lines of the log:"
    tail -n 80 "${ROOT_LOG}" || true
    exit 1
  fi

  MACHINE_ACLIC_SUFFIX="+"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${RUN_ID}" \
    "${SYS}" \
    "${LABEL}" \
    "${MC_FILE}" \
    "${OVR_FILE}" \
    "${SPLOG}" \
    "${ROOT_LOG}" \
    >> "${SUMMARY}"

  echo "[Machine done] ${RUN_ID} -> ${OUT_DIR}"
done

echo "Done."
echo "Summary:       ${SUMMARY}"
echo "Stability:     ${STABILITY}"
echo "sPlot logs:    ${SPLOT_LOG_DIR}"
echo "Machine runs:  ${RUN_DIR}"
echo "Spectra:       ${SCAN_DIR}/Output"
