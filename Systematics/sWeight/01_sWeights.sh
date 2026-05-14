#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# sWeight / sPlot systematic scan
#
# Run from the main project directory, i.e. the directory containing:
#   config.h, config_hist.h, sPlot/Simple_splot.C, Unfolding/Machine.C, Data/...
#
# Recommended run command:
#   bash Systematics/sWeight/AllIn_sWeight.sh
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
# sWeight / sPlot scan
# -------------------------
# 0 = nominal: Gaussian + Exponential
# 1 = Gaussian + Chebychev2 background
# 2 = Double Gaussian signal + Exponential
# 3 = Student-t signal + Exponential
# 4 = narrower fit range
# 5 = wider fit range
# 6 = keep negative bins in unfolded spectra (_CutOfNegative=false)
SYS_CODES=(0 1 2 3 4 5 6)

# -------------------------
# Simple_splot() and Machine() parameters
# -------------------------
MIN_PT_D0=1
MAX_PT_D0=10

FONLL_JET=1
CUT_NEG=1
CUT_NEG_KEEP=0
MIN_JET_PT_RECO_CUT=-30
SAVED_ITER=4

# sWeight scan does not use prior-shape weighting.
USE_PRIOR_SHAPE_WEIGHTING=0

# Safety limit for the Machine phase. 0 = no limit.
MAX_TOTAL_RUNS=0

# -------------------------
# Output location
# -------------------------
SCAN_DIR="${SCRIPT_DIR}/scanSWeight"
OVR_DIR="${SCAN_DIR}/overrides"
RUN_DIR="${SCAN_DIR}/runs"
SPLOT_LOG_DIR="${SCAN_DIR}/splot_logs"
SUMMARY="${SCAN_DIR}/summary.tsv"
STABILITY="${SCAN_DIR}/stability.tsv"

mkdir -p "${OVR_DIR}" "${RUN_DIR}" "${SPLOT_LOG_DIR}" "${SCAN_DIR}/Output" "${PROJECT_DIR}/Output" "${PROJECT_DIR}/OutputPdf"

# =========================
# Helper functions
# =========================
splot_label() {
  case "$1" in
    0) printf 'nominal' ;;
    1) printf 'cheby2_bkg' ;;
    2) printf 'double_gauss_sig' ;;
    3) printf 'student_t_sig' ;;
    4) printf 'narrow_fit' ;;
    5) printf 'wide_fit' ;;
    6) printf 'keep_negative' ;;
    *) printf 'unknown' ;;
  esac
}

splot_suffix() {
  case "$1" in
    0) printf '' ;;
    1) printf '_ChebysevBkg' ;;
    2) printf '_DoubleGaussSgn' ;;
    3) printf '_StudentTSignal' ;;
    4) printf '_NarrowFit' ;;
    5) printf '_WideFit' ;;
    6) printf '_KeepNegative' ;;
    *) printf '' ;;
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

splot_outputs_exist() {
  local sys="$1"
  local suffix
  suffix="$(splot_suffix "${sys}")"

  [[ -f "${PROJECT_DIR}/Output/Output_ICS${suffix}2.root" ]] && \
  [[ -f "${PROJECT_DIR}/Output/Output_AREA${suffix}2.root" ]] && \
  [[ -f "${PROJECT_DIR}/Output/Output_D0${suffix}2.root" ]]
}

write_sweight_override() {
  local file="$1"

  cat > "${file}" <<'EOF_OVR'
{
  // sWeight systematic scan:
  // keep the response matrix fixed and vary only the sWeights / negative-bin treatment.
  UseCachedRM    = true;
  FillStandardRM = false;
  FillCacheRM    = false;
}
EOF_OVR
}

# =========================
# Summary header
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\tsystematicSPlot\tsplot_label\tcut_negative\toverride_file\tsplot_log\tmachine_log\n" > "${SUMMARY}"
fi

# =========================
# Phase 1: prepare sWeights for all sPlot variations
# =========================
echo "========== PHASE 1: Simple_splot.C for all sWeight variations =========="

SPLOT_ACLIC_SUFFIX="++"

for SYS in "${SYS_CODES[@]}"; do
  LABEL="$(splot_label "${SYS}")"
  SPLOG="${SPLOT_LOG_DIR}/splot_sys${SYS}_${LABEL}.log"
  DONE="${SPLOT_LOG_DIR}/splot_sys${SYS}_${LABEL}.done"

  if [[ -f "${DONE}" ]] && splot_outputs_exist "${SYS}"; then
    echo "[sPlot resume] sys=${SYS} (${LABEL}) is already complete -> skip"
    continue
  fi

  rm -f "${DONE}"

  SPLOT_CMD="${SPLOT_MACRO}${SPLOT_ACLIC_SUFFIX}( \"${INPUT_FILE}\", \"Output\", \"Output2\", ${MIN_PT_D0}, ${MAX_PT_D0}, ${SYS} )"
  printf '%s\n' "root -l -b -q '${SPLOT_CMD}'" > "${SPLOT_LOG_DIR}/splot_sys${SYS}_${LABEL}_command.txt"

  echo "[sPlot] sys=${SYS} (${LABEL})"

  if ! root -l -b -q "${SPLOT_CMD}" > "${SPLOG}" 2>&1; then
    echo "[error] Simple_splot.C failed for sys=${SYS} (${LABEL}). Last 80 lines of the log:"
    tail -n 80 "${SPLOG}" || true
    exit 1
  fi

  # Only the first real Simple_splot() call should force recompilation.
  SPLOT_ACLIC_SUFFIX="+"

  if ! splot_outputs_exist "${SYS}"; then
    echo "[error] Simple_splot.C finished, but expected output files are missing for sys=${SYS} (${LABEL}). Log: ${SPLOG}"
    echo "[error] Expected files in ${PROJECT_DIR}/Output with suffix '$(splot_suffix "${SYS}")2.root'."
    exit 1
  fi

  touch "${DONE}"
  echo "[sPlot done] sys=${SYS} (${LABEL}) -> ${SPLOG}"
done

# =========================
# Phase 2: unfolding for all prepared sWeight variations
# =========================
echo "========== PHASE 2: Machine.C for all sWeight variations =========="

run_counter=0
runs_started=0
MACHINE_ACLIC_SUFFIX="++"

for SYS in "${SYS_CODES[@]}"; do
  ((run_counter += 1))
  printf -v RUN_ID "r%06d" "${run_counter}"

  LABEL="$(splot_label "${SYS}")"
  SPLOG="${SPLOT_LOG_DIR}/splot_sys${SYS}_${LABEL}.log"

  if run_is_complete "${RUN_ID}"; then
    echo "[Machine resume] ${RUN_ID} sys=${SYS} (${LABEL}) is already complete -> skip"
    continue
  fi

  if (( MAX_TOTAL_RUNS > 0 && runs_started >= MAX_TOTAL_RUNS )); then
    echo "Reached MAX_TOTAL_RUNS=${MAX_TOTAL_RUNS}, stopping."
    exit 0
  fi

  if ! splot_outputs_exist "${SYS}"; then
    echo "[error] Missing sPlot outputs for sys=${SYS} (${LABEL}). Run/fix Phase 1 first."
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

  write_sweight_override "${OVR_FILE}"

  CUT_NEG_THIS="${CUT_NEG}"
  if [[ "${SYS}" == "6" ]]; then
    CUT_NEG_THIS="${CUT_NEG_KEEP}"
  fi

  ROOT_CMD="${MACHINE_MACRO}${MACHINE_ACLIC_SUFFIX}( ${FONLL_JET}, ${CUT_NEG_THIS}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${RUN_ID}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\", \"${SCAN_DIR}\", ${USE_PRIOR_SHAPE_WEIGHTING}, ${SYS} )"
  printf '%s\n' "root -l -b -q '${ROOT_CMD}'" > "${OUT_DIR}/root_command.txt"

  echo "[Machine ${RUN_ID}] systematicSPlot=${SYS} (${LABEL}), CutOfNegative=${CUT_NEG_THIS}"

  if ! root -l -b -q "${ROOT_CMD}" > "${ROOT_LOG}" 2>&1; then
    echo "[error] Machine.C failed for ${RUN_ID}, sys=${SYS} (${LABEL}). Last 80 lines of the log:"
    tail -n 80 "${ROOT_LOG}" || true
    exit 1
  fi

  # Only the first real Machine() call should force recompilation.
  MACHINE_ACLIC_SUFFIX="+"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${RUN_ID}" \
    "${SYS}" \
    "${LABEL}" \
    "${CUT_NEG_THIS}" \
    "${OVR_FILE}" \
    "${SPLOG}" \
    "${ROOT_LOG}" \
    >> "${SUMMARY}"

  echo "[Machine done] ${RUN_ID} -> ${ROOT_LOG}"
done

echo "Done."
echo "Summary:    ${SUMMARY}"
echo "Stability:  ${STABILITY}"
echo "sPlot logs: ${SPLOT_LOG_DIR}"
echo "Spectra:    ${SCAN_DIR}/Output"
