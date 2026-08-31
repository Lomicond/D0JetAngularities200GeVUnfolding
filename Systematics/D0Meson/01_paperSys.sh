#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 MIN_PT_D0 MAX_PT_D0"
  echo "Example: $0 1 10"
  exit 1
fi

MIN_PT_D0="$1"
MAX_PT_D0="$2"
MIN_PT_D0_SPLOT=1 #For sPlot, it can be wider
MAX_PT_D0_SPLOT=10

if ! [[ "${MIN_PT_D0}" =~ ^[0-9]+$ && "${MAX_PT_D0}" =~ ^[0-9]+$ ]]; then
  echo "[error] MIN_PT_D0 and MAX_PT_D0 must be integers."
  exit 1
fi

if (( MIN_PT_D0 >= MAX_PT_D0 )); then
  echo "[error] MIN_PT_D0 must be smaller than MAX_PT_D0."
  exit 1
fi

# ============================================================
# Paper D0-reconstruction systematic scan
#
# Run from the main project directory, i.e. the directory containing:
#   config.h, config_hist.h, sPlot/Simple_splot.C, Unfolding/Machine.C, Data/...
#
# Recommended run command:
#   bash Systematics/D0Meson/AllIn_paperSys.sh
# ============================================================

# -------------------------
# Working-directory checks
# -------------------------
PROJECT_DIR="$(pwd -P)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

[[ -f "${PROJECT_DIR}/config.h" ]]             || { echo "[error] Run this script from the main project directory: missing ./config.h"; exit 1; }
[[ -f "${PROJECT_DIR}/config_hist.h" ]]        || { echo "[error] Run this script from the main project directory: missing ./config_hist.h"; exit 1; }
[[ -f "${PROJECT_DIR}/sPlot/Simple_splot.C" ]] || { echo "[error] Missing ./sPlot/Simple_splot.C"; exit 1; }
[[ -f "${PROJECT_DIR}/Unfolding/Machine.C" ]]  || { echo "[error] Missing ./Unfolding/Machine.C"; exit 1; }

command -v root >/dev/null 2>&1 || { echo "[error] root is not in PATH"; exit 1; }

# =========================
# User settings
# =========================
SPLOT_MACRO="${PROJECT_DIR}/sPlot/Simple_splot.C"
MACHINE_MACRO="${PROJECT_DIR}/Unfolding/Machine.C"
INPUT_FILE="${PROJECT_DIR}/Data/Output_real_final_01022026.root"

[[ -f "${INPUT_FILE}" ]] || { echo "[error] Missing input file: ${INPUT_FILE}"; exit 1; }

# -------------------------
# Paper D0-reconstruction systematics
# -------------------------
# 0  = nominal
# 7  = paper TPC tracking up
# 8  = paper TPC tracking down
# 9  = paper PID up
# 10 = paper PID down
# 11 = paper single-track pT up
# 12 = paper single-track pT down
# 13 = paper topological efficiency up
# 14 = paper topological efficiency down
# 15 = paper double-counting up
# 16 = paper double-counting down
# 17 = paper vertex correction up
# 18 = paper vertex correction down
# 19 = paper secondary-track up
# 20 = paper secondary-track down
SYS_CODES=(0 7 8 9 10 11 12 13 14 15 16 17 18 19 20)

# -------------------------
# Simple_splot() and Machine() parameters
# -------------------------
FONLL_JET=1
CUT_NEG=1
MIN_JET_PT_RECO_CUT=-30
SAVED_ITER=4

# Paper D0-reco systematics do not change prior-shape weighting.
USE_PRIOR_SHAPE=0

# Safety limit for the Machine phase. 0 = no limit.
MAX_TOTAL_RUNS=0

# -------------------------
# Output location
# -------------------------
SCAN_DIR="${SCRIPT_DIR}/scanPaperSys_${MIN_PT_D0}_${MAX_PT_D0}"
OVR_DIR="${SCAN_DIR}/overrides"
RUN_DIR="${SCAN_DIR}/runs"
SPLOT_LOG_DIR="${SCAN_DIR}/paperSys_logs"
SUMMARY="${SCAN_DIR}/summary.tsv"
STABILITY="${SCAN_DIR}/stability.tsv"

mkdir -p "${OVR_DIR}" "${RUN_DIR}" "${SPLOT_LOG_DIR}" "${SCAN_DIR}/Output" "${PROJECT_DIR}/Output" "${PROJECT_DIR}/OutputPdf"

# =========================
# Helper functions
# =========================
paper_label() {
  case "$1" in
    0)  printf 'nominal' ;;
    7)  printf 'paper_tpc_track_up' ;;
    8)  printf 'paper_tpc_track_down' ;;
    9)  printf 'paper_pid_up' ;;
    10) printf 'paper_pid_down' ;;
    11) printf 'paper_single_track_pt_up' ;;
    12) printf 'paper_single_track_pt_down' ;;
    13) printf 'paper_topo_eff_up' ;;
    14) printf 'paper_topo_eff_down' ;;
    15) printf 'paper_double_counting_up' ;;
    16) printf 'paper_double_counting_down' ;;
    17) printf 'paper_vertex_corr_up' ;;
    18) printf 'paper_vertex_corr_down' ;;
    19) printf 'paper_secondary_track_up' ;;
    20) printf 'paper_secondary_track_down' ;;
    *)  printf 'unknown' ;;
  esac
}

paper_suffix() {
  case "$1" in
    0)  printf '' ;;
    7)  printf '_paperTPCTrackUp' ;;
    8)  printf '_paperTPCTrackDown' ;;
    9)  printf '_paperPIDUp' ;;
    10) printf '_paperPIDDown' ;;
    11) printf '_paperSingleTrackPtUp' ;;
    12) printf '_paperSingleTrackPtDown' ;;
    13) printf '_paperTopoEffUp' ;;
    14) printf '_paperTopoEffDown' ;;
    15) printf '_paperDoubleCountingUp' ;;
    16) printf '_paperDoubleCountingDown' ;;
    17) printf '_paperVertexCorrUp' ;;
    18) printf '_paperVertexCorrDown' ;;
    19) printf '_paperSecondaryTrackUp' ;;
    20) printf '_paperSecondaryTrackDown' ;;
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

paper_outputs_exist() {
  local sys="$1"
  local suffix
  suffix="$(paper_suffix "${sys}")"

  [[ -f "${PROJECT_DIR}/Output/Output_ICS${suffix}2.root" ]] && \
  [[ -f "${PROJECT_DIR}/Output/Output_AREA${suffix}2.root" ]] && \
  [[ -f "${PROJECT_DIR}/Output/Output_D0${suffix}2.root" ]]
}

# =========================
# Summary header
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\tsystematicSPlot\tpaper_label\toverride_file\tpaperSys_log\tmachine_log\n" > "${SUMMARY}"
fi

# =========================
# PHASE 1: create paperSys sPlot outputs
# =========================
echo "========== PHASE 1: Simple_splot.C for paper D0-reco systematics =========="

SPLOT_ACLIC_SUFFIX="++"

for SYS in "${SYS_CODES[@]}"; do
  LABEL="$(paper_label "${SYS}")"
  SPLOG="${SPLOT_LOG_DIR}/paperSys_sys${SYS}_${LABEL}.log"
  DONE="${SPLOT_LOG_DIR}/paperSys_sys${SYS}_${LABEL}.done"

  if [[ -f "${DONE}" ]] && paper_outputs_exist "${SYS}"; then
    echo "[paperSys resume] sys=${SYS} (${LABEL}) is already complete -> skip"
    continue
  fi

  rm -f "${DONE}"

  SPLOT_CMD="${SPLOT_MACRO}${SPLOT_ACLIC_SUFFIX}( \"${INPUT_FILE}\", \"Output\", \"Output2\", ${MIN_PT_D0_SPLOT}, ${MAX_PT_D0_SPLOT}, ${SYS} )"
  printf '%s\n' "root -l -b -q '${SPLOT_CMD}'" > "${SPLOT_LOG_DIR}/paperSys_sys${SYS}_${LABEL}_command.txt"

  echo "[paperSys] sys=${SYS} (${LABEL})"

  if ! root -l -b -q "${SPLOT_CMD}" > "${SPLOG}" 2>&1; then
    echo "[error] Simple_splot.C failed for sys=${SYS} (${LABEL}). Last 80 lines of the log:"
    tail -n 80 "${SPLOG}" || true
    exit 1
  fi

  # Only the first actual Simple_splot() run should force recompilation.
  SPLOT_ACLIC_SUFFIX="+"

  if ! paper_outputs_exist "${SYS}"; then
    echo "[error] Simple_splot.C finished, but expected output is missing for sys=${SYS} (${LABEL})."
    echo "[error] Expected for example: ${PROJECT_DIR}/Output/Output_ICS$(paper_suffix "${SYS}")2.root"
    echo "[error] See log: ${SPLOG}"
    exit 1
  fi

  touch "${DONE}"
  echo "[paperSys done] sys=${SYS} (${LABEL}) -> ${SPLOG}"
done

# =========================
# PHASE 2: unfold all paperSys variants
# =========================
echo "========== PHASE 2: Machine.C for paper D0-reco systematics =========="

run_counter=0
runs_started=0
MACHINE_ACLIC_SUFFIX="++"

for SYS in "${SYS_CODES[@]}"; do
  ((run_counter += 1))
  printf -v RUN_ID "r%06d" "${run_counter}"

  LABEL="$(paper_label "${SYS}")"
  SPLOG="${SPLOT_LOG_DIR}/paperSys_sys${SYS}_${LABEL}.log"

  if run_is_complete "${RUN_ID}"; then
    echo "[Machine resume] ${RUN_ID} sys=${SYS} (${LABEL}) is already complete -> skip"
    continue
  fi

  if (( MAX_TOTAL_RUNS > 0 && runs_started >= MAX_TOTAL_RUNS )); then
    echo "Reached MAX_TOTAL_RUNS=${MAX_TOTAL_RUNS}, stopping."
    exit 0
  fi

  if ! paper_outputs_exist "${SYS}"; then
    echo "[error] Missing paperSys outputs for sys=${SYS} (${LABEL}). Run/fix Simple_splot.C first."
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
  // Keep the response matrix fixed for paper D0-reconstruction systematics.
  // Only the data-side D0-candidate weights / reconstructed input spectra change.
  UseCachedRM    = true;
  FillStandardRM = false;
  FillCacheRM    = false;
}
EOF_OVR

  ROOT_CMD="${MACHINE_MACRO}${MACHINE_ACLIC_SUFFIX}( ${FONLL_JET}, ${CUT_NEG}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${RUN_ID}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\", \"${SCAN_DIR}\", ${USE_PRIOR_SHAPE}, ${SYS} )"
  printf '%s\n' "root -l -b -q '${ROOT_CMD}'" > "${OUT_DIR}/root_command.txt"

  echo "[Machine ${RUN_ID}] systematicSPlot=${SYS} (${LABEL})"

  if ! root -l -b -q "${ROOT_CMD}" > "${ROOT_LOG}" 2>&1; then
    echo "[error] Machine.C failed for ${RUN_ID}, sys=${SYS} (${LABEL}). Last 80 lines of the log:"
    tail -n 80 "${ROOT_LOG}" || true
    exit 1
  fi

  # Only the first actual Machine() run should force recompilation.
  MACHINE_ACLIC_SUFFIX="+"

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${RUN_ID}" \
    "${SYS}" \
    "${LABEL}" \
    "${OVR_FILE}" \
    "${SPLOG}" \
    "${ROOT_LOG}" \
    >> "${SUMMARY}"

  echo "[Machine done] ${RUN_ID} -> ${OUT_DIR}"
done

echo "Done."
echo "Summary:       ${SUMMARY}"
echo "Stability:     ${STABILITY}"
echo "paperSys logs: ${SPLOT_LOG_DIR}"
echo "Outputs:       ${SCAN_DIR}/Output"
