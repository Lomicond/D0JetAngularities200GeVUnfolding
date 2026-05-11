#!/usr/bin/env bash
set -euo pipefail

# =========================
# Uživatelské nastavení
# =========================
SPLOT_MACRO="./sPlot/Simple_splot.C"
MACHINE_MACRO="./SuperIterace/Machine.C"
INPUT_FILE="./Data/Output_real_final_01022026.root"

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
# Parametry Simple_splot() a Machine()
# -------------------------
MIN_PT_D0=1
MAX_PT_D0=10

FONLL_JET=1
CUT_NEG=1
CUT_NEG_KEEP=0
MIN_JET_PT_RECO_CUT=-30
SAVED_ITER=4

# sWeight scan nemění prior-shape weighting.
USE_PRIOR_SHAPE=0

# bezpečnostní pojistka pro Machine fázi; 0 = bez limitu
MAX_TOTAL_RUNS=0

# -------------------------
# Kam ukládat
# -------------------------
OUT_BASE="scan_sWeight"
OVR_DIR="${OUT_BASE}/overrides"
RUN_DIR="${OUT_BASE}/runs"
SPLOT_LOG_DIR="${OUT_BASE}/splot_logs"
SUMMARY="${OUT_BASE}/summary.tsv"
STABILITY="${OUT_BASE}/stability.tsv"

mkdir -p "${OVR_DIR}" "${RUN_DIR}" "${SPLOT_LOG_DIR}" "${OUT_BASE}/Output" "Output" "OutputPdf"

# =========================
# Pomocné funkce
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

  [[ -f "./Output/Output_ICS${suffix}2.root" ]] && \
  [[ -f "./Output/Output_AREA${suffix}2.root" ]] && \
  [[ -f "./Output/Output_D0${suffix}2.root" ]]
}

# =========================
# Kontroly prostředí
# =========================
command -v root >/dev/null 2>&1 || { echo "root není v PATH"; exit 1; }
[[ -f "${SPLOT_MACRO}" ]]   || { echo "Nenalezeno SPLOT_MACRO=${SPLOT_MACRO}"; exit 1; }
[[ -f "${MACHINE_MACRO}" ]] || { echo "Nenalezeno MACHINE_MACRO=${MACHINE_MACRO}"; exit 1; }
[[ -f "${INPUT_FILE}" ]]    || { echo "Nenalezeno INPUT_FILE=${INPUT_FILE}"; exit 1; }

# =========================
# Hlavička summary
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\tsystematicSPlot\tsplot_label\toutput_dir\toverride_file\tsplot_log\tmachine_log\n" > "${SUMMARY}"
fi

# =========================
# FÁZE 1: připrav sWeights pro všechny systematiky
# =========================
echo "========== PHASE 1: Simple_splot.C pro všechny _systematicSPlot =========="

for SYS in "${SYS_CODES[@]}"; do
  LABEL="$(splot_label "${SYS}")"
  SPLOG="${SPLOT_LOG_DIR}/splot_sys${SYS}_${LABEL}.log"
  DONE="${SPLOT_LOG_DIR}/splot_sys${SYS}_${LABEL}.done"

  if [[ -f "${DONE}" ]] && splot_outputs_exist "${SYS}"; then
    echo "[sPlot resume] sys=${SYS} (${LABEL}) už je hotový -> skip"
    continue
  fi

  rm -f "${DONE}"

  SPLOT_CMD="${SPLOT_MACRO}+( \"${INPUT_FILE}\", \"Output\", \"Output2\", ${MIN_PT_D0}, ${MAX_PT_D0}, ${SYS} )"
  printf '%s\n' "root -l -b -q '${SPLOT_CMD}'" > "${SPLOT_LOG_DIR}/splot_sys${SYS}_${LABEL}_command.txt"

  echo "[sPlot] sys=${SYS} (${LABEL})"

  root -l -b -q "${SPLOT_CMD}" \
    > "${SPLOG}" 2>&1

  if ! splot_outputs_exist "${SYS}"; then
    echo "[ERROR] sPlot doběhl, ale chybí očekávaný výstup pro sys=${SYS}. Viz ${SPLOG}"
    exit 1
  fi

  touch "${DONE}"
  echo "[sPlot done] sys=${SYS} (${LABEL}) -> ${SPLOG}"
done

# =========================
# FÁZE 2: unfolding pro všechny připravené sWeight varianty
# =========================
echo "========== PHASE 2: Machine.C pro všechny _systematicSPlot =========="

run_counter=0
runs_started=0

for SYS in "${SYS_CODES[@]}"; do
  ((run_counter += 1))
  printf -v RUN_ID "r%06d" "${run_counter}"

  LABEL="$(splot_label "${SYS}")"
  SPLOG="${SPLOT_LOG_DIR}/splot_sys${SYS}_${LABEL}.log"

  if run_is_complete "${RUN_ID}"; then
    echo "[Machine resume] ${RUN_ID} sys=${SYS} (${LABEL}) už je hotový -> skip"
    continue
  fi

  if (( MAX_TOTAL_RUNS > 0 && runs_started >= MAX_TOTAL_RUNS )); then
    echo "Dosažen MAX_TOTAL_RUNS=${MAX_TOTAL_RUNS}, končím."
    exit 0
  fi

  if ! splot_outputs_exist "${SYS}"; then
    echo "[ERROR] Chybí sPlot výstupy pro sys=${SYS} (${LABEL}). Nejdřív oprav/spusť Simple_splot."
    exit 1
  fi

  remove_run_from_file "${RUN_ID}" "${SUMMARY}"
  remove_run_from_file "${RUN_ID}" "${STABILITY}"

  ((runs_started += 1))

  OVR_FILE="${OVR_DIR}/override_${RUN_ID}_${LABEL}.C"
  OUT_DIR="${RUN_DIR}/${RUN_ID}"
  rm -rf "${OUT_DIR}"
  mkdir -p "${OUT_DIR}"

  cat > "${OVR_FILE}" <<EOF2
{
  // Pro sWeight systematiku držíme response matrix stejnou.
  // Mění se pouze sWeights / případně ořez negativních hodnot.
  UseCachedRM     = true;
  FillStandardRM  = false;
  FillCacheRM     = false;
}
EOF2

  CUT_NEG_THIS="${CUT_NEG}"
  if [[ "${SYS}" == "6" ]]; then
    CUT_NEG_THIS="${CUT_NEG_KEEP}"
  fi

  ROOT_CMD="${MACHINE_MACRO}+( ${FONLL_JET}, ${CUT_NEG_THIS}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${OUT_DIR}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\", \"${OUT_BASE}\", ${USE_PRIOR_SHAPE}, ${SYS} )"
  printf '%s\n' "root -l -b -q '${ROOT_CMD}'" > "${OUT_DIR}/root_command.txt"

  echo "[Machine ${RUN_ID}] systematicSPlot=${SYS} (${LABEL})"

  root -l -b -q "${ROOT_CMD}" \
    2>&1 | awk '
      /pSWeight:/ { next }
      seen { print }
      /Loading RM from TTree\.\.\./ { seen=1; print }
      /Loading RM from cache histograms\.\.\./ { seen=1; print }
    ' > "${OUT_DIR}/root.log"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${RUN_ID}" \
    "${SYS}" \
    "${LABEL}" \
    "${OUT_DIR}" \
    "${OVR_FILE}" \
    "${SPLOG}" \
    "${OUT_DIR}/root.log" \
    >> "${SUMMARY}"

  echo "[Machine done] ${RUN_ID} -> ${OUT_DIR}"
done

echo "Hotovo."
echo "Summary:   ${SUMMARY}"
echo "Stability: ${STABILITY}"
echo "sPlot logs: ${SPLOT_LOG_DIR}"
echo "Outputs:   ${OUT_BASE}/Output"
