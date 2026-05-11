#!/usr/bin/env bash
set -euo pipefail

# =========================
# Uživatelské nastavení
# =========================
MACHINE_MACRO="./SuperIterace/Machine.C"
INPUT_FILE="./Data/Output_real_final_01022026.root"

# -------------------------
# Prior-shape scan
# -------------------------
# Podle aktuální logiky:
#   0  = no weighting
#   10 = jet pT +20%
#   11 = jet pT -20%
#   20 = second observable +20%  (z/lambda/pTD podle aktuální 2D RM)
#   21 = second observable -20%
PRIOR_CODES=(0 10 11 20 21)

# -------------------------
# Parametry Machine()
# -------------------------
FONLL_JET=1
CUT_NEG=1
MIN_JET_PT_RECO_CUT=-30
SAVED_ITER=4
MIN_PT_D0=1
MAX_PT_D0=10

# bezpečnostní pojistka; 0 = bez limitu
MAX_TOTAL_RUNS=0

# -------------------------
# Kam ukládat
# -------------------------
OUT_BASE="scan_priorShape"
OVR_DIR="${OUT_BASE}/overrides"
RUN_DIR="${OUT_BASE}/runs"
SUMMARY="${OUT_BASE}/summary.tsv"
STABILITY="${OUT_BASE}/stability.tsv"

mkdir -p "${OVR_DIR}" "${RUN_DIR}" "${OUT_BASE}/Output"

# =========================
# Pomocné funkce
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
# Kontroly prostředí
# =========================
command -v root >/dev/null 2>&1 || { echo "root není v PATH"; exit 1; }
[[ -f "${MACHINE_MACRO}" ]] || { echo "Nenalezeno MACHINE_MACRO=${MACHINE_MACRO}"; exit 1; }
[[ -f "${INPUT_FILE}" ]] || { echo "Nenalezeno INPUT_FILE=${INPUT_FILE}"; exit 1; }

# =========================
# Hlavička summary
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\tusePriorShapeWeighting\tprior_label\toutput_dir\toverride_file\n" > "${SUMMARY}"
fi

# =========================
# Hlavní smyčka
# =========================
run_counter=0
runs_started=0

for USE_PRIOR in "${PRIOR_CODES[@]}"; do
  ((run_counter += 1))
  printf -v RUN_ID "r%06d" "${run_counter}"

  LABEL="$(prior_label "${USE_PRIOR}")"

  if run_is_complete "${RUN_ID}"; then
    echo "[resume] ${RUN_ID} ${LABEL} už je hotový -> skip"
    continue
  fi

  if (( MAX_TOTAL_RUNS > 0 && runs_started >= MAX_TOTAL_RUNS )); then
    echo "Dosažen MAX_TOTAL_RUNS=${MAX_TOTAL_RUNS}, končím."
    exit 0
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
  // Prior-shape systematika má běžet z TTree, ne z cache.
  UseCachedRM  = false;
  FillStandardRM = true;
  FillCacheRM = false;

  // usePriorShapeWeighting je předán jako poslední argument Machine(...).
  // Tento override nechává binning beze změny.
}
EOF2

  ROOT_CMD="${MACHINE_MACRO}+( ${FONLL_JET}, ${CUT_NEG}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${OUT_DIR}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\", \"${OUT_BASE}\", ${USE_PRIOR} )"
  printf '%s\n' "root -l -b -q '${ROOT_CMD}'" > "${OUT_DIR}/root_command.txt"

  echo "[${RUN_ID}] usePriorShapeWeighting=${USE_PRIOR} (${LABEL})"

  root -l -b -q "${ROOT_CMD}" \
    2>&1 | awk '
      /pSWeight:/ { next }
      seen { print }
      /Loading RM from TTree\.\.\./ { seen=1; print }
      /Loading RM from cache histograms\.\.\./ { seen=1; print }
    ' > "${OUT_DIR}/root.log"

  printf "%s\t%s\t%s\t%s\t%s\n" \
    "${RUN_ID}" \
    "${USE_PRIOR}" \
    "${LABEL}" \
    "${OUT_DIR}" \
    "${OVR_FILE}" \
    >> "${SUMMARY}"

  echo "[done] ${RUN_ID} -> ${OUT_DIR}"
done

echo "Hotovo."
echo "Summary:   ${SUMMARY}"
echo "Stability: ${STABILITY}"
echo "Outputs:   ${OUT_BASE}/Output"
