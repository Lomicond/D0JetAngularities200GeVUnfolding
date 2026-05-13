#!/usr/bin/env bash
set -euo pipefail

# =========================
# Uživatelské nastavení
# =========================
MACHINE_MACRO="./SuperIterace/Machine.C"
INPUT_FILE="./Data/Output_real_final_01022026.root"

# -------------------------
# pTD scan only
# -------------------------
# Pozn.: index 5 = p_T^D / pTD.
# Override níže mění POUZE angRecoBinsVec[*][5] a angMcBinsVecCustom[*][5].
# Všechny ostatní binningy zůstávají tak, jak jsou v config.h.
PTD_INDEX=5

PTD_TRUE_START_LIST=(0 0.3 0.5)
PTD_RECO_START_LIST=(0 0.3 0.5)
PTD_RECO_ALLOWED_WIDTHS=(0.1 0.15 0.2 0.3)
PTD_RECO_MIN_WIDTH=0.1
PTD_RIGHT_EDGE=1.01
PTD_LAST_BIN_LEFT=0.9   # poslední bin je vždy 0.9 -> 1.01, efektivně jako 0.1

# "Ševelení" hran.
# Posun se aplikuje koherentně na vnitřní hrany reco i true binningu.
# Krajní hrany, tj. první hrana daného rozsahu a PTD_RIGHT_EDGE, se neposouvají.
# Pokud posun vytvoří neplatné hrany, daná varianta se přeskočí.
PTD_EDGE_SHIFTS=(0 -0.1 0.1)

# -------------------------
# Parametry Machine()
# -------------------------
FONLL_JET=1
CUT_NEG=1
MIN_JET_PT_RECO_CUT=-30
SAVED_ITER=4
MIN_PT_D0=5
MAX_PT_D0=10

# bezpečnostní pojistky
# 0 = bez limitu
MAX_TOTAL_RUNS=0

# -------------------------
# Kam ukládat
# -------------------------
OUT_BASE="scan_pTD"
OVR_DIR="${OUT_BASE}/overrides"
RUN_DIR="${OUT_BASE}/runs"
SUMMARY="${OUT_BASE}/summary.tsv"
STABILITY="${OUT_BASE}/stability.tsv"

mkdir -p "${OVR_DIR}" "${RUN_DIR}"

# =========================
# Pomocné funkce
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

# Posouváme jen vnitřní hrany. Levý okraj aktuálního rozsahu a 1.01 necháváme fixní.
if abs(shift) > 0:
    for i in range(1, len(edges) - 1):
        edges[i] = round(edges[i] + shift, 10)

# Bezpečnostní kontrola: hrany musí zůstat striktně rostoucí.
for a, b in zip(edges, edges[1:]):
    if not (b > a):
        raise SystemExit(f"Invalid shifted edges: {edges}")

if edges[0] < 0:
    raise SystemExit(f"Invalid shifted edges: {edges}")

# pTD má pravý okraj 1.01; po posunu vnitřních hran nesmí nic přeskočit za něj.
if edges[-1] > 1.01 + 1e-9:
    raise SystemExit(f"Invalid shifted edges: {edges}")

def fmt(x):
    s = f"{x:.10f}".rstrip('0').rstrip('.')
    return s if s else '0'

print(", ".join(fmt(x) for x in edges))
PY
}

# =========================
# Kontroly prostředí
# =========================
command -v root >/dev/null 2>&1 || { echo "root není v PATH"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "python3 není v PATH"; exit 1; }

# =========================
# Hlavička summary
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\tpTD_edge_shift\tpTD_true_start\tpTD_true_n_bins\tpTD_reco_start\tpTD_reco_n_bins\tpTD_reco_edges\tpTD_true_edges\n" > "${SUMMARY}"
fi

# =========================
# Resume počítadla
# =========================
run_counter=0
runs_started=0

# =========================
# Hlavní smyčka
# =========================
for PTD_TRUE_START in "${PTD_TRUE_START_LIST[@]}"; do
  PTD_TRUE_EDGES_BASE=$(format_true_edges_from_start "${PTD_TRUE_START}")

  for PTD_RECO_START in "${PTD_RECO_START_LIST[@]}"; do
    if ! python3 - "$PTD_RECO_START" "$PTD_TRUE_START" <<'PY' >/dev/null
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

        if ! PTD_TRUE_EDGES=$(shift_edges "${PTD_TRUE_EDGES_BASE}" "${PTD_EDGE_SHIFT}"); then
          echo "[skip] pTD shift=${PTD_EDGE_SHIFT}, true start=${PTD_TRUE_START} -> invalid true edges"
          continue
        fi

        PTD_RECO_N_BINS=$(count_edges_bins "${PTD_RECO_EDGES}")
        PTD_TRUE_N_BINS=$(count_edges_bins "${PTD_TRUE_EDGES}")

        if (( PTD_RECO_N_BINS < PTD_TRUE_N_BINS )); then
          echo "[skip] pTD shift=${PTD_EDGE_SHIFT}, true start=${PTD_TRUE_START}, reco start=${PTD_RECO_START} -> Nreco=${PTD_RECO_N_BINS} < Ntrue=${PTD_TRUE_N_BINS}"
          continue
        fi

        ((run_counter += 1))
        printf -v RUN_ID "r%06d" "${run_counter}"

        if run_is_complete "${RUN_ID}"; then
          echo "[resume] ${RUN_ID} už je hotový -> skip"
          continue
        fi

        if (( MAX_TOTAL_RUNS > 0 && runs_started >= MAX_TOTAL_RUNS )); then
          echo "Dosažen MAX_TOTAL_RUNS=${MAX_TOTAL_RUNS}, končím."
          exit 0
        fi

        remove_run_from_file "${RUN_ID}" "${SUMMARY}"
        remove_run_from_file "${RUN_ID}" "${STABILITY}"

        ((runs_started += 1))

        OVR_FILE="${OVR_DIR}/override_${RUN_ID}.C"
        OUT_DIR="${RUN_DIR}/${RUN_ID}"
        rm -rf "${OUT_DIR}"
        mkdir -p "${OUT_DIR}"

        cat > "${OVR_FILE}" <<EOF2
{
  // ---------- pTD scan only, index ${PTD_INDEX} ----------
  // All other binning variables stay exactly as defined in config.h.
  angRecoBinsVec[0][${PTD_INDEX}] = std::vector<double>{ ${PTD_RECO_EDGES} };
  angRecoBinsVec[1][${PTD_INDEX}] = std::vector<double>{ ${PTD_RECO_EDGES} };
  angRecoBinsVec[2][${PTD_INDEX}] = std::vector<double>{ ${PTD_RECO_EDGES} };

  angMcBinsVecCustom[0][${PTD_INDEX}] = std::vector<double>{ ${PTD_TRUE_EDGES} };
  angMcBinsVecCustom[1][${PTD_INDEX}] = std::vector<double>{ ${PTD_TRUE_EDGES} };
  angMcBinsVecCustom[2][${PTD_INDEX}] = std::vector<double>{ ${PTD_TRUE_EDGES} };
}
EOF2

        echo "[${RUN_ID}] pTD shift=${PTD_EDGE_SHIFT}, true start=${PTD_TRUE_START}, reco start=${PTD_RECO_START}, Nreco=${PTD_RECO_N_BINS}, Ntrue=${PTD_TRUE_N_BINS}"

        root -l -b -q \
          "${MACHINE_MACRO}+( ${FONLL_JET}, ${CUT_NEG}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${OUT_DIR}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\" )" \
          2>&1 | awk '
            seen { print }
            /Loading RM from cache histograms\.\.\./ { seen=1; print }
          ' > "${OUT_DIR}/root.log"

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
          "${RUN_ID}" \
          "${PTD_EDGE_SHIFT}" \
          "${PTD_TRUE_START}" "${PTD_TRUE_N_BINS}" \
          "${PTD_RECO_START}" "${PTD_RECO_N_BINS}" \
          "${PTD_RECO_EDGES}" "${PTD_TRUE_EDGES}" \
          >> "${SUMMARY}"

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

echo "Hotovo. Summary: ${SUMMARY}"
