#!/usr/bin/env bash
set -euo pipefail

# =========================
# Uživatelské nastavení
# =========================
MACHINE_MACRO="./SuperIterace/Machine.C"
INPUT_FILE="./Data/Output_real_final_01022026.root"

# -------------------------
# z scan only
# -------------------------
Z_TRUE_START_LIST=(0 0.2 0.4)
Z_RECO_START_LIST=(0 0.2 0.4)
Z_RECO_ALLOWED_WIDTHS=(0.1 0.2 0.3)
Z_RECO_MIN_WIDTH=0.1
Z_RIGHT_EDGE=1.01
Z_LAST_BIN_LEFT=0.9   # poslední bin je vždy 0.9 -> 1.01

# -------------------------
# Fixed pT binning
# -------------------------
PT_TRUE_EDGES=(1 2 3 4 5 7 9 11 13 15 20)
PT_RECO_EDGES_C0=(1 2 3 4 5 7 9 11 13 16 20)
PT_RECO_EDGES_C1=(1 1.5 2 2.5 3 3.5 5 7 11 15 20)
PT_RECO_EDGES_C2=(1 1.5 2 2.5 3 3.5 4 5 9 14 20)

# -------------------------
# Fixed angular binning (all observables stay baseline)
# -------------------------
RECO_ANG_C0_A0=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
RECO_ANG_C0_A1=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
RECO_ANG_C0_A2=(0 0.025 0.05 0.075 0.1 0.15 0.2 0.3 0.4)
RECO_ANG_C0_A3=(0 0.025 0.05 0.075 0.1 0.15 0.2 0.3)
RECO_ANG_C0_A4=(0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8)
RECO_ANG_C0_A5=(0 0.30 0.50 0.65 0.75 0.85 1.01)

RECO_ANG_C1_A0=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
RECO_ANG_C1_A1=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
RECO_ANG_C1_A2=(0 0.025 0.05 0.075 0.1 0.15 0.2 0.3 0.4)
RECO_ANG_C1_A3=(0 0.025 0.05 0.075 0.1 0.15 0.2 0.3)
RECO_ANG_C1_A4=(0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8)
RECO_ANG_C1_A5=(0 0.30 0.50 0.65 0.75 0.85 1.01)

RECO_ANG_C2_A0=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
RECO_ANG_C2_A1=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
RECO_ANG_C2_A2=(0 0.025 0.05 0.075 0.1 0.15 0.2 0.3 0.4)
RECO_ANG_C2_A3=(0 0.025 0.05 0.075 0.1 0.15 0.2 0.3)
RECO_ANG_C2_A4=(0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8)
RECO_ANG_C2_A5=(0 0.30 0.50 0.65 0.75 0.85 1.01)

TRUE_ANG_C0_A0=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
TRUE_ANG_C0_A1=(0 0.025 0.05 0.075 0.1 0.125 0.15 0.2 0.3 0.6)
TRUE_ANG_C0_A2=(0 0.025 0.05 0.075 0.1 0.15 0.2 0.3 0.6)
TRUE_ANG_C0_A3=(0 0.0125 0.025 0.0375 0.05 0.075 0.1 0.15 0.4)
TRUE_ANG_C0_A4=(0 0.1 0.2 0.3 0.4 0.5 0.6 0.9)
TRUE_ANG_C0_A5=(0 0.30 0.50 0.65 0.75 0.85 1.01)

TRUE_ANG_C1_A0=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
TRUE_ANG_C1_A1=(0 0.025 0.05 0.075 0.1 0.125 0.15 0.2 0.3 0.6)
TRUE_ANG_C1_A2=(0 0.025 0.05 0.075 0.1 0.15 0.2 0.3 0.6)
TRUE_ANG_C1_A3=(0 0.0125 0.025 0.0375 0.05 0.075 0.1 0.15 0.4)
TRUE_ANG_C1_A4=(0 0.1 0.2 0.3 0.4 0.5 0.6 0.9)
TRUE_ANG_C1_A5=(0 0.30 0.50 0.65 0.75 0.85 1.01)

TRUE_ANG_C2_A0=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
TRUE_ANG_C2_A1=(0 0.025 0.05 0.075 0.1 0.125 0.15 0.2 0.3 0.6)
TRUE_ANG_C2_A2=(0 0.025 0.05 0.075 0.1 0.15 0.2 0.3 0.6)
TRUE_ANG_C2_A3=(0 0.0125 0.025 0.0375 0.05 0.075 0.1 0.15 0.4)
TRUE_ANG_C2_A4=(0 0.1 0.2 0.3 0.4 0.5 0.6 0.9)
TRUE_ANG_C2_A5=(0 0.30 0.50 0.65 0.75 0.85 1.01)

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
MAX_TOTAL_RUNS=0

# -------------------------
# Kam ukládat
# -------------------------
OUT_BASE="scan_z"
OVR_DIR="${OUT_BASE}/overrides"
RUN_DIR="${OUT_BASE}/runs"
SUMMARY="${OUT_BASE}/summary.tsv"
STABILITY="${OUT_BASE}/stability.tsv"

mkdir -p "${OVR_DIR}" "${RUN_DIR}"

# =========================
# Pomocné funkce
# =========================
join_by_comma_space() {
  local out=""
  local x
  for x in "$@"; do
    out+="${out:+, }${x}"
  done
  printf '%s' "${out}"
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
command -v python3 >/dev/null 2>&1 || { echo "python3 není v PATH"; exit 1; }

# =========================
# Hlavička summary
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\tz_true_start\tz_true_n_bins\tz_reco_start\tz_reco_n_bins\tz_reco_edges\tz_true_edges\n" > "${SUMMARY}"
fi

# =========================
# Resume počítadla
# =========================
run_counter=0
runs_started=0

# =========================
# Hlavní smyčka
# =========================
for Z_TRUE_START in "${Z_TRUE_START_LIST[@]}"; do
  Z_TRUE_EDGES=$(python3 - "${Z_TRUE_START}" <<'PY'
import sys
start = float(sys.argv[1])
base = [0,0.2,0.4,0.6,0.7,0.8,0.9,1.01]
edges = [x for x in base if x >= start - 1e-9]
out=[]
for x in edges:
    s=f"{x:.10f}".rstrip('0').rstrip('.')
    out.append(s if s else '0')
print(", ".join(out))
print(len(edges)-1, file=sys.stderr)
PY
 2> /tmp/z_true_nbins.$$)
  Z_TRUE_N_BINS=$(python3 - "${Z_TRUE_START}" <<'PY'
import sys
start = float(sys.argv[1])
base = [0,0.2,0.4,0.6,0.7,0.8,0.9,1.01]
edges = [x for x in base if x >= start - 1e-9]
print(len(edges)-1)
PY
)

  for Z_RECO_START in "${Z_RECO_START_LIST[@]}"; do
    # reco musi byt stejne nebo sirsi => start <= true_start
    if ! python3 - "$Z_RECO_START" "$Z_TRUE_START" <<'PY' >/dev/null
import sys
reco = float(sys.argv[1])
tru  = float(sys.argv[2])
raise SystemExit(0 if reco <= tru + 1e-9 else 1)
PY
    then
      continue
    fi

    while IFS= read -r Z_RECO_EDGES; do
      [[ -n "${Z_RECO_EDGES}" ]] || continue

      Z_RECO_N_BINS=$(awk -F',' '{print NF-1}' <<< "${Z_RECO_EDGES}")
      if (( Z_RECO_N_BINS < Z_TRUE_N_BINS )); then
        echo "[skip] z true start=${Z_TRUE_START}, reco start=${Z_RECO_START} -> Nreco=${Z_RECO_N_BINS} < Ntrue=${Z_TRUE_N_BINS}"
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
  // ---------- fixed pT truth ----------
  ptMcBinsVecCustom[0] = std::vector<double>{$(join_by_comma_space "${PT_TRUE_EDGES[@]}")};
  ptMcBinsVecCustom[1] = std::vector<double>{$(join_by_comma_space "${PT_TRUE_EDGES[@]}")};
  ptMcBinsVecCustom[2] = std::vector<double>{$(join_by_comma_space "${PT_TRUE_EDGES[@]}")};

  // ---------- fixed pT reco ----------
  ptRecoBinsVec[0] = std::vector<double>{$(join_by_comma_space "${PT_RECO_EDGES_C0[@]}")};
  ptRecoBinsVec[1] = std::vector<double>{$(join_by_comma_space "${PT_RECO_EDGES_C1[@]}")};
  ptRecoBinsVec[2] = std::vector<double>{$(join_by_comma_space "${PT_RECO_EDGES_C2[@]}")};

  // ---------- z scan, same for all centralities ----------
  zRecoBinsVec[0] = std::vector<double>{ ${Z_RECO_EDGES} };
  zRecoBinsVec[1] = std::vector<double>{ ${Z_RECO_EDGES} };
  zRecoBinsVec[2] = std::vector<double>{ ${Z_RECO_EDGES} };

  zMcBinsVecCustom[0] = std::vector<double>{ ${Z_TRUE_EDGES} };
  zMcBinsVecCustom[1] = std::vector<double>{ ${Z_TRUE_EDGES} };
  zMcBinsVecCustom[2] = std::vector<double>{ ${Z_TRUE_EDGES} };

  // ---------- fixed angular bins ----------
  angRecoBinsVec[0][0] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C0_A0[@]}")};
  angRecoBinsVec[0][1] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C0_A1[@]}")};
  angRecoBinsVec[0][2] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C0_A2[@]}")};
  angRecoBinsVec[0][3] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C0_A3[@]}")};
  angRecoBinsVec[0][4] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C0_A4[@]}")};
  angRecoBinsVec[0][5] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C0_A5[@]}")};

  angRecoBinsVec[1][0] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C1_A0[@]}")};
  angRecoBinsVec[1][1] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C1_A1[@]}")};
  angRecoBinsVec[1][2] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C1_A2[@]}")};
  angRecoBinsVec[1][3] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C1_A3[@]}")};
  angRecoBinsVec[1][4] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C1_A4[@]}")};
  angRecoBinsVec[1][5] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C1_A5[@]}")};

  angRecoBinsVec[2][0] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C2_A0[@]}")};
  angRecoBinsVec[2][1] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C2_A1[@]}")};
  angRecoBinsVec[2][2] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C2_A2[@]}")};
  angRecoBinsVec[2][3] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C2_A3[@]}")};
  angRecoBinsVec[2][4] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C2_A4[@]}")};
  angRecoBinsVec[2][5] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C2_A5[@]}")};

  angMcBinsVecCustom[0][0] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C0_A0[@]}")};
  angMcBinsVecCustom[0][1] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C0_A1[@]}")};
  angMcBinsVecCustom[0][2] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C0_A2[@]}")};
  angMcBinsVecCustom[0][3] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C0_A3[@]}")};
  angMcBinsVecCustom[0][4] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C0_A4[@]}")};
  angMcBinsVecCustom[0][5] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C0_A5[@]}")};

  angMcBinsVecCustom[1][0] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C1_A0[@]}")};
  angMcBinsVecCustom[1][1] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C1_A1[@]}")};
  angMcBinsVecCustom[1][2] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C1_A2[@]}")};
  angMcBinsVecCustom[1][3] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C1_A3[@]}")};
  angMcBinsVecCustom[1][4] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C1_A4[@]}")};
  angMcBinsVecCustom[1][5] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C1_A5[@]}")};

  angMcBinsVecCustom[2][0] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C2_A0[@]}")};
  angMcBinsVecCustom[2][1] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C2_A1[@]}")};
  angMcBinsVecCustom[2][2] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C2_A2[@]}")};
  angMcBinsVecCustom[2][3] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C2_A3[@]}")};
  angMcBinsVecCustom[2][4] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C2_A4[@]}")};
  angMcBinsVecCustom[2][5] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C2_A5[@]}")};
}
EOF2

      echo "[${RUN_ID}] z true start=${Z_TRUE_START}, reco start=${Z_RECO_START}, Nreco=${Z_RECO_N_BINS}, Ntrue=${Z_TRUE_N_BINS}"

      root -l -b -q \
        "${MACHINE_MACRO}+( ${FONLL_JET}, ${CUT_NEG}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${OUT_DIR}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\" )" \
        2>&1 | awk '
          seen { print }
          /Loading RM from cache histograms\.\.\./ { seen=1; print }
        ' > "${OUT_DIR}/root.log"

      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${RUN_ID}" \
        "${Z_TRUE_START}" "${Z_TRUE_N_BINS}" \
        "${Z_RECO_START}" "${Z_RECO_N_BINS}" \
        "${Z_RECO_EDGES}" "${Z_TRUE_EDGES}" \
        >> "${SUMMARY}"

    done < <(
      python3 - "${Z_RECO_START}" "${Z_RIGHT_EDGE}" "${Z_LAST_BIN_LEFT}" "${Z_RECO_MIN_WIDTH}" "${Z_TRUE_N_BINS}" <<'PY'
import sys

start = float(sys.argv[1])
right_edge = float(sys.argv[2])
last_left = float(sys.argv[3])
min_width = float(sys.argv[4])
# min_n_bins is not enforced here; shell checks against true later
allowed = [0.1, 0.2, 0.3]
EPS = 1e-9

if start > last_left - EPS:
    raise SystemExit(0)

results = set()

# fixed last bin 0.9 -> 1.01

def rec(right_edge_current, prev_width_effective, acc_edges_desc):
    # acc_edges_desc stores edges from right to left, including right_edge and current left edge
    current_left = acc_edges_desc[-1]
    rem = current_left - start

    # close with one transition bin on the far left
    if rem >= prev_width_effective - EPS and rem >= min_width - EPS:
        edges_desc = acc_edges_desc + [start]
        edges = list(reversed(edges_desc))
        if abs(edges[-1] - right_edge) < 1e-8 and abs(edges[0] - start) < 1e-8:
            results.add(tuple(round(x, 10) for x in edges))

    # add another regular bin on the left
    for w in allowed:
        if w + EPS < prev_width_effective:
            continue
        new_left = current_left - w
        if new_left <= start + EPS:
            continue
        rec(right_edge_current, w, acc_edges_desc + [new_left])

# start from fixed last bin
rec(right_edge, 0.1, [right_edge, last_left])

for edges in sorted(results):
    out=[]
    for x in edges:
        s=f"{x:.10f}".rstrip('0').rstrip('.')
        out.append(s if s else '0')
    print(", ".join(out))
PY
)

  done
done

echo "Hotovo. Summary: ${SUMMARY}"
