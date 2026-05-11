#!/usr/bin/env bash
set -euo pipefail

# =========================
# Uživatelské nastavení
# =========================
MACHINE_MACRO="./SuperIterace/Machine.C"
INPUT_FILE="./Data/Output_real_final_01022026.root"

# -------------------------
# L12 scan only
# -------------------------
L12_INDEX=2
L12_TRUE_MAX_LIST=(0.4 0.6 0.8 1.0)
L12_RECO_MAX_LIST=(0.4 0.6 0.8 1.0)
L12_DENSE_BIN_COUNTS=(4 5 6)
L12_MIN_WIDTH=0.025
L12_TAIL_WIDTHS=(0.1 0.2 0.3 0.4)
L12_TRUE_N_BINS=8

# -------------------------
# Fixed pT binning
# -------------------------
PT_TRUE_EDGES=(1 2 3 4 5 7 9 11 13 15 20)
PT_RECO_EDGES_C0=(1 2 3 4 5 7 9 11 13 16 20)
PT_RECO_EDGES_C1=(1 1.5 2 2.5 3 3.5 5 7 11 15 20)
PT_RECO_EDGES_C2=(1 1.5 2 2.5 3 3.5 4 5 9 14 20)

# -------------------------
# Fixed angular binning for all observables except l12 (index 2)
# -------------------------
# Reco baseline (3 centralities x 6 observables)
RECO_ANG_C0_A0=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
RECO_ANG_C0_A1=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
RECO_ANG_C0_A3=(0 0.025 0.05 0.075 0.1 0.15 0.2 0.3)
RECO_ANG_C0_A4=(0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8)
RECO_ANG_C0_A5=(0 0.30 0.50 0.65 0.75 0.85 1.01)

RECO_ANG_C1_A0=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
RECO_ANG_C1_A1=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
RECO_ANG_C1_A3=(0 0.025 0.05 0.075 0.1 0.15 0.2 0.3)
RECO_ANG_C1_A4=(0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8)
RECO_ANG_C1_A5=(0 0.30 0.50 0.65 0.75 0.85 1.01)

RECO_ANG_C2_A0=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
RECO_ANG_C2_A1=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
RECO_ANG_C2_A3=(0 0.025 0.05 0.075 0.1 0.15 0.2 0.3)
RECO_ANG_C2_A4=(0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8)
RECO_ANG_C2_A5=(0 0.30 0.50 0.65 0.75 0.85 1.01)

# True baseline (3 centralities x 6 observables)
TRUE_ANG_C0_A0=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
TRUE_ANG_C0_A1=(0 0.025 0.05 0.075 0.1 0.125 0.15 0.2 0.3 0.6)
TRUE_ANG_C0_A3=(0 0.0125 0.025 0.0375 0.05 0.075 0.1 0.15 0.4)
TRUE_ANG_C0_A4=(0 0.1 0.2 0.3 0.4 0.5 0.6 0.9)
TRUE_ANG_C0_A5=(0 0.30 0.50 0.65 0.75 0.85 1.01)

TRUE_ANG_C1_A0=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
TRUE_ANG_C1_A1=(0 0.025 0.05 0.075 0.1 0.125 0.15 0.2 0.3 0.6)
TRUE_ANG_C1_A3=(0 0.0125 0.025 0.0375 0.05 0.075 0.1 0.15 0.4)
TRUE_ANG_C1_A4=(0 0.1 0.2 0.3 0.4 0.5 0.6 0.9)
TRUE_ANG_C1_A5=(0 0.30 0.50 0.65 0.75 0.85 1.01)

TRUE_ANG_C2_A0=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.4)
TRUE_ANG_C2_A1=(0 0.025 0.05 0.075 0.1 0.125 0.15 0.2 0.3 0.6)
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
OUT_BASE="scan_l12"
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
  printf "run_id\tl12_reco_max\tl12_dense_bins\tl12_reco_n_bins\tl12_true_max\tl12_true_n_bins\tl12_reco_edges\tl12_true_edges\n" > "${SUMMARY}"
fi

# =========================
# Resume počítadla
# =========================
run_counter=0
runs_started=0

# =========================
# Hlavní smyčka
# =========================
for L12_RECO_MAX in "${L12_RECO_MAX_LIST[@]}"; do
  for L12_DENSE_BINS in "${L12_DENSE_BIN_COUNTS[@]}"; do

    while IFS= read -r L12_RECO_EDGES; do
      [[ -n "${L12_RECO_EDGES}" ]] || continue

      L12_RECO_N_BINS=$(awk -F',' '{print NF-1}' <<< "${L12_RECO_EDGES}")

      for L12_TRUE_MAX in "${L12_TRUE_MAX_LIST[@]}"; do
        L12_TRUE_EDGES=$(python3 - "${L12_TRUE_MAX}" <<'PY'
import sys
xmax = float(sys.argv[1])
edges = [0,0.025,0.05,0.075,0.1,0.15,0.2,0.3,xmax]
print(", ".join(str(x).rstrip('0').rstrip('.') if '.' in str(x) else str(x) for x in edges))
PY
)
        L12_TRUE_N_BINS=8

        if (( L12_RECO_N_BINS < L12_TRUE_N_BINS )); then
          echo "[skip] l12 reco max=${L12_RECO_MAX}, dense=${L12_DENSE_BINS} -> Nreco=${L12_RECO_N_BINS} < Ntrue=${L12_TRUE_N_BINS}"
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

  // ---------- l12 scan (index 2), same for all centralities ----------
  angRecoBinsVec[0][${L12_INDEX}] = std::vector<double>{ ${L12_RECO_EDGES} };
  angRecoBinsVec[1][${L12_INDEX}] = std::vector<double>{ ${L12_RECO_EDGES} };
  angRecoBinsVec[2][${L12_INDEX}] = std::vector<double>{ ${L12_RECO_EDGES} };

  angMcBinsVecCustom[0][${L12_INDEX}] = std::vector<double>{ ${L12_TRUE_EDGES} };
  angMcBinsVecCustom[1][${L12_INDEX}] = std::vector<double>{ ${L12_TRUE_EDGES} };
  angMcBinsVecCustom[2][${L12_INDEX}] = std::vector<double>{ ${L12_TRUE_EDGES} };

  // ---------- fixed reco angular bins for remaining observables ----------
  angRecoBinsVec[0][0] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C0_A0[@]}")};
  angRecoBinsVec[0][1] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C0_A1[@]}")};
  angRecoBinsVec[0][3] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C0_A3[@]}")};
  angRecoBinsVec[0][4] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C0_A4[@]}")};
  angRecoBinsVec[0][5] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C0_A5[@]}")};

  angRecoBinsVec[1][0] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C1_A0[@]}")};
  angRecoBinsVec[1][1] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C1_A1[@]}")};
  angRecoBinsVec[1][3] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C1_A3[@]}")};
  angRecoBinsVec[1][4] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C1_A4[@]}")};
  angRecoBinsVec[1][5] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C1_A5[@]}")};

  angRecoBinsVec[2][0] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C2_A0[@]}")};
  angRecoBinsVec[2][1] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C2_A1[@]}")};
  angRecoBinsVec[2][3] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C2_A3[@]}")};
  angRecoBinsVec[2][4] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C2_A4[@]}")};
  angRecoBinsVec[2][5] = std::vector<double>{$(join_by_comma_space "${RECO_ANG_C2_A5[@]}")};

  // ---------- fixed true angular bins for remaining observables ----------
  angMcBinsVecCustom[0][0] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C0_A0[@]}")};
  angMcBinsVecCustom[0][1] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C0_A1[@]}")};
  angMcBinsVecCustom[0][3] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C0_A3[@]}")};
  angMcBinsVecCustom[0][4] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C0_A4[@]}")};
  angMcBinsVecCustom[0][5] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C0_A5[@]}")};

  angMcBinsVecCustom[1][0] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C1_A0[@]}")};
  angMcBinsVecCustom[1][1] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C1_A1[@]}")};
  angMcBinsVecCustom[1][3] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C1_A3[@]}")};
  angMcBinsVecCustom[1][4] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C1_A4[@]}")};
  angMcBinsVecCustom[1][5] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C1_A5[@]}")};

  angMcBinsVecCustom[2][0] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C2_A0[@]}")};
  angMcBinsVecCustom[2][1] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C2_A1[@]}")};
  angMcBinsVecCustom[2][3] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C2_A3[@]}")};
  angMcBinsVecCustom[2][4] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C2_A4[@]}")};
  angMcBinsVecCustom[2][5] = std::vector<double>{$(join_by_comma_space "${TRUE_ANG_C2_A5[@]}")};
}
EOF2

        echo "[${RUN_ID}] l12 reco max=${L12_RECO_MAX}, dense=${L12_DENSE_BINS}, Nreco=${L12_RECO_N_BINS} | l12 true max=${L12_TRUE_MAX}"

        root -l -b -q \
          "${MACHINE_MACRO}+( ${FONLL_JET}, ${CUT_NEG}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${OUT_DIR}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\" )" \
          2>&1 | awk '
            seen { print }
            /Loading RM from cache histograms\.\.\./ { seen=1; print }
          ' > "${OUT_DIR}/root.log"

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
          "${RUN_ID}" \
          "${L12_RECO_MAX}" "${L12_DENSE_BINS}" "${L12_RECO_N_BINS}" \
          "${L12_TRUE_MAX}" "${L12_TRUE_N_BINS}" \
          "${L12_RECO_EDGES}" "${L12_TRUE_EDGES}" \
          >> "${SUMMARY}"

      done
    done < <(
      python3 - "${L12_RECO_MAX}" "${L12_DENSE_BINS}" "${L12_MIN_WIDTH}" "${L12_TRUE_N_BINS}" <<'PY'
import sys

xmax = float(sys.argv[1])
dense_bins = int(sys.argv[2])
min_width = float(sys.argv[3])
min_n_bins = int(sys.argv[4])
allowed = [0.1, 0.2, 0.3, 0.4]
EPS = 1e-9

prefix = [0.0]
for _ in range(dense_bins):
    prefix.append(prefix[-1] + min_width)

fixed_sum = prefix[-1]
remaining = xmax - fixed_sum
if remaining < 0.1 - EPS:
    raise SystemExit(0)

results = set()

def rec(prev, used, widths):
    tail_count_so_far = len(widths)
    total_bins_so_far = dense_bins + tail_count_so_far + 1
    rem = remaining - used

    if total_bins_so_far >= min_n_bins:
        if rem + EPS >= max(prev, 0.1):
            tail = widths + [round(rem, 10)]
            edges = prefix[:]
            acc = fixed_sum
            for w in tail:
                acc += w
                edges.append(round(acc, 10))
            if abs(edges[-1] - xmax) < 1e-8:
                results.add(tuple(edges))

    for w in allowed:
        if w + EPS < max(prev, 0.1):
            continue
        rem_after = rem - w
        if rem_after < max(w, 0.1) - EPS:
            continue
        rec(w, used + w, widths + [w])

rec(0.1, 0.0, [])

for edges in sorted(results):
    out = []
    for x in edges:
        s = f"{x:.10f}".rstrip('0').rstrip('.')
        out.append(s if s else '0')
    print(", ".join(out))
PY
)

  done
done

echo "Hotovo. Summary: ${SUMMARY}"
