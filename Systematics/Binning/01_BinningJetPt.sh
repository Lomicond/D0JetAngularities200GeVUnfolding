#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Binning systematic scan: jet pT
#
# Run from the main project directory, i.e. the directory containing:
#   config.h, config_hist.h, Unfolding/Machine.C, Data/...
#
# Recommended run command:
#   bash Systematics/Binning/BinningJetPt.sh
# ============================================================

# -------------------------
# Working-directory checks
# -------------------------
PROJECT_DIR="$(pwd -P)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

[[ -f "${PROJECT_DIR}/config.h" ]]       || { echo "[error] Run this script from the main project directory: missing ./config.h"; exit 1; }
[[ -f "${PROJECT_DIR}/config_hist.h" ]]  || { echo "[error] Run this script from the main project directory: missing ./config_hist.h"; exit 1; }
[[ -f "${PROJECT_DIR}/Unfolding/Machine.C" ]] || { echo "[error] Missing ./Unfolding/Machine.C"; exit 1; }

command -v root    >/dev/null 2>&1 || { echo "[error] root is not in PATH"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "[error] python3 is not in PATH"; exit 1; }

# =========================
# User settings
# =========================
MACHINE_MACRO="${PROJECT_DIR}/Unfolding/Machine.C"
INPUT_FILE="${PROJECT_DIR}/Data/Output_real_final_01022026.root"

[[ -f "${INPUT_FILE}" ]] || { echo "[error] Missing input file: ${INPUT_FILE}"; exit 1; }

# -------------------------
# Reco-level jet-pT scan
# -------------------------
RECO_PT_MIN_LIST=(1 2 3 4 5)
RECO_PT_MAX_LIST=(20 25)

MIN_WIDTH=0.5
MAX_WIDTH=10
WIDTH_TREND="inc"            # "inc" or "dec"
STEP=""                      # empty => use MIN_WIDTH

# Automatic reco-bin count based on the interval length
AUTO_N_BINS=1
FIXED_N_BINS=12              # used only if AUTO_N_BINS=0
TARGET_AVG_WIDTH=2.0         # target average reco-bin width
N_BINS_OFFSETS=(0)           # e.g. (-1 0 1), but for the final scan it is better to control this explicitly
MIN_N_BINS_SCAN=5
MAX_N_BINS_SCAN=12

# Limit on the number of different reco binnings for one configuration
# (RECO_PT_MIN, RECO_PT_MAX, RECO_N_BINS). 0 = no limit.
MAX_RECO_PATTERNS_PER_CONFIG=0

# Edges above BREAK_PT must be integer values.
# This allows finer bins at low pT, but only integer edges at higher pT.
BREAK_PT=5.0

# -------------------------
# True-level jet-pT scan
# Only the lower and upper endpoints are varied;
# the internal edges are a subset of BASE_TRUE_EDGES.
# -------------------------
TRUE_PT_MIN_LIST=(1 3 5)
TRUE_PT_MAX_LIST=(20 25)
BASE_TRUE_EDGES=(1 2 3 4 5 7 9 11 13 15 20 25)

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
MIN_PT_D0=5
MAX_PT_D0=10
USE_PRIOR_SHAPE_WEIGHTING=0
SYSTEMATIC_SPLOT=0

# Safety limit for the whole scan. 0 = no limit.
MAX_TOTAL_RUNS=0

# -------------------------
# Output location
# -------------------------
SCAN_DIR="${SCRIPT_DIR}/scanJetPt"
OVR_DIR="${SCAN_DIR}/overrides"
RUN_DIR="${SCAN_DIR}/runs"
SUMMARY="${SCAN_DIR}/summary.tsv"

mkdir -p "${OVR_DIR}" "${RUN_DIR}" "${PROJECT_DIR}/OutputPdf"

# =========================
# Helper functions
# =========================
join_by_comma_space() {
  local out=""
  local x
  for x in "$@"; do
    out+="${out:+, }${x}"
  done
  printf '%s' "${out}"
}

build_true_edges() {
  local tmin="$1"
  local tmax="$2"

  local selected=()
  local e
  for e in "${BASE_TRUE_EDGES[@]}"; do
    if (( e >= tmin && e <= tmax )); then
      selected+=("${e}")
    fi
  done

  if (( ${#selected[@]} < 2 )); then
    return 1
  fi

  local last_idx=$(( ${#selected[@]} - 1 ))
  if [[ "${selected[0]}" != "${tmin}" ]]; then
    return 1
  fi
  if [[ "${selected[$last_idx]}" != "${tmax}" ]]; then
    return 1
  fi

  TRUE_N_BINS=$(( ${#selected[@]} - 1 ))
  TRUE_EDGES_STR="$(join_by_comma_space "${selected[@]}")"
  return 0
}

# STEP default
if [[ -z "${STEP}" ]]; then
  STEP="${MIN_WIDTH}"
fi

# =========================
# Summary header
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\treco_pt_min\treco_pt_max\treco_n_bins\ttrue_pt_min\ttrue_pt_max\ttrue_n_bins\tmin_width\tmax_width\tstep\ttrend\treco_edges\ttrue_edges\toverride_file\troot_log\n" > "${SUMMARY}"
fi

# =========================
# Determine the last run ID
# =========================
run_counter=0
shopt -s nullglob
for d in "${RUN_DIR}"/r*; do
  [[ -d "${d}" ]] || continue
  bn=$(basename "${d}")
  if [[ "${bn}" =~ ^r([0-9]+)$ ]]; then
    num=$((10#${BASH_REMATCH[1]}))
    (( num > run_counter )) && run_counter=${num}
  fi
done
shopt -u nullglob

runs_started=0

# =========================
# Main loop
# =========================
for RECO_PT_MIN in "${RECO_PT_MIN_LIST[@]}"; do
  for RECO_PT_MAX in "${RECO_PT_MAX_LIST[@]}"; do

    if (( RECO_PT_MAX <= RECO_PT_MIN )); then
      continue
    fi

    # -------------------------
    # Candidate reco-bin counts
    # -------------------------
    RECO_N_BINS_LIST=()

    if (( AUTO_N_BINS )); then
      NOMINAL_RECO_N_BINS=$(python3 - "${RECO_PT_MIN}" "${RECO_PT_MAX}" "${TARGET_AVG_WIDTH}" <<'PY'
import sys
pt_min = float(sys.argv[1])
pt_max = float(sys.argv[2])
target = float(sys.argv[3])
print(int(round((pt_max - pt_min) / target)))
PY
)

      declare -A seen_nbins=()
      for OFF in "${N_BINS_OFFSETS[@]}"; do
        NB=$((NOMINAL_RECO_N_BINS + OFF))
        (( NB < MIN_N_BINS_SCAN )) && continue
        (( NB > MAX_N_BINS_SCAN )) && continue
        if [[ -z "${seen_nbins[$NB]+x}" ]]; then
          RECO_N_BINS_LIST+=("${NB}")
          seen_nbins[$NB]=1
        fi
      done
      unset seen_nbins
    else
      RECO_N_BINS_LIST=("${FIXED_N_BINS}")
    fi

    (( ${#RECO_N_BINS_LIST[@]} > 0 )) || continue

    for RECO_N_BINS in "${RECO_N_BINS_LIST[@]}"; do

      while IFS= read -r RECO_EDGES; do
        [[ -n "${RECO_EDGES}" ]] || continue

        for TRUE_PT_MIN in "${TRUE_PT_MIN_LIST[@]}"; do
          for TRUE_PT_MAX in "${TRUE_PT_MAX_LIST[@]}"; do

            if (( TRUE_PT_MAX <= TRUE_PT_MIN )); then
              continue
            fi

            if ! build_true_edges "${TRUE_PT_MIN}" "${TRUE_PT_MAX}"; then
              continue
            fi

            # The true level must not have more bins than the reco level.
            if (( TRUE_N_BINS > RECO_N_BINS )); then
              echo "[skip] reco ${RECO_PT_MIN}-${RECO_PT_MAX} (N=${RECO_N_BINS}), true ${TRUE_PT_MIN}-${TRUE_PT_MAX} (N=${TRUE_N_BINS}) => N_true > N_reco"
              continue
            fi

            if (( MAX_TOTAL_RUNS > 0 && runs_started >= MAX_TOTAL_RUNS )); then
              echo "Reached MAX_TOTAL_RUNS=${MAX_TOTAL_RUNS}, stopping."
              exit 0
            fi

            ((run_counter += 1))
            ((runs_started += 1))
            printf -v RUN_ID "r%06d" "${run_counter}"

            OVR_FILE="${OVR_DIR}/override_${RUN_ID}.C"
            OUT_DIR="${RUN_DIR}/${RUN_ID}"
            ROOT_LOG="${OUT_DIR}/root.log"
            mkdir -p "${OUT_DIR}"

            # -------------------------
            # Override macro for Machine.C
            # Overrides the jet-pT reco and true binning for all centralities.
            # -------------------------
            cat > "${OVR_FILE}" <<EOF_OVR
{
  cout << "[override] jet-pT reco edges: ${RECO_EDGES}" << endl;
  cout << "[override] jet-pT true edges: ${TRUE_EDGES_STR}" << endl;

  for (int ic = 0; ic < nCentralityBins; ++ic) ptRecoBinsVec[ic].clear();
  double reco_edges[] = { ${RECO_EDGES} };
  int nReco = (int)(sizeof(reco_edges)/sizeof(double));
  for (int ic = 0; ic < nCentralityBins; ++ic) {
    for (int i = 0; i < nReco; ++i) ptRecoBinsVec[ic].push_back(reco_edges[i]);
  }

  for (int ic = 0; ic < nCentralityBins; ++ic) ptMcBinsVecCustom[ic].clear();
  double true_edges[] = { ${TRUE_EDGES_STR} };
  int nTrue = (int)(sizeof(true_edges)/sizeof(double));
  for (int ic = 0; ic < nCentralityBins; ++ic) {
    for (int i = 0; i < nTrue; ++i) ptMcBinsVecCustom[ic].push_back(true_edges[i]);
  }
}
EOF_OVR

            printf "%s\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
              "${RUN_ID}" \
              "${RECO_PT_MIN}" "${RECO_PT_MAX}" "${RECO_N_BINS}" \
              "${TRUE_PT_MIN}" "${TRUE_PT_MAX}" "${TRUE_N_BINS}" \
              "${MIN_WIDTH}" "${MAX_WIDTH}" "${STEP}" "${WIDTH_TREND}" \
              "${RECO_EDGES}" "${TRUE_EDGES_STR}" "${OVR_FILE}" "${ROOT_LOG}" \
              >> "${SUMMARY}"

            echo "[${RUN_ID}] reco: ${RECO_PT_MIN}-${RECO_PT_MAX}, N=${RECO_N_BINS} | true: ${TRUE_PT_MIN}-${TRUE_PT_MAX}, N=${TRUE_N_BINS}"

            # -------------------------
            # Run the ROOT macro
            # OutputFile = RUN_ID because Machine.C uses runId = BaseName(OutputFile)
            # ScanDir    = SCAN_DIR so stability.tsv and Output/OutputSpectra*.root go there
            # -------------------------
            if ! root -l -b -q \
              "${MACHINE_MACRO}+( ${FONLL_JET}, ${CUT_NEG}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${RUN_ID}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\", \"${SCAN_DIR}\", ${USE_PRIOR_SHAPE_WEIGHTING}, ${SYSTEMATIC_SPLOT} )" \
              > "${ROOT_LOG}" 2>&1; then
              echo "[error] ROOT failed for ${RUN_ID}. Last 60 lines of the log:"
              tail -n 60 "${ROOT_LOG}" || true
              exit 1
            fi

          done
        done

      done < <(
        python3 - "${RECO_PT_MIN}" "${RECO_PT_MAX}" "${RECO_N_BINS}" "${MIN_WIDTH}" "${MAX_WIDTH}" "${WIDTH_TREND}" "${STEP}" "${MAX_RECO_PATTERNS_PER_CONFIG}" "${BREAK_PT}" <<'PY'
import sys, math

pt_min        = float(sys.argv[1])
pt_max        = float(sys.argv[2])
n_bins        = int(sys.argv[3])
min_width     = float(sys.argv[4])
max_width_arg = sys.argv[5].strip()
trend         = sys.argv[6].strip().lower()
step          = float(sys.argv[7])
max_runs      = int(sys.argv[8])
break_pt      = float(sys.argv[9])

EPS = 1e-9

def skip(msg):
    print(f"[skip] reco config PT_MIN={pt_min}, PT_MAX={pt_max}, N_BINS={n_bins}: {msg}", file=sys.stderr)
    raise SystemExit(0)

if trend not in ("inc", "dec"):
    skip("WIDTH_TREND must be 'inc' or 'dec'")

L = pt_max - pt_min
if L <= 0:
    skip("PT_MAX must be > PT_MIN")

def close_to_int(x, eps=1e-9):
    return abs(x - round(x)) < eps

total_u = L / step
if not close_to_int(total_u):
    skip(f"range {L} is not divisible by STEP={step}")
total_u = int(round(total_u))

min_u = int(math.ceil(min_width / step - 1e-12))
if min_u <= 0:
    skip("MIN_WIDTH converted to STEP units is invalid")

max_u = None
if max_width_arg not in ("", "0", "0.0"):
    max_width = float(max_width_arg)
    if max_width < min_width:
        skip("MAX_WIDTH must be >= MIN_WIDTH")
    max_u = int(math.floor(max_width / step + 1e-12))
    if max_u < min_u:
        skip("MAX_WIDTH converted to STEP units is < MIN_WIDTH")

if n_bins * min_u > total_u:
    skip("N_BINS * MIN_WIDTH is larger than the interval")

if max_u is not None and n_bins * max_u < total_u:
    skip("N_BINS * MAX_WIDTH is smaller than the interval")

R = total_u - n_bins * min_u
count = 0

def gen_extra_inc(remaining, k, prev_e):
    if k == 0:
        if remaining == 0:
            yield []
        return
    min_rest = (k - 1) * prev_e
    max_e = remaining - min_rest
    for e in range(prev_e, max_e + 1):
        for tail in gen_extra_inc(remaining - e, k - 1, e):
            yield [e] + tail

def gen_extra_dec(remaining, k, prev_e):
    if k == 0:
        if remaining == 0:
            yield []
        return
    hi = min(prev_e, remaining)
    for e in range(hi, -1, -1):
        for tail in gen_extra_dec(remaining - e, k - 1, e):
            yield [e] + tail

gen = gen_extra_inc(R, n_bins, 0) if trend == "inc" else gen_extra_dec(R, n_bins, R)

for extra in gen:
    widths_u = [min_u + e for e in extra]

    if trend == "inc":
        ok = all(widths_u[i] <= widths_u[i+1] for i in range(n_bins - 1))
    else:
        ok = all(widths_u[i] >= widths_u[i+1] for i in range(n_bins - 1))
    if not ok:
        continue

    edges = [pt_min]
    acc_u = 0
    for w in widths_u:
        acc_u += w
        edges.append(pt_min + acc_u * step)
    edges[-1] = pt_max

    # Above BREAK_PT, allow only integer edges.
    if any((x > break_pt + EPS) and (abs(x - round(x)) > EPS) for x in edges):
        continue

    if max_u is not None and any(w > max_u for w in widths_u):
        continue

    if sum(widths_u) != total_u:
        continue

    nd = max(0, int(round(-math.log10(step))) + 2) if step < 1 else 6
    fmt = "{:" + f".{nd}f" + "}"

    out = []
    for x in edges:
        s = fmt.format(x)
        if "." in s:
            s = s.rstrip("0").rstrip(".")
        out.append(s)

    print(", ".join(out))
    count += 1

    if max_runs > 0 and count >= max_runs:
        break
PY
      )

    done
  done
done

echo "Done."
echo "Summary:   ${SUMMARY}"
echo "Stability: ${SCAN_DIR}/stability.tsv"
echo "Spectra:   ${SCAN_DIR}/Output/"
