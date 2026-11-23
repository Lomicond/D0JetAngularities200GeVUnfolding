#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Wide / diagnostic jet-pT binning scan
#
# Usage:
#   ./Systematics/Binning/01_BinningJetPtWide.sh 1 5
#
# Purpose:
#   Search a substantially wider set of jet-pT binnings than the final
#   binning-systematic scan, mainly to diagnose iteration stability.
#
# Philosophy:
#   - Machine.C receives MIN_PT_D0 and MAX_PT_D0 from the command line;
#   - every candidate jet-pT binning is applied to ALL centralities in a run,
#     so 0-10%, 10-40%, and 40-80% can be compared under exactly the same
#     response binning;
#   - reco candidates are generated from ALL THREE nominal centrality styles,
#     plus several explicit hybrid/coarse families;
#   - true candidates are controlled coarsenings of the nominal true binning;
#   - reco and true families are crossed, with the safety requirement
#       N_true_bins <= N_reco_bins;
#   - only the cached response is used: no TTree response rebuilding;
#   - this is a DIAGNOSTIC scan, not automatically part of the systematic
#     combination.
#
# Current plan:
#   27 unique reco candidates
#    8 true candidates
#   170 valid reco x true combinations after N_true <= N_reco.
#
# The scan is deterministic, so resume by run_id is safe.
#
# IMPORTANT:
#   Output goes to scanJetPtWide_MIN_MAX, NOT scanJetPt_MIN_MAX,
#   so the previous systematic scan is not overwritten.
# ============================================================

# -------------------------
# D0 pT range for Machine()
# -------------------------
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 MIN_PT_D0 MAX_PT_D0"
  echo "Example: $0 1 5"
  exit 1
fi

MIN_PT_D0="$1"
MAX_PT_D0="$2"

if ! [[ "${MIN_PT_D0}" =~ ^[0-9]+$ && "${MAX_PT_D0}" =~ ^[0-9]+$ ]]; then
  echo "[error] MIN_PT_D0 and MAX_PT_D0 must be integers."
  exit 1
fi

if (( MIN_PT_D0 >= MAX_PT_D0 )); then
  echo "[error] MIN_PT_D0 must be smaller than MAX_PT_D0."
  exit 1
fi

# -------------------------
# Working-directory checks
# -------------------------
PROJECT_DIR="$(pwd -P)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

[[ -f "${PROJECT_DIR}/config.h" ]]            || { echo "[error] Run this script from the main project directory: missing ./config.h"; exit 1; }
[[ -f "${PROJECT_DIR}/config_hist.h" ]]       || { echo "[error] Run this script from the main project directory: missing ./config_hist.h"; exit 1; }
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
# Machine() parameters
# -------------------------
# Machine(fonllJet, CutOfNegative, minJetPtRecoCut, savedIter,
#         InputFile, OutputFile, minPtD0Cut, maxPtD0Cut,
#         OverrideMacro, ScanDir, usePriorShapeWeighting, systematicSPlot)
FONLL_JET=1
CUT_NEG=1
MIN_JET_PT_RECO_CUT=-30
SAVED_ITER=4
USE_PRIOR_SHAPE_WEIGHTING=0
SYSTEMATIC_SPLOT=0

# Safety limit for a quick test. 0 = no limit.
# Example:
#   MAX_TOTAL_RUNS=5
MAX_TOTAL_RUNS=0

# This is only metadata for summary.tsv.  All candidate edges themselves
# are taken from cache-compatible half-GeV / existing bin boundaries.
EDGE_STEP=0.5

# -------------------------
# Output location
# -------------------------
SCAN_DIR="${SCRIPT_DIR}/scanJetPtWide_${MIN_PT_D0}_${MAX_PT_D0}"
OVR_DIR="${SCAN_DIR}/overrides"
RUN_DIR="${SCAN_DIR}/runs"
SUMMARY="${SCAN_DIR}/summary.tsv"
STABILITY="${SCAN_DIR}/stability.tsv"
PLAN="${SCAN_DIR}/planned_combinations.tsv"

mkdir -p \
  "${OVR_DIR}" \
  "${RUN_DIR}" \
  "${SCAN_DIR}/Output" \
  "${PROJECT_DIR}/OutputPdf"

# =========================
# Helper functions
# =========================
edge_stats() {
  local edges_in="$1"

  python3 - "${edges_in}" <<'PY'
import sys

edges = [float(x.strip()) for x in sys.argv[1].split(",") if x.strip()]
if len(edges) < 2:
    raise SystemExit("Need at least two edges")

widths = [b - a for a, b in zip(edges, edges[1:])]
if any(w <= 0 for w in widths):
    raise SystemExit(f"Non-increasing edges: {edges}")

def fmt(x):
    s = f"{x:.10f}".rstrip("0").rstrip(".")
    return s if s else "0"

print(
    fmt(edges[0]),
    fmt(edges[-1]),
    len(edges) - 1,
    fmt(min(widths)),
    fmt(max(widths)),
)
PY
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
# Build a complete deterministic plan
# =========================
#
# Reco reference styles are the three nominal centrality-dependent binnings:
#
#   0-10%:
#     1,2,3,4,5,6,7,8,9,10,13,16,25
#
#   10-40%:
#     1,1.5,2,2.5,3,3.5,4.5,7,10,13,16,19,25
#
#   40-80%:
#     1,1.5,2,2.5,3,4,5,7,9,11,13,18,25
#
# From EACH reference we generate:
#   nominal
#   low_x2 / low_x3    : coarsen only pT <= 5
#   high_x2 / high_x3  : coarsen only pT > 5
#   all_x2 / all_x3    : coarsen the whole internal edge sequence
#
# Then several explicit hybrid families are added.
#
# True-level families are all controlled coarsenings of:
#   1,2,3,4,5,7,9,11,13,15,20
#
# No "seveleni" is included here on purpose.  This scan is intended to answer
# whether COARSENING can stabilize the unfolding.  Edge nudges are already
# covered by the dedicated binning-systematic scan.
#
python3 - "${PLAN}" <<'PY'
import sys
from collections import OrderedDict

out_path = sys.argv[1]

reco_bases = OrderedDict([
    ("ref_0_10",  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 16, 25]),
    ("ref_10_40", [1, 1.5, 2, 2.5, 3, 3.5, 4.5, 7, 10, 13, 16, 19, 25]),
    ("ref_40_80", [1, 1.5, 2, 2.5, 3, 4, 5, 7, 9, 11, 13, 18, 25]),
])


def fmt_num(x):
    x = float(x)
    if x.is_integer():
        return str(int(x))
    return f"{x:.10f}".rstrip("0").rstrip(".")


def fmt_edges(edges):
    return ", ".join(fmt_num(x) for x in edges)


def coarsen_region(edges, region, stride):
    """
    Keep endpoints fixed.

    region="low":
        decimate internal edges <= 5; keep >5 unchanged.
    region="high":
        keep <=5 unchanged; decimate >5.
    region="all":
        decimate the complete internal sequence.

    The last edge of the decimated region is explicitly retained so the
    transition into the untouched region does not become accidental.
    """
    edges = list(map(float, edges))
    first, last = edges[0], edges[-1]
    internal = edges[1:-1]

    if region == "low":
        a = [x for x in internal if x <= 5.0]
        b = [x for x in internal if x > 5.0]
        kept = a[::stride]
        if a and a[-1] not in kept:
            kept.append(a[-1])
        result = [first] + kept + b + [last]

    elif region == "high":
        a = [x for x in internal if x <= 5.0]
        b = [x for x in internal if x > 5.0]
        kept = b[::stride]
        if b and b[-1] not in kept:
            kept.append(b[-1])
        result = [first] + a + kept + [last]

    elif region == "all":
        kept = internal[::stride]
        if internal and internal[-1] not in kept:
            kept.append(internal[-1])
        result = [first] + kept + [last]

    else:
        raise ValueError(region)

    # Deduplicate and sort defensively.
    return sorted(set(result))


# ------------------------------------------------------------
# Reco candidates
# ------------------------------------------------------------
# Use edge tuple as key to remove duplicate binnings generated from
# different source styles/families. Keep the first provenance label.
reco_candidates = OrderedDict()

for source, base in reco_bases.items():
    families = OrderedDict()
    families["nominal"] = base

    for stride in (2, 3):
        families[f"low_x{stride}"]  = coarsen_region(base, "low",  stride)
        families[f"high_x{stride}"] = coarsen_region(base, "high", stride)
        families[f"all_x{stride}"]  = coarsen_region(base, "all",  stride)

    for family, edges in families.items():
        reco_candidates.setdefault(
            tuple(map(float, edges)),
            (source, family, edges),
        )

# Explicit physically simple hybrid/coarse families.
manual_reco = OrderedDict([
    (
        "common_coarse_A",
        [1, 2, 3, 4, 5, 7, 9, 13, 18, 25],
    ),
    (
        "common_coarse_B",
        [1, 2, 3, 5, 7, 10, 13, 18, 25],
    ),
    (
        "common_coarse_C",
        [1, 2, 3, 5, 7, 10, 16, 25],
    ),
    (
        "very_coarse",
        [1, 2, 3, 5, 7, 10, 15, 25],
    ),
    (
        "lowfine_tailcoarse",
        [1, 1.5, 2, 2.5, 3, 4, 5, 7, 10, 15, 25],
    ),
    (
        "lowcoarse_tailfine",
        [1, 2, 3, 4, 5, 7, 9, 11, 13, 16, 25],
    ),
])

for family, edges in manual_reco.items():
    reco_candidates.setdefault(
        tuple(map(float, edges)),
        ("manual", family, edges),
    )


# ------------------------------------------------------------
# True candidates
# ------------------------------------------------------------
true_candidates = OrderedDict([
    (
        "true_nominal",
        [1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 20],
    ),
    (
        "true_low_coarse_A",
        [1, 2, 3, 5, 7, 9, 11, 13, 15, 20],
    ),
    (
        "true_low_coarse_B",
        [1, 3, 5, 7, 9, 11, 13, 15, 20],
    ),
    (
        "true_tail_coarse_A",
        [1, 2, 3, 4, 5, 7, 9, 13, 20],
    ),
    (
        "true_tail_coarse_B",
        [1, 2, 3, 4, 5, 9, 13, 20],
    ),
    (
        "true_balanced_A",
        [1, 2, 3, 5, 7, 9, 13, 20],
    ),
    (
        "true_balanced_B",
        [1, 3, 5, 7, 11, 15, 20],
    ),
    (
        "true_coarse",
        [1, 3, 5, 9, 13, 20],
    ),
])


# ------------------------------------------------------------
# Cross reco x true, subject to Ntrue <= Nreco
# ------------------------------------------------------------
rows = []
run_index = 0

for _rkey, (reco_source, reco_family, reco_edges) in reco_candidates.items():
    n_reco = len(reco_edges) - 1

    for true_family, true_edges in true_candidates.items():
        n_true = len(true_edges) - 1

        if n_true > n_reco:
            continue

        run_index += 1
        rows.append((
            run_index,
            reco_source,
            reco_family,
            fmt_edges(reco_edges),
            n_reco,
            true_family,
            fmt_edges(true_edges),
            n_true,
        ))

with open(out_path, "w", encoding="utf-8") as f:
    f.write(
        "plan_index\t"
        "reco_source\t"
        "reco_family\t"
        "reco_edges\t"
        "reco_n_bins\t"
        "true_family\t"
        "true_edges\t"
        "true_n_bins\n"
    )

    for row in rows:
        f.write("\t".join(map(str, row)) + "\n")

print(f"[plan] unique reco candidates: {len(reco_candidates)}")
print(f"[plan] true candidates:        {len(true_candidates)}")
print(f"[plan] valid combinations:     {len(rows)}")
PY

PLANNED_RUNS="$(awk 'NR>1 {n++} END {print n+0}' "${PLAN}")"

echo "========== Wide diagnostic jet-pT binning scan =========="
echo "D0 pT range:       ${MIN_PT_D0}-${MAX_PT_D0} GeV/c"
echo "Planned runs:      ${PLANNED_RUNS}"
echo "Scan directory:    ${SCAN_DIR}"
echo "Plan:              ${PLAN}"
echo
echo "Each candidate is applied to ALL centralities in the same Machine.C run."
echo "The scan uses cached response matrices only."
echo

# =========================
# Summary header
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\treco_pt_min\treco_pt_max\treco_n_bins\ttrue_pt_min\ttrue_pt_max\ttrue_n_bins\tmin_width\tmax_width\tstep\ttrend\treco_edges\ttrue_edges\toverride_file\troot_log\treco_reference\tedge_shift\treco_edges_base\treco_source\treco_family\ttrue_family\n" > "${SUMMARY}"
fi

# =========================
# Main loop
# =========================
runs_started=0
ACLIC_SUFFIX="++"

while IFS=$'\t' read -r \
  PLAN_INDEX \
  RECO_SOURCE \
  RECO_FAMILY \
  RECO_EDGES \
  RECO_N_BINS_PLAN \
  TRUE_FAMILY \
  TRUE_EDGES \
  TRUE_N_BINS_PLAN
do
  [[ "${PLAN_INDEX}" == "plan_index" ]] && continue

  printf -v RUN_ID "r%06d" "${PLAN_INDEX}"

  if run_is_complete "${RUN_ID}"; then
    echo "[resume] ${RUN_ID} ${RECO_SOURCE}/${RECO_FAMILY} x ${TRUE_FAMILY} is complete -> skip"
    continue
  fi

  if (( MAX_TOTAL_RUNS > 0 && runs_started >= MAX_TOTAL_RUNS )); then
    echo "Reached MAX_TOTAL_RUNS=${MAX_TOTAL_RUNS}, stopping."
    exit 0
  fi

  remove_run_from_file "${RUN_ID}" "${SUMMARY}"
  remove_run_from_file "${RUN_ID}" "${STABILITY}"

  read -r RECO_PT_MIN RECO_PT_MAX RECO_N_BINS RECO_MIN_WIDTH RECO_MAX_WIDTH \
    < <(edge_stats "${RECO_EDGES}")

  read -r TRUE_PT_MIN TRUE_PT_MAX TRUE_N_BINS TRUE_MIN_WIDTH TRUE_MAX_WIDTH \
    < <(edge_stats "${TRUE_EDGES}")

  if (( TRUE_N_BINS > RECO_N_BINS )); then
    echo "[skip] ${RUN_ID}: N_true=${TRUE_N_BINS} > N_reco=${RECO_N_BINS}"
    continue
  fi

  ((runs_started += 1))

  SAFE_RECO_FAMILY="${RECO_FAMILY//[^A-Za-z0-9_-]/_}"
  SAFE_TRUE_FAMILY="${TRUE_FAMILY//[^A-Za-z0-9_-]/_}"

  OVR_FILE="${OVR_DIR}/override_${RUN_ID}_${SAFE_RECO_FAMILY}_${SAFE_TRUE_FAMILY}.C"
  OUT_DIR="${RUN_DIR}/${RUN_ID}"
  ROOT_LOG="${OUT_DIR}/root.log"

  rm -rf "${OUT_DIR}"
  mkdir -p "${OUT_DIR}"

  # -------------------------
  # Override macro for Machine.C
  # -------------------------
  cat > "${OVR_FILE}" <<EOF_OVR
{
  cout << "[override] wide diagnostic jet-pT binning scan" << endl;
  cout << "[override] reco source: ${RECO_SOURCE}" << endl;
  cout << "[override] reco family: ${RECO_FAMILY}" << endl;
  cout << "[override] true family: ${TRUE_FAMILY}" << endl;
  cout << "[override] jet-pT reco edges: ${RECO_EDGES}" << endl;
  cout << "[override] jet-pT true edges: ${TRUE_EDGES}" << endl;

  // IMPORTANT:
  // Binning scans use the fine cached response. Do not rebuild from TTree.
  UseCachedRM = true;
  FillStandardRM = false;
  FillCacheRM = false;

  useCustomPtRecoBins = true;

  // Apply exactly the same candidate to all centralities.
  // This lets the stability diagnostics compare 0-10%, 10-40%, 40-80%
  // under the same reco/true jet-pT binning.
  for (int ic = 0; ic < nCentralityBins; ++ic) {
    ptRecoBinsVec[ic].clear();
  }

  double reco_edges[] = { ${RECO_EDGES} };
  int nReco = (int)(sizeof(reco_edges) / sizeof(double));

  for (int ic = 0; ic < nCentralityBins; ++ic) {
    for (int i = 0; i < nReco; ++i) {
      ptRecoBinsVec[ic].push_back(reco_edges[i]);
    }
  }

  for (int ic = 0; ic < nCentralityBins; ++ic) {
    ptMcBinsVecCustom[ic].clear();
  }

  double true_edges[] = { ${TRUE_EDGES} };
  int nTrue = (int)(sizeof(true_edges) / sizeof(double));

  for (int ic = 0; ic < nCentralityBins; ++ic) {
    for (int i = 0; i < nTrue; ++i) {
      ptMcBinsVecCustom[ic].push_back(true_edges[i]);
    }
  }
}
EOF_OVR

  TREND="${RECO_SOURCE}_${RECO_FAMILY}__${TRUE_FAMILY}"

  printf "%s\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${RUN_ID}" \
    "${RECO_PT_MIN}" "${RECO_PT_MAX}" "${RECO_N_BINS}" \
    "${TRUE_PT_MIN}" "${TRUE_PT_MAX}" "${TRUE_N_BINS}" \
    "${RECO_MIN_WIDTH}" "${RECO_MAX_WIDTH}" "${EDGE_STEP}" "${TREND}" \
    "${RECO_EDGES}" "${TRUE_EDGES}" "${OVR_FILE}" "${ROOT_LOG}" \
    "${RECO_SOURCE}" "0" "${RECO_EDGES}" \
    "${RECO_SOURCE}" "${RECO_FAMILY}" "${TRUE_FAMILY}" \
    >> "${SUMMARY}"

  ROOT_CMD="${MACHINE_MACRO}${ACLIC_SUFFIX}( ${FONLL_JET}, ${CUT_NEG}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${RUN_ID}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\", \"${SCAN_DIR}\", ${USE_PRIOR_SHAPE_WEIGHTING}, ${SYSTEMATIC_SPLOT} )"

  printf '%s\n' "root -l -b -q '${ROOT_CMD}'" > "${OUT_DIR}/root_command.txt"

  echo "[${RUN_ID}] reco=${RECO_SOURCE}/${RECO_FAMILY}  true=${TRUE_FAMILY}"
  echo "           reco (${RECO_N_BINS} bins): ${RECO_EDGES}"
  echo "           true (${TRUE_N_BINS} bins): ${TRUE_EDGES}"

  if ! root -l -b -q "${ROOT_CMD}" > "${ROOT_LOG}" 2>&1; then
    echo "[error] ROOT failed for ${RUN_ID}. Last 80 lines of the log:"
    tail -n 80 "${ROOT_LOG}" || true
    exit 1
  fi

  ACLIC_SUFFIX="+"

  echo "[done] ${RUN_ID} -> ${OUT_DIR}"

done < "${PLAN}"

echo
echo "Done."
echo "Planned combinations: ${PLAN}"
echo "Summary:              ${SUMMARY}"
echo "Stability:            ${STABILITY}"
echo "Machine runs:         ${RUN_DIR}"
echo "Spectra:              ${SCAN_DIR}/Output/"
