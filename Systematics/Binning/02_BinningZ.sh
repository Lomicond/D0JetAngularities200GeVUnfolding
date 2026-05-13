#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Binning systematic scan: z
#
# Run from the main project directory, i.e. the directory containing:
#   config.h, config_hist.h, Unfolding/Machine.C, Data/...
#
# Recommended run command:
#   bash Systematics/Binning/02_BinningZ.sh
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
# z scan
# -------------------------
Z_TRUE_START_LIST=(0 0.2 0.4)
Z_RECO_START_LIST=(0 0.2 0.4)

# True-level z edges are obtained by taking a suffix of this baseline list.
Z_TRUE_BASE_EDGES=(0 0.2 0.4 0.6 0.7 0.8 0.9 1.01)

# Reco-level z edges are generated from the right to the left.
# The last bin is fixed to Z_LAST_BIN_LEFT -> Z_RIGHT_EDGE.
Z_RECO_ALLOWED_WIDTHS=(0.1 0.2 0.3)
Z_RECO_MIN_WIDTH=0.1
Z_RIGHT_EDGE=1.01
Z_LAST_BIN_LEFT=0.9

# -------------------------
# Optional coherent z-edge nudging
# -------------------------
# The THnSparse z granularity is fine enough for 0.05 shifts.
# To avoid a combinatorial explosion, do not move each edge independently.
# Instead, use a small set of coherent variants:
#   nominal     : no edge shift
#   reco_minus  : move internal reco-z edges by -Z_EDGE_NUDGE_STEP
#   reco_plus   : move internal reco-z edges by +Z_EDGE_NUDGE_STEP
#
# Only reco-level edges are nudged. True-level edges are not shifted.
# Endpoints are kept fixed. The 0.9 edge is also kept fixed, so the last
# 0.9 -> 1.01 bin stays unchanged.
Z_EDGE_NUDGE_STEP=0.01
Z_EDGE_NUDGE_MIN_WIDTH=0.01
Z_NUDGE_MODES=(nominal reco_minus reco_plus)

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
MIN_PT_D0=1
MAX_PT_D0=10
USE_PRIOR_SHAPE_WEIGHTING=0
SYSTEMATIC_SPLOT=0

# Safety limit for the whole scan. 0 = no limit.
MAX_TOTAL_RUNS=0

# -------------------------
# Output location
# -------------------------
SCAN_DIR="${SCRIPT_DIR}/scanZ"
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

join_by_comma_plain() {
  local out=""
  local x
  for x in "$@"; do
    out+="${out:+,}${x}"
  done
  printf '%s' "${out}"
}

# =========================
# Summary header
# =========================
if [[ ! -f "${SUMMARY}" ]]; then
  printf "run_id\tz_true_start\tz_true_n_bins\tz_reco_start\tz_reco_n_bins\tz_nudge_mode\tz_reco_edges\tz_true_edges\toverride_file\troot_log\n" > "${SUMMARY}"
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
# Use forced ACLiC rebuild only for the first real Machine() call.
# After that, use the already compiled library.
ACLIC_SUFFIX="++"

# Convert bash arrays to compact strings for the embedded Python generator.
Z_TRUE_STARTS_CSV="$(join_by_comma_plain "${Z_TRUE_START_LIST[@]}")"
Z_RECO_STARTS_CSV="$(join_by_comma_plain "${Z_RECO_START_LIST[@]}")"
Z_TRUE_BASE_CSV="$(join_by_comma_plain "${Z_TRUE_BASE_EDGES[@]}")"
Z_RECO_WIDTHS_CSV="$(join_by_comma_plain "${Z_RECO_ALLOWED_WIDTHS[@]}")"
Z_NUDGE_MODES_CSV="$(join_by_comma_plain "${Z_NUDGE_MODES[@]}")"

# =========================
# Main loop
# =========================
while IFS=$'\t' read -r Z_TRUE_START Z_TRUE_N_BINS Z_RECO_START Z_RECO_N_BINS Z_NUDGE_MODE Z_RECO_EDGES Z_TRUE_EDGES; do
  [[ -n "${Z_RECO_EDGES}" ]] || continue

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
  # Only z binning is changed. All pT and angularity binnings stay exactly
  # as defined in config.h / config_hist.h.
  # -------------------------
  cat > "${OVR_FILE}" <<EOF_OVR
{
  cout << "[override] z reco edges: ${Z_RECO_EDGES}" << endl;
  cout << "[override] z true edges: ${Z_TRUE_EDGES}" << endl;
  cout << "[override] z nudge mode: ${Z_NUDGE_MODE}" << endl;

  for (int ic = 0; ic < nCentralityBins; ++ic) zRecoBinsVec[ic].clear();
  double z_reco_edges[] = { ${Z_RECO_EDGES} };
  int nZReco = (int)(sizeof(z_reco_edges)/sizeof(double));
  for (int ic = 0; ic < nCentralityBins; ++ic) {
    for (int i = 0; i < nZReco; ++i) zRecoBinsVec[ic].push_back(z_reco_edges[i]);
  }

  for (int ic = 0; ic < nCentralityBins; ++ic) zMcBinsVecCustom[ic].clear();
  double z_true_edges[] = { ${Z_TRUE_EDGES} };
  int nZTrue = (int)(sizeof(z_true_edges)/sizeof(double));
  for (int ic = 0; ic < nCentralityBins; ++ic) {
    for (int i = 0; i < nZTrue; ++i) zMcBinsVecCustom[ic].push_back(z_true_edges[i]);
  }
}
EOF_OVR

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${RUN_ID}" \
    "${Z_TRUE_START}" "${Z_TRUE_N_BINS}" \
    "${Z_RECO_START}" "${Z_RECO_N_BINS}" "${Z_NUDGE_MODE}" \
    "${Z_RECO_EDGES}" "${Z_TRUE_EDGES}" "${OVR_FILE}" "${ROOT_LOG}" \
    >> "${SUMMARY}"

  echo "[${RUN_ID}] z true start=${Z_TRUE_START}, reco start=${Z_RECO_START}, Nreco=${Z_RECO_N_BINS}, Ntrue=${Z_TRUE_N_BINS}, mode=${Z_NUDGE_MODE}"

  # -------------------------
  # Run the ROOT macro
  # OutputFile = RUN_ID because Machine.C uses runId = BaseName(OutputFile)
  # ScanDir    = SCAN_DIR so stability.tsv and Output/OutputSpectra*.root go there
  # -------------------------
  if ! root -l -b -q \
    "${MACHINE_MACRO}${ACLIC_SUFFIX}( ${FONLL_JET}, ${CUT_NEG}, ${MIN_JET_PT_RECO_CUT}, ${SAVED_ITER}, \"${INPUT_FILE}\", \"${RUN_ID}\", ${MIN_PT_D0}, ${MAX_PT_D0}, \"${OVR_FILE}\", \"${SCAN_DIR}\", ${USE_PRIOR_SHAPE_WEIGHTING}, ${SYSTEMATIC_SPLOT} )" \
    > "${ROOT_LOG}" 2>&1; then
    echo "[error] ROOT failed for ${RUN_ID}. Last 60 lines of the log:"
    tail -n 60 "${ROOT_LOG}" || true
    exit 1
  fi

  # Only the first run should force recompilation.
  ACLIC_SUFFIX="+"

done < <(
  python3 - \
    "${Z_TRUE_STARTS_CSV}" \
    "${Z_RECO_STARTS_CSV}" \
    "${Z_TRUE_BASE_CSV}" \
    "${Z_RECO_WIDTHS_CSV}" \
    "${Z_RECO_MIN_WIDTH}" \
    "${Z_RIGHT_EDGE}" \
    "${Z_LAST_BIN_LEFT}" \
    "${Z_EDGE_NUDGE_STEP}" \
    "${Z_EDGE_NUDGE_MIN_WIDTH}" \
    "${Z_NUDGE_MODES_CSV}" <<'PY'
import sys

true_starts = [float(x) for x in sys.argv[1].split(',') if x]
reco_starts = [float(x) for x in sys.argv[2].split(',') if x]
true_base = [float(x) for x in sys.argv[3].split(',') if x]
reco_allowed_widths = [float(x) for x in sys.argv[4].split(',') if x]
reco_min_width = float(sys.argv[5])
right_edge = float(sys.argv[6])
last_bin_left = float(sys.argv[7])
nudge_step = float(sys.argv[8])
nudge_min_width = float(sys.argv[9])
nudge_modes = [x.strip() for x in sys.argv[10].split(',') if x.strip()]

EPS = 1e-9


def fmt_edges(edges):
    out = []
    for x in edges:
        s = f"{x:.10f}".rstrip('0').rstrip('.')
        out.append(s if s else '0')
    return ", ".join(out)


def build_true_edges(start):
    edges = [x for x in true_base if x >= start - EPS]
    if len(edges) < 2:
        return None
    if abs(edges[0] - start) > EPS:
        return None
    return tuple(round(x, 10) for x in edges)


def generate_reco_edges(start):
    if start > last_bin_left - EPS:
        return []

    results = set()

    # acc_edges_desc stores edges from right to left, including right_edge
    # and the current left edge. The last bin is fixed to last_bin_left -> right_edge.
    def rec(prev_width_effective, acc_edges_desc):
        current_left = acc_edges_desc[-1]
        rem = current_left - start

        # Close with one transition bin on the far left.
        if rem >= prev_width_effective - EPS and rem >= reco_min_width - EPS:
            edges_desc = acc_edges_desc + [start]
            edges = list(reversed(edges_desc))
            if abs(edges[-1] - right_edge) < 1e-8 and abs(edges[0] - start) < 1e-8:
                results.add(tuple(round(x, 10) for x in edges))

        # Add another regular bin on the left.
        for w in reco_allowed_widths:
            if w + EPS < prev_width_effective:
                continue
            new_left = current_left - w
            if new_left <= start + EPS:
                continue
            rec(w, acc_edges_desc + [new_left])

    rec(reco_min_width, [right_edge, last_bin_left])
    return sorted(results)


def shifted_edges(edges, shift):
    """Move internal edges coherently by shift.

    Keep the first endpoint, the final right endpoint, and the fixed 0.9 edge.
    This preserves the special last bin 0.9 -> 1.01.
    """
    shifted = list(edges)
    for i in range(1, len(edges) - 1):
        if abs(edges[i] - last_bin_left) < EPS:
            continue
        shifted[i] = round(edges[i] + shift, 10)

    if any(x < -EPS or x > right_edge + EPS for x in shifted):
        return None

    for a, b in zip(shifted, shifted[1:]):
        if b - a < nudge_min_width - EPS:
            return None

    return tuple(round(x, 10) for x in shifted)


def apply_mode(reco_edges, true_edges, mode):
    # Nudge only the reco-level z edges.
    # True-level z edges are never shifted here; they are changed only by Z_TRUE_START_LIST.
    if mode == "nominal":
        return reco_edges, true_edges
    if mode == "reco_minus":
        return shifted_edges(reco_edges, -nudge_step), true_edges
    if mode == "reco_plus":
        return shifted_edges(reco_edges, +nudge_step), true_edges

    print(f"[warning] Unknown Z_NUDGE_MODE='{mode}', skipping", file=sys.stderr)
    return None, None


seen = set()

for true_start in true_starts:
    true_edges_base = build_true_edges(true_start)
    if true_edges_base is None:
        continue

    for reco_start in reco_starts:
        # Reco range must be equal or wider than true range.
        if reco_start > true_start + EPS:
            continue

        for reco_edges_base in generate_reco_edges(reco_start):
            if len(reco_edges_base) - 1 < len(true_edges_base) - 1:
                print(
                    f"[skip] z true start={true_start:g}, reco start={reco_start:g} -> "
                    f"Nreco={len(reco_edges_base)-1} < Ntrue={len(true_edges_base)-1}",
                    file=sys.stderr,
                )
                continue

            for mode in nudge_modes:
                reco_edges, true_edges = apply_mode(reco_edges_base, true_edges_base, mode)
                if reco_edges is None or true_edges is None:
                    continue
                if len(reco_edges) - 1 < len(true_edges) - 1:
                    continue

                key = (reco_edges, true_edges)
                if key in seen:
                    continue
                seen.add(key)

                print(
                    f"{true_start:g}\t{len(true_edges)-1}\t"
                    f"{reco_start:g}\t{len(reco_edges)-1}\t"
                    f"{mode}\t{fmt_edges(reco_edges)}\t{fmt_edges(true_edges)}"
                )
PY
)

echo "Done."
echo "Summary:   ${SUMMARY}"
echo "Stability: ${SCAN_DIR}/stability.tsv"
echo "Spectra:   ${SCAN_DIR}/Output/"
