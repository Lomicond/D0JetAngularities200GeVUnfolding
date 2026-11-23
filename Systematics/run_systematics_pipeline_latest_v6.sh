#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Range-aware systematics pipeline
#
# Keeps different D0-pT selections completely separated:
#   - reads only scan..._<MIN>_<MAX> directories
#   - writes range-specific long/wide tables
#   - uses a range-specific combiner config
#   - writes a range-specific final ROOT file
#   - does NOT run finalPlot_from_systematics_root_*.py
#
# Usage from the project root:
#   ./Systematics/run_systematics_pipeline_latest_v6.sh 1 5 10
#   ./Systematics/run_systematics_pipeline_latest_v6.sh 2 5 10
#
# Mode:
#   1 = GUI config editor + non-GUI combiner
#   2 = non-GUI combiner only
# ============================================================

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 MODE MIN_PT_D0 MAX_PT_D0"
  echo "Example: $0 1 5 10"
  echo
  echo "MODE:"
  echo "  1 = GUI config editor + non-GUI combiner"
  echo "  2 = non-GUI combiner only"
  exit 1
fi

MODE="$1"
MIN_PT_D0="$2"
MAX_PT_D0="$3"

if ! [[ "${MIN_PT_D0}" =~ ^[0-9]+$ && "${MAX_PT_D0}" =~ ^[0-9]+$ ]]; then
  echo "[error] MIN_PT_D0 and MAX_PT_D0 must be integers." >&2
  exit 1
fi

if (( MIN_PT_D0 >= MAX_PT_D0 )); then
  echo "[error] MIN_PT_D0 must be smaller than MAX_PT_D0." >&2
  exit 1
fi

D0_TAG="${MIN_PT_D0}_${MAX_PT_D0}"

SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
SCRIPT_BASE="$(basename "$SCRIPT_DIR")"

if [[ "$SCRIPT_BASE" == "Systematics" ]]; then
  PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
else
  PROJECT_DIR="$(pwd)"
  if [[ -d "$PROJECT_DIR/Systematics" ]]; then
    SCRIPT_DIR="$PROJECT_DIR/Systematics"
  fi
fi

cd "$PROJECT_DIR"

need_file() {
  local path="$1"
  if [[ ! -e "$path" ]]; then
    echo "[error] Missing file: $path" >&2
    exit 1
  fi
}

need_dir() {
  local path="$1"
  if [[ ! -d "$path" ]]; then
    echo "[error] Missing range-specific scan directory:" >&2
    echo "        $path" >&2
    echo "[error] Refusing to fall back to an unsuffixed/older scan." >&2
    exit 1
  fi
}

need_nonempty_file() {
  local path="$1"
  if [[ ! -s "$path" ]]; then
    echo "[error] Expected non-empty file was not created: $path" >&2
    exit 1
  fi
}

latest_versioned_file() {
  local pattern="$1"
  local best=""
  local bestv=-1
  local f base v

  shopt -s nullglob
  for f in $pattern; do
    base="$(basename "$f")"
    if [[ "$base" =~ _v([0-9]+)\. ]]; then
      v="${BASH_REMATCH[1]}"
      if (( v > bestv )); then
        bestv="$v"
        best="$f"
      fi
    fi
  done
  shopt -u nullglob

  if [[ -z "$best" ]]; then
    echo ""
    return 1
  fi
  echo "$best"
}

version_of() {
  local base
  base="$(basename "$1")"
  if [[ "$base" =~ _v([0-9]+)\. ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo ""
  fi
}

choose_config() {
  local real_pattern="$1"
  local template_pattern="$2"
  local cfg=""

  cfg="$(latest_versioned_file "$real_pattern" 2>/dev/null || true)"
  if [[ -z "$cfg" ]]; then
    cfg="$(latest_versioned_file "$template_pattern" 2>/dev/null || true)"
  fi
  echo "$cfg"
}

run_cmd() {
  echo "+ $*"
  "$@"
}

config_root_file() {
  local cfg="$1"
  python3 - "$cfg" "$PROJECT_DIR" <<'PY'
import json
import os
import sys
from pathlib import Path

cfg_path = Path(sys.argv[1])
project_dir_from_wrapper = Path(sys.argv[2]).resolve()

with cfg_path.open("r", encoding="utf-8") as f:
    cfg = json.load(f)

project_dir = Path(os.path.expanduser(os.path.expandvars(str(cfg.get("project_dir", ".")))))
if not project_dir.is_absolute():
    project_dir = (project_dir_from_wrapper / project_dir).resolve()

root_path = Path(os.path.expanduser(os.path.expandvars(
    str(cfg.get("output_root", "Systematics/final_systematics_results.root"))
)))
if not root_path.is_absolute():
    root_path = (project_dir / root_path).resolve()

print(root_path)
PY
}

# ------------------------------------------------------------------
# Select the newest analysis scripts and the normal/base JSON configs.
# The range-specific JSONs are generated below.
# ------------------------------------------------------------------
COMP_BUILDER="$(latest_versioned_file "$SCRIPT_DIR/systematics_components_builder_v*.py")"
BASE_COMP_CONFIG="$(choose_config \
  "$SCRIPT_DIR/systematics_components_config_v*.json" \
  "$SCRIPT_DIR/systematics_components_config_template_v*.json")"
WIDE_BUILDER="$(latest_versioned_file "$SCRIPT_DIR/systematics_components_wide_builder_v*.py")"
NON_GUI_COMBINER="$(latest_versioned_file "$SCRIPT_DIR/systematics_final_combiner_v*.py")"
GUI_COMBINER="$(latest_versioned_file "$SCRIPT_DIR/systematics_final_combiner_gui_v*.py" 2>/dev/null || true)"

need_file "$COMP_BUILDER"
need_file "$BASE_COMP_CONFIG"
need_file "$WIDE_BUILDER"
need_file "$NON_GUI_COMBINER"

case "$MODE" in
  1|gui|GUI)
    MODE_NAME="gui"
    need_file "$GUI_COMBINER"
    BASE_COMBINER_CONFIG="$(choose_config \
      "$SCRIPT_DIR/systematics_final_combiner_gui_config_v*.json" \
      "$SCRIPT_DIR/systematics_final_combiner_gui_config_template_v*.json")"
    ;;
  2|nogui|non-gui|NOGUI)
    MODE_NAME="nogui"
    BASE_COMBINER_CONFIG="$(choose_config \
      "$SCRIPT_DIR/systematics_final_combiner_config_v*.json" \
      "$SCRIPT_DIR/systematics_final_combiner_config_template_v*.json")"
    ;;
  *)
    echo "[error] Unknown mode: $MODE" >&2
    exit 1
    ;;
esac

need_file "$BASE_COMBINER_CONFIG"

COMP_V="$(version_of "$COMP_BUILDER")"
WIDE_V="$(version_of "$WIDE_BUILDER")"

COMP_TSV="Systematics/systematics_components_v${COMP_V}_${D0_TAG}.tsv"
COMP_CSV="Systematics/systematics_components_v${COMP_V}_${D0_TAG}.csv"
WIDE_TSV="Systematics/systematics_components_wide_v${WIDE_V}_${D0_TAG}.tsv"
WIDE_CSV="Systematics/systematics_components_wide_v${WIDE_V}_${D0_TAG}.csv"

RANGE_CONFIG_DIR="$SCRIPT_DIR/range_configs/${D0_TAG}"
mkdir -p "$RANGE_CONFIG_DIR"

BASE_COMP_NAME="$(basename "$BASE_COMP_CONFIG" .json)"
BASE_COMB_NAME="$(basename "$BASE_COMBINER_CONFIG" .json)"

RANGE_COMP_CONFIG="$RANGE_CONFIG_DIR/${BASE_COMP_NAME}_${D0_TAG}.json"
RANGE_COMBINER_CONFIG="$RANGE_CONFIG_DIR/${BASE_COMB_NAME}_${D0_TAG}.json"

# ------------------------------------------------------------------
# Hard preflight: all inputs must have the requested D0-pT suffix.
# This deliberately prevents accidental mixing with the old 1-10 scans.
# ------------------------------------------------------------------
RANGE_SCAN_DIRS=(
  "Systematics/Nominal/scanNominal_${D0_TAG}"
  "Systematics/JetsReco/scanJetsReco_${D0_TAG}"
  "Systematics/sWeight/scanSWeight_${D0_TAG}"
  "Systematics/PriorShape/scanPriorShape_${D0_TAG}"
  "Systematics/D0Meson/scanPaperSys_${D0_TAG}"
  "Systematics/Binning/scanJetPt_${D0_TAG}"
  "Systematics/Binning/scanZ_${D0_TAG}"
  "Systematics/Binning/scanL11_${D0_TAG}"
  "Systematics/Binning/scanL11p5_${D0_TAG}"
  "Systematics/Binning/scanL12_${D0_TAG}"
  "Systematics/Binning/scanL13_${D0_TAG}"
  "Systematics/Binning/scanL10p5_${D0_TAG}"
  "Systematics/Binning/scanPTD_${D0_TAG}"
)

for d in "${RANGE_SCAN_DIRS[@]}"; do
  need_dir "$d"
done

# ------------------------------------------------------------------
# Build the range-aware COMPONENTS config.
#
# Important detail:
# The old config has a jet-pT endpoint filter fixed to min=1.
# Our new scanJetPt_<range> has already been pre-filtered at generation
# time, including the selected old binning style and the 5-GeV lower
# region.  Therefore the obsolete downstream min=1 filter is disabled
# for binning_first_variable; otherwise a 5-10 scan would be rejected.
#
# The second-variable endpoint filters are left unchanged.
# ------------------------------------------------------------------
python3 - \
  "$BASE_COMP_CONFIG" \
  "$RANGE_COMP_CONFIG" \
  "$D0_TAG" \
  "$COMP_TSV" \
  "$COMP_CSV" <<'PY'
import json
import sys
from pathlib import Path

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
tag = sys.argv[3]
comp_tsv = sys.argv[4]
comp_csv = sys.argv[5]

with src.open("r", encoding="utf-8") as f:
    cfg = json.load(f)

cfg["output_tsv"] = comp_tsv
cfg["output_csv"] = comp_csv

# Nominal
nom = cfg.get("nominal", {})
if isinstance(nom, dict) and "file_glob" in nom:
    nom["file_glob"] = f"Systematics/Nominal/scanNominal_{tag}/Output/OutputSpectra*.root"

# Ordinary systematic groups
group_scan_map = {
    "JetsReco":   f"Systematics/JetsReco/scanJetsReco_{tag}",
    "sWeight":    f"Systematics/sWeight/scanSWeight_{tag}",
    "PriorShape": f"Systematics/PriorShape/scanPriorShape_{tag}",
    "D0Meson":    f"Systematics/D0Meson/scanPaperSys_{tag}",
}
for g in cfg.get("groups", []):
    source_group = str(g.get("source_group", ""))
    if source_group in group_scan_map:
        g["scan_dir"] = group_scan_map[source_group]

# Binning sources
binning = cfg.get("binning", {})
for source in binning.get("sources", []):
    source_name = source.get("source_name", "")

    if source_name == "binning_first_variable":
        source["scan_dir"] = f"Systematics/Binning/scanJetPt_{tag}"

        # The new jet-pT scan is already filtered when generated.
        if isinstance(source.get("true_range_filter"), dict):
            source["true_range_filter"]["enabled"] = False
            source["true_range_filter"]["mode"] = "ignore"
        if isinstance(source.get("reco_range_filter"), dict):
            source["reco_range_filter"]["enabled"] = False
            source["reco_range_filter"]["mode"] = "ignore"

        old_notes = str(source.get("notes", ""))
        source["notes"] = (
            old_notes
            + f" Range-aware wrapper: scanJetPt_{tag} is pre-filtered; "
              "legacy endpoint filters disabled."
        ).strip()

    elif source_name == "binning_second_variable":
        source["scan_dir_map"] = {
            "z":       f"Systematics/Binning/scanZ_{tag}",
            "l11":     f"Systematics/Binning/scanL11_{tag}",
            "l15":     f"Systematics/Binning/scanL11p5_{tag}",
            "l21":     f"Systematics/Binning/scanL12_{tag}",
            "l31":     f"Systematics/Binning/scanL13_{tag}",
            "l051":    f"Systematics/Binning/scanL10p5_{tag}",
            "pTD":     f"Systematics/Binning/scanPTD_{tag}",
            "default": f"Systematics/Binning/scanZ_{tag}",
        }

dst.parent.mkdir(parents=True, exist_ok=True)
with dst.open("w", encoding="utf-8") as f:
    json.dump(cfg, f, indent=2)
    f.write("\n")

print(f"[info] Wrote range-specific components config: {dst}")
PY

# ------------------------------------------------------------------
# Build/update the range-aware FINAL-COMBINER config.
#
# First run for a range:
#   clone the normal/base config, then suffix all outputs.
#
# Later runs:
#   reuse the existing range config so GUI choices are preserved,
#   but reassert the range-specific input/output paths.
# ------------------------------------------------------------------
python3 - \
  "$BASE_COMBINER_CONFIG" \
  "$RANGE_COMBINER_CONFIG" \
  "$D0_TAG" \
  "$COMP_TSV" <<'PY'
import json
import sys
from pathlib import Path

base_path = Path(sys.argv[1])
range_path = Path(sys.argv[2])
tag = sys.argv[3]
comp_tsv = sys.argv[4]

source_path = range_path if range_path.exists() else base_path

with source_path.open("r", encoding="utf-8") as f:
    cfg = json.load(f)

def suffix_path(value, tag):
    p = Path(str(value))
    if not p.name:
        return str(p)
    if p.suffix:
        return str(p.with_name(f"{p.stem}_{tag}{p.suffix}"))
    return str(p.with_name(f"{p.name}_{tag}"))

cfg["input"] = comp_tsv

for key in ("output_tsv", "output_csv", "details_tsv", "output_root"):
    if key in cfg:
        current = str(cfg[key])

        # Avoid repeatedly appending the same tag on subsequent runs.
        stemmed_tag = f"_{tag}"
        p = Path(current)
        if p.suffix:
            if not p.stem.endswith(stemmed_tag):
                cfg[key] = str(p.with_name(f"{p.stem}_{tag}{p.suffix}"))
        else:
            if not p.name.endswith(stemmed_tag):
                cfg[key] = str(p.with_name(f"{p.name}_{tag}"))

plots = cfg.get("plots")
if isinstance(plots, dict):
    # Give each D0 range a fully independent plot directory.
    plots["output_dir"] = f"Systematics/FinalPlots_{tag}"

range_path.parent.mkdir(parents=True, exist_ok=True)
with range_path.open("w", encoding="utf-8") as f:
    json.dump(cfg, f, indent=2)
    f.write("\n")

print(f"[info] Using range-specific combiner config: {range_path}")
PY

cat <<EOM
[info] Working directory : $PROJECT_DIR
[info] D0 pT range       : ${MIN_PT_D0}-${MAX_PT_D0} GeV/c
[info] Range tag         : ${D0_TAG}
[info] Combiner mode     : ${MODE_NAME}

[info] Selected files:
  components builder : $(basename "$COMP_BUILDER")
  base comp. config  : $(basename "$BASE_COMP_CONFIG")
  range comp. config : $RANGE_COMP_CONFIG
  wide builder       : $(basename "$WIDE_BUILDER")
  combiner script    : $(basename "$NON_GUI_COMBINER")
  base comb. config  : $(basename "$BASE_COMBINER_CONFIG")
  range comb. config : $RANGE_COMBINER_CONFIG

[info] Range-specific outputs:
  long TSV           : $COMP_TSV
  long CSV           : $COMP_CSV
  wide TSV           : $WIDE_TSV
  wide CSV           : $WIDE_CSV
EOM

echo
printf '========== Build long systematics table ==========\n'
run_cmd python3 "$COMP_BUILDER" --config "$RANGE_COMP_CONFIG"

need_nonempty_file "$COMP_TSV"
rows=$(( $(wc -l < "$COMP_TSV") - 1 ))
echo "[info] Long table rows: $rows"
if (( rows <= 0 )); then
  echo "[error] Long table has zero data rows: $COMP_TSV" >&2
  exit 1
fi

echo
printf '========== Build wide systematics table ==========\n'
run_cmd python3 "$WIDE_BUILDER" \
  --input "$COMP_TSV" \
  --output-tsv "$WIDE_TSV" \
  --output-csv "$WIDE_CSV"

need_nonempty_file "$WIDE_TSV"

if [[ "$MODE_NAME" == "gui" ]]; then
  echo
  printf '========== Edit range-specific final-combiner config in GUI ==========\n'
  echo "[info] This GUI edits ONLY:"
  echo "       $RANGE_COMBINER_CONFIG"
  echo "[info] Your normal/old config is not modified."
  echo "[info] After editing, click Save config and close the GUI."

  run_cmd python3 "$GUI_COMBINER" \
    --config "$RANGE_COMBINER_CONFIG" \
    --input "$COMP_TSV" \
    --combiner "$NON_GUI_COMBINER"

  need_file "$RANGE_COMBINER_CONFIG"
fi

ROOT_FILE="$(config_root_file "$RANGE_COMBINER_CONFIG")"
echo
echo "[info] Final ROOT file requested by range config:"
echo "       $ROOT_FILE"

# Because the ROOT filename itself contains the D0 range, replacing it is safe:
# no unsuffixed/old result can be overwritten.
if [[ -e "$ROOT_FILE" ]]; then
  stamp="$(date +%Y%m%d_%H%M%S)"
  backup="${ROOT_FILE}.bak_${stamp}"
  echo "[info] Existing SAME-RANGE ROOT file will be moved aside:"
  echo "       $backup"
  mv "$ROOT_FILE" "$backup"
fi

echo
printf '========== Combine final systematics and write fresh ROOT ==========\n'
run_cmd python3 "$NON_GUI_COMBINER" \
  --config "$RANGE_COMBINER_CONFIG" \
  --force-root \
  --output-root "$ROOT_FILE"

need_nonempty_file "$ROOT_FILE"

echo
printf '========== Done ==========\n'
echo "[info] Final ROOT: $ROOT_FILE"
echo "[info] Final plotting was intentionally NOT run."
echo "[info] Different D0-pT ranges remain in separate files/directories."
