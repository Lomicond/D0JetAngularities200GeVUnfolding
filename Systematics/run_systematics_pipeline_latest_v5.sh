#!/usr/bin/env bash
set -euo pipefail

# Run the full systematics pipeline using the newest versioned scripts.
# The final ROOT file is always produced by systematics_final_combiner_v*.py
# from the currently selected config.  No pre-existing write_final_*.C macro is
# executed by this wrapper, so stale generated ROOT macros cannot overwrite the
# current GUI/non-GUI choices.
#
# Usage from project root:
#   ./Systematics/run_systematics_pipeline_latest_v5.sh 1   # GUI config editor, then run combiner
#   ./Systematics/run_systematics_pipeline_latest_v5.sh 2   # non-GUI combiner
#
# Mode:
#   1 = GUI config editor + non-GUI combiner
#   2 = non-GUI combiner only

MODE="${1:-}"

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
  local base="$(basename "$1")"
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

choose_final_plotter() {
  local plotter=""
  plotter="$(latest_versioned_file "$SCRIPT_DIR/finalPlot_from_systematics_root_v*.py" 2>/dev/null || true)"
  if [[ -z "$plotter" && -e "$SCRIPT_DIR/finalPlot_from_systematics_root.py" ]]; then
    plotter="$SCRIPT_DIR/finalPlot_from_systematics_root.py"
  fi
  echo "$plotter"
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
if not cfg_path.exists():
    raise SystemExit(f"missing config: {cfg_path}")

with cfg_path.open("r", encoding="utf-8") as f:
    cfg = json.load(f)

project_dir = Path(os.path.expanduser(os.path.expandvars(str(cfg.get("project_dir", ".")))))
if not project_dir.is_absolute():
    project_dir = (project_dir_from_wrapper / project_dir).resolve()

root_path = Path(os.path.expanduser(os.path.expandvars(str(cfg.get("output_root", "Systematics/final_systematics_results_v6.root")))))
if not root_path.is_absolute():
    root_path = (project_dir / root_path).resolve()

print(root_path)
PY
}

backup_existing_file() {
  local path="$1"
  if [[ -e "$path" ]]; then
    local stamp backup
    stamp="$(date +%Y%m%d_%H%M%S)"
    backup="${path}.bak_${stamp}"
    echo "[info] Existing ROOT file will be moved aside to prevent stale plotting:"
    echo "       $backup"
    mv "$path" "$backup"
  fi
}

# Files live in Systematics/ even though we run from project root.
COMP_BUILDER="$(latest_versioned_file "$SCRIPT_DIR/systematics_components_builder_v*.py")"
COMP_CONFIG="$(choose_config "$SCRIPT_DIR/systematics_components_config_v*.json" "$SCRIPT_DIR/systematics_components_config_template_v*.json")"
WIDE_BUILDER="$(latest_versioned_file "$SCRIPT_DIR/systematics_components_wide_builder_v*.py")"
NON_GUI_COMBINER="$(latest_versioned_file "$SCRIPT_DIR/systematics_final_combiner_v*.py")"
GUI_COMBINER="$(latest_versioned_file "$SCRIPT_DIR/systematics_final_combiner_gui_v*.py" 2>/dev/null || true)"
FINAL_PLOTTER="$(choose_final_plotter)"

need_file "$COMP_BUILDER"
need_file "$COMP_CONFIG"
need_file "$WIDE_BUILDER"
need_file "$NON_GUI_COMBINER"
need_file "$FINAL_PLOTTER"

if [[ -z "$MODE" ]]; then
  echo "Choose final combiner mode:"
  echo "  1) GUI config editor, then non-GUI combiner"
  echo "  2) non-GUI combiner"
  read -r -p "Selection [1/2]: " MODE
fi

case "$MODE" in
  1|gui|GUI)
    MODE_NAME="gui"
    need_file "$GUI_COMBINER"
    COMBINER_CONFIG="$(choose_config "$SCRIPT_DIR/systematics_final_combiner_gui_config_v*.json" "$SCRIPT_DIR/systematics_final_combiner_gui_config_template_v*.json")"
    if [[ -z "$COMBINER_CONFIG" ]]; then
      COMBINER_CONFIG="$SCRIPT_DIR/systematics_final_combiner_gui_config_v5.json"
      echo "[info] No GUI combiner config found yet. The GUI will save a new one here:"
      echo "       $COMBINER_CONFIG"
    fi
    ;;
  2|nogui|non-gui|NOGUI)
    MODE_NAME="nogui"
    COMBINER_CONFIG="$(choose_config "$SCRIPT_DIR/systematics_final_combiner_config_v*.json" "$SCRIPT_DIR/systematics_final_combiner_config_template_v*.json")"
    need_file "$COMBINER_CONFIG"
    ;;
  *)
    echo "[error] Unknown mode: $MODE" >&2
    exit 1
    ;;
esac

COMP_V="$(version_of "$COMP_BUILDER")"
WIDE_V="$(version_of "$WIDE_BUILDER")"
COMB_V="$(version_of "$NON_GUI_COMBINER")"
PLOT_V="$(version_of "$FINAL_PLOTTER")"

COMP_TSV="Systematics/systematics_components_v${COMP_V}.tsv"
COMP_CSV="Systematics/systematics_components_v${COMP_V}.csv"
WIDE_TSV="Systematics/systematics_components_wide_v${WIDE_V}.tsv"
WIDE_CSV="Systematics/systematics_components_wide_v${WIDE_V}.csv"

cat <<EOM
[info] Working directory : $PROJECT_DIR
[info] Script directory  : $SCRIPT_DIR
[info] Selected files:
  components builder : $(basename "$COMP_BUILDER")
  components config  : $(basename "$COMP_CONFIG")
  wide builder       : $(basename "$WIDE_BUILDER")
  combiner mode      : $MODE_NAME
  GUI editor         : ${GUI_COMBINER:+$(basename "$GUI_COMBINER")}
  combiner script    : $(basename "$NON_GUI_COMBINER")
  combiner config    : $(basename "$COMBINER_CONFIG")
  final plotter      : $(basename "$FINAL_PLOTTER")
EOM

echo
printf '========== Build long systematics table ==========%s' ""; echo
run_cmd python3 "$COMP_BUILDER" --config "$COMP_CONFIG"

need_nonempty_file "$COMP_TSV"
rows=$(( $(wc -l < "$COMP_TSV") - 1 ))
echo "[info] Long table rows: $rows"
if (( rows <= 0 )); then
  echo "[error] Long table has zero data rows: $COMP_TSV" >&2
  exit 1
fi

echo
printf '========== Build wide systematics table ==========%s' ""; echo
run_cmd python3 "$WIDE_BUILDER" --input "$COMP_TSV" --output-tsv "$WIDE_TSV" --output-csv "$WIDE_CSV"
need_nonempty_file "$WIDE_TSV"

if [[ "$MODE_NAME" == "gui" ]]; then
  echo
  printf '========== Edit final-combiner config in GUI ==========%s' ""; echo
  echo "[info] After editing, click Save config and close the GUI."
  echo "[info] The wrapper will then run the non-GUI combiner using that saved config."
  run_cmd python3 "$GUI_COMBINER" \
    --config "$COMBINER_CONFIG" \
    --input "$COMP_TSV" \
    --combiner "$NON_GUI_COMBINER"
  need_file "$COMBINER_CONFIG"
fi

ROOT_FILE="$(config_root_file "$COMBINER_CONFIG")"
echo "[info] Final ROOT file requested by config: $ROOT_FILE"
backup_existing_file "$ROOT_FILE"

echo
printf '========== Combine final systematics and write fresh ROOT ==========%s' ""; echo
run_cmd python3 "$NON_GUI_COMBINER" --config "$COMBINER_CONFIG" --force-root --output-root "$ROOT_FILE"
need_nonempty_file "$ROOT_FILE"

echo
printf '========== Draw final plots from fresh ROOT ==========%s' ""; echo
run_cmd python3 "$FINAL_PLOTTER" --root-file "$ROOT_FILE"

echo
printf '========== Done ==========%s' ""; echo
