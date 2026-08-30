#!/usr/bin/env bash
set -Eeuo pipefail

DEPLOY_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="$DEPLOY_ROOT/venv/bin/python"
PROJECT_ID="6a581fe4b219fa804843a951"
INPUT_DIR="$DEPLOY_ROOT/data_root/extracted/$PROJECT_ID"
AC_RUN_DIR="$DEPLOY_ROOT/data_root/ac_runs/$PROJECT_ID"
STATUS_FILE="$DEPLOY_ROOT/results/run_status.txt"

export AA_DATA_REPO="$DEPLOY_ROOT/reference_data"
export PYTHONPATH="$DEPLOY_ROOT/AmpliconClassifier:$DEPLOY_ROOT/algorithm_revised_src:$DEPLOY_ROOT/original_algorithm_src:$DEPLOY_ROOT/runner"
export MPLBACKEND="Agg"
export PYTHONHASHSEED="0"
export OMP_NUM_THREADS="1"
export OPENBLAS_NUM_THREADS="1"
export MKL_NUM_THREADS="1"

test -x "$PYTHON_BIN"
test -d "$INPUT_DIR/results"
mkdir -p "$AC_RUN_DIR" "$DEPLOY_ROOT/data_root/lp_runs" "$DEPLOY_ROOT/results"
printf 'RUNNING\nstarted=%s\n' "$(date --iso-8601=seconds)" > "$STATUS_FILE"

on_exit() {
  exit_code=$?
  trap - EXIT
  if [ "$exit_code" -eq 0 ]; then
    state="COMPLETE"
  else
    state="FAILED_OR_TIMEOUT"
  fi
  printf '%s\nfinished=%s\nexit_code=%s\n' "$state" "$(date --iso-8601=seconds)" "$exit_code" > "$STATUS_FILE"
  exit "$exit_code"
}
trap on_exit EXIT

"$PYTHON_BIN" "$DEPLOY_ROOT/AmpliconClassifier/amplicon_classifier.py" \
  --ref GRCh38 \
  --AA_results "$INPUT_DIR" \
  -o "$AC_RUN_DIR/$PROJECT_ID" \
  --jobs 4 \
  --bfb_threads 1 \
  --no_bfbarchitect \
  --no_results_table

"$PYTHON_BIN" "$DEPLOY_ROOT/runner/run_lwcn_and_merge.py" \
  --dataset-manifest "$DEPLOY_ROOT/runner/coral_manifest.json" \
  --data-root "$DEPLOY_ROOT/data_root" \
  --output-dir "$DEPLOY_ROOT/results" \
  --jobs 4 \
  --force

"$PYTHON_BIN" "$DEPLOY_ROOT/runner/verify_server_task2.py" \
  --deploy-root "$DEPLOY_ROOT"
