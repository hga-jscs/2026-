#!/usr/bin/env bash
set -Eeuo pipefail

DEPLOY_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="$DEPLOY_ROOT/venv/bin/python"
PROJECT_ID="6a581fe4b219fa804843a951"
INPUT_DIR="$DEPLOY_ROOT/data_root/extracted/$PROJECT_ID"
STATUS_FILE="$DEPLOY_ROOT/results/run_status.txt"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"

export AA_DATA_REPO="$DEPLOY_ROOT/reference_data"
export PYTHONPATH="$DEPLOY_ROOT/AmpliconClassifier:$DEPLOY_ROOT/algorithm_revised_src:$DEPLOY_ROOT/original_algorithm_src:$DEPLOY_ROOT/runner"
export MPLBACKEND="Agg"
export PYTHONHASHSEED="0"
export OMP_NUM_THREADS="1"
export OPENBLAS_NUM_THREADS="1"
export MKL_NUM_THREADS="1"

test -x "$PYTHON_BIN"
test -d "$INPUT_DIR/results"
mkdir -p "$DEPLOY_ROOT/data_root/ac_runs/$PROJECT_ID" "$DEPLOY_ROOT/data_root/lp_runs" "$DEPLOY_ROOT/results"
printf 'RUNNING\nrun_id=%s\nhostname=%s\nstarted=%s\n' \
  "$RUN_ID" "$(hostname)" "$(date --iso-8601=seconds)" > "$STATUS_FILE"

on_exit() {
  exit_code=$?
  trap - EXIT
  if [ "$exit_code" -eq 0 ]; then
    state="COMPLETE"
  else
    state="FAILED_OR_TIMEOUT"
  fi
  printf '%s\nrun_id=%s\nhostname=%s\nfinished=%s\nexit_code=%s\n' \
    "$state" "$RUN_ID" "$(hostname)" "$(date --iso-8601=seconds)" "$exit_code" > "$STATUS_FILE"
  exit "$exit_code"
}
trap on_exit EXIT

chmod u+x "$DEPLOY_ROOT/AmpliconClassifier/ampclasslib/make_input.sh"

"$PYTHON_BIN" "$DEPLOY_ROOT/runner/run_ac_on_linux.py" \
  --deploy-root "$DEPLOY_ROOT" \
  --project-id "$PROJECT_ID" \
  --reference GRCh38 \
  --jobs 4

"$PYTHON_BIN" "$DEPLOY_ROOT/runner/run_lwcn_and_merge.py" \
  --dataset-manifest "$DEPLOY_ROOT/runner/coral_manifest.json" \
  --data-root "$DEPLOY_ROOT/data_root" \
  --output-dir "$DEPLOY_ROOT/results" \
  --jobs 4 \
  --force

cp "$DEPLOY_ROOT/results/全量AC与环状LWCN结果.csv" \
  "$DEPLOY_ROOT/results/server_coral_112_ac_lwcn_ratio.csv"

"$PYTHON_BIN" "$DEPLOY_ROOT/runner/verify_server_task2.py" \
  --deploy-root "$DEPLOY_ROOT"
