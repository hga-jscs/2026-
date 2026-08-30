#!/usr/bin/env bash
set -Eeuo pipefail

# 只在 task2 私有目录中建目录和符号链接；不使用 sudo，也不修改共享环境。
DEPLOY_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXPECTED_ROOT="${TASK2_ROOT:-$HOME/task2_ac_lwcn_20260830}"
BASE_ROOT="${TASK1_ROOT:-$HOME/task1_coral_20260825}"

test "$DEPLOY_ROOT" = "$EXPECTED_ROOT"
test -x "$BASE_ROOT/venv/bin/python"
test -d "$BASE_ROOT/data_root/extracted/6a581fe4b219fa804843a951/results"
test -d "$BASE_ROOT/reference_data/GRCh38"

mkdir -p "$DEPLOY_ROOT/data_root/ac_runs" \
  "$DEPLOY_ROOT/data_root/lp_runs" \
  "$DEPLOY_ROOT/results" \
  "$DEPLOY_ROOT/logs"

test ! -e "$DEPLOY_ROOT/venv" || test -L "$DEPLOY_ROOT/venv"
test ! -e "$DEPLOY_ROOT/reference_data" || test -L "$DEPLOY_ROOT/reference_data"
test ! -e "$DEPLOY_ROOT/data_root/extracted" || test -L "$DEPLOY_ROOT/data_root/extracted"

ln -sfn "$BASE_ROOT/venv" "$DEPLOY_ROOT/venv"
ln -sfn "$BASE_ROOT/reference_data" "$DEPLOY_ROOT/reference_data"
ln -sfn "$BASE_ROOT/data_root/extracted" "$DEPLOY_ROOT/data_root/extracted"

printf 'TASK2_SERVER_SETUP_PASSED root=%s\n' "$DEPLOY_ROOT"
