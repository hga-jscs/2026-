#!/usr/bin/env bash
set -Eeuo pipefail

DEPLOY_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_FILE="$DEPLOY_ROOT/logs/run_$(date -u +%Y%m%dT%H%M%SZ).log"
mkdir -p "$DEPLOY_ROOT/logs"

set +e
timeout --signal=TERM --kill-after=2m 30m \
  nice -n 10 bash "$DEPLOY_ROOT/runner/run_coral_task2.sh" \
  2>&1 | tee "$LOG_FILE"
status=${PIPESTATUS[0]}
set -e

printf 'wrapper_exit_code=%s\nlog=%s\n' "$status" "$LOG_FILE"
exit "$status"
