#!/usr/bin/env bash
# Stop web app and worker. Redis is stopped only with --all.
# A running workflow is interrupted by stopping the worker; it will be marked
# FAILED (resumable) the next time the worker starts.
set -uo pipefail
APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_DIR="${NGS_CONDA_ENV:-$HOME/miniconda3/envs/ngs_platform}"
cd "$APP_DIR"

stop_pid() {
  local name=$1 file=$2
  if [[ -f $file ]] && kill -0 "$(cat "$file")" 2>/dev/null; then
    echo "Stopping $name (pid $(cat "$file"))"
    kill -TERM "$(cat "$file")"
    for _ in $(seq 1 30); do kill -0 "$(cat "$file")" 2>/dev/null || break; sleep 1; done
  else
    echo "$name not running"
  fi
  rm -f "$file"
}

stop_pid "web app" instance/web.pid
stop_pid "celery worker" instance/worker.pid
if [[ "${1:-}" == "--all" ]]; then
  "$ENV_DIR/bin/redis-cli" -p "${NGS_REDIS_PORT:-6379}" shutdown 2>/dev/null && echo "Redis stopped"
fi
