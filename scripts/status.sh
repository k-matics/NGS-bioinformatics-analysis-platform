#!/usr/bin/env bash
APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_DIR="${NGS_CONDA_ENV:-$HOME/miniconda3/envs/ngs_platform}"
cd "$APP_DIR"
check() { [[ -f $2 ]] && kill -0 "$(cat "$2")" 2>/dev/null && echo "  $1: running (pid $(cat "$2"))" || echo "  $1: STOPPED"; }
echo "NGS platform services:"
"$ENV_DIR/bin/redis-cli" -p "${NGS_REDIS_PORT:-6379}" ping >/dev/null 2>&1 && echo "  redis: running" || echo "  redis: STOPPED"
check "celery worker" instance/worker.pid
check "web app" instance/web.pid
if [[ -f "$APP_DIR/config/local.env" ]]; then set -a; source "$APP_DIR/config/local.env"; set +a; fi
RESULTS="${NGS_RESULTS_DIR:-$APP_DIR/results}"
echo "Results directory: $RESULTS"
echo "Disk free on results filesystem: $(df -h "$RESULTS" | awk 'NR==2{print $4" ("$5" used)"}')"
