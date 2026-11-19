#!/usr/bin/env bash
# Start Redis (job queue), the Celery worker (runs pipelines) and the web app.
# Everything runs as the current user, bound to 127.0.0.1 only.
# Access from Windows:  ssh -L 5000:127.0.0.1:5000 <user>@<server>   then open http://localhost:5000
set -euo pipefail

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Server-specific settings (e.g. NGS_RESULTS_DIR, NGS_DATA_ROOTS) live in config/local.env,
# which is not in git. See config/local.env.example.
if [[ -f "$APP_DIR/config/local.env" ]]; then
  set -a; source "$APP_DIR/config/local.env"; set +a
fi
ENV_DIR="${NGS_CONDA_ENV:-$HOME/miniconda3/envs/ngs_platform}"
PORT="${NGS_PORT:-5000}"
REDIS_PORT="${NGS_REDIS_PORT:-6379}"
CONCURRENCY="${NGS_WORKER_CONCURRENCY:-1}"

export PATH="$ENV_DIR/bin:$PATH"
export NGS_REDIS_URL="redis://127.0.0.1:${REDIS_PORT}"
cd "$APP_DIR"
mkdir -p logs instance database/redis

if ! redis-cli -p "$REDIS_PORT" ping >/dev/null 2>&1; then
  echo "Starting Redis on 127.0.0.1:${REDIS_PORT}"
  redis-server --bind 127.0.0.1 --port "$REDIS_PORT" --protected-mode yes \
    --daemonize yes --dir "$APP_DIR/database/redis" --appendonly yes \
    --pidfile "$APP_DIR/instance/redis.pid" --logfile "$APP_DIR/logs/redis.log"
  sleep 1
fi
redis-cli -p "$REDIS_PORT" ping

if [[ -f instance/worker.pid ]] && kill -0 "$(cat instance/worker.pid)" 2>/dev/null; then
  echo "Celery worker already running (pid $(cat instance/worker.pid))"
else
  echo "Starting Celery worker (concurrency=${CONCURRENCY})"
  celery -A celery_worker.celery worker -Q workflows --concurrency="$CONCURRENCY" \
    --hostname "ngs-worker@%h" --loglevel=INFO --logfile "$APP_DIR/logs/worker.log" \
    --pidfile "$APP_DIR/instance/worker.pid" --detach
fi

if [[ -f instance/web.pid ]] && kill -0 "$(cat instance/web.pid)" 2>/dev/null; then
  echo "Web app already running (pid $(cat instance/web.pid))"
else
  echo "Starting web app on http://127.0.0.1:${PORT}"
  # One process + threads: SSE streams each hold a thread; SQLite prefers one writer process.
  gunicorn -k gthread -w 1 --threads 32 -b "127.0.0.1:${PORT}" --timeout 120 \
    --pid "$APP_DIR/instance/web.pid" --access-logfile "$APP_DIR/logs/access.log" \
    --error-logfile "$APP_DIR/logs/web.log" --daemon run:app
fi
sleep 2
"$APP_DIR/scripts/status.sh"
