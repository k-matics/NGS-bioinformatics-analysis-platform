"""JSON + Server-Sent-Events endpoints used by the browser for live updates."""
import json
import time
from pathlib import Path

from flask import Blueprint, Response, current_app, jsonify, request, stream_with_context, url_for

from ..extensions import db
from ..models import WorkflowRun
from ..models import status as S
from ..services.file_service import list_directory
from ..services.system_service import system_snapshot
from ..utils.validators import ValidationError
from .dashboard import run_counts

bp = Blueprint("api", __name__, url_prefix="/api")

LOG_CHUNK = 64 * 1024


def run_state(run: WorkflowRun) -> dict:
    current = run.current_step
    last_cmd = None
    if current and current.commands:
        last_cmd = current.commands[-1].display
    return {
        "run_id": run.run_id,
        "status": run.status,
        "progress": run.progress,
        "elapsed": run.elapsed_seconds,
        "error": run.error_message,
        "warnings": run.warnings,
        "control": run.control,
        "current_step": current.name if current else None,
        "current_command": last_cmd,
        "steps": [
            {"key": s.key, "name": s.name, "status": s.status, "message": s.message,
             "duration": s.duration_seconds, "warnings": len(s.warnings or []),
             "has_log": bool(s.log_path)}
            for s in run.steps
        ],
    }


@bp.route("/system")
def system():
    return jsonify(system_snapshot())


@bp.route("/dashboard")
def dashboard():
    active = (WorkflowRun.query.filter(WorkflowRun.status.in_(list(S.ACTIVE_STATUSES)))
              .order_by(WorkflowRun.created_at).all())
    return jsonify({
        "counts": run_counts(),
        "active": [
            {**run_state(r), "sample": r.sample.sample_id, "workflow": r.workflow_type,
             "url": url_for("workflows.detail", run_id=r.run_id)}
            for r in active
        ],
    })


@bp.route("/browse")
def browse():
    try:
        roots = current_app.config["DATA_ROOTS"]
        return jsonify(list_directory(request.args.get("path") or str(roots[0]), roots,
                                      show_all=request.args.get("all") == "1"))
    except ValidationError as exc:
        return jsonify({"error": str(exc)}), 400


@bp.route("/runs/<run_id>/state")
def state(run_id):
    run = WorkflowRun.query.filter_by(run_id=run_id).first_or_404()
    return jsonify(run_state(run))


def _read_new(path: Path, offset: int) -> tuple[str, int]:
    try:
        size = path.stat().st_size
    except OSError:
        return "", offset
    if offset < 0:  # first read: only the tail of a long log
        offset = max(0, size - LOG_CHUNK)
    if size <= offset:
        return "", offset
    with open(path, "rb") as fh:
        fh.seek(offset)
        data = fh.read(LOG_CHUNK)
    return data.decode(errors="replace"), offset + len(data)


@bp.route("/runs/<run_id>/events")
def events(run_id):
    """SSE stream: 'state' events when the run changes, 'log' events with new log text."""

    def generate():
        last_payload, offsets, last_ping = None, {}, time.time()
        yield "retry: 3000\n\n"
        while True:
            run = WorkflowRun.query.filter_by(run_id=run_id).first()
            if run is None:
                yield "event: error\ndata: {\"error\": \"run not found\"}\n\n"
                return
            payload = json.dumps(run_state(run), default=str)
            if payload != last_payload:
                yield f"event: state\ndata: {payload}\n\n"
                last_payload = payload
            for s in run.steps:
                if not s.log_path:
                    continue
                text, offsets[s.key] = _read_new(Path(s.log_path), offsets.get(s.key, -1))
                if text:
                    yield f"event: log\ndata: {json.dumps({'step': s.key, 'text': text})}\n\n"
            terminal = run.status in S.TERMINAL_STATUSES
            db.session.close()  # release the read snapshot so the next poll sees new data
            if terminal:
                yield "event: end\ndata: {}\n\n"
                return
            if time.time() - last_ping > 15:
                yield ": ping\n\n"
                last_ping = time.time()
            time.sleep(1)

    return Response(stream_with_context(generate()), mimetype="text/event-stream",
                    headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})
