"""Executes a WorkflowRun step by step inside the Celery worker.

Responsibilities
* claim the run atomically (prevents two workers executing the same run)
* verify inputs have not changed since submission (size + mtime)
* run each step, recording status, timings, commands, exit codes, versions
* skip steps already completed in an earlier attempt **only if their outputs
  still validate** (resume)
* honour cancel / pause requests from the UI
* write a run manifest (provenance) at the end, whatever the outcome
"""
from __future__ import annotations

import json
import logging
import socket
import traceback
from pathlib import Path

from sqlalchemy import select, update

from ..extensions import db
from ..models import CommandExecution, SampleFile, StepRun, WorkflowEvent, WorkflowRun
from ..models import status as S
from ..models.base import utcnow
from ..services.provenance import write_manifest
from ..services.settings_service import get_setting
from ..services.tool_service import tool_version
from .registry import get_workflow
from .runner import CommandFailed, CommandRunner, StepCancelled
from .step import StepContext, StepFailure

log = logging.getLogger(__name__)


def add_event(run: WorkflowRun, message: str, level: str = "INFO") -> None:
    db.session.add(WorkflowEvent(run_pk=run.id, message=message, level=level))


def _read_control(run_pk: int) -> str | None:
    """Read the control flag with a fresh connection (not the session's snapshot)."""
    with db.engine.connect() as conn:
        return conn.execute(select(WorkflowRun.control).where(WorkflowRun.id == run_pk)).scalar()


def claim_run(run_pk: int) -> bool:
    result = db.session.execute(
        update(WorkflowRun)
        .where(WorkflowRun.id == run_pk, WorkflowRun.status == S.QUEUED)
        .values(status=S.RUNNING, started_at=utcnow(), finished_at=None, error_message=None,
                worker_host=socket.gethostname())
    )
    db.session.commit()
    return result.rowcount == 1


def _inputs_changed(files: list[SampleFile]) -> list[str]:
    changed = []
    for f in files:
        p = Path(f.path)
        if not p.exists():
            changed.append(f"{p} no longer exists")
            continue
        st = p.stat()
        if st.st_size != f.size_bytes or abs(st.st_mtime - (f.mtime or 0)) > 1e-3:
            changed.append(f"{p} was modified after registration (size/mtime changed)")
    return changed


def execute_run(run_pk: int) -> str:
    if not claim_run(run_pk):
        log.warning("Run %s not in queued state — not executing (duplicate or cancelled).", run_pk)
        return "skipped"

    run = db.session.get(WorkflowRun, run_pk)
    definition = get_workflow(run.workflow_type)
    files = [db.session.get(SampleFile, fid) for fid in run.input_file_ids]
    files = [f for f in files if f is not None]
    output_dir = Path(run.output_dir)
    (output_dir / "logs").mkdir(parents=True, exist_ok=True)
    add_event(run, f"Execution started on {run.worker_host}.")
    db.session.commit()

    step_impls = {s.key: s for s in definition.build_steps(run.params)}
    ctx = StepContext(
        run=run, sample=run.sample, reference=run.reference, params=run.params,
        output_dir=output_dir, threads=int(run.params.get("threads") or get_setting("threads_per_job", 4)),
        input_files=files, stop_requested=lambda: _read_control(run_pk) == "cancel",
    )
    current: StepRun | None = None
    try:
        problems = _inputs_changed(files)
        if problems:
            raise StepFailure("Input files changed since this run was created; start a new run. "
                              + "; ".join(problems))
        for sr in run.steps:
            impl = step_impls[sr.key]
            if sr.status == S.SKIPPED:
                continue
            if sr.status == S.COMPLETED:
                if impl.outputs_still_valid(ctx, sr.outputs):
                    add_event(run, f"Resume: '{sr.name}' already completed, outputs re-validated — not re-run.")
                    db.session.commit()
                    continue
                add_event(run, f"Resume: outputs of '{sr.name}' no longer valid — re-running.", "WARNING")
            current = sr
            _run_step(run, sr, impl, ctx)
            current = None
        run.status = S.COMPLETED
        add_event(run, "Workflow completed.")
    except StepCancelled:
        run.status = S.CANCELLED
        if current:
            current.status, current.finished_at = S.CANCELLED, utcnow()
        add_event(run, "Workflow cancelled by user. Intermediate files preserved.", "WARNING")
    except (StepFailure, CommandFailed, FileNotFoundError) as exc:
        _fail(run, current, str(exc))
    except Exception as exc:  # unexpected bug: keep full traceback in the step log
        if current and current.log_path:
            with open(current.log_path, "a") as fh:
                fh.write(traceback.format_exc())
        _fail(run, current, f"Internal error: {exc!r}")
    finally:
        run.finished_at = utcnow()
        run.control = None
        db.session.commit()
        try:
            write_manifest(run)
        except Exception:
            log.exception("Could not write manifest for %s", run.run_id)
    return run.status


def _fail(run: WorkflowRun, step: StepRun | None, message: str) -> None:
    run.status = S.FAILED
    run.error_message = message
    if step:
        step.status, step.finished_at, step.message = S.FAILED, utcnow(), message
    add_event(run, f"Workflow failed{' at ' + step.name if step else ''}: {message}", "ERROR")


def _run_step(run: WorkflowRun, sr: StepRun, impl, ctx: StepContext) -> None:
    sr.status = S.RUNNING
    sr.started_at, sr.finished_at = utcnow(), None
    sr.attempts += 1
    sr.message, sr.warnings, sr.summary, sr.outputs = None, [], {}, {}
    sr.log_path = str(ctx.output_dir / "logs" / f"{sr.order:02d}_{sr.key}.log")
    ctx.log_path = Path(sr.log_path)
    add_event(run, f"Step started: {sr.name}")
    db.session.commit()
    ctx.log(f"=== {sr.name} (attempt {sr.attempts}) — run {run.run_id}, sample {run.sample.sample_id} ===")

    ok, reason = impl.applicable(ctx)
    if not ok:
        sr.status, sr.message, sr.finished_at = S.SKIPPED, reason, utcnow()
        ctx.log(f"Skipped: {reason}")
        add_event(run, f"Step skipped: {sr.name} — {reason}")
        db.session.commit()
        return

    def on_start(argv, display, pid):
        versions = {}
        for cmd in argv:
            name = Path(cmd[0]).name
            versions[name] = tool_version(name, cmd[0]) or "unknown"
        rec = CommandExecution(step_pk=sr.id, argv=argv, display=display, tool_versions=versions, pid=pid)
        db.session.add(rec)
        db.session.commit()
        return rec

    def on_finish(rec, codes):
        rec.exit_codes, rec.finished_at = codes, utcnow()
        db.session.commit()

    def on_pause_change(paused: bool):
        run.status = S.PAUSED if paused else S.RUNNING
        sr.status = S.PAUSED if paused else S.RUNNING
        add_event(run, "Paused." if paused else "Resumed.")
        db.session.commit()

    runner = CommandRunner(
        Path(sr.log_path), control=lambda: _read_control(run.id),
        on_start=on_start, on_finish=on_finish, on_pause_change=on_pause_change,
    )
    result = impl.execute(ctx, runner)
    ctx.check_cancel()
    impl.validate(ctx, result)
    sr.outputs, sr.summary, sr.warnings = result.outputs, result.summary, result.warnings
    if result.warnings:
        run.warnings = list(run.warnings or []) + [f"{sr.name}: {w}" for w in result.warnings]
        for w in result.warnings:
            ctx.log(f"WARNING: {w}")
    sr.status, sr.finished_at = S.COMPLETED, utcnow()
    ctx.log(f"Step completed. Outputs: {json.dumps(result.outputs)}")
    add_event(run, f"Step completed: {sr.name}" + (f" ({len(result.warnings)} warning(s))" if result.warnings else ""))
    db.session.commit()


def mark_interrupted_runs() -> int:
    """Called when a worker starts: runs left RUNNING/PAUSED by a dead worker
    are marked FAILED (resumable), so they never look 'running' forever."""
    stale = WorkflowRun.query.filter(WorkflowRun.status.in_([S.RUNNING, S.PAUSED])).all()
    for run in stale:
        run.status = S.FAILED
        run.error_message = "Worker stopped while this run was executing (server restart or crash). Use Resume."
        run.finished_at = utcnow()
        for sr in run.steps:
            if sr.status in (S.RUNNING, S.PAUSED):
                sr.status, sr.message = S.FAILED, "Interrupted"
        add_event(run, run.error_message, "ERROR")
    db.session.commit()
    return len(stale)
