"""Create, submit, cancel, pause and resume workflow runs (called by the web routes)."""
from __future__ import annotations

import copy
import getpass
import hashlib
import json
import secrets
from datetime import datetime
from pathlib import Path

from flask import current_app

from ..engine.registry import get_workflow
from ..extensions import db
from ..models import Sample, SampleFile, StepRun, WorkflowRun
from ..models import status as S
from ..models.base import utcnow
from .file_service import disk_free_gb
from .settings_service import get_setting


class WorkflowError(Exception):
    """User-facing error (shown as a flash message)."""


def new_run_id() -> str:
    return f"WF-{datetime.now():%Y%m%d-%H%M%S}-{secrets.token_hex(3)}"


def fingerprint(workflow_type: str, sample: Sample, files: list[SampleFile], params: dict, reference_pk) -> str:
    payload = {
        "workflow": workflow_type,
        "sample": sample.sample_id,
        "inputs": sorted(f.path for f in files),
        "params": params,
        "reference": reference_pk,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


def create_run(sample: Sample, workflow_type: str, files: list[SampleFile], reference,
               params: dict | None = None, submitted_by: str | None = None, enqueue: bool = True) -> WorkflowRun:
    definition = get_workflow(workflow_type)
    merged = copy.deepcopy(definition.default_params)
    merged.update(copy.deepcopy(params or {}))
    merged.setdefault("threads", int(get_setting("threads_per_job", 4)))

    if definition.input_roles is not None:
        # only the file types this workflow consumes become run inputs
        files = [f for f in files if f.role in definition.input_roles]
    problems = definition.check_inputs(sample, files)
    if problems:
        raise WorkflowError("; ".join(problems))
    if reference is None:
        raise WorkflowError("No reference genome selected.")
    problems = definition.validate_params(merged, reference)
    if problems:
        raise WorkflowError("; ".join(problems))
    if reference.status == "error":
        raise WorkflowError(f"Reference '{reference.name}' failed its integrity check — fix it in Settings.")

    fp = fingerprint(workflow_type, sample, files, merged, reference.id)
    duplicate = WorkflowRun.query.filter(
        WorkflowRun.fingerprint == fp, WorkflowRun.status.in_(S.ACTIVE_STATUSES)
    ).first()
    if duplicate:
        raise WorkflowError(f"An identical run is already active: {duplicate.run_id}. "
                            "Duplicate execution was prevented.")

    results_dir: Path = current_app.config["RESULTS_DIR"]
    needed = definition.estimate_disk_gb(sample, files) + float(get_setting("min_free_disk_gb", 3))
    free = disk_free_gb(results_dir)
    if free < needed:
        raise WorkflowError(f"Not enough disk space: {free:.1f} GB free, about {needed:.1f} GB needed "
                            "(estimate + safety margin). Free space or change the margin in Settings.")

    run_id = new_run_id()
    output_dir = results_dir / sample.sample_id / run_id
    output_dir.mkdir(parents=True, exist_ok=False)  # never reuse/overwrite a run directory

    run = WorkflowRun(
        run_id=run_id, sample=sample, reference=reference, workflow_type=workflow_type,
        workflow_version=definition.version, status=S.PENDING, params=merged,
        input_file_ids=[f.id for f in files], fingerprint=fp, output_dir=str(output_dir),
        submitted_by=submitted_by or getpass.getuser(),
    )
    db.session.add(run)
    db.session.flush()
    for i, step in enumerate(definition.build_steps(merged), start=1):
        db.session.add(StepRun(run_pk=run.id, order=i, key=step.key, name=step.name))
    from ..engine.executor import add_event
    add_event(run, f"Run created: {definition.name} for sample {sample.sample_id}.")
    db.session.commit()
    if enqueue:
        enqueue_run(run)
    return run


def enqueue_run(run: WorkflowRun) -> None:
    from ..tasks import run_workflow
    from ..engine.executor import add_event

    run.status = S.QUEUED
    run.control = None
    add_event(run, "Queued for execution.")
    db.session.commit()
    try:
        result = run_workflow.delay(run.id)
    except Exception as exc:  # Redis down etc.: do not pretend the job is queued
        run.status = S.FAILED
        run.error_message = f"Could not reach the job queue (is Redis running?): {exc}"
        add_event(run, run.error_message, "ERROR")
        db.session.commit()
        raise WorkflowError(run.error_message) from exc
    db.session.refresh(run)
    if run.celery_task_id is None and result is not None:
        run.celery_task_id = result.id
        db.session.commit()


def request_cancel(run: WorkflowRun) -> None:
    from ..engine.executor import add_event

    if run.status in (S.PENDING, S.QUEUED):
        # Not started yet: the executor's atomic claim will refuse to run it.
        run.status = S.CANCELLED
        run.finished_at = utcnow()
        for sr in run.steps:
            if sr.status == S.PENDING:
                sr.status = S.CANCELLED
        add_event(run, "Cancelled before execution started.", "WARNING")
    elif run.status in (S.RUNNING, S.PAUSED):
        run.control = "cancel"
        add_event(run, "Cancellation requested.", "WARNING")
    else:
        raise WorkflowError(f"Run is {run.status}; nothing to cancel.")
    db.session.commit()


def request_pause(run: WorkflowRun) -> None:
    if run.status != S.RUNNING:
        raise WorkflowError("Only a running workflow can be paused.")
    run.control = "pause"
    db.session.commit()


def request_unpause(run: WorkflowRun) -> None:
    if run.status != S.PAUSED and run.control != "pause":
        raise WorkflowError("Run is not paused.")
    run.control = None
    db.session.commit()


def resume_run(run: WorkflowRun) -> None:
    """Re-queue a failed/cancelled run. Completed steps are kept (and re-validated
    by the executor); failed/cancelled/pending steps run again."""
    from ..engine.executor import add_event

    if not run.can_resume:
        raise WorkflowError(f"Run is {run.status}; only failed or cancelled runs can be resumed.")
    current_version = get_workflow(run.workflow_type).version
    if current_version != run.workflow_version:
        raise WorkflowError(
            f"The workflow changed since this run was created (v{run.workflow_version} → v{current_version}). "
            "Resuming would mix results of two pipeline versions; start a new run instead.")
    duplicate = WorkflowRun.query.filter(
        WorkflowRun.fingerprint == run.fingerprint, WorkflowRun.status.in_(S.ACTIVE_STATUSES),
        WorkflowRun.id != run.id,
    ).first()
    if duplicate:
        raise WorkflowError(f"An identical run is already active: {duplicate.run_id}.")
    for sr in run.steps:
        if sr.status not in (S.COMPLETED, S.SKIPPED):
            sr.status, sr.message = S.PENDING, None
    run.error_message = None
    add_event(run, "Resume requested: completed steps will be re-validated, remaining steps re-run.")
    db.session.commit()
    enqueue_run(run)
