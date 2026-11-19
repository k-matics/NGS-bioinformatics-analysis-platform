"""Workflow provenance model.

WorkflowRun     one execution of a workflow for one sample (unique run_id)
  StepRun       one processing step inside the run (FASTQ QC, alignment, ...)
    CommandExecution  every external command the step ran: argv, exit code, timing, tool version
  WorkflowEvent human-readable audit trail (submitted, started, cancelled, ...)
"""
from ..extensions import db
from . import status as S
from .base import utcnow


class WorkflowRun(db.Model):
    __tablename__ = "workflow_runs"

    id = db.Column(db.Integer, primary_key=True)
    run_id = db.Column(db.String(40), unique=True, nullable=False, index=True)
    sample_pk = db.Column(db.Integer, db.ForeignKey("samples.id"), nullable=False)
    reference_pk = db.Column(db.Integer, db.ForeignKey("reference_genomes.id"))
    workflow_type = db.Column(db.String(64), nullable=False)
    workflow_version = db.Column(db.String(16), nullable=False)
    status = db.Column(db.String(16), nullable=False, default=S.PENDING, index=True)
    params = db.Column(db.JSON, nullable=False, default=dict)
    input_file_ids = db.Column(db.JSON, nullable=False, default=list)
    fingerprint = db.Column(db.String(64), index=True)  # duplicate-execution guard
    output_dir = db.Column(db.String(1024), nullable=False)
    celery_task_id = db.Column(db.String(64))
    control = db.Column(db.String(16))  # cancel / pause / resume requests from the UI
    submitted_by = db.Column(db.String(64))
    worker_host = db.Column(db.String(128))
    created_at = db.Column(db.DateTime, default=utcnow, nullable=False)
    started_at = db.Column(db.DateTime)
    finished_at = db.Column(db.DateTime)
    error_message = db.Column(db.Text)
    warnings = db.Column(db.JSON, nullable=False, default=list)

    sample = db.relationship("Sample", back_populates="runs")
    reference = db.relationship("ReferenceGenome")
    steps = db.relationship(
        "StepRun", back_populates="run", cascade="all, delete-orphan", order_by="StepRun.order"
    )
    events = db.relationship(
        "WorkflowEvent", back_populates="run", cascade="all, delete-orphan",
        order_by="WorkflowEvent.id",
    )

    # --- derived values (computed from real step states, never simulated) ---
    @property
    def progress(self) -> int:
        """Percentage of steps finished (completed or skipped)."""
        applicable = [s for s in self.steps]
        if not applicable:
            return 0
        done = sum(1 for s in applicable if s.status in (S.COMPLETED, S.SKIPPED))
        return round(100 * done / len(applicable))

    @property
    def current_step(self):
        for s in self.steps:
            if s.status in (S.RUNNING, S.PAUSED):
                return s
        return None

    @property
    def elapsed_seconds(self):
        if not self.started_at:
            return None
        end = self.finished_at or utcnow()
        return (end - self.started_at).total_seconds()

    @property
    def is_active(self) -> bool:
        return self.status in S.ACTIVE_STATUSES

    @property
    def can_resume(self) -> bool:
        return self.status in (S.FAILED, S.CANCELLED)


class StepRun(db.Model):
    __tablename__ = "step_runs"

    id = db.Column(db.Integer, primary_key=True)
    run_pk = db.Column(db.Integer, db.ForeignKey("workflow_runs.id"), nullable=False)
    order = db.Column(db.Integer, nullable=False)
    key = db.Column(db.String(64), nullable=False)
    name = db.Column(db.String(128), nullable=False)
    status = db.Column(db.String(16), nullable=False, default=S.PENDING)
    started_at = db.Column(db.DateTime)
    finished_at = db.Column(db.DateTime)
    log_path = db.Column(db.String(1024))
    outputs = db.Column(db.JSON, nullable=False, default=dict)   # name -> path
    summary = db.Column(db.JSON, nullable=False, default=dict)   # key metrics for the UI
    warnings = db.Column(db.JSON, nullable=False, default=list)
    message = db.Column(db.Text)                                 # failure reason / skip reason
    attempts = db.Column(db.Integer, nullable=False, default=0)

    run = db.relationship("WorkflowRun", back_populates="steps")
    commands = db.relationship(
        "CommandExecution", back_populates="step", cascade="all, delete-orphan",
        order_by="CommandExecution.id",
    )

    @property
    def duration_seconds(self):
        if not self.started_at:
            return None
        end = self.finished_at or utcnow()
        return (end - self.started_at).total_seconds()


class CommandExecution(db.Model):
    """One external command (or pipe of commands) run by a step."""

    __tablename__ = "command_executions"

    id = db.Column(db.Integer, primary_key=True)
    step_pk = db.Column(db.Integer, db.ForeignKey("step_runs.id"), nullable=False)
    argv = db.Column(db.JSON, nullable=False)          # list of argument lists (pipeline)
    display = db.Column(db.Text, nullable=False)       # shell-quoted, for humans only
    tool_versions = db.Column(db.JSON, nullable=False, default=dict)
    started_at = db.Column(db.DateTime, default=utcnow, nullable=False)
    finished_at = db.Column(db.DateTime)
    exit_codes = db.Column(db.JSON)                    # one per process in the pipeline
    pid = db.Column(db.Integer)

    step = db.relationship("StepRun", back_populates="commands")

    @property
    def duration_seconds(self):
        if not self.finished_at:
            return None
        return (self.finished_at - self.started_at).total_seconds()


class WorkflowEvent(db.Model):
    __tablename__ = "workflow_events"

    id = db.Column(db.Integer, primary_key=True)
    run_pk = db.Column(db.Integer, db.ForeignKey("workflow_runs.id"), nullable=False, index=True)
    timestamp = db.Column(db.DateTime, default=utcnow, nullable=False)
    level = db.Column(db.String(8), nullable=False, default="INFO")
    message = db.Column(db.Text, nullable=False)

    run = db.relationship("WorkflowRun", back_populates="events")
