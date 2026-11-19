"""Workflow / step status vocabulary (single source of truth)."""

PENDING = "pending"      # created, not yet sent to the queue
QUEUED = "queued"        # in the job queue, waiting for a worker
RUNNING = "running"
PAUSED = "paused"        # processes stopped with SIGSTOP, resumable
COMPLETED = "completed"
FAILED = "failed"
CANCELLED = "cancelled"
SKIPPED = "skipped"      # step not applicable (e.g. no BAM registered)

ACTIVE_STATUSES = {PENDING, QUEUED, RUNNING, PAUSED}
TERMINAL_STATUSES = {COMPLETED, FAILED, CANCELLED}
RUN_STATUSES = [PENDING, QUEUED, RUNNING, PAUSED, COMPLETED, FAILED, CANCELLED]

STATUS_BADGE = {
    PENDING: "secondary",
    QUEUED: "info",
    RUNNING: "primary",
    PAUSED: "warning",
    COMPLETED: "success",
    FAILED: "danger",
    CANCELLED: "dark",
    SKIPPED: "light",
}
