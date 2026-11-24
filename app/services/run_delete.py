"""Remove a finished workflow run from the platform, optionally with its output folder.

Two modes, both chosen explicitly by the user on a confirmation page:

* ``records``  - delete the run's database records (run, steps, commands, events). Files on disk
  are kept, and so are the files the run registered to the sample (they still exist).
* ``outputs``  - also delete the run's output folder and un-register the sample files inside it.

Safety rules (any failure -> nothing is deleted):
* only finished runs (completed / failed / cancelled) - never a run a worker may still write to;
* the folder deleted must be the run's own folder: inside RESULTS_DIR, named after the run ID,
  not a symlink. Symlinks inside it are removed as links; their targets (e.g. an external BAM)
  are never touched;
* output files are not deleted while another run still uses one of them as input
  (input_file_ids or baseline_file_id) - delete those runs first, or keep the files;
* every deletion is appended to LOG_DIR/deletions.jsonl (who, when, what, how much) so the
  removal itself stays traceable.
"""
from __future__ import annotations

import getpass
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

from flask import current_app

from ..extensions import db
from ..models import SampleFile, WorkflowRun
from ..models import status as S

MODES = ("records", "outputs")


class DeleteError(Exception):
    pass


def _inside(path: Path, folder: Path) -> bool:
    path, folder = Path(os.path.abspath(path)), Path(os.path.abspath(folder))
    return path == folder or folder in path.parents


def run_folder(run: WorkflowRun) -> Path:
    return Path(os.path.abspath(run.output_dir))


def folder_problem(run: WorkflowRun) -> str | None:
    """Why the run's output folder must NOT be deleted, or None if it is safe."""
    folder = run_folder(run)
    results = Path(os.path.abspath(current_app.config["RESULTS_DIR"]))
    if folder.name != run.run_id:
        return f"The output folder name ({folder.name}) is not the run ID - refusing to delete it."
    if not _inside(folder, results) or folder == results:
        return f"The output folder is not inside the results folder ({results}) - refusing to delete it."
    if folder.is_symlink():
        return "The output folder is a symbolic link - refusing to delete it."
    return None


def folder_usage(folder: Path) -> tuple[int, int, int]:
    """(files, symlinks, bytes) without following symlinks."""
    files = links = size = 0
    if not folder.is_dir():
        return 0, 0, 0
    for root, dirs, names in os.walk(folder, followlinks=False):
        for d in dirs:
            if os.path.islink(os.path.join(root, d)):
                links += 1
        for n in names:
            p = os.path.join(root, n)
            if os.path.islink(p):
                links += 1
            else:
                files += 1
                size += os.lstat(p).st_size
    return files, links, size


def registered_files(run: WorkflowRun) -> list[SampleFile]:
    """Sample files whose path lies inside this run's output folder."""
    folder = run_folder(run)
    return [f for f in SampleFile.query.filter_by(sample_pk=run.sample_pk).all()
            if _inside(Path(f.path), folder)]


def dependent_runs(run: WorkflowRun, files: list[SampleFile]) -> list[tuple[WorkflowRun, list[SampleFile]]]:
    """Other runs that use one of `files` as input."""
    ids = {f.id for f in files}
    by_id = {f.id: f for f in files}
    out = []
    if not ids:
        return out
    for other in WorkflowRun.query.filter(WorkflowRun.id != run.id).all():
        used = set(other.input_file_ids or [])
        bid = (other.params or {}).get("baseline_file_id")
        if bid is not None:
            used.add(bid)
        hit = sorted(used & ids)
        if hit:
            out.append((other, [by_id[i] for i in hit]))
    return out


def deletion_plan(run: WorkflowRun) -> dict:
    """Everything the confirmation page shows. Nothing is changed."""
    folder = run_folder(run)
    files = registered_files(run)
    deps = dependent_runs(run, files)
    n_files, n_links, size = folder_usage(folder)
    blockers = []
    if run.status not in S.TERMINAL_STATUSES:
        blockers.append(f"The run is {run.status}. Cancel it and wait until it has stopped before removing it.")
    output_blockers = []
    problem = folder_problem(run)
    if problem:
        output_blockers.append(problem)
    if deps:
        output_blockers.append(f"{len(deps)} other run(s) use output files of this run as input. "
                               "Remove those runs first, or remove only the record and keep the files.")
    return {"folder": str(folder), "folder_exists": folder.is_dir(), "files": n_files, "symlinks": n_links,
            "bytes": size, "registered": files, "dependents": deps, "blockers": blockers,
            "output_blockers": output_blockers,
            "steps": len(run.steps), "commands": sum(len(s.commands) for s in run.steps),
            "events": len(run.events)}


def _audit(entry: dict) -> None:
    log_dir = Path(current_app.config["LOG_DIR"])
    log_dir.mkdir(parents=True, exist_ok=True)
    with open(log_dir / "deletions.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry, default=str) + "\n")


def delete_run(run: WorkflowRun, mode: str, user: str | None = None) -> dict:
    """Delete the run (mode 'records' or 'outputs'). Raises DeleteError and changes nothing
    if any safety rule fails. Returns what was removed."""
    if mode not in MODES:
        raise DeleteError(f"Unknown mode {mode!r}.")
    plan = deletion_plan(run)
    if plan["blockers"]:
        raise DeleteError(" ".join(plan["blockers"]))
    if mode == "outputs" and plan["output_blockers"]:
        raise DeleteError(" ".join(plan["output_blockers"]))

    folder = run_folder(run)
    entry = {"time": datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "user": user or getpass.getuser(), "mode": mode, "run_id": run.run_id,
             "sample": run.sample.sample_id, "workflow": f"{run.workflow_type} v{run.workflow_version}",
             "status": run.status, "output_dir": str(folder),
             "files_removed": 0, "bytes_removed": 0, "unregistered_files": []}

    # 1. remove the folder first: if that fails, the database still describes what is on disk
    if mode == "outputs":
        if plan["folder_exists"]:
            try:
                shutil.rmtree(folder)          # unlinks symlinks, never follows them
            except OSError as exc:
                entry["error"] = str(exc)
                _audit(entry)
                raise DeleteError(f"Could not delete {folder}: {exc}. Some files may already be gone; "
                                  "the run record was kept so you can see what it was.") from exc
        entry["files_removed"], entry["bytes_removed"] = plan["files"], plan["bytes"]
        for f in plan["registered"]:
            entry["unregistered_files"].append({"id": f.id, "role": f.role, "path": f.path})
            db.session.delete(f)

    # 2. database records (steps, commands and events cascade)
    db.session.delete(run)
    db.session.commit()
    _audit(entry)
    return entry
