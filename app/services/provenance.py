"""Run manifest: a self-contained JSON record of what was run, on what, with what.

It is written into every run directory so results remain traceable even
without the database.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

from ..extensions import db
from ..models import SampleFile


def _iso(dt):
    return dt.isoformat() + "Z" if dt else None


def build_manifest(run) -> dict:
    files = [db.session.get(SampleFile, fid) for fid in run.input_file_ids]
    ref = run.reference
    return {
        "manifest_version": 1,
        "run_id": run.run_id,
        "workflow": {"type": run.workflow_type, "version": run.workflow_version},
        "sample_id": run.sample.sample_id,
        "status": run.status,
        "error_message": run.error_message,
        "warnings": run.warnings,
        "submitted_by": run.submitted_by,
        "worker_host": run.worker_host,
        "created_at": _iso(run.created_at),
        "started_at": _iso(run.started_at),
        "finished_at": _iso(run.finished_at),
        "parameters": run.params,
        "reference": None if ref is None else {
            "name": ref.name, "fasta": ref.fasta_path, "build": ref.build,
            "contigs": ref.contig_count, "status": ref.status,
        },
        "inputs": [
            {"role": f.role, "origin": f.origin, "path": f.path, "size_bytes": f.size_bytes,
             "sha256": f.sha256, "detected_build": f.detected_build}
            for f in files if f is not None
        ],
        "steps": [
            {
                "order": s.order, "key": s.key, "name": s.name, "status": s.status,
                "started_at": _iso(s.started_at), "finished_at": _iso(s.finished_at),
                "attempts": s.attempts, "message": s.message, "warnings": s.warnings,
                "outputs": s.outputs, "summary": s.summary, "log": s.log_path,
                "commands": [
                    {"command": c.display, "argv": c.argv, "tool_versions": c.tool_versions,
                     "exit_codes": c.exit_codes, "started_at": _iso(c.started_at),
                     "finished_at": _iso(c.finished_at)}
                    for c in s.commands
                ],
            }
            for s in run.steps
        ],
    }


def write_manifest(run) -> Path:
    path = Path(run.output_dir) / "run_manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(build_manifest(run), indent=2, default=str))
    os.replace(tmp, path)  # atomic: never leaves a half-written manifest
    return path
