"""Integration tests: full intake workflow through the engine (Celery in eager mode)."""
import gzip
import json
from pathlib import Path

import pytest

from app.engine import executor
from app.extensions import db
from app.models import SampleFile, WorkflowRun
from app.models import status as S
from app.services.workflow_service import WorkflowError, create_run, request_cancel, resume_run
from conftest import GRCH37, GRCH38, write_bam, write_fastq_pair


def _run(sample, reference, **kw):
    run = create_run(sample, "intake_validation", list(sample.files), reference, **kw)
    db.session.expire_all()
    return db.session.get(WorkflowRun, run.id)


def test_full_intake_success(make_sample, reference, data_dir):
    r1, r2 = write_fastq_pair(data_dir)
    bam = write_bam(data_dir / "ours.bam", GRCH38)
    sample = make_sample(files=[("fastq_r1", r1, "internal"), ("fastq_r2", r2, "internal"),
                                ("bam", bam, "internal")])
    run = _run(sample, reference)
    assert run.status == S.COMPLETED, run.error_message
    assert run.progress == 100
    statuses = {s.key: s.status for s in run.steps}
    assert statuses["vcf_validation"] == S.SKIPPED  # no VCF registered
    assert statuses["fastq_integrity"] == S.COMPLETED
    # checksums recorded and a sha256sum-compatible file written
    for f in sample.files:
        assert f.sha256 and len(f.sha256) == 64
    out = Path(run.output_dir)
    assert (out / "intake" / "checksums.sha256").exists()
    assert (out / "reports" / "intake_report.html").exists()
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["status"] == "completed" and manifest["inputs"][0]["sha256"]
    # every external command recorded with exit codes + versions
    cmds = [c for s in run.steps for c in s.commands]
    assert cmds and all(c.exit_codes == [0] for c in cmds)
    assert any("samtools" in c.tool_versions for c in cmds)
    # reference ALT warning surfaces as a run warning; BAM without index warns
    assert any("ALT" in w for w in run.warnings)
    assert any("no index" in w for w in run.warnings)


def test_grch37_bam_flagged_incompatible(make_sample, reference, data_dir):
    bam = write_bam(data_dir / "external_hg19.bam", GRCH37)
    sample = make_sample(files=[("bam", bam, "external")])
    run = _run(sample, reference)
    assert run.status == S.COMPLETED
    f = SampleFile.query.one()
    assert f.detected_build == "GRCh37" and f.validation_status == "error"
    assert any("GRCh37" in w for w in run.warnings)


def test_bam_without_read_group_warns(make_sample, reference, data_dir):
    bam = write_bam(data_dir / "norg.bam", GRCH38, read_group=False)
    run = _run(make_sample(files=[("bam", bam, "internal")]), reference)
    assert any("@RG" in w for w in run.warnings)


def test_truncated_gzip_fails_and_resume_reruns_only_failed(make_sample, reference, data_dir):
    r1, r2 = write_fastq_pair(data_dir)
    data = r2.read_bytes()
    r2.write_bytes(data[: len(data) // 2])  # truncate
    sample = make_sample(files=[("fastq_r1", r1, "internal"), ("fastq_r2", r2, "internal")])
    run = _run(sample, reference)
    assert run.status == S.FAILED
    failed = [s for s in run.steps if s.status == S.FAILED]
    assert failed and failed[0].key == "fastq_integrity"
    assert "pigz" in run.error_message or "exit codes" in run.error_message

    # Resume without fixing: inputs are unchanged, so fails again at the same step,
    # but the completed steps are NOT re-executed.
    checksum_step = next(s for s in run.steps if s.key == "checksums")
    attempts_before = checksum_step.attempts
    resume_run(run)
    db.session.expire_all()
    run = db.session.get(WorkflowRun, run.id)
    assert run.status == S.FAILED
    assert next(s for s in run.steps if s.key == "checksums").attempts == attempts_before


def test_modified_input_blocks_resume(make_sample, reference, data_dir):
    r1, r2 = write_fastq_pair(data_dir)
    r2.write_bytes(r2.read_bytes()[:50])
    sample = make_sample(files=[("fastq_r1", r1, "internal"), ("fastq_r2", r2, "internal")])
    run = _run(sample, reference)
    assert run.status == S.FAILED
    write_fastq_pair(data_dir)  # "fix" the file in place -> content changed since registration
    resume_run(run)
    db.session.expire_all()
    run = db.session.get(WorkflowRun, run.id)
    assert run.status == S.FAILED and "changed" in run.error_message


def test_unpaired_reads_detected(make_sample, reference, data_dir):
    r1, r2 = write_fastq_pair(data_dir, broken_r2=True)
    run = _run(make_sample(files=[("fastq_r1", r1, "internal"), ("fastq_r2", r2, "internal")]), reference)
    assert run.status == S.FAILED and "out of sync" in run.error_message


def test_malformed_fastq(make_sample, reference, data_dir):
    bad = data_dir / "bad_R1.fastq.gz"
    with gzip.open(bad, "wt") as fh:
        fh.write("@r1\nACGT\n+\nII\n")
    run = _run(make_sample(files=[("fastq_r1", bad, "internal")]), reference)
    assert run.status == S.FAILED and "length" in run.error_message


def test_duplicate_active_run_prevented(make_sample, reference, data_dir):
    bam = write_bam(data_dir / "x.bam", GRCH38)
    sample = make_sample(files=[("bam", bam, "internal")])
    create_run(sample, "intake_validation", list(sample.files), reference, enqueue=False)
    with pytest.raises(WorkflowError, match="already active"):
        create_run(sample, "intake_validation", list(sample.files), reference, enqueue=False)


def test_cancel_before_start_is_never_executed(make_sample, reference, data_dir):
    bam = write_bam(data_dir / "x.bam", GRCH38)
    sample = make_sample(files=[("bam", bam, "internal")])
    run = create_run(sample, "intake_validation", list(sample.files), reference, enqueue=False)
    run.status = S.QUEUED
    db.session.commit()
    request_cancel(run)
    assert executor.execute_run(run.id) == "skipped"
    assert db.session.get(WorkflowRun, run.id).status == S.CANCELLED


def test_interrupted_runs_marked_failed(make_sample, reference, data_dir):
    bam = write_bam(data_dir / "x.bam", GRCH38)
    sample = make_sample(files=[("bam", bam, "internal")])
    run = create_run(sample, "intake_validation", list(sample.files), reference, enqueue=False)
    run.status = S.RUNNING
    db.session.commit()
    assert executor.mark_interrupted_runs() == 1
    assert run.status == S.FAILED and run.can_resume


def test_no_disk_space_refused(app, make_sample, reference, data_dir):
    bam = write_bam(data_dir / "x.bam", GRCH38)
    sample = make_sample(files=[("bam", bam, "internal")])
    app.config["DEFAULTS"]["min_free_disk_gb"] = 10**6
    with pytest.raises(WorkflowError, match="disk space"):
        create_run(sample, "intake_validation", list(sample.files), reference, enqueue=False)


def test_run_directories_unique(make_sample, reference, data_dir):
    bam = write_bam(data_dir / "x.bam", GRCH38)
    sample = make_sample(files=[("bam", bam, "internal")])
    a = _run(sample, reference)
    b = _run(sample, reference)  # first is finished, so a new run is allowed
    assert a.output_dir != b.output_dir and Path(a.output_dir).exists()
