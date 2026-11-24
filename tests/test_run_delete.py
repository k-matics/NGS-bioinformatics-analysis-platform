"""Removing a workflow run: record only, record + outputs, and the safety rules."""
import json
from pathlib import Path

from app.extensions import db
from app.models import SampleFile, StepRun, WorkflowEvent, WorkflowRun
from app.services.run_delete import DeleteError, delete_run, deletion_plan


def _run(app, sample, reference, run_id, status="completed", inputs=(), params=None, folder=None):
    out = Path(folder) if folder else Path(app.config["RESULTS_DIR"]) / sample.sample_id / run_id
    out.mkdir(parents=True, exist_ok=True)
    run = WorkflowRun(run_id=run_id, sample=sample, reference=reference, workflow_type="germline",
                      workflow_version="1.3", status=status, params=params or {},
                      input_file_ids=list(inputs), output_dir=str(out))
    run.steps.append(StepRun(order=1, key="x", name="x", status="completed"))
    run.events.append(WorkflowEvent(message="submitted"))
    db.session.add(run)
    db.session.commit()
    return run, out


def _with_vcf(app, make_sample, reference, data_dir, run_id="WF-1"):
    sample = make_sample()
    run, out = _run(app, sample, reference, run_id)
    (out / "vcf").mkdir()
    vcf = out / "vcf" / "S1.vcf.gz"
    vcf.write_bytes(b"x" * 100)
    external = data_dir / "external.bam"
    external.write_bytes(b"bam")
    (out / "bam_link.bam").symlink_to(external)            # e.g. an external BAM linked into the run
    f = SampleFile(sample=sample, role="vcf", origin="pipeline", path=str(vcf))
    db.session.add(f)
    db.session.commit()
    return sample, run, out, f, external


def test_delete_with_outputs(app, make_sample, reference, data_dir):
    sample, run, out, f, external = _with_vcf(app, make_sample, reference, data_dir)
    plan = deletion_plan(run)
    assert plan["files"] == 1 and plan["symlinks"] == 1 and plan["registered"] == [f]
    entry = delete_run(run, "outputs", user="tester")
    assert not out.exists()
    assert external.exists()                               # link target untouched
    assert WorkflowRun.query.count() == 0 and StepRun.query.count() == 0 and WorkflowEvent.query.count() == 0
    assert SampleFile.query.filter_by(role="vcf").count() == 0
    assert entry["bytes_removed"] == 100
    log = Path(app.config["LOG_DIR"]) / "deletions.jsonl"
    assert json.loads(log.read_text().splitlines()[-1])["run_id"] == "WF-1"


def test_delete_record_only_keeps_files(app, make_sample, reference, data_dir):
    sample, run, out, f, _ = _with_vcf(app, make_sample, reference, data_dir)
    delete_run(run, "records")
    assert (out / "vcf" / "S1.vcf.gz").exists()
    assert WorkflowRun.query.count() == 0
    assert SampleFile.query.filter_by(role="vcf").count() == 1


def test_active_run_refused(app, make_sample, reference):
    run, out = _run(app, make_sample(), reference, "WF-2", status="running")
    for mode in ("records", "outputs"):
        try:
            delete_run(run, mode)
            raise AssertionError("should refuse")
        except DeleteError:
            pass
    assert out.exists() and WorkflowRun.query.count() == 1


def test_outputs_used_by_other_run_refused(app, make_sample, reference, data_dir):
    sample, run, out, f, _ = _with_vcf(app, make_sample, reference, data_dir)
    _run(app, sample, reference, "WF-3", params={"baseline_file_id": f.id})
    assert deletion_plan(run)["dependents"]
    try:
        delete_run(run, "outputs")
        raise AssertionError("should refuse")
    except DeleteError:
        pass
    assert out.exists()
    delete_run(run, "records")                              # record-only is still allowed
    assert WorkflowRun.query.filter_by(run_id="WF-1").first() is None


def test_folder_outside_results_refused(app, make_sample, reference, data_dir):
    run, out = _run(app, make_sample(), reference, "WF-4", folder=data_dir / "WF-4")
    try:
        delete_run(run, "outputs")
        raise AssertionError("should refuse")
    except DeleteError:
        pass
    assert out.exists()


def test_remove_routes(client, app, make_sample, reference, data_dir):
    sample, run, out, f, _ = _with_vcf(app, make_sample, reference, data_dir)
    assert b"delete the output files" in client.get("/workflows/WF-1/remove").data
    assert b"Remove run" in client.get("/workflows/WF-1").data
    r = client.post("/workflows/WF-1/remove", data={"mode": "outputs", "confirm": "yes", "typed_run_id": "WF-x"})
    assert r.status_code == 302 and out.exists()            # wrong run ID typed -> nothing deleted
    r = client.post("/workflows/WF-1/remove", data={"mode": "outputs", "confirm": "yes", "typed_run_id": "WF-1"})
    assert r.status_code == 302 and not out.exists()
    assert client.get("/workflows/WF-1").status_code == 404
