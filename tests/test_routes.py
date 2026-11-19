from app.models import Sample, WorkflowRun
from conftest import GRCH38, write_bam


def test_pages_load(client):
    for url in ["/", "/samples/", "/samples/new", "/workflows/", "/settings/", "/module/fastq",
                "/api/system", "/api/dashboard"]:
        assert client.get(url).status_code == 200, url


def test_create_sample_register_file_and_run(client, data_dir):
    r = client.post("/samples/new", data={"sample_id": "SAMPLE_A", "sequencing_type": "WES",
                                          "library_layout": "paired"}, follow_redirects=True)
    assert b"SAMPLE_A" in r.data
    bam = write_bam(data_dir / "sample_a.bam", GRCH38)
    client.post("/samples/SAMPLE_A/files", data={"path": str(bam), "role": "", "origin": "internal"})
    sample = Sample.query.filter_by(sample_id="SAMPLE_A").one()
    assert sample.files[0].role == "bam"
    ref_id = 1
    r = client.post("/samples/SAMPLE_A/run", data={"workflow_type": "intake_validation", "reference_id": ref_id,
                                                  "file_ids": [sample.files[0].id]})
    assert r.status_code == 302
    run = WorkflowRun.query.one()
    assert client.get(f"/workflows/{run.run_id}").status_code == 200
    state = client.get(f"/api/runs/{run.run_id}/state").get_json()
    assert state["status"] == "completed" and state["progress"] == 100
    # SSE stream terminates for finished runs and contains state + end events
    body = client.get(f"/api/runs/{run.run_id}/events").get_data(as_text=True)
    assert "event: state" in body and "event: end" in body and "event: log" in body
    # report download within run dir works
    assert client.get(f"/workflows/{run.run_id}/files/reports/intake_report.html").status_code == 200


def test_bad_sample_id_rejected(client):
    r = client.post("/samples/new", data={"sample_id": "../../x", "sequencing_type": "WES",
                                          "library_layout": "paired"}, follow_redirects=True)
    assert b"Sample ID must be" in r.data
    assert Sample.query.count() == 0


def test_file_outside_roots_rejected(client):
    client.post("/samples/new", data={"sample_id": "S2", "sequencing_type": "WES", "library_layout": "paired"})
    r = client.post("/samples/S2/files", data={"path": "/etc/passwd", "origin": "internal"}, follow_redirects=True)
    assert b"outside the allowed" in r.data


def test_download_path_traversal_blocked(client, data_dir):
    client.post("/samples/new", data={"sample_id": "S3", "sequencing_type": "WES", "library_layout": "paired"})
    bam = write_bam(data_dir / "s3.bam", GRCH38)
    client.post("/samples/S3/files", data={"path": str(bam), "origin": "internal"})
    s = Sample.query.filter_by(sample_id="S3").one()
    client.post("/samples/S3/run", data={"workflow_type": "intake_validation", "reference_id": 1,
                                         "file_ids": [s.files[0].id]})
    run = WorkflowRun.query.one()
    assert client.get(f"/workflows/{run.run_id}/files/../../../../etc/passwd").status_code in (403, 404)


def test_browse_outside_roots(client):
    assert client.get("/api/browse?path=/etc").status_code == 400
