"""Collects report files and metrics across runs (Reports module, CSV exports)."""
from __future__ import annotations

import os
from pathlib import Path

REPORT_STEPS = ("compare_callers", "vcf_compare", "bam_compare", "read_evidence", "report")


def run_reports(run) -> list[dict]:
    """HTML reports produced by a run, as paths relative to its output folder."""
    out = []
    base = Path(run.output_dir)
    for s in run.steps:
        if s.status != "completed" or s.key not in REPORT_STEPS:
            continue
        html = (s.outputs or {}).get("report_html")
        if html and Path(html).exists():
            out.append({"name": "Run report" if s.key == "report" else s.name,
                        "rel": os.path.relpath(html, base), "key": s.key})
    mq = base / "reports" / "multiqc" / "multiqc_report.html"
    if mq.exists():
        out.append({"name": "MultiQC", "rel": os.path.relpath(mq, base), "key": "multiqc"})
    return out


def _s(run, key):
    for s in run.steps:
        if s.key == key and s.status == "completed":
            return s.summary or {}
    return {}


def qc_rows(runs) -> list[dict]:
    rows = []
    for r in runs:
        fp, md, qc = _s(r, "fastp"), _s(r, "markdup"), _s(r, "bam_qc")
        pct = qc.get("pct_target_bases") or {}
        rows.append({
            "run_id": r.run_id, "sample": r.sample.sample_id, "workflow": r.workflow_type,
            "version": r.workflow_version, "status": r.status, "reference": r.reference.name if r.reference else "",
            "created": r.created_at.isoformat() if r.created_at else "",
            "reads_raw": fp.get("reads_before"), "reads_removed_pct": fp.get("reads_removed_pct"),
            "q30_after_pct": fp.get("q30_after_pct"), "primary_reads": qc.get("primary") or md.get("primary"),
            "mapped_pct": qc.get("primary_mapped_pct"), "duplicate_pct": qc.get("duplicate_pct"),
            "properly_paired_pct": qc.get("properly_paired_pct"), "insert_size_mean": qc.get("insert_size_mean"),
            "mean_target_coverage": qc.get("mean_target_coverage"),
            **{f"pct_ge_{k}": v for k, v in pct.items()},
            "warnings": len(r.warnings or []),
        })
    return rows


def variant_rows(runs) -> list[dict]:
    rows = []
    for r in runs:
        for s in r.steps:
            if not s.key.startswith("call_") or s.status != "completed":
                continue
            k = s.key[5:]
            c, f = s.summary or {}, _s(r, f"filter_{k}")
            rows.append({"run_id": r.run_id, "sample": r.sample.sample_id, "caller": k,
                         "records": c.get("records"), "snvs": c.get("snvs"), "indels": c.get("indels"),
                         "ts_tv": c.get("ts_tv"), "het": c.get("het"), "hom_alt": c.get("hom_alt"),
                         "pass": f.get("pass"), "pass_pct": f.get("pass_pct"), "pass_snvs": f.get("pass_snvs"),
                         "pass_indels": f.get("pass_indels"), "pass_ts_tv": f.get("pass_ts_tv")})
    return rows


def comparison_rows(runs) -> list[dict]:
    rows = []
    for r in runs:
        for key in ("compare_callers", "vcf_compare"):
            s = _s(r, key)
            for p in s.get("pairwise", []):
                rows.append({"run_id": r.run_id, "sample": r.sample.sample_id, "comparison": key, **p})
    return rows
