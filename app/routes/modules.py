"""Sections that are planned for later phases. They state honestly what is
coming, instead of showing non-functional controls."""
from flask import Blueprint, abort, render_template

bp = Blueprint("modules", __name__)

# Modules that are working (Phase 3): key -> (title, workflow modules, intro)
ACTIVE = {
    "fastq": ("FASTQ Analysis", ("fastq", "bam"),
              "Read quality control, adapter/quality trimming and alignment to the reference."),
    "bam": ("BAM Analysis", ("bam",),
            "Alignment, duplicate marking, indexing and BAM QC / target coverage."),
    "variants": ("Variant Calling", ("vcf",),
                 "Germline short-variant calling with one or more callers, then normalisation."),
}

ACTIVE.update({
    "caller-comparison": ("Caller Comparison", ("comparison",),
                          "Compare variant callers on the same BAM. Select several callers in a germline workflow "
                          "(the run then contains a caller-comparison step), or compare any registered VCFs of a "
                          "sample with the VCF-vs-VCF workflow. Optional truth set (e.g. GIAB) for precision/recall."),
    "external": ("External Data Comparison", ("external", "comparison"),
                 "Compare our results with files from external laboratories or providers: BAM vs BAM QC, VCF vs VCF "
                 "concordance and read-level evidence of a VCF in a BAM. Genome builds are checked and never mixed."),
})
COMPARISON_STEPS = {"caller-comparison": ("compare_callers", "vcf_compare"),
                    "external": ("bam_compare", "vcf_compare", "read_evidence")}

PLANNED: dict = {}


@bp.route("/module/<key>")
def planned(key):
    if key == "reports":
        return reports_page()
    if key in ACTIVE:
        return _active(key)
    if key not in PLANNED:
        abort(404)
    title, phase, items = PLANNED[key]
    return render_template("planned.html", title=title, phase=phase, items=items)


def _active(key):
    from ..engine.registry import all_workflows
    from ..models import WorkflowRun
    from ..pipelines.variant_calling.callers import all_callers

    title, modules, intro = ACTIVE[key]
    workflows = [w for w in all_workflows() if w.module in modules]
    runs = (WorkflowRun.query.filter(WorkflowRun.workflow_type.in_([w.key for w in workflows]))
            .order_by(WorkflowRun.id.desc()).limit(25).all())
    callers = [(c, *c.available()) for c in all_callers()] if key in ("variants", "caller-comparison") else []
    comparisons = []
    if key in COMPARISON_STEPS:
        from ..models import StepRun
        q = (StepRun.query.filter(StepRun.key.in_(COMPARISON_STEPS[key]), StepRun.status == "completed")
             .order_by(StepRun.id.desc()).limit(30).all())
        for s in q:
            html = (s.outputs or {}).get("report_html")
            if html:
                import os
                comparisons.append({"run": s.run, "step": s, "summary": s.summary or {},
                                    "rel": os.path.relpath(html, s.run.output_dir)})
    return render_template("module.html", title=title, intro=intro, workflows=workflows, runs=runs,
                           callers=callers, comparisons=comparisons)


def reports_page():
    from ..models import WorkflowRun
    from ..services.report_service import run_reports
    runs = (WorkflowRun.query.filter(WorkflowRun.status == "completed").order_by(WorkflowRun.id.desc())
            .limit(200).all())
    items = [(r, run_reports(r)) for r in runs]
    return render_template("reports/index.html", items=items)


@bp.route("/reports/<kind>.csv")
def reports_csv(kind):
    import csv
    import io

    from flask import Response

    from ..models import WorkflowRun
    from ..services import report_service as rs
    fn = {"qc": rs.qc_rows, "variants": rs.variant_rows, "comparisons": rs.comparison_rows}.get(kind)
    if fn is None:
        abort(404)
    runs = WorkflowRun.query.filter(WorkflowRun.status == "completed").order_by(WorkflowRun.id).all()
    rows = fn(runs)
    buf = io.StringIO()
    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    w = csv.DictWriter(buf, fieldnames=fields or ["no_rows"])
    w.writeheader()
    w.writerows(rows)
    return Response(buf.getvalue(), mimetype="text/csv",
                    headers={"Content-Disposition": f"attachment; filename=ngs_{kind}.csv"})
