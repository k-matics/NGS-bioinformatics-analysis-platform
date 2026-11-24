import os
from pathlib import Path

from flask import Blueprint, abort, flash, redirect, render_template, request, send_file, url_for

from ..engine.registry import get_workflow
from ..models import WorkflowRun
from ..models import status as S
from ..services.file_service import safe_child
from ..services.limitations import run_limitations
from ..services.report_service import run_reports
from ..services.run_delete import DeleteError, delete_run, deletion_plan
from ..services.workflow_service import (
    WorkflowError,
    request_cancel,
    request_pause,
    request_unpause,
    resume_run,
)

bp = Blueprint("workflows", __name__, url_prefix="/workflows")


def _get_run(run_id: str) -> WorkflowRun:
    run = WorkflowRun.query.filter_by(run_id=run_id).first()
    if run is None:
        abort(404)
    return run


@bp.route("/")
def history():
    status = request.args.get("status", "")
    query = WorkflowRun.query
    if status in S.RUN_STATUSES:
        query = query.filter_by(status=status)
    runs = query.order_by(WorkflowRun.created_at.desc()).limit(500).all()
    return render_template("workflows/list.html", runs=runs, status=status, statuses=S.RUN_STATUSES)


@bp.route("/<run_id>")
def detail(run_id):
    run = _get_run(run_id)
    try:
        definition = get_workflow(run.workflow_type)
    except KeyError:
        definition = None
    return render_template("workflows/detail.html", run=run, definition=definition,
        reports=run_reports(run), limitations=run_limitations(run.reference, run.params, run.sample))


@bp.route("/<run_id>/remove", methods=["GET", "POST"])
def remove(run_id):
    """Confirmation page (GET) and removal (POST) of a finished run."""
    run = _get_run(run_id)
    if request.method == "GET":
        return render_template("workflows/remove.html", run=run, plan=deletion_plan(run))
    mode = request.form.get("mode", "")
    if request.form.get("confirm") != "yes":
        flash("Tick the confirmation box to remove the run.", "warning")
        return redirect(url_for("workflows.remove", run_id=run_id))
    if mode == "outputs" and request.form.get("typed_run_id", "").strip() != run.run_id:
        flash("To delete the output files, type the run ID exactly.", "warning")
        return redirect(url_for("workflows.remove", run_id=run_id))
    try:
        entry = delete_run(run, mode)
    except DeleteError as exc:
        flash(str(exc), "danger")
        if WorkflowRun.query.filter_by(run_id=run_id).first() is None:
            return redirect(url_for("workflows.history"))
        return redirect(url_for("workflows.remove", run_id=run_id))
    _VC_CACHE.clear()
    if mode == "outputs":
        flash(f"Removed {run_id}: run record and output folder ({entry['files_removed']} files, "
              f"{entry['bytes_removed'] / 1e9:.2f} GB, {len(entry['unregistered_files'])} registered file(s) "
              "un-registered). Logged in logs/deletions.jsonl.", "success")
    else:
        flash(f"Removed the record of {run_id}. Its files were kept on disk in {entry['output_dir']}. "
              "Logged in logs/deletions.jsonl.", "success")
    return redirect(url_for("workflows.history"))


@bp.route("/<run_id>/<action>", methods=["POST"])
def control(run_id, action):
    run = _get_run(run_id)
    actions = {"cancel": request_cancel, "pause": request_pause,
               "unpause": request_unpause, "resume": resume_run}
    if action not in actions:
        abort(404)
    try:
        actions[action](run)
        flash(f"{action.capitalize()} requested for {run_id}.", "info")
    except WorkflowError as exc:
        flash(str(exc), "danger")
    return redirect(url_for("workflows.detail", run_id=run_id))


@bp.route("/<run_id>/files/")
@bp.route("/<run_id>/files/<path:relpath>")
def files(run_id, relpath=""):
    """Browse / download files inside the run's output directory only."""
    run = _get_run(run_id)
    base = Path(run.output_dir)
    try:
        target = safe_child(base, relpath)
    except PermissionError:
        abort(403)
    if not target.exists():
        abort(404)
    if target.is_file():
        inline = target.suffix in (".html", ".txt", ".log", ".json") and request.args.get("download") != "1"
        return send_file(target, as_attachment=not inline,
                         mimetype="text/plain" if target.suffix in (".log", ".sha256") else None)
    entries = []
    for entry in sorted(os.scandir(target), key=lambda e: (not e.is_dir(), e.name)):
        entries.append({"name": entry.name, "is_dir": entry.is_dir(),
                        "size": None if entry.is_dir() else entry.stat().st_size,
                        "rel": str(Path(relpath) / entry.name) if relpath else entry.name})
    return render_template("workflows/files.html", run=run, relpath=relpath, entries=entries, base=base)


# ------------------------------------------------------------------ variant-by-variant viewer
_VC_CACHE: dict = {}
VC_FILTER_FIELDS = ("min_dp", "min_gq", "min_vaf", "min_qual", "b_min_qual", "b_min_dp", "b_min_mq", "b_max_mq")


def _vc_rows(run):
    from ..pipelines.comparison.variant_check import read_rows
    step = next((s for s in run.steps if s.key == "variant_check" and s.status == S.COMPLETED), None)
    if step is None or not (step.outputs or {}).get("variant_tsv"):
        abort(404)
    path = Path(step.outputs["variant_tsv"])
    key = (str(path), path.stat().st_mtime_ns)
    if key not in _VC_CACHE:
        _VC_CACHE.clear()                       # keep only one run's table in memory
        _VC_CACHE[key] = read_rows(path)
    return step, _VC_CACHE[key]


def _float_arg(name):
    raw = request.args.get(name, "").strip()
    if not raw:
        return None
    try:
        return float(raw)
    except ValueError:
        abort(400, f"{name} must be a number")


def _vc_search():
    """The viewer's variant search box: (locus or None, window bp, error message or None)."""
    from ..pipelines.comparison.variant_check import parse_locus
    q = request.args.get("q", "").strip()
    window = max(0, min(request.args.get("near", 0, type=int) or 0, 100_000))
    if not q:
        return None, window, None
    try:
        return parse_locus(q), window, None
    except ValueError as exc:
        return None, window, str(exc)


def _vc_filtered(rows):
    """Apply the viewer's filters. Returns (rows in view, thresholds, what-if result)."""
    from ..pipelines.comparison.variant_check import match_locus, what_if
    locus, window, _ = _vc_search()
    if locus:
        rows = match_locus(rows, locus, window)
    statuses = set(request.args.getlist("status"))
    types = set(request.args.getlist("type"))
    chrom = request.args.get("chrom", "").strip()
    th = {k: _float_arg(k) for k in VC_FILTER_FIELDS}
    ignore = request.args.get("ignore_filters") == "1"
    view = request.args.get("view", "all")

    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    base = []
    for r in rows:
        if types and r["TYPE"] not in types:
            continue
        if chrom and r["CHROM"] != chrom:
            continue
        if th["b_min_qual"] is not None and (num(r["B_QUAL"]) or 0) < th["b_min_qual"]:
            continue
        if th["b_min_dp"] is not None and (num(r["B_DP"]) or 0) < th["b_min_dp"]:
            continue
        # baseline mapping quality (INFO/MQ); rows without MQ (older tables, callers without MQ) are
        # excluded by an MQ filter rather than silently counted as passing it
        if th["b_min_mq"] is not None and (num(r.get("B_MQ")) is None or num(r.get("B_MQ")) < th["b_min_mq"]):
            continue
        if th["b_max_mq"] is not None and (num(r.get("B_MQ")) is None or num(r.get("B_MQ")) > th["b_max_mq"]):
            continue
        base.append(r)
    kw = dict(min_dp=th["min_dp"], min_gq=th["min_gq"], min_vaf=th["min_vaf"], min_qual=th["min_qual"],
              ignore_our_filters=ignore)
    result = what_if(base, **kw)
    out = []
    for r in base:
        if statuses and r["STATUS"] not in statuses:
            continue
        if view != "all":
            found = what_if([r], **kw)["found"] == 1
            if (view == "found") != found:
                continue
        out.append(r)
    return out, th, ignore, view, result


@bp.route("/<run_id>/variant-check")
def variant_check(run_id):
    from ..pipelines.comparison.variant_check import COLUMNS, STATUS_LABELS
    run = _get_run(run_id)
    step, rows = _vc_rows(run)
    view_rows, th, ignore, view, result = _vc_filtered(rows)
    page = max(1, request.args.get("page", 1, type=int))
    per = 100
    pages = max(1, (len(view_rows) + per - 1) // per)
    page = min(page, pages)
    chroms = sorted({r["CHROM"] for r in rows}, key=lambda c: (len(c), c))
    args = request.args.to_dict(flat=False)
    args.pop("page", None)
    locus, window, search_error = _vc_search()
    ours_only = None
    if locus:                       # also look in 'our PASS variants absent from the baseline'
        from ..pipelines.comparison.variant_check import match_locus, read_rows
        qtsv = (step.outputs or {}).get("query_only_tsv")
        if qtsv and Path(qtsv).is_file():
            ours_only = match_locus(read_rows(Path(qtsv)), locus, window)
    return render_template("workflows/variant_check.html", run=run, s=step.summary, rows=view_rows[(page - 1) * per:page * per],
                           q=request.args.get("q", ""), near=window, locus=locus, search_error=search_error,
                           ours_only=ours_only,
                           n=len(view_rows), page=page, pages=pages, columns=COLUMNS, labels=STATUS_LABELS,
                           th=th, ignore=ignore, view=view, result=result, chroms=chroms,
                           sel_status=set(request.args.getlist("status")), sel_type=set(request.args.getlist("type")),
                           chrom=request.args.get("chrom", ""), args=args)


@bp.route("/<run_id>/variant-check.csv")
def variant_check_csv(run_id):
    import csv
    import io

    from flask import Response

    from ..pipelines.comparison.variant_check import COLUMNS
    run = _get_run(run_id)
    _, rows = _vc_rows(run)
    view_rows, th, ignore, view, result = _vc_filtered(rows)
    buf = io.StringIO()
    buf.write(f"# run {run_id}; filters: {request.query_string.decode()[:500]}; "
              f"found under thresholds: {result['found']}/{result['total']}\n")
    w = csv.DictWriter(buf, fieldnames=COLUMNS)
    w.writeheader()
    w.writerows(view_rows)
    return Response(buf.getvalue(), mimetype="text/csv",
                    headers={"Content-Disposition": f"attachment; filename={run_id}_variant_check.csv"})


@bp.route("/<run_id>/variant-check/site")
def variant_check_site(run_id):
    """'Check reads': the reads of our BAM at one variant, and why a caller may have missed it."""
    import re

    from ..models import SampleFile
    from ..services.site_inspect import bam_for_vcf, inspect_site
    run = _get_run(run_id)
    step, _ = _vc_rows(run)
    chrom = request.args.get("chrom", "")
    pos = request.args.get("pos", type=int)
    ref, alt = request.args.get("ref", "").upper(), request.args.get("alt", "").upper()
    if not re.fullmatch(r"[A-Za-z0-9_.]{1,64}", chrom) or not pos or pos < 1 \
            or not re.fullmatch(r"[ACGTN]{1,500}", ref) or not re.fullmatch(r"[ACGTN]{1,500}", alt):
        abort(400, "chrom, pos, ref and alt are required (ref/alt as A/C/G/T/N).")
    bams = [f for f in run.sample.files if f.role == "bam"]
    default = bam_for_vcf(run.sample, (step.summary or {}).get("query", {}).get("source", ""))
    bam = default
    if request.args.get("bam_id", type=int):
        bam = next((f for f in bams if f.id == request.args.get("bam_id", type=int)), None)
        if bam is None:
            abort(404)
    min_mq = max(0, min(request.args.get("min_mq", 20, type=int), 60))
    min_bq = max(0, min(request.args.get("min_bq", 10, type=int), 60))
    result = error = None
    if bam is not None:
        try:
            result = inspect_site(Path(bam.path), chrom, pos, ref, alt, min_mq=min_mq, min_bq=min_bq)
        except (ValueError, OSError) as exc:
            error = str(exc)
    row = next((r for r in _vc_rows(run)[1] if r["CHROM"] == chrom and r["POS"] == str(pos)
                and r["REF"] == ref and r["ALT"] == alt), None)
    return render_template("workflows/variant_check_site.html", run=run, s=step.summary, chrom=chrom, pos=pos,
                           ref=ref, alt=alt, bams=bams, bam=bam, default_bam=default, result=result, error=error,
                           row=row, min_mq=min_mq, min_bq=min_bq)
