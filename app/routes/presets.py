"""Parameter presets: list, view, edit, duplicate, delete, import a filter profile from a VCF,
save the Start form or a past run as a preset."""
from __future__ import annotations

import json
import re

from flask import (Blueprint, Response, current_app, flash, redirect, render_template, request,
                   url_for)
from sqlalchemy.exc import IntegrityError

from ..engine.params import ParamError, parse_form
from ..engine.registry import get_workflow
from ..extensions import db
from ..models import ParamPreset, WorkflowRun
from ..services import preset_service as ps
from ..utils.validators import ValidationError, resolve_within_roots

bp = Blueprint("presets", __name__, url_prefix="/presets")
NAME_RE = re.compile(r"^[A-Za-z0-9 _().,+-]{1,80}$")


def _name(raw: str) -> str:
    name = (raw or "").strip()
    if not NAME_RE.match(name):
        raise ValidationError("Preset name: 1-80 characters, letters, digits, spaces and _ ( ) . , + - only.")
    return name


def _save(preset: ParamPreset, next_url: str):
    try:
        db.session.add(preset)
        db.session.commit()
    except IntegrityError:
        db.session.rollback()
        flash(f"A preset named '{preset.name}' already exists.", "danger")
        return redirect(next_url)
    flash(f"Preset '{preset.name}' saved.", "success")
    return redirect(url_for("presets.detail", preset_id=preset.id))


@bp.route("/")
def index():
    ps.ensure_builtin_presets()
    presets = ParamPreset.query.order_by(ParamPreset.locked.desc(), ParamPreset.name).all()
    return render_template("presets/index.html", presets=presets, origins=ParamPreset.ORIGINS,
                           data_roots=current_app.config["DATA_ROOTS"])


@bp.route("/<int:preset_id>")
def detail(preset_id):
    p = db.get_or_404(ParamPreset, preset_id)
    flat = sorted(ps.flatten_for_form(p.params or {}).items())
    runs = [r for r in WorkflowRun.query.order_by(WorkflowRun.created_at.desc()).limit(500)
            if ((r.params or {}).get("preset") or {}).get("id") == p.id]
    return render_template("presets/detail.html", p=p, flat=flat, origins=ParamPreset.ORIGINS, runs=runs,
                           params_json=json.dumps(p.params, indent=2, sort_keys=True))


@bp.route("/<int:preset_id>.json")
def download(preset_id):
    p = db.get_or_404(ParamPreset, preset_id)
    body = json.dumps({"name": p.name, "description": p.description, "origin": p.origin, "source": p.source,
                       "params": p.params, "notes": p.notes}, indent=2)
    return Response(body, mimetype="application/json",
                    headers={"Content-Disposition": f"attachment; filename=preset_{p.id}.json"})


@bp.route("/<int:preset_id>/edit", methods=["POST"])
def edit(preset_id):
    p = db.get_or_404(ParamPreset, preset_id)
    back = url_for("presets.detail", preset_id=p.id)
    if p.locked:
        flash("Built-in presets are locked. Duplicate it to make an editable copy.", "warning")
        return redirect(back)
    try:
        name = _name(request.form.get("name"))
        params = json.loads(request.form.get("params_json") or "{}")
    except ValidationError as exc:
        flash(str(exc), "danger")
        return redirect(back)
    except json.JSONDecodeError as exc:
        flash(f"Parameters are not valid JSON: {exc}", "danger")
        return redirect(back)
    clean, problems = ps.validate_preset_params(params)
    if problems:
        flash("Not saved — " + "; ".join(problems[:8]), "danger")
        return redirect(back)
    p.name, p.params = name, clean
    p.description = (request.form.get("description") or "").strip()[:2000]
    return _save(p, back)


@bp.route("/<int:preset_id>/duplicate", methods=["POST"])
def duplicate(preset_id):
    p = db.get_or_404(ParamPreset, preset_id)
    base = f"Copy of {p.name}"[:74]
    name, i = base, 2
    while ParamPreset.query.filter_by(name=name).first():
        name, i = f"{base} ({i})", i + 1
    return _save(ParamPreset(name=name, description=p.description, origin="user", source=f"copy of preset {p.id}",
                             params=dict(p.params or {}), notes=p.notes, profile=p.profile),
                 url_for("presets.detail", preset_id=p.id))


@bp.route("/<int:preset_id>/delete", methods=["POST"])
def delete(preset_id):
    p = db.get_or_404(ParamPreset, preset_id)
    if p.locked:
        flash("Built-in presets cannot be deleted.", "warning")
        return redirect(url_for("presets.detail", preset_id=p.id))
    db.session.delete(p)
    db.session.commit()
    flash(f"Preset '{p.name}' deleted. Runs that used it keep their full parameters.", "info")
    return redirect(url_for("presets.index"))


@bp.route("/import", methods=["POST"])
def import_vcf():
    """Grab the filtering profile of a VCF and make it a preset."""
    from ..services.filter_profile import extract_profile, profile_to_params
    back = url_for("presets.index")
    try:
        path = resolve_within_roots(request.form.get("path", "").strip(), current_app.config["DATA_ROOTS"])
        if not path.is_file() or not (path.name.endswith(".vcf") or path.name.endswith(".vcf.gz")):
            raise ValidationError(f"Not a VCF file: {path}")
        name = _name(request.form.get("name") or f"{path.name.split('.')[0]}-like")
        profile = extract_profile(path)
    except ValidationError as exc:
        flash(str(exc), "danger")
        return redirect(back)
    except (ValueError, OSError) as exc:
        flash(f"Could not read the VCF: {exc}", "danger")
        return redirect(back)
    params, notes = profile_to_params(profile)
    clean, problems = ps.validate_preset_params(params)
    if problems:
        flash("Mapping produced invalid values — " + "; ".join(problems), "danger")
        return redirect(back)
    desc = (request.form.get("description") or "").strip() or (
        f"Filter thresholds derived from {path.name} ({profile['caller']}, {profile['records']:,} records). "
        "Only filtering and GATK's emission threshold are set; everything else keeps the form's values.")
    return _save(ParamPreset(name=name, description=desc[:2000], origin="imported", source=str(path),
                             params=clean, notes=notes, profile=profile), back)


@bp.route("/from-form", methods=["POST"])
def from_form():
    """'Save as preset' button on the Start-workflow form: the same parsing as Start."""
    from ..pipelines.germline import form_specs
    back = request.form.get("next") or url_for("presets.index")
    wf = request.form.get("workflow_type", "")
    try:
        definition = get_workflow(wf)
        name = _name(request.form.get("preset_name"))
        params = parse_form(form_specs(definition), request.form, f"p.{wf}.")
        if definition.parse_extra:
            params.update(definition.parse_extra(request.form, wf))
        params.pop("baseline_file_id", None)        # file choices are not parameters
    except (KeyError, ParamError, ValidationError) as exc:
        flash(f"Preset not saved — {exc}", "danger")
        return redirect(back)
    clean, problems = ps.validate_preset_params(params)
    if problems:
        flash("Preset not saved — " + "; ".join(problems[:8]), "danger")
        return redirect(back)
    return _save(ParamPreset(name=name, origin="user", source=f"Start form ({definition.name})", params=clean,
                             description=(request.form.get("preset_description") or "").strip()[:2000]), back)


@bp.route("/from-run/<run_id>", methods=["POST"])
def from_run(run_id):
    run = WorkflowRun.query.filter_by(run_id=run_id).first_or_404()
    back = url_for("workflows.detail", run_id=run_id)
    try:
        name = _name(request.form.get("preset_name") or f"Params of {run_id}")
    except ValidationError as exc:
        flash(str(exc), "danger")
        return redirect(back)
    params = {k: v for k, v in (run.params or {}).items() if k not in ("preset", "threads", "baseline_file_id")}
    top = ps._top_level_specs()
    dropped = [k for k in params if k != "caller_params" and k not in top]
    for k in dropped:
        params.pop(k)
    clean, problems = ps.validate_preset_params(params)
    if problems:
        flash("Preset not saved — " + "; ".join(problems[:8]), "danger")
        return redirect(back)
    notes = [{"param": k, "value": None, "status": "not copied",
              "reason": "Internal run setting, not a form parameter."} for k in dropped]
    desc = f"Parameters of run {run_id} ({run.workflow_type}, workflow v{run.workflow_version})."
    return _save(ParamPreset(name=name, origin="run", source=run_id, params=clean, notes=notes or None,
                             description=desc), back)


def presets_for_form() -> list[dict]:
    """Presets as the Start form's JS needs them.
    The 'raw command' switch is a mode, not a threshold: a preset that does not set it means
    'off', so switching from the raw preset to another one never leaves raw mode switched on."""
    from ..pipelines.variant_calling.callers import all_callers
    ps.ensure_builtin_presets()
    out = []
    for p in ParamPreset.query.order_by(ParamPreset.locked.desc(), ParamPreset.name):
        values = ps.flatten_for_form(p.params or {})
        for c in all_callers("germline"):
            values.setdefault(f"c:{c.key}.tool_defaults", False)
        out.append({"id": p.id, "name": p.name, "origin": p.origin, "description": p.description or "",
                    "values": values})
    return out

