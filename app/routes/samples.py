from pathlib import Path

from flask import Blueprint, abort, current_app, flash, redirect, render_template, request, url_for
from sqlalchemy.exc import IntegrityError

from ..engine.params import ParamError, parse_form
from ..engine.registry import all_workflows, get_workflow
from ..extensions import db
from ..models import ParamPreset, ReferenceGenome, Sample, SampleFile
from ..services.preset_service import record_on_params
from .presets import presets_for_form
from ..models.sample import FILE_ORIGINS, FILE_ROLES
from ..services.file_service import file_fingerprint
from ..services.workflow_service import WorkflowError, create_run
from ..utils.validators import (
    ValidationError,
    check_role_matches_extension,
    guess_role,
    resolve_within_roots,
    validate_sample_id,
)

bp = Blueprint("samples", __name__, url_prefix="/samples")

SEQ_TYPES = ["WES", "WGS", "Targeted panel"]


def _get_sample(sample_id: str) -> Sample:
    sample = Sample.query.filter_by(sample_id=sample_id).first()
    if sample is None:
        abort(404)
    return sample


@bp.route("/")
def list_samples():
    q = request.args.get("q", "").strip()
    query = Sample.query
    if q:
        query = query.filter(Sample.sample_id.ilike(f"%{q}%"))
    return render_template("samples/list.html", samples=query.order_by(Sample.created_at.desc()).all(), q=q)


@bp.route("/new", methods=["GET", "POST"])
def new_sample():
    if request.method == "POST":
        try:
            sample_id = validate_sample_id(request.form.get("sample_id", ""))
            seq_type = request.form.get("sequencing_type", "WES")
            layout = request.form.get("library_layout", "paired")
            if seq_type not in SEQ_TYPES or layout not in ("paired", "single"):
                raise ValidationError("Invalid sequencing type or layout.")
            sample = Sample(sample_id=sample_id, sequencing_type=seq_type, library_layout=layout,
                            description=request.form.get("description", "").strip()[:2000])
            db.session.add(sample)
            db.session.commit()
            flash(f"Sample {sample_id} created. Now register its input files.", "success")
            return redirect(url_for("samples.detail", sample_id=sample_id))
        except ValidationError as exc:
            flash(str(exc), "danger")
        except IntegrityError:
            db.session.rollback()
            flash("A sample with this ID already exists.", "danger")
    return render_template("samples/new.html", seq_types=SEQ_TYPES)


@bp.route("/<sample_id>")
def detail(sample_id):
    from ..pipelines.germline import form_specs
    from ..pipelines.variant_calling.callers import all_callers

    sample = _get_sample(sample_id)
    workflows = all_workflows()
    return render_template(
        "samples/detail.html", sample=sample, roles=FILE_ROLES, origins=FILE_ORIGINS,
        workflows=workflows, specs={w.key: form_specs(w) for w in workflows},
        callers=all_callers("germline"), references=ReferenceGenome.query.all(),
        data_roots=current_app.config["DATA_ROOTS"], presets=presets_for_form(),
    )


@bp.route("/<sample_id>/files", methods=["POST"])
def add_file(sample_id):
    sample = _get_sample(sample_id)
    try:
        path = resolve_within_roots(request.form.get("path", "").strip(), current_app.config["DATA_ROOTS"])
        if not path.is_file():
            raise ValidationError(f"File not found: {path}")
        role = request.form.get("role") or guess_role(path)
        if not role:
            raise ValidationError("Could not determine the file type; choose a role.")
        check_role_matches_extension(path, role)
        origin = request.form.get("origin", "internal")
        if origin not in FILE_ORIGINS or origin == "pipeline":
            raise ValidationError("Invalid origin.")
        if role in ("fastq_r1", "fastq_r2", "fastq_se") and sample.files_by_role(role):
            raise ValidationError(f"A {FILE_ROLES[role]} file is already registered for this sample.")
        fp = file_fingerprint(path)
        sf = SampleFile(sample=sample, role=role, origin=origin, path=str(path),
                        label=request.form.get("label", "").strip()[:128] or None, **fp)
        db.session.add(sf)
        db.session.commit()
        flash(f"Registered {path.name} as {FILE_ROLES[role]} ({FILE_ORIGINS[origin]}). "
              "Run intake validation to checksum and check it.", "success")
    except ValidationError as exc:
        flash(str(exc), "danger")
    except IntegrityError:
        db.session.rollback()
        flash("This file is already registered for the sample.", "warning")
    return redirect(url_for("samples.detail", sample_id=sample_id))


@bp.route("/<sample_id>/files/<int:file_id>/remove", methods=["POST"])
def remove_file(sample_id, file_id):
    sample = _get_sample(sample_id)
    sf = db.session.get(SampleFile, file_id)
    if sf is None or sf.sample_pk != sample.id:
        abort(404)
    used = [r.run_id for r in sample.runs if file_id in (r.input_file_ids or [])]
    if used:
        flash("File is referenced by workflow runs and cannot be unregistered (traceability). "
              f"Runs: {', '.join(used[:3])}", "warning")
    else:
        db.session.delete(sf)
        db.session.commit()
        flash(f"Unregistered {Path(sf.path).name} (the file itself was not touched).", "info")
    return redirect(url_for("samples.detail", sample_id=sample_id))


@bp.route("/<sample_id>/run", methods=["POST"])
def start_workflow(sample_id):
    sample = _get_sample(sample_id)
    workflow_type = request.form.get("workflow_type", "")
    ref = db.session.get(ReferenceGenome, request.form.get("reference_id", type=int) or -1)
    selected = set(request.form.getlist("file_ids", type=int))
    files = [f for f in sample.files if f.id in selected]
    try:
        from ..pipelines.germline import form_specs

        definition = get_workflow(workflow_type)
        params = parse_form(form_specs(definition), request.form, f"p.{workflow_type}.")
        if definition.parse_extra:
            params.update(definition.parse_extra(request.form, workflow_type))
        preset = db.session.get(ParamPreset, request.form.get("preset_id", type=int) or -1)
        if preset is not None:
            record_on_params(preset, {s.key for s in form_specs(definition)}, params)
        run = create_run(sample, workflow_type, files, ref, params=params)
        flash(f"Workflow {run.run_id} queued.", "success")
        return redirect(url_for("workflows.detail", run_id=run.run_id))
    except ParamError as exc:
        flash(f"Invalid parameter — {exc}", "danger")
        return redirect(url_for("samples.detail", sample_id=sample_id))
    except (WorkflowError, KeyError) as exc:
        flash(str(exc), "danger")
        return redirect(url_for("samples.detail", sample_id=sample_id))
