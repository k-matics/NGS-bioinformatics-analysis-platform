from pathlib import Path

from flask import Blueprint, current_app, flash, redirect, render_template, request, url_for
from sqlalchemy.exc import IntegrityError

from ..extensions import db
from ..models import ReferenceGenome
from ..services.reference_service import apply_report, inspect_reference
from ..services.settings_service import all_settings, set_setting
from ..services.tool_service import find_gatk_java_home, tool_inventory, tool_version
from ..utils.validators import ValidationError, resolve_within_roots

bp = Blueprint("settings", __name__, url_prefix="/settings")


@bp.route("/")
def index():
    values = all_settings()
    tool_version.cache_clear()
    return render_template(
        "settings/index.html", settings=values, references=ReferenceGenome.query.all(),
        tools=tool_inventory(values.get("gatk_path")), config=current_app.config,
        gatk_java=find_gatk_java_home(values.get("gatk_java_home") or None),
    )


@bp.route("/runtime", methods=["POST"])
def save_runtime():
    try:
        threads = int(request.form["threads_per_job"])
        min_free = float(request.form["min_free_disk_gb"])
        records = int(request.form["fastq_check_records"])
        if not 1 <= threads <= (current_app.config.get("MAX_THREADS") or 64):
            raise ValueError("threads out of range")
        if min_free < 0 or not 1000 <= records <= 10_000_000:
            raise ValueError("value out of range")
        gatk = request.form.get("gatk_path", "").strip()
        set_setting("threads_per_job", threads)
        set_setting("min_free_disk_gb", min_free)
        set_setting("fastq_check_records", records)
        set_setting("retain_intermediates", request.form.get("retain_intermediates") == "on")
        if gatk:
            set_setting("gatk_path", gatk)
        java_home = request.form.get("gatk_java_home", "").strip()
        if java_home and not Path(java_home, "bin", "java").is_file():
            raise ValueError(f"no bin/java under {java_home}")
        set_setting("gatk_java_home", java_home)
        find_gatk_java_home.cache_clear()
        import re
        image = request.form.get("deepvariant_image", "").strip()
        if image:
            if not re.fullmatch(r"[A-Za-z0-9._/-]+(:[A-Za-z0-9._-]+)?(@sha256:[0-9a-f]{64})?", image):
                raise ValueError("invalid Docker image name")
            set_setting("deepvariant_image", image)
        rtg_mem = int(request.form.get("rtg_mem_gb") or 8)
        if not 2 <= rtg_mem <= 48:
            raise ValueError("RTG memory must be 2-48 GB")
        set_setting("rtg_mem_gb", rtg_mem)
        flash("Settings saved. They apply to newly submitted runs.", "success")
    except (KeyError, ValueError) as exc:
        flash(f"Invalid settings: {exc}", "danger")
    return redirect(url_for("settings.index"))


@bp.route("/references", methods=["POST"])
def add_reference():
    try:
        name = request.form.get("name", "").strip()
        if not name or len(name) > 64:
            raise ValidationError("Reference name required (max 64 chars).")
        roots = current_app.config["DATA_ROOTS"]
        fasta = resolve_within_roots(request.form.get("fasta_path", "").strip(), roots)
        ref = ReferenceGenome(name=name, fasta_path=str(fasta), notes=request.form.get("notes", "")[:1000])
        apply_report(ref, inspect_reference(fasta))
        db.session.add(ref)
        db.session.commit()
        flash(f"Reference {name} registered: build {ref.build}, status {ref.status}.", "success")
    except ValidationError as exc:
        flash(str(exc), "danger")
    except IntegrityError:
        db.session.rollback()
        flash("A reference with this name already exists.", "danger")
    return redirect(url_for("settings.index"))


@bp.route("/references/<int:ref_id>/check", methods=["POST"])
def check_reference(ref_id):
    ref = db.get_or_404(ReferenceGenome, ref_id)
    apply_report(ref, inspect_reference(Path(ref.fasta_path)))
    db.session.commit()
    flash(f"Reference {ref.name} re-checked: {ref.status}.", "info")
    return redirect(url_for("settings.index"))


@bp.route("/references/<int:ref_id>/default", methods=["POST"])
def make_default(ref_id):
    ref = db.get_or_404(ReferenceGenome, ref_id)
    for other in ReferenceGenome.query.all():
        other.is_default = other.id == ref.id
    db.session.commit()
    flash(f"{ref.name} is now the default reference.", "info")
    return redirect(url_for("settings.index"))
