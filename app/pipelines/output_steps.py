"""End-of-run steps: register pipeline outputs as sample files, optional cleanup of intermediates."""
from __future__ import annotations

import shutil
from pathlib import Path

from ..engine.step import Step, StepContext, StepResult
from ..extensions import db
from ..models.base import utcnow
from ..services.file_service import file_fingerprint, sha256_file
from .common import prior_output, write_json


class RegisterOutputsStep(Step):
    """Makes the final BAM and filtered VCFs selectable as inputs of comparison workflows
    (registered with origin 'pipeline', a SHA-256 checksum and the run ID in the label)."""
    key = "register_outputs"
    name = "Register outputs"
    description = "Registers the analysis BAM and filtered VCFs as sample files (with checksums) for comparisons."

    def execute(self, ctx: StepContext, runner) -> StepResult:
        from ..models import SampleFile
        from .variant_calling.callers import get_caller
        items = []
        pr = (ctx.params.get("preset") or {})
        preset = f" · {pr['name']}{' (modified)' if pr.get('changed_after_applying') else ''}" if pr.get("name") else ""
        src = next((f for f in ctx.files("bam") if f.origin != "pipeline"), None)
        if src is not None:                      # calls made on a provider's / external BAM
            preset += f" · on {src.label or Path(src.path).name}"
        bam = prior_output(ctx, "index_bam", "bam", required=False)
        if bam:
            items.append(("bam", bam, f"{ctx.run.run_id} · pipeline BAM{preset}"))
        for k in ctx.params.get("callers", []):
            vcf = prior_output(ctx, f"filter_{k}", "filtered_vcf", required=False)
            if vcf:
                items.append(("vcf", vcf, f"{ctx.run.run_id} · {get_caller(k).name} (filtered){preset}"))
        registered = []
        for role, path, label in items:
            p = Path(path)
            if SampleFile.query.filter_by(sample_pk=ctx.sample.id, path=str(p)).first():
                continue
            ctx.log(f"Checksumming {p.name} ...")
            f = SampleFile(sample_pk=ctx.sample.id, role=role, origin="pipeline", label=label[:128],
                           path=str(p), sha256=sha256_file(p), checksum_at=utcnow(),
                           detected_build=ctx.reference.build, validation_status="ok",
                           validation_message=f"Produced by {ctx.run.run_id} ({ctx.run.workflow_type} "
                                              f"v{ctx.run.workflow_version}) on {ctx.reference.name}.",
                           **file_fingerprint(p))
            db.session.add(f)
            registered.append({"role": role, "path": str(p), "label": label})
        db.session.commit()
        return StepResult(summary={"registered": registered})

    def validate(self, ctx, result):
        return None


INTERMEDIATE_KEYS = [("fastp", "trimmed_r1"), ("fastp", "trimmed_r2")]


class CleanupStep(Step):
    """Opt-in removal of large re-creatable intermediates after a SUCCESSFUL run.
    Removed reads are never deleted (failed_reads.fastq.gz is kept). Every removed file is
    listed (path, size) in reports/cleanup_manifest.json."""
    key = "cleanup"
    name = "Clean up intermediate files"
    description = "Only when 'Keep intermediate files' is off: removes trimmed FASTQs and caller temp folders."

    def applicable(self, ctx):
        keep = ctx.params.get("retain_intermediates", True)
        return (not keep, "Intermediate files are kept (default)." if keep else "")

    def execute(self, ctx: StepContext, runner) -> StepResult:
        removed = []
        for step, name in INTERMEDIATE_KEYS:
            p = prior_output(ctx, step, name, required=False)
            if p and Path(p).exists():
                removed.append({"path": p, "bytes": Path(p).stat().st_size})
                Path(p).unlink()
        for k in ctx.params.get("callers", []):
            tmp = ctx.output_dir / "vcf" / k / "tmp"
            if tmp.exists():
                size = sum(f.stat().st_size for f in tmp.rglob("*") if f.is_file())
                shutil.rmtree(tmp)
                removed.append({"path": str(tmp), "bytes": size})
        for r in removed:
            ctx.log(f"Removed intermediate: {r['path']} ({r['bytes'] / 1e9:.2f} GB)")
        manifest = write_json(ctx.step_dir("reports") / "cleanup_manifest.json",
                              {"removed": removed, "at": utcnow().isoformat() + "Z"})
        return StepResult(outputs={"manifest": manifest},
                          summary={"removed_files": len(removed),
                                   "freed_gb": round(sum(r["bytes"] for r in removed) / 1e9, 2)})
