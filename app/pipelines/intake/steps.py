"""Steps of the intake validation workflow."""
from __future__ import annotations

import gzip
import json
import re
from collections import Counter
from pathlib import Path

from flask import render_template

from ...engine.step import Step, StepContext, StepFailure, StepResult
from ...extensions import db
from ...models.base import utcnow
from ...services.file_service import sha256_file
from ...services.genome_build import compare_to_reference
from ...services.hts_inspect import bam_header_summary, vcf_header_summary
from ...services.reference_service import apply_report, inspect_reference, reference_contigs


def _write_json(path: Path, data) -> str:
    path.write_text(json.dumps(data, indent=2, default=str))
    return str(path)


# ---------------------------------------------------------------------------
class ReferenceCheckStep(Step):
    key = "reference_check"
    name = "Reference genome check"
    description = "FASTA, .fai, .dict, BWA index presence/consistency; build and naming detection."

    def execute(self, ctx: StepContext, runner) -> StepResult:
        ref = ctx.reference
        ctx.log(f"Inspecting reference {ref.name}: {ref.fasta_path}")
        report = inspect_reference(Path(ref.fasta_path))
        apply_report(ref, report)
        db.session.commit()
        out = _write_json(ctx.step_dir("intake") / "reference_check.json", report)
        for e in report["errors"]:
            ctx.log(f"ERROR: {e}")
        if report["errors"]:
            raise StepFailure("Reference check failed: " + "; ".join(report["errors"]))
        ctx.shared["reference_contigs"] = reference_contigs(ref)
        return StepResult(
            outputs={"reference_check": out},
            summary={"build": report.get("build"), "naming": report.get("naming"),
                     "contigs": report.get("contig_count"), "classes": report.get("contig_classes")},
            warnings=report["warnings"],
        )


# ---------------------------------------------------------------------------
class ChecksumStep(Step):
    key = "checksums"
    name = "Input checksums (SHA-256)"
    description = "Fingerprints every input file so later results can be tied to exact inputs."

    def execute(self, ctx: StepContext, runner) -> StepResult:
        lines, summary = [], {}
        for f in ctx.input_files:
            ctx.log(f"SHA-256: {f.path}")
            digest = sha256_file(Path(f.path), should_stop=ctx.stop_requested)
            ctx.check_cancel()
            if f.sha256 and f.sha256 != digest:
                raise StepFailure(
                    f"Checksum of {f.path} differs from the value recorded earlier "
                    f"({f.sha256[:12]}… vs {digest[:12]}…). The file content changed — do not use it."
                )
            f.sha256, f.checksum_at = digest, utcnow()
            db.session.commit()
            lines.append(f"{digest}  {f.path}")
            summary[Path(f.path).name] = digest
            ctx.log(f"  {digest}")
        out = ctx.step_dir("intake") / "checksums.sha256"
        out.write_text("\n".join(lines) + "\n")
        ctx.log(f"Verify later with: sha256sum -c {out}")
        return StepResult(outputs={"checksums": str(out)}, summary=summary)


# ---------------------------------------------------------------------------
FASTQ_SEQ_RE = re.compile(rb"^[ACGTNacgtn.]*$")


def _read_name(header: bytes) -> bytes:
    name = header[1:].split()[0]
    if name.endswith((b"/1", b"/2")):
        name = name[:-2]
    return name


def scan_fastq(path: Path, max_records: int, check_cancel) -> dict:
    """Validate the 4-line FASTQ structure on the first `max_records` reads."""
    lengths, min_q, max_q, names, errors = Counter(), 255, 0, [], []
    n = 0
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rb") as fh:
        while n < max_records:
            header = fh.readline()
            if not header:
                break
            seq, plus, qual = fh.readline().rstrip(b"\r\n"), fh.readline(), fh.readline().rstrip(b"\r\n")
            n += 1
            line_no = 4 * (n - 1) + 1
            if not header.startswith(b"@"):
                errors.append(f"record {n} (line {line_no}): header does not start with '@'")
            elif not plus.startswith(b"+"):
                errors.append(f"record {n}: third line does not start with '+'")
            elif len(seq) != len(qual):
                errors.append(f"record {n}: sequence length {len(seq)} != quality length {len(qual)}")
            elif not FASTQ_SEQ_RE.match(seq):
                errors.append(f"record {n}: unexpected characters in sequence")
            if errors and len(errors) >= 5:
                break
            lengths[len(seq)] += 1
            if qual:
                min_q, max_q = min(min_q, min(qual)), max(max_q, max(qual))
            names.append(_read_name(header))
            if n % 20000 == 0:
                check_cancel()
    if n == 0:
        errors.append("file contains no reads")
    encoding = "Phred+33 (Sanger/Illumina 1.8+)" if min_q < 59 else ("Phred+64 (old Illumina)" if min_q >= 64 else "ambiguous")
    return {
        "records_checked": n, "errors": errors, "read_length_min": min(lengths) if lengths else 0,
        "read_length_max": max(lengths) if lengths else 0,
        "read_length_mode": lengths.most_common(1)[0][0] if lengths else 0,
        "quality_chars": f"{chr(min_q)}..{chr(max_q)}" if n else "",
        "quality_encoding": encoding,
        "names": names,
    }


class FastqIntegrityStep(Step):
    key = "fastq_integrity"
    name = "FASTQ integrity & format"
    description = "Full gzip integrity test (pigz -t), FASTQ record structure, quality encoding, R1/R2 pairing."
    tools = ("pigz",)

    def applicable(self, ctx):
        return (bool(ctx.files("fastq_r1", "fastq_r2", "fastq_se")), "No FASTQ files registered.")

    def execute(self, ctx: StepContext, runner) -> StepResult:
        fastqs = ctx.files("fastq_r1", "fastq_r2", "fastq_se")
        max_records = int(ctx.params.get("fastq_check_records", 100_000))
        warnings, summary, scans = [], {}, {}
        for f in fastqs:
            p = Path(f.path)
            if p.name.endswith(".gz"):
                ctx.log(f"Testing gzip integrity of the whole file: {p.name}")
                runner.run(["pigz", "-t", "-p", str(min(ctx.threads, 4)), p])
            ctx.log(f"Checking FASTQ structure of the first {max_records:,} records: {p.name}")
            scan = scan_fastq(p, max_records, ctx.check_cancel)
            if scan["errors"]:
                f.validation_status, f.validation_message = "error", "; ".join(scan["errors"])
                db.session.commit()
                raise StepFailure(f"{p.name}: malformed FASTQ — " + "; ".join(scan["errors"]))
            if "Phred+33" not in scan["quality_encoding"]:
                warnings.append(f"{p.name}: quality encoding looks {scan['quality_encoding']}.")
            scans[f.role] = scan
            f.validation_status = "ok"
            f.validation_message = (f"gzip OK; {scan['records_checked']:,} records checked; "
                                    f"read length {scan['read_length_min']}-{scan['read_length_max']}")
            db.session.commit()
            summary[p.name] = {k: v for k, v in scan.items() if k != "names"}
            ctx.log(f"  OK: {summary[p.name]}")

        if "fastq_r1" in scans and "fastq_r2" in scans:
            r1, r2 = scans["fastq_r1"]["names"], scans["fastq_r2"]["names"]
            n = min(len(r1), len(r2))
            mismatches = sum(1 for a, b in zip(r1[:n], r2[:n]) if a != b)
            summary["pairing"] = {"records_compared": n, "name_mismatches": mismatches}
            ctx.log(f"R1/R2 read-name pairing: {mismatches} mismatches in {n:,} records")
            if mismatches:
                raise StepFailure(
                    f"R1 and R2 are out of sync: {mismatches} of the first {n:,} read names differ. "
                    "Files may be from different samples/lanes or one was truncated/re-sorted."
                )
            if len(r1) != len(r2):
                warnings.append("R1 and R2 have different numbers of records in the checked region.")
        elif ctx.sample.library_layout == "paired" and ("fastq_r1" in scans) != ("fastq_r2" in scans):
            warnings.append("Sample is paired-end but only one of R1/R2 is registered.")

        out = _write_json(ctx.step_dir("intake") / "fastq_integrity.json", summary)
        return StepResult(outputs={"fastq_report": out}, summary=summary, warnings=warnings)


# ---------------------------------------------------------------------------
def _contig_check(f, contigs: dict, ref_contigs: dict, kind: str) -> tuple[dict, list[str]]:
    comparison = compare_to_reference(contigs, ref_contigs)
    warnings = []
    f.detected_build = comparison["file_build"]
    if comparison["verdict"] == "incompatible":
        f.validation_status = "error"
        f.validation_message = "INCOMPATIBLE with reference: " + comparison["reason"]
        warnings.append(f"{Path(f.path).name}: {comparison['reason']} "
                        f"This {kind} must not be compared directly with results on our reference.")
    elif comparison["verdict"] == "primary_compatible":
        f.validation_status = "warning"
        f.validation_message = comparison["reason"]
        warnings.append(f"{Path(f.path).name}: {comparison['reason']}")
    else:
        f.validation_status = "ok"
        f.validation_message = comparison["reason"]
    return comparison, warnings


class BamValidationStep(Step):
    key = "bam_validation"
    name = "BAM validation & build check"
    description = ("samtools quickcheck (header + EOF marker), sort order, read groups, index, "
                   "aligner provenance, contig/genome-build compatibility with the reference.")
    tools = ("samtools",)

    def applicable(self, ctx):
        return (bool(ctx.files("bam")), "No BAM files registered.")

    def execute(self, ctx: StepContext, runner) -> StepResult:
        ref_contigs = ctx.shared.get("reference_contigs") or reference_contigs(ctx.reference)
        warnings, summary = [], {}
        for f in ctx.files("bam"):
            p = Path(f.path)
            ctx.log(f"samtools quickcheck: {p.name}")
            try:
                runner.run(["samtools", "quickcheck", "-v", p])
            except Exception:
                f.validation_status, f.validation_message = "error", "samtools quickcheck failed (truncated or corrupt BAM)"
                db.session.commit()
                raise
            hdr = bam_header_summary(p)
            comparison, w = _contig_check(f, hdr["contigs"], ref_contigs, "BAM")
            warnings += w
            name = p.name
            if hdr["sort_order"] != "coordinate":
                warnings.append(f"{name}: not coordinate-sorted (SO:{hdr['sort_order']}).")
            if not hdr["read_groups"]:
                warnings.append(f"{name}: no @RG read group — GATK will refuse this BAM; "
                                "duplicate marking cannot be done per library.")
            if not hdr["index"]:
                warnings.append(f"{name}: no index (.bai) — region queries and most callers need one "
                                "(samtools index).")
            if len(hdr["samples"]) > 1:
                warnings.append(f"{name}: contains several samples {hdr['samples']}.")
            programs = [f"{pg['id']} {pg['version']}".strip() for pg in hdr["programs"] if pg["id"] != "samtools"]
            summary[name] = {
                "origin": f.origin, "verdict": comparison["verdict"], "build": comparison["file_build"],
                "contigs": comparison["file_contigs"], "contig_classes": comparison["file_contig_classes"],
                "missing_in_reference": comparison["missing_in_reference_count"],
                "reason": comparison["reason"], "sort_order": hdr["sort_order"],
                "read_groups": hdr["read_groups"], "samples": hdr["samples"],
                "indexed": bool(hdr["index"]), "programs": programs[:6],
            }
            db.session.commit()
            ctx.log(f"  {name}: build={comparison['file_build']} verdict={comparison['verdict']} "
                    f"RG={len(hdr['read_groups'])} indexed={bool(hdr['index'])}")
        out = _write_json(ctx.step_dir("intake") / "bam_validation.json", summary)
        return StepResult(outputs={"bam_report": out}, summary=summary, warnings=warnings)


# ---------------------------------------------------------------------------
class VcfValidationStep(Step):
    key = "vcf_validation"
    name = "VCF validation & build check"
    description = "bcftools header parse, sample names, ##contig build compatibility, index presence."
    tools = ("bcftools",)

    def applicable(self, ctx):
        return (bool(ctx.files("vcf")), "No VCF files registered.")

    def execute(self, ctx: StepContext, runner) -> StepResult:
        ref_contigs = ctx.shared.get("reference_contigs") or reference_contigs(ctx.reference)
        warnings, summary = [], {}
        out_dir = ctx.step_dir("intake")
        for f in ctx.files("vcf"):
            p = Path(f.path)
            header_out = out_dir / f"{p.name}.header.txt"
            runner.run(["bcftools", "view", "-h", p], stdout_path=header_out)
            hdr = vcf_header_summary(p)
            if hdr["contigs"]:
                comparison, w = _contig_check(f, hdr["contigs"], ref_contigs, "VCF")
                warnings += w
            else:
                comparison = {"verdict": "unknown", "file_build": "unknown",
                              "reason": "VCF has no ##contig lengths; genome build cannot be verified."}
                f.validation_status, f.validation_message = "warning", comparison["reason"]
                warnings.append(f"{p.name}: {comparison['reason']}")
            if not hdr["index"]:
                warnings.append(f"{p.name}: no .tbi/.csi index.")
            summary[p.name] = {"origin": f.origin, "verdict": comparison["verdict"],
                               "build": comparison["file_build"], "samples": hdr["samples"],
                               "source": hdr["source"], "indexed": bool(hdr["index"]),
                               "reason": comparison["reason"]}
            db.session.commit()
        out = _write_json(out_dir / "vcf_validation.json", summary)
        return StepResult(outputs={"vcf_report": out}, summary=summary, warnings=warnings)


# ---------------------------------------------------------------------------
class IntakeReportStep(Step):
    key = "intake_report"
    name = "Intake report"
    description = "Combined HTML + JSON report of all checks."

    def execute(self, ctx: StepContext, runner) -> StepResult:
        run = ctx.run
        data = {
            "run_id": run.run_id, "sample": ctx.sample.sample_id,
            "reference": {"name": ctx.reference.name, "path": ctx.reference.fasta_path,
                          "build": ctx.reference.build},
            "files": [
                {"role": f.role_label, "origin": f.origin, "path": f.path, "size": f.size_bytes,
                 "sha256": f.sha256, "build": f.detected_build, "status": f.validation_status,
                 "message": f.validation_message}
                for f in ctx.input_files
            ],
            "steps": [{"name": s.name, "status": s.status, "summary": s.summary, "warnings": s.warnings,
                       "message": s.message} for s in run.steps if s.key != self.key],
            "generated_at": utcnow().isoformat() + "Z",
        }
        reports = ctx.step_dir("reports")
        json_out = _write_json(reports / "intake_report.json", data)
        html_out = reports / "intake_report.html"
        html_out.write_text(render_template("reports/intake_report.html", data=data))
        return StepResult(outputs={"report_json": json_out, "report_html": str(html_out)})
