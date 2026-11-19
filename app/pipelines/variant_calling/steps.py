"""Germline variant calling (one step per selected caller), normalisation and reporting.

For every caller the step produces, in ``vcf/<caller>/``:
  <stem>.<caller>.raw.vcf.gz       exactly what the caller wrote (+ .tbi)
  <stem>.<caller>.norm.vcf.gz      left-aligned, multi-allelic sites split (bcftools norm) (+ .tbi)
  <stem>.<caller>.stats.txt        bcftools stats on the normalised VCF

<stem> is the sample ID, plus the preset when the run was started from one
(e.g. SAMPLE_A.provider-like-dragen-4-3.gatk_hc.filtered.vcf.gz) so files from
different parameter sets can be told apart after download.

Why normalise? The same indel can be written in several equivalent ways
(e.g. in a repeat, at different positions). Left-alignment + splitting
multi-allelic records gives one canonical representation, which is required
before VCFs from different callers can be compared (Phase 4).

These are RAW, UNFILTERED calls for research / pipeline development. They are
not clinically validated.
"""
from __future__ import annotations

import os
from pathlib import Path

from flask import render_template

from ...engine.step import Step, StepContext, StepFailure, StepResult
from ...models import CommandExecution
from ...models.base import utcnow
from ..bam_processing.steps import analysis_bam
from ..common import check_vcf, output_stem, parse_bcftools_stats, prior_output, write_json
from .callers import CallJob, get_caller


def _psc(stats_path: Path) -> dict:
    """Per-sample counts line (PSC) from bcftools stats -s -."""
    with open(stats_path) as fh:
        for line in fh:
            if line.startswith("PSC\t"):
                p = line.rstrip("\n").split("\t")
                return {"hom_ref": int(p[3]), "hom_alt": int(p[4]), "het": int(p[5]),
                        "mean_depth": float(p[9])}
    return {}


class VariantCallStep(Step):
    def __init__(self, caller_key: str):
        self.caller = get_caller(caller_key)
        self.key = f"call_{caller_key}"
        self.name = f"Variant calling: {self.caller.name}"
        self.description = self.caller.description
        self.tools = (self.caller.tool, "bcftools")

    def execute(self, ctx: StepContext, runner) -> StepResult:
        ok, why = self.caller.available()
        if not ok:
            raise StepFailure(f"{self.caller.name} is not available: {why}")
        sid = ctx.sample.sample_id          # VCF sample name
        stem = output_stem(ctx)             # file names
        k = self.caller.key
        d = ctx.step_dir(f"vcf/{k}")
        ref = Path(ctx.reference.fasta_path)
        raw = d / f"{stem}.{k}.raw.vcf.gz"
        norm = d / f"{stem}.{k}.norm.vcf.gz"
        stats = d / f"{stem}.{k}.stats.txt"
        cparams = (ctx.params.get("caller_params") or {}).get(k, {})
        job = CallJob(bam=analysis_bam(ctx), reference=ref,
                      regions_bed=Path(prior_output(ctx, "regions", "calling_bed")),
                      out_vcf=raw, work_dir=ctx.step_dir(f"vcf/{k}/tmp"), threads=ctx.threads,
                      sample_id=sid, params=cparams, label=ctx.run.run_id)
        ctx.log(f"Caller parameters: {cparams}")
        for cmd in self.caller.build_commands(job):
            try:
                runner.run(cmd.commands, stdout_path=cmd.stdout_path, env=cmd.env)
            except BaseException:
                if cmd.cleanup:   # e.g. remove a container left behind by a cancelled run
                    try:
                        runner.run(cmd.cleanup, ok_exit_codes=tuple(range(256)))
                    except Exception:  # noqa: BLE001 - cleanup is best effort, failure already raised
                        pass
                raise
        t = str(min(ctx.threads, 4))
        runner.run(["bcftools", "index", "-f", "-t", raw])
        runner.run(["bcftools", "norm", "--threads", t, "-f", ref, "-m", "-any", "--check-ref", "w",
                    "-Oz", "-o", norm, raw])
        runner.run(["bcftools", "index", "-t", norm])
        runner.run(["bcftools", "stats", "-F", ref, "-s", "-", norm], stdout_path=stats)
        try:
            job.work_dir.rmdir()
        except OSError:
            pass
        st = parse_bcftools_stats(stats)
        psc = _psc(stats)
        summary = {
            "caller": self.caller.name, "params": cparams,
            "records": st.get("records", 0), "snvs": st.get("SNPs", 0), "indels": st.get("indels", 0),
            "mnps": st.get("MNPs", 0), "others": st.get("others", 0),
            "multiallelic_sites": st.get("multiallelic_sites", 0),
            "ts_tv": st.get("ts_tv"), **psc,
            "filter_state": "raw/unfiltered (see the filtering step)",
            **({"provenance": self.caller.version_info()} if self.caller.version_info() else {}),
        }
        if psc.get("hom_alt"):
            summary["het_hom_ratio"] = round(psc["het"] / psc["hom_alt"], 2)
        warnings = []
        if summary["records"] == 0:
            warnings.append(f"{self.caller.name} produced no variants in the calling regions.")
        tstv = summary.get("ts_tv")
        if tstv and summary["snvs"] > 500 and not 1.8 <= tstv <= 3.6:
            warnings.append(f"{self.caller.name}: Ti/Tv = {tstv} is outside the expected range "
                            "(~2.0-2.1 genome-wide, ~2.6-3.3 in exons) — possible excess of false positives.")
        write_json(d / f"{stem}.{k}.summary.json", summary)
        return StepResult(outputs={"raw_vcf": str(raw), "norm_vcf": str(norm), "stats": str(stats)},
                          summary=summary, warnings=warnings)

    def validate(self, ctx, result):
        super().validate(ctx, result)
        check_vcf(Path(result.outputs["raw_vcf"]))
        check_vcf(Path(result.outputs["norm_vcf"]))


# ---------------------------------------------------------------------------
class VariantFilterStep(Step):
    """Soft filtering: every record is kept; failing records get a FILTER label.

    Why soft? A hard filter deletes evidence. With soft filters the FILTER column
    shows *why* a call is doubtful, comparisons can use PASS-only or all calls, and
    a reviewer can still see filtered variants in IGV.
    """

    def __init__(self, caller_key: str):
        self.caller = get_caller(caller_key)
        self.key = f"filter_{caller_key}"
        self.name = f"Variant filtering: {self.caller.name}"
        self.description = ("Caller-specific soft filters (QUAL, depth, genotype quality, allele fraction, "
                            "caller-specific rules). Records are labelled, never removed.")
        self.tools = ("bcftools",)

    def execute(self, ctx: StepContext, runner) -> StepResult:
        k = self.caller.key
        stem = output_stem(ctx)
        d = ctx.step_dir(f"vcf/{k}")
        norm = prior_output(ctx, f"call_{k}", "norm_vcf")
        cparams = (ctx.params.get("caller_params") or {}).get(k, {})
        rules = self.caller.filter_rules(cparams)
        out = d / f"{stem}.{k}.filtered.vcf.gz"
        pipe = []
        for i, r in enumerate(rules):
            argv = ["bcftools", "filter", "-m", "+", "-s", r.name, "-e", r.expr]
            argv += ["-Oz", "-o", out] if i == len(rules) - 1 else ["-Ou"]
            if i == 0:
                argv.append(norm)
            pipe.append(argv)
        for r in rules:
            ctx.log(f"Filter {r.name}: {r.description}   [{r.expr}]")
        runner.run(pipe)
        runner.run(["bcftools", "index", "-t", out])
        ref = ctx.reference.fasta_path
        stats = d / f"{stem}.{k}.pass.stats.txt"
        runner.run(["bcftools", "stats", "-f", "PASS", "-F", ref, "-s", "-", out], stdout_path=stats)
        import pysam
        counts: dict = {}
        total = passed = 0
        with pysam.VariantFile(str(out)) as fh:
            for rec in fh:
                total += 1
                f = [x for x in rec.filter.keys() if x != "PASS"]
                if not f:
                    passed += 1
                for x in f:
                    counts[x] = counts.get(x, 0) + 1
        st = parse_bcftools_stats(stats)
        psc = _psc(stats)
        summary = {
            "caller": self.caller.name,
            "rules": [{"name": r.name, "expr": r.expr, "description": r.description} for r in rules],
            "records": total, "pass": passed, "pass_pct": round(100 * passed / total, 2) if total else None,
            "filtered_by": dict(sorted(counts.items(), key=lambda x: -x[1])),
            "pass_snvs": st.get("SNPs", 0), "pass_indels": st.get("indels", 0),
            "pass_ts_tv": st.get("ts_tv"), "pass_het": psc.get("het"), "pass_hom_alt": psc.get("hom_alt"),
            "filter_state": "soft-filtered (FILTER column set; no record removed)",
        }
        warnings = []
        if total and passed / total < 0.5:
            warnings.append(f"{self.caller.name}: only {summary['pass_pct']}% of records PASS the filters "
                            "— check the filter thresholds or the data quality.")
        write_json(d / f"{stem}.{k}.filter_summary.json", summary)
        return StepResult(outputs={"filtered_vcf": str(out), "pass_stats": str(stats)},
                          summary=summary, warnings=warnings)

    def validate(self, ctx, result):
        super().validate(ctx, result)
        check_vcf(Path(result.outputs["filtered_vcf"]))


# ---------------------------------------------------------------------------
class MultiQCStep(Step):
    key = "multiqc"
    name = "MultiQC summary"
    description = "Aggregates FastQC, fastp, samtools, mosdepth and bcftools reports into one HTML page."
    tools = ("multiqc",)

    def execute(self, ctx: StepContext, runner) -> StepResult:
        out = ctx.step_dir("reports/multiqc")
        runner.run(["multiqc", "--no-ansi", "--outdir", out, "--filename", "multiqc_report.html",
                    "--title", f"{ctx.sample.sample_id} — {ctx.run.run_id}",
                    "--ignore", "reports/*", ctx.output_dir])
        return StepResult(outputs={"html": str(out / "multiqc_report.html")})


class RunReportStep(Step):
    key = "report"
    name = "Run report"
    description = "HTML + JSON report: QC metrics, read accounting, variant counts, tool versions, warnings."

    def execute(self, ctx: StepContext, runner) -> StepResult:
        from ...engine.registry import get_workflow

        run = ctx.run
        steps = [s for s in run.steps if s.key != self.key]
        versions = {}
        for s in steps:
            for c in CommandExecution.query.filter_by(step_pk=s.id).all():
                versions.update(c.tool_versions or {})
        by_key = {s.key: s for s in steps}
        d_reports = ctx.step_dir("reports")
        accounting = []
        if "fastp" in by_key and by_key["fastp"].summary:
            fp = by_key["fastp"].summary
            accounting.append(("Raw reads (R1+R2)", fp.get("reads_before")))
            accounting.append(("Removed by fastp (kept in failed_reads.fastq.gz)", fp.get("reads_removed")))
            accounting.append(("Reads entering alignment", fp.get("reads_after")))
        if "bwa_align" in by_key and by_key["bwa_align"].summary:
            accounting.append(("Primary alignments in BAM", by_key["bwa_align"].summary.get("primary")))
        if "markdup" in by_key and by_key["markdup"].summary:
            md = by_key["markdup"].summary
            accounting.append(("Primary records after duplicate marking", md.get("primary")))
            accounting.append(("Of which flagged duplicate", md.get("duplicates")))
        callers = []
        for s in steps:
            if s.key.startswith("call_") and s.summary:
                f = by_key.get("filter_" + s.key[5:])
                callers.append(s.summary | {"key": s.key,
                                            "filter": f.summary if f is not None and f.summary else None})
        from ...services.limitations import run_limitations
        comparison = None
        report_links = []
        for key in ("compare_callers", "vcf_compare", "bam_compare", "read_evidence", "variant_check"):
            if key in by_key and by_key[key].status == "completed":
                if key in ("compare_callers", "vcf_compare"):
                    comparison = by_key[key].summary
                html_out = (by_key[key].outputs or {}).get("report_html")
                if html_out:
                    report_links.append((by_key[key].name, os.path.relpath(html_out, d_reports)))
        data = {
            "run_id": run.run_id, "workflow": get_workflow(run.workflow_type).name,
            "workflow_version": run.workflow_version, "sample": ctx.sample.sample_id,
            "reference": {"name": ctx.reference.name, "path": ctx.reference.fasta_path,
                          "build": ctx.reference.build},
            "params": run.params,
            "inputs": [{"role": f.role_label, "path": f.path, "sha256": f.sha256} for f in ctx.input_files],
            "steps": [{"key": s.key, "name": s.name, "status": s.status, "summary": s.summary,
                       "warnings": s.warnings, "message": s.message,
                       "duration": (s.finished_at - s.started_at).total_seconds()
                       if s.finished_at and s.started_at else None} for s in steps],
            "accounting": accounting, "callers": callers, "tool_versions": versions,
            "warnings": list(run.warnings or []),
            "comparison": comparison, "report_links": report_links,
            "limitations": run_limitations(ctx.reference, run.params, ctx.sample),
            "generated_at": utcnow().isoformat() + "Z",
        }
        d = ctx.step_dir("reports")
        json_out = write_json(d / "run_report.json", data)
        html = d / "run_report.html"
        html.write_text(render_template("reports/run_report.html", data=data))
        return StepResult(outputs={"report_json": json_out, "report_html": str(html)})
