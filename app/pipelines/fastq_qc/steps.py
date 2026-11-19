"""FASTQ quality control and adapter/quality trimming.

* FastQC   - per-file QC report (base quality, GC, adapters, duplication...).
             Run on raw reads and again on trimmed reads to show the effect.
* fastp    - adapter + quality trimming and read filtering. Reads that are
             removed are NOT thrown away silently: they are written to
             ``failed_reads.fastq.gz`` together with the reason, and the
             counts are reported and reconciled later against the BAM.
"""
from __future__ import annotations

from pathlib import Path

from ...engine.step import Step, StepContext, StepFailure, StepResult
from ..common import check_fastq_gz, prior_output, read_json, write_json


def _pair(ctx: StepContext) -> tuple[Path, Path]:
    r1, r2 = ctx.files("fastq_r1"), ctx.files("fastq_r2")
    if len(r1) != 1 or len(r2) != 1:
        raise StepFailure("Exactly one R1 and one R2 FASTQ are required.")
    return Path(r1[0].path), Path(r2[0].path)


def analysis_fastqs(ctx: StepContext) -> tuple[Path, Path]:
    """Trimmed reads if trimming ran, otherwise the raw input reads."""
    if ctx.params.get("trimming", True):
        return (Path(prior_output(ctx, "fastp", "trimmed_r1")), Path(prior_output(ctx, "fastp", "trimmed_r2")))
    return _pair(ctx)


def _fastqc_summary(zip_path: Path) -> tuple[dict, int | None]:
    """Read PASS/WARN/FAIL per FastQC module and the total read count from the zip."""
    import zipfile
    modules, total = {}, None
    try:
        with zipfile.ZipFile(zip_path) as zf:
            names = zf.namelist()
            summ = next(n for n in names if n.endswith("/summary.txt"))
            for line in zf.read(summ).decode().splitlines():
                status, module, _ = line.split("\t")
                modules[module] = status
            data = next(n for n in names if n.endswith("/fastqc_data.txt"))
            for line in zf.read(data).decode().splitlines():
                if line.startswith("Total Sequences"):
                    total = int(line.split("\t")[1])
                    break
    except (OSError, StopIteration, ValueError, zipfile.BadZipFile):
        pass
    return modules, total


class FastQCStep(Step):
    tools = ("fastqc",)

    def __init__(self, stage: str):
        assert stage in ("raw", "trimmed")
        self.stage = stage
        self.key = f"fastqc_{stage}"
        self.name = f"FastQC ({stage} reads)"
        self.description = "Per-file read quality report."

    def applicable(self, ctx):
        if self.stage == "trimmed" and not ctx.params.get("trimming", True):
            return False, "Trimming disabled — no trimmed reads to check."
        return True, ""

    def execute(self, ctx: StepContext, runner) -> StepResult:
        files = list(_pair(ctx)) if self.stage == "raw" else list(analysis_fastqs(ctx))
        out_dir = ctx.step_dir(f"qc/fastqc_{self.stage}")
        runner.run(["fastqc", "--threads", str(min(ctx.threads, len(files))),
                    "--outdir", out_dir, *files])
        outputs, summary, warnings = {}, {}, []
        for i, f in enumerate(files, 1):
            base = f.name.removesuffix(".gz").removesuffix(".fastq").removesuffix(".fq")
            html, zp = out_dir / f"{base}_fastqc.html", out_dir / f"{base}_fastqc.zip"
            outputs[f"r{i}_html"], outputs[f"r{i}_zip"] = str(html), str(zp)
            mods, total = _fastqc_summary(zp)
            summary[f.name] = {"total_sequences": total, "modules": mods}
            fails = [m for m, s in mods.items() if s == "FAIL"]
            if fails:
                # FAIL flags are heuristics tuned for WGS libraries; for exomes
                # "Sequence Duplication" and "Per base sequence content" often FAIL normally.
                warnings.append(f"{f.name}: FastQC FAIL in {', '.join(fails)} (review the report; some "
                                "modules routinely fail for exome/amplicon libraries).")
        return StepResult(outputs=outputs, summary=summary, warnings=warnings)


class FastpStep(Step):
    key = "fastp"
    name = "Adapter & quality trimming (fastp)"
    description = "Adapter removal, poly-G trimming (NovaSeq), low-quality/short read filtering."
    tools = ("fastp",)

    def applicable(self, ctx):
        return (bool(ctx.params.get("trimming", True)), "Trimming disabled by parameter.")

    def execute(self, ctx: StepContext, runner) -> StepResult:
        r1, r2 = _pair(ctx)
        d = ctx.step_dir("fastq")
        sid = ctx.sample.sample_id
        out1, out2 = d / f"{sid}.trimmed_R1.fastq.gz", d / f"{sid}.trimmed_R2.fastq.gz"
        failed = d / f"{sid}.failed_reads.fastq.gz"
        js, html = d / f"{sid}.fastp.json", d / f"{sid}.fastp.html"
        p = ctx.params
        argv = ["fastp", "--in1", r1, "--in2", r2, "--out1", out1, "--out2", out2,
                "--failed_out", failed, "--json", js, "--html", html,
                "--report_title", f"fastp {sid}", "--thread", str(min(ctx.threads, 16)),
                "--qualified_quality_phred", str(p.get("fastp_qualified_quality", 15)),
                "--unqualified_percent_limit", str(p.get("fastp_unqualified_percent", 40)),
                "--n_base_limit", str(p.get("fastp_n_base_limit", 5)),
                "--length_required", str(p.get("fastp_min_length", 36))]
        if p.get("fastp_adapter_trimming", True):
            argv.append("--detect_adapter_for_pe")
        else:
            argv.append("--disable_adapter_trimming")
        runner.run(argv)

        rep = read_json(js)
        before, after = rep["summary"]["before_filtering"], rep["summary"]["after_filtering"]
        filt = rep.get("filtering_result", {})
        removed = before["total_reads"] - after["total_reads"]
        summary = {
            "reads_before": before["total_reads"], "reads_after": after["total_reads"],
            "reads_removed": removed,
            "reads_removed_pct": round(100 * removed / before["total_reads"], 3) if before["total_reads"] else 0,
            "removal_reasons": filt,
            "q30_before_pct": round(100 * before["q30_rate"], 2), "q30_after_pct": round(100 * after["q30_rate"], 2),
            "gc_pct": round(100 * after["gc_content"], 2),
            "mean_length_after": [after.get("read1_mean_length"), after.get("read2_mean_length")],
            "adapter_trimmed_reads": rep.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0),
            "duplication_rate_pct": round(100 * rep.get("duplication", {}).get("rate", 0), 2),
            "insert_size_peak": rep.get("insert_size", {}).get("peak"),
        }
        warnings = []
        if summary["reads_removed_pct"] > 10:
            warnings.append(f"fastp removed {summary['reads_removed_pct']}% of reads (> 10%) — check "
                            "library quality; removed reads are kept in failed_reads.fastq.gz.")
        if after["total_reads"] == 0:
            raise StepFailure("No reads passed fastp filtering.")
        write_json(d / "fastp_summary.json", summary)
        return StepResult(
            outputs={"trimmed_r1": str(out1), "trimmed_r2": str(out2), "failed_reads": str(failed),
                     "json": str(js), "html": str(html)},
            summary=summary, warnings=warnings,
        )

    def validate(self, ctx, result):
        # failed_reads may legitimately be empty (gzip header only) -> only check existence
        for name in ("trimmed_r1", "trimmed_r2", "json"):
            p = Path(result.outputs.get(name, ""))
            if not p.is_file() or p.stat().st_size == 0:
                raise StepFailure(f"fastp output '{name}' missing or empty: {p}")
        if not Path(result.outputs.get("failed_reads", "")).exists():
            raise StepFailure("fastp failed-reads file missing.")
        check_fastq_gz(Path(result.outputs["trimmed_r1"]))
        check_fastq_gz(Path(result.outputs["trimmed_r2"]))
