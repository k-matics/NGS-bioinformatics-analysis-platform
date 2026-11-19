"""Read alignment and BAM post-processing.

Order (GATK/clinical convention):
  1. BWA-MEM alignment -> samtools fixmate -> coordinate sort   (one streamed step,
     so no huge unsorted BAM is written to disk)
  2. Duplicate MARKING (flag 0x400). Reads are never removed here.
  3. Optional duplicate REMOVAL -> a *separate* BAM; the marked BAM is kept.
  4. BAM index (.bai)

Read accounting: after alignment and after duplicate marking we compare the
number of primary records with the number of reads that entered the step.
Any difference fails the run - reads must never be lost silently.
"""
from __future__ import annotations

from pathlib import Path

from ...engine.step import Step, StepContext, StepFailure, StepResult
from ...services.settings_service import get_setting
from ...services.tool_service import gatk_env
from ..common import (
    check_index_fresh,
    parse_flagstat_json,
    prior_output,
    prior_summary,
    read_group_from_fastq,
    rg_line,
    samtools_quickcheck,
    write_json,
)
from ..fastq_qc.steps import analysis_fastqs


def flagstat(runner, ctx: StepContext, bam: Path, out: Path) -> dict:
    runner.run(["samtools", "flagstat", "-@", str(min(ctx.threads, 4)), "-O", "json", bam], stdout_path=out)
    return parse_flagstat_json(out)


def expected_input_reads(ctx: StepContext) -> int | None:
    """Reads entering alignment: fastp output if trimming ran, else FastQC raw totals."""
    if ctx.params.get("trimming", True):
        return prior_summary(ctx, "fastp").get("reads_after")
    totals = [v.get("total_sequences") for v in prior_summary(ctx, "fastqc_raw").values()
              if isinstance(v, dict)]
    return sum(totals) if totals and all(t is not None for t in totals) else None


class BwaAlignStep(Step):
    key = "bwa_align"
    name = "Alignment (BWA-MEM) + coordinate sort"
    description = "bwa mem with read group -> samtools fixmate -m -> samtools sort (streamed, no shell)."
    tools = ("bwa", "samtools")

    def execute(self, ctx: StepContext, runner) -> StepResult:
        ref = Path(ctx.reference.fasta_path)
        r1, r2 = analysis_fastqs(ctx)
        sid = ctx.sample.sample_id
        rg = read_group_from_fastq(r1, sid, ctx.params.get("library_id") or None)
        d = ctx.step_dir("bam")
        tmp = ctx.step_dir("bam/tmp_sort")
        out = d / f"{sid}.aligned.sorted.bam"
        sort_threads = max(1, min(4, ctx.threads // 2))
        ctx.log(f"Read group: {rg}")
        runner.run([
            ["bwa", "mem", "-t", str(ctx.threads), "-K", "100000000", "-Y", "-R", rg_line(rg), ref, r1, r2],
            ["samtools", "fixmate", "-m", "-u", "-@", "2", "-", "-"],
            ["samtools", "sort", "-@", str(sort_threads), "-m", "1G", "-T", tmp / "sort", "-o", out, "-"],
        ])
        try:
            tmp.rmdir()   # empty after a successful sort
        except OSError:
            pass
        samtools_quickcheck(out)
        fs = flagstat(runner, ctx, out, d / f"{sid}.aligned.flagstat.json")
        warnings = []
        expected = expected_input_reads(ctx)
        if expected is not None and fs["primary"] != expected:
            raise StepFailure(f"Read accounting mismatch: {expected:,} reads entered alignment but the BAM "
                              f"has {fs['primary']:,} primary records. Reads were lost.")
        if expected is None:
            warnings.append("Input read count unknown — read accounting could not be verified.")
        if fs["primary_mapped_pct"] < 90:
            warnings.append(f"Only {fs['primary_mapped_pct']}% of reads mapped (< 90%): possible contamination, "
                            "wrong reference/species or poor library.")
        summary = {"read_group": rg, "reads_in": expected, **fs}
        write_json(d / "alignment_summary.json", summary)
        return StepResult(outputs={"bam": str(out), "flagstat": str(d / f"{sid}.aligned.flagstat.json")},
                          summary=summary, warnings=warnings)

    def validate(self, ctx, result):
        super().validate(ctx, result)
        samtools_quickcheck(Path(result.outputs["bam"]))


class MarkDuplicatesStep(Step):
    key = "markdup"
    name = "Duplicate marking"
    description = "Flags PCR/optical duplicates (SAM flag 0x400). No reads are removed."

    @property
    def tools(self):
        return ("samtools", "gatk")

    def execute(self, ctx: StepContext, runner) -> StepResult:
        sid = ctx.sample.sample_id
        aligned = Path(prior_output(ctx, "bwa_align", "bam"))
        d = ctx.step_dir("bam")
        out = d / f"{sid}.markdup.bam"
        metrics = d / f"{sid}.markdup.metrics.txt"
        tool = ctx.params.get("markdup_tool", "samtools")
        pixel = str(ctx.params.get("optical_pixel_distance", 2500))
        if tool == "gatk":
            gatk = get_setting("gatk_path")
            runner.run([gatk, "--java-options", f"-Xmx{int(ctx.params.get('java_xmx_gb', 8))}g",
                        "MarkDuplicates", "-I", aligned, "-O", out, "-M", metrics,
                        "--OPTICAL_DUPLICATE_PIXEL_DISTANCE", pixel, "--CREATE_INDEX", "false",
                        "--TMP_DIR", ctx.step_dir("bam/tmp_markdup")],
                       env=gatk_env(get_setting("gatk_java_home", "")))
        else:
            runner.run(["samtools", "markdup", "-@", str(min(ctx.threads, 8)), "-d", pixel,
                        "-f", metrics, "-T", d / "markdup_tmp", aligned, out])
        samtools_quickcheck(out)
        fs = flagstat(runner, ctx, out, d / f"{sid}.markdup.flagstat.json")
        before = prior_summary(ctx, "bwa_align").get("primary")
        if before is not None and fs["primary"] != before:
            raise StepFailure(f"Duplicate marking changed the number of primary records "
                              f"({before:,} -> {fs['primary']:,}). Reads must not be removed here.")
        warnings = []
        if fs.get("duplicate_pct", 0) > 30:
            warnings.append(f"High duplicate rate {fs['duplicate_pct']}% (> 30%): low library complexity.")
        summary = {"tool": tool, "optical_pixel_distance": int(pixel), **fs}
        return StepResult(outputs={"bam": str(out), "metrics": str(metrics),
                                   "flagstat": str(d / f"{sid}.markdup.flagstat.json")},
                          summary=summary, warnings=warnings)

    def validate(self, ctx, result):
        super().validate(ctx, result)
        samtools_quickcheck(Path(result.outputs["bam"]))


class IndexBamStep(Step):
    key = "index_bam"
    name = "BAM index"
    description = "samtools index on the duplicate-marked BAM (needed for region access and callers)."
    tools = ("samtools",)

    def execute(self, ctx: StepContext, runner) -> StepResult:
        bam = Path(prior_output(ctx, "markdup", "bam"))
        runner.run(["samtools", "index", "-@", str(min(ctx.threads, 8)), bam])
        return StepResult(outputs={"bam": str(bam), "bai": str(bam) + ".bai"})

    def validate(self, ctx, result):
        super().validate(ctx, result)
        check_index_fresh(Path(result.outputs["bam"]))


class RemoveDuplicatesStep(Step):
    key = "remove_duplicates"
    name = "Duplicate removal (separate BAM)"
    description = "Optional: writes a second BAM without duplicate-flagged reads. The marked BAM is kept."
    tools = ("samtools",)

    def applicable(self, ctx):
        return (bool(ctx.params.get("remove_duplicates")), "Not requested (callers ignore flagged duplicates).")

    def execute(self, ctx: StepContext, runner) -> StepResult:
        sid = ctx.sample.sample_id
        marked = Path(prior_output(ctx, "markdup", "bam"))
        d = ctx.step_dir("bam")
        out = d / f"{sid}.dedup.bam"
        runner.run(["samtools", "view", "-@", str(min(ctx.threads, 8)), "-b", "-F", "0x400",
                    "-o", out, marked])
        runner.run(["samtools", "index", "-@", str(min(ctx.threads, 8)), out])
        fs = flagstat(runner, ctx, out, d / f"{sid}.dedup.flagstat.json")
        marked_fs = prior_summary(ctx, "markdup")
        expected = marked_fs.get("total_records", 0) - marked_fs.get("duplicates", 0)
        summary = {**fs, "removed_duplicates": marked_fs.get("duplicates")}
        warnings = []
        if fs["total_records"] != expected:
            warnings.append(f"Dedup BAM has {fs['total_records']:,} records, expected {expected:,} "
                            "(supplementary/secondary duplicates are counted differently).")
        return StepResult(outputs={"bam": str(out), "bai": str(out) + ".bai"}, summary=summary,
                          warnings=warnings)

    def validate(self, ctx, result):
        super().validate(ctx, result)
        samtools_quickcheck(Path(result.outputs["bam"]))
