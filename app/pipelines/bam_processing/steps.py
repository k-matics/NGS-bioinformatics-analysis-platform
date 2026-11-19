"""Analysis regions, input-BAM preparation and BAM quality control."""
from __future__ import annotations

import os
from pathlib import Path

from ...engine.step import Step, StepContext, StepFailure, StepResult
from ...services.genome_build import compare_to_reference
from ...services.hts_inspect import bam_header_summary
from ...services.reference_service import reference_contigs
from ..common import (
    bed_size,
    check_index_fresh,
    merge_intervals,
    parse_contig_selection,
    parse_flagstat_json,
    parse_mosdepth,
    prior_output,
    read_bed,
    samtools_quickcheck,
    write_bed,
    write_json,
)

MOSDEPTH_THRESHOLDS = [1, 10, 20, 30, 50, 100]
DUP_PROGRAMS = ("markdup", "markduplicates", "dragen", "sambamba", "biobambam", "bammarkduplicates")


def analysis_bam(ctx: StepContext) -> Path:
    """The BAM that QC and callers should use in this run."""
    keys = {s.key for s in ctx.run.steps}
    if "index_bam" in keys:
        return Path(prior_output(ctx, "index_bam", "bam"))
    return Path(prior_output(ctx, "bam_prepare", "bam"))


# ---------------------------------------------------------------------------
class AnalysisRegionsStep(Step):
    key = "regions"
    name = "Analysis regions"
    description = ("Builds the target BED (for coverage) and the padded calling-region BED, restricted to "
                   "the selected contigs (default: primary chromosomes, excluding ALT/decoy contigs).")

    def execute(self, ctx: StepContext, runner) -> StepResult:
        ref_contigs = reference_contigs(ctx.reference)
        if not ref_contigs:
            raise StepFailure("Reference .fai missing or empty.")
        contigs = parse_contig_selection(ctx.params.get("contigs", "primary"), ref_contigs)
        padding = int(ctx.params.get("interval_padding", 100))
        d = ctx.step_dir("regions")
        warnings, summary, outputs = [], {"contigs": contigs if len(contigs) <= 30 else f"{len(contigs)} contigs"}, {}

        beds = ctx.files("bed")
        if len(beds) > 1:
            raise StepFailure("Select at most one target BED file for a run.")
        if beds:
            src = Path(beds[0].path)
            raw = read_bed(src)
            not_in_ref = sorted({c for c, _, _ in raw if c not in ref_contigs})
            if not_in_ref and len(not_in_ref) == len({c for c, _, _ in raw}):
                raise StepFailure(f"No contig of {src.name} exists in the reference "
                                  f"(BED uses e.g. {not_in_ref[0]!r}). Wrong build or naming (chr1 vs 1)?")
            if not_in_ref:
                warnings.append(f"{src.name}: {len(not_in_ref)} contig(s) not in the reference ignored: "
                                f"{not_in_ref[:5]}")
            selected = set(contigs)
            raw = [iv for iv in raw if iv[0] in selected]
            if not raw:
                raise StepFailure(f"The target BED has no intervals on the selected contigs {contigs[:5]}.")
            targets = merge_intervals(raw, ref_contigs, 0)
            calling = merge_intervals(raw, ref_contigs, padding)
            outputs["targets_bed"] = write_bed(d / "targets.bed", targets)
            summary.update({"target_bed_source": str(src), "target_intervals": len(targets),
                            "target_bp": bed_size(targets), "padding": padding,
                            "calling_intervals": len(calling), "calling_bp": bed_size(calling)})
        else:
            calling = [(c, 0, ref_contigs[c]) for c in contigs]
            summary.update({"target_bed_source": None, "calling_bp": bed_size(calling)})
            if ctx.sample.sequencing_type != "WGS":
                warnings.append(f"No target BED selected for a {ctx.sample.sequencing_type} sample: coverage "
                                "and variant calling cover whole chromosomes (slower; off-target calls "
                                "will be included).")
        outputs["calling_bed"] = write_bed(d / "calling_regions.bed", calling)
        write_json(d / "regions_summary.json", summary)
        return StepResult(outputs=outputs, summary=summary, warnings=warnings)


def prepare_bam(ctx: StepContext, runner, f, d: Path, need_rg: bool = False):
    """Integrity, build, sort-order and duplicate checks of a registered BAM; links it into ``d``
    and creates an index (and, if needed, a read-group copy) there. The input is never modified.
    Returns (bam, bai, summary, warnings)."""
    src = Path(f.path)
    runner.run(["samtools", "quickcheck", "-v", src])
    hdr = bam_header_summary(src)
    cmp = compare_to_reference(hdr["contigs"], reference_contigs(ctx.reference))
    warnings = []
    if cmp["verdict"] == "incompatible":
        raise StepFailure(f"{src.name}: {cmp['reason']} Using it with this reference would mix "
                          "genome builds. Re-align the FASTQ to this reference instead.")
    if cmp["verdict"] == "primary_compatible":
        warnings.append(f"{src.name}: {cmp['reason']} GATK may refuse BAMs whose contigs differ "
                        "from the reference dictionary.")
    if hdr["sort_order"] != "coordinate":
        raise StepFailure(f"{src.name} is not coordinate-sorted (SO:{hdr['sort_order']}). Sort it first.")
    progs = " ".join(f"{p['id']} {p['name']} {p['command'][:100]}" for p in hdr["programs"]).lower()
    dups_marked = any(t in progs for t in DUP_PROGRAMS)
    if not dups_marked:
        warnings.append(f"{src.name}: no duplicate-marking program in the header (@PG). PCR duplicates "
                        "will be counted as independent evidence → inflated confidence/false positives.")
    link = d / src.name
    if not link.exists():
        os.symlink(src, link)
    bam = link
    rg_added = False
    if not hdr["read_groups"]:
        if need_rg:
            sid = ctx.sample.sample_id
            bam = d / f"{sid}.with_rg.bam"
            ctx.log("No @RG in the BAM; GATK requires one. Writing a copy with a read group into the run "
                    "folder (the input BAM is not modified).")
            runner.run(["samtools", "addreplacerg", "-@", str(min(ctx.threads, 8)),
                        "-r", f"ID:{sid}", "-r", f"SM:{sid}", "-r", f"LB:{sid}", "-r", "PL:ILLUMINA",
                        "-o", bam, link])
            rg_added = True
            warnings.append(f"{src.name}: had no read group; a copy with RG ID/SM={sid} was created.")
        else:
            warnings.append(f"{src.name}: no read group; VCF sample name will be generic.")
    elif hdr["samples"] and ctx.sample.sample_id not in hdr["samples"]:
        warnings.append(f"BAM sample name(s) {hdr['samples']} differ from sample ID "
                        f"{ctx.sample.sample_id}; the VCF will use the BAM's name.")
    bai = Path(str(bam) + ".bai")
    if not rg_added and hdr["index"] and hdr["index"].endswith(".bai") and not bai.exists():
        os.symlink(hdr["index"], bai)
        if bai.stat().st_mtime < src.stat().st_mtime:
            warnings.append("The existing BAM index is older than the BAM — re-indexing in the run folder.")
            bai.unlink()
    if not bai.exists():
        runner.run(["samtools", "index", "-@", str(min(ctx.threads, 8)), bam])
    summary = {"source": str(src), "verdict": cmp["verdict"], "build": cmp["file_build"],
               "read_groups": hdr["read_groups"], "read_group_added": rg_added,
               "duplicates_marked": dups_marked, "indexed_in_run_folder": not hdr["index"] or rg_added,
               "size_bytes": src.stat().st_size, "programs": [p["name"] for p in hdr["programs"]][:10]}
    return bam, bai, summary, warnings


# ---------------------------------------------------------------------------
class BamPrepareStep(Step):
    key = "bam_prepare"
    name = "Input BAM preparation"
    description = ("Integrity + genome-build check of the registered BAM, links it into the run folder, "
                   "creates an index / read group in the run folder if missing (input is never modified).")
    tools = ("samtools",)

    def execute(self, ctx: StepContext, runner) -> StepResult:
        from ..variant_calling.callers import get_caller

        bams = ctx.files("bam")
        if len(bams) != 1:
            raise StepFailure("Select exactly one BAM file for BAM → VCF.")
        callers = [get_caller(k) for k in ctx.params.get("callers", [])]
        need_rg = any(c.requires_read_group for c in callers)
        bam, bai, summary, warnings = prepare_bam(ctx, runner, bams[0], ctx.step_dir("bam"), need_rg)
        return StepResult(outputs={"bam": str(bam), "bai": str(bai)}, summary=summary, warnings=warnings)

    def validate(self, ctx, result):
        bam = Path(result.outputs.get("bam", ""))
        if not bam.exists():
            raise StepFailure(f"Prepared BAM missing: {bam}")
        samtools_quickcheck(bam)
        if not Path(result.outputs["bai"]).exists():
            raise StepFailure("BAM index missing.")


# ---------------------------------------------------------------------------
class BamQcStep(Step):
    key = "bam_qc"
    name = "BAM QC (flagstat, stats, coverage)"
    description = ("samtools flagstat + stats (mapping, pairing, insert size, error rate) and mosdepth "
                   "coverage over the targets (duplicates and QC-fail reads excluded).")
    tools = ("samtools", "mosdepth")

    def execute(self, ctx: StepContext, runner) -> StepResult:
        bam = analysis_bam(ctx)
        sid = ctx.sample.sample_id
        d = ctx.step_dir("qc/bam")
        t = str(min(ctx.threads, 4))
        fs_path = d / f"{sid}.flagstat.json"
        runner.run(["samtools", "flagstat", "-@", t, "-O", "json", bam], stdout_path=fs_path)
        stats_path = d / f"{sid}.samtools_stats.txt"
        # Restrict stats to the analysis contigs: (1) external BAMs may contain contigs our reference
        # lacks (e.g. DRAGEN decoys) -> "-r" cannot fetch their sequence; (2) QC of different BAMs is
        # then computed over the same chromosomes and is comparable.
        ref_contigs = reference_contigs(ctx.reference)
        contigs = parse_contig_selection(ctx.params.get("contigs", "primary"), ref_contigs)
        runner.run(["samtools", "stats", "-@", t, "-r", ctx.reference.fasta_path, bam, *contigs],
                   stdout_path=stats_path)
        targets = prior_output(ctx, "regions", "targets_bed", required=False)
        by = targets or prior_output(ctx, "regions", "calling_bed")
        prefix = d / sid
        runner.run(["mosdepth", "-t", t, "-n", "-x", "-Q", str(ctx.params.get("coverage_min_mapq", 20)),
                    "--by", by, "-T", ",".join(map(str, MOSDEPTH_THRESHOLDS)), prefix, bam])
        fs = parse_flagstat_json(fs_path)
        cov = parse_mosdepth(prefix, MOSDEPTH_THRESHOLDS)
        sn = {}
        with open(stats_path) as fh:
            for line in fh:
                if line.startswith("SN\t"):
                    k, v = line.split("\t")[1:3]
                    sn[k.rstrip(":")] = v
        summary = {
            **fs,
            "insert_size_mean": float(sn.get("insert size average", 0)),
            "insert_size_sd": float(sn.get("insert size standard deviation", 0)),
            "error_rate": float(sn.get("error rate", 0)),
            "coverage_region": "targets" if targets else "calling regions",
            "stats_scope": f"{len(contigs)} analysis contig(s)" if len(contigs) > 3 else ", ".join(contigs),
            **cov,
        }
        warnings = []
        mean = cov.get("mean_target_coverage")
        pct20 = (cov.get("pct_target_bases") or {}).get("20x")
        if mean is not None and mean < 30:
            warnings.append(f"Mean coverage over {summary['coverage_region']} is {mean:.1f}x (< 30x).")
        if pct20 is not None and pct20 < 90:
            warnings.append(f"Only {pct20}% of {summary['coverage_region']} bases reach 20x (< 90%).")
        write_json(d / "bam_qc_summary.json", summary)
        return StepResult(
            outputs={"flagstat": str(fs_path), "stats": str(stats_path),
                     "mosdepth_summary": f"{prefix}.mosdepth.summary.txt",
                     "mosdepth_regions": f"{prefix}.regions.bed.gz"},
            summary=summary, warnings=warnings)

    def validate(self, ctx, result):
        super().validate(ctx, result)
        check_index_fresh(analysis_bam(ctx))
