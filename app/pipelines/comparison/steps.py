"""Comparison steps: caller concordance (inside a germline run) and VCF-vs-VCF (standalone)."""
from __future__ import annotations

import itertools
import re
from pathlib import Path

from flask import render_template

from ...engine.step import Step, StepContext, StepFailure, StepResult
from ...models.base import utcnow
from ...services.genome_build import compare_to_reference
from ...services.hts_inspect import vcf_header_summary
from ...services.reference_service import reference_contigs
from ..common import prior_output, write_json
from . import engine
from ...services.limitations import run_limitations

ORIGIN_LABELS = {"internal": "In-house", "external": "External", "pipeline": "Platform pipeline"}


def safe_label(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")[:60] or "set"


def run_comparison(ctx: StepContext, runner, inputs: list[dict], *, out_dir: Path, regions: Path | None,
                   title: str) -> tuple[dict, dict, list[str]]:
    """inputs: [{label, vcf (comparison-ready), source, origin}]. Returns (summary, outputs, warnings)."""
    pass_only = bool(ctx.params.get("compare_pass_only", True))
    warnings: list[str] = []
    sets, skipped = {}, {}
    for inp in inputs:
        try:
            v, sk = engine.load_variants(Path(inp["vcf"]), pass_only=pass_only)
        except ValueError as exc:
            raise StepFailure(str(exc)) from None
        sets[inp["label"]] = v
        skipped[inp["label"]] = dict(sk)
        ctx.log(f"{inp['label']}: {len(v)} variants used; not used: {dict(sk)}")
        ctx.check_cancel()
    result = engine.compare_sets(sets)
    for inp in inputs:
        ps = result["per_set"][inp["label"]]
        ps["not_used"] = skipped[inp["label"]]
        if inp.get("raw"):
            ps["multiallelic_raw"] = engine.count_multiallelic(Path(inp["raw"]))
    files = engine.write_tables(sets, result, out_dir)

    # haplotype-aware pairwise comparison (and truth-set benchmarking when provided)
    vcfeval = {"pairwise": [], "truth": []}
    if ctx.params.get("compare_vcfeval", True) and len(inputs) >= 1:
        sdf = engine.ensure_sdf(runner, ctx.reference, ctx.log)
        t = min(ctx.threads, 8)
        exact = {(p["a"], p["b"]): p["shared"] for p in result["pairwise"]}
        for a, b in itertools.combinations(inputs, 2):
            row = {"a": a["label"], "b": b["label"]}
            for mode, squash in (("genotype", False), ("allele", True)):
                od = out_dir / "vcfeval" / f"{a['label']}__vs__{b['label']}.{mode}"
                r = engine.run_vcfeval(runner, sdf, Path(a["vcf"]), Path(b["vcf"]), od, regions=regions,
                                       threads=t, squash_ploidy=squash, all_records=not pass_only)
                row[mode] = r["overall"]
                ctx.check_cancel()
            shared_exact = exact.get((a["label"], b["label"]), 0)
            if row["allele"]:
                row["exact_shared"] = shared_exact
                row["representation_only"] = row["allele"]["tp_baseline"] - shared_exact
            vcfeval["pairwise"].append(row)
        truth = [f for f in ctx.input_files if f.role == "truth_vcf"]
        truth_bed = [f for f in ctx.input_files if f.role == "truth_bed"]
        if truth:
            tdir = ctx.step_dir("comparison/truth")
            tv = engine.normalise_for_comparison(runner, Path(truth[0].path), Path(ctx.reference.fasta_path),
                                                 tdir / "truth.cmp.vcf.gz", threads=t)
            for inp in inputs:
                od = out_dir / "vcfeval" / f"truth__vs__{inp['label']}"
                r = engine.run_vcfeval(runner, sdf, tv, Path(inp["vcf"]), od, regions=regions,
                                       eval_regions=Path(truth_bed[0].path) if truth_bed else None,
                                       threads=t, all_records=not pass_only)
                vcfeval["truth"].append({"set": inp["label"], **r})
            if not truth_bed:
                warnings.append("Truth VCF given without a confident-regions BED: variants outside the truth "
                                "set's confident regions are counted as false positives.")
    rows = []
    for r in vcfeval["pairwise"]:
        for mode in ("genotype", "allele"):
            if r.get(mode):
                rows.append({"baseline": r["a"], "calls": r["b"], "match": mode, **r[mode]})
    if rows:
        files["vcfeval_csv"] = engine._csv(out_dir / "vcfeval_pairwise.csv", rows)
    if vcfeval["truth"]:
        files["truth_csv"] = engine._csv(out_dir / "truth_benchmark.csv", [
            {"set": r["set"], "class": c, **r[c]} for r in vcfeval["truth"] for c in ("overall", "snp", "non_snp")
            if r.get(c)])

    summary = {
        "title": title, "pass_only": pass_only,
        "inputs": [{k: inp.get(k) for k in ("label", "source", "origin", "vcf")} for inp in inputs],
        **result, "vcfeval": vcfeval,
        "regions": str(regions) if regions else None,
    }
    data = {"summary": summary, "run_id": ctx.run.run_id, "sample": ctx.sample.sample_id,
            "reference": ctx.reference.name, "files": files, "generated_at": utcnow().isoformat() + "Z",
            "origin_labels": ORIGIN_LABELS,
            "limitations": run_limitations(ctx.reference, ctx.params, ctx.sample)}
    html = out_dir / "comparison_report.html"
    html.write_text(render_template("reports/comparison_report.html", **data))
    files["report_html"] = str(html)
    files["summary_json"] = write_json(out_dir / "comparison_summary.json", summary)
    # the DB summary keeps the small tables only
    small = {k: summary[k] for k in ("title", "pass_only", "inputs", "sets", "union", "in_all", "per_set",
                                     "pairwise", "vcfeval")}
    small["intersections"] = result["intersections"][:40]
    return small, files, warnings


# ---------------------------------------------------------------------------
class CompareCallersStep(Step):
    key = "compare_callers"
    name = "Caller comparison (concordance)"
    description = ("Compares the callers of this run on the same BAM: normalised allele matching "
                   "(concordance matrix, intersections, genotype concordance) + RTG vcfeval.")
    tools = ("bcftools", "rtg")

    def applicable(self, ctx):
        n = len(ctx.params.get("callers", []))
        return (n >= 2, "" if n >= 2 else "Only one caller selected — nothing to compare.")

    def execute(self, ctx: StepContext, runner) -> StepResult:
        from ..variant_calling.callers import get_caller
        out = ctx.step_dir("comparison/callers")
        cmp_dir = ctx.step_dir("comparison/callers/vcf")
        ref = Path(ctx.reference.fasta_path)
        inputs = []
        for k in ctx.params["callers"]:
            vcf = Path(prior_output(ctx, f"filter_{k}", "filtered_vcf"))
            cv = engine.normalise_for_comparison(runner, vcf, ref, cmp_dir / f"{k}.cmp.vcf.gz",
                                                 threads=min(ctx.threads, 4))
            inputs.append({"label": k, "vcf": str(cv), "source": str(vcf), "origin": "pipeline",
                           "name": get_caller(k).name,
                           "raw": prior_output(ctx, f"call_{k}", "raw_vcf")})
        regions = Path(prior_output(ctx, "regions", "calling_bed"))
        summary, files, warnings = run_comparison(ctx, runner, inputs, out_dir=out, regions=regions,
                                                  title=f"Caller comparison — {ctx.sample.sample_id}")
        return StepResult(outputs=files, summary=summary, warnings=warnings)


# ---------------------------------------------------------------------------
def _file_label(f, used: set) -> str:
    base = f.label or f"{ORIGIN_LABELS.get(f.origin, f.origin)}_{Path(f.path).name.split('.')[0]}"
    if f.origin == "pipeline" and f.label and " · " in f.label:
        run, what = f.label.split(" · ", 1)            # "WF-…-336ce2 · GATK HaplotypeCaller (filtered)"
        base = f"{run[-6:]}_{what.split(' ')[0]}"         # -> "336ce2_GATK"
    lab = safe_label(base)
    i = 2
    while lab in used:
        lab = f"{safe_label(base)}_{i}"
        i += 1
    used.add(lab)
    return lab


class VcfInputsStep(Step):
    """Checks each selected VCF against the reference and makes it comparable.

    Build guard, two independent checks:
      * ##contig header lines (names + lengths) vs the reference (when present);
      * REF alleles of the first 5,000 records vs the reference sequence. A VCF from
        another build (e.g. GRCh37 on GRCh38) has mostly wrong REF bases -> refused.
    """
    key = "vcf_prepare"
    name = "VCF inputs: build check + normalisation"
    description = ("Genome-build check (contig header + REF alleles), sample selection, then split/"
                   "atomise/left-align against our reference so representations are comparable.")
    tools = ("bcftools",)

    def __init__(self, min_files: int = 2):
        self.min_files = min_files

    def execute(self, ctx: StepContext, runner) -> StepResult:
        files = ctx.files("vcf")
        if len(files) < self.min_files:
            raise StepFailure(f"Select at least {self.min_files} VCF file(s) (selected: {len(files)}).")
        ref = Path(ctx.reference.fasta_path)
        ref_contigs = reference_contigs(ctx.reference)
        d = ctx.step_dir("comparison/vcf_inputs")
        regions = prior_output(ctx, "regions", "calling_bed", required=False)
        used: set = set()
        prepared, warnings, checks = [], [], []
        for f in files:
            src = Path(f.path)
            label = _file_label(f, used)
            hdr = vcf_header_summary(src)
            check = {"label": label, "source": str(src), "origin": f.origin, "samples": hdr["samples"]}
            if hdr["contigs"]:
                cmp = compare_to_reference(hdr["contigs"], ref_contigs)
                check["header_verdict"] = cmp["verdict"]
                if cmp["verdict"] == "incompatible":
                    raise StepFailure(f"{src.name}: {cmp['reason']} Comparing it against {ctx.reference.name} "
                                      "would mix genome builds.")
            else:
                check["header_verdict"] = "no ##contig lines"
                warnings.append(f"{src.name}: no ##contig header lines; the build is checked from REF alleles only.")
            n, bad = engine.ref_mismatch_rate(src, ref)
            check["ref_checked"], check["ref_mismatch"] = n, bad
            if n == 0:
                raise StepFailure(f"{src.name}: none of its records lie on contigs of {ctx.reference.name} "
                                  "(different chromosome naming, e.g. '1' vs 'chr1', or another build).")
            if bad / n > 0.01:
                raise StepFailure(f"{src.name}: {bad} of {n} REF alleles ({100 * bad / n:.1f}%) do not match "
                                  f"{ctx.reference.name} — the VCF is from another genome build or reference. "
                                  "Lift it over or re-call; builds are never mixed.")
            if bad:
                warnings.append(f"{src.name}: {bad}/{n} REF alleles differ from the reference (≤1%, tolerated; "
                                "these records are reported by bcftools norm in the log).")
            sample = None
            if not hdr["samples"]:
                raise StepFailure(f"{src.name} has no sample column (sites-only VCF): genotypes cannot be compared.")
            if len(hdr["samples"]) > 1:
                wanted = ctx.sample.sample_id if ctx.sample.sample_id in hdr["samples"] else hdr["samples"][0]
                sample = wanted
                warnings.append(f"{src.name}: multi-sample VCF; using sample '{wanted}'.")
            out = engine.normalise_for_comparison(runner, src, ref, d / f"{label}.cmp.vcf.gz",
                                                  regions=Path(regions) if regions else None, sample=sample,
                                                  threads=min(ctx.threads, 4))
            prepared.append({"label": label, "vcf": str(out), "source": str(src), "origin": f.origin,
                             "raw": str(src), "sample": sample or (hdr["samples"] or [None])[0]})
            checks.append(check)
        write_json(d / "inputs.json", prepared)
        return StepResult(outputs={"inputs_json": str(d / "inputs.json")},
                          summary={"inputs": checks, "restricted_to": regions}, warnings=warnings)


class VcfCompareStep(Step):
    key = "vcf_compare"
    name = "VCF comparison"
    description = ("Concordance of the selected VCFs: shared / unique variants, genotype and allele "
                   "concordance, representation differences (RTG vcfeval), optional truth-set benchmark.")
    tools = ("rtg",)

    def execute(self, ctx: StepContext, runner) -> StepResult:
        from ..common import read_json
        inputs = read_json(Path(prior_output(ctx, "vcf_prepare", "inputs_json")))
        regions = prior_output(ctx, "regions", "calling_bed", required=False)
        out = ctx.step_dir("comparison/vcf")
        summary, files, warnings = run_comparison(ctx, runner, inputs, out_dir=out,
                                                  regions=Path(regions) if regions else None,
                                                  title=f"VCF comparison — {ctx.sample.sample_id}")
        return StepResult(outputs=files, summary=summary, warnings=warnings)
