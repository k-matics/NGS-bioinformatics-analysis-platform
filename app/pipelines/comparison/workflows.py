"""Standalone comparison workflows (Caller Comparison / External Data Comparison modules)."""
from __future__ import annotations

from ...engine.params import ParamSpec, defaults
from ...engine.registry import WorkflowDefinition, register
from ..bam_processing.steps import AnalysisRegionsStep, BamPrepareStep
from ..germline import COMPARE_PARAMS, REGION_PARAMS, validate_germline_params
from ..intake.steps import ReferenceCheckStep
from ..variant_calling.steps import RunReportStep
from .bam_steps import BamCompareStep, BamMetricsStep
from .evidence import ReadEvidenceStep
from .steps import VcfCompareStep, VcfInputsStep

VERSION = "1.0"
REGIONS = [p for p in REGION_PARAMS if p.key in ("contigs", "interval_padding")]


def _count(files, role):
    return len([f for f in files if f.role == role])


def _bad_build(files, roles):
    return [f"{f.path} failed intake validation (e.g. different genome build): {f.validation_message}"
            for f in files if f.role in roles and f.validation_status == "error"]


def check_vcf_compare(sample, files):
    p = _bad_build(files, ("vcf",))
    if _count(files, "vcf") < 2:
        p.append("Select at least 2 VCF files.")
    if _count(files, "vcf") > 8:
        p.append("Select at most 8 VCF files.")
    for r in ("bed", "truth_vcf", "truth_bed"):
        if _count(files, r) > 1:
            p.append(f"Select at most one {r} file.")
    return p


def check_bam_compare(sample, files):
    p = _bad_build(files, ("bam",))
    if not 2 <= _count(files, "bam") <= 6:
        p.append("Select 2-6 BAM files.")
    if _count(files, "bed") > 1:
        p.append("Select at most one target BED.")
    return p


def check_evidence(sample, files):
    p = _bad_build(files, ("bam", "vcf"))
    if _count(files, "bam") != 1 or _count(files, "vcf") != 1:
        p.append("Select exactly one BAM and one VCF.")
    if _count(files, "bed") > 1:
        p.append("Select at most one BED.")
    return p


BAM_COMPARE_PARAMS = [
    ParamSpec("coverage_thresholds", "Coverage thresholds (x)", "text", "1,10,20,30,50,100",
              "Comma-separated depths for '% target bases ≥ N x'.", group="Coverage"),
    ParamSpec("coverage_min_mapq", "Coverage: min mapping quality", "int", 20,
              "Reads below this MAPQ are not counted as coverage.", min=0, max=60, group="Coverage"),
    ParamSpec("target_diff_threshold", "Per-target difference threshold (x)", "int", 20,
              "List targets whose mean depth reaches this in one BAM but not in the other.",
              min=1, max=1000, group="Coverage"),
]
EVIDENCE_PARAMS = [
    ParamSpec("evidence_pass_only", "PASS variants only", "bool", False,
              "Off = every record of the VCF is checked (filtered ones too).", group="Evidence"),
    ParamSpec("evidence_min_bq", "Min base quality (SNV bases)", "int", 13,
              "ALT/REF bases below this quality are not counted (as bcftools mpileup default).",
              min=0, max=60, group="Evidence"),
    ParamSpec("evidence_max_sites", "Max variants", "int", 200000,
              "Safety limit; the run stops (no silent sampling) if the VCF has more.",
              min=1, max=5000000, group="Evidence", advanced=True),
]


def _reg(key, name, description, build, specs, roles, check, module, disk=0.5):
    register(WorkflowDefinition(
        key=key, name=name, version=VERSION, description=description, build_steps=build,
        default_params=defaults(specs), check_inputs=check, estimate_disk_gb=lambda s, f: disk,
        param_specs=specs, input_roles=roles, module=module, validate_params=validate_germline_params))


_reg("vcf_compare", "VCF vs VCF comparison",
     "Compare 2-8 VCFs of one sample (callers, in-house vs external provider): build check, normalisation, "
     "concordance matrix, genotype concordance, RTG vcfeval; optional truth-set benchmark (GIAB).",
     lambda p: [ReferenceCheckStep(), AnalysisRegionsStep(), VcfInputsStep(2), VcfCompareStep(), RunReportStep()],
     REGIONS + COMPARE_PARAMS, ("vcf", "bed", "truth_vcf", "truth_bed"), check_vcf_compare, "comparison")

_reg("bam_compare", "BAM vs BAM comparison",
     "Side-by-side QC of 2-6 BAMs over the same targets: read categories, MAPQ distribution, reads per contig "
     "class, insert size, error rate, target coverage, targets covered differently.",
     lambda p: [ReferenceCheckStep(), AnalysisRegionsStep(), BamMetricsStep(), BamCompareStep(), RunReportStep()],
     [p for p in REGIONS if p.key == "contigs"] + BAM_COMPARE_PARAMS, ("bam", "bed"), check_bam_compare,
     "external", disk=1.0)

_reg("bam_vcf_evidence", "BAM vs VCF read-level validation",
     "For every variant of a VCF, measure the supporting reads in a BAM (REF/ALT counts, VAF, strand, MAPQ, "
     "soft-clips, read-end position) and flag variants whose evidence needs review; IGV batch file.",
     lambda p: [ReferenceCheckStep(), AnalysisRegionsStep(), BamPrepareStep(), VcfInputsStep(1),
                ReadEvidenceStep(), RunReportStep()],
     REGIONS + EVIDENCE_PARAMS, ("bam", "vcf", "bed"), check_evidence, "external")


# ------------------------------------------------------------------ variant-by-variant check (new, 2026-11-19)
# A separate workflow: the existing vcf_compare above is unchanged.
from .variant_check import VariantCheckStep  # noqa: E402

VARIANT_CHECK_PARAMS = [
    ParamSpec("baseline_pass_only", "Baseline: PASS variants only", "bool", True,
              "Check only the baseline's PASS variants (recommended). Off = every baseline record.",
              group="Variant check"),
    ParamSpec("query_all_records", "Our VCF: include filtered records", "bool", True,
              "On = a variant counts as present even if our filters labelled it (shown as 'found, filtered', "
              "with the filter name). Lets you check a permissive VCF and filter afterwards in the viewer.",
              group="Variant check"),
    ParamSpec("nearby_window", "Nearby window (bp)", "int", 10,
              "For variants not found: list our variants within this distance (different allele, shifted indel).",
              min=0, max=1000, group="Variant check"),
    ParamSpec("check_haplotype", "Rescue representation differences (RTG vcfeval)", "bool", True,
              "Variants written differently but meaning the same sequence count as found.",
              group="Variant check", advanced=True),
    ParamSpec("report_query_only", "Also list our PASS variants missing from the baseline", "bool", True,
              "", group="Variant check", advanced=True),
]


def check_variant_check(sample, files):
    p = _bad_build(files, ("vcf",))
    if _count(files, "vcf") != 2:
        p.append("Select exactly two VCF files: the baseline (e.g. the provider's) and ours.")
    if _count(files, "bed") > 1:
        p.append("Select at most one BED.")
    return p


def parse_variant_check(form, wf_key):
    fid = form.get(f"x.{wf_key}.baseline_file_id", type=int)
    return {"baseline_file_id": fid} if fid else {}


register(WorkflowDefinition(
    key="vcf_variant_check", name="VCF variant-by-variant check (baseline → ours)", version=VERSION,
    description=("For every variant of a baseline VCF (e.g. the provider's), look it up in our VCF: found and "
                 "PASS, found but filtered (by which filter), genotyped as reference, written differently, a "
                 "different variant nearby, or not called. Side-by-side QUAL/DP/AD/VAF/GQ; interactive viewer "
                 "with threshold what-ifs."),
    build_steps=lambda p: [ReferenceCheckStep(), AnalysisRegionsStep(), VcfInputsStep(2), VariantCheckStep(),
                           RunReportStep()],
    default_params=defaults(REGIONS + VARIANT_CHECK_PARAMS), check_inputs=check_variant_check,
    estimate_disk_gb=lambda s, f: 0.5, param_specs=REGIONS + VARIANT_CHECK_PARAMS, input_roles=("vcf", "bed"),
    module="comparison", parse_extra=parse_variant_check, validate_params=validate_germline_params))
