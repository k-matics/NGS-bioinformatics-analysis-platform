"""Germline short-variant workflows (Phase 3).

Registered workflows (all share the same step classes):

  fastq_qc       FASTQ QC & trimming        FastQC raw -> fastp -> FastQC trimmed -> MultiQC
  fastq_to_bam   FASTQ -> analysis-ready BAM  + BWA-MEM/sort -> mark duplicates -> index -> BAM QC
  bam_to_vcf     BAM -> VCF                 prepare BAM -> BAM QC -> callers -> normalise -> report
  fastq_to_vcf   FASTQ -> BAM -> VCF        everything above in one run

Germline only. Somatic (tumour/normal) calling will be separate workflows with
their own callers and inputs - the two are never mixed.
"""
from __future__ import annotations

from ...engine.params import ParamError, ParamSpec, defaults, parse_form, validate_dict
from ...engine.registry import WorkflowDefinition, register
from ..alignment.steps import BwaAlignStep, IndexBamStep, MarkDuplicatesStep, RemoveDuplicatesStep
from ..bam_processing.steps import AnalysisRegionsStep, BamPrepareStep, BamQcStep
from ..common import sum_size_gb
from ..fastq_qc.steps import FastpStep, FastQCStep
from ..intake.steps import FastqIntegrityStep, ReferenceCheckStep
from ..variant_calling.callers import all_callers, get_caller
from ..comparison.steps import CompareCallersStep
from ..output_steps import CleanupStep, RegisterOutputsStep
from ..variant_calling.steps import MultiQCStep, RunReportStep, VariantCallStep, VariantFilterStep

FASTQ_ROLES = ("fastq_r1", "fastq_r2")

# Bump when the scientific behaviour of any step changes (recorded with every run).
# 1.1: samtools stats restricted to the analysis contigs (external BAMs with decoy contigs).
# 1.2: soft filtering per caller, bcftools call -a GQ, DeepVariant, caller comparison,
#      output registration, optional cleanup.
VERSION = "1.2"

# --------------------------------------------------------------------------- parameters
TRIM_PARAMS = [
    ParamSpec("trimming", "Adapter/quality trimming (fastp)", "bool", True,
              "Recommended. Untrimmed adapters cause soft-clipping and false indels near read ends.",
              group="Trimming"),
    ParamSpec("fastp_adapter_trimming", "Adapter trimming", "bool", True,
              "Auto-detect adapters from read-pair overlap (--detect_adapter_for_pe).", group="Trimming"),
    ParamSpec("fastp_min_length", "Minimum read length after trimming", "int", 36,
              "Shorter reads are removed (kept in failed_reads.fastq.gz).", min=15, max=300, group="Trimming"),
    ParamSpec("fastp_qualified_quality", "Qualified base quality", "int", 15,
              "A base is 'low quality' below this Phred score.", min=0, max=40, group="Trimming", advanced=True),
    ParamSpec("fastp_unqualified_percent", "Max % low-quality bases per read", "int", 40,
              "Reads with more low-quality bases are removed.", min=0, max=100, group="Trimming", advanced=True),
    ParamSpec("fastp_n_base_limit", "Max N bases per read", "int", 5, "", min=0, max=50,
              group="Trimming", advanced=True),
]
ALIGN_PARAMS = [
    ParamSpec("library_id", "Library ID (@RG LB)", "text", "",
              "Defaults to the sample ID. Duplicate marking is done per library.", group="Alignment"),
    ParamSpec("markdup_tool", "Duplicate marking tool", "choice", "samtools",
              "samtools markdup (fast) or GATK/Picard MarkDuplicates (GATK Best Practices; slower).",
              choices=[("samtools", "samtools markdup"), ("gatk", "GATK MarkDuplicates (Picard)")],
              group="Alignment"),
    ParamSpec("optical_pixel_distance", "Optical duplicate pixel distance", "int", 2500,
              "2500 for patterned flow cells (NovaSeq/NextSeq 2000); 100 for older non-patterned.",
              min=0, max=10000, group="Alignment", advanced=True),
    ParamSpec("remove_duplicates", "Also write a duplicate-removed BAM", "bool", False,
              "Creates an extra BAM; the duplicate-MARKED BAM is always kept and used for calling.",
              group="Alignment"),
]
OUTPUT_PARAMS = [
    ParamSpec("retain_intermediates", "Keep intermediate files", "bool", True,
              "Off = after a successful run, trimmed FASTQs and caller temp folders are deleted (listed in "
              "reports/cleanup_manifest.json). Removed reads (failed_reads.fastq.gz) are always kept.",
              group="Outputs"),
]
COMPARE_PARAMS = [
    ParamSpec("compare_pass_only", "Compare PASS variants only", "bool", True,
              "Use only records that pass the caller's filters (recommended). Off = all records.",
              group="Comparison"),
    ParamSpec("compare_vcfeval", "Haplotype-aware comparison (RTG vcfeval)", "bool", True,
              "Also compare with RTG vcfeval, which matches equivalent variant representations.",
              group="Comparison", advanced=True),
]
REGION_PARAMS = [
    ParamSpec("contigs", "Contigs to analyse", "text", "primary",
              "'primary' = chr1-22, X, Y, M (ALT/decoy contigs excluded). Or a list, e.g. chr22 or chr1,chr2.",
              group="Regions"),
    ParamSpec("interval_padding", "Target padding (bp)", "int", 100,
              "Added to each side of the target intervals for calling (captures splice-site flanks).",
              min=0, max=1000, group="Regions"),
    ParamSpec("coverage_min_mapq", "Coverage: min mapping quality", "int", 20,
              "Reads below this MAPQ are not counted as coverage (mosdepth -Q).", min=0, max=60,
              group="Regions", advanced=True),
]


def caller_param(enabled_default=("bcftools",)):
    choices = []
    for c in all_callers("germline"):
        ok, why = c.available()
        choices.append((c.key, c.name if ok else f"{c.name} — unavailable: {why}"))
    return ParamSpec("callers", "Variant callers", "multichoice", list(enabled_default),
                     "One output VCF per caller. Each caller has its own parameters below.",
                     choices=choices, group="Variant calling")


def caller_defaults() -> dict:
    return {c.key: defaults(c.params) for c in all_callers("germline")}


def parse_caller_params(form, wf_key: str) -> dict:
    """Form fields named  c.<workflow>.<caller>.<param>  ->  {"caller_params": {caller: {...}}}."""
    selected = form.getlist(f"p.{wf_key}.callers")
    out = {}
    for key in selected:
        try:
            c = get_caller(key)
        except KeyError:
            raise ParamError(f"Unknown caller {key}") from None
        out[key] = parse_form(c.params, form, f"c.{wf_key}.{key}.")
    return {"caller_params": out}


def validate_germline_params(params: dict, reference) -> list[str]:
    problems = []
    callers = params.get("callers")
    if callers is not None:
        if not callers:
            problems.append("Select at least one variant caller.")
        cp = params.setdefault("caller_params", {})
        for key in callers:
            try:
                c = get_caller(key)
            except KeyError:
                problems.append(f"Unknown caller {key}.")
                continue
            ok, why = c.available()
            if not ok:
                problems.append(f"{c.name} is not available: {why}")
            try:
                cp[key] = validate_dict(c.params, cp.get(key, {}), f"{c.name}: ")
            except (ParamError, ValueError, TypeError) as exc:
                problems.append(str(exc))
        for key in list(cp):
            if key not in callers:
                cp.pop(key)   # keep only parameters of callers that will actually run
    contigs = params.get("contigs")
    if contigs and contigs not in ("primary", "all") and reference is not None:
        from ...services.reference_service import reference_contigs
        known = reference_contigs(reference)
        missing = [c for c in contigs.split(",") if c.strip() and c.strip() not in known]
        if missing:
            problems.append(f"Contig(s) not in reference {reference.name}: {missing}")
    return problems


# --------------------------------------------------------------------------- input checks
def check_fastq_inputs(sample, files):
    problems = []
    r1 = [f for f in files if f.role == "fastq_r1"]
    r2 = [f for f in files if f.role == "fastq_r2"]
    if len(r1) != 1 or len(r2) != 1:
        problems.append("Select exactly one FASTQ R1 and one FASTQ R2 (single-end is not supported yet).")
    for f in r1 + r2:
        if f.validation_status == "error":
            problems.append(f"{f.path} failed intake validation: {f.validation_message}")
    if len([f for f in files if f.role == "bed"]) > 1:
        problems.append("Select at most one target BED.")
    return problems


def check_bam_inputs(sample, files):
    problems = []
    bams = [f for f in files if f.role == "bam"]
    if len(bams) != 1:
        problems.append(f"Select exactly one BAM (selected: {len(bams)}).")
    for f in bams:
        if f.validation_status == "error":
            problems.append(f"{f.path} failed intake validation (e.g. different genome build): "
                            f"{f.validation_message}")
    if len([f for f in files if f.role == "bed"]) > 1:
        problems.append("Select at most one target BED.")
    return problems


def _bam_copy_gb(files) -> float:
    """A BAM without read group needs a copy with @RG for GATK (worst case: assume GATK is chosen)."""
    from ...services.hts_inspect import bam_header_summary
    from pathlib import Path
    total = 0.0
    for f in files:
        if f.role == "bam":
            try:
                if not bam_header_summary(Path(f.path))["read_groups"]:
                    total += (f.size_bytes or 0) / 1e9
            except (OSError, ValueError):
                pass
    return total


# --------------------------------------------------------------------------- step lists
def qc_steps():
    return [ReferenceCheckStep(), FastqIntegrityStep(), FastQCStep("raw"), FastpStep(), FastQCStep("trimmed")]


def bam_steps():
    return [BwaAlignStep(), MarkDuplicatesStep(), IndexBamStep(), RemoveDuplicatesStep()]


def call_steps(params):
    steps = []
    for k in params.get("callers", []):
        steps += [VariantCallStep(k), VariantFilterStep(k)]
    return steps + [CompareCallersStep()]


def tail_steps():
    """MultiQC, register outputs, report, then (opt-in) cleanup as the very last step."""
    return [MultiQCStep(), RegisterOutputsStep(), RunReportStep(), CleanupStep()]


TRUTH_ROLES = ("truth_vcf", "truth_bed")


def _register(key, name, description, build, specs, roles, check, disk, module, with_callers):
    default_params = defaults(specs)
    default_params["fastq_check_records"] = 100_000
    if with_callers:
        default_params["callers"] = ["bcftools"]
        default_params["caller_params"] = {"bcftools": defaults(get_caller("bcftools").params)}
    register(WorkflowDefinition(
        key=key, name=name, version=VERSION, description=description, build_steps=build,
        default_params=default_params, check_inputs=check, estimate_disk_gb=disk,
        param_specs=specs, input_roles=roles, module=module,
        parse_extra=parse_caller_params if with_callers else None,
        validate_params=validate_germline_params,
    ))


_register(
    "fastq_qc", "FASTQ QC & trimming",
    "FastQC on raw reads, fastp adapter/quality trimming, FastQC on trimmed reads, MultiQC summary.",
    lambda p: qc_steps() + [MultiQCStep(), RunReportStep()],
    TRIM_PARAMS, FASTQ_ROLES, check_fastq_inputs,
    lambda s, files: 1.1 * sum_size_gb(files, FASTQ_ROLES) + 0.2, "fastq", False,
)
_register(
    "fastq_to_bam", "FASTQ → analysis-ready BAM",
    "QC + trimming, BWA-MEM alignment with read group, sort, duplicate marking, index, BAM QC/coverage.",
    lambda p: qc_steps() + [AnalysisRegionsStep()] + bam_steps() + [BamQcStep()] + tail_steps(),
    TRIM_PARAMS + ALIGN_PARAMS + REGION_PARAMS + OUTPUT_PARAMS, FASTQ_ROLES + ("bed",), check_fastq_inputs,
    lambda s, files: 2.8 * sum_size_gb(files, FASTQ_ROLES) + 0.5, "bam", False,
)
_register(
    "bam_to_vcf", "BAM → VCF (germline)",
    "Validate/prepare an existing BAM, BAM QC, call variants with one or more callers, normalise, report.",
    lambda p: [ReferenceCheckStep(), AnalysisRegionsStep(), BamPrepareStep(), BamQcStep()]
    + call_steps(p) + tail_steps(),
    REGION_PARAMS + COMPARE_PARAMS, ("bam", "bed") + TRUTH_ROLES, check_bam_inputs,
    lambda s, files: 1.0 + _bam_copy_gb(files), "vcf", True,
)
_register(
    "fastq_to_vcf", "FASTQ → BAM → VCF (germline, full workflow)",
    "Complete germline pipeline: QC, trimming, alignment, duplicate marking, BAM QC, multi-caller "
    "variant calling and normalisation, MultiQC + run report.",
    lambda p: qc_steps() + [AnalysisRegionsStep()] + bam_steps() + [BamQcStep()] + call_steps(p)
    + tail_steps(),
    TRIM_PARAMS + ALIGN_PARAMS + REGION_PARAMS + COMPARE_PARAMS + OUTPUT_PARAMS,
    FASTQ_ROLES + ("bed",) + TRUTH_ROLES, check_fastq_inputs,
    lambda s, files: 2.8 * sum_size_gb(files, FASTQ_ROLES) + 4.0, "vcf", True,
)


def form_specs(definition) -> list[ParamSpec]:
    """Workflow param specs plus the (dynamic) caller selector when the workflow calls variants."""
    import dataclasses

    from ...services.settings_service import get_setting
    specs = [dataclasses.replace(p, default=bool(get_setting("retain_intermediates", True)))
             if p.key == "retain_intermediates" else p for p in definition.param_specs]
    if definition.parse_extra is parse_caller_params:
        specs.append(caller_param())
    return specs
