"""Known limitations, shown as warnings on the dashboard, run pages and in every report.

These are deliberate, documented gaps — not hidden. Each item says what is
affected and where it is explained, so results are interpreted correctly.
"""
from __future__ import annotations

from pathlib import Path


def reference_limitations(ref) -> list[dict]:
    out = []
    if ref is None:
        return out
    alt = ref.alt_contig_count or 0
    has_alt_file = Path(str(ref.fasta_path) + ".alt").is_file()
    rep = ref.check_report or {}
    masked = bool((rep.get("alt_masking") or {}).get("verified"))
    if alt and masked:
        out.append({
            "key": "alt_masked", "severity": "info",
            "title": "ALT contigs masked (DRAGEN-style)",
            "text": (f"{ref.name}: the ALT regions listed in the DRAGEN mask BED are replaced by N, so reads map "
                     "to the primary chromosomes. This matches DRAGEN's default hg38 reference; it is not BWA's "
                     "ALT-aware mode, and reads are not assigned to ALT haplotypes."),
            "doc": "BAM_PROCESSING.md §5",
        })
    elif alt and not has_alt_file:
        out.append({
            "key": "non_alt_aware", "severity": "warning",
            "title": "Reference is not ALT-aware",
            "text": (f"{ref.name} contains {alt} ALT contigs but no BWA .alt file. Reads from regions with ALT "
                     "copies (e.g. the MHC/HLA region on chr6, parts of chr15-17, 19, 22) get MAPQ 0 and are "
                     "ignored by coverage and callers. SAMPLE_A: 4.6% of target bases had no usable coverage. "
                     "Fix planned: ALT-masked / no-ALT analysis-set reference."),
            "doc": "BAM_PROCESSING.md §5",
        })
    if (ref.build or "").startswith("GRCh38") and alt and rep.get("par_masked") is not True:
        out.append({
            "key": "par_unmasked", "severity": "warning",
            "title": "chrY PAR not masked",
            "text": ("The pseudo-autosomal regions exist on both chrX and chrY in this reference, so PAR reads "
                     "get MAPQ 0 (part of the chrX coverage gap)."),
            "doc": "BAM_PROCESSING.md §5",
        })
    return out


STATIC = [
    {"key": "ploidy", "severity": "info", "title": "Sex chromosomes called as diploid",
     "text": ("All callers assume two copies of chrX/chrY. For a male sample, non-PAR chrX/chrY genotypes should "
              "be haploid; no sex check is performed yet."), "doc": "VARIANT_CALLING.md §8"},
    {"key": "no_bqsr", "severity": "info", "title": "No base-quality recalibration (BQSR)",
     "text": "Optional in current GATK practice; not applied.", "doc": "VARIANT_CALLING.md §8"},
    {"key": "difficult_genes", "severity": "info", "title": "Paralogous / pseudogene regions",
     "text": ("Genes with near-identical copies (PMS2, SMN1/2, CYP21A2, GBA, CYP2D6, …) are unreliable with "
              "short-read alignment and need dedicated methods."), "doc": "VARIANT_CALLING.md §3"},
    {"key": "not_validated", "severity": "danger", "title": "Not clinically validated",
     "text": ("Research use only. No truth-set benchmark of this pipeline has been performed on this server; "
              "producing a VCF is not evidence of analytical accuracy."), "doc": "VARIANT_CALLER_COMPARISON.md"},
]


def run_limitations(ref, params: dict | None = None, sample=None) -> list[dict]:
    items = reference_limitations(ref) + [dict(x) for x in STATIC]
    params = params or {}
    if "deepvariant" in (params.get("callers") or []):
        items.append({"key": "dv_model", "severity": "info", "title": "DeepVariant model",
                      "text": (f"Model {((params.get('caller_params') or {}).get('deepvariant') or {}).get('model_type', 'WES')} "
                               "was trained on specific library types; performance on other capture kits is "
                               "not established here."), "doc": "VARIANT_CALLING.md §2"})
    return items


def platform_limitations() -> list[dict]:
    from .reference_service import get_default_reference
    try:
        ref = get_default_reference()
    except Exception:  # noqa: BLE001 - dashboard must render even if DB is not initialised
        ref = None
    return reference_limitations(ref) + [dict(x) for x in STATIC]
