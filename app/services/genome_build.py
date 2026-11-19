"""Genome-build detection and contig compatibility.

Why this exists
---------------
A BAM/VCF only makes sense relative to the reference it was aligned to.
GRCh37 (hg19) and GRCh38 (hg38) have different chromosome lengths and
different coordinates for the same gene, so mixing them silently produces
nonsense comparisons. We identify the build from the *lengths* of the main
chromosomes (names alone are not reliable: both builds use "chr1" in UCSC style).
"""
from __future__ import annotations

import re

# Lengths of primary chromosomes that differ between builds.
BUILD_SIGNATURES = {
    "GRCh38": {"1": 248_956_422, "2": 242_193_529, "X": 156_040_895, "22": 50_818_468},
    "GRCh37": {"1": 249_250_621, "2": 243_199_373, "X": 155_270_560, "22": 51_304_566},
}
PRIMARY_RE = re.compile(r"^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$")


def strip_chr(name: str) -> str:
    n = name[3:] if name.startswith("chr") else name
    return "M" if n == "MT" else n


def is_primary(name: str) -> bool:
    return bool(PRIMARY_RE.match(name))


def naming_style(contigs: dict[str, int]) -> str:
    primary = [c for c in contigs if is_primary(c)]
    if not primary:
        return "unknown"
    with_chr = sum(c.startswith("chr") for c in primary)
    if with_chr == len(primary):
        return "ucsc"
    if with_chr == 0:
        return "ensembl"
    return "mixed"


def detect_build(contigs: dict[str, int]) -> str:
    """Return GRCh38 / GRCh37 / unknown from a {contig: length} mapping."""
    norm = {strip_chr(k): v for k, v in contigs.items() if is_primary(k)}
    for build, sig in BUILD_SIGNATURES.items():
        hits = [norm.get(c) == length for c, length in sig.items() if c in norm]
        if hits and all(hits):
            return build
    return "unknown"


def classify_contigs(contigs: dict[str, int]) -> dict[str, int]:
    counts = {"primary": 0, "alt": 0, "random": 0, "unplaced": 0, "decoy": 0, "hla": 0, "other": 0}
    for name in contigs:
        if is_primary(name):
            counts["primary"] += 1
        elif name.endswith("_alt"):
            counts["alt"] += 1
        elif name.endswith("_random"):
            counts["random"] += 1
        # decoys before chrUn: 1000G/DRAGEN decoys are named chrUn_JTFH01000001v1_decoy
        elif "decoy" in name or name == "chrEBV" or name.startswith("hs37d5"):
            counts["decoy"] += 1
        elif name.startswith("chrUn") or name.startswith("GL") or name.startswith("KI"):
            counts["unplaced"] += 1
        elif name.startswith("HLA"):
            counts["hla"] += 1
        else:
            counts["other"] += 1
    return counts


def compare_to_reference(other: dict[str, int], reference: dict[str, int]) -> dict:
    """Compare a file's contigs (BAM @SQ / VCF ##contig) to the reference.

    verdict:
      identical        same contig names and lengths
      compatible       every contig in the file exists in the reference with the same length
      primary_compatible  primary chromosomes match, but the file has extra contigs
                          (decoys/HLA/spike-ins) the reference lacks -> usable, with caveats
      incompatible     a primary chromosome length differs (different build) or naming differs
    """
    other_build = detect_build(other)
    ref_build = detect_build(reference)
    mismatched = sorted(
        c for c in other if c in reference and other[c] != reference[c]
    )
    missing = sorted(c for c in other if c not in reference)
    primary_other = {c: l for c, l in other.items() if is_primary(c)}
    primary_mismatch = [c for c in primary_other if c in reference and reference[c] != primary_other[c]]
    primary_missing = [c for c in primary_other if c not in reference]

    if other_build != "unknown" and ref_build != "unknown" and other_build != ref_build:
        verdict = "incompatible"
        reason = f"Different genome build: file is {other_build}, reference is {ref_build}."
    elif primary_mismatch or (primary_other and len(primary_missing) == len(primary_other)):
        verdict = "incompatible"
        reason = ("Primary chromosome names/lengths do not match the reference "
                  f"(naming: file={naming_style(other)}, reference={naming_style(reference)}).")
    elif not mismatched and not missing and set(other) == set(reference):
        verdict, reason = "identical", "Contigs identical to the reference."
    elif not mismatched and not missing:
        verdict, reason = "compatible", "All contigs present in the reference with identical lengths."
    else:
        verdict = "primary_compatible"
        reason = (f"Primary chromosomes match, but {len(missing)} contig(s) in the file are not in the "
                  "reference (e.g. decoy/HLA/spike-in sequences). Reads on those contigs have no "
                  "equivalent in our reference — expect coverage/variant differences in related regions.")
    return {
        "verdict": verdict,
        "reason": reason,
        "file_build": other_build,
        "reference_build": ref_build,
        "file_naming": naming_style(other),
        "reference_naming": naming_style(reference),
        "file_contigs": len(other),
        "file_contig_classes": classify_contigs(other),
        "length_mismatches": mismatched[:50],
        "missing_in_reference": missing[:50],
        "missing_in_reference_count": len(missing),
    }
