"""Input validation. Everything coming from a browser goes through here
before it can reach the filesystem or a command line."""
from __future__ import annotations

import re
from pathlib import Path

SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
BAM_SUFFIXES = (".bam",)
CRAM_SUFFIXES = (".cram",)
VCF_SUFFIXES = (".vcf.gz", ".vcf", ".bcf")
BED_SUFFIXES = (".bed",)
KNOWN_SUFFIXES = FASTQ_SUFFIXES + BAM_SUFFIXES + CRAM_SUFFIXES + VCF_SUFFIXES + BED_SUFFIXES


class ValidationError(ValueError):
    pass


def validate_sample_id(value: str) -> str:
    value = (value or "").strip()
    if not SAMPLE_ID_RE.match(value):
        raise ValidationError(
            "Sample ID must be 1-64 characters: letters, digits, '_', '-', '.', "
            "starting with a letter or digit."
        )
    return value


def resolve_within_roots(raw_path: str, roots: list[Path]) -> Path:
    """Resolve a user-supplied path (following symlinks) and ensure it lies
    inside one of the allowed data roots. Raises ValidationError otherwise."""
    if not raw_path or "\x00" in raw_path:
        raise ValidationError("Empty or invalid path.")
    path = Path(raw_path).expanduser()
    if not path.is_absolute():
        raise ValidationError("Please give an absolute server path (starting with /).")
    resolved = path.resolve()
    for root in roots:
        if resolved == root or root in resolved.parents:
            return resolved
    raise ValidationError(
        f"Path is outside the allowed data folders ({', '.join(map(str, roots))})."
    )


def guess_role(path: Path) -> str | None:
    """Guess a file role from its name (the user can always override)."""
    name = path.name.lower()
    if name.endswith(FASTQ_SUFFIXES):
        if re.search(r"(_|\.)r?2(_\d{3})?\.f(ast)?q", name):
            return "fastq_r2"
        if re.search(r"(_|\.)r?1(_\d{3})?\.f(ast)?q", name):
            return "fastq_r1"
        return "fastq_se"
    if name.endswith(BAM_SUFFIXES):
        return "bam"
    if name.endswith(CRAM_SUFFIXES):
        return "cram"
    if name.endswith(VCF_SUFFIXES):
        return "vcf"
    if name.endswith(BED_SUFFIXES):
        return "bed"
    return None


def check_role_matches_extension(path: Path, role: str) -> None:
    name = path.name.lower()
    expected = {
        "fastq_r1": FASTQ_SUFFIXES,
        "fastq_r2": FASTQ_SUFFIXES,
        "fastq_se": FASTQ_SUFFIXES,
        "bam": BAM_SUFFIXES,
        "cram": CRAM_SUFFIXES,
        "vcf": VCF_SUFFIXES,
        "bed": BED_SUFFIXES,
        "truth_vcf": VCF_SUFFIXES,
        "truth_bed": BED_SUFFIXES,
    }.get(role)
    if expected is None:
        raise ValidationError(f"Unknown file role '{role}'.")
    if not name.endswith(expected):
        raise ValidationError(f"{path.name} does not look like a {role} file ({', '.join(expected)}).")
