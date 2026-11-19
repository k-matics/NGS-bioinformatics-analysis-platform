"""Reference genome registry and integrity checks.

The reference FASTA is *never* modified. We only read the FASTA header line
and its index files.
"""
from __future__ import annotations

import logging
from pathlib import Path

from flask import current_app

from ..extensions import db
from ..models import ReferenceGenome
from ..models.base import utcnow
from .genome_build import classify_contigs, detect_build, naming_style

log = logging.getLogger(__name__)

BWA_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")
BWA_MEM2_SUFFIXES = (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")


def read_fai(fai_path: Path) -> dict[str, int]:
    contigs = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                contigs[parts[0]] = int(parts[1])
    return contigs


def read_dict(dict_path: Path) -> dict[str, int]:
    contigs = {}
    with open(dict_path) as fh:
        for line in fh:
            if not line.startswith("@SQ"):
                continue
            fields = dict(f.split(":", 1) for f in line.rstrip("\n").split("\t")[1:] if ":" in f)
            contigs[fields["SN"]] = int(fields["LN"])
    return contigs


def inspect_reference(fasta: Path) -> dict:
    """Check a reference FASTA and its companion files. Returns a report dict
    with 'status' = ok / warning / error and lists of errors/warnings."""
    fasta = Path(fasta).expanduser()
    report: dict = {"fasta": str(fasta), "errors": [], "warnings": [], "files": {}}

    if not fasta.is_file():
        report["errors"].append(f"Reference FASTA not found: {fasta}")
        report["status"] = "error"
        return report
    report["files"]["fasta"] = {"exists": True, "size_bytes": fasta.stat().st_size}

    if fasta.suffix == ".gz":
        gzi = Path(str(fasta) + ".gzi")
        if not gzi.exists():
            report["errors"].append("Gzipped FASTA must be bgzip-compressed with a .gzi index.")

    with open(fasta, "rb") as fh:
        first = fh.readline(200)
    if not first.startswith(b">"):
        report["errors"].append("File does not start with a FASTA header ('>').")

    fai = Path(str(fasta) + ".fai")
    contigs: dict[str, int] = {}
    if fai.is_file():
        contigs = read_fai(fai)
        report["files"]["fai"] = {"exists": True, "contigs": len(contigs)}
        if fai.stat().st_mtime < fasta.stat().st_mtime:
            report["warnings"].append(".fai is older than the FASTA — it may be stale.")
    else:
        report["files"]["fai"] = {"exists": False}
        report["errors"].append("FASTA index (.fai) missing — create with: samtools faidx <fasta>")

    dict_path = fasta.with_suffix(".dict") if fasta.suffix != ".gz" else Path(str(fasta)[:-3]).with_suffix(".dict")
    if dict_path.is_file():
        dict_contigs = read_dict(dict_path)
        same = dict_contigs == contigs
        report["files"]["dict"] = {"exists": True, "contigs": len(dict_contigs), "matches_fai": same}
        if contigs and not same:
            report["errors"].append(".dict contigs/lengths do not match the .fai.")
    else:
        report["files"]["dict"] = {"exists": False}
        report["warnings"].append("Sequence dictionary (.dict) missing — GATK needs it "
                                  "(gatk CreateSequenceDictionary -R <fasta>).")

    bwa_ok = all(Path(str(fasta) + s).is_file() for s in BWA_SUFFIXES)
    report["files"]["bwa_index"] = {"exists": bwa_ok}
    if bwa_ok:
        ann = Path(str(fasta) + ".ann")
        with open(ann) as fh:
            n_seqs = int(fh.readline().split()[1])
        report["files"]["bwa_index"]["contigs"] = n_seqs
        if contigs and n_seqs != len(contigs):
            report["errors"].append(f"BWA index has {n_seqs} contigs but .fai has {len(contigs)} — "
                                    "index built from a different FASTA.")
    else:
        report["warnings"].append("BWA index missing — alignment with BWA will not be possible.")
    report["files"]["bwa_mem2_index"] = {"exists": all(Path(str(fasta) + s).is_file() for s in BWA_MEM2_SUFFIXES)}

    alt_file = Path(str(fasta) + ".alt")
    report["files"]["alt_file"] = {"exists": alt_file.is_file()}

    # DRAGEN-style ALT masking: a sidecar <fasta>.masking.json written by
    # scripts/build_masked_reference.sh describes which regions were replaced by N.
    # We do not trust the sidecar alone: a sample of the masked intervals is read
    # back from the FASTA and must be all N.
    masking = None
    if contigs and fasta.suffix != ".gz":
        masking = verify_masking(fasta, contigs)
        if masking is not None:
            report["alt_masking"] = masking
            if not masking["verified"]:
                report["errors"].append("ALT-masking sidecar present but the FASTA is not masked as described: "
                                        + "; ".join(masking["problems"][:3]))
        par = check_par_masked(fasta, contigs)
        if par is not None:
            report["par_masked"] = par

    if contigs:
        classes = classify_contigs(contigs)
        report["build"] = detect_build(contigs)
        report["naming"] = naming_style(contigs)
        report["contig_count"] = len(contigs)
        report["contig_classes"] = classes
        if report["build"] == "unknown":
            report["warnings"].append("Could not recognise the genome build from chromosome lengths.")
        if classes["alt"] and masking and masking["verified"]:
            pass  # ALT copies masked DRAGEN-style: reads map to the primary assembly
        elif classes["alt"] and not alt_file.is_file():
            report["warnings"].append(
                f"Reference contains {classes['alt']} ALT contigs but no BWA .alt file: alignment is "
                "NOT ALT-aware. Reads in regions with ALT haplotypes (e.g. HLA, KIR, some clinically "
                "relevant genes) may receive MAPQ 0 and be ignored by variant callers."
            )
    report["status"] = "error" if report["errors"] else ("warning" if report["warnings"] else "ok")
    return report


def _read_region(fasta: Path, fai: dict[str, tuple[int, int, int, int]], contig: str,
                 start: int, end: int) -> bytes:
    """Read 0-based [start, end) of a contig straight from an uncompressed FASTA via its .fai."""
    length, offset, line_bases, line_bytes = fai[contig]
    end = min(end, length)
    if start >= end:
        return b""
    first = offset + (start // line_bases) * line_bytes + start % line_bases
    last = offset + ((end - 1) // line_bases) * line_bytes + (end - 1) % line_bases
    with open(fasta, "rb") as fh:
        fh.seek(first)
        raw = fh.read(last - first + 1)
    return raw.replace(b"\n", b"").replace(b"\r", b"")


def _read_fai_full(fai_path: Path) -> dict[str, tuple[int, int, int, int]]:
    out = {}
    with open(fai_path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                out[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
    return out


def verify_masking(fasta: Path, contigs: dict[str, int], sample: int = 60) -> dict | None:
    """Return None if the FASTA has no masking sidecar, else a verification dict."""
    import json
    side = Path(str(fasta) + ".masking.json")
    if not side.is_file():
        return None
    info = json.loads(side.read_text())
    out = {"sidecar": str(side), "method": info.get("method"), "mask_bed": info.get("mask_bed"),
           "masked_bp": info.get("masked_bp"), "intervals": info.get("intervals"),
           "source_fasta": info.get("source_fasta"), "problems": [], "checked": 0}
    bed = Path(info.get("mask_bed_copy") or "")
    if not bed.is_file():
        out["problems"].append(f"mask BED copy not found: {bed}")
        out["verified"] = False
        return out
    ivs = []
    with open(bed) as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) >= 3 and not line.startswith(("#", "track")):
                ivs.append((p[0], int(p[1]), int(p[2])))
    fai = _read_fai_full(Path(str(fasta) + ".fai"))
    step = max(1, len(ivs) // sample)
    for chrom, s, e in ivs[::step]:
        if chrom not in fai:
            out["problems"].append(f"{chrom} not in FASTA")
            continue
        # check the first and last 2 kb of each sampled interval (fast, and catches off-by-one shifts)
        for a, b in ((s, min(e, s + 2000)), (max(s, e - 2000), e)):
            seq = _read_region(fasta, fai, chrom, a, b).upper()
            if seq.strip(b"N"):
                out["problems"].append(f"{chrom}:{a + 1}-{b} is not all N")
        out["checked"] += 1
    out["verified"] = not out["problems"]
    return out


# chrY pseudo-autosomal regions (GRCh38, 1-based inclusive)
PAR_Y_GRCH38 = [("chrY", 10_001, 2_781_479), ("chrY", 56_887_903, 57_217_415)]


def check_par_masked(fasta: Path, contigs: dict[str, int]) -> bool | None:
    """True if the chrY PARs are hard-masked (N), as in analysis-set references."""
    if contigs.get("chrY") != 57_227_415:
        return None
    fai = _read_fai_full(Path(str(fasta) + ".fai"))
    for chrom, s, e in PAR_Y_GRCH38:
        mid = (s + e) // 2
        for a in (s - 1, mid, e - 1000):
            if _read_region(fasta, fai, chrom, a, a + 1000).upper().strip(b"N"):
                return False
    return True


def apply_report(ref: ReferenceGenome, report: dict) -> None:
    ref.check_report = report
    ref.status = report["status"]
    ref.build = report.get("build", "unknown")
    ref.naming = report.get("naming")
    ref.contig_count = report.get("contig_count")
    ref.alt_contig_count = (report.get("contig_classes") or {}).get("alt")
    ref.last_checked = utcnow()


def ensure_default_reference() -> None:
    """Register the configured default reference on first start."""
    if ReferenceGenome.query.count():
        return
    fasta = current_app.config["DEFAULT_REFERENCE_FASTA"]
    ref = ReferenceGenome(name=current_app.config["DEFAULT_REFERENCE_NAME"],
                          fasta_path=str(fasta), is_default=True)
    try:
        apply_report(ref, inspect_reference(fasta))
    except Exception as exc:  # never block app start-up on this
        log.exception("Reference inspection failed")
        ref.status, ref.notes = "error", str(exc)
    db.session.add(ref)
    db.session.commit()


def get_default_reference() -> ReferenceGenome | None:
    return ReferenceGenome.query.filter_by(is_default=True).first()


def reference_contigs(ref: ReferenceGenome) -> dict[str, int]:
    return read_fai(Path(ref.fasta_path + ".fai"))
