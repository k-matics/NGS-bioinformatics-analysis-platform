"""Read BAM/CRAM/VCF headers (no reads are loaded)."""
from __future__ import annotations

from pathlib import Path

import pysam


def bam_header_summary(path: Path) -> dict:
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as fh:
        header = fh.header.to_dict()
    contigs = {sq["SN"]: int(sq["LN"]) for sq in header.get("SQ", [])}
    programs = [
        {"id": pg.get("ID", "").strip(), "name": pg.get("PN", ""), "version": str(pg.get("VN", "")),
         "command": (pg.get("CL", "") or "")[:500]}
        for pg in header.get("PG", [])
    ]
    read_groups = [
        {k: rg.get(k) for k in ("ID", "SM", "LB", "PL", "PU", "CN")} for rg in header.get("RG", [])
    ]
    index_candidates = [Path(str(path) + ".bai"), path.with_suffix(".bai"), Path(str(path) + ".csi")]
    return {
        "sort_order": header.get("HD", {}).get("SO", "unknown"),
        "format_version": header.get("HD", {}).get("VN"),
        "contigs": contigs,
        "read_groups": read_groups,
        "samples": sorted({rg.get("SM") for rg in read_groups if rg.get("SM")}),
        "programs": programs,
        "index": next((str(p) for p in index_candidates if p.exists()), None),
    }


def vcf_header_summary(path: Path) -> dict:
    with pysam.VariantFile(str(path)) as fh:
        contigs = {name: rec.length for name, rec in fh.header.contigs.items() if rec.length}
        samples = list(fh.header.samples)
        sources = [r.value for r in fh.header.records if r.key in ("source", "reference")]
    index_candidates = [Path(str(path) + ".tbi"), Path(str(path) + ".csi")]
    return {
        "contigs": contigs,
        "samples": samples,
        "source": sources,
        "index": next((str(p) for p in index_candidates if p.exists()), None),
    }
