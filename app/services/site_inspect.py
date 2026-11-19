"""Look at the reads behind one variant: why did (or didn't) a caller call it?

Used by the variant-by-variant viewer ("Check reads"). Reads the BAM at one position
(read-only, pysam) and reports, per read: mapping quality (MAPQ), base, base quality,
duplicate flag and where else the aligner says the read fits (BWA's XA tag). From that
it explains the usual reasons a variant visible in another pipeline is missing here:

* MAPQ 0 / low MAPQ - the read fits equally well somewhere else (paralogue, unplaced
  contig, segmental duplication); callers ignore such reads (GATK/bcftools/FreeBayes
  default MAPQ >= 20 here, DeepVariant >= 5).
* low base quality, duplicates, too few ALT reads, or simply no reads.

Nothing is written; no read is changed or discarded.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

MAX_READS = 5000


def bam_for_vcf(sample, vcf_path: str):
    """The BAM that produced a pipeline VCF: a BAM registered inside the same run folder, or the
    BAM that run used as input. None for external VCFs (the user picks a BAM)."""
    from ..models import SampleFile, WorkflowRun
    vcf = str(Path(vcf_path).resolve())
    for run in WorkflowRun.query.filter_by(sample_pk=sample.id).all():
        out = str(Path(run.output_dir).resolve())
        if not vcf.startswith(out + "/"):
            continue
        for f in sample.files:
            if f.role == "bam" and str(Path(f.path).resolve()).startswith(out + "/"):
                return f
        for fid in run.input_file_ids or []:
            f = SampleFile.query.get(fid)
            if f is not None and f.role == "bam":
                return f
    return None


def _has_index(bam: Path) -> bool:
    return any(Path(str(bam) + ext).exists() for ext in (".bai", ".csi")) or bam.with_suffix(".bai").exists()


def inspect_site(bam: Path, chrom: str, pos: int, ref: str, alt: str, min_mq: int = 20, min_bq: int = 10,
                 index: Path | None = None) -> dict:
    """pos is 1-based (VCF). For SNVs the base at pos is compared with REF/ALT; for indels the
    reads with an insertion/deletion starting at pos are counted as ALT."""
    import pysam
    bam = Path(bam)
    if not bam.is_file():
        raise ValueError(f"BAM not found: {bam}")
    if index is None and not _has_index(bam):
        raise ValueError(f"{bam.name} has no index (.bai/.csi); random access is not possible. "
                         "Index it with 'samtools index' (writes a new .bai next to it) to inspect it here.")
    ref, alt = ref.upper(), alt.upper()
    snv = len(ref) == 1 and len(alt) == 1
    ins = len(alt) > len(ref)
    indel_len = abs(len(alt) - len(ref))
    reads, truncated = [], False
    with pysam.AlignmentFile(str(bam), index_filename=str(index) if index else None) as fh:
        if chrom not in fh.references:
            raise ValueError(f"{chrom} is not a contig of {bam.name}")
        for col in fh.pileup(chrom, pos - 1, pos, truncate=True, stepper="nofilter", ignore_overlaps=False,
                             min_base_quality=0, max_depth=MAX_READS + 1, ignore_orphans=False):
            for pr in col.pileups:
                if len(reads) >= MAX_READS:
                    truncated = True
                    break
                a = pr.alignment
                if a.is_secondary or a.is_supplementary:
                    continue
                if pr.is_refskip:
                    continue
                if snv:
                    if pr.is_del:
                        allele, bq = "del", None
                    else:
                        b = a.query_sequence[pr.query_position]
                        bq = a.query_qualities[pr.query_position] if a.query_qualities is not None else None
                        allele = "ALT" if b == alt else ("REF" if b == ref else b)
                else:
                    bq = (a.query_qualities[pr.query_position]
                          if pr.query_position is not None and a.query_qualities is not None else None)
                    if ins:
                        allele = "ALT" if pr.indel == indel_len else ("REF" if pr.indel == 0 else "other")
                    else:
                        allele = "ALT" if pr.indel == -indel_len else ("REF" if pr.indel == 0 else "other")
                xa = a.get_tag("XA") if a.has_tag("XA") else ""
                reads.append({"mapq": a.mapping_quality, "allele": allele, "bq": bq, "dup": a.is_duplicate,
                              "qcfail": a.is_qcfail, "reverse": a.is_reverse,
                              "xa": [x.split(",")[0] for x in xa.split(";") if x],
                              "nm": a.get_tag("NM") if a.has_tag("NM") else None,
                              "as": a.get_tag("AS") if a.has_tag("AS") else None,
                              "xs": a.get_tag("XS") if a.has_tag("XS") else None})
    return summarise_site(reads, min_mq, min_bq, truncated, bam, chrom, pos, ref, alt)


def _mq_bin(q):
    return "0" if q == 0 else ("1-19" if q < 20 else ("20-39" if q < 40 else "≥40"))


def summarise_site(reads, min_mq, min_bq, truncated, bam, chrom, pos, ref, alt) -> dict:
    n = len(reads)
    usable = [r for r in reads if not r["dup"] and not r["qcfail"]]
    good = [r for r in usable if r["mapq"] >= min_mq and (r["bq"] is None or r["bq"] >= min_bq)]
    alleles_all = Counter(r["allele"] for r in usable)
    alleles_good = Counter(r["allele"] for r in good)
    mq_bins = Counter(_mq_bin(r["mapq"]) for r in usable)
    alt_mq = Counter(_mq_bin(r["mapq"]) for r in usable if r["allele"] == "ALT")
    xa = Counter(c for r in usable for c in r["xa"])
    alt_strand = Counter("reverse" if r["reverse"] else "forward" for r in good if r["allele"] == "ALT")
    reasons = []
    if n == 0:
        reasons.append("No reads cover this position in our BAM (outside the capture, or reads were not "
                       "aligned here).")
    else:
        mq0 = sum(1 for r in usable if r["mapq"] == 0)
        if usable and mq0 / len(usable) >= 0.5:
            where = ", ".join(f"{c} ({k})" for c, k in xa.most_common(3)) or "several places (no XA tag)"
            reasons.append(f"{mq0} of {len(usable)} reads have MAPQ 0: the aligner found an equally good "
                           f"alignment elsewhere — {where}. Variant callers ignore reads below MAPQ {min_mq} "
                           "(DeepVariant: 5), so the site has (almost) no usable reads.")
        elif usable and sum(1 for r in usable if r["mapq"] < min_mq) / len(usable) >= 0.5:
            reasons.append(f"Most reads have MAPQ below {min_mq} (uncertain placement); callers ignore them.")
        if alleles_all.get("ALT", 0) and alleles_good.get("ALT", 0) < 2:
            reasons.append(f"Only {alleles_good.get('ALT', 0)} ALT read(s) pass MAPQ ≥ {min_mq} and base "
                           f"quality ≥ {min_bq} (of {alleles_all['ALT']} ALT reads in total).")
        if not alleles_all.get("ALT"):
            reasons.append("No read in our BAM carries the ALT allele here.")
        dups = n - len(usable)
        if n and dups / n >= 0.5:
            reasons.append(f"{dups} of {n} reads are marked duplicate (or QC-fail) and are ignored by callers.")
        if len(alt_strand) == 1 and sum(alt_strand.values()) >= 4:
            reasons.append(f"All usable ALT reads are on the {next(iter(alt_strand))} strand (strand bias).")
    if not reasons and alleles_good.get("ALT", 0) >= 2:
        reasons.append("Our reads support the ALT allele with good mapping and base quality — the caller should "
                       "be able to see it; check the caller's raw VCF and its own filters/thresholds.")
    return {"bam": str(bam), "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "min_mq": min_mq,
            "min_bq": min_bq, "reads": n, "truncated": truncated, "duplicates_or_qcfail": n - len(usable),
            "alleles_all": dict(alleles_all), "alleles_good": dict(alleles_good),
            "mapq_bins": [(b, mq_bins.get(b, 0)) for b in ("0", "1-19", "20-39", "≥40")],
            "alt_mapq_bins": [(b, alt_mq.get(b, 0)) for b in ("0", "1-19", "20-39", "≥40")],
            "alt_strand": dict(alt_strand), "xa_contigs": xa.most_common(8),
            "examples": reads[:30], "reasons": reasons}
