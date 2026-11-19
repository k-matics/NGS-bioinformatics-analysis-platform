"""BAM vs VCF: read-level evidence for every variant of a VCF.

For each (normalised, split) variant the reads of the BAM covering the position
are inspected with pysam:

  ref / alt / other read counts, depth, observed VAF (alt / (ref+alt)),
  forward / reverse strand counts and a Fisher strand-bias p-value,
  mean MAPQ of ref and alt reads and the fraction of MAPQ-0 reads,
  mean base quality of alt bases (SNVs), fraction of alt reads that are soft-clipped,
  median distance of the variant from the nearest read end in alt reads.

Duplicates, QC-fail, secondary and supplementary reads are excluded (as callers
do). Evidence *notes* are attached (e.g. "VAF inconsistent with GT"), but a
note is NOT a verdict: interpretation needs the context (region, IGV, other
callers). Nothing is classified as an error automatically.
"""
from __future__ import annotations

import gzip
import math
from pathlib import Path
from statistics import median

import pysam

from ...engine.step import Step, StepContext, StepFailure, StepResult
from ..common import prior_output, read_json, write_json
from .engine import _csv

EXCLUDE = 0x4 | 0x100 | 0x200 | 0x400 | 0x800   # unmapped, secondary, qc-fail, duplicate, supplementary


def fisher_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact test for [[a, b], [c, d]] (ref fwd/rev, alt fwd/rev)."""
    n = a + b + c + d
    if n == 0:
        return 1.0
    r1, c1 = a + b, a + c

    def lp(x):
        return (math.lgamma(r1 + 1) + math.lgamma(n - r1 + 1) + math.lgamma(c1 + 1) + math.lgamma(n - c1 + 1)
                - math.lgamma(n + 1) - math.lgamma(x + 1) - math.lgamma(r1 - x + 1) - math.lgamma(c1 - x + 1)
                - math.lgamma(n - r1 - c1 + x + 1))
    lo, hi = max(0, r1 + c1 - n), min(r1, c1)
    p0 = lp(a)
    p = sum(math.exp(lp(x)) for x in range(lo, hi + 1) if lp(x) <= p0 + 1e-7)
    return min(1.0, p)


def _allele(pr, ref: str, alt: str):
    """Classify one pileup read at the VCF anchor base: 'ref' | 'alt' | 'other' | None (skip)."""
    if pr.is_refskip:
        return None
    if pr.is_del:
        return "other"                      # read has a deletion spanning the anchor base
    seq = pr.alignment.query_sequence or ""
    qp = pr.query_position
    if len(ref) == 1 and len(alt) == 1:
        b = seq[qp].upper()
        return "alt" if b == alt else "ref" if b == ref else "other"
    if len(alt) > len(ref):                 # insertion after the anchor base
        if pr.indel > 0:
            return "alt" if seq[qp + 1:qp + 1 + pr.indel].upper() == alt[1:] else "other"
        return "other" if pr.indel < 0 else "ref"
    if len(ref) > len(alt):                 # deletion after the anchor base
        if pr.indel < 0:
            return "alt" if -pr.indel == len(ref) - len(alt) else "other"
        return "other" if pr.indel > 0 else "ref"
    return None                             # MNP/complex: not evaluated


def site_evidence(bam: pysam.AlignmentFile, chrom: str, pos: int, ref: str, alt: str, min_bq: int) -> dict:
    pos0 = pos - 1
    snv = len(ref) == 1 and len(alt) == 1
    c = {"ref": 0, "alt": 0, "other": 0}
    strand = {"ref": [0, 0], "alt": [0, 0]}
    mq = {"ref": [], "alt": []}
    mq0 = 0
    alt_bq, alt_clip, alt_end = [], 0, []
    for col in bam.pileup(chrom, pos0, pos0 + 1, truncate=True, stepper="nofilter", ignore_overlaps=False,
                          min_base_quality=0, max_depth=1_000_000):
        for pr in col.pileups:
            read = pr.alignment
            if read.flag & EXCLUDE:
                continue
            allele = _allele(pr, ref, alt)
            if allele is None:
                continue
            q = pr.query_position
            if snv and q is not None and read.query_qualities is not None and read.query_qualities[q] < min_bq:
                continue
            c[allele] += 1
            mq0 += read.mapping_quality == 0
            if allele in strand:
                strand[allele][1 if read.is_reverse else 0] += 1
                mq[allele].append(read.mapping_quality)
            if allele == "alt":
                if snv and q is not None and read.query_qualities is not None:
                    alt_bq.append(read.query_qualities[q])
                if any(op == 4 for op, _ in (read.cigartuples or [])):
                    alt_clip += 1
                if q is not None:
                    alt_end.append(min(q, read.query_length - 1 - q))
    depth = c["ref"] + c["alt"] + c["other"]
    ra = c["ref"] + c["alt"]
    p_sb = fisher_two_sided(strand["ref"][0], strand["ref"][1], strand["alt"][0], strand["alt"][1])
    return {
        "depth": depth, "ref_reads": c["ref"], "alt_reads": c["alt"], "other_reads": c["other"],
        "vaf": round(c["alt"] / ra, 3) if ra else None,
        "alt_fwd": strand["alt"][0], "alt_rev": strand["alt"][1],
        "ref_fwd": strand["ref"][0], "ref_rev": strand["ref"][1],
        "strand_bias_p": round(p_sb, 6),
        "mean_mq_ref": round(sum(mq["ref"]) / len(mq["ref"]), 1) if mq["ref"] else None,
        "mean_mq_alt": round(sum(mq["alt"]) / len(mq["alt"]), 1) if mq["alt"] else None,
        "mq0_frac": round(mq0 / depth, 3) if depth else None,
        "mean_bq_alt": round(sum(alt_bq) / len(alt_bq), 1) if alt_bq else None,
        "alt_softclip_frac": round(alt_clip / c["alt"], 3) if c["alt"] else None,
        "alt_median_dist_to_end": median(alt_end) if alt_end else None,
    }


NOTE_LABELS = {
    "no_reads": "no usable reads at the position", "no_alt": "no ALT reads in BAM",
    "low_alt": "low ALT support (<3 reads)", "vaf_het": "VAF unusual for a heterozygous call",
    "vaf_hom": "VAF low for a homozygous-ALT call", "strand": "ALT reads on one strand only",
    "low_mq": "low mapping quality of ALT reads (<30)", "mq0": ">20% MAPQ-0 reads (repetitive/ALT region)",
    "softclip": "most ALT reads soft-clipped", "read_end": "ALT mostly near read ends (<10 bp)",
    "third_allele": ">20% reads with a third allele / other indel",
}


def notes_for(ev: dict, gt_zyg: str | None) -> list[str]:
    """Evidence note codes (see NOTE_LABELS). Notes are observations, not verdicts."""
    if ev["depth"] == 0:
        return ["no_reads"]
    n = []
    if ev["alt_reads"] == 0:
        n.append("no_alt")
    elif ev["alt_reads"] < 3:
        n.append("low_alt")
    v = ev["vaf"]
    if v is not None and ev["alt_reads"]:
        if gt_zyg == "het" and not 0.2 <= v <= 0.8:
            n.append("vaf_het")
        if gt_zyg == "hom" and v < 0.8:
            n.append("vaf_hom")
    if ev["alt_reads"] >= 4 and (ev["alt_fwd"] == 0 or ev["alt_rev"] == 0) and ev["strand_bias_p"] < 0.01:
        n.append("strand")
    if ev["mean_mq_alt"] is not None and ev["mean_mq_alt"] < 30:
        n.append("low_mq")
    if ev["mq0_frac"] is not None and ev["mq0_frac"] > 0.2:
        n.append("mq0")
    if ev["alt_softclip_frac"] is not None and ev["alt_reads"] >= 3 and ev["alt_softclip_frac"] > 0.5:
        n.append("softclip")
    if ev["alt_median_dist_to_end"] is not None and ev["alt_reads"] >= 3 and ev["alt_median_dist_to_end"] < 10:
        n.append("read_end")
    if ev["other_reads"] and ev["other_reads"] / ev["depth"] > 0.2:
        n.append("third_allele")
    return n


class ReadEvidenceStep(Step):
    key = "read_evidence"
    name = "Read-level evidence (BAM vs VCF)"
    description = ("For every variant: REF/ALT read counts, VAF, strand balance, mapping and base quality, "
                   "soft-clipping and read-end position from the BAM, plus evidence notes and IGV files.")
    tools = ("pysam",)

    def execute(self, ctx: StepContext, runner) -> StepResult:
        inputs = read_json(Path(prior_output(ctx, "vcf_prepare", "inputs_json")))
        if len(inputs) != 1:
            raise StepFailure("Select exactly one VCF for BAM-vs-VCF evidence.")
        inp = inputs[0]
        bam_path = prior_output(ctx, "bam_prepare", "bam")
        pass_only = bool(ctx.params.get("evidence_pass_only", False))
        max_sites = int(ctx.params.get("evidence_max_sites", 200_000))
        min_bq = int(ctx.params.get("evidence_min_bq", 13))
        d = ctx.step_dir("comparison/evidence")
        rows, note_counts = [], {}
        n_seen = 0
        with pysam.VariantFile(inp["vcf"]) as vf, pysam.AlignmentFile(bam_path) as bam:
            bam_contigs = set(bam.references)
            for rec in vf:
                if not rec.alts or rec.alts[0].startswith("<") or rec.alts[0] == "*":
                    continue
                filt = [f for f in rec.filter.keys() if f != "PASS"]
                if pass_only and filt:
                    continue
                n_seen += 1
                if n_seen > max_sites:
                    raise StepFailure(f"The VCF has more than {max_sites} variants to evaluate. Restrict it with a "
                                      "BED, use PASS-only, or raise 'Max variants' (nothing was sampled silently).")
                if rec.chrom not in bam_contigs:
                    continue
                gt = rec.samples[0].get("GT") or ()
                al = sorted(a for a in gt if a is not None)
                zyg = None if not al or max(al) == 0 else ("hom" if len(set(al)) == 1 else "het")
                vcf_ad = rec.samples[0].get("AD") if "AD" in rec.format else None
                ev = site_evidence(bam, rec.chrom, rec.pos, rec.ref.upper(), rec.alts[0].upper(), min_bq)
                notes = notes_for(ev, zyg)
                for x in notes:
                    note_counts[x] = note_counts.get(x, 0) + 1
                rows.append({"chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref, "alt": rec.alts[0],
                             "filter": ";".join(filt) or "PASS", "gt": "/".join(map(str, al)) or ".",
                             "vcf_qual": rec.qual,
                             "vcf_ad": ",".join(str(x) for x in vcf_ad) if vcf_ad and None not in vcf_ad else "",
                             **ev, "notes": ";".join(notes)})
                if n_seen % 2000 == 0:
                    ctx.check_cancel()
                    ctx.log(f"... {n_seen} variants evaluated")
        tsv = d / "read_evidence.tsv.gz"
        cols = list(rows[0]) if rows else ["chrom", "pos", "ref", "alt", "notes"]
        with gzip.open(tsv, "wt") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")
        flagged = [r for r in rows if r["notes"]]
        _csv(d / "flagged_variants.csv", flagged[:20000], cols)
        bed = d / "igv_flagged_sites.bed"
        with open(bed, "w") as fh:
            fh.write(f'track name="flagged_{ctx.sample.sample_id}" description="Variants with evidence notes '
                     f'({len(flagged)})"\n')
            for r in flagged:
                fh.write(f"{r['chrom']}\t{max(0, r['pos'] - 1)}\t{r['pos'] - 1 + len(r['ref'])}\t"
                         f"{r['ref']}>{r['alt']}|{r['notes'][:80]}\n")
        batch = d / "igv_batch.txt"
        with open(batch, "w") as fh:
            fh.write(f"new\ngenome {ctx.reference.fasta_path}\nload {bam_path}\nload {inp['source']}\n"
                     f"snapshotDirectory {d / 'igv_snapshots'}\n")
            for r in flagged[:500]:
                fh.write(f"goto {r['chrom']}:{max(1, r['pos'] - 50)}-{r['pos'] + 50}\nsnapshot\n")
        vafs = [r["vaf"] for r in rows if r["vaf"] is not None and r["alt_reads"]]
        summary = {"vcf": inp["source"], "bam": bam_path, "evaluated": len(rows), "flagged": len(flagged),
                   "pass_only": pass_only, "note_counts": dict(sorted(note_counts.items(), key=lambda x: -x[1])),
                   "note_labels": NOTE_LABELS,
                   "no_alt_reads": sum(1 for r in rows if r["alt_reads"] == 0),
                   "median_vaf_het": _median([r["vaf"] for r in rows if r["gt"] == "0/1" and r["vaf"] is not None]),
                   "median_vaf_hom": _median([r["vaf"] for r in rows if r["gt"] == "1/1" and r["vaf"] is not None]),
                   "vaf_hist": _hist(vafs)}
        write_json(d / "read_evidence_summary.json", summary)
        from flask import render_template
        from ...models.base import utcnow
        from ...services.limitations import run_limitations
        html = d / "read_evidence_report.html"
        html.write_text(render_template("reports/evidence_report.html", s=summary, top=flagged[:300],
                                        run_id=ctx.run.run_id, sample=ctx.sample.sample_id,
                                        reference=ctx.reference.name, generated_at=utcnow().isoformat() + "Z",
                                        limitations=run_limitations(ctx.reference, ctx.params)))
        return StepResult(outputs={"evidence_tsv": str(tsv), "flagged_csv": str(d / "flagged_variants.csv"),
                                   "igv_bed": str(bed), "igv_batch": str(batch), "report_html": str(html)},
                          summary=summary)


def _median(x):
    return round(median(x), 3) if x else None


def _hist(vafs, bins=20):
    h = [0] * bins
    for v in vafs:
        h[min(bins - 1, int(v * bins))] += 1
    return h
