"""Variant-set comparison engine (used by the in-run caller comparison and by VCF-vs-VCF).

Two complementary methods are used, because they answer different questions:

1. **Normalised allele matching** (this module, Python/pysam). Every VCF is first
   made comparable: multi-allelic records split, complex/MNP records decomposed
   into their SNVs/indels ("atomised"), indels left-aligned against OUR reference.
   A variant is then the key (CHROM, POS, REF, ALT). This gives the concordance
   matrix (which variant was found by which caller), intersections, and genotype
   concordance for shared variants.

2. **Haplotype-aware matching** with RTG ``vcfeval``. Some differences are only
   *representation* differences (e.g. one caller writes an indel + nearby SNV as one
   complex allele): vcfeval replays the variants onto the reference sequence and
   compares the resulting haplotypes, so equivalent representations match. The
   difference between (1) and (2) is reported as "representation differences".

Genotype comparison: GT alleles are sorted (0|1 == 1|0 == 0/1). Records whose
sample genotype has no ALT allele (0/0, ./.) are not variants for that caller.
"""
from __future__ import annotations

import csv
import gzip
import itertools
import os
from collections import Counter, defaultdict
from pathlib import Path

import pysam

from ...services.settings_service import get_setting

VARIANT_TYPES = ("SNV", "indel", "MNP/complex")


# --------------------------------------------------------------------------- preparation
def normalise_for_comparison(runner, vcf: Path, ref: Path, out: Path, *, regions: Path | None = None,
                             sample: str | None = None, threads: int = 2) -> Path:
    """Make a VCF comparable: (subset sample) -> restrict to regions -> atomise -> split ->
    left-align -> sort -> bgzip + index. Streams (no index needed on the input)."""
    view = ["bcftools", "view", "-Ou"]
    if sample:
        view += ["-s", sample, "--trim-alt-alleles"]
    if regions:
        view += ["-T", regions]          # -T streams (works on unindexed external files)
    view.append(vcf)
    atom = ["bcftools", "norm", "-a", "--atom-overlaps", ".", "-Ou"]
    split = ["bcftools", "norm", "-f", ref, "-m", "-any", "--check-ref", "w", "-Ou"]
    tmp = out.parent / (out.name + ".sort_tmp")
    sort = ["bcftools", "sort", "-T", tmp, "-Oz", "-o", out]
    runner.run([view, atom, split, sort])
    runner.run(["bcftools", "index", "-f", "-t", out])
    return out


def ref_mismatch_rate(vcf: Path, fasta: Path, n: int = 5000) -> tuple[int, int]:
    """(checked, mismatched) REF alleles vs the reference for the first n records.
    A high mismatch rate means the VCF is from another genome build."""
    checked = bad = 0
    with pysam.VariantFile(str(vcf)) as fh, pysam.FastaFile(str(fasta)) as fa:
        refs = set(fa.references)
        for rec in fh:
            if rec.chrom not in refs or not rec.ref or len(rec.ref) > 50 or "N" in rec.ref.upper():
                continue
            seq = fa.fetch(rec.chrom, rec.pos - 1, rec.pos - 1 + len(rec.ref)).upper()
            checked += 1
            bad += seq != rec.ref.upper()
            if checked >= n:
                break
    return checked, bad


# --------------------------------------------------------------------------- loading
def variant_type(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) != len(alt):
        return "indel"
    return "MNP/complex"


def load_variants(vcf: Path, pass_only: bool = True) -> tuple[dict, Counter]:
    """{(chrom,pos,ref,alt): (gt, zygosity, type)} plus counts of skipped records by reason."""
    out: dict = {}
    skipped: Counter = Counter()
    with pysam.VariantFile(str(vcf)) as fh:
        if not fh.header.samples:
            raise ValueError(f"{vcf.name} has no sample column (sites-only VCF).")
        for rec in fh:
            alt = rec.alts[0] if rec.alts else None
            if alt is None or alt in ("*", "<*>", "<NON_REF>") or alt.startswith("<"):
                skipped["no/symbolic ALT"] += 1
                continue
            filters = [f for f in rec.filter.keys() if f != "PASS"]
            if pass_only and filters:
                skipped["not PASS"] += 1
                continue
            gt = rec.samples[0].get("GT") or ()
            alleles = sorted(a for a in gt if a is not None)
            if not alleles or max(alleles) == 0:
                skipped["no ALT in genotype"] += 1
                continue
            gts = "/".join(str(a) for a in alleles)
            zyg = "hom" if len(set(alleles)) == 1 and len(alleles) > 1 else (
                "hemi" if len(alleles) == 1 else "het")
            key = (rec.chrom, rec.pos, rec.ref.upper(), alt.upper())
            if key in out:
                skipped["duplicate record"] += 1
                continue
            out[key] = (gts, zyg, variant_type(rec.ref, alt))
    return out, skipped


def count_multiallelic(vcf: Path) -> int:
    n = 0
    with pysam.VariantFile(str(vcf)) as fh:
        for rec in fh:
            n += bool(rec.alts) and len(rec.alts) > 1
    return n


# --------------------------------------------------------------------------- comparison
def compare_sets(sets: dict[str, dict]) -> dict:
    names = list(sets)
    union = set().union(*[set(v) for v in sets.values()]) if sets else set()
    per = {}
    for n in names:
        v = sets[n]
        c = Counter(t for _, _, t in v.values())
        z = Counter(zg for _, zg, _ in v.values())
        others = set().union(*[set(sets[m]) for m in names if m != n]) if len(names) > 1 else set()
        per[n] = {"total": len(v), "snv": c["SNV"], "indel": c["indel"], "mnp_complex": c["MNP/complex"],
                  "het": z["het"], "hom": z["hom"], "hemi": z["hemi"],
                  "unique": len(set(v) - others)}
    # intersections (UpSet-style): which exact combination of sets contains each variant
    inter: dict = defaultdict(Counter)
    for key in union:
        pattern = tuple(n for n in names if key in sets[n])
        t = next(sets[n][key][2] for n in pattern)
        inter[pattern]["total"] += 1
        inter[pattern][t] += 1
    intersections = sorted(
        ({"sets": list(p), "n_sets": len(p), "total": c["total"], "snv": c["SNV"], "indel": c["indel"],
          "mnp_complex": c["MNP/complex"]} for p, c in inter.items()),
        key=lambda r: -r["total"])
    pairwise = []
    for a, b in itertools.combinations(names, 2):
        sa, sb = set(sets[a]), set(sets[b])
        shared = sa & sb
        same_gt = sum(sets[a][k][0] == sets[b][k][0] for k in shared)
        pairwise.append({
            "a": a, "b": b, "shared": len(shared), "a_only": len(sa - sb), "b_only": len(sb - sa),
            "jaccard": round(len(shared) / len(sa | sb), 4) if sa | sb else None,
            "gt_concordant": same_gt, "gt_discordant": len(shared) - same_gt,
            "gt_concordance_pct": round(100 * same_gt / len(shared), 2) if shared else None,
        })
    in_all = sum(1 for k in union if all(k in sets[n] for n in names))
    return {"sets": names, "union": len(union), "in_all": in_all, "per_set": per,
            "intersections": intersections, "pairwise": pairwise}


def write_tables(sets: dict[str, dict], result: dict, out_dir: Path) -> dict:
    """matrix.tsv.gz (every variant x every set), discordant genotypes, CSV summaries."""
    names = result["sets"]
    union = sorted(set().union(*[set(v) for v in sets.values()]),
                   key=lambda k: (k[0], k[1], k[2], k[3]))
    files = {}
    matrix = out_dir / "variant_matrix.tsv.gz"
    disc = out_dir / "discordant_genotypes.tsv"
    with gzip.open(matrix, "wt") as fm, open(disc, "w") as fd:
        hdr = ["CHROM", "POS", "REF", "ALT", "TYPE", "N_SETS", *[f"GT_{n}" for n in names]]
        fm.write("\t".join(hdr) + "\n")
        fd.write("\t".join(hdr) + "\n")
        for k in union:
            gts = [sets[n][k][0] if k in sets[n] else "." for n in names]
            t = next(sets[n][k][2] for n in names if k in sets[n])
            row = [k[0], str(k[1]), k[2], k[3], t, str(sum(g != "." for g in gts)), *gts]
            fm.write("\t".join(row) + "\n")
            present = [g for g in gts if g != "."]
            if len(set(present)) > 1:
                fd.write("\t".join(row) + "\n")
    files["matrix"] = str(matrix)
    files["discordant_genotypes"] = str(disc)
    files["per_set_csv"] = _csv(out_dir / "per_set_counts.csv",
                                [{"set": n, **result["per_set"][n]} for n in names])
    files["intersections_csv"] = _csv(out_dir / "intersections.csv",
                                      [{**r, "sets": "+".join(r["sets"])} for r in result["intersections"]])
    files["pairwise_csv"] = _csv(out_dir / "pairwise.csv", result["pairwise"])
    return files


def _csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> str:
    """Write rows as CSV. An empty table still gets its header (an empty result is a result)."""
    fields = list(rows[0]) if rows else (fields or ["no_rows"])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    return str(path)


# --------------------------------------------------------------------------- RTG vcfeval
def sdf_dir_for(reference) -> Path:
    base = get_setting("rtg_sdf_dir") or ""
    if not base:
        from flask import current_app
        base = Path(current_app.instance_path) / "rtg_sdf"
    return Path(base) / f"ref{reference.id}_{Path(reference.fasta_path).stem}"


def rtg_env() -> dict:
    env = dict(os.environ)
    env["RTG_MEM"] = f"{int(get_setting('rtg_mem_gb', 8) or 8)}G"
    return env


def ensure_sdf(runner, reference, log) -> Path:
    """RTG works on its own binary copy of the reference (SDF). It is built once per reference,
    outside the reference directory (the reference itself is never modified)."""
    sdf = sdf_dir_for(reference)
    if (sdf / "done").exists():
        return sdf
    sdf.parent.mkdir(parents=True, exist_ok=True)
    tmp = sdf.parent / (sdf.name + f".building.{os.getpid()}")
    log(f"Building RTG SDF of the reference (one-off) in {sdf} ...")
    runner.run(["rtg", "format", "-o", tmp, reference.fasta_path], env=rtg_env())
    (tmp / "done").write_text(reference.fasta_path + "\n")
    tmp.rename(sdf)
    return sdf


def _roc_last(path: Path) -> dict:
    """Last data row of a vcfeval ROC file = all calls (no score threshold)."""
    last = None
    cols = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") and "true_positives_baseline" in line:
                cols = line.lstrip("#").rstrip("\n").split("\t")
            elif not line.startswith("#") and line.strip():
                last = line.rstrip("\n").split("\t")
    if not last or not cols:
        return {}
    d = dict(zip(cols, last))
    return {"tp_baseline": int(float(d.get("true_positives_baseline", 0))),
            "tp_call": int(float(d.get("true_positives_call", 0))),
            "fp": int(float(d.get("false_positives", 0))),
            "fn": int(float(d.get("false_negatives", 0))),
            "precision": _f(d.get("precision")), "sensitivity": _f(d.get("sensitivity")),
            "f_measure": _f(d.get("f_measure"))}


def _f(x):
    try:
        return round(float(x), 4)
    except (TypeError, ValueError):
        return None


def run_vcfeval(runner, sdf: Path, baseline: Path, calls: Path, out_dir: Path, *, regions: Path | None,
                eval_regions: Path | None = None, threads: int = 4, squash_ploidy: bool = False,
                all_records: bool = False) -> dict:
    """baseline vs calls. Returns overall + SNP / non-SNP metrics (all calls, no score threshold)."""
    argv = ["rtg", "vcfeval", "-b", baseline, "-c", calls, "-t", sdf, "-o", out_dir,
            "--threads", str(threads), "--ref-overlap", "--vcf-score-field", "QUAL"]
    if regions:
        argv += ["--bed-regions", regions]
    if eval_regions:
        argv += ["--evaluation-regions", eval_regions]
    if squash_ploidy:
        argv.append("--squash-ploidy")
    if all_records:
        argv.append("--all-records")
    runner.run(argv, env=rtg_env())
    res = {"overall": _roc_last(out_dir / "weighted_roc.tsv.gz"),
           "snp": _roc_last(out_dir / "snp_roc.tsv.gz"),
           "non_snp": _roc_last(out_dir / "non_snp_roc.tsv.gz"),
           "dir": str(out_dir)}
    return res
