"""Variant-by-variant check: is every variant of a *baseline* VCF (e.g. the provider's) in our VCF?

This is a new, separate workflow. The existing VCF-vs-VCF comparison is unchanged; this one
answers a different question, for one baseline variant at a time:

* **found_pass**        – our VCF has the same allele and it PASSes our filters
* **found_filtered**    – our VCF has it, but one of our filters labelled it (which one is listed)
* **found_ref_gt**      – our caller looked at the site but genotyped it as reference / no-call
                          (e.g. DeepVariant RefCall)
* **found_other_repr**  – not the same CHROM/POS/REF/ALT, but RTG vcfeval shows the same haplotype
                          (the same variant written differently)
* **nearby**            – none of the above, but our VCF has another variant within N bp
                          (a different allele, or a shifted indel)
* **not_called**        – nothing in our VCF

Our VCF is read with **all** records (filtered ones too) by default, so a permissive VCF can be
checked and then filtered interactively in the viewer ("what if our min depth were 5?") without
re-running anything. Both VCFs have been normalised identically by the VcfInputsStep before this
step (split, atomised, left-aligned against our reference, restricted to the same regions).
"""
from __future__ import annotations

import bisect
import gzip
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import pysam

from ...engine.step import Step, StepContext, StepFailure, StepResult
from ..common import prior_output, read_json, vcf_stem, write_json
from . import engine

STATUS_LABELS = {
    "found_pass": "Found in our VCF, PASS",
    "found_filtered": "Found, but filtered by our filters",
    "found_ref_gt": "Site in our VCF, genotyped as reference / no-call",
    "found_other_repr": "Found, written differently (haplotype match, RTG vcfeval)",
    "nearby": "Not found; our VCF has a different variant nearby",
    "not_called": "Not in our VCF",
}
FOUND = ("found_pass", "found_filtered", "found_ref_gt", "found_other_repr")
COLUMNS = ["CHROM", "POS", "REF", "ALT", "TYPE", "STATUS", "GT_MATCH",
           "B_GT", "B_QUAL", "B_FILTER", "B_DP", "B_AD_REF", "B_AD_ALT", "B_VAF", "B_GQ", "B_MQ",
           "Q_GT", "Q_QUAL", "Q_FILTER", "Q_DP", "Q_AD_REF", "Q_AD_ALT", "Q_VAF", "Q_GQ", "Q_MQ", "NEARBY"]
Q_COLUMNS = ["CHROM", "POS", "REF", "ALT", "TYPE", "Q_GT", "Q_QUAL", "Q_FILTER", "Q_DP", "Q_AD_REF",
             "Q_AD_ALT", "Q_VAF", "Q_GQ", "Q_MQ"]


# ------------------------------------------------------------------------------ record reading
def _num(x, nd=2):
    if x is None:
        return None
    try:
        return round(float(x), nd)
    except (TypeError, ValueError):
        return None


def record_info(rec) -> dict:
    """GT / QUAL / FILTER / DP / AD / VAF / GQ of the first sample of a (split) record, and the
    site mapping quality (INFO/MQ: GATK, bcftools, DRAGEN; INFO/MQM: FreeBayes; DeepVariant has none).
    MQ scales differ: DRAGEN's VCF MQ goes above 60 (e.g. 250), BWA/GATK's stops at 60."""
    s = rec.samples[0] if rec.samples else None
    gt = (s.get("GT") if s is not None else None) or ()
    alleles = sorted(a for a in gt if a is not None)
    fl = list(rec.filter.keys())
    dp = ad_r = ad_a = vaf = gq = None
    if s is not None:
        fmt = rec.format
        if "DP" in fmt:
            dp = s.get("DP")
        if "AD" in fmt:
            ad = s.get("AD")
            if ad and len(ad) >= 2 and ad[0] is not None and ad[1] is not None:
                ad_r, ad_a = ad[0], ad[1]
                if ad_r + ad_a > 0:
                    vaf = ad_a / (ad_r + ad_a)
        if vaf is None:
            for k in ("VAF", "AF"):
                if k in fmt and s.get(k) is not None:
                    v = s.get(k)
                    vaf = v[0] if isinstance(v, tuple) else v
                    break
        if "GQ" in fmt:
            gq = s.get("GQ")
    if dp is None and "DP" in rec.info:
        dp = rec.info["DP"]
    mq = None
    for k in ("MQ", "MQM"):
        if k in rec.info:
            v = rec.info[k]
            mq = v[0] if isinstance(v, tuple) else v
            break
    return {"gt": "/".join(str(a) for a in alleles) if alleles else "./.",
            "has_alt": bool(alleles) and max(alleles) > 0,
            "qual": _num(rec.qual), "filter": ";".join(fl) if fl else ".",
            "is_pass": fl == ["PASS"] or not fl,
            "dp": dp, "ad_ref": ad_r, "ad_alt": ad_a, "vaf": _num(vaf, 3), "gq": gq, "mq": _num(mq, 1)}


def _key(rec):
    alt = rec.alts[0] if rec.alts else None
    if alt is None or alt.startswith("<") or alt == "*":
        return None
    return (rec.chrom, rec.pos, rec.ref.upper(), alt.upper())


def load_query(vcf: Path, all_records: bool) -> tuple[dict, dict]:
    recs: dict = {}
    by_chrom: dict = defaultdict(list)
    with pysam.VariantFile(str(vcf)) as fh:
        for rec in fh:
            k = _key(rec)
            if k is None or k in recs:
                continue
            info = record_info(rec)
            if not all_records and not info["is_pass"]:
                continue
            recs[k] = info
            by_chrom[k[0]].append((k[1], k))
    for v in by_chrom.values():
        v.sort()
    return recs, by_chrom


def _nearby(by_chrom: dict, k, window: int) -> list:
    lst = by_chrom.get(k[0]) or []
    i = bisect.bisect_left(lst, (k[1] - window,))
    out = []
    while i < len(lst) and lst[i][0] <= k[1] + window:
        out.append(lst[i][1])
        i += 1
    return out


def _vcfeval_matches(path: Path) -> set:
    keys = set()
    if path.exists():
        with pysam.VariantFile(str(path)) as fh:
            for rec in fh:
                for alt in rec.alts or ():
                    keys.add((rec.chrom, rec.pos, rec.ref.upper(), alt.upper()))
    return keys


# ------------------------------------------------------------------------------ the check
def classify(baseline: Path, query: dict, by_chrom: dict, *, baseline_pass_only: bool, window: int,
             haplotype_matches: set, check_cancel=lambda: None) -> tuple[list[dict], Counter]:
    rows, skipped = [], Counter()
    with pysam.VariantFile(str(baseline)) as fh:
        for i, rec in enumerate(fh):
            if i % 20000 == 0:
                check_cancel()
            k = _key(rec)
            if k is None:
                skipped["symbolic / no ALT"] += 1
                continue
            b = record_info(rec)
            if baseline_pass_only and not b["is_pass"]:
                skipped["baseline not PASS"] += 1
                continue
            if not b["has_alt"]:
                skipped["baseline genotype has no ALT"] += 1
                continue
            q = query.get(k)
            near = ""
            if q is not None and q["has_alt"]:
                status = "found_pass" if q["is_pass"] else "found_filtered"
                gt_match = "same" if q["gt"] == b["gt"] else "different"
            elif q is not None:
                status, gt_match = "found_ref_gt", "different"
            elif k in haplotype_matches:
                status, gt_match = "found_other_repr", ""
            else:
                nb = [n for n in _nearby(by_chrom, k, window) if query[n]["has_alt"]]
                status = "nearby" if nb else "not_called"
                gt_match = ""
                near = ",".join(f"{n[1]}:{n[2]}>{n[3]}({query[n]['filter']})" for n in nb[:5])
            q = q or {}
            rows.append({
                "CHROM": k[0], "POS": k[1], "REF": k[2], "ALT": k[3], "TYPE": engine.variant_type(k[2], k[3]),
                "STATUS": status, "GT_MATCH": gt_match,
                "B_GT": b["gt"], "B_QUAL": b["qual"], "B_FILTER": b["filter"], "B_DP": b["dp"],
                "B_AD_REF": b["ad_ref"], "B_AD_ALT": b["ad_alt"], "B_VAF": b["vaf"], "B_GQ": b["gq"], "B_MQ": b["mq"],
                "Q_GT": q.get("gt"), "Q_QUAL": q.get("qual"), "Q_FILTER": q.get("filter"), "Q_DP": q.get("dp"),
                "Q_AD_REF": q.get("ad_ref"), "Q_AD_ALT": q.get("ad_alt"), "Q_VAF": q.get("vaf"),
                "Q_GQ": q.get("gq"), "Q_MQ": q.get("mq"), "NEARBY": near})
    return rows, skipped


def what_if(rows: list[dict], *, min_dp=None, min_gq=None, min_vaf=None, min_qual=None,
            ignore_our_filters=False) -> dict:
    """How many baseline variants would count as found if our VCF were filtered with these thresholds.
    A variant is 'found' when our VCF has the same allele with an ALT genotype that meets every
    given threshold (and PASSes our own filters unless ignore_our_filters)."""
    def ok(r):
        if r["STATUS"] not in ("found_pass", "found_filtered"):
            return False
        if not ignore_our_filters and r["STATUS"] != "found_pass":
            return False
        for col, thr in (("Q_DP", min_dp), ("Q_GQ", min_gq), ("Q_VAF", min_vaf), ("Q_QUAL", min_qual)):
            if thr is not None:
                v = r.get(col)
                if v in (None, "") or float(v) < thr:
                    return False
        return True
    n = len(rows)
    found = sum(1 for r in rows if ok(r))
    return {"total": n, "found": found, "missing": n - found,
            "found_pct": round(100 * found / n, 2) if n else None}


def write_rows(path: Path, rows: list[dict], cols: list[str]) -> None:
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")


LOCUS_RE = re.compile(
    r"^\s*(?P<chrom>[A-Za-z0-9_.]+?)\s*[:\s\-_]\s*(?P<start>[\d,]+)"
    r"(?:\s*-\s*(?P<end>[\d,]+)(?=\s|$))?"
    r"(?:\s*[:\s\-_]?\s*(?P<ref>[ACGTNacgtn*]+)\s*(?:>|:|/|-|_|\s)\s*(?P<alt>[ACGTNacgtn*,]+))?\s*$")


def parse_locus(text: str) -> dict:
    """Parse a variant search: 'chr1:100000 A>G', 'chr1:100000', 'chr1:99,000-101,000', '1-100000-A-G',
    'chr1 100000 A G'. Returns {chrom, start, end, ref, alt} (ref/alt None when not given)."""
    m = LOCUS_RE.match(text or "")
    if not m:
        raise ValueError("Use chrom:position, chrom:start-end or chrom:position REF>ALT "
                         "(e.g. chr1:100000 A>G).")
    start = int(m["start"].replace(",", ""))
    end = int(m["end"].replace(",", "")) if m["end"] else start
    if end < start:
        raise ValueError("End position is before the start position.")
    return {"chrom": m["chrom"], "start": start, "end": end,
            "ref": m["ref"].upper() if m["ref"] else None, "alt": m["alt"].upper() if m["alt"] else None}


def match_locus(rows: list[dict], locus: dict, window: int = 0) -> list[dict]:
    """Rows at the locus (± window bp). Chromosome names are matched with or without 'chr'.
    REF/ALT, when given, must match exactly (after normalisation both files are left-aligned)."""
    c = locus["chrom"]
    names = {c, c[3:] if c.lower().startswith("chr") else "chr" + c}
    lo, hi = locus["start"] - window, locus["end"] + window
    out = []
    for r in rows:
        if r["CHROM"] not in names:
            continue
        pos = int(r["POS"])
        if not lo <= pos <= hi:
            continue
        if locus["ref"] and r["REF"].upper() != locus["ref"]:
            continue
        if locus["alt"] and locus["alt"] not in r["ALT"].upper().split(","):
            continue
        out.append(r)
    return out


def read_rows(path: Path) -> list[dict]:
    with gzip.open(path, "rt") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(cols, line.rstrip("\n").split("\t"))) for line in fh]


def _bins(values, edges):
    out = Counter()
    for v in values:
        if v is None:
            out["unknown"] += 1
            continue
        label = f"≥ {edges[-1]}"
        for lo, hi in zip(edges, edges[1:]):
            if lo <= v < hi:
                label = f"{lo}–{hi - 1 if isinstance(hi, int) else hi}"
                break
        out[label] += 1
    order = [f"{lo}–{hi - 1 if isinstance(hi, int) else hi}" for lo, hi in zip(edges, edges[1:])] + \
            [f"≥ {edges[-1]}", "unknown"]
    return [(k, out[k]) for k in order if out[k]]


def summarise(rows: list[dict], skipped: Counter, query_only: int | None, meta: dict) -> dict:
    by_status = Counter(r["STATUS"] for r in rows)
    by_type = defaultdict(Counter)
    for r in rows:
        by_type[r["TYPE"]][r["STATUS"]] += 1
    filt_combo = Counter(r["Q_FILTER"] for r in rows if r["STATUS"] == "found_filtered")
    filt_label = Counter()
    for combo, n in filt_combo.items():
        for lab in (combo or "").split(";"):
            filt_label[lab] += n
    n = len(rows)
    found_any = sum(by_status[s] for s in FOUND)
    missing = [r for r in rows if r["STATUS"] in ("not_called", "nearby")]
    return {
        **meta, "baseline_variants": n, "skipped": dict(skipped), "by_status": dict(by_status),
        "status_labels": STATUS_LABELS,
        "pct_found_pass": round(100 * by_status["found_pass"] / n, 2) if n else None,
        "pct_found_any": round(100 * found_any / n, 2) if n else None,
        "gt_different": sum(1 for r in rows if r["GT_MATCH"] == "different" and r["STATUS"] != "found_ref_gt"),
        "by_type": {t: dict(c) for t, c in by_type.items()},
        "filtered_by_combination": dict(filt_combo.most_common()),
        "filtered_by_label": dict(filt_label.most_common()),
        "missing_baseline_qual_bins": _bins([r["B_QUAL"] for r in missing], [0, 10, 30, 100]),
        "missing_baseline_dp_bins": _bins([r["B_DP"] for r in missing], [0, 5, 10, 20, 50]),
        "missing_baseline_mq_bins": _bins([r.get("B_MQ") for r in missing], [0, 10, 40, 60]),
        "missing_by_chrom": dict(Counter(r["CHROM"] for r in missing).most_common(10)),
        "query_only_pass": query_only,
    }


class VariantCheckStep(Step):
    key = "variant_check"
    name = "Variant-by-variant check (baseline → our VCF)"
    description = ("For every variant of the baseline VCF: is it in our VCF (PASS, filtered, genotyped as "
                   "reference, written differently, a different variant nearby, or not called)? Side-by-side "
                   "QUAL/DP/AD/VAF/GQ of both files.")
    tools = ("bcftools", "rtg")

    def execute(self, ctx: StepContext, runner) -> StepResult:
        inputs = read_json(Path(prior_output(ctx, "vcf_prepare", "inputs_json")))
        if len(inputs) != 2:
            raise StepFailure(f"Exactly two VCFs are needed (baseline + ours); got {len(inputs)}.")
        base = self._baseline(ctx, inputs)
        query = next(i for i in inputs if i is not base)
        out = ctx.step_dir("comparison/variant_check")
        p = ctx.params
        warnings = []
        ctx.log(f"Baseline: {base['source']} ({base['origin']}); our VCF: {query['source']} ({query['origin']})")

        q_recs, by_chrom = load_query(Path(query["vcf"]), bool(p.get("query_all_records", True)))
        ctx.log(f"Our VCF: {len(q_recs):,} records loaded "
                f"({'all records, filtered ones included' if p.get('query_all_records', True) else 'PASS only'}).")

        hap = set()
        if p.get("check_haplotype", True):
            if shutil.which("rtg"):
                sdf = engine.ensure_sdf(runner, ctx.reference, ctx.log)
                vdir = out / "vcfeval"
                if vdir.exists():                     # resume: vcfeval refuses an existing output folder
                    shutil.move(str(vdir), str(out / f"vcfeval.previous.{vdir.stat().st_mtime_ns}"))
                engine.run_vcfeval(runner, sdf, Path(base["vcf"]), Path(query["vcf"]), vdir, regions=None,
                                   threads=min(ctx.threads, 8), squash_ploidy=True, all_records=True)
                hap = _vcfeval_matches(vdir / "tp-baseline.vcf.gz")
            else:
                warnings.append("RTG Tools not found: representation differences are not rescued "
                                "(they appear as 'nearby' or 'not called').")

        rows, skipped = classify(Path(base["vcf"]), q_recs, by_chrom,
                                 baseline_pass_only=bool(p.get("baseline_pass_only", True)),
                                 window=int(p.get("nearby_window", 10)), haplotype_matches=hap,
                                 check_cancel=ctx.check_cancel)
        if not rows:
            warnings.append("The baseline VCF has no variants in the analysis regions (after the PASS/genotype "
                            "selection).")
        # Identifiable names: <our VCF>_vs_<baseline VCF>.<what>, e.g.
        # SAMPLE_A.provider-like-dragen-4-3.gatk_hc.filtered_vs_SAMPLE_A_provider.hard-filtered.variant_check.tsv.gz
        name = f"{vcf_stem(query['source'])}_vs_{vcf_stem(base['source'])}"
        ctx.log(f"Output file prefix: {name}")
        tsv = out / f"{name}.variant_check.tsv.gz"
        write_rows(tsv, rows, COLUMNS)
        outputs = {"variant_tsv": str(tsv)}

        n_qonly = None
        if p.get("report_query_only", True):
            bkeys = {(r["CHROM"], int(r["POS"]), r["REF"], r["ALT"]) for r in rows}
            tp_call = _vcfeval_matches(out / "vcfeval" / "tp.vcf.gz") if hap else set()
            qonly = []
            for k, q in q_recs.items():
                if q["is_pass"] and q["has_alt"] and k not in bkeys and k not in tp_call:
                    qonly.append({"CHROM": k[0], "POS": k[1], "REF": k[2], "ALT": k[3],
                                  "TYPE": engine.variant_type(k[2], k[3]),
                                  **{"Q_" + c.upper(): q[c] for c in ("gt", "qual", "filter", "dp", "ad_ref", "mq",
                                                                      "ad_alt", "vaf", "gq")}})
            qonly.sort(key=lambda r: (r["CHROM"], r["POS"]))
            qtsv = out / f"{name}.ours_pass_not_in_baseline.tsv.gz"
            write_rows(qtsv, qonly, Q_COLUMNS)
            outputs["query_only_tsv"] = str(qtsv)
            n_qonly = len(qonly)

        meta = {"baseline": {"label": base["label"], "source": base["source"], "origin": base["origin"]},
                "query": {"label": query["label"], "source": query["source"], "origin": query["origin"]},
                "baseline_pass_only": bool(p.get("baseline_pass_only", True)),
                "query_all_records": bool(p.get("query_all_records", True)),
                "nearby_window": int(p.get("nearby_window", 10)), "haplotype_check": bool(hap) or None,
                "regions": prior_output(ctx, "regions", "calling_bed", required=False)}
        summary = summarise(rows, skipped, n_qonly, meta)
        summary["what_if_examples"] = [
            {"label": "As our VCF is (our FILTER column)", **what_if(rows)},
            {"label": "Ignore our filters (any record with the allele)", **what_if(rows, ignore_our_filters=True)},
            {"label": "Ignore our filters, depth ≥ 2 (DRAGEN-like)", **what_if(rows, ignore_our_filters=True, min_dp=2)},
            {"label": "Ignore our filters, depth ≥ 10", **what_if(rows, ignore_our_filters=True, min_dp=10)},
        ]
        write_json(out / f"{name}.variant_check_summary.json", summary)

        from flask import render_template

        from ...models.base import utcnow
        from ...services.limitations import run_limitations
        examples = {s: [r for r in rows if r["STATUS"] == s][:150] for s in ("not_called", "nearby", "found_filtered")}
        html = out / f"{name}.variant_check_report.html"
        html.write_text(render_template("reports/variant_check_report.html", s=summary, examples=examples,
                                        run_id=ctx.run.run_id, sample=ctx.sample.sample_id, file_prefix=name,
                                        reference=ctx.reference.name, generated_at=utcnow().isoformat() + "Z",
                                        limitations=run_limitations(ctx.reference, ctx.params)))
        outputs["report_html"] = str(html)
        outputs["summary_json"] = str(out / f"{name}.variant_check_summary.json")
        return StepResult(outputs=outputs, summary=summary, warnings=warnings)

    def validate(self, ctx, result):
        for name in ("variant_tsv", "report_html", "summary_json"):
            p = Path(result.outputs.get(name, ""))
            if not p.is_file() or p.stat().st_size == 0:
                raise StepFailure(f"Expected output '{name}' missing or empty: {p}")

    @staticmethod
    def _baseline(ctx, inputs):
        fid = ctx.params.get("baseline_file_id")
        if fid:
            f = next((f for f in ctx.input_files if f.id == fid), None)
            if f is None:
                raise StepFailure("The chosen baseline VCF is not one of the selected input files.")
            b = next((i for i in inputs if i["source"] == f.path), None)
            if b is None:
                raise StepFailure("The chosen baseline VCF was not prepared (internal error).")
            return b
        ext = [i for i in inputs if i["origin"] == "external"]
        if len(ext) == 1:
            return ext[0]
        raise StepFailure("Choose which VCF is the baseline (exactly one 'External provider' VCF, or select it "
                          "on the form).")
