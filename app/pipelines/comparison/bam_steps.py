"""BAM vs BAM comparison.

The same QC is computed for every BAM, over the SAME regions with the SAME
settings (otherwise numbers are not comparable):

* samtools flagstat  - read counts by category (primary, secondary, supplementary,
                       duplicates, properly paired, singletons, mate on another chromosome)
* samtools stats     - insert size, error rate, MAPQ distribution (analysis contigs only)
* samtools idxstats  - mapped reads per contig (shows e.g. reads on decoy/ALT/HLA contigs)
* mosdepth           - target coverage at configurable thresholds, per target region

Then the per-target coverage is compared: targets that reach the threshold in one
BAM but not in the other are listed — this is where alignment differences
(reference version, ALT handling, MAPQ) become visible.
"""
from __future__ import annotations

import gzip
from pathlib import Path

from flask import render_template

from ...engine.step import Step, StepContext, StepFailure, StepResult
from ...models.base import utcnow
from ...services.reference_service import reference_contigs
from ..bam_processing.steps import prepare_bam
from ..common import parse_contig_selection, parse_flagstat_json, parse_mosdepth, prior_output, write_json
from .engine import _csv
from ...services.limitations import run_limitations
from .steps import ORIGIN_LABELS, _file_label

MAPQ_BINS = [(0, 0, "0"), (1, 9, "1-9"), (10, 19, "10-19"), (20, 39, "20-39"), (40, 59, "40-59"), (60, 255, "60+")]


def parse_thresholds(text: str) -> list[int]:
    try:
        vals = sorted({int(x) for x in str(text).split(",") if x.strip()})
    except ValueError:
        raise StepFailure(f"Coverage thresholds must be comma-separated integers: {text!r}") from None
    if not vals or vals[0] < 1 or len(vals) > 12:
        raise StepFailure("Give 1-12 coverage thresholds ≥ 1, e.g. 1,10,20,30,50,100.")
    return vals


def _samtools_stats(path: Path) -> tuple[dict, dict]:
    sn, mapq = {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("SN\t"):
                k, v = line.split("\t")[1:3]
                sn[k.rstrip(":")] = v
            elif line.startswith("MAPQ\t"):
                _, q, n = line.rstrip("\n").split("\t")[:3]
                mapq[int(q)] = int(n)
    return sn, mapq


def _mapq_bins(mapq: dict) -> dict:
    total = sum(mapq.values()) or 1
    out = {}
    for lo, hi, label in MAPQ_BINS:
        n = sum(c for q, c in mapq.items() if lo <= q <= hi)
        out[label] = {"reads": n, "pct": round(100 * n / total, 2)}
    return out


def _idxstats(path: Path) -> dict:
    groups = {"primary": 0, "alt": 0, "decoy": 0, "hla": 0, "random/unplaced": 0, "other": 0, "unmapped": 0}
    per = []
    with open(path) as fh:
        for line in fh:
            name, length, mapped, unmapped = line.rstrip("\n").split("\t")
            mapped, unmapped = int(mapped), int(unmapped)
            if name == "*":
                groups["unmapped"] += unmapped
                continue
            n = name.lower()
            g = ("alt" if n.endswith("_alt") else "decoy" if "decoy" in n else "hla" if n.startswith("hla")
                 else "random/unplaced" if ("random" in n or n.startswith("chrun")) else
                 "primary" if n.replace("chr", "") in [str(i) for i in range(1, 23)] + ["x", "y", "m", "mt"]
                 else "other")
            groups[g] += mapped
            per.append({"contig": name, "length": int(length), "mapped": mapped})
    return {"groups": groups, "per_contig": per}


def _region_means(path: Path) -> dict:
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            out[(p[0], int(p[1]), int(p[2]))] = float(p[-1])
    return out


class BamMetricsStep(Step):
    key = "bam_metrics"
    name = "BAM metrics (each BAM)"
    description = ("For every selected BAM: integrity/build check, flagstat, stats (MAPQ, insert size, error "
                   "rate), idxstats and target coverage — same regions and settings for all.")
    tools = ("samtools", "mosdepth")

    def execute(self, ctx: StepContext, runner) -> StepResult:
        files = ctx.files("bam")
        if len(files) < 2:
            raise StepFailure(f"Select at least 2 BAM files to compare (selected: {len(files)}).")
        thresholds = parse_thresholds(ctx.params.get("coverage_thresholds", "1,10,20,30,50,100"))
        ref_contigs = reference_contigs(ctx.reference)
        contigs = parse_contig_selection(ctx.params.get("contigs", "primary"), ref_contigs)
        targets = prior_output(ctx, "regions", "targets_bed", required=False)
        by = targets or prior_output(ctx, "regions", "calling_bed")
        t = str(min(ctx.threads, 4))
        used: set = set()
        results, warnings = [], []
        for f in files:
            label = _file_label(f, used)
            d = ctx.step_dir(f"qc/compare/{label}")
            ctx.log(f"--- {label}: {f.path}")
            bam, _, prep, w = prepare_bam(ctx, runner, f, d)
            warnings += [f"{label}: {x}" for x in w]
            fs = d / "flagstat.json"
            runner.run(["samtools", "flagstat", "-@", t, "-O", "json", bam], stdout_path=fs)
            st = d / "samtools_stats.txt"
            runner.run(["samtools", "stats", "-@", t, "-r", ctx.reference.fasta_path, bam, *contigs],
                       stdout_path=st)
            ix = d / "idxstats.tsv"
            runner.run(["samtools", "idxstats", "-@", t, bam], stdout_path=ix)
            prefix = d / label
            runner.run(["mosdepth", "-t", t, "-n", "-x", "-Q", str(ctx.params.get("coverage_min_mapq", 20)),
                        "--by", by, "-T", ",".join(map(str, thresholds)), prefix, bam])
            sn, mapq = _samtools_stats(st)
            results.append({
                "label": label, "source": f.path, "origin": f.origin, "prep": prep,
                "flagstat": parse_flagstat_json(fs),
                "stats": {"insert_size_mean": float(sn.get("insert size average", 0)),
                          "insert_size_sd": float(sn.get("insert size standard deviation", 0)),
                          "error_rate": float(sn.get("error rate", 0)),
                          "average_length": float(sn.get("average length", 0)),
                          "reads_mq0": int(sn.get("reads MQ0", 0))},
                "mapq": _mapq_bins(mapq),
                "idx_groups": _idxstats(ix)["groups"],
                "coverage": parse_mosdepth(prefix, thresholds),
                "files": {"flagstat": str(fs), "stats": str(st), "idxstats": str(ix),
                          "regions": f"{prefix}.regions.bed.gz", "summary": f"{prefix}.mosdepth.summary.txt"},
            })
            ctx.check_cancel()
        out = ctx.step_dir("comparison/bam")
        path = write_json(out / "bam_metrics.json", {"thresholds": thresholds, "bams": results,
                                                       "coverage_region": "targets" if targets else "calling regions"})
        return StepResult(outputs={"metrics_json": path},
                          summary={"bams": [r["label"] for r in results], "thresholds": thresholds},
                          warnings=warnings)


METRIC_ROWS = [
    ("File size (GB)", lambda r: round(r["prep"]["size_bytes"] / 1e9, 2)),
    ("Total records", lambda r: r["flagstat"]["total_records"]),
    ("Primary reads", lambda r: r["flagstat"]["primary"]),
    ("Secondary", lambda r: r["flagstat"]["secondary"]),
    ("Supplementary", lambda r: r["flagstat"]["supplementary"]),
    ("Mapped (all records)", lambda r: r["flagstat"]["mapped"]),
    ("Primary mapped %", lambda r: r["flagstat"]["primary_mapped_pct"]),
    ("Duplicates (primary)", lambda r: r["flagstat"]["duplicates"]),
    ("Duplicate %", lambda r: r["flagstat"].get("duplicate_pct")),
    ("Properly paired %", lambda r: r["flagstat"]["properly_paired_pct"]),
    ("Singletons", lambda r: r["flagstat"]["singletons"]),
    ("Mate on other chromosome", lambda r: r["flagstat"]["mate_diff_chr"]),
    ("  … with MAPQ ≥ 5", lambda r: r["flagstat"]["mate_diff_chr_mapq5"]),
    ("Reads MAPQ 0 (analysis contigs)", lambda r: r["stats"]["reads_mq0"]),
    ("Insert size mean", lambda r: r["stats"]["insert_size_mean"]),
    ("Error rate", lambda r: r["stats"]["error_rate"]),
    ("Average read length", lambda r: r["stats"]["average_length"]),
    ("Mean target coverage", lambda r: r["coverage"].get("mean_target_coverage")),
]


class BamCompareStep(Step):
    key = "bam_compare"
    name = "BAM comparison"
    description = "Side-by-side QC table, MAPQ distribution, reads per contig class, per-target coverage differences."

    def execute(self, ctx: StepContext, runner) -> StepResult:
        from ..common import read_json
        m = read_json(Path(prior_output(ctx, "bam_metrics", "metrics_json")))
        bams, thr = m["bams"], m["thresholds"]
        labels = [b["label"] for b in bams]
        rows = []
        for name, fn in METRIC_ROWS:
            vals = [fn(b) for b in bams]
            rows.append({"metric": name, **dict(zip(labels, vals))})
        for t in thr:
            rows.append({"metric": f"% target bases ≥ {t}x",
                         **{b["label"]: (b["coverage"].get("pct_target_bases") or {}).get(f"{t}x") for b in bams}})
        for _, _, lab in MAPQ_BINS:
            rows.append({"metric": f"MAPQ {lab} (% reads)", **{b["label"]: b["mapq"][lab]["pct"] for b in bams}})
        for g in bams[0]["idx_groups"]:
            rows.append({"metric": f"Reads on {g} contigs", **{b["label"]: b["idx_groups"][g] for b in bams}})
        out = ctx.step_dir("comparison/bam")
        files = {"metrics_csv": _csv(out / "bam_comparison.csv", rows)}

        # per-target: covered (mean >= threshold) in one BAM but not the other (first two BAMs)
        diff_thr = int(ctx.params.get("target_diff_threshold", 20))
        a, b = bams[0], bams[1]
        ra, rb = _region_means(Path(a["files"]["regions"])), _region_means(Path(b["files"]["regions"]))
        only_a, only_b, diffs = 0, 0, []
        for reg in ra.keys() & rb.keys():
            x, y = ra[reg], rb[reg]
            if (x >= diff_thr) != (y >= diff_thr):
                only_a += x >= diff_thr
                only_b += y >= diff_thr
                diffs.append({"chrom": reg[0], "start": reg[1], "end": reg[2], a["label"]: x, b["label"]: y})
        diffs.sort(key=lambda r: (r["chrom"], r["start"]))
        files["target_diff_tsv"] = _csv(out / "targets_differently_covered.csv", diffs,
                                        ["chrom", "start", "end", a["label"], b["label"]])
        by_chrom: dict = {}
        for r in diffs:
            by_chrom[r["chrom"]] = by_chrom.get(r["chrom"], 0) + 1
        summary = {
            "bams": [{"label": x["label"], "source": x["source"], "origin": x["origin"],
                      "build": x["prep"].get("build"), "verdict": x["prep"].get("verdict")} for x in bams],
            "table": rows, "thresholds": thr, "coverage_region": m["coverage_region"],
            "target_diff": {"threshold": diff_thr, "a": a["label"], "b": b["label"],
                            "covered_only_a": only_a, "covered_only_b": only_b,
                            "targets_compared": len(ra.keys() & rb.keys()),
                            "by_chrom": dict(sorted(by_chrom.items(), key=lambda x: -x[1])[:15])},
        }
        html = out / "bam_comparison_report.html"
        html.write_text(render_template("reports/bam_comparison_report.html", s=summary, labels=labels,
                                        run_id=ctx.run.run_id, sample=ctx.sample.sample_id,
                                        reference=ctx.reference.name, files=files,
                                        origin_labels=ORIGIN_LABELS, top_diffs=diffs[:200],
                                        limitations=run_limitations(ctx.reference, ctx.params),
                                        generated_at=utcnow().isoformat() + "Z"))
        files["report_html"] = str(html)
        write_json(out / "bam_comparison_summary.json", summary)
        return StepResult(outputs=files, summary=summary)
