"""Helpers shared by the FASTQ, BAM and variant-calling pipelines."""
from __future__ import annotations

import gzip
import json
import os
import re
import subprocess
from pathlib import Path

from ..engine.step import StepContext, StepFailure
from ..models import status as S
from ..services.genome_build import is_primary
from ..services.reference_service import reference_contigs

FLOWCELL_HEADER_RE = re.compile(
    r"^@(?P<instrument>[^:\s]+):(?P<run>\d+):(?P<flowcell>[^:\s]+):(?P<lane>\d+):")
RG_SAFE_RE = re.compile(r"[^A-Za-z0-9_.-]")


SLUG_RE = re.compile(r"[^A-Za-z0-9]+")


def slug(text: str, max_len: int = 40) -> str:
    """'Provider-like (DRAGEN 4.3)' -> 'provider-like-dragen-4-3' (safe in a file name)."""
    s = SLUG_RE.sub("-", text or "").strip("-").lower()[:max_len].strip("-")
    return s or "x"


def preset_tag(params: dict | None) -> str:
    """File-name tag of the preset a run was started from ('' when none).
    '-modified' is added when fields were changed after applying the preset, so a file
    name never claims a preset whose values were not used as they are."""
    pr = (params or {}).get("preset") or {}
    if not pr.get("name"):
        return ""
    return slug(pr["name"]) + ("-modified" if pr.get("changed_after_applying") else "")


def input_bam_tag(ctx: StepContext) -> str:
    """Label tag of an input BAM that was NOT produced by this platform (e.g. a provider's DRAGEN
    BAM), so its calls are never confused with calls from our own alignment. '' otherwise."""
    for f in getattr(ctx, "input_files", None) or []:
        if getattr(f, "role", None) == "bam" and getattr(f, "origin", "pipeline") != "pipeline":
            return slug(f.label or Path(f.path).stem, 30) + "-bam"
    return ""


def output_stem(ctx: StepContext) -> str:
    """Identifiable output name: <sample>[.<input BAM>][.<preset>], e.g.
    SAMPLE_A.provider-like-dragen-4-3 or SAMPLE_A.dragen-4-3-provider-bam.raw-tool-defaults-no-filters.
    Read groups / VCF sample names keep the plain sample ID."""
    parts = [ctx.sample.sample_id, input_bam_tag(ctx), preset_tag(ctx.params)]
    return ".".join(x for x in parts if x)


def vcf_stem(path) -> str:
    """File name without .vcf / .vcf.gz."""
    name = Path(path).name
    for ext in (".vcf.gz", ".vcf.bgz", ".vcf"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def write_json(path: Path, data) -> str:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, default=str))
    os.replace(tmp, path)
    return str(path)


def read_json(path) -> dict:
    with open(path) as fh:
        return json.load(fh)


def prior_output(ctx: StepContext, step_key: str, name: str, required: bool = True) -> str | None:
    """Output of an earlier step of this run (works across resume, read from the DB)."""
    for sr in ctx.run.steps:
        if sr.key == step_key and sr.status == S.COMPLETED:
            value = (sr.outputs or {}).get(name)
            if value:
                return value
    if required:
        raise StepFailure(f"Required output '{name}' of step '{step_key}' is not available — "
                          "that step did not complete.")
    return None


def prior_summary(ctx: StepContext, step_key: str) -> dict:
    for sr in ctx.run.steps:
        if sr.key == step_key and sr.status == S.COMPLETED:
            return sr.summary or {}
    return {}


# --------------------------------------------------------------------------- read groups
def read_group_from_fastq(fastq: Path, sample_id: str, library: str | None, platform: str = "ILLUMINA") -> dict:
    """Build the @RG fields from the first Illumina read header.

    ID = flowcell.lane identifies the sequencing *lane* (duplicate marking and
    base-quality models are per read group). SM = our sample ID, LB = library.
    """
    opener = gzip.open if fastq.name.endswith(".gz") else open
    with opener(fastq, "rt", errors="replace") as fh:
        header = fh.readline().strip()
    m = FLOWCELL_HEADER_RE.match(header)
    if m:
        flowcell, lane = m["flowcell"], m["lane"]
        rg_id = f"{flowcell}.{lane}"
        pu = f"{flowcell}.{lane}.{sample_id}"
    else:
        rg_id, pu = f"{sample_id}.1", f"{sample_id}.1"
    rg = {"ID": rg_id, "SM": sample_id, "LB": library or sample_id, "PL": platform, "PU": pu}
    return {k: RG_SAFE_RE.sub("_", v) for k, v in rg.items()}


def rg_line(rg: dict) -> str:
    # bwa interprets the two-character sequence "\t" as a tab
    return "@RG\\t" + "\\t".join(f"{k}:{rg[k]}" for k in ("ID", "SM", "LB", "PL", "PU") if k in rg)


# --------------------------------------------------------------------------- regions
def parse_contig_selection(value: str, ref_contigs: dict[str, int]) -> list[str]:
    """'primary' -> chr1..22,X,Y,M (reference order); else a comma list validated against the .fai."""
    value = (value or "primary").strip()
    if value == "primary":
        return [c for c in ref_contigs if is_primary(c)]
    if value == "all":
        return list(ref_contigs)
    wanted = [c.strip() for c in value.split(",") if c.strip()]
    missing = [c for c in wanted if c not in ref_contigs]
    if missing:
        raise StepFailure(f"Contig(s) not in the reference: {missing}. "
                          f"Reference uses names like {next(iter(ref_contigs))!r}.")
    order = {c: i for i, c in enumerate(ref_contigs)}
    return sorted(set(wanted), key=order.get)


def read_bed(path: Path) -> list[tuple[str, int, int]]:
    intervals = []
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise StepFailure(f"{path.name} line {n}: BED needs at least 3 tab-separated columns.")
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                raise StepFailure(f"{path.name} line {n}: start/end are not integers.") from None
            if start < 0 or end <= start:
                raise StepFailure(f"{path.name} line {n}: invalid interval {start}-{end}.")
            intervals.append((parts[0], start, end))
    return intervals


def merge_intervals(intervals, ref_contigs: dict[str, int], padding: int = 0):
    """Pad, clip to contig ends, sort in reference order and merge overlaps."""
    order = {c: i for i, c in enumerate(ref_contigs)}
    padded = sorted(
        ((c, max(0, s - padding), min(ref_contigs[c], e + padding)) for c, s, e in intervals),
        key=lambda iv: (order[iv[0]], iv[1]),
    )
    merged: list[list] = []
    for c, s, e in padded:
        if merged and merged[-1][0] == c and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([c, s, e])
    return [tuple(m) for m in merged]


def write_bed(path: Path, intervals) -> str:
    with open(path, "w") as fh:
        for c, s, e in intervals:
            fh.write(f"{c}\t{s}\t{e}\n")
    return str(path)


def bed_size(intervals) -> int:
    return sum(e - s for _, s, e in intervals)


# --------------------------------------------------------------------------- validation helpers
def samtools_quickcheck(path: Path) -> None:
    try:
        proc = subprocess.run(["samtools", "quickcheck", "-v", str(path)], capture_output=True, text=True,
                              errors="replace", timeout=300)
    except FileNotFoundError as exc:
        raise StepFailure("samtools not found for output validation") from exc
    if proc.returncode != 0:
        raise StepFailure(f"BAM failed samtools quickcheck (truncated/corrupt): {path}")


def check_vcf(path: Path) -> None:
    proc = subprocess.run(["bcftools", "view", "-h", str(path)], capture_output=True, text=True,
                          errors="replace", timeout=300)
    if proc.returncode != 0 or "#CHROM" not in proc.stdout:
        raise StepFailure(f"VCF header could not be read: {path}")
    if not Path(str(path) + ".tbi").exists() and not Path(str(path) + ".csi").exists():
        raise StepFailure(f"VCF index missing: {path}")


def check_fastq_gz(path: Path) -> None:
    try:
        with gzip.open(path, "rb") as fh:
            first = fh.readline()
    except (OSError, EOFError) as exc:
        raise StepFailure(f"{path.name} is not a readable gzip file: {exc}") from exc
    if first and not first.startswith(b"@"):
        raise StepFailure(f"{path.name} does not start with a FASTQ record.")


def check_index_fresh(bam: Path) -> None:
    bai = Path(str(bam) + ".bai")
    if not bai.exists():
        raise StepFailure(f"BAM index missing: {bai}")
    if bai.stat().st_mtime < bam.stat().st_mtime:
        raise StepFailure(f"BAM index is older than the BAM (stale): {bai}")


# --------------------------------------------------------------------------- tool output parsers
def parse_flagstat_json(path: Path) -> dict:
    data = read_json(path)["QC-passed reads"]
    total = data.get("total", 0)
    out = {
        "total_records": total,
        "primary": data.get("primary", 0),
        "secondary": data.get("secondary", 0),
        "supplementary": data.get("supplementary", 0),
        "duplicates": data.get("primary duplicates", data.get("duplicates", 0)),
        "primary_mapped": data.get("primary mapped", 0),
        "primary_mapped_pct": data.get("primary mapped %", 0),
        "properly_paired_pct": data.get("properly paired %", 0),
        "singletons_pct": data.get("singletons %", 0),
        "mapped": data.get("mapped", 0),
        "mapped_pct": data.get("mapped %", 0),
        "paired": data.get("paired in sequencing", 0),
        "properly_paired": data.get("properly paired", 0),
        "singletons": data.get("singletons", 0),
        "mate_diff_chr": data.get("with mate mapped to a different chr", 0),
        "mate_diff_chr_mapq5": data.get("with mate mapped to a different chr (mapQ >= 5)", 0),
    }
    if out["primary"]:
        out["duplicate_pct"] = round(100 * out["duplicates"] / out["primary"], 2)
    return out


def parse_bcftools_stats(path: Path) -> dict:
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "SN":
                key = parts[2].rstrip(":").replace("number of ", "").replace(" ", "_")
                out[key] = int(parts[3])
            elif parts[0] == "TSTV" and "ts_tv" not in out:
                out["ts"], out["tv"], out["ts_tv"] = int(parts[2]), int(parts[3]), float(parts[4])
    return out


def parse_mosdepth(prefix: Path, thresholds: list[int]) -> dict:
    summary = Path(f"{prefix}.mosdepth.summary.txt")
    out: dict = {}
    with open(summary) as fh:
        next(fh)
        for line in fh:
            chrom, length, bases, mean, mn, mx = line.split()[:6]
            # "total" = whole genome (meaningless for an exome); "total_region" = the --by regions
            if chrom == "total_region":
                out["mean_target_coverage"] = float(mean)
                out["target_bases"] = int(length)
    thr_file = Path(f"{prefix}.thresholds.bed.gz")
    if thr_file.exists():
        covered = [0] * len(thresholds)
        total = 0
        with gzip.open(thr_file, "rt") as fh:
            next(fh)
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                total += int(parts[2]) - int(parts[1])
                for i in range(len(thresholds)):
                    covered[i] += int(parts[4 + i])
        if total:
            out["pct_target_bases"] = {f"{t}x": round(100 * c / total, 2) for t, c in zip(thresholds, covered)}
    return out


def sum_size_gb(files, roles) -> float:
    return sum((f.size_bytes or 0) for f in files if f.role in roles) / 1e9
