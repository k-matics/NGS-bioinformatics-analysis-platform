"""Phase 4-5: filtering, DeepVariant, caller comparison, truth benchmark, VCF/BAM comparison, read evidence.

Uses the simulated mini-genome of test_germline.py, whose variants are KNOWN, so
comparison results can be checked against the truth (not just "a file exists").
"""
import gzip
import shutil
import subprocess
from pathlib import Path

import pysam
import pytest

from app.extensions import db
from app.models import SampleFile, WorkflowRun
from app.models import status as S
from app.pipelines.comparison.engine import compare_sets
from app.pipelines.comparison.evidence import fisher_two_sided, notes_for
from app.pipelines.variant_calling.callers import get_caller
from app.pipelines.variant_calling.callers.deepvariant import DEFAULT_IMAGE, _mounts, image_id
from app.services.workflow_service import WorkflowError, create_run
from test_germline import HET_DEL, HET_SNV1, HET_SNV2, HOM_SNV, TOOLS, _sample, mini, mini_ref  # noqa: F401

pytestmark = pytest.mark.skipif(any(shutil.which(t) is None for t in TOOLS + ("rtg",)),
                                reason="bioinformatics tools not on PATH")


def _reg(sample, role, origin, path, label=None):
    from app.services.file_service import file_fingerprint
    f = SampleFile(sample=sample, role=role, origin=origin, label=label, path=str(path),
                   **file_fingerprint(Path(path)))
    db.session.add(f)
    db.session.commit()
    return f


def _reload(run):
    db.session.expire_all()
    return db.session.get(WorkflowRun, run.id)


def _truth_vcf(mini, path: Path) -> Path:
    seq, alt = mini["chr22"], mini["alt"]
    rows = [(HOM_SNV, seq[HOM_SNV - 1], alt[HOM_SNV], "1/1"), (HET_SNV1, seq[HET_SNV1 - 1], alt[HET_SNV1], "0/1"),
            (HET_SNV2, seq[HET_SNV2 - 1], alt[HET_SNV2], "0/1"),
            (HET_DEL, seq[HET_DEL - 1:HET_DEL + 2], seq[HET_DEL - 1], "0/1")]
    txt = path.with_suffix("")
    with open(txt, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=5000>\n##contig=<ID=chr22,length=30000>\n"
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\n")
        for pos, r, a, gt in rows:
            fh.write(f"chr22\t{pos}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t{gt}\n")
    pysam.tabix_index(str(txt), preset="vcf", force=True)
    return Path(str(txt) + ".gz")


# --------------------------------------------------------------------------- unit tests
def test_fisher_exact_matches_known_values():
    assert abs(fisher_two_sided(3, 1, 1, 3) - 0.4857) < 1e-3
    assert abs(fisher_two_sided(10, 10, 10, 10) - 1.0) < 1e-9
    assert abs(fisher_two_sided(20, 20, 0, 15) - 3.4622e-4) < 1e-7   # all ALT reads on one strand


def test_compare_sets_counts_and_genotype_concordance():
    a = {("chr1", 10, "A", "G"): ("0/1", "het", "SNV"), ("chr1", 20, "AT", "A"): ("1/1", "hom", "indel"),
         ("chr1", 30, "C", "T"): ("0/1", "het", "SNV")}
    b = {("chr1", 10, "A", "G"): ("1/1", "hom", "SNV"), ("chr1", 20, "AT", "A"): ("1/1", "hom", "indel"),
         ("chr1", 40, "G", "A"): ("0/1", "het", "SNV")}
    r = compare_sets({"a": a, "b": b})
    p = r["pairwise"][0]
    assert (p["shared"], p["a_only"], p["b_only"]) == (2, 1, 1)
    assert p["gt_concordant"] == 1 and p["gt_discordant"] == 1 and p["gt_concordance_pct"] == 50.0
    assert r["in_all"] == 2 and r["union"] == 4
    assert {tuple(i["sets"]): i["total"] for i in r["intersections"]} == {("a", "b"): 2, ("a",): 1, ("b",): 1}


def test_evidence_notes_are_observations():
    ev = {"depth": 40, "ref_reads": 36, "alt_reads": 4, "other_reads": 0, "vaf": 0.1, "alt_fwd": 4, "alt_rev": 0,
          "strand_bias_p": 0.001, "mean_mq_alt": 60, "mq0_frac": 0, "alt_softclip_frac": 0,
          "alt_median_dist_to_end": 50}
    assert set(notes_for(ev, "het")) == {"vaf_het", "strand"}
    assert notes_for({**ev, "depth": 0}, "het") == ["no_reads"]


def test_filter_rules_are_caller_specific():
    gatk = [r.name for r in get_caller("gatk_hc").filter_rules({"filter_gatk_hard": True, "filter_min_depth": 10})]
    fb = [r.name for r in get_caller("freebayes").filter_rules({"filter_strand": True, "filter_min_qual": 20})]
    assert "HardFilterSNP" in gatk and "HardFilterIndel" in gatk and "StrandBias" not in gatk
    assert "StrandBias" in fb and "LowQual" in fb and "LowGQ" not in fb     # FreeBayes has no GQ


def test_deepvariant_mounts_inputs_read_only(tmp_path):
    ref = tmp_path / "ref" / "g.fa"
    ref.parent.mkdir()
    ref.write_text(">c\nA\n")
    run = tmp_path / "run"
    (run / "bam").mkdir(parents=True)
    bam = run / "bam" / "x.bam"
    bam.write_text("x")
    args = _mounts([ref, bam], [run])
    assert f"{ref.parent}:{ref.parent}:ro" in args and f"{run}:{run}" in args
    assert not any(a.startswith(f"{run / 'bam'}") for a in args)       # inside the writable mount


# --------------------------------------------------------------------------- integration
def test_filter_compare_truth_register_and_comparison_workflows(app, make_sample, mini, mini_ref, tmp_path):
    sample = _sample(make_sample, mini)
    truth = _truth_vcf(mini, tmp_path / "truth.vcf.gz")
    _reg(sample, "truth_vcf", "external", truth)
    truth_bed = tmp_path / "truth_confident.bed"
    shutil.copy(mini["bed"], truth_bed)
    _reg(sample, "truth_bed", "external", truth_bed)
    run = create_run(sample, "fastq_to_vcf", list(sample.files), mini_ref,
                     params={"contigs": "chr22", "callers": ["bcftools", "freebayes"], "retain_intermediates": False,
                             "caller_params": {"freebayes": {"min_coverage": 3}}})
    run = _reload(run)
    assert run.status == S.COMPLETED, run.error_message
    st = {s.key: s for s in run.steps}

    # --- soft filtering: nothing removed, FILTER set, planted variants PASS
    for k in ("bcftools", "freebayes"):
        f = st[f"filter_{k}"]
        assert f.status == S.COMPLETED
        with pysam.VariantFile(f.outputs["filtered_vcf"]) as fh:
            recs = list(fh)
            assert "NoAltGT" in fh.header.filters and "LowDP" in fh.header.filters
        n_norm = sum(1 for _ in pysam.VariantFile(st[f"call_{k}"].outputs["norm_vcf"]))
        assert len(recs) == n_norm == f.summary["records"]            # soft filter: same record count
        passed = {(r.pos) for r in recs if list(r.filter.keys()) in ([], ["PASS"])}
        assert {HOM_SNV, HET_SNV1, HET_SNV2} <= passed

    # --- caller comparison: planted variants shared, vcfeval ran, truth benchmark perfect recall
    c = st["compare_callers"].summary
    assert c["sets"] == ["bcftools", "freebayes"]
    assert c["pairwise"][0]["shared"] >= 4 and c["pairwise"][0]["gt_concordance_pct"] == 100.0
    assert c["vcfeval"]["pairwise"][0]["genotype"]["tp_baseline"] >= 4
    truth_rows = {r["set"]: r for r in c["vcfeval"]["truth"]}
    for k in ("bcftools", "freebayes"):
        assert truth_rows[k]["overall"]["sensitivity"] == 1.0, truth_rows[k]
        assert truth_rows[k]["overall"]["fn"] == 0
    out = Path(run.output_dir)
    assert (out / "comparison" / "callers" / "comparison_report.html").exists()
    with gzip.open(out / "comparison" / "callers" / "variant_matrix.tsv.gz", "rt") as fh:
        assert fh.readline().startswith("CHROM\tPOS\tREF\tALT\tTYPE\tN_SETS\tGT_bcftools\tGT_freebayes")
    report = (out / "reports" / "run_report.html").read_text()
    assert "Known limitations" in report and "Filtering (soft filters" in report

    # --- cleanup (opt-in): trimmed FASTQs gone and listed; removed reads kept
    assert st["cleanup"].status == S.COMPLETED
    assert not Path(st["fastp"].outputs["trimmed_r1"]).exists()
    assert Path(st["fastp"].outputs["failed_reads"]).exists()
    assert (out / "reports" / "cleanup_manifest.json").exists()

    # --- outputs registered as sample files (origin pipeline, checksummed)
    db.session.expire_all()
    reg = [f for f in sample.files if f.origin == "pipeline"]
    assert sorted(f.role for f in reg) == ["bam", "vcf", "vcf"] and all(f.sha256 for f in reg)
    bam = next(f for f in reg if f.role == "bam")
    vcfs = [f for f in reg if f.role == "vcf"]
    bed = next(f for f in sample.files if f.role == "bed")

    # --- VCF vs VCF (the two registered VCFs + truth) ---------------------------------
    run2 = _reload(create_run(sample, "vcf_compare", vcfs + [bed] + [f for f in sample.files if f.role.startswith("truth")],
                              mini_ref, params={"contigs": "chr22"}))
    assert run2.status == S.COMPLETED, run2.error_message
    s2 = {s.key: s for s in run2.steps}["vcf_compare"].summary
    assert len(s2["sets"]) == 2 and s2["in_all"] >= 4
    assert all(r["overall"]["sensitivity"] == 1.0 for r in s2["vcfeval"]["truth"])

    # --- BAM vs BAM (pipeline BAM vs a 50 % subsample registered as 'external') -------
    sub = tmp_path / "subsample.bam"
    subprocess.run(["samtools", "view", "-b", "-s", "7.5", "-o", sub, bam.path], check=True)
    subprocess.run(["samtools", "index", sub], check=True)
    ext = _reg(sample, "bam", "external", sub, label="provider")
    run3 = _reload(create_run(sample, "bam_compare", [bam, ext, bed], mini_ref,
                              params={"contigs": "chr22", "coverage_thresholds": "1,10,20"}))
    assert run3.status == S.COMPLETED, run3.error_message
    s3 = {s.key: s for s in run3.steps}["bam_compare"].summary
    rows = {r["metric"]: r for r in s3["table"]}
    labels = [b["label"] for b in s3["bams"]]
    prim = rows["Primary reads"]
    assert 0.4 < prim["provider"] / prim[labels[0]] < 0.6               # subsample halves the reads
    assert "% target bases ≥ 20x" in rows and "MAPQ 60+ (% reads)" in rows
    assert Path({s.key: s for s in run3.steps}["bam_compare"].outputs["report_html"]).exists()

    # --- BAM vs VCF evidence ------------------------------------------------------------
    run4 = _reload(create_run(sample, "bam_vcf_evidence", [bam, vcfs[0], bed], mini_ref, params={"contigs": "chr22"}))
    assert run4.status == S.COMPLETED, run4.error_message
    ev = {s.key: s for s in run4.steps}["read_evidence"]
    rows = {}
    with gzip.open(ev.outputs["evidence_tsv"], "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            rows[int(r["pos"])] = r
    assert float(rows[HOM_SNV]["vaf"]) > 0.95
    assert 0.3 < float(rows[HET_SNV1]["vaf"]) < 0.7
    assert ev.summary["evaluated"] >= 4
    assert Path(ev.outputs["igv_batch"]).read_text().startswith("new\ngenome ")


def test_vcf_from_other_build_is_refused(app, make_sample, mini, mini_ref, tmp_path):
    """A VCF whose REF alleles do not match the reference (= other build) must never be compared."""
    sample = _sample(make_sample, mini)
    seq = mini["chr22"]
    paths = []
    for i, shift in enumerate((0, 1)):
        txt = tmp_path / f"v{i}.vcf"
        with open(txt, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr22,length=30000>\n"
                     '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\n")
            for pos in range(100, 29000, 250):
                ref = seq[pos - 1 + shift]          # shift=1: REF taken from the wrong position
                alt = "A" if ref != "A" else "C"
                fh.write(f"chr22\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t0/1\n")
        pysam.tabix_index(str(txt), preset="vcf", force=True)
        paths.append(_reg(sample, "vcf", "external", str(txt) + ".gz"))
    run = _reload(create_run(sample, "vcf_compare", paths, mini_ref, params={"contigs": "chr22"}))
    assert run.status == S.FAILED
    assert "do not match" in (run.error_message or "") or any(
        "do not match" in (s.message or "") for s in run.steps)


def test_comparison_input_checks(app, make_sample, mini, mini_ref):
    sample = _sample(make_sample, mini)
    with pytest.raises(WorkflowError, match="at least 2 VCF"):
        create_run(sample, "vcf_compare", list(sample.files), mini_ref, enqueue=False)
    with pytest.raises(WorkflowError, match="2-6 BAM"):
        create_run(sample, "bam_compare", list(sample.files), mini_ref, enqueue=False)


@pytest.mark.skipif(not image_id(DEFAULT_IMAGE) if shutil.which("docker") else True,
                    reason="DeepVariant Docker image not available")
def test_deepvariant_calls_planted_variants(app, make_sample, mini, mini_ref):
    sample = _sample(make_sample, mini)
    r1 = _reload(create_run(sample, "fastq_to_bam", list(sample.files), mini_ref, params={"contigs": "chr22"}))
    assert r1.status == S.COMPLETED, r1.error_message
    db.session.expire_all()
    bam = next(f for f in sample.files if f.origin == "pipeline" and f.role == "bam")
    bed = next(f for f in sample.files if f.role == "bed")
    r2 = _reload(create_run(sample, "bam_to_vcf", [bam, bed], mini_ref,
                            params={"contigs": "chr22", "callers": ["deepvariant"], "threads": 2}))
    assert r2.status == S.COMPLETED, r2.error_message
    st = {s.key: s for s in r2.steps}
    assert st["compare_callers"].status == S.SKIPPED
    with pysam.VariantFile(st["filter_deepvariant"].outputs["filtered_vcf"]) as fh:
        calls = {r.pos: tuple(sorted(r.samples[0]["GT"])) for r in fh if "PASS" in r.filter.keys()}
    assert calls.get(HOM_SNV) == (1, 1) and calls.get(HET_SNV1) == (0, 1) and calls.get(HET_SNV2) == (0, 1)
    assert st["call_deepvariant"].summary["provenance"]["image_id"].startswith("sha256:")
