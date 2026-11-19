"""End-to-end tests of the Phase 3 germline workflows on a simulated mini-genome.

A 30 kb random "chr22" (plus a 5 kb "chr1") is generated, BWA-indexed, and
paired-end reads are simulated from a diploid genome carrying KNOWN variants
(homozygous SNV, heterozygous SNVs, a heterozygous 2-bp deletion). Adapter
read-through and PCR duplicates are included on purpose so that trimming
and duplicate marking have something to do.

The tests check the science, not just that files exist:
  * the planted variants are called with the right genotype by each caller,
  * no reads are lost between FASTQ, trimming, alignment and duplicate marking,
  * duplicates are flagged, not removed,
  * build-incompatible BAMs are refused.
"""
import gzip
import random
import shutil
import subprocess
from pathlib import Path

import pysam
import pytest

from app.engine.params import ParamError
from app.extensions import db
from app.models import ReferenceGenome, WorkflowRun
from app.models import status as S
from app.services.reference_service import apply_report, inspect_reference
from app.services.tool_service import find_gatk_java_home
from app.services.workflow_service import WorkflowError, create_run
from conftest import GRCH37, write_bam

TOOLS = ("bwa", "samtools", "bcftools", "fastp", "fastqc", "multiqc", "mosdepth", "freebayes", "bgzip")
pytestmark = pytest.mark.skipif(any(shutil.which(t) is None for t in TOOLS),
                                reason="bioinformatics tools not on PATH")

ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
COMP = str.maketrans("ACGT", "TGCA")
# 1-based VCF positions of planted variants on chr22
HOM_SNV, HET_SNV1, HET_SNV2, HET_DEL = 5001, 10001, 15001, 20001


def revcomp(s):
    return s.translate(COMP)[::-1]


@pytest.fixture(scope="module")
def mini(tmp_path_factory):
    """Reference + simulated reads, built once for the module."""
    rnd = random.Random(42)
    d = tmp_path_factory.mktemp("mini")
    chr22 = "".join(rnd.choice("ACGT") for _ in range(30_000))
    chr1 = "".join(rnd.choice("ACGT") for _ in range(5_000))
    fasta = d / "mini.fa"
    with open(fasta, "w") as fh:
        for name, seq in (("chr1", chr1), ("chr22", chr22)):
            fh.write(f">{name}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    subprocess.run(["samtools", "faidx", fasta], check=True)
    subprocess.run(["samtools", "dict", "-o", d / "mini.dict", fasta], check=True)
    subprocess.run(["bwa", "index", fasta], check=True, capture_output=True)

    def mutate(seq, pos1, alt):
        i = pos1 - 1
        return seq[:i] + alt + seq[i + 1:]

    ref_base = {p: chr22[p - 1] for p in (HOM_SNV, HET_SNV1, HET_SNV2)}
    alt_base = {p: next(b for b in "ACGT" if b != ref_base[p]) for p in ref_base}
    h1 = mutate(chr22, HOM_SNV, alt_base[HOM_SNV])
    h1 = mutate(h1, HET_SNV1, alt_base[HET_SNV1])
    h2 = mutate(chr22, HOM_SNV, alt_base[HOM_SNV])
    h2 = mutate(h2, HET_SNV2, alt_base[HET_SNV2])
    # heterozygous 2-bp deletion on haplotype 2 (deletes the 2 bases after HET_DEL)
    h2 = h2[:HET_DEL] + h2[HET_DEL + 2:]

    r1p, r2p = d / "SIM_R1_001.fastq.gz", d / "SIM_R2_001.fastq.gz"
    qual = "I" * 100
    n = 0
    with gzip.open(r1p, "wt") as f1, gzip.open(r2p, "wt") as f2:
        def emit(frag):
            nonlocal n
            n += 1
            a = (frag + ADAPTER + "A" * 100)[:100]
            b = (revcomp(frag) + ADAPTER2 + "A" * 100)[:100]
            name = f"@SIM:1:FCSIM01:1:1101:{1000 + n * 7 % 30000}:{1000 + n}"
            f1.write(f"{name} 1:N:0:ACGTACGT\n{a}\n+\n{qual}\n")
            f2.write(f"{name} 2:N:0:ACGTACGT\n{b}\n+\n{qual}\n")
        for i in range(7000):
            hap = h1 if rnd.random() < 0.5 else h2
            size = 70 if i % 20 == 0 else rnd.randint(250, 350)     # 5 % adapter read-through
            start = rnd.randint(0, len(hap) - size)
            frag = hap[start:start + size]
            if rnd.random() < 0.5:
                frag = revcomp(frag)
            emit(frag)
            if i % 10 == 0:          # 10 % PCR duplicates (same fragment, new cluster)
                emit(frag)
    bed = d / "targets.bed"
    bed.write_text("chr22\t4000\t6000\nchr22\t9000\t16000\nchr22\t19000\t21000\n")
    return {"fasta": fasta, "r1": r1p, "r2": r2p, "bed": bed, "pairs": n, "alt": alt_base, "chr22": chr22}


@pytest.fixture
def mini_ref(app, mini):
    ref = ReferenceGenome(name="mini", fasta_path=str(mini["fasta"]))
    apply_report(ref, inspect_reference(mini["fasta"]))
    db.session.add(ref)
    db.session.commit()
    assert ref.status != "error", ref.check_report
    return ref


def _sample(make_sample, mini, with_bed=True):
    files = [("fastq_r1", mini["r1"], "internal"), ("fastq_r2", mini["r2"], "internal")]
    if with_bed:
        files.append(("bed", mini["bed"], "internal"))
    return make_sample(files=files)


def _reload(run):
    db.session.expire_all()
    return db.session.get(WorkflowRun, run.id)


def _genotypes(vcf):
    out = {}
    with pysam.VariantFile(str(vcf)) as fh:
        for rec in fh:
            gt = tuple(sorted(a for a in rec.samples[0]["GT"] if a is not None))
            out[(rec.chrom, rec.pos, rec.ref, rec.alts[0])] = gt
    return out


def _check_planted(gts, mini, caller):
    alt = mini["alt"]
    seq = mini["chr22"]
    expected = {
        ("chr22", HOM_SNV, seq[HOM_SNV - 1], alt[HOM_SNV]): (1, 1),
        ("chr22", HET_SNV1, seq[HET_SNV1 - 1], alt[HET_SNV1]): (0, 1),
        ("chr22", HET_SNV2, seq[HET_SNV2 - 1], alt[HET_SNV2]): (0, 1),
    }
    for key, gt in expected.items():
        assert gts.get(key) == gt, f"{caller}: expected {key} {gt}, got {gts.get(key)}"
    # the 2-bp deletion, left-aligned by bcftools norm (position may shift left in a repeat)
    dels = [k for k, g in gts.items() if len(k[2]) == 3 and len(k[3]) == 1 and abs(k[1] - HET_DEL) <= 5]
    assert dels and gts[dels[0]] == (0, 1), f"{caller}: het deletion near {HET_DEL} not found: {dels}"


# ---------------------------------------------------------------------------
def test_fastq_to_vcf_end_to_end(make_sample, mini, mini_ref):
    sample = _sample(make_sample, mini)
    run = create_run(sample, "fastq_to_vcf", list(sample.files), mini_ref,
                     params={"contigs": "chr22", "callers": ["bcftools", "freebayes"],
                             "caller_params": {"freebayes": {"min_coverage": 3}}})
    run = _reload(run)
    assert run.status == S.COMPLETED, run.error_message
    st = {s.key: s for s in run.steps}
    assert st["remove_duplicates"].status == S.SKIPPED

    # --- read accounting: nothing lost between steps
    fp, al, md = st["fastp"].summary, st["bwa_align"].summary, st["markdup"].summary
    assert fp["reads_before"] == 2 * mini["pairs"]
    assert fp["reads_after"] + fp["reads_removed"] == fp["reads_before"]
    assert fp["adapter_trimmed_reads"] > 0
    assert al["primary"] == fp["reads_after"] == md["primary"]
    assert al["read_group"]["ID"] == "FCSIM01.1" and al["read_group"]["SM"] == "S1"
    # duplicates flagged (about 10 % planted), never removed
    assert 5 < md["duplicate_pct"] < 20
    assert Path(st["fastp"].outputs["failed_reads"]).exists()

    # --- coverage only over the BED targets
    qc = st["bam_qc"].summary
    assert qc["coverage_region"] == "targets" and qc["target_bases"] == 11_000
    assert qc["mean_target_coverage"] > 20

    # --- variant calls: planted genotypes found by both callers
    for caller in ("bcftools", "freebayes"):
        s = st[f"call_{caller}"]
        assert s.status == S.COMPLETED
        _check_planted(_genotypes(s.outputs["norm_vcf"]), mini, caller)
        assert Path(s.outputs["norm_vcf"] + ".tbi").exists()
    # calls stay inside padded targets
    for key in _genotypes(st["call_bcftools"].outputs["norm_vcf"]):
        assert any(a - 100 < key[1] <= b + 100 for a, b in ((4000, 6000), (9000, 16000), (19000, 21000)))

    out = Path(run.output_dir)
    assert (out / "reports" / "run_report.html").read_text().count("Research use only") >= 1
    assert (out / "reports" / "multiqc" / "multiqc_report.html").exists()
    assert (out / "run_manifest.json").exists()
    # the BAM in the run folder is coordinate-sorted, indexed, carries the read group
    bam = st["index_bam"].outputs["bam"]
    with pysam.AlignmentFile(bam) as fh:
        assert fh.header.to_dict()["HD"]["SO"] == "coordinate"
        assert fh.header.to_dict()["RG"][0]["SM"] == "S1"


@pytest.mark.skipif(find_gatk_java_home() is None, reason="no Java 17 for GATK")
def test_bam_to_vcf_gatk_on_pipeline_bam(app, make_sample, mini, mini_ref, tmp_path):
    from app.services.settings_service import set_setting
    import os
    gatk = os.environ.get("NGS_GATK_PATH") or _local_env("NGS_GATK_PATH") or shutil.which("gatk") or ""
    if not gatk or not Path(gatk).exists():
        pytest.skip("GATK not installed")
    set_setting("gatk_path", gatk)
    # produce a BAM with the FASTQ -> BAM workflow, then register it as a BAM input
    sample = _sample(make_sample, mini)
    run = _reload(create_run(sample, "fastq_to_bam", list(sample.files), mini_ref,
                             params={"contigs": "chr22", "markdup_tool": "gatk"}))
    assert run.status == S.COMPLETED, run.error_message
    bam = Path({s.key: s for s in run.steps}["index_bam"].outputs["bam"])
    s2 = make_sample("S2", files=[("bam", bam, "internal"), ("bed", mini["bed"], "internal")])
    run2 = _reload(create_run(s2, "bam_to_vcf", list(s2.files), mini_ref,
                              params={"contigs": "chr22", "callers": ["gatk_hc"]}))
    assert run2.status == S.COMPLETED, run2.error_message
    st = {s.key: s for s in run2.steps}
    assert st["bam_prepare"].summary["duplicates_marked"] is True
    assert Path(st["bam_prepare"].outputs["bam"]).is_symlink()      # input linked, not copied
    _check_planted(_genotypes(st["call_gatk_hc"].outputs["norm_vcf"]), mini, "gatk_hc")


def test_bam_to_vcf_refuses_other_build(make_sample, mini_ref, data_dir):
    bam = write_bam(data_dir / "grch37.bam", GRCH37)
    pysam.index(str(bam))
    sample = make_sample(files=[("bam", bam, "external")])
    run = _reload(create_run(sample, "bam_to_vcf", list(sample.files), mini_ref,
                             params={"contigs": "chr22"}))
    assert run.status == S.FAILED
    assert "mix" in run.error_message and "genome builds" in run.error_message


def test_bam_without_read_group_gets_rg_copy_for_gatk(app, make_sample, mini_ref, data_dir):
    """GATK needs @RG: a copy with RG is written into the run folder; the input is untouched."""
    from app.pipelines.bam_processing.steps import BamPrepareStep
    from app.engine.step import StepContext
    from app.engine.runner import CommandRunner
    contigs = {"chr1": 5000, "chr22": 30000}
    bam = write_bam(data_dir / "norg.bam", contigs, read_group=False)
    before = bam.stat().st_mtime
    sample = make_sample(files=[("bam", bam, "external")])
    run = create_run(sample, "bam_to_vcf", list(sample.files), mini_ref,
                     params={"contigs": "chr22", "callers": ["bcftools"]}, enqueue=False)
    run.params = dict(run.params, callers=["gatk_hc"])
    ctx = StepContext(run=run, sample=sample, reference=mini_ref, params=run.params,
                      output_dir=Path(run.output_dir), threads=1, input_files=list(sample.files))
    log = Path(run.output_dir) / "t.log"
    res = BamPrepareStep().execute(ctx, CommandRunner(log))
    assert res.summary["read_group_added"] is True
    assert pysam.AlignmentFile(res.outputs["bam"]).header.to_dict()["RG"][0]["SM"] == "S1"
    assert bam.stat().st_mtime == before
    assert not Path(str(bam) + ".bai").exists()         # nothing written next to the input


def test_read_accounting_mismatch_fails(make_sample, mini, mini_ref, monkeypatch):
    from app.pipelines.alignment import steps as align_steps
    monkeypatch.setattr(align_steps, "expected_input_reads", lambda ctx: 999_999)
    sample = _sample(make_sample, mini)
    run = _reload(create_run(sample, "fastq_to_bam", list(sample.files), mini_ref,
                             params={"contigs": "chr22", "trimming": False}))
    assert run.status == S.FAILED
    assert "Read accounting mismatch" in run.error_message


# --------------------------------------------------------------------------- parameter validation
def test_unknown_contig_rejected(make_sample, mini, mini_ref):
    sample = _sample(make_sample, mini)
    with pytest.raises(WorkflowError, match="not in reference"):
        create_run(sample, "fastq_to_vcf", list(sample.files), mini_ref, params={"contigs": "22"},
                   enqueue=False)


def test_unavailable_caller_rejected(make_sample, mini, mini_ref):
    from app.services.settings_service import set_setting
    set_setting("deepvariant_image", "google/deepvariant:no-such-tag")   # image not present locally
    sample = _sample(make_sample, mini)
    with pytest.raises(WorkflowError, match="DeepVariant"):
        create_run(sample, "fastq_to_vcf", list(sample.files), mini_ref,
                   params={"callers": ["deepvariant"]}, enqueue=False)


def test_caller_param_out_of_range_rejected(make_sample, mini, mini_ref):
    sample = _sample(make_sample, mini)
    with pytest.raises(WorkflowError, match="Min mapping quality"):
        create_run(sample, "bam_to_vcf" if False else "fastq_to_vcf", list(sample.files), mini_ref,
                   params={"callers": ["bcftools"], "caller_params": {"bcftools": {"min_MQ": 500}}},
                   enqueue=False)


def test_caller_params_are_caller_specific(make_sample, mini, mini_ref):
    sample = _sample(make_sample, mini)
    run = create_run(sample, "fastq_to_vcf", list(sample.files), mini_ref,
                     params={"callers": ["freebayes"]}, enqueue=False)
    cp = run.params["caller_params"]
    assert set(cp) == {"freebayes"}                    # bcftools defaults dropped: not selected
    assert cp["freebayes"]["min_alt_fraction"] == 0.2


def test_fastq_workflow_ignores_unrelated_files(make_sample, mini, mini_ref, data_dir):
    bam = write_bam(data_dir / "x.bam", {"chr22": 30000})
    sample = make_sample(files=[("fastq_r1", mini["r1"], "internal"), ("fastq_r2", mini["r2"], "internal"),
                                ("bam", bam, "internal")])
    run = create_run(sample, "fastq_qc", list(sample.files), mini_ref, enqueue=False)
    roles = {f.role for f in sample.files if f.id in run.input_file_ids}
    assert roles == {"fastq_r1", "fastq_r2"}


def test_param_spec_rejects_shell_metacharacters():
    from app.pipelines.germline import REGION_PARAMS
    contigs = next(s for s in REGION_PARAMS if s.key == "contigs")
    with pytest.raises(ParamError):
        contigs.parse("chr22; rm -rf /", True)


def test_start_form_posts_params(client, make_sample, mini, mini_ref):
    """The web form sends p.<workflow>.<param> fields; values reach the run."""
    sample = _sample(make_sample, mini, with_bed=False)
    resp = client.post(f"/samples/{sample.sample_id}/run", data={
        "workflow_type": "fastq_qc", "reference_id": mini_ref.id,
        "file_ids": [f.id for f in sample.files],
        "p.fastq_qc.trimming": "on", "p.fastq_qc.fastp_adapter_trimming": "on",
        "p.fastq_qc.fastp_min_length": "50",
    })
    assert resp.status_code == 302
    run = WorkflowRun.query.filter_by(workflow_type="fastq_qc").one()
    assert run.params["fastp_min_length"] == 50
    assert run.status == S.COMPLETED, run.error_message
    assert {s.key: s for s in run.steps}["fastp"].summary["reads_after"] > 0


def test_bam_with_contigs_missing_from_reference(make_sample, mini, mini_ref, data_dir):
    """External BAMs (e.g. DRAGEN) can carry decoy contigs our reference lacks. QC and calling must
    work on the shared contigs instead of failing on reference look-ups (seen on real data)."""
    contigs = {"chr1": 5000, "chr22": 30000, "chrUn_KN707606v1_decoy": 2000}
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": c, "LN": l} for c, l in contigs.items()],
              "RG": [{"ID": "rg1", "SM": "S1", "LB": "lib1", "PL": "ILLUMINA"}],
              "PG": [{"ID": "markdup", "PN": "samtools"}]}
    bam = data_dir / "dragen_like.bam"
    seq = mini["chr22"]
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for tid, starts in ((1, range(4000, 5900, 20)), (2, range(100, 1000, 50))):
            for i, start in enumerate(starts):
                seg = pysam.AlignedSegment(out.header)
                seg.query_name = f"r{tid}_{i}"
                seg.query_sequence = seq[start:start + 100] if tid == 1 else "A" * 100
                seg.flag, seg.reference_id, seg.reference_start = 0, tid, start
                seg.mapping_quality, seg.cigarstring = 60, "100M"
                seg.query_qualities = pysam.qualitystring_to_array("I" * 100)
                out.write(seg)
    pysam.index(str(bam))
    sample = make_sample(files=[("bam", bam, "external"), ("bed", mini["bed"], "internal")])
    run = _reload(create_run(sample, "bam_to_vcf", list(sample.files), mini_ref,
                             params={"contigs": "chr22", "callers": ["bcftools"]}))
    assert run.status == S.COMPLETED, run.error_message
    st = {s.key: s for s in run.steps}
    assert st["bam_prepare"].summary["verdict"] == "primary_compatible"
    assert st["bam_qc"].summary["stats_scope"] == "chr22"


def test_resume_refused_after_workflow_version_change(make_sample, mini, mini_ref):
    from app.services.workflow_service import resume_run
    sample = _sample(make_sample, mini)
    run = create_run(sample, "fastq_qc", list(sample.files), mini_ref, enqueue=False)
    run.status, run.workflow_version = S.FAILED, "0.9"
    db.session.commit()
    with pytest.raises(WorkflowError, match="workflow changed"):
        resume_run(run)


def _local_env(key):
    """A value from this server's config/local.env (tests do not load that file globally)."""
    f = Path(__file__).resolve().parents[1] / "config" / "local.env"
    if f.is_file():
        for line in f.read_text().splitlines():
            if line.startswith(key + "="):
                return line.split("=", 1)[1].strip()
    return None
