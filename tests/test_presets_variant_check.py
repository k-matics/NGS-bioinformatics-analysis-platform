"""Parameter presets, filter-profile import and the variant-by-variant check workflow."""
import gzip
import shutil
from pathlib import Path

import pysam
import pytest

from app.extensions import db
from app.models import ParamPreset, SampleFile, WorkflowRun
from app.models import status as S
from app.pipelines.comparison.variant_check import what_if
from app.services import preset_service as ps
from app.services.filter_profile import extract_profile, profile_to_params
from app.services.workflow_service import WorkflowError, create_run
from test_germline import TOOLS, _sample, mini, mini_ref  # noqa: F401

needs_tools = pytest.mark.skipif(any(shutil.which(t) is None for t in TOOLS + ("rtg",)),
                                 reason="bioinformatics tools not on PATH")


def _vcf(path: Path, header_extra: str, rows: list[str], sample="S") -> Path:
    txt = path.with_suffix("")
    with open(txt, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=5000>\n##contig=<ID=chr22,length=30000>\n"
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                 '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
                 '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
                 '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n' + header_extra +
                 f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        fh.write("".join(r + "\n" for r in rows))
    pysam.tabix_index(str(txt), preset="vcf", force=True)
    return Path(str(txt) + ".gz")


DRAGEN_HDR = ('##FILTER=<ID=DRAGENSnpHardQUAL,Description="Set if true:QUAL < 0.4">\n'
              '##FILTER=<ID=LowDepth,Description="Set if true:DP <= 1">\n'
              '##FILTER=<ID=MosaicLowAF,Description="Set if true:AF < 0.2">\n'
              '##FILTER=<ID=TargetedLowGQ,Description="Set if call has GQ < 3 and JGQ is not present.">\n'
              '##FILTER=<ID=weak_evidence,Description="Somatic variant score does not meet the threshold of 3">\n'
              '##DRAGENCommandLine=<ID=dragen,Version="SW: 4.3.6",Date="x",CommandLineOptions="-r /ref '
              '--enable-variant-caller true --cnv-target-bed /t.bed">\n##source=DRAGEN_SNV\n')


# ----------------------------------------------------------------------------- presets
def test_builtin_preset_is_valid_and_locked(app):
    ps.ensure_builtin_presets()
    p = ParamPreset.query.filter_by(name=ps.DEFAULT_PRESET).one()
    assert p.locked and p.origin == "builtin"
    assert set(p.params["callers"]) >= {"bcftools", "gatk_hc"}
    assert p.params["caller_params"]["gatk_hc"]["filter_min_depth"] == 10
    clean, problems = ps.validate_preset_params(p.params)
    assert not problems and clean == p.params


def test_preset_validation_rejects_bad_values(app):
    _, problems = ps.validate_preset_params({"nope": 1, "fastp_min_length": 5, "trimming": "yes",
                                             "caller_params": {"gatk_hc": {"filter_min_depth": -1, "x": 1},
                                                               "ghost": {}}})
    text = " | ".join(problems)
    assert "Unknown parameter 'nope'" in text and "unknown caller 'ghost'" in text
    assert "≥ 15" in text and "true/false" in text and "unknown parameter 'x'" in text
    assert "≥ 0" in text
    clean, problems = ps.validate_preset_params({"caller_params": {"gatk_hc": {"filter_min_depth": 2}}})
    assert not problems and clean == {"caller_params": {"gatk_hc": {"filter_min_depth": 2}}}


def test_preset_diff_lists_changed_values():
    preset = {"trimming": True, "caller_params": {"gatk_hc": {"filter_min_depth": 2}, "freebayes": {"x": 1}}}
    final = {"trimming": False, "callers": ["gatk_hc"], "caller_params": {"gatk_hc": {"filter_min_depth": 2}}}
    assert ps.preset_diff(preset, final) == ["trimming"]   # freebayes not run -> not a change


# ----------------------------------------------------------------------------- filter profile
def test_dragen_profile_is_mapped_explicitly(app, tmp_path):
    rows = ["chr22\t100\t.\tA\tG\t30\tPASS\t.\tGT:AD:DP:GQ\t0/1:5,6:11:40",
            "chr22\t200\t.\tC\tT\t5\tLowDepth\t.\tGT:AD:DP:GQ\t1/1:0,1:1:3",
            "chr22\t300\t.\tG\tA\t2\tPASS\t.\tGT:AD:DP:GQ\t0/1:3,1:4:2",
            "chr22\t400\t.\tT\tC\t9\tMosaicLowAF\t.\tGT:AD:DP:GQ\t0/1:30,3:33:9"]
    prof = extract_profile(_vcf(tmp_path / "d.vcf.gz", DRAGEN_HDR, rows))
    assert prof["caller"] == "DRAGEN" and prof["records"] == 4
    assert prof["filter_counts"] == {"PASS": 2, "LowDepth": 1, "MosaicLowAF": 1}
    assert prof["observed_pass"]["DP"]["min"] == 4 and prof["observed_pass"]["VAF"]["min"] == 0.25
    assert any(f["key"] == "--vc-target-bed" and f["value"] == "not set" for f in prof["facts"])
    lowdepth = next(f for f in prof["filters"] if f["id"] == "LowDepth")
    assert lowdepth["used"] and lowdepth["rules"] == [{"field": "DP", "op": "<=", "value": 1.0}]
    assert not next(f for f in prof["filters"] if f["id"] == "weak_evidence")["used"]

    params, notes = profile_to_params(prof)
    cp = params["caller_params"]
    assert all(cp[c]["filter_min_depth"] == 2 for c in cp)            # DP <= 1 fails -> keep DP >= 2
    assert cp["gatk_hc"]["filter_min_gq"] == 0                       # TargetedLowGQ is not germline
    assert cp["gatk_hc"]["filter_min_vaf"] == 0.0                    # MosaicLowAF is mosaic-only
    assert cp["bcftools"]["filter_min_qual"] == 0.4 and cp["gatk_hc"]["stand_call_conf"] == 1.0
    assert cp["gatk_hc"]["filter_gatk_hard"] is False and cp["freebayes"]["filter_strand"] is False
    status = {n["param"]: n["status"] for n in notes}
    assert status["filter_min_depth"] == "mapped" and status["gatk_hc.stand_call_conf"] == "approximate"
    assert status["MosaicLowAF"] == "not mapped"
    _, problems = ps.validate_preset_params(params)
    assert not problems


def test_preset_routes_import_and_start_form(app, tmp_path, make_sample):
    client = app.test_client()
    vcf = _vcf(tmp_path / "prov.vcf.gz", DRAGEN_HDR, ["chr22\t100\t.\tA\tG\t30\tPASS\t.\tGT:AD:DP:GQ\t0/1:5,6:11:40"])
    r = client.post("/presets/import", data={"path": str(vcf), "name": "Provider-like"})
    assert r.status_code == 302
    p = ParamPreset.query.filter_by(name="Provider-like").one()
    assert p.origin == "imported" and p.profile["caller"] == "DRAGEN" and p.notes
    assert client.get(f"/presets/{p.id}").status_code == 200
    assert b"Provider-like" in client.get("/presets/").data
    # outside the allowed folders -> refused
    r = client.post("/presets/import", data={"path": "/etc/passwd", "name": "x"}, follow_redirects=True)
    assert ParamPreset.query.filter_by(name="x").first() is None
    # edit: invalid JSON value refused, valid saved
    r = client.post(f"/presets/{p.id}/edit", data={"name": "Provider-like", "params_json":
                                                    '{"caller_params": {"gatk_hc": {"filter_min_depth": -5}}}'},
                    follow_redirects=True)
    assert b"Not saved" in r.data
    builtin = ParamPreset.query.filter_by(locked=True).first()
    client.post(f"/presets/{builtin.id}/delete")
    assert db.session.get(ParamPreset, builtin.id) is not None       # locked

    sample = make_sample(sample_id="S2")
    page = client.get(f"/samples/{sample.sample_id}")
    assert page.status_code == 200 and b"preset-select" in page.data and b"Provider-like" in page.data


# ----------------------------------------------------------------------------- variant check
def test_what_if_thresholds():
    rows = [{"STATUS": "found_pass", "Q_DP": "30", "Q_GQ": "50", "Q_VAF": "0.5", "Q_QUAL": "100"},
            {"STATUS": "found_filtered", "Q_DP": "4", "Q_GQ": "10", "Q_VAF": "0.4", "Q_QUAL": "20"},
            {"STATUS": "not_called", "Q_DP": "", "Q_GQ": "", "Q_VAF": "", "Q_QUAL": ""}]
    assert what_if(rows)["found"] == 1
    assert what_if(rows, ignore_our_filters=True)["found"] == 2
    assert what_if(rows, ignore_our_filters=True, min_dp=5)["found"] == 1
    assert what_if(rows, ignore_our_filters=True, min_dp=2)["found_pct"] == 66.67


@needs_tools
def test_variant_check_workflow(app, make_sample, mini, mini_ref, tmp_path):
    seq = mini["chr22"]
    def alt(p):
        return "A" if seq[p - 1] != "A" else "C"
    P = [1000, 2000, 3000, 4000, 5000, 6000]
    base_rows = [f"chr22\t{p}\t.\t{seq[p - 1]}\t{alt(p)}\t50\tPASS\t.\tGT:AD:DP:GQ\t0/1:10,10:20:60" for p in P]
    base_rows.append(f"chr22\t7000\t.\t{seq[6999]}\t{alt(7000)}\t3\tLowDepth\t.\tGT:AD:DP:GQ\t1/1:0,1:1:3")
    ours = [
        f"chr22\t1000\t.\t{seq[999]}\t{alt(1000)}\t60\tPASS\t.\tGT:AD:DP:GQ\t0/1:12,11:23:70",       # found_pass
        f"chr22\t2000\t.\t{seq[1999]}\t{alt(2000)}\t60\tPASS\t.\tGT:AD:DP:GQ\t1/1:0,20:20:60",       # GT differs
        f"chr22\t3000\t.\t{seq[2999]}\t{alt(3000)}\t15\tLowDP\t.\tGT:AD:DP:GQ\t0/1:3,3:6:30",        # filtered
        f"chr22\t4000\t.\t{seq[3999]}\t{alt(4000)}\t0\tRefCall\t.\tGT:AD:DP:GQ\t0/0:20,2:22:20",     # ref GT
        f"chr22\t5003\t.\t{seq[5002]}\t{alt(5003)}\t40\tPASS\t.\tGT:AD:DP:GQ\t0/1:10,9:19:50",       # nearby
        f"chr22\t9000\t.\t{seq[8999]}\t{alt(9000)}\t40\tPASS\t.\tGT:AD:DP:GQ\t0/1:10,9:19:50",       # ours only
    ]                                                                                               # 6000: missing
    hdr_f = ('##FILTER=<ID=LowDepth,Description="DP <= 1">\n##FILTER=<ID=LowDP,Description="d">\n'
             '##FILTER=<ID=RefCall,Description="r">\n')
    bvcf = _vcf(tmp_path / "provider.vcf.gz", hdr_f, base_rows, sample="PROV")
    qvcf = _vcf(tmp_path / "ours.vcf.gz", hdr_f, ours, sample="OURS")
    sample = _sample(make_sample, mini)
    from app.services.file_service import file_fingerprint
    fb = SampleFile(sample=sample, role="vcf", origin="external", path=str(bvcf), **file_fingerprint(bvcf))
    fq = SampleFile(sample=sample, role="vcf", origin="internal", path=str(qvcf), **file_fingerprint(qvcf))
    db.session.add_all([fb, fq])
    db.session.commit()

    with pytest.raises(WorkflowError, match="exactly two VCF"):
        create_run(sample, "vcf_variant_check", [fb], mini_ref, enqueue=False)

    run = create_run(sample, "vcf_variant_check", [fb, fq], mini_ref,
                     params={"contigs": "chr22", "baseline_file_id": fb.id})
    db.session.expire_all()
    run = db.session.get(WorkflowRun, run.id)
    assert run.status == S.COMPLETED, run.error_message
    step = {s.key: s for s in run.steps}["variant_check"]
    s = step.summary
    assert s["baseline"]["origin"] == "external" and s["baseline_variants"] == 6
    assert s["skipped"] == {"baseline not PASS": 1}
    assert s["by_status"] == {"found_pass": 2, "found_filtered": 1, "found_ref_gt": 1, "nearby": 1, "not_called": 1}
    assert s["gt_different"] == 1 and s["filtered_by_label"] == {"LowDP": 1} and s["query_only_pass"] == 2  # 5003 + 9000
    with gzip.open(step.outputs["variant_tsv"], "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        rows = {int(r["POS"]): r for r in (dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh)}
    assert rows[3000]["Q_FILTER"] == "LowDP" and rows[3000]["Q_DP"] == "6"
    assert rows[5000]["NEARBY"].startswith("5003:")
    assert rows[1000]["B_VAF"] == "0.5" and rows[1000]["Q_GQ"] == "70"

    client = app.test_client()
    page = client.get(f"/workflows/{run.run_id}/variant-check?ignore_filters=1&min_dp=5")
    assert page.status_code == 200 and b"found_filtered" in page.data
    csv = client.get(f"/workflows/{run.run_id}/variant-check.csv?status=not_called").data.decode()
    assert csv.count("\nchr22,") == 1 and "6000" in csv
    assert Path(step.outputs["report_html"]).read_text().count("Not in our VCF") >= 1


def test_identifiable_output_names():
    from types import SimpleNamespace

    from app.pipelines.common import output_stem, preset_tag, vcf_stem
    assert preset_tag({}) == ""
    assert preset_tag({"preset": {"name": "Provider-like (DRAGEN 4.3)", "changed_after_applying": []}}) \
        == "provider-like-dragen-4-3"
    assert preset_tag({"preset": {"name": "In-house default", "changed_after_applying": ["gatk_hc.x"]}}) \
        == "in-house-default-modified"
    ctx = SimpleNamespace(sample=SimpleNamespace(sample_id="SAMPLE_A"),
                          params={"preset": {"name": "../evil name"}})
    assert output_stem(ctx) == "SAMPLE_A.evil-name"          # never a path separator
    ctx.params = {}
    assert output_stem(ctx) == "SAMPLE_A"
    ctx.input_files = [SimpleNamespace(role="bam", origin="external", label="DRAGEN 4.3 provider", path="/x/a.bam")]
    ctx.params = {"preset": {"name": "Raw tool defaults (no filters)"}}
    assert output_stem(ctx) == "SAMPLE_A.dragen-4-3-provider-bam.raw-tool-defaults-no-filters"
    ctx.input_files = [SimpleNamespace(role="bam", origin="pipeline", label="x", path="/x/a.bam")]
    assert output_stem(ctx) == "SAMPLE_A.raw-tool-defaults-no-filters"
    assert vcf_stem("/x/SAMPLE_A_provider.hard-filtered.vcf.gz") == "SAMPLE_A_provider.hard-filtered"


def test_variant_search_parsing_and_matching():
    import pytest

    from app.pipelines.comparison.variant_check import match_locus, parse_locus
    assert parse_locus("chr1:100000 A>G") == {"chrom": "chr1", "start": 100000, "end": 100000, "ref": "A", "alt": "G"}
    assert parse_locus("chr1:100000A>G")["alt"] == "G"
    assert parse_locus("1-100000-A-G") == {"chrom": "1", "start": 100000, "end": 100000, "ref": "A", "alt": "G"}
    assert parse_locus("chr1 100000 a g")["ref"] == "A"
    assert parse_locus("chr1:99,000-101,000") == {"chrom": "chr1", "start": 99000, "end": 101000, "ref": None, "alt": None}
    assert parse_locus("chrX:5")["ref"] is None
    assert parse_locus("chr1:100 AT>A")["ref"] == "AT"
    for bad in ("", "chr1", "chr1:abc", "chr1:200-100", "chr1:5 T>C; rm -rf /"):
        with pytest.raises(ValueError):
            parse_locus(bad)
    rows = [{"CHROM": "chr1", "POS": "100000", "REF": "A", "ALT": "G"},
            {"CHROM": "chr1", "POS": "100000", "REF": "A", "ALT": "T"},
            {"CHROM": "chr1", "POS": "100007", "REF": "TA", "ALT": "T"},
            {"CHROM": "chr2", "POS": "100000", "REF": "A", "ALT": "G"}]
    assert len(match_locus(rows, parse_locus("chr1:100000 A>G"))) == 1
    assert len(match_locus(rows, parse_locus("1:100000"))) == 2            # 'chr' optional
    assert len(match_locus(rows, parse_locus("chr1:100000"), window=10)) == 3
    assert match_locus(rows, parse_locus("chr1:100001 T>C")) == []


def test_raw_tool_defaults_commands(app):
    from pathlib import Path

    from app.pipelines.variant_calling.callers import get_caller
    from app.pipelines.variant_calling.callers.base import CallJob
    from app.services import preset_service as ps
    with app.app_context():
        ps.ensure_builtin_presets()
        from app.models import ParamPreset
        raw = ParamPreset.query.filter_by(name=ps.RAW_PRESET).one()
        assert raw.locked
        assert raw.params["caller_params"]["gatk_hc"]["tool_defaults"] is True
        assert raw.params["caller_params"]["gatk_hc"]["filter_min_depth"] == 0
        assert raw.params["caller_params"]["gatk_hc"]["filter_gatk_hard"] is False
        for key in ("gatk_hc", "bcftools", "freebayes", "deepvariant"):
            c = get_caller(key)
            prm = {s.key: s.default for s in c.params}
            job = CallJob(bam=Path("/x/s.bam"), reference=Path("/x/r.fa"), regions_bed=Path("/x/r.bed"),
                          out_vcf=Path("/x/run/vcf/k/o.vcf.gz"), work_dir=Path("/x/tmp"), threads=4,
                          sample_id="S", params=prm)
            normal = " ".join(map(str, _flat(c.build_commands(job)[0].commands)))
            job.params = {**prm, "tool_defaults": True}
            rawcmd = " ".join(map(str, _flat(c.build_commands(job)[0].commands)))
            for opt in ("--standard-min-confidence-threshold-for-calling", "--minimum-mapping-quality",
                        "--pcr-indel-model", "min_mapping_quality=", " -q ", " -Q ", " -d ", " -F ", " -C ",
                        "--min-coverage"):
                assert opt not in rawcmd, (key, opt)
            assert "/x/s.bam" in rawcmd and "/x/r.fa" in rawcmd and "/x/r.bed" in rawcmd
            assert len(rawcmd) < len(normal)


def _flat(cmds):
    if cmds and isinstance(cmds[0], list):
        return [x for c in cmds for x in c]
    return cmds


def test_site_summary_reasons():
    from app.services.site_inspect import summarise_site
    reads = [{"mapq": 0, "allele": "ALT", "bq": 30, "dup": False, "qcfail": False, "reverse": i % 2 == 0,
              "xa": ["chr11_KI270721v1_random"], "nm": 1, "as": 148, "xs": 148} for i in range(40)]
    s = summarise_site(reads, 20, 10, False, "/x.bam", "chr11", 2000000, "T", "A")
    assert s["alleles_all"] == {"ALT": 40} and s["alleles_good"] == {}
    assert any("MAPQ 0" in r and "chr11_KI270721v1_random" in r for r in s["reasons"])
    good = [dict(r, mapq=60, xa=[]) for r in reads]
    s = summarise_site(good, 20, 10, False, "/x.bam", "chr11", 2000000, "T", "A")
    assert s["alleles_good"] == {"ALT": 40} and "should be able" in s["reasons"][0]


def test_gatk_dragen_mode_commands(app, tmp_path):
    from pathlib import Path

    from app.pipelines.variant_calling.callers import get_caller
    from app.pipelines.variant_calling.callers.base import CallJob
    from app.pipelines.variant_calling.callers.gatk import dragstr_table_path
    c = get_caller("gatk_hc")
    prm = {s.key: s.default for s in c.params}
    prm.update(dragen_mode=True, min_MQ=0, stand_call_conf=1.0)
    with app.app_context():
        app.config  # instance_path used for the cache
        job = CallJob(bam=Path("/x/s.bam"), reference=tmp_path / "ref.fa", regions_bed=Path("/x/r.bed"),
                      out_vcf=tmp_path / "o.vcf.gz", work_dir=tmp_path, threads=4, sample_id="S", params=prm)
        table = dragstr_table_path(job.reference)
        assert str(tmp_path) not in str(table.parent)           # never written next to the reference
        cmds = [list(map(str, x.commands)) for x in c.build_commands(job)]
        tools = [next(t for t in cmd if t in ("ComposeSTRTableFile", "CalibrateDragstrModel", "HaplotypeCaller", "mv"))
                 for cmd in cmds]
        assert tools == ["ComposeSTRTableFile", "mv", "CalibrateDragstrModel", "HaplotypeCaller"]
        hc = cmds[-1]
        assert hc.index("--dragen-mode") < hc.index("--minimum-mapping-quality")   # explicit values win
        assert hc[hc.index("--minimum-mapping-quality") + 1] == "0"
        assert "--dragstr-params-path" in hc and "--pcr-indel-model" not in hc
        table.parent.mkdir(parents=True, exist_ok=True)
        table.write_text("x")                                     # cached: no rebuild
        try:
            tools2 = [next(t for t in map(str, x.commands) if t in ("ComposeSTRTableFile", "CalibrateDragstrModel",
                                                                      "HaplotypeCaller")) for x in c.build_commands(job)]
            assert tools2 == ["CalibrateDragstrModel", "HaplotypeCaller"]
            job.params = {**prm, "dragstr_calibrate": False}
            hc = list(map(str, c.build_commands(job)[-1].commands))
            assert "--dragen-mode" in hc and "--dragstr-params-path" not in hc and "--pcr-indel-model" in hc
        finally:
            table.unlink()
