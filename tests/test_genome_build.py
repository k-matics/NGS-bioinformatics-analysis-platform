from app.services.genome_build import compare_to_reference, detect_build, naming_style
from app.services.reference_service import inspect_reference
from conftest import GRCH37, GRCH38, write_reference


def test_detect_build():
    assert detect_build(GRCH38) == "GRCh38"
    assert detect_build(GRCH37) == "GRCh37"
    assert detect_build({k.replace("chr", ""): v for k, v in GRCH38.items()}) == "GRCh38"
    assert detect_build({"contigA": 100}) == "unknown"


def test_naming_style():
    assert naming_style(GRCH38) == "ucsc"
    assert naming_style({"1": 1, "2": 2}) == "ensembl"


def test_grch37_file_vs_grch38_reference_is_incompatible():
    result = compare_to_reference(GRCH37, GRCH38)
    assert result["verdict"] == "incompatible"
    assert "GRCh37" in result["reason"]


def test_naming_mismatch_is_incompatible():
    ensembl = {k.replace("chr", ""): v for k, v in GRCH38.items()}
    assert compare_to_reference(ensembl, GRCH38)["verdict"] == "incompatible"


def test_decoys_give_primary_compatible():
    dragen_like = {**GRCH38, "chrUn_JTFH01000001v1_decoy": 25139, "HLA-A*01:01:01:01": 3503}
    result = compare_to_reference(dragen_like, GRCH38)
    assert result["verdict"] == "primary_compatible"
    assert result["missing_in_reference_count"] == 2


def test_identical():
    assert compare_to_reference(GRCH38, GRCH38)["verdict"] == "identical"


def test_inspect_reference_ok_and_alt_warning(tmp_path):
    fasta = write_reference(tmp_path, GRCH38, alt=2)
    rep = inspect_reference(fasta)
    assert rep["build"] == "GRCh38" and rep["status"] == "warning"
    assert rep["files"]["dict"]["matches_fai"]
    assert any("ALT" in w for w in rep["warnings"])


def test_inspect_reference_missing_fai_is_error(tmp_path):
    fasta = write_reference(tmp_path, GRCH38)
    (tmp_path / "ref.fa.fai").unlink()
    rep = inspect_reference(fasta)
    assert rep["status"] == "error"


def test_inspect_reference_bwa_index_from_other_fasta(tmp_path):
    fasta = write_reference(tmp_path, GRCH38)
    (tmp_path / "ref.fa.ann").write_text("1000 999 11\n")
    assert any("BWA index" in e for e in inspect_reference(fasta)["errors"])


def test_missing_fasta(tmp_path):
    assert inspect_reference(tmp_path / "nope.fa")["status"] == "error"
