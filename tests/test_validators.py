from pathlib import Path

import pytest

from app.utils.validators import (
    ValidationError,
    check_role_matches_extension,
    guess_role,
    resolve_within_roots,
    validate_sample_id,
)


@pytest.mark.parametrize("value", ["SAMPLE_A", "SAMPLE_B", "a.b-c_1"])
def test_valid_sample_ids(value):
    assert validate_sample_id(value) == value


@pytest.mark.parametrize("value", ["", "../etc", "a b", "x;rm -rf /", "-start", "a" * 65, "$(id)"])
def test_invalid_sample_ids(value):
    with pytest.raises(ValidationError):
        validate_sample_id(value)


def test_path_must_be_inside_roots(tmp_path):
    inside = tmp_path / "data" / "x.bam"
    inside.parent.mkdir()
    inside.write_text("x")
    assert resolve_within_roots(str(inside), [tmp_path]) == inside.resolve()
    with pytest.raises(ValidationError):
        resolve_within_roots("/etc/passwd", [tmp_path])
    with pytest.raises(ValidationError):  # traversal
        resolve_within_roots(str(tmp_path / "data" / ".." / ".." / "etc"), [tmp_path])
    with pytest.raises(ValidationError):  # relative
        resolve_within_roots("data/x.bam", [tmp_path])


def test_symlink_escape_rejected(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    (root / "link").symlink_to("/etc")
    with pytest.raises(ValidationError):
        resolve_within_roots(str(root / "link" / "passwd"), [root])


@pytest.mark.parametrize("name,role", [
    ("SAMPLE_A_S1_L001_R1_001.fastq.gz", "fastq_r1"),
    ("SAMPLE_A_S1_L001_R2_001.fastq.gz", "fastq_r2"),
    ("sample_1.fq.gz", "fastq_r1"),
    ("sample.fastq.gz", "fastq_se"),
    ("x.sorted.bam", "bam"),
    ("calls.vcf.gz", "vcf"),
    ("notes.txt", None),
])
def test_guess_role(name, role):
    assert guess_role(Path(name)) == role


def test_role_extension_mismatch():
    with pytest.raises(ValidationError):
        check_role_matches_extension(Path("x.bam"), "vcf")
