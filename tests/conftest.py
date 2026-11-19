"""Shared fixtures: an isolated app (temp DB + temp results), a synthetic
reference (metadata files only) and tiny FASTQ/BAM/VCF test files."""
import gzip
import os
import sys
from pathlib import Path

import pysam
import pytest

os.environ["NGS_SKIP_LOCAL_ENV"] = "1"      # never read this server's config/local.env in tests
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app import create_app  # noqa: E402
from app.extensions import db  # noqa: E402
from app.models import ReferenceGenome, Sample, SampleFile  # noqa: E402
from app.services.file_service import file_fingerprint  # noqa: E402

GRCH38 = {"chr1": 248_956_422, "chr2": 242_193_529, "chr22": 50_818_468, "chrX": 156_040_895, "chrM": 16_569}
GRCH37 = {"chr1": 249_250_621, "chr2": 243_199_373, "chr22": 51_304_566, "chrX": 155_270_560, "chrM": 16_569}


def write_reference(directory: Path, contigs: dict, alt: int = 0) -> Path:
    """Write a FASTA stub with .fai/.dict/BWA index stubs describing `contigs`.
    Only metadata is realistic; sequences are tiny (enough for header checks)."""
    directory.mkdir(parents=True, exist_ok=True)
    contigs = dict(contigs)
    for i in range(alt):
        contigs[f"chr1_KI27076{i}v1_alt"] = 1000 + i
    fasta = directory / "ref.fa"
    fasta.write_text("".join(f">{c}\nACGTACGTNN\n" for c in contigs))
    (directory / "ref.fa.fai").write_text("".join(f"{c}\t{l}\t0\t60\t61\n" for c, l in contigs.items()))
    (directory / "ref.dict").write_text("@HD\tVN:1.6\n" + "".join(f"@SQ\tSN:{c}\tLN:{l}\n" for c, l in contigs.items()))
    for suffix in (".amb", ".bwt", ".pac", ".sa"):
        (directory / f"ref.fa{suffix}").write_text("x")
    (directory / "ref.fa.ann").write_text(f"1000 {len(contigs)} 11\n")
    return fasta


def write_fastq_pair(directory: Path, n: int = 200, broken_r2: bool = False) -> tuple[Path, Path]:
    r1, r2 = directory / "S1_R1_001.fastq.gz", directory / "S1_R2_001.fastq.gz"
    with gzip.open(r1, "wt") as a, gzip.open(r2, "wt") as b:
        for i in range(n):
            a.write(f"@read{i} 1:N:0:ACGT\nACGTACGTAC\n+\nIIIIIIIII9\n")
            name = f"read{i + 1 if broken_r2 and i > 50 else i}"
            b.write(f"@{name} 2:N:0:ACGT\nTTGCATGCAA\n+\nIIIII9IIII\n")
    return r1, r2


def write_bam(path: Path, contigs: dict, read_group: bool = True) -> Path:
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": c, "LN": l} for c, l in contigs.items()]}
    if read_group:
        header["RG"] = [{"ID": "rg1", "SM": "S1", "LB": "lib1", "PL": "ILLUMINA"}]
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i in range(5):
            seg = pysam.AlignedSegment(out.header)
            seg.query_name, seg.query_sequence = f"r{i}", "ACGTACGTAC"
            seg.flag, seg.reference_id, seg.reference_start = 0, 0, 1000 + i * 10
            seg.mapping_quality, seg.cigarstring = 60, "10M"
            seg.query_qualities = pysam.qualitystring_to_array("IIIIIIIIII")
            out.write(seg)
    return path


@pytest.fixture
def app(tmp_path):
    ref = write_reference(tmp_path / "ref", GRCH38, alt=3)
    application = create_app("config.TestConfig", overrides={
        "SQLALCHEMY_DATABASE_URI": f"sqlite:///{tmp_path / 'test.db'}",
        "RESULTS_DIR": tmp_path / "results",
        "UPLOAD_DIR": tmp_path / "uploads",
        "LOG_DIR": tmp_path / "logs",
        "DATA_ROOTS": [tmp_path],
        "DEFAULT_REFERENCE_FASTA": ref,
        "DEFAULTS": {"threads_per_job": 2, "min_free_disk_gb": 0, "retain_intermediates": True,
                     "gatk_path": "", "fastq_check_records": 1000},
    })
    with application.app_context():
        yield application
        db.session.remove()


@pytest.fixture
def client(app):
    return app.test_client()


@pytest.fixture
def data_dir(tmp_path):
    d = tmp_path / "data"
    d.mkdir()
    return d


@pytest.fixture
def make_sample(app):
    def _make(sample_id="S1", files=()):
        sample = Sample(sample_id=sample_id, sequencing_type="WES", library_layout="paired")
        db.session.add(sample)
        for role, path, origin in files:
            db.session.add(SampleFile(sample=sample, role=role, origin=origin, path=str(path),
                                      **file_fingerprint(Path(path))))
        db.session.commit()
        return sample
    return _make


@pytest.fixture
def reference(app):
    return ReferenceGenome.query.filter_by(is_default=True).one()
