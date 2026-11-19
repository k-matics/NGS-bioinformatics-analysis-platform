from ..extensions import db
from .base import utcnow

FILE_ROLES = {
    "fastq_r1": "FASTQ R1",
    "fastq_r2": "FASTQ R2",
    "fastq_se": "FASTQ (single-end)",
    "bam": "BAM",
    "cram": "CRAM",
    "vcf": "VCF",
    "bed": "Target regions BED",
    "truth_vcf": "Truth VCF (benchmark, e.g. GIAB)",
    "truth_bed": "Truth confident-regions BED",
}
FILE_ORIGINS = {"internal": "In-house", "external": "External provider", "pipeline": "Platform pipeline output"}


class Sample(db.Model):
    __tablename__ = "samples"

    id = db.Column(db.Integer, primary_key=True)
    sample_id = db.Column(db.String(64), unique=True, nullable=False, index=True)
    sequencing_type = db.Column(db.String(16), nullable=False, default="WES")
    library_layout = db.Column(db.String(16), nullable=False, default="paired")
    description = db.Column(db.Text)
    created_at = db.Column(db.DateTime, default=utcnow, nullable=False)

    files = db.relationship(
        "SampleFile", back_populates="sample", cascade="all, delete-orphan", order_by="SampleFile.id"
    )
    runs = db.relationship("WorkflowRun", back_populates="sample", order_by="WorkflowRun.id.desc()")

    def files_by_role(self, *roles):
        return [f for f in self.files if f.role in roles]


class SampleFile(db.Model):
    """An input file registered for a sample. Files are referenced in place, never copied or moved."""

    __tablename__ = "sample_files"
    __table_args__ = (db.UniqueConstraint("sample_pk", "path", name="uq_sample_path"),)

    id = db.Column(db.Integer, primary_key=True)
    sample_pk = db.Column(db.Integer, db.ForeignKey("samples.id"), nullable=False)
    role = db.Column(db.String(16), nullable=False)
    origin = db.Column(db.String(16), nullable=False, default="internal")
    label = db.Column(db.String(128))                 # e.g. "DRAGEN 4.3 (provider X)"
    path = db.Column(db.String(1024), nullable=False)
    size_bytes = db.Column(db.BigInteger)
    mtime = db.Column(db.Float)                       # modification time at registration
    sha256 = db.Column(db.String(64))
    checksum_at = db.Column(db.DateTime)
    detected_build = db.Column(db.String(32))
    validation_status = db.Column(db.String(16), default="unchecked")  # ok/warning/error/unchecked
    validation_message = db.Column(db.Text)
    registered_at = db.Column(db.DateTime, default=utcnow, nullable=False)

    sample = db.relationship("Sample", back_populates="files")

    @property
    def role_label(self):
        return FILE_ROLES.get(self.role, self.role)
