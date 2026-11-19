from ..extensions import db
from .base import utcnow


class ReferenceGenome(db.Model):
    """A registered reference FASTA plus the result of its last integrity check."""

    __tablename__ = "reference_genomes"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(64), unique=True, nullable=False)
    fasta_path = db.Column(db.String(1024), nullable=False)
    build = db.Column(db.String(32))             # detected: GRCh38 / GRCh37 / unknown
    naming = db.Column(db.String(16))            # "ucsc" (chr1) / "ensembl" (1)
    contig_count = db.Column(db.Integer)
    alt_contig_count = db.Column(db.Integer)
    is_default = db.Column(db.Boolean, default=False, nullable=False)
    status = db.Column(db.String(16), default="unchecked")  # ok / warning / error / unchecked
    check_report = db.Column(db.JSON)            # full report from reference_service
    last_checked = db.Column(db.DateTime)
    notes = db.Column(db.Text)
    created_at = db.Column(db.DateTime, default=utcnow, nullable=False)

    def __repr__(self):
        return f"<Reference {self.name} {self.build}>"
