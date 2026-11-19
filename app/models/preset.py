from ..extensions import db
from .base import utcnow


class ParamPreset(db.Model):
    """A named, reusable set of workflow parameters (e.g. "In-house default", "Provider-like").

    ``params`` may be partial: only the keys it contains are applied to the Start-workflow
    form. Every run records which preset it started from and which values were changed
    afterwards (run params["preset"]).
    """

    __tablename__ = "param_presets"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(80), unique=True, nullable=False)
    description = db.Column(db.Text)
    origin = db.Column(db.String(16), nullable=False, default="user")   # builtin / user / imported / run
    source = db.Column(db.String(1024))          # imported: VCF path; run: run ID
    locked = db.Column(db.Boolean, nullable=False, default=False)      # built-ins cannot be edited/deleted
    params = db.Column(db.JSON, nullable=False, default=dict)
    notes = db.Column(db.JSON)                   # imported: mapping notes (how each value was derived)
    profile = db.Column(db.JSON)                 # imported: the filter profile extracted from the VCF
    created_at = db.Column(db.DateTime, default=utcnow, nullable=False)
    updated_at = db.Column(db.DateTime, default=utcnow, onupdate=utcnow, nullable=False)

    ORIGINS = {"builtin": "Built-in", "user": "User", "imported": "Imported from VCF", "run": "From a run"}

    def __repr__(self):
        return f"<Preset {self.name}>"
