from ..extensions import db
from .base import utcnow


class AppSetting(db.Model):
    __tablename__ = "app_settings"

    key = db.Column(db.String(64), primary_key=True)
    value = db.Column(db.JSON)
    updated_at = db.Column(db.DateTime, default=utcnow, onupdate=utcnow)
