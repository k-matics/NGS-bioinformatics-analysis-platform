"""Runtime settings: config defaults overridden by values saved from the Settings page."""
from __future__ import annotations

from flask import current_app

from ..extensions import db
from ..models import AppSetting


def get_setting(key: str, default=None):
    row = db.session.get(AppSetting, key)
    if row is not None:
        return row.value
    return current_app.config["DEFAULTS"].get(key, default)


def all_settings() -> dict:
    values = dict(current_app.config["DEFAULTS"])
    for row in AppSetting.query.all():
        values[row.key] = row.value
    return values


def set_setting(key: str, value) -> None:
    row = db.session.get(AppSetting, key)
    if row is None:
        row = AppSetting(key=key, value=value)
        db.session.add(row)
    else:
        row.value = value
    db.session.commit()
