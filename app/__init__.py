"""Application factory."""
from __future__ import annotations

from datetime import datetime, timezone

from flask import Flask

from .extensions import db, init_celery
from .utils.logging_setup import configure_logging


def create_app(config_object: str | type = "config.Config", overrides: dict | None = None) -> Flask:
    app = Flask(__name__, instance_relative_config=False)
    app.config.from_object(config_object)
    if overrides:
        app.config.update(overrides)

    for key in ("RESULTS_DIR", "UPLOAD_DIR", "LOG_DIR"):
        app.config[key].mkdir(parents=True, exist_ok=True)
    if app.config["SQLALCHEMY_DATABASE_URI"].startswith("sqlite:///"):
        app.config["DATABASE_DIR"].mkdir(parents=True, exist_ok=True)

    configure_logging(app)
    db.init_app(app)
    init_celery(app)

    from . import models  # noqa: F401  (register tables)
    from .routes import register_blueprints

    register_blueprints(app)
    _register_template_helpers(app)

    with app.app_context():
        db.create_all()
        from .services.reference_service import ensure_default_reference

        ensure_default_reference()
        try:
            from .services.preset_service import ensure_builtin_presets
            ensure_builtin_presets()
        except Exception:  # never block start-up; the Presets page retries and shows the error
            app.logger.exception("Could not create the built-in parameter preset")

    return app


def _register_template_helpers(app: Flask) -> None:
    from .models.status import STATUS_BADGE

    @app.template_filter("dt")
    def format_dt(value):
        if not value:
            return "—"
        if value.tzinfo is None:
            value = value.replace(tzinfo=timezone.utc)
        return value.astimezone().strftime("%Y-%m-%d %H:%M:%S")

    @app.template_filter("duration")
    def format_duration(seconds):
        if seconds is None:
            return "—"
        seconds = int(seconds)
        h, rem = divmod(seconds, 3600)
        m, s = divmod(rem, 60)
        return f"{h}h {m:02d}m {s:02d}s" if h else f"{m}m {s:02d}s"

    @app.template_filter("filesize")
    def format_size(num):
        if num is None:
            return "—"
        for unit in ("B", "KB", "MB", "GB", "TB"):
            if abs(num) < 1024:
                return f"{num:.1f} {unit}" if unit != "B" else f"{num} B"
            num /= 1024
        return f"{num:.1f} PB"

    @app.context_processor
    def inject_globals():
        return {"STATUS_BADGE": STATUS_BADGE, "now": datetime.now(timezone.utc)}
