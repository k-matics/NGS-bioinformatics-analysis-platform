"""Application-level logging (per-workflow logs are written by the engine)."""
import logging
from logging.handlers import RotatingFileHandler


def configure_logging(app) -> None:
    if app.config.get("TESTING"):
        return
    log_dir = app.config["LOG_DIR"]
    handler = RotatingFileHandler(log_dir / "app.log", maxBytes=10_000_000, backupCount=10)
    handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)s [%(process)d] %(name)s: %(message)s")
    )
    handler.setLevel(logging.INFO)
    root = logging.getLogger()
    if not any(isinstance(h, RotatingFileHandler) for h in root.handlers):
        root.addHandler(handler)
    root.setLevel(logging.INFO)
