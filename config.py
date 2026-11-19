"""Central configuration for the NGS platform.

Values here are *defaults*. Every important value can be overridden with an
environment variable (prefix ``NGS_``) and the runtime-tunable ones (threads,
retention policy, tool paths...) can also be changed from the Settings page,
which stores them in the database (see ``app/services/settings_service.py``).

No credentials (SSH or otherwise) are ever stored here.
"""
from __future__ import annotations

import os
import secrets
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent


def _load_local_env(path: Path = BASE_DIR / "config" / "local.env") -> None:
    """Server-specific settings (config/local.env, not in git) for EVERY entry point: web app, worker,
    scripts and shells. Real environment variables win; tests set NGS_SKIP_LOCAL_ENV=1."""
    if os.environ.get("NGS_SKIP_LOCAL_ENV") or not path.is_file():
        return
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        os.environ.setdefault(key.strip(), value.strip().strip('"').strip("'"))


_load_local_env()


def _path(env_name: str, default: str | Path) -> Path:
    return Path(os.environ.get(env_name, str(default))).expanduser().resolve()


def _secret_key() -> str:
    """Read (or create once) a random secret key stored outside the code."""
    env_key = os.environ.get("NGS_SECRET_KEY")
    if env_key:
        return env_key
    key_file = BASE_DIR / "instance" / "secret_key"
    key_file.parent.mkdir(parents=True, exist_ok=True)
    if not key_file.exists():
        key_file.write_text(secrets.token_hex(32))
        key_file.chmod(0o600)
    return key_file.read_text().strip()


class Config:
    BASE_DIR = BASE_DIR
    SECRET_KEY = _secret_key()

    # --- Database -----------------------------------------------------------
    # SQLite now; set NGS_DATABASE_URL=postgresql+psycopg2://... to switch.
    DATABASE_DIR = _path("NGS_DATABASE_DIR", BASE_DIR / "database")
    SQLALCHEMY_DATABASE_URI = os.environ.get(
        "NGS_DATABASE_URL", f"sqlite:///{DATABASE_DIR / 'ngs_platform.db'}"
    )
    SQLALCHEMY_TRACK_MODIFICATIONS = False

    # --- Directories --------------------------------------------------------
    RESULTS_DIR = _path("NGS_RESULTS_DIR", BASE_DIR / "results")
    UPLOAD_DIR = _path("NGS_UPLOAD_DIR", BASE_DIR / "uploads")
    LOG_DIR = _path("NGS_LOG_DIR", BASE_DIR / "logs")

    # Only files below these roots can be registered/browsed (colon separated).
    DATA_ROOTS = [
        Path(p).expanduser().resolve()
        for p in os.environ.get("NGS_DATA_ROOTS", "~").split(":")
        if p
    ]

    # --- Reference genome (the only hard-coded data path, per requirements) --
    DEFAULT_REFERENCE_NAME = "hg38_ucsc"
    DEFAULT_REFERENCE_FASTA = _path("NGS_REFERENCE_FASTA", "~/reference_genome/hg38.fa")

    # --- Job queue ----------------------------------------------------------
    REDIS_URL = os.environ.get("NGS_REDIS_URL", "redis://127.0.0.1:6379")
    CELERY = {
        "broker_url": f"{REDIS_URL}/0",
        "result_backend": f"{REDIS_URL}/1",
        "task_ignore_result": True,
        # Long bioinformatics jobs: take one task at a time, never re-deliver
        # a task that is still running (that would cause duplicate execution).
        "worker_prefetch_multiplier": 1,
        "task_acks_late": False,
        "broker_transport_options": {"visibility_timeout": 7 * 24 * 3600},
        "task_default_queue": "workflows",
        "broker_connection_retry_on_startup": True,
    }

    # --- Runtime defaults (overridable from Settings page) --------------------
    DEFAULTS = {
        "threads_per_job": 8,
        "min_free_disk_gb": 3,
        "retain_intermediates": True,
        "gatk_path": os.environ.get("NGS_GATK_PATH", "gatk"),   # launcher path, or "gatk" on PATH
        # empty = auto-detect a Java 17-21 under /usr/lib/jvm (GATK 4.6 does not support Java 25)
        "gatk_java_home": "",
        "fastq_check_records": 100_000,
        # Phase 4/5
        "deepvariant_image": "google/deepvariant:1.10.0",
        "rtg_sdf_dir": "",          # empty = <instance>/rtg_sdf (RTG's binary copy of the reference)
        "rtg_mem_gb": 8,
    }


class TestConfig(Config):
    TESTING = True
    SQLALCHEMY_DATABASE_URI = "sqlite://"
    CELERY = dict(Config.CELERY, task_always_eager=True, task_eager_propagates=True)
