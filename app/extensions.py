"""Shared extension objects (created once, bound to the app in create_app)."""
from celery import Celery, Task
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import event
from sqlalchemy.engine import Engine

db = SQLAlchemy()
_state = {"flask_app": None}


class FlaskTask(Task):
    """Run every task inside an application context of the *current* Flask app."""

    def __call__(self, *args, **kwargs):
        with _state["flask_app"].app_context():
            return self.run(*args, **kwargs)


celery_app = Celery("ngs_platform", task_cls=FlaskTask)


@event.listens_for(Engine, "connect")
def _sqlite_pragmas(dbapi_connection, _record):
    """WAL mode lets the web app read while the worker writes."""
    if dbapi_connection.__class__.__module__.startswith("sqlite3"):
        cur = dbapi_connection.cursor()
        cur.execute("PRAGMA journal_mode=WAL")
        cur.execute("PRAGMA busy_timeout=30000")
        cur.execute("PRAGMA foreign_keys=ON")
        cur.close()


def init_celery(flask_app) -> Celery:
    _state["flask_app"] = flask_app
    celery_app.config_from_object(flask_app.config["CELERY"])
    celery_app.set_default()
    flask_app.extensions["celery"] = celery_app
    return celery_app
