"""Celery worker entry point:

    celery -A celery_worker.celery worker --concurrency=1 -Q workflows
"""
import logging

from celery.signals import worker_ready

from app import create_app
from app.extensions import celery_app as celery

flask_app = create_app()


@worker_ready.connect
def _recover_interrupted_runs(**_kwargs):
    with flask_app.app_context():
        from app.engine.executor import mark_interrupted_runs

        n = mark_interrupted_runs()
        if n:
            logging.getLogger(__name__).warning("Marked %d interrupted run(s) as failed (resumable).", n)


import app.tasks  # noqa: E402,F401  (register tasks)
