"""Celery tasks (executed by the worker process, not the web server)."""
from .extensions import celery_app


@celery_app.task(name="ngs.run_workflow", bind=True)
def run_workflow(self, run_pk: int) -> str:
    from .engine.executor import execute_run
    from .extensions import db
    from .models import WorkflowRun

    run = db.session.get(WorkflowRun, run_pk)
    if run is not None and self.request.id and not run.celery_task_id:
        run.celery_task_id = self.request.id
        db.session.commit()
    return execute_run(run_pk)
