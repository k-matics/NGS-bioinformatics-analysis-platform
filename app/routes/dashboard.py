from flask import Blueprint, render_template
from sqlalchemy import func

from ..extensions import db
from ..models import Sample, WorkflowEvent, WorkflowRun
from ..models import status as S
from ..services.limitations import platform_limitations

bp = Blueprint("dashboard", __name__)


def run_counts() -> dict:
    rows = db.session.query(WorkflowRun.status, func.count()).group_by(WorkflowRun.status).all()
    counts = {s: 0 for s in S.RUN_STATUSES}
    counts.update(dict(rows))
    return counts


@bp.route("/")
def index():
    counts = run_counts()
    active = (WorkflowRun.query.filter(WorkflowRun.status.in_([S.RUNNING, S.PAUSED, S.QUEUED, S.PENDING]))
              .order_by(WorkflowRun.created_at).all())
    events = WorkflowEvent.query.order_by(WorkflowEvent.id.desc()).limit(15).all()
    samples_processed = (db.session.query(func.count(func.distinct(WorkflowRun.sample_pk)))
                         .filter(WorkflowRun.status == S.COMPLETED).scalar())
    return render_template(
        "dashboard.html", counts=counts, active=active, events=events,
        total_samples=Sample.query.count(), samples_processed=samples_processed,
        limitations=platform_limitations(),
    )
