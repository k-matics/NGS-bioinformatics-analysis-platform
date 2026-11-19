from .api import bp as api_bp
from .dashboard import bp as dashboard_bp
from .modules import bp as modules_bp
from .presets import bp as presets_bp
from .samples import bp as samples_bp
from .settings import bp as settings_bp
from .workflows import bp as workflows_bp


def register_blueprints(app):
    for bp in (dashboard_bp, samples_bp, workflows_bp, settings_bp, modules_bp, presets_bp, api_bp):
        app.register_blueprint(bp)
