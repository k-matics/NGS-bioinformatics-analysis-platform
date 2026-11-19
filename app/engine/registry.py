"""Registry of available workflows. Each pipeline package registers itself here."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from .params import ParamSpec
from .step import Step


@dataclass
class WorkflowDefinition:
    key: str
    name: str
    version: str
    description: str
    build_steps: Callable[[dict], list[Step]]
    default_params: dict = field(default_factory=dict)
    # returns a list of human-readable reasons why the sample cannot run this workflow
    check_inputs: Callable[[object, list], list[str]] = lambda sample, files: []
    # rough disk estimate in GB for the pre-flight check
    estimate_disk_gb: Callable[[object, list], float] = lambda sample, files: 0.1
    # typed parameters shown on the Start-workflow form
    param_specs: list[ParamSpec] = field(default_factory=list)
    # file roles this workflow consumes (None = all registered files)
    input_roles: tuple[str, ...] | None = None
    # UI module this workflow belongs to (dashboard / module pages)
    module: str = "intake"
    # parses extra, workflow-specific form fields (e.g. per-caller parameters)
    parse_extra: Callable[[object], dict] | None = None
    # semantic checks needing the reference etc.: (params, reference) -> problems
    validate_params: Callable[[dict, object], list[str]] = lambda params, reference: []


_REGISTRY: dict[str, WorkflowDefinition] = {}


def register(defn: WorkflowDefinition) -> WorkflowDefinition:
    _REGISTRY[defn.key] = defn
    return defn


def get_workflow(key: str) -> WorkflowDefinition:
    _load_pipelines()
    if key not in _REGISTRY:
        raise KeyError(f"Unknown workflow '{key}'")
    return _REGISTRY[key]


def all_workflows() -> list[WorkflowDefinition]:
    _load_pipelines()
    return list(_REGISTRY.values())


def _load_pipelines() -> None:
    # Importing a pipeline package registers its workflows.
    from ..pipelines import germline, intake  # noqa: F401
    from ..pipelines.comparison import workflows  # noqa: F401
