"""Parameter presets: validation, built-ins, applying a preset, recording it with a run.

A preset is a (possibly partial) parameter dict in exactly the shape a run stores:

    {"trimming": true, "callers": ["gatk_hc", "deepvariant"],
     "caller_params": {"gatk_hc": {"filter_min_depth": 2, ...}}}

Every value is validated against the same ParamSpecs that validate the Start form, so a
preset can never put a value into a command line that the form itself would refuse.
"""
from __future__ import annotations

import copy

from ..engine.params import ParamError, ParamSpec
from ..extensions import db
from ..models import ParamPreset

DEFAULT_PRESET = "In-house default"
RAW_PRESET = "Raw tool defaults (no filters)"


def _top_level_specs() -> dict[str, ParamSpec]:
    """Union of the parameters of every registered workflow (same key = same meaning)."""
    from ..engine.registry import all_workflows
    from ..pipelines.germline import form_specs
    out: dict[str, ParamSpec] = {}
    for w in all_workflows():
        for s in form_specs(w):
            out.setdefault(s.key, s)
    return out


def _caller_specs() -> dict[str, dict[str, ParamSpec]]:
    from ..pipelines.variant_calling.callers import all_callers
    return {c.key: {s.key: s for s in c.params} for c in all_callers("germline")}


def _check(spec: ParamSpec, value, where: str):
    try:
        if spec.type == "bool":
            if not isinstance(value, bool):
                raise ParamError(f"{spec.label}: must be true/false")
            return value
        if spec.type == "multichoice":
            if not isinstance(value, list):
                raise ParamError(f"{spec.label}: must be a list")
            return spec.parse([str(v) for v in value], True)
        if spec.type in ("int", "float"):
            if isinstance(value, bool) or not isinstance(value, (int, float, str)):
                raise ParamError(f"{spec.label}: must be a number")
            return spec.parse(str(value), True)
        return spec.parse(str(value), True)
    except ParamError as exc:
        raise ParamError(f"{where}{exc}") from None


def validate_preset_params(params: dict) -> tuple[dict, list[str]]:
    """Returns (cleaned params, problems). Unknown keys are problems, never silently dropped."""
    problems: list[str] = []
    clean: dict = {}
    if not isinstance(params, dict):
        return {}, ["Preset parameters must be a JSON object."]
    top = _top_level_specs()
    callers = _caller_specs()
    for key, value in params.items():
        if key == "caller_params":
            if not isinstance(value, dict):
                problems.append("caller_params must be an object {caller: {param: value}}.")
                continue
            cp = {}
            for ck, cvals in value.items():
                if ck not in callers:
                    problems.append(f"caller_params: unknown caller '{ck}'.")
                    continue
                if not isinstance(cvals, dict):
                    problems.append(f"caller_params.{ck} must be an object.")
                    continue
                cp[ck] = {}
                for pk, pv in cvals.items():
                    spec = callers[ck].get(pk)
                    if spec is None:
                        problems.append(f"caller_params.{ck}: unknown parameter '{pk}'.")
                        continue
                    try:
                        cp[ck][pk] = _check(spec, pv, f"{ck}: ")
                    except ParamError as exc:
                        problems.append(str(exc))
            clean["caller_params"] = cp
            continue
        spec = top.get(key)
        if spec is None:
            problems.append(f"Unknown parameter '{key}'.")
            continue
        try:
            clean[key] = _check(spec, value, "")
        except ParamError as exc:
            problems.append(str(exc))
    return clean, problems


def in_house_defaults() -> dict:
    """The platform's own defaults for a full germline run with all four callers."""
    from ..engine.registry import get_workflow
    from ..pipelines.germline import form_specs
    from ..pipelines.variant_calling.callers import all_callers
    w = get_workflow("fastq_to_vcf")
    params = {s.key: copy.deepcopy(s.default) for s in form_specs(w) if s.key != "callers"}
    keys = [c.key for c in all_callers("germline")]
    params["callers"] = keys
    params["caller_params"] = {c.key: {s.key: s.default for s in c.params} for c in all_callers("germline")}
    return params


def raw_tool_defaults() -> dict:
    """Every caller with its own defaults (bare command) and every soft filter off.
    Only caller settings are set; trimming / duplicates / regions keep the form's values."""
    from ..pipelines.variant_calling.callers import all_callers
    cp = {}
    for c in all_callers("germline"):
        vals = {}
        for s in c.params:
            if s.key == "tool_defaults":
                vals[s.key] = True
            elif s.key.startswith("filter_"):
                vals[s.key] = False if s.type == "bool" else 0
        cp[c.key] = vals
    return {"caller_params": cp}


BUILTINS = [
    (DEFAULT_PRESET, in_house_defaults,
     "The platform's standard settings (all four callers, per-caller filters: depth ≥ 10, GQ ≥ 20, "
     "VAF ≥ 0.15, GATK hard filters, FreeBayes strand filter). Locked; duplicate it to make a variant."),
    (RAW_PRESET, raw_tool_defaults,
     "Each caller runs its bare command with the tool's own defaults (e.g. GATK: gatk HaplotypeCaller "
     "-R ref -I bam -O vcf, plus the calling regions and thread count) and no soft filter is applied "
     "(only NoAltGT labels 0/0 records). Use it to see what a tool calls out of the box. Note: tool "
     "defaults still ignore low-MAPQ reads (GATK MAPQ < 20, DeepVariant < 5, FreeBayes < 1). Locked."),
]


def ensure_builtin_presets() -> None:
    added = False
    for name, build, desc in BUILTINS:
        if ParamPreset.query.filter_by(name=name).first():
            continue
        params, problems = validate_preset_params(build())
        if problems:   # a bug in the defaults themselves - do not hide it
            raise RuntimeError(f"Built-in preset '{name}' invalid: " + "; ".join(problems))
        db.session.add(ParamPreset(name=name, origin="builtin", locked=True, params=params, description=desc))
        added = True
    if added:
        db.session.commit()


def applicable_params(preset_params: dict, workflow_keys: set[str]) -> dict:
    """The part of a preset that a given workflow understands."""
    out = {k: v for k, v in preset_params.items() if k in workflow_keys}
    if "callers" in workflow_keys and "caller_params" in preset_params:
        out["caller_params"] = preset_params["caller_params"]
    return out


def preset_diff(preset_params: dict, final: dict) -> list[str]:
    """Keys whose submitted value differs from the preset (the user changed them after applying)."""
    changed = []
    for k, v in preset_params.items():
        if k == "caller_params":
            for ck, cv in (v or {}).items():
                if ck not in (final.get("callers") or []):
                    continue
                for pk, pv in cv.items():
                    if (final.get("caller_params") or {}).get(ck, {}).get(pk) != pv:
                        changed.append(f"{ck}.{pk}")
        elif k in final and final[k] != v:
            changed.append(k)
    return changed


def record_on_params(preset: ParamPreset, workflow_keys: set[str], params: dict) -> None:
    applied = applicable_params(preset.params or {}, workflow_keys)
    params["preset"] = {"id": preset.id, "name": preset.name, "origin": preset.origin,
                        "preset_updated": preset.updated_at.isoformat() if preset.updated_at else None,
                        "changed_after_applying": preset_diff(applied, params)}


def flatten_for_form(preset_params: dict) -> dict:
    """{form field suffix: value} as used by the Start form JS: 'key' and 'caller.param'."""
    out = {}
    for k, v in preset_params.items():
        if k == "caller_params":
            for ck, cv in v.items():
                for pk, pv in cv.items():
                    out[f"c:{ck}.{pk}"] = pv
        else:
            out[f"p:{k}"] = v
    return out
