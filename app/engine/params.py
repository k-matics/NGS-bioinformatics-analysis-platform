"""Typed, validated workflow parameters.

Each workflow (and each variant-caller plug-in) declares its parameters as
``ParamSpec`` objects. The Start-workflow form is generated from them, and
every submitted value is type-checked and range-checked here before it can
reach a command line. Values are passed to tools as separate argv items,
never interpolated into a shell string.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

SAFE_TEXT_RE = re.compile(r"^[A-Za-z0-9_.,:+-]*$")


class ParamError(ValueError):
    pass


@dataclass
class ParamSpec:
    key: str
    label: str
    type: str                      # bool | int | float | choice | multichoice | text
    default: object = None
    help: str = ""
    min: float | None = None
    max: float | None = None
    choices: list[tuple[str, str]] = field(default_factory=list)   # (value, label)
    group: str = "General"
    advanced: bool = False

    def parse(self, raw, present: bool):
        """Convert a raw form value (str / list[str] / None) into a typed value."""
        if self.type == "bool":
            return bool(raw) if present else False
        if self.type == "multichoice":
            values = raw or []
            allowed = {c[0] for c in self.choices}
            bad = [v for v in values if v not in allowed]
            if bad:
                raise ParamError(f"{self.label}: invalid choice(s) {bad}")
            return values
        if raw is None or raw == "":
            return self.default
        if self.type == "int":
            try:
                value = int(raw)
            except (TypeError, ValueError):
                raise ParamError(f"{self.label}: '{raw}' is not a whole number") from None
            return self._check_range(value)
        if self.type == "float":
            try:
                value = float(raw)
            except (TypeError, ValueError):
                raise ParamError(f"{self.label}: '{raw}' is not a number") from None
            return self._check_range(value)
        if self.type == "choice":
            if raw not in {c[0] for c in self.choices}:
                raise ParamError(f"{self.label}: invalid choice '{raw}'")
            return raw
        if self.type == "text":
            raw = str(raw).strip()
            if len(raw) > 2000 or not SAFE_TEXT_RE.match(raw):
                raise ParamError(f"{self.label}: only letters, digits and _ . , : + - are allowed")
            return raw
        raise ParamError(f"Unknown parameter type {self.type}")

    def _check_range(self, value):
        if self.min is not None and value < self.min:
            raise ParamError(f"{self.label}: must be ≥ {self.min}")
        if self.max is not None and value > self.max:
            raise ParamError(f"{self.label}: must be ≤ {self.max}")
        return value


def defaults(specs: list[ParamSpec]) -> dict:
    return {s.key: s.default for s in specs}


def parse_form(specs: list[ParamSpec], form, prefix: str) -> dict:
    """Read ``prefix + spec.key`` fields from a Flask form (MultiDict)."""
    out = {}
    for s in specs:
        name = prefix + s.key
        if s.type == "multichoice":
            out[s.key] = s.parse(form.getlist(name), name in form)
        else:
            out[s.key] = s.parse(form.get(name), name in form)
    return out


def validate_dict(specs: list[ParamSpec], values: dict, label: str = "") -> dict:
    """Validate a dict (e.g. from an API/test) against specs, filling defaults."""
    out = {}
    known = {s.key for s in specs}
    unknown = set(values) - known
    if unknown:
        raise ParamError(f"{label}unknown parameter(s): {sorted(unknown)}")
    for s in specs:
        if s.key not in values:
            out[s.key] = s.default
            continue
        v = values[s.key]
        if s.type == "bool":
            out[s.key] = bool(v)
        elif s.type == "multichoice":
            out[s.key] = s.parse(list(v), True)
        else:
            out[s.key] = s.parse(v if isinstance(v, str) else v, True) if s.type in ("choice", "text") \
                else s._check_range(int(v) if s.type == "int" else float(v))
    return out
