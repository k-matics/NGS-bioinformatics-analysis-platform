"""Building blocks for workflows.

A *Step* is one scientifically meaningful unit of work (e.g. "FASTQ integrity
check", "Duplicate marking"). It declares:

* ``applicable(ctx)``  - should it run for this sample/inputs? (else SKIPPED)
* ``execute(ctx, runner)`` - do the work; external tools go through ``runner.run``
* ``validate(ctx, result)`` - check the outputs are technically valid.
  A command exiting with code 0 is NOT enough (e.g. a truncated BAM).

Steps must write only inside ``ctx.step_dir(...)`` so every run is isolated
and nothing is overwritten between runs.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Callable

from .runner import CommandRunner, StepCancelled


class StepFailure(Exception):
    """Raised by a step for a scientific/technical failure with a clear message."""


@dataclass
class StepResult:
    outputs: dict[str, str] = field(default_factory=dict)   # name -> absolute path
    summary: dict = field(default_factory=dict)              # small metrics for the UI/report
    warnings: list[str] = field(default_factory=list)


@dataclass
class StepContext:
    run: object            # WorkflowRun
    sample: object         # Sample
    reference: object      # ReferenceGenome | None
    params: dict
    output_dir: Path
    threads: int
    input_files: list      # SampleFile records selected for this run
    log_path: Path | None = None
    stop_requested: Callable[[], bool] = lambda: False
    shared: dict = field(default_factory=dict)  # results passed between steps in one execution

    def step_dir(self, name: str) -> Path:
        d = self.output_dir / name
        d.mkdir(parents=True, exist_ok=True)
        return d

    def log(self, message: str) -> None:
        if self.log_path:
            with open(self.log_path, "a") as fh:
                fh.write(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {message}\n")

    def check_cancel(self) -> None:
        """Call regularly inside long pure-Python loops."""
        if self.stop_requested():
            self.log("Cancelled by user.")
            raise StepCancelled()

    def files(self, *roles: str) -> list:
        return [f for f in self.input_files if f.role in roles]


class Step:
    key: str = ""
    name: str = ""
    description: str = ""
    tools: tuple[str, ...] = ()   # recorded for provenance

    def applicable(self, ctx: StepContext) -> tuple[bool, str]:
        return True, ""

    def execute(self, ctx: StepContext, runner: CommandRunner) -> StepResult:
        raise NotImplementedError

    def validate(self, ctx: StepContext, result: StepResult) -> None:
        """Default: every declared output must exist and be non-empty."""
        for name, path in result.outputs.items():
            p = Path(path)
            if not p.exists() or (p.is_file() and p.stat().st_size == 0):
                raise StepFailure(f"Expected output '{name}' missing or empty: {p}")

    def outputs_still_valid(self, ctx: StepContext, recorded_outputs: dict) -> bool:
        """Used on resume: can a previously completed step be trusted?"""
        try:
            self.validate(ctx, StepResult(outputs=recorded_outputs))
            return True
        except StepFailure:
            return False
