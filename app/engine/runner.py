"""Safe execution of external commands.

Rules enforced here (clinical traceability + security):
* Commands are argument lists, never shell strings -> no shell injection.
  A pipe (``bwa mem ... | samtools sort ...``) is expressed as a list of
  argument lists and connected with OS pipes, still without a shell.
* The executable is resolved to an absolute path and logged.
* stdout/stderr go to the step log file; nothing is hidden.
* Every process of a pipeline gets its own exit code; *any* non-zero exit
  fails the step (like ``set -o pipefail``).
* All processes run in one new process group, so Cancel/Pause can signal
  the whole pipeline (SIGTERM -> SIGKILL, SIGSTOP / SIGCONT).
"""
from __future__ import annotations

import os
import shlex
import shutil
import signal
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import Callable, Sequence

Argv = Sequence[str | os.PathLike]


class StepCancelled(Exception):
    """Raised when the user cancelled the workflow."""


class CommandFailed(Exception):
    def __init__(self, display: str, exit_codes: list[int], log_tail: str = ""):
        self.display = display
        self.exit_codes = exit_codes
        self.log_tail = log_tail
        super().__init__(f"Command failed (exit codes {exit_codes}): {display}")


def _normalise(commands) -> list[list[str]]:
    if isinstance(commands, (str, bytes)):
        raise TypeError("Commands must be argument lists, not shell strings.")
    if commands and isinstance(commands[0], (str, os.PathLike)):
        commands = [commands]  # a single command
    pipeline = []
    for argv in commands:
        if isinstance(argv, (str, bytes)) or not argv:
            raise TypeError("Each command must be a non-empty argument list.")
        if not all(isinstance(a, (str, os.PathLike, int, float)) for a in argv):
            raise TypeError(f"Invalid argument type in {argv!r}")
        pipeline.append([str(a) for a in argv])
    return pipeline


def display_command(pipeline: list[list[str]], stdout_path=None, stdin_path=None) -> str:
    text = " | ".join(shlex.join(c) for c in pipeline)
    if stdin_path:
        text += f" < {shlex.quote(str(stdin_path))}"
    if stdout_path:
        text += f" > {shlex.quote(str(stdout_path))}"
    return text


def tail_file(path: Path, max_bytes: int = 4000) -> str:
    try:
        with open(path, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            fh.seek(max(0, size - max_bytes))
            return fh.read().decode(errors="replace")
    except OSError:
        return ""


class CommandRunner:
    POLL_SECONDS = 0.5
    TERM_GRACE_SECONDS = 15

    def __init__(
        self,
        log_path: Path,
        control: Callable[[], str | None] | None = None,
        on_start: Callable[[list[list[str]], str, int], object] | None = None,
        on_finish: Callable[[object, list[int]], None] | None = None,
        on_pause_change: Callable[[bool], None] | None = None,
    ):
        """
        control()        -> None (keep running) | "pause" | "cancel"
        on_start(argv, display, pid) -> handle (e.g. DB record)
        on_finish(handle, exit_codes)
        """
        self.log_path = Path(log_path)
        self.control = control
        self.on_start = on_start
        self.on_finish = on_finish
        self.on_pause_change = on_pause_change

    def log(self, message: str) -> None:
        with open(self.log_path, "a") as fh:
            fh.write(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {message}\n")

    def run(
        self,
        commands,
        *,
        stdout_path: Path | None = None,
        stdin_path: Path | None = None,
        cwd: Path | None = None,
        env: dict | None = None,
        allow_sigpipe: bool = True,
        ok_exit_codes: Sequence[int] = (0,),
    ) -> list[int]:
        pipeline = _normalise(commands)
        for argv in pipeline:
            exe = shutil.which(argv[0])
            if exe is None:
                raise FileNotFoundError(f"Tool not found on PATH: {argv[0]}")
            argv[0] = exe
        display = display_command(pipeline, stdout_path, stdin_path)
        self.log(f"$ {display}")

        procs: list[subprocess.Popen] = []
        stdin_fh = open(stdin_path, "rb") if stdin_path else subprocess.DEVNULL
        out_fh = open(stdout_path, "wb") if stdout_path else None
        log_fh = open(self.log_path, "ab")
        handle = None
        try:
            pgid = 0  # 0 = start a new process group led by the first process
            prev = stdin_fh
            for i, argv in enumerate(pipeline):
                last = i == len(pipeline) - 1
                proc = subprocess.Popen(
                    argv,
                    stdin=prev,
                    stdout=(out_fh or log_fh) if last else subprocess.PIPE,
                    stderr=log_fh,
                    cwd=cwd,
                    env=env,
                    process_group=pgid,
                )
                if i == 0:
                    pgid = proc.pid
                else:
                    prev.close()  # parent's copy of the pipe; lets SIGPIPE propagate
                prev = proc.stdout
                procs.append(proc)
            if self.on_start:
                handle = self.on_start(pipeline, display, pgid)
            codes = self._wait(procs, pgid)
        finally:
            if stdin_fh is not subprocess.DEVNULL:
                stdin_fh.close()
            if out_fh:
                out_fh.close()
            log_fh.close()

        if self.on_finish:
            self.on_finish(handle, codes)
        self.log(f"exit codes: {codes}")
        bad = [
            c for i, c in enumerate(codes)
            if c not in ok_exit_codes
            and not (allow_sigpipe and c == -signal.SIGPIPE and i < len(codes) - 1)
        ]
        if bad:
            raise CommandFailed(display, codes, tail_file(self.log_path))
        return codes

    # ------------------------------------------------------------------
    def _wait(self, procs: list[subprocess.Popen], pgid: int) -> list[int]:
        paused = False
        while True:
            codes = [p.poll() for p in procs]
            if all(c is not None for c in codes):
                return codes
            request = self.control() if self.control else None
            if request == "cancel":
                self._terminate(procs, pgid)
                self.log("Cancelled by user: process group terminated.")
                raise StepCancelled()
            if request == "pause" and not paused:
                os.killpg(pgid, signal.SIGSTOP)
                paused = True
                self.log("Paused (SIGSTOP sent to process group).")
                if self.on_pause_change:
                    self.on_pause_change(True)
            elif request != "pause" and paused:
                os.killpg(pgid, signal.SIGCONT)
                paused = False
                self.log("Resumed (SIGCONT sent to process group).")
                if self.on_pause_change:
                    self.on_pause_change(False)
            time.sleep(self.POLL_SECONDS)

    def _terminate(self, procs, pgid) -> None:
        try:
            os.killpg(pgid, signal.SIGCONT)  # a stopped process cannot handle SIGTERM
            os.killpg(pgid, signal.SIGTERM)
        except ProcessLookupError:
            return
        deadline = time.time() + self.TERM_GRACE_SECONDS
        while time.time() < deadline and any(p.poll() is None for p in procs):
            time.sleep(0.2)
        if any(p.poll() is None for p in procs):
            try:
                os.killpg(pgid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        for p in procs:
            p.wait()
