"""The command runner is the security- and integrity-critical core."""
import threading
import time

import pytest

from app.engine.runner import CommandFailed, CommandRunner, StepCancelled


def test_success_and_log(tmp_path):
    log = tmp_path / "step.log"
    starts, finishes = [], []
    runner = CommandRunner(log, on_start=lambda a, d, p: starts.append(d) or "h",
                           on_finish=lambda h, c: finishes.append(c))
    assert runner.run(["echo", "hello"]) == [0]
    text = log.read_text()
    assert "hello" in text and "$ " in text and "exit codes: [0]" in text
    assert finishes == [[0]] and "echo hello" in starts[0]


def test_failure_records_exit_code(tmp_path):
    runner = CommandRunner(tmp_path / "s.log")
    with pytest.raises(CommandFailed) as exc:
        runner.run(["sh", "-c", "echo boom >&2; exit 3"])
    assert exc.value.exit_codes == [3]
    assert "boom" in exc.value.log_tail


def test_pipeline_and_stdout_file(tmp_path):
    out = tmp_path / "out.txt"
    CommandRunner(tmp_path / "s.log").run([["printf", "a\\nb\\nc\\n"], ["grep", "b"]], stdout_path=out)
    assert out.read_text() == "b\n"


def test_pipefail_semantics(tmp_path):
    """A failure in the *first* command of a pipe must not be hidden."""
    with pytest.raises(CommandFailed) as exc:
        CommandRunner(tmp_path / "s.log").run([["sh", "-c", "exit 2"], ["cat"]])
    assert exc.value.exit_codes[0] == 2


def test_shell_strings_rejected(tmp_path):
    with pytest.raises(TypeError):
        CommandRunner(tmp_path / "s.log").run("echo hi; rm -rf /")


def test_no_shell_interpretation(tmp_path):
    marker = tmp_path / "pwned"
    CommandRunner(tmp_path / "s.log").run(["echo", f"x; touch {marker}"])
    assert not marker.exists()


def test_missing_tool(tmp_path):
    with pytest.raises(FileNotFoundError):
        CommandRunner(tmp_path / "s.log").run(["definitely-not-a-tool-xyz"])


def test_cancel_kills_whole_process_group(tmp_path):
    flag = {"v": None}
    runner = CommandRunner(tmp_path / "s.log", control=lambda: flag["v"])
    runner.POLL_SECONDS = 0.1
    threading.Timer(0.5, lambda: flag.update(v="cancel")).start()
    t0 = time.time()
    with pytest.raises(StepCancelled):
        runner.run([["sleep", "30"], ["cat"]])
    assert time.time() - t0 < 10


def test_pause_and_resume(tmp_path):
    flag = {"v": None}
    changes = []
    runner = CommandRunner(tmp_path / "s.log", control=lambda: flag["v"], on_pause_change=changes.append)
    runner.POLL_SECONDS = 0.1
    threading.Timer(0.3, lambda: flag.update(v="pause")).start()
    threading.Timer(1.0, lambda: flag.update(v=None)).start()
    assert runner.run(["sleep", "1.5"]) == [0]
    assert changes == [True, False]
