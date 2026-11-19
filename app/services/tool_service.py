"""Locate external tools and record their versions (for provenance).

The platform runs inside the `ngs_platform` conda environment, so tools are
resolved from PATH (conda env first). The resolved *absolute* path and the
version are stored with every command that is executed.
"""
from __future__ import annotations

import re
import shutil
import subprocess
from functools import lru_cache

# tool -> (version command, regex capturing the version)
TOOL_SPECS: dict[str, tuple[list[str], str]] = {
    "samtools": (["samtools", "--version"], r"samtools ([\d.]+)"),
    "bcftools": (["bcftools", "--version"], r"bcftools ([\d.]+)"),
    "bgzip": (["bgzip", "--version"], r"bgzip \(htslib\) ([\d.]+)"),
    "tabix": (["tabix", "--version"], r"tabix \(htslib\) ([\d.]+)"),
    "pigz": (["pigz", "--version"], r"pigz ([\d.]+)"),
    "bwa": (["bwa"], r"Version: (\S+)"),
    "fastqc": (["fastqc", "--version"], r"FastQC v([\d.]+)"),
    "multiqc": (["multiqc", "--version"], r"version ([\d.]+)"),
    "fastp": (["fastp", "--version"], r"fastp ([\d.]+)"),
    "mosdepth": (["mosdepth", "--version"], r"mosdepth ([\d.]+)"),
    "bedtools": (["bedtools", "--version"], r"bedtools v([\d.]+)"),
    "freebayes": (["freebayes", "--version"], r"version:\s+v?([\d.]+)"),
    "rtg": (["rtg", "version"], r"Product: RTG Tools ([\d.]+)"),
    "java": (["java", "-version"], r'version "([^"]+)"'),
    "docker": (["docker", "--version"], r"Docker version ([\d.]+)"),
}
# Tools listed on the Settings page, with the phase that first needs them.
TOOL_PHASE = {
    "samtools": 2, "bcftools": 2, "pigz": 2, "bgzip": 2, "tabix": 2,
    "fastqc": 3, "multiqc": 3, "fastp": 3, "bwa": 3, "mosdepth": 3, "bedtools": 3,
    "gatk": 3, "java": 3, "freebayes": 4, "docker": 4, "rtg": 5,
}


def resolve(tool: str, override: str | None = None) -> str | None:
    if override:
        return override if shutil.which(override) else None
    return shutil.which(tool)


@lru_cache(maxsize=64)
def tool_version(tool: str, path: str | None = None) -> str | None:
    """Return the version string for a tool, or None if not found."""
    exe = path or resolve(tool)
    if not exe:
        return None
    if tool == "gatk":
        cmd, pattern = [exe, "--version"], r"\(GATK\) v([\d.]+)"
    elif tool in TOOL_SPECS:
        cmd, pattern = TOOL_SPECS[tool]
        cmd = [exe, *cmd[1:]]
    else:
        return "unknown"
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, errors="replace", timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return None
    m = re.search(pattern, proc.stdout + proc.stderr)
    return m.group(1) if m else "unknown"


def tool_inventory(gatk_path: str | None = None) -> list[dict]:
    rows = []
    for tool, phase in TOOL_PHASE.items():
        path = resolve(tool, gatk_path if tool == "gatk" else None)
        rows.append({
            "tool": tool,
            "phase": phase,
            "path": path,
            "version": tool_version(tool, path) if path else None,
        })
    return rows


# --------------------------------------------------------------------------- Java for GATK
GATK_JAVA_MAJOR = (17, 18, 19, 20, 21)   # GATK 4.6 is built and tested on Java 17


def _java_major(java: str) -> int | None:
    try:
        proc = subprocess.run([java, "-version"], capture_output=True, text=True, errors="replace", timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return None
    m = re.search(r'version "(\d+)', proc.stderr + proc.stdout)
    return int(m.group(1)) if m else None


@lru_cache(maxsize=4)
def find_gatk_java_home(configured: str | None = None) -> str | None:
    """Locate a Java 17-21 installation for GATK.

    The conda environment ships Java 25 (needed by nothing else in particular),
    which GATK 4.6 does not support. Order: the 'gatk_java_home' setting, then
    standard system JVM folders. Nothing is installed or changed.
    """
    from pathlib import Path
    candidates = [Path(configured)] if configured else []
    jvm_root = Path("/usr/lib/jvm")
    if jvm_root.is_dir():
        candidates += sorted(p for p in jvm_root.iterdir() if p.is_dir() and not p.is_symlink())
    for home in candidates:
        java = home / "bin" / "java"
        if java.is_file() and _java_major(str(java)) in GATK_JAVA_MAJOR:
            return str(home)
    return None


def gatk_env(configured_java_home: str | None = None) -> dict:
    """Environment for running GATK with a supported Java (PATH + JAVA_HOME)."""
    import os
    home = find_gatk_java_home(configured_java_home or None)
    if home is None:
        raise RuntimeError("No Java 17-21 found for GATK (set 'gatk_java_home' in Settings).")
    env = dict(os.environ)
    env["JAVA_HOME"] = home
    env["PATH"] = f"{home}/bin:{env.get('PATH', '')}"
    return env
