"""Filesystem helpers: checksums, safe directory listing, disk space."""
from __future__ import annotations

import hashlib
import os
import shutil
from pathlib import Path
from typing import Callable

from ..utils.validators import KNOWN_SUFFIXES, resolve_within_roots


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024,
                should_stop: Callable[[], bool] | None = None) -> str:
    """Stream a file through SHA-256 (works for multi-GB files, constant memory)."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for i, chunk in enumerate(iter(lambda: fh.read(chunk_size), b"")):
            digest.update(chunk)
            if should_stop and i % 16 == 0 and should_stop():
                raise InterruptedError("checksum interrupted")
    return digest.hexdigest()


def file_fingerprint(path: Path) -> dict:
    st = path.stat()
    return {"size_bytes": st.st_size, "mtime": st.st_mtime}


def list_directory(raw_path: str, roots: list[Path], show_all: bool = False) -> dict:
    """List a directory for the server-side file picker (restricted to data roots)."""
    path = resolve_within_roots(raw_path, roots)
    if not path.is_dir():
        path = path.parent
    entries = []
    try:
        with os.scandir(path) as it:
            for entry in it:
                if entry.name.startswith("."):
                    continue
                try:
                    is_dir = entry.is_dir()
                    if not is_dir and not show_all and not entry.name.lower().endswith(KNOWN_SUFFIXES):
                        continue
                    size = None if is_dir else entry.stat().st_size
                except OSError:
                    continue
                entries.append({"name": entry.name, "path": str(Path(entry.path)), "is_dir": is_dir, "size": size})
    except PermissionError:
        entries = []
    entries.sort(key=lambda e: (not e["is_dir"], e["name"].lower()))
    parent = str(path.parent) if any(path.parent == r or r in path.parent.parents for r in roots) else None
    return {"path": str(path), "parent": parent, "entries": entries}


def disk_free_gb(path: Path) -> float:
    path = Path(path)
    while not path.exists():
        path = path.parent
    return shutil.disk_usage(path).free / 1024**3


def safe_child(base: Path, relative: str) -> Path:
    """Join base + relative and make sure the result stays inside base."""
    base = base.resolve()
    target = (base / relative).resolve()
    if target != base and base not in target.parents:
        raise PermissionError("Path escapes the output directory.")
    return target
