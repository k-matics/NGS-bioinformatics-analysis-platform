"""Live server resource snapshot for the dashboard."""
from __future__ import annotations

import os
import shutil

import psutil
from flask import current_app


def system_snapshot() -> dict:
    mem = psutil.virtual_memory()
    disks = []
    seen = set()
    for label, path in (("Root filesystem (/)", "/"), ("Results directory", current_app.config["RESULTS_DIR"])):
        try:
            dev = os.stat(path).st_dev
        except OSError:
            continue
        usage = shutil.disk_usage(path)
        disks.append({
            "label": label,
            "path": str(path),
            "same_device_as_previous": dev in seen,
            "total_gb": round(usage.total / 1024**3, 1),
            "free_gb": round(usage.free / 1024**3, 1),
            "used_percent": round(100 * (usage.total - usage.free) / usage.total, 1),
        })
        seen.add(dev)
    load1, load5, load15 = os.getloadavg()
    return {
        "cpu_percent": psutil.cpu_percent(interval=0.3),
        "cpu_count": psutil.cpu_count(),
        "load_avg": [round(load1, 2), round(load5, 2), round(load15, 2)],
        "memory": {
            "total_gb": round(mem.total / 1024**3, 1),
            "available_gb": round(mem.available / 1024**3, 1),
            "used_percent": mem.percent,
        },
        "disks": disks,
    }
