#!/usr/bin/env python
"""Create a small, realistic paired-end FASTQ test set from an existing BAM.

Why: a full WES run takes hours and ~15-25 GB. For development and quick
demos we extract only the read pairs where BOTH mates aligned to one
chromosome (default chr22). Because the reads come from the real sample,
they have real error profiles, duplicates and adapters - unlike simulated
reads. The original FASTQ sequences/qualities are restored by
`samtools fastq` (reverse-strand reads are reverse-complemented back).

Usage:
    python scripts/make_test_subset.py IN.bam OUT_DIR PREFIX [--contig chr22] [--threads 4]

The BAM does not need an index (it is streamed). Nothing is overwritten:
OUT_DIR must not already contain the output files.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("prefix")
    ap.add_argument("--contig", default="chr22")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()

    if not re.fullmatch(r"[A-Za-z0-9_.-]+", a.contig) or not re.fullmatch(r"[A-Za-z0-9_.-]+", a.prefix):
        sys.exit("contig/prefix may only contain letters, digits, '_', '-', '.'")
    a.out_dir.mkdir(parents=True, exist_ok=True)
    r1 = a.out_dir / f"{a.prefix}_R1.fastq.gz"
    r2 = a.out_dir / f"{a.prefix}_R2.fastq.gz"
    for p in (r1, r2):
        if p.exists():
            sys.exit(f"Refusing to overwrite existing file: {p}")
    samtools = shutil.which("samtools") or sys.exit("samtools not found on PATH")

    expr = f'rname=="{a.contig}" && mrname=="{a.contig}"'
    view = [samtools, "view", "-u", "-f", "1", "-F", "0x900", "-e", expr, "-@", str(a.threads), str(a.bam)]
    collate = [samtools, "collate", "-u", "-O", "-@", str(a.threads), "-", str(a.out_dir / "collate_tmp")]
    fastq = [samtools, "fastq", "-n", "-@", str(a.threads), "-1", str(r1), "-2", str(r2),
             "-0", "/dev/null", "-s", "/dev/null", "-"]
    print("$", " ".join(view), "|", " ".join(collate), "|", " ".join(fastq), flush=True)
    p1 = subprocess.Popen(view, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(collate, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    p3 = subprocess.Popen(fastq, stdin=p2.stdout)
    p2.stdout.close()
    codes = [p3.wait(), p2.wait(), p1.wait()][::-1]
    print("exit codes:", codes)
    if any(codes):
        return 1
    for p in (r1, r2):
        print(f"{p}  {p.stat().st_size / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
