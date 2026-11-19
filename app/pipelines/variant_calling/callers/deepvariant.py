"""DeepVariant (Google) — deep-learning germline caller, run from its official Docker image.

How it works: ``make_examples`` turns the reads around each candidate site into
a multi-channel "pileup image" (bases, base quality, mapping quality, strand,
...); a convolutional neural network (``call_variants``) classifies each image
as hom-ref / het / hom-alt; ``postprocess_variants`` writes the VCF.
The WES model is trained on exome data, the WGS model on genomes — use the one
that matches the library.

Container rules followed here:
* runs as the calling user (``--user uid:gid``) so output files are not root-owned;
* inputs (reference, BAM) are mounted READ-ONLY at their own paths;
* only the run's output folder is writable;
* a fixed container name lets a cancelled/failed run clean up its container.
Pause (SIGSTOP) is not supported for containers; cancel is.
"""
from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path

from ....engine.params import ParamSpec
from ....services.settings_service import get_setting
from ....services.tool_service import resolve
from .base import CallCommand, CallJob, VariantCaller, common_filter_params, register_caller, tool_defaults_param

DEFAULT_IMAGE = "google/deepvariant:1.10.0"
_cache: dict = {}


def image_name() -> str:
    return get_setting("deepvariant_image", DEFAULT_IMAGE) or DEFAULT_IMAGE


def image_id(image: str) -> str | None:
    """Local image ID (content digest) or None; cached 60 s (called on page renders)."""
    hit = _cache.get(image)
    if hit and time.time() - hit[1] < 60:
        return hit[0]
    docker = resolve("docker")
    iid = None
    if docker:
        try:
            r = subprocess.run([docker, "image", "inspect", "--format", "{{.Id}}", image],
                               capture_output=True, text=True, timeout=20)
            iid = r.stdout.strip() if r.returncode == 0 else None
        except (OSError, subprocess.TimeoutExpired):
            iid = None
    _cache[image] = (iid, time.time())
    return iid


def _mounts(read_only: list[Path], writable: list[Path]) -> list[str]:
    """-v arguments: each directory mounted at the same path inside the container."""
    rw = {p.resolve() for p in writable}
    ro = set()
    for p in read_only:
        d = p.resolve().parent if p.resolve().is_file() or not p.resolve().exists() else p.resolve()
        if not any(d == w or w in d.parents for w in rw):
            ro.add(d)
    args = []
    for d in sorted(rw):
        args += ["-v", f"{d}:{d}"]
    for d in sorted(ro):
        args += ["-v", f"{d}:{d}:ro"]
    return args


@register_caller
class DeepVariantCaller(VariantCaller):
    key = "deepvariant"
    name = "DeepVariant"
    tool = "docker"
    description = ("Deep-learning caller: reads around each candidate site are encoded as an image and "
                   "classified by a neural network. Top accuracy in public benchmarks. Runs in Docker.")
    params = [
        ParamSpec("model_type", "Model", "choice", "WES",
                  "Must match the library: WES for exomes/large panels, WGS for genomes.",
                  choices=[("WES", "WES (exome)"), ("WGS", "WGS (genome)")]),
        tool_defaults_param("make_examples MAPQ ≥ 5, base quality ≥ 10 (the model type above is still "
                            "used: run_deepvariant requires it)."),
        ParamSpec("min_MQ", "Min mapping quality", "int", 5,
                  "make_examples min_mapping_quality (DeepVariant default 5).", min=0, max=60),
        ParamSpec("min_BQ", "Min base quality", "int", 10,
                  "make_examples min_base_quality (DeepVariant default 10).", min=0, max=60),
    ] + common_filter_params(min_qual=0, min_depth=10, min_vaf=0.15, min_gq=20)

    def available(self):
        if not resolve("docker"):
            return False, "Docker not available."
        if not image_id(image_name()):
            return False, f"Docker image {image_name()} not present (docker pull {image_name()})."
        return True, ""

    def version_info(self):
        return {"image": image_name(), "image_id": image_id(image_name())}

    def build_commands(self, job: CallJob):
        p = job.params
        run_root = job.out_vcf.parent.parent.parent        # <run>/vcf/deepvariant/x.vcf.gz -> <run>
        ref = Path(job.reference)
        bam = Path(job.bam)
        inputs = [ref, Path(f"{ref}.fai"), bam, Path(f"{bam}.bai")]
        if bam.is_symlink():
            inputs.append(bam.resolve())
        bai = Path(f"{bam}.bai")
        if bai.is_symlink():
            inputs.append(bai.resolve())
        name = f"ngs-dv-{job.label or os.getpid()}".lower()
        argv = ["docker", "run", "--rm", "--init", "--name", name,
                "--user", f"{os.getuid()}:{os.getgid()}",
                *_mounts(inputs, [run_root]), image_name(),
                "/opt/deepvariant/bin/run_deepvariant",
                f"--model_type={p['model_type']}", f"--ref={ref}", f"--reads={bam}",
                f"--output_vcf={job.out_vcf}", f"--num_shards={max(1, job.threads)}",
                f"--intermediate_results_dir={job.work_dir}",
                f"--logging_dir={job.work_dir}/logs"]
        if not p.get("tool_defaults"):
            argv.append(f"--make_examples_extra_args=min_mapping_quality={p['min_MQ']},"
                        f"min_base_quality={p['min_BQ']}")
        if job.regions_bed:
            argv.append(f"--regions={job.regions_bed}")
        return [CallCommand(argv, cleanup=["docker", "rm", "-f", name])]
