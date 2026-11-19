# Developer Guide

## Adding a workflow step

```python
# app/pipelines/<module>/steps.py
from ...engine.step import Step, StepResult, StepFailure

class SamtoolsFlagstat(Step):
    key = "flagstat"                      # stable identifier (stored in DB)
    name = "samtools flagstat"            # shown in the UI
    tools = ("samtools",)

    def applicable(self, ctx):            # optional: return (False, reason) to SKIP
        return bool(ctx.files("bam")), "No BAM registered."

    def execute(self, ctx, runner):
        bam = ctx.files("bam")[0].path
        out = ctx.step_dir("qc") / "flagstat.txt"          # always write inside the run dir
        runner.run(["samtools", "flagstat", "-@", ctx.threads, bam], stdout_path=out)
        return StepResult(outputs={"flagstat": str(out)}, summary={...}, warnings=[...])

    def validate(self, ctx, result):      # optional: stronger than "file exists and non-empty"
        super().validate(ctx, result)
```

Rules:
* **Never** build shell strings. `runner.run([...])` for one command, `runner.run([[cmd1...], [cmd2...]])`
  for a pipe. Use `stdout_path=` instead of `>`.
* Raise `StepFailure("clear, actionable message")` for scientific/technical failures.
* Long pure-Python loops must call `ctx.check_cancel()` regularly.
* `validate()` must check outputs are *technically valid* (e.g. `samtools quickcheck`), not just present.
  Resume relies on it.

## Registering a workflow

```python
# app/pipelines/<module>/__init__.py
register(WorkflowDefinition(key="...", name="...", version="1.0", description="...",
    build_steps=lambda params: [StepA(), StepB()], default_params={...},
    check_inputs=lambda sample, files: [...problems...],
    estimate_disk_gb=lambda sample, files: ...))
```
and import the package in `app/engine/registry.py::_load_pipelines`.
Bump `version` whenever the scientific behaviour of a workflow changes (it is recorded per run).

## Adding a variant caller (plug-in)

```python
# app/pipelines/variant_calling/callers/mycaller.py
from ....engine.params import ParamSpec
from .base import CallCommand, CallJob, VariantCaller, register_caller

@register_caller
class MyCaller(VariantCaller):
    key, name, tool, mode = "mycaller", "My Caller", "mycaller", "germline"
    requires_read_group = False
    description = "One-line explanation shown in the UI."
    params = [ParamSpec("min_MQ", "Min mapping quality", "int", 20, "help text", min=0, max=60)]

    def available(self):                       # (False, reason) greys it out in the form
        return True, ""

    def build_commands(self, job: CallJob):    # must produce bgzipped job.out_vcf
        argv = ["mycaller", "--ref", job.reference, "--bam", job.bam, "--min-mq", str(job.params["min_MQ"])]
        if job.regions_bed:
            argv += ["--bed", job.regions_bed]
        return [CallCommand([argv, ["bgzip", "-c"]], stdout_path=job.out_vcf)]
```
Then import the module in `callers/__init__.py`. Indexing, normalisation, stats, validation, the
form fields (`c.<workflow>.<caller>.<param>`) and the report entry come for free.

**Filtering:** add `common_filter_params(min_qual=…, min_gq=None if the caller writes no GQ)` to `params`.
For caller-specific rules, override `filter_rules(params)`, which returns `FilterRule(name, bcftools_expr,
description)`. Only build expressions from validated numeric parameters, never from free text.
**Containers:** return `CallCommand(argv, cleanup=["docker", "rm", "-f", name])`. The cleanup runs if the command
fails or is cancelled. `version_info()` adds provenance (e.g. the image digest) to the step summary. Add an end-to-end
check to `tests/test_germline.py::_check_planted`.

## Workflow parameters

Declare `ParamSpec`s in `param_specs`; the form, type/range validation, defaults and storage in
`run.params` are automatic. `text` parameters only allow `[A-Za-z0-9_.,:+-]`. Semantic checks that
need the reference (e.g. "contig exists") go in `validate_params(params, reference)`.

## Parameter presets

`ParamPreset` (`app/models/preset.py`) stores a possibly partial params dict in the same shape as `run.params`
(top-level keys + `caller_params`). `services/preset_service.py`:
* `validate_preset_params()` checks every value against the union of all workflows' `form_specs()` and the caller
  plug-ins' `params` — the same `ParamSpec`s as the form, so a preset cannot carry a value the form would refuse;
* `ensure_builtin_presets()` creates the locked *In-house default* (called at start-up and on the Presets page);
* `record_on_params()` writes `params["preset"] = {id, name, changed_after_applying}` when a run starts.
Adding a new parameter to a workflow or caller needs no preset change; old presets simply don't set it.

`services/filter_profile.py`: `extract_profile(vcf)` (header + record scan, no sampling) and
`profile_to_params(profile)` (explicit mapping + notes). To support another caller's header format, extend
`RULE_RE` / `_caller_from_header` and add a test with a header snippet.

## Variant-by-variant check

`pipelines/comparison/variant_check.py` (`VariantCheckStep`, workflow `vcf_variant_check` registered at the end of
`comparison/workflows.py`). It reuses `VcfInputsStep` (build guard + normalisation) and `engine.run_vcfeval`; it does
not modify `vcf_compare`. Viewer: `routes/workflows.py` (`variant_check`, `variant_check_csv`), filtering in
`_vc_filtered`, counting in `variant_check.what_if()`.

## Running tests

```bash
~/miniconda3/envs/ngs_platform/bin/python -m pytest -q tests
```
Tests use temporary databases/directories and synthetic files (a reference *stub* with real GRCh37/38
chromosome lengths, tiny FASTQ and BAM files). They never touch real data or `results/`.

| Test file | Covers |
|---|---|
| `test_validators.py` | sample ID rules, path-root confinement incl. symlink escape, role detection |
| `test_genome_build.py` | build detection, GRCh37-vs-38 incompatibility, decoy handling, reference integrity |
| `test_runner.py` | no-shell execution, pipefail semantics, exit codes, cancel kills process group, pause/resume |
| `test_workflow.py` | full intake runs, truncated gzip, unpaired R1/R2, malformed FASTQ, resume, changed inputs, duplicate prevention, crash recovery, disk pre-flight |
| `test_routes.py` | pages, sample/file/run flow via HTTP, SSE stream, path traversal |
| `test_comparison.py` | Fisher test against exact values, allele/genotype concordance arithmetic, caller-specific filter rules, read-only Docker mounts; end to end on the mini-genome: soft filtering keeps every record, planted variants PASS, caller comparison + vcfeval, **truth benchmark sensitivity 1.0**, opt-in cleanup, output registration, VCF vs VCF, BAM vs BAM (50 % subsample halves reads), read evidence (VAF ≈ 1.0 hom / ≈ 0.5 het), other-build VCF refused, **DeepVariant in Docker calls the planted genotypes** |
| `test_germline.py` | **simulated mini-genome with planted variants**: FASTQ→VCF end to end with bcftools + FreeBayes (+ GATK HC via BAM→VCF), genotype correctness, read accounting, duplicates flagged not removed, regions, build refusal, RG copy, parameter validation, form posting |
