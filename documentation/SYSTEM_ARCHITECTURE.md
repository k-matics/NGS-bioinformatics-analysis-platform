# System Architecture

## Components

```
 Windows browser ──SSH tunnel──► gunicorn (Flask web app, 127.0.0.1:5000)
                                   │  reads/writes
                                   ▼
                              SQLite database (database/ngs_platform.db, WAL mode)
                                   ▲  reads/writes
 Flask ──enqueue run id──► Redis (127.0.0.1:6379) ──► Celery worker ──► workflow engine
                                                                          │ subprocess (no shell)
                                                                          ▼
                                                     samtools, bcftools, pigz, bwa, gatk …
                                                                          │
                                                                          ▼
                                                     results/<SAMPLE>/<RUN_ID>/…
```

* **Web app** never runs bioinformatics tools itself, so the UI stays responsive.
* **Worker** executes one workflow at a time by default (`NGS_WORKER_CONCURRENCY`), each with
  *threads_per_job* threads (Settings).
* **Live updates**: the run page opens a Server-Sent-Events stream (`/api/runs/<id>/events`) that pushes
  step status and new log lines once per second. The dashboard polls `/api/dashboard` and `/api/system`.

## Code layout

```
ngs_platform/
├── config.py                 defaults; every value overridable with NGS_* env vars
├── run.py                    web entry point (gunicorn run:app)
├── celery_worker.py          worker entry point (+ recovery of interrupted runs)
├── app/
│   ├── __init__.py           application factory
│   ├── extensions.py         SQLAlchemy + Celery objects
│   ├── models/               database tables (see below)
│   ├── engine/               workflow engine
│   │   ├── runner.py         CommandRunner: safe subprocess execution, pipes, cancel/pause
│   │   ├── step.py           Step base class, StepContext, StepResult
│   │   ├── executor.py       runs a WorkflowRun step by step; resume; provenance
│   │   ├── registry.py       workflow registry (WorkflowDefinition: steps, params, input roles)
│   │   └── params.py         ParamSpec: typed, range-checked workflow/caller parameters → form + validation
│   ├── pipelines/            one package per module
│   │   ├── intake/           Phase 2: sample intake validation
│   │   ├── fastq_qc/         FastQC, fastp steps
│   │   ├── alignment/        BWA-MEM+sort, duplicate marking, index, optional duplicate removal
│   │   ├── bam_processing/   analysis regions, input-BAM preparation, BAM QC (flagstat/stats/mosdepth)
│   │   ├── variant_calling/  VariantCallStep, VariantFilterStep (soft filters), MultiQC, run report
│   │   │   └── callers/      caller PLUG-INS: bcftools, gatk (HaplotypeCaller), freebayes, deepvariant (Docker)
│   │   ├── germline/         registers fastq_qc, fastq_to_bam, bam_to_vcf, fastq_to_vcf from the steps above
│   │   ├── comparison/       engine.py (normalise, allele matching, RTG vcfeval), steps.py (caller comparison,
│   │   │                     VCF inputs + build guard, VCF vs VCF), bam_steps.py (BAM vs BAM),
│   │   │                     evidence.py (BAM vs VCF read evidence), workflows.py (registers the 3 workflows)
│   │   ├── output_steps.py   register outputs as sample files; opt-in cleanup of intermediates
│   │   └── common.py         shared helpers (read groups, BED/regions, parsers, output validation)
│   ├── services/             reusable logic (reference checks, genome build, checksums, limitations.py =
│   │                         known-limitations warnings, report_service.py = report index + CSV exports, …)
│   ├── routes/               Flask blueprints (pages + JSON/SSE API)
│   ├── templates/, static/   Bootstrap 5 UI (served locally, no internet needed)
│   └── utils/                validators, logging
├── scripts/                  start_services.sh, stop_services.sh, status.sh, make_test_subset.py
├── tests/                    pytest suite
├── documentation/
├── results/  logs/  database/  uploads/  instance/ (secret key, pid files)
```

## Data model (provenance)

| Table | Holds |
|---|---|
| `samples` | Sample ID, sequencing type, layout |
| `sample_files` | Registered input files (path, role, origin internal/external, size, mtime, SHA-256, detected build, validation status). Files are referenced in place, never copied/modified. |
| `reference_genomes` | Registered references + last integrity report |
| `workflow_runs` | Unique run ID, workflow type+version, parameters, inputs, reference, status, times, output dir, warnings, errors, who submitted |
| `step_runs` | Per step: status, times, attempts, outputs, summary metrics, warnings, log path |
| `command_executions` | Every external command: argv, exit code(s), tool versions, times |
| `workflow_events` | Audit trail |
| `app_settings` | Settings changed from the UI |

Each run directory also gets a self-contained **`run_manifest.json`** with the same information.

## Run lifecycle

```
PENDING ─► QUEUED ─► RUNNING ⇄ PAUSED
                       │
          ┌────────────┼──────────────┐
          ▼            ▼              ▼
      COMPLETED      FAILED       CANCELLED
                       └── Resume ──┘  (completed steps re-validated, others re-run)
```

* **Duplicate prevention**: a fingerprint (workflow + inputs + parameters + reference) is stored; an
  identical run cannot be submitted while one is active. The worker also *claims* a run atomically
  (`UPDATE … WHERE status='queued'`), so a run can never execute twice.
* **Cancel**: SIGTERM to the whole process group (then SIGKILL after 15 s). Files are kept.
* **Pause**: SIGSTOP / SIGCONT to the process group (applies to external tools, not pure-Python checks).
* **Resume**: allowed for failed/cancelled runs. Inputs must be unchanged (size + mtime); completed
  steps are skipped only if their outputs still validate.
* **Worker crash / server restart**: on start, the worker marks runs left RUNNING as FAILED (resumable).
* **Disk pre-flight**: a run is refused if free space < workflow estimate + safety margin.

## Security

* Commands are argument lists executed without a shell; user input never becomes shell code.
* Sample IDs are restricted to `[A-Za-z0-9_.-]`; file paths must resolve (symlinks followed) inside the
  configured data roots (`NGS_DATA_ROOTS`, default `/data:~`).
* Downloads are restricted to the run's own output directory.
* Web app, Redis and worker bind to 127.0.0.1 only; access via SSH tunnel. No credentials are stored.
* No user authentication yet — acceptable while access is via SSH tunnel only; to be added before any
  wider (LAN) exposure.

## Switching to PostgreSQL later

Set `NGS_DATABASE_URL=postgresql+psycopg2://user@host/ngs` (and install `psycopg2`). Models use only
portable column types (JSON, DateTime, BigInteger). A migration tool (Alembic/Flask-Migrate) will be
added before the schema is used for production data.

## Workflow composition (Phase 3)

Workflows are lists of reusable steps. The four germline workflows share the same step classes, so a
change to e.g. duplicate marking applies everywhere, and the workflow version is bumped.

```
fastq_qc      : ref check → FASTQ integrity → FastQC raw → fastp → FastQC trimmed → MultiQC → report
fastq_to_bam  : … fastq_qc steps … → regions → BWA+sort → markdup → index → [dedup] → BAM QC → MultiQC → report
bam_to_vcf    : ref check → regions → BAM prepare → BAM QC → call_<caller> × N → MultiQC → report
fastq_to_vcf  : fastq_to_bam steps + call_<caller> × N
```

Steps pass results to later steps via the **database** (`StepRun.outputs` / `summary`), never via
memory, so a resumed run (new worker process) finds the outputs of steps that completed earlier.

### Variant-caller plug-ins
`app/pipelines/variant_calling/callers/<name>.py` — a class with `key`, `mode` (germline/somatic),
`params` (its own ParamSpec list), `requires_read_group`, `available()` and `build_commands(job)`.
`VariantCallStep` handles everything common: calling regions, running the commands, indexing,
`bcftools norm`, `bcftools stats`, validation. See DEVELOPER_GUIDE.md.

### GATK and Java
The conda environment's Java is 25; GATK 4.6 needs Java 17. GATK commands run with an environment
whose `JAVA_HOME`/`PATH` point to a Java 17-21 found under `/usr/lib/jvm` (or the *gatk_java_home*
setting). Nothing system-wide is changed.

## Phase 4-6 additions

* **Germline run (v1.2):** `… → call_<caller> → filter_<caller> (per caller) → compare_callers (≥2 callers) →
  multiqc → register_outputs → report → cleanup (opt-in)`.
* **Comparison workflows** (module *comparison* / *external*):

  | Workflow | Steps |
  |---|---|
  | `vcf_compare` | reference_check → regions → vcf_prepare → vcf_compare → report |
  | `bam_compare` | reference_check → regions → bam_metrics → bam_compare → report |
  | `bam_vcf_evidence` | reference_check → regions → bam_prepare → vcf_prepare → read_evidence → report |

  They reuse the germline building blocks: `AnalysisRegionsStep`, `prepare_bam()` and the build guard.
* **External tools** added:
  - **Docker** (DeepVariant image pinned in Settings; its image ID is recorded per run).
  - **RTG Tools** (vcfeval, with an SDF cache in `instance/rtg_sdf/`).
* **Pipeline outputs as inputs:** `SampleFile.origin = "pipeline"` marks files produced by the platform (checksummed,
  with the run ID in the label). Users cannot register files with this origin by hand.
* **Limitations:** `services/limitations.py` is the single list of known limitations. It is computed from the
  reference (e.g. ALT contigs without a `.alt` file) plus fixed items, and shown on the dashboard, the run page and
  in every HTML report.
