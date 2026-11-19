# NGS Platform

A modular web platform for processing and validating next-generation sequencing (NGS) data,
built for a clinical bioinformatics team. Runs are tracked end to end: every command, tool
version, parameter and file checksum is recorded.

**What it does**
- **Germline pipeline:** FASTQ QC and trimming (fastp), alignment (BWA-MEM), duplicate marking,
  BAM QC and target coverage, then variant calling with bcftools, GATK HaplotypeCaller
  (optionally in DRAGEN-GATK mode), FreeBayes and DeepVariant. Calls are normalised and
  soft-filtered: records are labelled, never removed.
- **Comparison and validation:** caller concordance (exact-allele matching and RTG vcfeval),
  BAM vs BAM, VCF vs VCF, and a variant-by-variant check of an external VCF (e.g. from DRAGEN)
  against ours. A viewer lets you search a variant, change thresholds and see the effect
  instantly, and check the reads behind any call.
- **Parameter presets:** built-in, user-defined, or derived from another pipeline's VCF header
  and filters. Each run records which preset it used and any changes made after applying it.
- **Safety:** reference build checks (GRCh37/GRCh38 never mixed silently), ALT-masked GRCh38
  reference support, no shell string building from user input, and identifiable output names.

**Stack:** Flask, Bootstrap, SQLite, Celery + Redis, server-sent events for live progress;
tools run from a conda environment and Docker (DeepVariant).

**Status:** research use only — not clinically validated. Known limitations are shown as
warnings in the UI and in every report.

## Quick start

```bash
conda env create -f environment.yml              # tools + Python packages (see INSTALLATION_GUIDE)
cp config/local.env.example config/local.env     # server-specific paths (results folder, GATK, data roots)
scripts/start_services.sh                        # Redis + Celery worker + web app on 127.0.0.1:5000
ssh -L 5000:127.0.0.1:5000 <user>@<server>       # from your PC, then open http://localhost:5000
```

Tests:

```bash
conda run -n ngs_platform python -m pytest -q tests
```

## Documentation

In `documentation/`: start with `USER_GUIDE.md` and `INSTALLATION_GUIDE.md`, then
`SYSTEM_ARCHITECTURE.md`, `BIOINFORMATICS_WORKFLOW.md`, `FASTQ_PROCESSING.md`,
`BAM_PROCESSING.md`, `VARIANT_CALLING.md`, `VARIANT_CALLER_COMPARISON.md`,
`BAM_VCF_VALIDATION.md`, `TROUBLESHOOTING.md` and `DEVELOPER_GUIDE.md`.

Sample names in the documentation (SAMPLE_A, SAMPLE_B, "Provider") are pseudonyms. No
sequencing data, results or databases are part of this repository.
