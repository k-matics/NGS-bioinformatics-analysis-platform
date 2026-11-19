# Installation Guide

Everything is installed in user space — no `sudo`, no system packages, existing conda environments untouched.

## 1. Conda environment

```bash
~/miniconda3/bin/conda env create -f environment.yml          # creates "ngs_platform"
```
(This is what was run on 2026-11-05; it uses ~2.5 GB in `~/miniconda3/envs/ngs_platform`
plus the conda package cache in `~/miniconda3/pkgs`.)

Installed versions (2026-11-05): Python 3.12.16, Flask 3.1.3, SQLAlchemy 2.1.8, Celery 5.6.3,
Redis 8.10.2, pysam 0.24.1, cyvcf2 0.31.1, samtools/bcftools/htslib 1.21, bwa 0.7.19, fastp 1.4.0,
FastQC 0.13.0, MultiQC 1.35, mosdepth 0.3.10, FreeBayes 1.3.9, RTG Tools 3.13, bedtools 2.31.1, pigz 2.8.

> **Note on Java:** the conda env brings OpenJDK 25 (needed by FastQC/RTG). GATK 4.6.2.0 is
> validated on Java 17 (`/usr/bin/java`). GATK commands are run with `JAVA_HOME` set to a Java 17-21 found under
> `/usr/lib/jvm` (override: Settings → *Java home for GATK*). Nothing system-wide is changed.

> **DeepVariant (Phase 4, installed 2026-11-09):** `docker pull google/deepvariant:1.10.0`, the CPU image, 4.9 GB
> on disk, image ID `sha256:96e40015…`. This needs no system change, because user `common` is already in the `docker`
> group. The image is set in Settings and pinned to a version; never use `latest`, so results stay reproducible.
>
> **RTG vcfeval reference cache:** on the first comparison the platform builds `instance/rtg_sdf/ref<id>_<name>/`
> (`rtg format`, about 1.2 GB for hg38). It is a derived copy, and the reference FASTA is not touched.

> **Note on BWA:** the hg38 BWA index was built with bwa 0.7.17; bwa 0.7.19 reads the same index format.

To reclaim disk space after installation (safe — only removes downloaded archives):
```bash
~/miniconda3/bin/conda clean --tarballs
```

## 2. Start / stop

```bash
cd <app_dir>
scripts/start_services.sh      # Redis + Celery worker + web app (all 127.0.0.1)
scripts/status.sh
scripts/stop_services.sh       # web + worker;  --all also stops Redis
```
Logs: `logs/app.log`, `logs/worker.log`, `logs/web.log`, `logs/redis.log`.

## 3. Open the UI from Windows

In PowerShell / Command Prompt on the Windows workstation:
```bash
ssh -L 5000:127.0.0.1:5000 <user>@<server-address>
```
Keep that window open and browse to **http://localhost:5000**.

## 4. Configuration (environment variables)

| Variable | Default | Purpose |
|---|---|---|
| `NGS_REFERENCE_FASTA` | `~/reference_genome/hg38.fa` | Default reference registered on first start |
| `NGS_DATA_ROOTS` | `/data:~` | Folders users may register files from |
| `NGS_RESULTS_DIR` | `ngs_platform/results` | Where run outputs go (move to a bigger disk when available) |
| `NGS_DATABASE_URL` | SQLite in `database/` | Switch to PostgreSQL later |
| `NGS_WORKER_CONCURRENCY` | `1` | Workflows executed in parallel |
| `NGS_PORT` / `NGS_REDIS_PORT` | 5000 / 6379 | Ports (localhost only) |

### Server-specific settings: `config/local.env`

Settings that differ per server go in `config/local.env` (not in git; template:
`config/local.env.example`). `config.py` reads it for every entry point — web app, worker,
scripts, `flask shell` — and real environment variables still win. Tests ignore it.

On the analysis server (since 2026-11-19):

```
NGS_RESULTS_DIR=/data/ngs_platform_results
```

### Moving the results folder

`scripts/move_results.sh <old> <new> [--delete-old]` (services stopped first): copies with rsync,
verifies every file by size + SHA-256 (manifest written into the new folder), re-points internal
symlinks, backs up the database and rewrites stored paths, sets `NGS_RESULTS_DIR` in
`config/local.env`, and deletes the old folder only with `--delete-old` and only after every check
passed. It refuses a non-empty target or running workflows. The 2026-11-19 move: 3,204 files
(34 GB) to `/data/ngs_platform_results`, 0 mismatches; manifest
`MIGRATION_MANIFEST_20261119-105700.tsv`, database backup
`database/ngs_platform.db.before-results-move-20261119-105700`.

## 4a. Building the DRAGEN-style masked reference (done 2026-11-18)

Needs ~12 GB disk (source FASTA 3.3 GB + masked FASTA 3.3 GB + BWA index 5.4 GB) and ~5 GB RAM for `bwa index`.
Nothing existing is modified; the script refuses to overwrite.

```bash
D=~/reference_genome/hg38_dragen_masked; mkdir -p $D/source; cd $D/source
B=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome
curl -fL -C - -O $B/GRCh38_full_analysis_set_plus_decoy_hla.fa      # 3.3 GB; -C - resumes
curl -fL -O $B/GRCh38_full_analysis_set_plus_decoy_hla.dict         # for MD5 verification
curl -fL -O https://raw.githubusercontent.com/Illumina/DRAGMAP/master/fasta_mask/hg38_alt_mask.bed
bash <app_dir>/scripts/build_masked_reference.sh \
     $D/source/GRCh38_full_analysis_set_plus_decoy_hla.fa $D/source/hg38_alt_mask.bed \
     $D/GRCh38_analysis_set_dragen_masked.fa
```
Then register it in **Settings → References** (name `hg38_dragen_masked`) and press *Check*. The check must
report "ALT masking verified".

## 5. Tests

```bash
cd <app_dir>
~/miniconda3/envs/ngs_platform/bin/python -m pytest -q tests
```
