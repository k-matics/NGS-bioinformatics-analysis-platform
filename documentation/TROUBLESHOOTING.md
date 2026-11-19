# Troubleshooting

## Platform

| Symptom | Cause / fix |
|---|---|
| Run stays **QUEUED** | Worker not running → `scripts/status.sh`, then `scripts/start_services.sh`. Check `logs/worker.log`. |
| "Could not reach the job queue (is Redis running?)" | Redis stopped → `scripts/start_services.sh`. |
| Run became **FAILED: Worker stopped while this run was executing** | Server/worker restarted mid-run. Click **Resume**. |
| "Not enough disk space" | The root disk is ~99 % full. Free space, or lower the safety margin in Settings (not recommended below 2 GB — a full root disk breaks the OS and database). |
| "An identical run is already active" | Duplicate protection: the same workflow, inputs and parameters are already queued/running. |
| "Input files changed since this run was created" on Resume | A file was modified/replaced after registration. Start a **new** run so results are tied to the new file. |
| Page doesn't load from Windows | SSH tunnel not open, or web app stopped. |

## Intake validation findings

| Message | Meaning and action |
|---|---|
| `pigz -t` fails (exit code 1/2) | The `.gz` file is truncated or corrupt — typically an incomplete transfer. Re-copy it and compare checksums with the provider. |
| "R1 and R2 are out of sync" | Read names differ at the same position: files from different samples/lanes, or one was truncated or re-sorted. Do not align. |
| "sequence length ≠ quality length" | Malformed FASTQ (corruption or bad conversion). |
| Quality encoding "Phred+64" | Very old Illumina data — must be converted before alignment. |
| BAM: "Different genome build: file is GRCh37" | The BAM was aligned to hg19. Do not compare with GRCh38 results: re-align the FASTQ (preferred) or lift over the VCF. |
| BAM: "Primary chromosomes match, but N contigs are not in the reference" | e.g. DRAGEN decoy/HLA contigs. Comparable on main chromosomes, but expect differences in decoy/HLA-related regions. |
| BAM: "no @RG read group" | GATK refuses it. Regenerate the BAM with read groups (Phase 3 pipeline does this). |
| BAM: "no index (.bai)" | `samtools index file.bam` (writes a new `.bai` next to it — needs write permission and a little space). |
| `samtools quickcheck` fails | BAM truncated (missing EOF block) or not a BAM. |
| Reference: "NOT ALT-aware" | The run used `hg38.fa`. Select reference `hg38_dragen_masked` when starting the run (BIOINFORMATICS_WORKFLOW.md §3). |
| Reference check error "ALT-masking sidecar present but the FASTA is not masked" | The FASTA next to a `.masking.json` was replaced or edited. Rebuild it with `scripts/build_masked_reference.sh` |

## Processing workflows (Phase 3)

| Message | Meaning and action |
|---|---|
| "Read accounting mismatch" | The BAM does not contain one primary record per input read. Never ignore: check disk space (truncated output) and the step log. |
| "Duplicate marking changed the number of primary records" | The marking tool dropped reads. Check the markdup log. |
| fastp "removed > 10 % of reads" | Library/sequencing quality problem. Inspect `failed_reads.fastq.gz` and the fastp report (reasons per read). |
| "Only N % of reads mapped" | Contamination, wrong species/reference, or adapter dimers. Check FastQC adapter/overrepresented sequences. |
| "Mean coverage … < 30x" / "… % bases reach 20x" | Under-sequenced or failed capture. ~4.6 % of SAMPLE_A target bases have no MAPQ ≥ 20 coverage because the reference is not ALT-aware (MHC, chr15-17, 19, 22): a known limitation, see BAM_PROCESSING.md §5. |
| GATK "No Java 17-21 found" | Settings → *Java home for GATK*. GATK 4.6 does not run on the conda Java 25. |
| GATK error about incompatible contigs / sequence dictionary | The BAM's @SQ lines differ from the reference (e.g. DRAGEN decoys). Use bcftools/FreeBayes, or re-align the FASTQ. |
| "Contig(s) not in the reference" | Contig names are case/style sensitive: `chr22` for hg38 UCSC, not `22`. |
| "… is not available" for a caller | Shown greyed out in the form with the reason (e.g. DeepVariant image not installed). |
| Ti/Tv outside 1.8-3.6 | Excess false positives (low-quality calls, contamination, wrong regions). Check the PASS Ti/Tv after filtering. |

## Filtering, DeepVariant, comparisons (Phases 4-5)

| Message | Meaning and action |
|---|---|
| "only N % of records PASS the filters" | Thresholds too strict for this data (e.g. low depth), or poor data. Check *Filtered by* in the run report. |
| DeepVariant "Docker image … not present" | `docker pull google/deepvariant:1.10.0`, or set the image in Settings. |
| DeepVariant permission denied writing output | The run folder must be writable by your user. The container runs as your user, not as root. |
| DeepVariant cannot be paused | Containers don't support pausing. Cancel works and removes the container. |
| "REF alleles do not match … another genome build" | The VCF is from another build or reference. It is refused; re-call it or lift it over (outside the platform, documented). |
| "has no sample column (sites-only VCF)" | Genotypes can't be compared. Request the sample VCF. |
| RTG "SDF" built on first comparison | One-off, ~3 min and ~1.2 GB in `instance/rtg_sdf/`. |
| vcfeval out of memory | Settings → RTG vcfeval memory. |
| "more than N variants to evaluate" (read evidence) | Safety limit, no silent sampling: restrict with a BED, use *PASS only*, or raise *Max variants*. |
| "Select 2-6 BAM files" / "at least 2 VCF files" | Tick the files to compare on the sample page. Pipeline outputs are registered automatically after a run. |
| Input files changed since this run was created | A registered file was modified after registration. Re-register it; results must stay tied to exact inputs. |

## "A high-quality baseline variant is not in our VCF at all"

Open the variant check, search the variant and click **Check reads**. If our reads carry the ALT
allele but all have MAPQ 0 with an `XA` hit on another contig (e.g. `chr11_KI270721v1_random`),
the reads are ambiguously placed on our reference and every caller ignores them. Lowering
filters or using raw tool defaults will not help. The difference is in the reference mask:
DRAGEN 4.3 uses Illumina's v2 mask (`hg38_v2_mask`), which also masks such duplicated unplaced
sequence; our `hg38_dragen_masked` uses the public DRAGMAP v1 mask (ALT contigs only).
