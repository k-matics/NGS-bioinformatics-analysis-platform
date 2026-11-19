# FASTQ Processing

This guide covers what happens to raw reads before alignment, why each step exists, and how to read the results.
These steps are used by these workflows: **FASTQ QC & trimming**, **FASTQ → analysis-ready BAM** and
**FASTQ → BAM → VCF**.

```
R1/R2 .fastq.gz ─► integrity check ─► FastQC (raw) ─► fastp trimming ─► FastQC (trimmed) ─► alignment
                     (pigz -t, format,                 │
                      R1/R2 pairing)                   └─► failed_reads.fastq.gz  (removed reads are KEPT here)
```

## 1. What is in a FASTQ file

Each read takes up four lines:

```
@LH00263:26:23HCHGLT4:3:1391:4705:6997 1:N:0:ACGTACGT   ← read name + comment
ATATCGCTCGTAATGATG...                                  ← bases
+
IIIIIIII9IIII...                                       ← one quality character per base
```

* **Read name.** It tells you where the read came from: instrument `LH00263`, run `26`, **flow cell `23HCHGLT4`**,
  **lane `3`**, then the tile and the x/y position of the cluster. R1 and R2 of the same fragment share the name.
  That is how the pairing check knows they belong together.
* **Quality characters** encode Phred scores: Q = ASCII value − 33.
  * `I` = Q40, which means a 1 in 10,000 chance the base is wrong.
  * `9` = Q24, `#` = Q2.
  * NovaSeq X uses **binned** qualities: only a few distinct values such as 2, 12, 24 and 40. That is normal.

## 2. Integrity check (reused from intake)

Before any work, the platform:

* tests the entire gzip stream (`pigz -t`);
* checks the 4-line structure and quality encoding on the first 100,000 records;
* confirms R1 and R2 read names match record by record.

A truncated download is caught here, not three hours later.

## 3. FastQC: quality *description*

FastQC is run twice: on the raw reads and on the trimmed reads. Comparing the two shows what trimming fixed.
The report has a PASS/WARN/FAIL flag per module. **These flags were designed for whole-genome libraries, so read them
with care:**

| Module | Typical for good exome data | Worry when… |
|---|---|---|
| Per base sequence quality | High (Q30+), slight drop at the 3′ end | Median falls below Q20 over many cycles |
| Per tile sequence quality | Often **FAIL** on NovaSeq X: a few tiles have bubbles or edge effects | Many tiles affected, or quality is low overall |
| Per base sequence content | Often WARN/FAIL in the first ~10 bases (random-priming bias) | Imbalance along the whole read |
| Per sequence GC content | Exome GC (~50–60 %) is higher than genome GC (~41 %), so it often WARNs | Two peaks, which suggests contamination |
| Sequence duplication | **High is expected** for exomes (captured regions are re-sequenced deeply) | Above 50 % at moderate depth |
| Adapter content | FAIL on raw reads when inserts are shorter than the read length | Still present **after** trimming |

The platform turns each FAIL into a *warning*, never an error, and asks you to review the report.

## 4. fastp: trimming and filtering

### Why trim?

* **Adapters.** If the DNA fragment (the *insert*) is shorter than the read length (151 bp), the sequencer reads past
  the end of the fragment into the adapter. The aligner cannot match adapter sequence to the genome. It soft-clips
  it at best, or reports false mismatches and indels near read ends at worst.
  * SAMPLE_A has a median insert of ~160 bp, so about **half of the reads contained adapter**.
    fastp trimmed 835,581 of 1.64 M reads.
* **Poly-G tails.** On 2-colour instruments (NovaSeq, NextSeq), "no signal" is read as G. A cluster that has
  faded produces GGGG… tails. fastp detects NovaSeq/NextSeq from the read names and trims these automatically.
* **Low-quality reads.** fastp removes reads that have too many low-quality bases, too many N bases, or are too short
  after trimming.

### Parameters (Start-workflow form → Trimming)

| Parameter | Default | Meaning |
|---|---|---|
| Adapter/quality trimming | on | Turns fastp on or off. If off, alignment uses the raw reads |
| Adapter trimming | on | `--detect_adapter_for_pe`: finds adapters from the overlap between R1 and R2 |
| Minimum read length | 36 | Reads shorter than this after trimming are removed. Very short reads map ambiguously |
| Qualified base quality | 15 | A base is "low quality" below Q15 |
| Max % low-quality bases | 40 | Reads with more low-quality bases than this are removed |
| Max N bases | 5 | Reads with more N bases than this are removed |

fastp works on **pairs**: if one mate fails, both are removed. This keeps R1 and R2 in sync.

### Never silently discarding reads

Removed reads are written to `fastq/<sample>.failed_reads.fastq.gz`. The reason for each removal is in the
read comment and in `fastp.json → filtering_result`. The run report contains a **read accounting** table, and the run
**fails** if the numbers don't reconcile:

```
raw reads (R1+R2)            1,637,776
removed by fastp                 8,482   (0.52 %)   → failed_reads.fastq.gz
entering alignment           1,629,294
primary alignments in BAM    1,629,294   ✔ equal
after duplicate marking      1,629,294   ✔ equal (duplicates are flagged, not removed)
```

If more than 10 % of reads are removed, the run shows a warning. That usually indicates a library problem.

## 5. Outputs (inside `results/<sample>/<run_id>/`)

| File | What it is |
|---|---|
| `qc/fastqc_raw/*_fastqc.html` and `qc/fastqc_trimmed/*_fastqc.html` | FastQC reports |
| `fastq/<sample>.trimmed_R1/R2.fastq.gz` | Trimmed reads, used for alignment |
| `fastq/<sample>.failed_reads.fastq.gz` | Every read fastp removed |
| `fastq/<sample>.fastp.html` and `.json` | fastp report, before/after statistics |
| `reports/multiqc/multiqc_report.html` | All QC in one page |

## 6. Example: SAMPLE_A, chr22 subset

| Metric | Value | Interpretation |
|---|---|---|
| Q30 before → after | 94.9 % → 96.6 % | Excellent run quality |
| Reads removed | 0.52 % | Normal |
| Adapter-trimmed reads | 51 % | Short inserts (peak 113 bp, mean ~160 bp) |
| Duplication (fastp estimate) | 6.7 % | fastp estimates from sequence only; the real figure comes from duplicate marking after alignment (10.7 %) |
| FastQC FAIL | Per tile quality, Adapter content (raw only) | Adapter FAIL disappears after trimming, as it should. Per-tile FAIL is common on NovaSeq X |

## 7. Test subsets for development

`scripts/make_test_subset.py` extracts read pairs where **both mates mapped to one chromosome** from an existing BAM
and turns them back into FASTQ, restoring the original sequence and quality orientation:

```bash
python scripts/make_test_subset.py SAMPLE_A/sample_a.sorted.bam test_subsets/SAMPLE_A_chr22 SAMPLE_A_chr22 --contig chr22
```

Real reads keep real error profiles, adapters and duplicates, unlike simulated reads.

**Caveat:** read pairs whose mate mapped to another chromosome are not in the subset.

> **Correction (full-exome run WF-20261106-082225-5eeddc):** an earlier version of this guide blamed the ~6 % of chr22
> target bases with no coverage on the subset extraction. The **full exome shows the same gap on chr22 (5.9 %)**, so
> the cause is the **non-ALT-aware reference** (see BAM_PROCESSING.md §5), not the subset.
