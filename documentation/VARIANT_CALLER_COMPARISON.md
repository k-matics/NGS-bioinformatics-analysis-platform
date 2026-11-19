# Variant Caller Comparison

This guide covers how the platform compares variant callers, and VCFs in general, and how to read the results.

> **Agreement is not accuracy.** If two callers agree, they may simply share the same error, for example a
> mis-mapped read in a duplicated gene that both of them trust. Only a **truth set** measures accuracy (§6).

## 1. Where comparisons happen

| Situation | How |
|---|---|
| Several callers on the **same BAM** | Select 2 or more callers in *FASTQ → BAM → VCF* or *BAM → VCF*. The run then contains a **Caller comparison** step |
| Any VCFs of one sample: in-house vs provider, old vs new pipeline, two runs | Use the **VCF vs VCF comparison** workflow (2–8 VCFs) |
| Accuracy against a known answer | Register a **Truth VCF** and **Truth confident-regions BED** for the sample. Both comparisons above then add precision and recall |

Every completed germline run **registers its BAM and filtered VCFs as sample files** (origin *Platform pipeline*,
with a SHA-256 checksum and the run ID in the label). They can be selected directly for comparisons.

## 2. Step 1 — filtering (per caller, soft)

The callers report *candidate* variants with very different levels of confidence. Before comparing them, each
caller's VCF is **soft-filtered**: a failing record gets a label in the FILTER column and **is never deleted**.

| FILTER | Rule | Callers |
|---|---|---|
| `NoAltGT` | Genotype is 0/0 or missing: no variant in this sample | all |
| `LowQual` | QUAL below threshold (default 20; off for GATK/DeepVariant, which already threshold) | bcftools, FreeBayes |
| `LowDP` | FORMAT/DP < 10 | all |
| `LowGQ` | Genotype quality < 20 | bcftools (`call -a GQ`), GATK, DeepVariant |
| `LowVAF` | ALT reads / (REF+ALT) < 0.15 | all |
| `HardFilterSNP` / `HardFilterIndel` | GATK Best-Practices hard filters (QD, FS, MQ, SOR, MQRankSum, ReadPosRankSum) | GATK |
| `StrandBias` | All ALT reads on one strand (SAF = 0 or SAR = 0) | FreeBayes |
| `RefCall`, `LowQual` | DeepVariant's own labels (kept) | DeepVariant |

**Why per caller?** QUAL has a different scale in every caller, and only some callers write GQ, SAF/SAR or QD.
A single shared filter would treat callers unfairly. All thresholds are on the start form (*Filtering* group of each
caller) and are stored with the run.

**Why soft?**
- Filtered calls remain visible, so you can still check them in IGV.
- The comparison can use *PASS only* (the default) or *all records*.
- The report shows how many records each rule labelled.

Output: `vcf/<caller>/<sample>.<caller>.filtered.vcf.gz`, with PASS statistics in `*.pass.stats.txt`.

## 3. Step 2 — making VCFs comparable (normalisation)

The same biological variant can be written in several ways. Before matching, every VCF goes through
`comparison/.../<set>.cmp.vcf.gz`:

```
bcftools view  [-s sample] [-T regions]      one sample, analysis regions only
bcftools norm  -a --atom-overlaps .         complex/MNP alleles → separate SNVs/indels
bcftools norm  -f hg38.fa -m -any           split multi-allelic sites, left-align indels
bcftools sort  + index
```

Examples:
- FreeBayes may report `ACG → TCA` as one complex allele. After *atomising* it becomes two SNVs (A→T, G→A), which
  is how GATK writes it.
- A deletion inside `AAAA` may be written at any of the four A's. *Left-alignment* moves it to the left-most one.

## 4. Step 3 — two matching methods

### 4a. Normalised allele matching (concordance matrix)

After normalisation, a variant is identified by **CHROM, POS, REF, ALT**.

| Result | Meaning |
|---|---|
| **Per set** | Variants, SNVs, indels, het, hom-alt, multi-allelic sites (counted in the raw VCF, before splitting), variants unique to this set, records not used (not PASS, no ALT in the genotype, …) |
| **Intersections** (UpSet-style) | For every exact combination of sets, how many variants are in exactly that combination. Example row: "GATK + DeepVariant, not bcftools, not FreeBayes" |
| **Pairwise** | Shared variants, A-only, B-only, Jaccard index (shared / union) |
| **Allele concordance** | Same CHROM/POS/REF/ALT |
| **Genotype concordance** | Of the shared variants, the fraction with the same genotype (0/1 vs 1/1). Genotypes are compared unphased: 0\|1 = 1\|0 = 0/1 |
| `variant_matrix.tsv.gz` | **One row per variant, one GT column per set** (`.` = not called). This is the concordance matrix |
| `discordant_genotypes.tsv` | Variants called by several sets with different genotypes |

### 4b. Haplotype-aware matching (RTG vcfeval)

Some differences remain after normalisation. Example: one caller writes a deletion plus a nearby SNV as one complex
event, and another writes them separately at slightly different positions.

**RTG vcfeval** applies each call set to the reference sequence, builds the two resulting haplotypes, and checks
whether they are the same DNA sequence. If they are, the calls match, whatever their representation.

- **Genotype mode:** alleles and genotype must match.
- **Allele mode** (`--squash-ploidy`): the allele must match; zygosity is ignored.
- **Representation-only** = matches found by vcfeval (allele mode) minus matches found by exact allele matching:
  the same variant, written differently.

**hap.py** (Illumina) is the other common benchmarking tool. It uses vcfeval as its matching engine and adds
GA4GH-style stratification. vcfeval was chosen because it is installed natively (RTG Tools 3.13) and hap.py
needs a legacy Python 2 build or Docker. Stratification (by region type) is a possible extension.

vcfeval needs the reference in its own binary format (an *SDF*). This is built once per reference in
`instance/rtg_sdf/` (≈ 1.2 GB for hg38, ~3 minutes). The reference itself is not touched.

## 5. Reading a caller comparison

1. **Start with the "found by all" row.** On real exome data, about 85–90 % of PASS variants are found by all
   modern callers. These are the easy, high-confidence variants.
2. **Look at the single-caller rows.** Variants found by only one caller are enriched for false positives of that
   caller, but they can also be real variants that the other callers missed.
3. **For a discordant variant, the next step is evidence, not a vote.**
   - Run the *BAM vs VCF read-level validation* workflow on it (BAM_VCF_VALIDATION.md).
   - Or open it in IGV.
   - Typical patterns: low MAPQ (duplicated region), strand bias, VAF around 0.1 (artefact or mosaic), a homopolymer
     indel.
4. **Genotype discordance** (0/1 vs 1/1) usually means low depth or allelic imbalance at the site.
5. **Multi-allelic sites** (e.g. 1/2) are counted from the raw VCF. They are split for matching.

## 6. Truth-set benchmarking (accuracy)

A **truth set** is a sample whose variants are known with high confidence, e.g. **GIAB HG001 (NA12878) / HG002**.
It is published as:
- a truth VCF;
- a **confident-regions BED**: only inside these regions is "not in the truth VCF" reliable evidence of "no variant".

How to use it:
1. Sequence the GIAB reference DNA in your own lab, with your own capture kit and sequencer.
2. Register the FASTQ, the GIAB truth VCF (*Truth VCF*) and the confident BED (*Truth confident-regions BED*).
3. Run *FASTQ → BAM → VCF* with the callers you want to evaluate.

For each caller, vcfeval reports:

| Metric | Formula | Question |
|---|---|---|
| **Precision** (PPV) | TP / (TP + FP) | Of what we called, how much is real? |
| **Sensitivity** (recall) | TP / (TP + FN) | Of what is real, how much did we find? |
| **F1** | 2·P·R / (P + R) | Both in one number |

These are reported for all variants, SNVs and indels/other. Indels are always harder.

This is the analytical-validation step a clinical lab needs. **No truth sample has been processed on this server
yet**, so no accuracy claim is made for this pipeline. The automated tests use a simulated genome with planted
variants, where sensitivity is 1.0 for every caller.

## 7. Output files (`comparison/callers/` or `comparison/vcf/`)

| File | Content |
|---|---|
| `comparison_report.html` | Everything above, readable |
| `comparison_summary.json` | Machine-readable summary |
| `per_set_counts.csv`, `intersections.csv`, `pairwise.csv` | Tables |
| `variant_matrix.tsv.gz` | Variant × set genotype matrix |
| `discordant_genotypes.tsv` | Same variant, different genotype |
| `vcfeval/<A>__vs__<B>.<mode>/` | Complete vcfeval output (tp, fp, fn VCFs and ROC curves). Open `fp.vcf.gz`/`fn.vcf.gz` in IGV |
| `vcfeval_pairwise.csv`, `truth_benchmark.csv` | vcfeval metrics |

## 8. Common problems

| Symptom | Cause / fix |
|---|---|
| "REF alleles do not match" in VCF vs VCF | The VCF comes from another build (e.g. GRCh37) or reference. It is refused: builds are never mixed. Re-call it from the FASTQ, or lift it over outside the platform and record that |
| "no sample column (sites-only VCF)" | Genotypes can't be compared. Ask the provider for the sample VCF |
| Multi-sample VCF | The sample with the platform sample ID is used, otherwise the first one. A warning is shown |
| Very low overlap with a provider VCF | Check that the regions match (provider VCFs are often whole-exome plus flanks, or filtered differently), whether PASS-only was used, and the sample identity |
| vcfeval: "Not enough memory" | Settings → RTG vcfeval memory |

## 9. Real-data examples (SAMPLE_A, 2026-11-09)

### 9a. Four callers on the same BAM, chr22 only (run WF-20261109-040542-336ce2, in-house BWA BAM, 5 min)

| Caller | Records | PASS | PASS Ti/Tv | Main filter labels |
|---|---|---|---|---|
| BCFtools | 1,487 | 1,355 (91 %) | 2.70 | LowDP 101, LowQual 48, LowVAF 28 |
| GATK HC | 1,468 | 1,342 (91 %) | 2.67 | HardFilterSNP 74, LowDP 49, LowGQ 26 |
| FreeBayes | 1,391 | 1,297 (93 %) | 2.71 | StrandBias 65, LowDP 50 |
| DeepVariant | 1,818 | 1,335 (73 %) | 2.66 | RefCall 384, LowVAF 246 |

Ti/Tv rose after filtering for every caller (e.g. BCFtools 2.49 → 2.70). The removed calls were enriched for errors.
DeepVariant's lower PASS fraction is expected: it also writes *RefCall* records, which are candidates it decided are
not variants.

* **1,189 of 1,472 distinct PASS variants (81 %) were found by all four callers.**
* The next-largest groups are variants found by 3 of 4 callers (24–38 each).
* 21–31 variants were found by only one caller.
* Pairwise genotype concordance on shared variants: **98.9–99.5 %**.
* vcfeval found 5–20 extra matches per pair that exact matching missed (*representation-only*). FreeBayes vs
  DeepVariant had the most (20), because FreeBayes writes complex alleles.

### 9b. Same callers, different BAM, chr22 (in-house BWA vs provider DRAGEN BAM), VCF vs VCF (WF-20261109-042715-b7ff86)

| Pair | Shared | Only first | Only second | Genotype concordance |
|---|---|---|---|---|
| GATK (our BAM) vs GATK (DRAGEN BAM) | 1,311 | 30 | **138** | 99.85 % |
| DeepVariant (our BAM) vs DeepVariant (DRAGEN BAM) | 1,306 | 29 | **141** | 99.77 % |
| GATK vs DeepVariant (same DRAGEN BAM) | 1,374 | 75 | 73 | 99.49 % |

**Interpretation:**
* With the same caller, the **BAM** changes about 140 calls on chr22, more than the difference between two callers
  on the same BAM.
* The extra calls are mostly on the DRAGEN side. Read evidence (BAM_VCF_VALIDATION.md §4b) shows they sit in regions
  where our alignment has MAPQ 0, which is the non-ALT-aware reference limitation.

**Note on a negative "representation-only" value (−5):** vcfeval matched slightly *fewer* variants than exact
matching for one pair. vcfeval evaluates whole haplotype regions and can count a cluster of nearby variants as one
mismatch when one record in it differs, or skip a region that is too complex. Small negative values are therefore
possible. The exact-match counts and vcfeval's own TP/FP/FN VCFs are both kept for inspection.

### 9c. Four callers on the same BAM, whole exome (run WF-20261109-043041-7bc18d, in-house BWA BAM)

All primary chromosomes, target BED + 100 bp, 10 threads. Total run time **1 h 57 min**; the run folder is 814 MB.

| Step | Time |
|---|---|
| BAM QC (mosdepth, stats) | 4 min |
| BCFtools | 14 min |
| GATK HaplotypeCaller | 42 min |
| FreeBayes | 17 min |
| DeepVariant 1.10 (Docker, WES model) | 38 min |
| Comparison incl. 12 vcfeval runs | 2 min |

| Caller | Records | PASS | Ti/Tv raw → PASS | Main filter labels |
|---|---|---|---|---|
| BCFtools | 74,003 | 66,181 (89 %) | 2.32 → 2.47 | LowDP 5,969, LowQual 3,034, LowVAF 1,452 |
| GATK HC | 72,740 | 66,036 (91 %) | 2.42 → 2.48 | LowDP 3,197, HardFilterSNP 2,940, LowGQ 1,094 |
| FreeBayes | 68,764 | 63,046 (92 %) | 2.44 → 2.52 | StrandBias 3,840, LowDP 3,204, LowQual 1,446 |
| DeepVariant | 91,203 | 64,245 (70 %) | 2.10 → 2.48 | NoAltGT 22,403, RefCall 21,417, LowVAF 13,056 |

A PASS exome Ti/Tv of about 2.5 is in the expected range for coding regions (≈ 2.8–3.3 for pure coding SNVs,
lower when flanks and UTRs are included). Filtering raised Ti/Tv for every caller, as on chr22.

* **58,358 of 71,346 distinct PASS variants (82 %) were found by all four callers**, the same proportion as chr22.
* Variants found by **only one caller** are mostly **indels**: GATK-only 1,387 (1,209 indels), BCFtools-only 1,269
  (882 indels), FreeBayes-only 1,045 (338 indels), DeepVariant-only 467 (245 indels). Indels in repeats and
  homopolymers are where callers disagree most. DeepVariant has the fewest single-caller calls.
* Pairwise genotype concordance on shared variants: **98.8–99.4 %**; Jaccard 0.89–0.90.
* vcfeval allele mode found up to 418 extra matches per pair that exact matching missed (FreeBayes vs DeepVariant;
  FreeBayes pairs again have the most because of its complex alleles). BCFtools vs GATK had −18 (see the note above).

These are **agreement figures, not accuracy**. The run (on `hg38.fa`) carries the non-ALT-aware reference warning; the ~5 % of
target bases without MAPQ ≥ 20 coverage (BAM_PROCESSING.md §5) are largely invisible to all four callers alike, so
"found by all four" does not mean those regions are covered.

### 9d. After the reference fix: agreement with the provider (Provider DRAGEN 4.3) VCF

Run **WF-20261119-001642-08530a**: same FASTQ, trimming, duplicate marking, callers and filters as before; only
the reference changed to `hg38_dragen_masked` (BIOINFORMATICS_WORKFLOW.md §3). Run time 2 h 37 min.

The provider VCF was restricted to the same target regions (xGen targets + 100 bp), PASS only, both sides
normalised with the same reference. (Their VCF also has 132,765 records outside the targets because DRAGEN's
small-variant caller was run without `--vc-target-bed`; those are not comparable and are excluded.)

| Our caller | Reference | Shared with DRAGEN | DRAGEN-only | Ours-only | % of DRAGEN calls found |
|---|---|---|---|---|---|
| GATK HC | `hg38.fa` | 63,310 | 12,485 | 2,726 | 83.5 % |
| GATK HC | **`hg38_dragen_masked`** | **66,746** | **9,049** | 2,846 | **88.0 %** |
| DeepVariant | `hg38.fa` | 63,519 | 12,276 | 726 | 83.8 % |
| DeepVariant | **`hg38_dragen_masked`** | **66,814** | **8,981** | 762 | **88.1 %** |

Provider in-target PASS calls: 75,795.

* **DRAGEN-only calls dropped by about 3,400 (−28 %).** On chr6 they fell from 1,931 to 725 (MHC region).
* Of the 9,049 still DRAGEN-only for GATK:
  * **4,841 are in our VCF but were removed by our filters** (LowDP 1,733, HardFilterSNP 1,506, LowDP+LowGQ 546, …).
    These differ because DRAGEN's own PASS threshold is much more permissive (QUAL ≥ 0.4, DP ≥ 2).
  * **4,208 were not called by GATK at all**; only 355 of them have DRAGEN QUAL ≥ 30. These are mostly
    low-confidence calls where callers legitimately differ.
* Within our run, the four callers still agree as before: 61,357 of 75,261 PASS variants (82 %) found by all four,
  genotype concordance 98.7–99.4 %. Every caller gained about 3,500 PASS variants.

**Remaining differences, and what would close them:** most of the gap is now filtering policy, not alignment.
A "DRAGEN-comparable" filter profile would show how much is threshold choice. Matching DRAGEN's counts exactly is
not the goal; a GIAB truth sample is needed to decide which pipeline is more accurate.
