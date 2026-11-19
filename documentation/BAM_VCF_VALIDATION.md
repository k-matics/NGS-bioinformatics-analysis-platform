# External BAM and VCF Validation

We often receive BAM and VCF files from external laboratories or providers, e.g. DRAGEN BAMs. This module answers
three questions:

| Question | Workflow |
|---|---|
| Are the two **alignments** comparable, and where do they differ? | **BAM vs BAM comparison** |
| Do the two **variant lists** agree? | **VCF vs VCF comparison** (see VARIANT_CALLER_COMPARISON.md) |
| Do the **reads** actually support a variant? | **BAM vs VCF read-level validation** |

All three workflows start from a sample's page, and their files must be registered for that sample. Mark provider files
with origin **External provider**: reports label each set as *In-house*, *External* or *Platform pipeline*.

> **A difference is a finding to explain, not automatically an error.** The aim is to find out whether a difference
> comes from alignment, preprocessing, variant calling or variant representation.

## 0. Genome-build guard (all three workflows)

Before anything is compared:

* **BAM:** `@SQ` names and lengths are checked against our reference.
  * A GRCh37 BAM against GRCh38 is refused.
  * Extra contigs, such as DRAGEN decoys and HLA, give a warning (*primary_compatible*).
* **VCF:** two independent checks.
  1. `##contig` lines are checked, when the file has them.
  2. The **REF alleles of the first 5,000 records are checked against our reference sequence.** A VCF from another
     build has mostly wrong REF bases, so if more than 1 % differ it is refused. This also catches VCFs that have
     no contig header lines.
* The original external file is **never modified.** Links, indexes and normalised copies are written in the run
  folder.

## 1. BAM vs BAM

Every BAM (2–6) gets **the same tools, settings and regions**. Comparing QC that was computed differently would be
meaningless.

| Tool | Metrics |
|---|---|
| `samtools flagstat` | Total, primary, **secondary**, **supplementary**, mapped, duplicates, properly paired, **singletons**, **mate on another chromosome** (all, and with MAPQ ≥ 5) |
| `samtools stats` (analysis contigs) | Insert size, error rate, read length, **MAPQ distribution** (0, 1–9, 10–19, 20–39, 40–59, 60+) and reads with MAPQ 0 |
| `samtools idxstats` | Reads per contig, grouped as primary, ALT, decoy, HLA, random/unplaced and unmapped |
| `mosdepth` | Mean target coverage and % of target bases ≥ N× (**thresholds configurable**, default 1,10,20,30,50,100), per target region |

Then **per-target differences** are listed. These are targets whose mean depth reaches the threshold (default 20×)
in one BAM but not in the other, in `targets_differently_covered.csv`, with counts per chromosome.

How to read the main differences:

| Difference | Typical cause |
|---|---|
| Different total/primary reads | Different trimming (reads removed), or a different FASTQ (check the checksums) |
| Supplementary count | Aligner behaviour with chimeric or split reads (`-Y` soft- vs hard-clipping) |
| Duplicate % | Duplicate-marking tool and optical-distance setting |
| **MAPQ 0 fraction** | **Reference version and ALT handling.** A non-ALT-aware alignment gives MAPQ 0 in ALT-duplicated regions |
| Reads on decoy/HLA contigs | Only references that contain these contigs, e.g. DRAGEN's, can place reads there. Otherwise those reads are forced onto similar primary sequence |
| Targets covered in only one BAM, clustered on chr6 (MHC), 17, 19, 22 | The ALT-contig effect (BAM_PROCESSING.md §5) |

The original outputs (flagstat JSON, stats, idxstats, mosdepth) are kept per BAM in `qc/compare/<label>/`.

## 2. VCF vs VCF

This uses the same engine as the caller comparison. See VARIANT_CALLER_COMPARISON.md §3–§5.

Points specific to external VCFs:
* **Restrict both VCFs to the same regions.** The target BED padded by 100 bp is used by default. Provider VCFs
  often extend further, which would otherwise create apparent "external-only" calls.
* **Filtering:** providers apply their own filters. *PASS only* compares what each side reports as final. Untick it
  to compare all records.
* **Sample names** may differ (e.g. `SAMPLE_A_provider` vs `SAMPLE_A`). Matching is by position and allele, not by name.

## 3. BAM vs VCF read-level validation

For **every variant of a VCF**, the workflow measures the reads of a BAM at that position. Duplicates, QC-fail,
secondary and supplementary reads are excluded, as the callers do.

| Column | Meaning | Why it matters |
|---|---|---|
| `depth`, `ref_reads`, `alt_reads`, `other_reads` | Reads showing REF, ALT, or something else (a third base or a different indel) | Is there any ALT support at all? |
| `vaf` | ALT / (REF + ALT) | Germline expectation: about 0.5 (het) or about 1.0 (hom). Around 0.1 suggests an artefact, contamination or mosaicism |
| `alt_fwd` / `alt_rev`, `strand_bias_p` | ALT reads per strand, with a Fisher exact test against the REF strand ratio | Real variants appear on both strands. Strand-specific errors are a classic artefact |
| `mean_mq_ref` / `mean_mq_alt`, `mq0_frac` | Mapping quality | Low ALT MAPQ or many MAPQ-0 reads suggest a paralogous or ALT region, where the reads may belong elsewhere |
| `mean_bq_alt` | Base quality of ALT bases (SNVs) | Low values suggest sequencing errors |
| `alt_softclip_frac` | ALT reads that are soft-clipped | Local misalignment, a nearby indel, or a structural variant |
| `alt_median_dist_to_end` | Position of the ALT within the read | Errors cluster at read ends, e.g. leftover adapter or end-repair artefacts |
| `vcf_ad`, `gt`, `filter` | What the VCF itself claims | Compare the claim with the observation |

**Evidence notes** are observations that need review. They are never automatic verdicts.

`no_reads`, `no_alt`, `low_alt` (fewer than 3 ALT reads), `vaf_het`, `vaf_hom`, `strand`, `low_mq`, `mq0`,
`softclip`, `read_end`, `third_allele`.

A note can have an innocent explanation. A real heterozygous variant at the end of a capture target can show strand
imbalance, and a mosaic variant has a low VAF.

**IGV-compatible outputs:**
- `igv_flagged_sites.bed`: load it as a track.
- `igv_batch.txt`: *Tools → Run Batch Script*. It loads the reference, BAM and VCF, jumps to each flagged variant
  (up to 500) and saves a snapshot. Paths are server paths, so adjust them if IGV runs on Windows with a mapped drive.

**Typical uses:**
* Does our BAM contain the reads for the provider's variants? Use **external VCF + our BAM**. A provider-only variant
  with ALT reads in our BAM was a calling or filtering difference. One without ALT reads is an alignment difference
  (MAPQ, ALT region) or a provider artefact.
* Does the provider's BAM support our calls? Use **our VCF + the external BAM**.
* QC of one VCF: **our VCF + our BAM** shows the VAF distribution and flags artefact-like calls.

Limits: complex/MNP alleles are decomposed before evaluation. Very long indels and structural variants are not
evaluated reliably from the pileup, so use IGV for those.

## 4. Real-data examples (SAMPLE_A, 2026-11-09)

### 4a. BAM vs BAM: in-house BWA vs provider DRAGEN 4.3, all primary chromosomes (WF-20261109-040542-2cbd80, 7 min)

| Metric | In-house BWA | DRAGEN | Explanation |
|---|---|---|---|
| Primary reads | 69,884,036 | 69,785,495 | Trimming differs (DRAGEN: quality 30, min length 50) |
| Average read length | 132 | 148 | Our fastp adapter trimming shortens reads more; DRAGEN trims less aggressively |
| Duplicates | 10.86 % | 11.32 % | Different duplicate-marking implementations |
| Singletons | 1,481 | 384,480 | DRAGEN leaves more mates unmapped (255,232 unmapped records vs 126) |
| Mate on another chromosome | 463,372 | 221,690 | Without decoys and ALT handling, more pairs are split across contigs in our alignment |
| **Reads on ALT contigs** | **1,709,728** | 10,856 | **Non-ALT-aware alignment.** DRAGEN masks the ALT copies |
| Reads on decoy / HLA contigs | 0 / 0 | 16,436 / 687 | Our reference has no decoys or HLA |
| MAPQ 0 (% reads) | 5.70 % | 4.43 % | Same cause |
| Mean target coverage | 139.6× | 143.0× | |
| **% target bases ≥ 20×** | **94.61 %** | **97.33 %** | |

**Targets reaching 20× in only one BAM:**
- 6,052 only in DRAGEN and 600 only in ours, out of 197,676 compared.
- Top chromosomes: **chr6 1,461** (MHC), chr17 805, chr19 478, chr16 427, chr15 414.
- This is exactly the ALT-contig pattern. The tool found it without being told where to look.

### 4b. BAM vs VCF: DRAGEN-BAM GATK calls checked in our BAM (WF-20261109-042715-84b482, chr22)

The 1,601 GATK calls made on the DRAGEN BAM were checked against the reads in our BAM:

| Note | Variants |
|---|---|
| >20 % MAPQ-0 reads at the site (`mq0`) | 168 |
| Low MAPQ of ALT reads (`low_mq`) | 131 |
| No ALT read in our BAM (`no_alt`) | 81 (95 incl. sites with no usable reads at all) |
| VAF unusual for het / low for hom | 63 / 36 |
| Third allele / other indel | 35 |

Median observed VAF: **0.455 for 0/1 calls, 1.0 for 1/1 calls**, as expected for germline variants.

The provider-only calls are largely explained by **mapping quality in our alignment**: the reads are there but sit at
MAPQ 0, or they were placed on an ALT contig. This is an alignment difference, not a calling difference. Open
`igv_flagged_sites.bed` together with both BAMs in IGV to see it.

### 4c. Variant-by-variant: is every Provider variant in our VCF? (2026-11-19)

Workflow **VCF variant-by-variant check**. Baseline: the Provider DRAGEN 4.3 VCF (PASS, inside xGen targets + 100 bp,
75,795 variants). Ours: run WF-20261119-001642-08530a on the masked reference, **all records** of our filtered VCF.

| Status | GATK (WF-20261119-043156-62ad76) | DeepVariant (WF-20261119-043156-26c341) |
|---|---|---|
| Found, PASS | 66,746 (88.1 %) | 66,814 (88.2 %) |
| Found, but filtered by our filters | 4,841 | 4,078 |
| Site seen, genotyped as reference (DeepVariant RefCall) | — | 1,060 |
| Written differently (vcfeval haplotype match) | 16 | 104 |
| Different variant within 10 bp | 487 | 469 |
| Not in our VCF | 3,705 | 3,270 |
| **Present in our VCF in any form** | **94.5 %** | **95.1 %** |

What the numbers say:
* **Our filters account for 4,841 of the 9,049 missing GATK calls** (LowDP 2,962, HardFilterSNP 2,294, LowGQ 1,063;
  a call can carry several labels). Ignoring our FILTER column, 94.5 % of the Provider variants are present.
  Requiring only depth ≥ 2 (DRAGEN's own rule) gives the same 94.4 %; depth ≥ 10 gives 90.5 %.
* **The missing ones are mostly low-confidence in DRAGEN too**: of the 4,192 not called / nearby, 3,840 have DRAGEN
  QUAL < 30 and only 352 have QUAL ≥ 30.
* Restricting the baseline to DRAGEN QUAL ≥ 30 (viewer: *baseline QUAL ≥ 30*): 58,066 variants, of which
  **54,935 (94.6 %) are found PASS** by our GATK and only **276 are not called at all**. These 276 are the list to
  review in IGV (CSV from the viewer).
* Most remaining differences are therefore a *filtering policy* question (how much low-quality evidence to report),
  not missing data. A truth sample is needed to decide which policy is more accurate.
