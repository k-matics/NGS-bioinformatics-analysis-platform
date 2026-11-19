# Variant Calling (germline)

This guide covers how variants are called from the analysis-ready BAM, how the caller plug-ins work, and how to read
the output.

> **Scope:** the calling covers germline short variants (SNVs and small indels) in a single sample. Each caller's
> VCF is kept raw (as the caller wrote it), normalised, and **soft-filtered** with caller-specific technical filters
> (§8). The output is intended for research and pipeline development. A VCF is not a validated clinical result.

## 1. Germline vs somatic

* **Germline** variants are inherited and present in every cell. For a diploid genome the expected allele fractions
  are about **50 % (heterozygous, 0/1)** or **100 % (homozygous alternate, 1/1)**.
* **Somatic** variants (tumours) can occur at any fraction, e.g. 3 %. They need different statistical models and
  tumour/normal inputs (e.g. Mutect2).
* The platform keeps these apart. All callers here declare `mode = "germline"`. Somatic callers will be separate
  plug-ins and separate workflows.

## 2. Callers available

| Caller | How it works | Strengths | Weaknesses |
|---|---|---|---|
| **BCFtools** mpileup/call (`-m`) | Looks at the pile of bases at each position and computes genotype likelihoods | Fast, transparent, low memory | Indels in repeats; no local re-assembly |
| **GATK HaplotypeCaller** | Finds *active regions*, **re-assembles** the reads into candidate haplotypes, then genotypes | GATK Best Practices; strong for indels | Slower; needs read groups and Java 17 |
| **FreeBayes** | Bayesian, haplotype-based over short windows | Sensitive; reports complex alleles (MNPs) | Single-threaded here; more low-quality calls without filtering |
| **DeepVariant** 1.10 (Docker) | Encodes the reads around each candidate site as a multi-channel *image* (bases, qualities, MAPQ, strand) and classifies it with a convolutional neural network (hom-ref / het / hom-alt) | Top accuracy in public benchmarks (GIAB/precisionFDA); well calibrated GQ | Needs the model that matches the data (**WES** vs WGS); runs in a container (pause is not supported, cancel is); slower on CPU |

Select one or more callers in the start-workflow form. **Each caller gets its own VCF.**

## 3. Caller-specific parameters

Parameters are **not shared** between callers. "Minimum base quality" is implemented differently in each tool, so
each caller has its own settings. The values you choose are stored with the run (`run.params.caller_params`) and
printed in the report.

| Caller | Parameter | Default | Option |
|---|---|---|---|
| BCFtools | Min mapping quality | 20 | `mpileup -q` |
| | Min base quality | 20 | `mpileup -Q` |
| | Max depth | 1000 | `mpileup -d`. Reads beyond this per position are ignored; raise it for very deep panels |
| | Variant sites only | on | `call -v` |
| GATK HC | Min call confidence | 30 | `--standard-min-confidence-threshold-for-calling` |
| | Min mapping / base quality | 20 / 10 | |
| | PCR indel model | CONSERVATIVE | Use **NONE for PCR-free** libraries. Exomes are PCR-amplified |
| | Java heap | 8 GB | |
| FreeBayes | Min mapping / base quality | 20 / 20 | |
| | Min alternate fraction | 0.2 | An allele must be in ≥ 20 % of reads to be considered |
| | Min alternate count | 2 | |
| | Min coverage | 5 | |
| DeepVariant | Model | WES | `--model_type`: must match the library |
| | Min mapping / base quality | 5 / 10 | DeepVariant defaults (`make_examples`) |

Every caller also has its own **Filtering** parameters (§8). How DeepVariant runs:
- `docker run --user <you> -v <reference dir>:ro -v <BAM dir>:ro -v <run folder> google/deepvariant:1.10.0 run_deepvariant …`
- Inputs are mounted read-only and only the run folder is writable.
- The **image ID** (a content digest) is recorded with the run, so the exact model version is traceable.

**Mapping quality (MAPQ)** is the aligner's confidence that a read is placed correctly.

* MAPQ 60 means unique.
* **MAPQ 0** means the read fits equally well in several places, as in segmental duplications. It also happens in
  **ALT-contig regions with our non-ALT-aware reference** (see BIOINFORMATICS_WORKFLOW.md §3).
* Reads below the minimum MAPQ are ignored. This is why genes such as *PMS2*, *SMN1/2* and *CYP2D6* need special
  handling in clinical practice.

### GATK: DRAGEN-GATK mode

`dragen_mode` adds `--dragen-mode`: DRAGEN's genotyping model inside GATK (base-quality dropout,
foreign-read detection, DRAGEN MAPQ transform, posterior-based QUAL, soft-clipping of low-quality
read ends, 150 bp indel padding). Values set on the form (min MAPQ, min call confidence, base
quality) are passed after it and win. With `dragstr_calibrate` (default on) the STR indel-error
model is fitted to the BAM:

1. `ComposeSTRTableFile` — once per reference, cached in `instance/dragstr/` (outside the
   reference folder);
2. `CalibrateDragstrModel` on the run's BAM and calling regions → `dragstr_params.txt`;
3. `HaplotypeCaller --dragen-mode --dragstr-params-path …` (the PCR indel model is then not passed:
   DragSTR replaces it).

Why: on SAMPLE_A, 57 % of the Provider variants that plain GATK missed on DRAGEN's own BAM sit on
reads with MAPQ < 20, which GATK ignores by default and DRAGEN uses. DRAGEN mode + MAPQ ≥ 0 +
QUAL ≥ 1 recovered 706 of 2,619 missed sites near those positions with few extra calls. The preset
"DRAGEN-GATK, MAPQ 0 (Provider-like)" sets this up. Calls from low-MAPQ reads have uncertain
placement — review them before clinical use.

### Raw mode (tool defaults)

Each caller has a `tool_defaults` switch ("Raw command"). When on, only the inputs/outputs,
calling regions and speed-only options are passed; every other option of that caller is left to
the tool's own default. Filters are separate and still apply (the built-in preset
"Raw tool defaults (no filters)" also turns them all off). Use it to see what a tool does out of
the box, or to reproduce a command from a paper or another lab.

## 4. Calling regions

All callers are restricted to `regions/calling_regions.bed`: the target BED on the selected contigs, padded by 100 bp.

* Calling off-target would add slow, low-coverage, unreliable calls.
* If no BED is registered, whole chromosomes are called, with a warning for WES samples.

## 5. Normalisation

The same variant can be written in several valid ways. For example, deleting one `A` from `AAAA` can be written at
any of the four positions. Callers also differ in how they write multi-allelic sites.

```
bcftools norm -f hg38.fa -m -any --check-ref w
```

* **Left-aligns** indels to the left-most equivalent position and makes them parsimonious.
* **Splits multi-allelic records** (`-m -any`): a site with genotype `1/2` becomes two records.
* **Checks the REF allele** against the reference. Mismatches are reported, which is a strong sign of a build or
  reference mismatch.

Both files are kept for every caller:

* `*.raw.vcf.gz`: exactly what the caller produced, for provenance.
* `*.norm.vcf.gz`: canonical, ready for comparison (Phase 4).

FreeBayes reports **MNPs and complex alleles**, e.g. two adjacent SNVs as one record. Normalisation does not split
these. Comparison tools such as RTG vcfeval and hap.py handle different representations, which is why Phase 4 uses
them instead of naive position matching.

## 6. Output files (`vcf/<caller>/`)

| File | Content |
|---|---|
| `<sample>.<caller>.raw.vcf.gz` + `.tbi` | Caller output (bgzip + tabix index) |
| `<sample>.<caller>.norm.vcf.gz` + `.tbi` | Normalised |
| `<sample>.<caller>.stats.txt` | `bcftools stats` (also shown in MultiQC) |
| `<sample>.<caller>.summary.json` | Counts shown in the report |

**File names are identifiable.** When a run is started from a parameter preset, the preset is part
of every caller file name, so files from different parameter sets cannot be confused after
download — e.g. `SAMPLE_A.provider-like-dragen-4-3.gatk_hc.filtered.vcf.gz` instead of
`SAMPLE_A.gatk_hc.filtered.vcf.gz`. If you changed fields after applying the preset, `-modified` is
added (`SAMPLE_A.provider-like-dragen-4-3-modified.gatk_hc...`), so a name never claims values that
were not used. Without a preset the name is `<sample>.<caller>...` as before. The VCF's *sample
column* is always the plain sample ID. The registered file label also shows the preset.

## 7. Sanity metrics

These metrics catch gross problems. They **do not** measure accuracy; that needs a truth set such as GIAB.

| Metric | Expected (germline) | If off… |
|---|---|---|
| **Ti/Tv** (transitions A↔G, C↔T / transversions) | about 2.0–2.1 genome-wide, **about 2.6–3.3 in exome targets** | Low Ti/Tv means excess false positives, since random errors give Ti/Tv ≈ 0.5. A warning is raised outside 1.8–3.6 |
| **Het/Hom-alt ratio** | about 1.5–2.0 | Very high suggests contamination or a sample mix |
| Records | about 20–25k SNVs per exome, about 2–3k per Mb… | Far fewer suggests coverage or filter problems |

**SAMPLE_A chr22 subset (fastq_to_vcf, three callers):**

| Caller | Records | SNVs | Indels | Ti/Tv | Het/Hom |
|---|---|---|---|---|---|
| BCFtools | 1,487 | 1,323 | 164 | 2.49 | 1.68 |
| GATK HC | 1,467 | 1,283 | 183 | 2.61 | 1.69 |
| FreeBayes | 1,392 | 1,238 | 127 (+27 MNP/complex) | 2.67 | 1.65 |

All three agree on **1,314 sites** (exact allele match after normalisation). Each caller also has 31–53 calls the
others lack. Phase 4 will analyse why, by variant type, depth and region.

## 8. Filtering (soft, per caller)

After normalisation each VCF is **soft-filtered**: records failing a rule get a FILTER label (e.g. `LowDP`,
`HardFilterSNP`) and **are kept**. The output is `<sample>.<caller>.filtered.vcf.gz`.

Filtering is **per caller** for two reasons:
- QUAL is on a different scale in every caller.
- Only some callers annotate GQ, QD or SAF/SAR.

| Caller | Default rules |
|---|---|
| All | `NoAltGT` (genotype 0/0 or missing), `LowDP` (FORMAT/DP < 10), `LowVAF` (ALT fraction < 0.15) |
| BCFtools | `LowQual` (QUAL < 20), `LowGQ` (GQ < 20; `bcftools call -a GQ` writes GQ) |
| GATK HC | `LowGQ`, plus **GATK Best-Practices hard filters**. SNP: QD < 2, FS > 60, MQ < 40, SOR > 3, MQRankSum < −12.5, ReadPosRankSum < −8. Indel: QD < 2, FS > 200, SOR > 10, ReadPosRankSum < −20 |
| FreeBayes | `LowQual` (QUAL < 20), `StrandBias` (all ALT reads on one strand). FreeBayes writes no GQ |
| DeepVariant | Its own `RefCall`/`LowQual` labels are kept, plus `LowGQ` |

- **VQSR** (GATK's machine-learning filter) needs many samples or a whole genome, so it is not used for a single
  exome. That is why hard filters are the GATK fallback.
- Thresholds are editable per caller on the start form. The run report shows how many records each rule labelled.
- The report also shows PASS counts and PASS Ti/Tv. Ti/Tv should **rise** after filtering, because false positives
  have a low Ti/Tv.

## 9. Not included yet (deliberately)

* **BQSR** (base-quality recalibration): optional in modern GATK practice and of limited benefit on NovaSeq binned
  qualities. It can be added as an optional step; the known-sites files (dbSNP, Mills) are present.
* **Sex chromosomes:** all callers assume diploid. A male sample's chrX/chrY calls (outside the PARs) should be
  haploid. This needs a sex check and ploidy settings (planned).
* **gVCF / joint genotyping** for cohorts.
