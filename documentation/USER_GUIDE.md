# User Guide (Phase 3)

## 1. Open the platform
Start the SSH tunnel (see INSTALLATION_GUIDE.md §3) and open http://localhost:5000.

## 2. Create a sample
**Sample Management → New sample.** Enter the Sample ID (e.g. `SAMPLE_A`), sequencing type and layout.
The ID becomes a folder name, so only letters, digits, `_ - .` are allowed.

## 3. Register input files
On the sample page, under **Input files**:
1. Click **Browse** to pick a file on the server (or paste its absolute path).
2. Role: leave *Auto-detect* (R1/R2 detected from `_R1_`/`_R2_` in the name) or choose explicitly.
3. Origin: **In-house** for our own data, **External provider** for files received from another lab.
4. Optional label, e.g. "DRAGEN 4.3 provider X".
5. For exomes/panels also register the capture kit's **target BED** (role *Target regions BED*), e.g.
   `/data/xgen-exome-hyb-panel-v2-targets-hg38.bed`. It defines where coverage is measured and
   where variants are called.

Files are **referenced where they are** — never copied, moved or modified. A file used by a run cannot be
unregistered (traceability).

## 4. Run intake validation
Under **Start a workflow**: choose *Sample intake validation*, the reference, tick the files, **Start**.
You are taken to the run page.

## 5. The run page (monitoring)
* **Status badge**, **overall progress** (= finished steps ÷ total steps — never simulated),
  **elapsed time**, **current stage** and the **exact command** being executed.
* **Processing steps** table with per-step status and duration. Click a step to show only its log.
* **Live log** of all tool output (stdout + stderr).
* Buttons: **Pause / Continue** (external tools only), **Cancel** (terminates tools, keeps files),
  **Resume** (after failure/cancel: finished steps are re-validated and kept, the rest re-run),
  **Output files** (browse/download).
* **Warnings** are shown in yellow. A warning is something to *review*, not automatically an error.
* At the bottom: every executed command with its exit code and tool versions, and the audit trail.

Estimated remaining time is deliberately not shown yet — it will appear once enough run history exists
to estimate it reliably.

## 6. Intake results
On the sample page each file now shows its **SHA-256**, **detected build** and **validation status**:

| Status | Meaning |
|---|---|
| OK | Intact and compatible with the reference |
| WARNING | Usable with caveats (e.g. external BAM with extra decoy/HLA contigs) |
| ERROR | Corrupt, malformed, or **different genome build** — must not be compared directly |

The full report is `reports/intake_report.html` in the run's output files.
The checksum file can be re-verified at any time on the server:
```bash
sha256sum -c results/<SAMPLE>/<RUN_ID>/intake/checksums.sha256
```

## 7. Processing workflows
Choose the workflow under **Start a workflow** on the sample page. The input files are pre-selected
by type (FASTQ + BED, or one BAM + BED); untick what you don't want. Parameters for the chosen workflow
appear below, grouped (Trimming, Alignment, Regions, Variant calling); rarely changed ones are under
*Advanced*. Every value is validated before the run is created and stored with the run.

| Workflow | Use it to… | Typical time (full WES / chr22 subset) |
|---|---|---|
| FASTQ QC & trimming | check a new FASTQ pair | ~15 min / 30 s |
| FASTQ → analysis-ready BAM | produce a duplicate-marked, indexed BAM + coverage QC | ~1.5-2 h / 1.5 min |
| BAM → VCF (germline) | call variants on an existing BAM (in-house or provider) | depends on callers |
| FASTQ → BAM → VCF | the complete pipeline | ~3-5 h / 4 min |

**Quick test run:** set *Contigs to analyse* = `chr22`. **Variant callers:** tick one or more; each
caller's own parameters appear when it is ticked, including its **Filtering** thresholds. Callers that
cannot run on this server are shown greyed out with the reason. With two or more callers, the run also
compares them (*Caller comparison* step). *Keep intermediate files* (Outputs group) is on by default;
off deletes trimmed FASTQs and caller temp folders after a successful run (listed in
`reports/cleanup_manifest.json`).

**Results** (run page → *Output files*):
* `reports/run_report.html` — read accounting, trimming, alignment/coverage, variant counts per caller,
  warnings, tool versions, parameters.
* `reports/multiqc/multiqc_report.html` — all QC plots in one page.
* `bam/<sample>.markdup.bam` (+ `.bai`) — the analysis-ready BAM (open in IGV).
* `vcf/<caller>/<sample>.<caller>.filtered.vcf.gz` — **the VCF to use**: normalised and soft-filtered
  (FILTER column; nothing removed). The raw and normalised VCFs are kept beside it.
  With a preset the preset is in the name: `<sample>.<preset>.<caller>.filtered.vcf.gz`.
* `comparison/callers/comparison_report.html` — caller concordance (2+ callers).

The run page shows links to all reports at the top, and the **Known limitations** box.
At the end of a germline run its BAM and filtered VCFs are **registered as sample files**
(origin *Platform pipeline*), so they can be picked for comparisons straight away.

Concepts and how to judge the numbers: FASTQ_PROCESSING.md, BAM_PROCESSING.md, VARIANT_CALLING.md.
The module pages (FASTQ Analysis, BAM Analysis, Variant Calling) list the workflows, their steps,
caller availability and recent runs.

## 8. Comparison workflows (Caller Comparison / External Data Comparison)
Started from the sample page like the others; register the provider's files with origin **External provider**.

| Workflow | Select | Result |
|---|---|---|
| VCF vs VCF comparison | 2–8 VCFs (+ BED; + Truth VCF and Truth BED for accuracy) | concordance matrix, intersections, genotype concordance, RTG vcfeval, precision/recall if a truth set is given |
| BAM vs BAM comparison | 2–6 BAMs (+ BED) | side-by-side QC table, MAPQ distribution, reads per contig class, targets covered differently |
| BAM vs VCF read-level validation | 1 BAM + 1 VCF (+ BED) | read evidence per variant (VAF, strand, MAPQ, soft-clips), evidence notes, IGV batch file |
| **VCF variant-by-variant check** (new) | exactly 2 VCFs (+ BED); choose the **baseline** | for every baseline variant: found & PASS / found but filtered (which filter) / genotyped as reference / written differently / different variant nearby / not called, with both files' QUAL·DP·AD·VAF·GQ side by side; interactive viewer |

Files from another genome build are refused with an explanation. How to read the results:
VARIANT_CALLER_COMPARISON.md and BAM_VCF_VALIDATION.md. The module pages list recent comparison
results with their key number and a link to the report.

### 8a. Variant-by-variant check: "Is every Provider variant in our VCF?"

1. On the sample page choose **VCF variant-by-variant check**, tick exactly two VCFs (the provider's and ours) and the
   target BED, and pick the provider's VCF as **Baseline**.
2. Keep *Our VCF: include filtered records* on (default). Then a variant our filters removed is still found and shown as
   **found, filtered** with the filter's name, instead of disappearing into "not called". This is why our VCF may be
   made with *less stringent* parameters: nothing is lost, you filter afterwards.
3. When the run has finished, open **Open variant-by-variant viewer** on the run page:
   * left box: thresholds for *our* VCF (depth, GQ, VAF, QUAL; optionally ignore our FILTER column). The big number
     at the top shows how many baseline variants would count as found with exactly these thresholds — a what-if,
     without re-running anything;
   * middle box: which baseline variants to include (e.g. only baseline QUAL ≥ 30, only SNVs, one chromosome);
   * right box: which rows to list (by status, or only the ones missing under your thresholds);
   * **CSV of these rows** downloads exactly what is shown (all pages).
4. The run report (HTML) has the same summary plus the first 150 variants of each "missing" category.

The existing *VCF vs VCF comparison* is unchanged; this is a separate workflow for the one-directional question.

Variant-check files are named after both VCFs, e.g.
`SAMPLE_A.provider-like-dragen-4-3.gatk_hc.filtered_vs_SAMPLE_A_provider.hard-filtered.variant_check.tsv.gz`
(plus `.ours_pass_not_in_baseline.tsv.gz`, `.variant_check_summary.json`, `.variant_check_report.html`).

#### Searching, MQ and "Check reads"

* **Find a variant**: type `chr1:100000 A>G`, `chr1:100000` or a range; ± bp catches shifted indels.
* **MQ columns**: the site mapping quality (INFO/MQ) of both VCFs. A red baseline MQ (< 20) means
  the baseline's reads at that site fit equally well somewhere else. The *baseline MQ* min/max
  boxes select such sites (e.g. max 10 → "variants the baseline called from ambiguously placed
  reads"). DRAGEN's VCF MQ is on its own scale (values up to 250); BWA/GATK's stops at 60.
  DeepVariant writes no MQ.
* **Check reads** (every row): reads our BAM at that position, read-only, and shows how many reads
  carry REF/ALT, their MAPQ, how many are usable by the callers (MAPQ ≥ 20, base quality ≥ 10)
  and where else the aligner placed them (BWA `XA` tag), with the likely reason a variant is
  missing. You can pick another BAM of the sample (it must be indexed).

**Worked example — a variant inside the chr11 region duplicated by `chr11_KI270721v1_random`**
(called homozygous by the provider with good quality; not in our VCF): every read in our BAM
carries the ALT base, but every one has **MAPQ 0** because it aligns equally well to
`chr11_KI270721v1_random`, an unplaced GRCh38 contig that repeats this part of chr11 (*H19*
region). Callers ignore MAPQ-0 reads, so no caller saw a usable read. DRAGEN placed the same reads
with MAPQ 60: its 4.3 reference was built with Illumina's newer mask
(`hg38_alt_masked_v3/hg38_v2_mask…bed`), which also masks such duplicated unplaced sequence; the
public DRAGMAP mask we used (v1, ALT contigs only) does not. This is a reference difference, not a
filter: no threshold in our pipeline can recover these sites.

**CSV in Excel**: Excel turns a genotype like `1/1` into a date ("1-Jan"). Import the CSV with the
GT columns as *Text*, or open it in LibreOffice/a text editor.

## 8b. Parameter presets

*Parameter Presets* (sidebar) holds named parameter sets, e.g. the built-in **In-house default** (locked) and
**Provider-like (DRAGEN 4.3)**.

* **Use one:** on the sample page, choose a preset above the file list. It fills the fields of the selected
  workflow; fields the preset does not contain keep their values. You can still change anything. The run page shows
  the preset and lists every value you changed after applying it, so a run is always reproducible.
* **Create one:**
  * from the Start form: set the fields, type a name next to *Save as preset*, click it (nothing is started);
  * from a finished run: *Save parameters as preset* on the run page;
  * from a provider VCF: *Grab the filtering profile of a VCF* on the Presets page (next section);
  * *Duplicate* any preset and edit the copy (JSON, checked against the same rules as the form).
* Presets are validated: an unknown parameter or an out-of-range value is refused with the reason.

### Built-in "Raw tool defaults (no filters)"

Every caller runs its **bare command** with the tool's own defaults — for GATK
`gatk HaplotypeCaller -R ref -I bam -O vcf` — and all soft filters are off. Each caller has a
*Raw command: use the tool's own defaults* switch (that is all this preset turns on); its help text
lists the tool defaults. Still passed: the calling regions (`-L`/`-R`/`-t`/`--regions`; to call
the whole BAM, start the run without a BED), thread count, temp folder and Java heap — they
change speed, not calls. bcftools gets `call -m -v` (it needs a model) and still writes AD/DP/GQ
annotations. The exact command is in the run's command log.

Tool defaults still ignore low-MAPQ reads (GATK < 20, DeepVariant < 5, FreeBayes < 1), so "raw"
does **not** recover MAPQ-0 sites such as the chr11 example above.

## 8c. Grabbing the filtering profile of a VCF

On the Presets page, give the path of a VCF (e.g. `SAMPLE_A/vcf/SAMPLE_A_provider.hard-filtered.vcf.gz`). The platform reads:

| From | What | Example (Provider DRAGEN 4.3) |
|---|---|---|
| Header | Caller and command line | DRAGEN 4.3.6; `--vc-target-bed` **not set** → whole genome called |
| Header | Every `##FILTER` rule that states its threshold | `LowDepth: DP <= 1`, `DRAGENSnpHardQUAL: QUAL < 0.4`, `MosaicLowAF: AF < 0.2` |
| Records | Which filters were really applied (headers define many unused ones) | PASS 137,725 · LowDepth 69,580 · MosaicLowAF 4,241 · 18 somatic filters never used |
| Records | Lowest / 1 % / 5 % / median QUAL, DP, GQ, VAF of PASS calls = the *effective* thresholds | PASS calls with DP 0–1 and GQ 0 exist |

It then proposes values for our filters, each labelled **mapped** (same meaning), **approximate** (QUAL: every
caller computes it differently) or **not mapped / not mappable** (e.g. DRAGEN's mosaic and ploidy filters). For
Provider: depth ≥ 2, GQ/VAF filters off, GATK hard filters and FreeBayes strand filter off, GATK emission threshold
QUAL ≥ 1 (in-house 30). The preset changes **filtering only**; alignment and calling settings stay as the form has
them, because the VCF does not record them.

> A "Provider-like" preset makes our output *comparable* to theirs. It does not make it more accurate — lenient
> thresholds also pass more false positives. Use it for comparison, keep the in-house preset for reporting, and decide
> with a truth sample (GIAB).

## 9. Dashboard
Counts by status, active jobs with live progress, recent activity, CPU/RAM/disk (refreshes every 5 s),
and the **Known limitations** of the current setup (e.g. non-ALT-aware reference) — these are also
printed in every report.

## 10. Reports
*Reports* lists every completed run with links to its run report, MultiQC and comparison reports, and
offers CSV exports across all runs: QC metrics per run, variant counts per caller (before/after
filtering), pairwise concordance.

## 11. Settings
Reference genomes (integrity check results, register another reference), runtime defaults (threads,
disk safety margin, default for intermediate-file retention, GATK launcher and its Java 17 home — empty =
auto-detect, DeepVariant Docker image (pinned version), RTG vcfeval memory), detected tools and versions, paths.
