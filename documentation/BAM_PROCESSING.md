# BAM Processing

This guide covers how trimmed reads become an analysis-ready BAM, and how an *existing* BAM (in-house or from a
provider) is prepared for variant calling.

```
trimmed FASTQ ─► bwa mem (+@RG) ─► samtools fixmate -m ─► samtools sort ─► aligned.sorted.bam
                                                                               │
                         duplicate MARKING (flag 0x400, nothing removed) ◄─────┘
                                   │
                                   ├─► markdup.bam + .bai   ← used by QC and all callers
                                   └─► (optional) dedup.bam ← a separate copy without duplicates
```

## 1. Alignment: BWA-MEM

`bwa mem` finds the best position of each read in the reference. Options used and why:

| Option | Why |
|---|---|
| `-R '@RG\tID:…\tSM:…\tLB:…\tPL:ILLUMINA\tPU:…'` | Adds a **read group** (see below). GATK refuses BAMs without one |
| `-K 100000000` | Processes a fixed number of bases per batch, so **results are identical whatever the thread count**. This matters for reproducibility |
| `-Y` | Soft-clips supplementary alignments (GATK convention). Keeps the full read sequence in the BAM |
| `-t <threads>` | Speed |

**BWA-MEM2** would be 2–3× faster, but its index needs about 87 GB of RAM to build and this server has 62 GB.
We use BWA-MEM 0.7.19. Its alignments are equivalent to BWA-MEM2, which was designed to give identical output.

### Read groups (`@RG`)

A read group identifies *which sequencing lane produced which reads*. The platform derives it from the first read name:

| Field | Value (SAMPLE_A) | Purpose |
|---|---|---|
| `ID` | `23HCHGLT4.3` (flow cell.lane) | Unique per lane. Base-quality models and optical duplicates are per lane |
| `SM` | `SAMPLE_A_chr22` (sample ID) | **Becomes the sample name in the VCF** |
| `LB` | sample ID (or *Library ID* parameter) | Duplicate marking is done **per library** |
| `PL` | `ILLUMINA` | Platform (error model) |
| `PU` | flowcell.lane.sample | Platform unit |

### Streaming: no shell, no huge temporary file

`bwa mem | samtools fixmate -m | samtools sort` runs as **one step**: three processes connected by OS pipes, not
a shell. The large unsorted BAM is never written to disk. If any of the three fails, the step fails
(`pipefail` behaviour), and the exit code of each process is recorded.

* `fixmate -m` adds mate information (MC/ms tags) that `samtools markdup` needs.
* `sort` orders reads by chromosome and position (*coordinate-sorted*). Indexing and callers need that order.

### Integrity checks after alignment

* `samtools quickcheck`: the BAM has a valid header and the end-of-file marker, so it is not truncated.
* **Read accounting:** the number of *primary* alignments must equal the number of reads that went in. Every input
  read produces exactly one primary record, mapped or not. A mismatch fails the run.
* A warning is raised if fewer than 90 % of reads map. Possible causes are contamination, the wrong species or
  reference, or a poor library.

## 2. Duplicate marking (not removal)

**What duplicates are:**

* **PCR duplicates:** the same original DNA fragment copied during library amplification and sequenced several times.
* **Optical duplicates:** one cluster read twice by the camera.

Duplicates are *not independent evidence*. Ten copies of a fragment carrying a PCR error look like a confident
variant.

**How they are detected:** read pairs with the same start positions and orientation for both mates. The copy with the
best base qualities is kept as the original. The others get SAM flag **0x400 (1024)**.

**Why mark, not remove:**

* Callers (GATK, bcftools, FreeBayes, DeepVariant) and coverage tools (mosdepth) **ignore flagged duplicates
  automatically**.
* The reads stay in the BAM, so the data is complete and auditable, and duplicate rate can be measured.
* If you need a BAM without duplicates (some external tools do), tick *Also write a duplicate-removed BAM*. A
  **separate** `dedup.bam` is created and the marked BAM is kept.

| Parameter | Default | Notes |
|---|---|---|
| Duplicate marking tool | samtools markdup | Fast. **GATK MarkDuplicates** (Picard) is the GATK Best Practices tool and is available as an option; results are very similar |
| Optical pixel distance | 2500 | Use 2500 for patterned flow cells (NovaSeq, NextSeq 2000) and 100 for older instruments |

Duplicate rate in exomes: 5–20 % is typical, and **above 30 %** gives a warning (low library complexity).
SAMPLE_A chr22: **10.7 %**.

After marking, the number of primary records is checked again: it must be unchanged.

## 3. Index

`samtools index` writes `markdup.bam.bai`. Without an index, a program has to read the whole file to find chr22.
The index is checked to be **newer than the BAM**. A stale index after a BAM is replaced can give wrong results
silently.

## 4. Analysis regions

These are built from the sample's registered **target BED** (the capture kit design, e.g. xGen Exome v2):

| File | Content | Used by |
|---|---|---|
| `regions/targets.bed` | Targets on the selected contigs, merged | Coverage (mosdepth) |
| `regions/calling_regions.bed` | Targets ± **padding** (default 100 bp), merged | Variant callers |

* **Padding** captures splice-site regions and the flanks of exons, where clinically relevant variants occur.
* **Contigs** defaults to `primary` (chr1–22, X, Y, M). ALT, decoy and unplaced contigs are excluded from calling.
* Use `chr22` (or a comma list) for quick tests.
* BED contigs that aren't in the reference (wrong naming, e.g. `1` vs `chr1`) give a warning. If none match, the step
  fails.

## 5. BAM QC

| Tool | Metrics |
|---|---|
| `samtools flagstat` | Total, primary, mapped %, properly paired %, duplicates, supplementary |
| `samtools stats` | Insert size (mean ± SD), error rate (mismatches vs reference), read lengths |
| `mosdepth` (targets, MAPQ ≥ 20, duplicates excluded) | **Mean target coverage** and **% target bases ≥ 1/10/20/30/50/100×** |

**How to judge an exome:**

* Mean target coverage ≥ 100× is typical for clinical WES.
* **% bases ≥ 20× is the key metric.** Many labs require ≥ 95 %. Below 90 % gives a warning, because variants in
  poorly covered bases can be missed.
* Properly paired above 95 %, and an error rate below about 0.5 %.

**SAMPLE_A full exome (run WF-20261106-082225-5eeddc, 34.1 Mb of targets):**

* Mean target coverage 140×, duplicates 10.9 %, properly paired 99.1 %, error rate 0.23 %.
* **94.6 % of target bases ≥ 20×**, but **4.6 % have no usable coverage at all** (MAPQ ≥ 20).

**Where the uncovered bases are, and why (a real limitation, not a sequencing problem):**

| Region | Target bases with no MAPQ ≥ 20 reads |
|---|---|
| chr6, of which the **MHC/HLA region** (28.5–33.5 Mb) | 15.8 % of chr6 targets; 94 % of these are in the MHC |
| chr17, chr15, chr16, chr19, chr22 | 5–9 % each |
| chrY | 99.8 %. Expected if the sample is female; needs a sex check |

Our `hg38.fa` contains 261 **ALT contigs** (alternative versions of regions such as the MHC), but BWA was given no
`.alt` file, so the alignment is **not ALT-aware**. A read from such a region matches the primary chromosome and
its ALT copy equally well. BWA therefore gives it **MAPQ 0**, and coverage and all callers ignore it. DRAGEN handles
ALT contigs, which is why the provider BAM reaches 96.7 % ≥ 20× on chr22 against our 93.4 %.

**After the fix (run WF-20261119-001642-08530a, same FASTQ and settings, reference `hg38_dragen_masked`):**

| Metric | `hg38.fa` (not ALT-aware) | `hg38_dragen_masked` | Provider DRAGEN BAM |
|---|---|---|---|
| Mean target coverage | 139.6× | **144.5×** | 143.0× |
| % target bases ≥ 1× | 95.37 % | **98.36 %** | |
| **% target bases ≥ 20×** | 94.61 % | **97.66 %** | 97.33 % |
| % target bases ≥ 100× | 76.65 % | 79.24 % | |

Only the reference changed, so the +3 percentage points are entirely the ALT effect. About 1 Mb of target
sequence (mostly the MHC on chr6, and chr15–17, 19) became usable. The remaining ~1.6 % of target bases
without coverage include chrY (expected for a female sample) and genuinely repetitive regions.

> An earlier version of this guide blamed the chr22 gap on the test subset. The full exome shows the same 5.9 % gap on
> chr22, so that explanation was wrong.

## 6. Using an existing BAM (BAM → VCF workflow)

The **Input BAM preparation** step never modifies the registered file:

| Check / action | Behaviour |
|---|---|
| `samtools quickcheck` | Truncated or corrupt BAM, so the run fails |
| Genome build vs reference | **GRCh37 BAM on a GRCh38 reference makes the run fail.** Builds are never mixed |
| Extra contigs (e.g. DRAGEN decoys/HLA) | Warning. GATK may refuse such BAMs |
| Coordinate-sorted? | Fails if not |
| Duplicates marked? (looks for markdup / MarkDuplicates / DRAGEN in `@PG`) | Warning if not. Duplicates would count as evidence |
| Read group | If missing and GATK is selected, a **copy with an @RG** is written into the run folder. This needs disk space about equal to the BAM size |
| Index | The BAM is **symlinked** into the run folder. An existing `.bai` is linked too, or a new index is created **in the run folder**. Nothing is written next to the original file |

### Example: the provider's DRAGEN BAM for SAMPLE_A (`bam_to_vcf`, chr22)

| Check | Result |
|---|---|
| Build | GRCh38. Verdict **primary_compatible**: 2,911 extra DRAGEN decoy/HLA contigs |
| Index | None was provided. A `.bai` was created **in the run folder**; the BAM itself is only symlinked |
| Duplicates | Already marked by DRAGEN (11.3 %) |
| Read group SM | `SAMPLE_A_provider`, not `SAMPLE_A`. A warning says the VCF will carry the BAM's sample name |
| Coverage, chr22 targets | Mean 160×, **96.7 % ≥ 20×** (our alignment: 93.4 %, lower because it is not ALT-aware, see §5) |
| GATK on a BAM with extra contigs | Accepted |

**Lesson from real data:** `samtools stats -r hg38.fa` **failed** on this BAM because it has reads on
`chrUn_KN707606v1_decoy`, which is not in our reference. The failure was reported, not hidden. Since workflow v1.1,
`samtools stats` is restricted to the analysis contigs, so external BAMs work and QC is comparable between BAMs.
