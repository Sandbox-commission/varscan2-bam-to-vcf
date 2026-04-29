# VarScan2 Somatic Variant and CNV Analysis Pipeline

Paired case/control somatic variant calling (SNP, INDEL, Germline, LOH) and
copy number variation (CNV) analysis using VarScan2, samtools, bam-readcount,
and fpfilter. Implemented as a single-binary Rust pipeline. Designed for whole
exome sequencing (WES) with full WGS compatibility.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Stages](#pipeline-stages)
- [Prerequisites](#prerequisites)
- [Directory Structure](#directory-structure)
- [Input File Format](#input-file-format)
- [Configuration](#configuration)
- [Parameter Reference](#parameter-reference)
  - [Somatic Calling Parameters](#1-somatic-calling-parameters-varscan-somatic)
  - [processSomatic Parameters](#2-processsomatic-parameters)
  - [CNV Parameters](#3-cnv-parameters)
  - [copyCaller Parameters](#4-copycaller-parameters)
  - [Data-Ratio Calculation](#5-data-ratio-calculation)
  - [bam-readcount Parameters](#6-bam-readcount-parameters)
  - [fpfilter Parameters](#7-fpfilter-parameters)
  - [Summary Table](#summary-table)
- [Usage](#usage)
- [Output Files](#output-files)
- [Key Design Decisions](#key-design-decisions)
- [Post-Pipeline Analysis](#post-pipeline-analysis)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)
- [License](#license)

---

---

## Overview

```
BAM files
    |
    v
[Stage 1]  samtools flagstat       — alignment QC; read counts for data-ratio
    |
    v
[Stage 2]  samtools mpileup        — paired pileup  (normal FIRST, tumour SECOND)
    |
    v
[Stage 3]  VarScan somatic         — raw variant calls (.snp.vcf, .indel.vcf)
    |
    v
[Stage 4]  VarScan processSomatic  — classify: Somatic / Germline / LOH (.hc.vcf)
    |
    v
[Stage 5]  VarScan copynumber      — per-window log2 depth ratios (.copynumber)
    |
    v
[Stage 6]  VarScan copyCaller      — segment CNV calls (.called, .homdel)
    |
    v
[Stage 7]  awk                     — convert .hc.vcf to position lists (.var)
    |
    v
[Stage 8]  bam-readcount           — per-position allele read support
    |
    v
[Stage 9]  fpfilter.pl             — false positive removal (.fpfilter.vcf)
    |
    v
[Stage 10] Summary report
```

---

## Pipeline Stages

| Stage | Tool | Input | Output |
|-------|------|-------|--------|
| 1 | `samtools flagstat` | BAM files | `flagstats/*.flagstats` |
| 2 | `samtools mpileup` | BAM pairs + reference | `mpileup/*.mpileup` |
| 3 | `VarScan somatic` | paired mpileup | `somatic/*.snp.vcf`, `somatic/*.indel.vcf` |
| 4 | `VarScan processSomatic` | raw VCFs from Stage 3 | `somatic/*.{Somatic,Germline,LOH}.hc.vcf` |
| 5 | `VarScan copynumber` | paired mpileup + flagstats | `copynumber/*.copynumber` |
| 6 | `VarScan copyCaller` | `.copynumber` files | `*.copynumber.called`, `*.copynumber.homdel` |
| 7 | `awk` | `.hc.vcf` from Stage 4 | `snp-VAR/*.var` (all SNP classes), `indel-VAR/*.var` (Somatic INDELs only) |
| 8 | `bam-readcount` | `.var` files + BAMs | `readcount/*.readcount` |
| 9 | `fpfilter.pl` | readcounts + VCFs | `filtered/*.fpfilter.vcf` |
| 10 | — | all outputs | `varscan_pipeline_summary.txt` |

---

## Prerequisites

| Tool | Minimum version | Notes |
|------|----------------|-------|
| bash | 4.0 | Requires `(( ))` arithmetic and `[[ ]]` conditionals |
| Java JRE | 8 | Required to run VarScan2 .jar |
| VarScan2 | 2.3.9 | `VarScan.v2.3.9.jar` placed in `$SOFTWAREDIR` |
| samtools | 1.13 | Must be on `$PATH` |
| bam-readcount | 0.8 | Must be on `$PATH` |
| Perl | 5.10 | Required for fpfilter.pl |
| fpfilter.pl | — | VarScan2 companion script; place in `$SCRIPTSDIR` |
| bc | any | Floating-point data-ratio calculation |

All BAM files must be:
- **Coordinate-sorted** — required by samtools mpileup and bam-readcount
- **Indexed** — `.bai` file present alongside each BAM (checked at Stage 8)
- **Duplicate-marked** — PCR/optical duplicates flagged by `MarkDuplicates`; the
  pipeline excludes them via `-F 0x400` in mpileup
- **Chromosome naming consistent with the reference FASTA** — Ensembl uses `1, 2,
  ..., X, Y, MT`; UCSC uses `chr1, chr2, ..., chrX, chrY, chrM`. A mismatch
  produces silent empty output at every stage. Verify with:
  ```bash
  samtools view -H sample.bam | grep '^@SQ' | head -3
  ```

---

## Build

Install Rust (if not already installed):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

Build the pipeline binary:

```bash
cd /path/to/varscan
cargo build --release
```

Binary location: `target/release/varscan2_pipeline`

> The legacy `varscan2_pipeline.sh` is preserved for reference but the Rust
> binary is canonical. All new development targets the Rust binary.

---

## Obtaining VarScan2 and fpfilter.pl

**VarScan.v2.3.9.jar** — place in `software/`:

```bash
wget -O software/VarScan.v2.3.9.jar \
  https://github.com/dkoboldt/varscan/releases/download/2.3.9/VarScan.v2.3.9.jar
```

**fpfilter.pl** — place in `scripts/`:

```bash
wget -O scripts/fpfilter.pl \
  https://raw.githubusercontent.com/genome/fpfilter-tool/master/fpfilter.pl
chmod +x scripts/fpfilter.pl
```

Java JRE 8+ is required for VarScan (see Prerequisites table above).

---

## Directory Structure

```
<working_directory>/
|-- varscan_pipeline.sh        # this script
|-- sample_pairs.csv           # pairs file (tumour col 1, normal col 2)
|
|-- flagstats/                 # samtools flagstat outputs
|-- mpileup/                   # paired .mpileup files
|-- somatic/                   # VarScan somatic and processSomatic outputs
|-- copynumber/                # VarScan copynumber and copyCaller outputs
|-- snp-VAR/                   # SNP position lists (.var) for bam-readcount
|-- indel-VAR/                 # INDEL position lists (.var) for bam-readcount
|-- readcount/                 # bam-readcount outputs
|-- filtered/                  # fpfilter.pl output VCFs
|-- varscan_pipeline_summary.txt
```

---

## Input File Format

The pairs file is a plain CSV with **one tumour/normal pair per line**.
Column 1 is the **tumour** identifier; column 2 is the **normal** identifier.

> **Column order is critical.** Swapping tumour and normal silently reverses all
> Somatic / Germline / LOH classifications because VarScan assigns sample roles
> by column position in the mpileup (normal must be first in the mpileup call,
> which is determined by which entry is treated as normal here).

The sample name is derived by stripping `PAIRS_SUFFIX` from each identifier.
The BAM file is then located as `$BAM_DIR/<sample_name>$BAM_SUFFIX`.

**Incomplete pairs** (empty normal field) are automatically skipped with a
warning. Verify all pairs before running; do not omit unpaired tumour-only samples.

### Example A — typical `sample_pairs.csv`

```
case01_final.bam,control01_final.bam
case02_final.bam,control02_final.bam
```

```bash
FILE_PAIRS_LIST="sample_pairs.csv"
PAIRS_SUFFIX="_final.bam"
BAM_SUFFIX="_final.bam"
BAM_DIR="$PWD"
```

> Samples with no matched control are skipped automatically with a WARNING message.

### Example B — pairs file lists flagstat filenames; BAMs are in a subdirectory

```
# file_pairs.txt  (tumour first, normal second)
tumour1.flagstats,normal1.flagstats
tumour2.flagstats,normal2.flagstats
```

```bash
FILE_PAIRS_LIST="file_pairs.txt"
PAIRS_SUFFIX=".flagstats"
BAM_SUFFIX="_sm.bam"
BAM_DIR="$PWD/sm-bam"
```

The two suffix variables decouple the pairs file format from the BAM naming
convention, so the same script handles any combination without modification
to the stage logic.

---

## Configuration

> **Rust binary:** Parameters are compile-time constants at the top of
> `varscan2_pipeline.rs`. Edit the constant, then rebuild with
> `cargo build --release`. This ensures the binary and its configuration are
> always in sync and avoids runtime config-file parsing errors.
>
> **`GENOMEIDX1` is mandatory.** The binary will refuse to start with a clear
> error message if it is left empty. Set it to the absolute path of your
> indexed GRCh38 reference FASTA before building.
>
> The legacy bash script (`varscan_pipeline.sh`) uses shell variables in its
> `CONFIGURATION` section — the table below applies to both.

### Path settings

| Variable | Description |
|----------|-------------|
| `GENOMEIDX1` | Full path to the indexed reference FASTA (**must be set**) |
| `SOFTWAREDIR` | Directory containing `VarScan.v2.3.9.jar` |
| `SCRIPTSDIR` | Directory containing `fpfilter.pl` |
| `BAM_DIR` | Directory containing all BAM files |
| `BAM_SUFFIX` | Filename suffix of BAM files (e.g. `.bqsr.bam`, `_sm.bam`) |
| `FILE_PAIRS_LIST` | Path to the tumour/normal pairs CSV |
| `PAIRS_SUFFIX` | Extension to strip from pairs-file entries to get sample name |
| `TARGET_BED` | **WES:** Path to capture kit BED file (e.g. Agilent SureSelect, Twist Exome). Passed as `-l` to `samtools mpileup` — restricts pileup to captured regions, reduces mpileup file sizes by ~95%, and naturally excludes alt contigs and decoy sequences. Leave empty for WGS. |
| `VCF_SAMPLE_LIST` | Optional: plain-text file with sample names (normal first, tumour second, one per line). Provides correct column labels in VCF `##SAMPLE` headers instead of generic `NORMAL`/`TUMOR`. Leave empty to use VarScan defaults. |

For detailed justification of every numeric parameter see the
[Parameter Reference](#parameter-reference) section below.

---

## Parameter Reference

---

## 1. Somatic Calling Parameters (`VarScan somatic`)

---

### `MIN_COVERAGE = 20`

**What it does:** Minimum total read depth (normal + tumour combined) at a
position for it to be considered for variant calling.

**Justification:** At fewer than 20× coverage, allele frequency estimates become
statistically unreliable. For a heterozygous somatic mutation at 50% VAF, 10×
coverage gives a binomial standard error of ~16%, making it nearly impossible to
distinguish a true 10% somatic mutation from noise. At 20×, there are enough
observations to apply a meaningful statistical test. This is the widely accepted
minimum for somatic calling in whole-exome sequencing (WES). For whole-genome
sequencing (WGS) at shallower depths (10–15×), this may be lowered, but at the
cost of sensitivity and specificity.

---

### `MIN_COVERAGE_NORMAL = 10` and `MIN_COVERAGE_TUMOR = 20`

**What they do:** Minimum depth required independently in each sample before a
site is considered.

**Justification:** These are set deliberately lower than `MIN_COVERAGE` because
the combined threshold already enforces adequate total depth. Setting separate
per-sample floors prevents the edge case where one sample dominates coverage —
for example, a site with 25× in normal and 2× in tumour would pass
`MIN_COVERAGE=20` (combined) but the tumour VAF estimate from 2 reads is
meaningless.

`MIN_COVERAGE_TUMOR=20` is set to 20 for WES, where exon capture uniformity is
high and 20× tumour depth is routinely achieved. At `MIN_VAR_FREQ=0.10`, a
20× minimum guarantees at least 2 variant-supporting reads for any called
variant — the absolute minimum for meaningful somatic evidence. At lower depths
(e.g. 8×), a single miscalled read can satisfy the 10% VAF threshold, making
the call entirely unreliable.

`MIN_COVERAGE_NORMAL=10` is set lower than the tumour floor because the normal
sample requires sufficient depth to assess germline status, but false negatives
in the normal (missing a germline variant due to low depth) produce a more
recoverable error (germline called as somatic) than false negatives in the
tumour (missing a true somatic mutation entirely).

For shallow WGS (10–20×), both values can be reduced proportionally.

---

### `MIN_BASE_QUAL = 20`

**What it does:** Minimum Phred-scaled base quality applied at three points in
the pipeline:
- `samtools mpileup -Q` — excludes bases below this quality from the pileup
- `VarScan somatic --min-base-qual` — applied internally per mpileup column
- `VarScan copynumber --min-base-qual` — applied per CNV window

**Justification:** 20 is the VarScan2 official default and corresponds to a
per-base error probability of 1% (`10^(-20/10)`). Applying it consistently
across all three points ensures the read sets used by mpileup, somatic calling,
and CNV analysis are on exactly the same quality-filtered basis. Without the
`-Q` flag in mpileup, bases with quality 13–19 enter the pileup but are then
silently excluded by VarScan's internal filter — producing a discrepancy
between the total depth VarScan reports and what a user would count in IGV.

---

### `MIN_VAR_FREQ = 0.10`

**What it does:** Minimum variant allele frequency (VAF) in the tumour for a
variant to be reported at all. A variant with VAF below this threshold is not
emitted by VarScan.

**Justification:** 0.10 (10%) is the **VarScan2 official default** for the
somatic command (`--min-var-freq [0.10]` per the manual). It represents a
validated lower bound for reliable somatic detection at standard WES depths
(20–100×). Illumina base error rates are ~0.1–0.5% per base; at 20× coverage,
random errors can produce apparent VAFs of 2–5%. The 10% floor provides a
2–5× safety margin above the noise ceiling. The VarScan2 publication
(Koboldt et al. 2012) benchmarks sensitivity and specificity at this threshold.

For ultra-deep sequencing (>500×, liquid biopsy), lower to 0.01–0.02 where
statistical power at high depth is sufficient to distinguish true low-VAF
variants from noise.

---

### `MIN_FREQ_FOR_HOM = 0.75`

**What it does:** Minimum VAF required to call a variant homozygous (as opposed
to heterozygous).

**Justification:** For a truly homozygous variant in a diploid cell, the
expected VAF is 1.0. In practice, read sampling, alignment errors, and residual
normal contamination in tumour samples push this below 1.0. Setting the
threshold at 0.75 means VarScan calls homozygous only when at least 75% of
reads support the alt allele — accounting for up to ~25% normal cell
contamination or sequencing noise while still being distinct from the 0.50
expected for heterozygous variants. Setting this too low (e.g. 0.60) risks
calling heterozygous LOH events as homozygous.

---

### `NORMAL_PURITY = 1.0` and `TUMOR_PURITY = 1.0`

**What they do:** Inform VarScan of the expected proportion of tumour and normal
cells in each sample. Used internally to adjust VAF expectations during
classification.

**Justification:** `NORMAL_PURITY=1.0` is almost always correct — blood or
adjacent tissue normals are essentially pure diploid cells. `TUMOR_PURITY=1.0`
is the default starting assumption; however, most solid tumours have 20–80%
tumour cell content (the rest being stromal and immune cells). Setting it to
1.0 means VarScan does not adjust its somatic/germline decision boundary for
tumour purity.

If tumour purity is known from pathology estimates or computational tools
(ABSOLUTE, PURPLE, TitanCNA), set `TUMOR_PURITY` accordingly. For a 60% pure
tumour, a true heterozygous somatic mutation will appear at ~30% VAF rather
than 50%. With `TUMOR_PURITY=1.0`, VarScan may under-classify such variants.
The `processSomatic` step with appropriate `MIN_TUMOR_FREQ` compensates for
this partially.

---

### `P_VALUE = 0.99`

**What it does:** The p-value threshold for calling a position as a variant at
all (the initial detection gate). A site passes if the probability of being
non-reference exceeds this threshold — equivalently, variants with calling
p-value ≤ 0.99 are emitted.

**Justification:** This is the VarScan2 default and is **intentionally
permissive**. At this stage, VarScan is asking: *"Is there any statistically
meaningful deviation from the reference?"* A permissive gate casts a wide net,
capturing both true variants and potential false positives. The critical filters
are applied downstream:

- `SOMATIC_P_VALUE` controls whether a detected variant is classified as somatic
- `processSomatic` separates high-confidence from low-confidence classifications
- `fpfilter` removes strand-bias and mapping artefacts

If `P_VALUE` were set to a stringent value (e.g. 0.005), true somatic mutations
with moderate read support would be discarded before they even reach
classification — a false negative that cannot be recovered downstream.

---

### `SOMATIC_P_VALUE = 0.05`

**What it does:** The significance threshold for classifying a detected variant
as **somatic** (present in tumour, absent in normal). This is the Fisher's
exact test p-value comparing tumour vs normal allele counts.

**Justification:** 0.05 is the conventional significance level and appropriate
here because:

1. It is applied per-variant, not genome-wide — VarScan somatic is not
   performing a single genome-wide test but rather a targeted comparison at
   called sites
2. The subsequent `processSomatic` step provides a second layer of filtering
   with `MIN_TUMOR_FREQ` and `MAX_NORMAL_FREQ`, so the combined filter is more
   stringent than 0.05 alone
3. A more stringent threshold (e.g. 0.01) would reduce sensitivity for
   subclonal mutations with low read support

For projects requiring strict FDR control, Bonferroni correction or a
genome-wide FDR analysis on the output VCF is recommended rather than
tightening this threshold.

---

### `STRAND_FILTER = 1`

**What it does:** When enabled, VarScan excludes variants where more than 90%
of supporting reads come from a single strand.

**Justification:** Strand bias is one of the strongest indicators of a
sequencing artefact, particularly:

- **FFPE artefacts:** C→T / G→A transitions caused by cytosine deamination
  appear predominantly on one strand
- **PCR amplification bias:** Low-complexity regions can show strand-specific
  amplification errors
- **Library preparation artefacts:** End-repair and A-tailing introduce
  strand-specific errors near read ends

Disabling strand filter (`STRAND_FILTER=0`) is only appropriate for amplicon
sequencing where by design all reads originate from the same strand relative
to the amplicon.

---

## 2. processSomatic Parameters

---

### `MIN_TUMOR_FREQ = 0.10`

**What it does:** Minimum VAF in the tumour for a somatic variant to be
classified as **high-confidence** (`.hc.vcf`).

**Justification:** After the initial calling (Stage 3), all variants above
`MIN_VAR_FREQ=0.10` are emitted. The `MIN_TUMOR_FREQ=0.10` threshold for
high-confidence classification matches the calling floor, ensuring that only
variants meeting both the emission threshold and the classification threshold
receive the `.hc` label. Acknowledging that:

- Very low VAF calls near the 10% threshold in standard WES carry higher false positive rates
- The `.hc.vcf` files feed directly into bam-readcount and fpfilter — a cleaner
  input set reduces unnecessary computation
- For low-purity tumours, lowering this to 0.05 is biologically justified; for
  high-purity cases, 0.10–0.20 is appropriate

The non-high-confidence variants (`.Somatic.vcf` without `.hc`) are still
available for review but require additional validation.

---

### `MAX_NORMAL_FREQ = 0.05`

**What it does:** Maximum allowable VAF in the normal sample for a somatic call
to be classified high-confidence.

**Justification:** A true somatic mutation should be absent from the germline
normal. A VAF of up to 5% in the normal is tolerated rather than requiring 0%
for two reasons:

1. **Tumour contamination in normal:** If the normal sample is adjacent tissue
   rather than blood, low-level tumour DNA may be present
2. **Sequencing noise floor:** At 20–30× depth, 1–3 miscalled reads can produce
   apparent VAFs of 5–10%

Setting this too low (e.g. 0.01) causes false negatives when there is minimal
cross-contamination or sequencing noise in the normal. Setting it too high
(e.g. 0.15) risks classifying rare germline variants or pre-existing mosaic
variants as somatic.

---

### `PROCESS_P_VALUE = 0.05`

**What it does:** p-value threshold applied by `processSomatic` when classifying
variants into high-confidence categories.

**Justification:** This operates on the already-called variants from Stage 3 and
uses a Fisher's exact test comparing normal and tumour allele counts. At 0.05,
only variants where the normal-tumour difference is statistically significant at
the 5% level receive the `.hc` designation. This directly gates all downstream
analysis (VAR files, bam-readcount, fpfilter), making it the primary
sensitivity/specificity control for the final output.

Setting this above 0.05 admits marginally significant calls into high-confidence
output — increasing the false positive burden on bam-readcount and fpfilter
unnecessarily.

---

## 3. CNV Parameters

---

### `CNV_MIN_COVERAGE = 20`

**What it does:** Minimum read depth in both normal and tumour at a genomic
window for it to be included in copy number analysis.

**Justification:** Low-coverage windows produce unreliable log2 ratios due to
Poisson sampling variance. At 10× depth, the coefficient of variation for read
count is ~31%, producing log2 ratio noise of ±0.4 — large enough to mimic copy
number gains or losses. At 20×, noise drops to ±0.22 log2, which is below the
threshold for a single-copy gain (log2 = +0.58). Excluding under-covered windows
prevents noisy centromeric, repetitive, and low-complexity regions from
corrupting segmentation.

---

### `CNV_P_VALUE = 0.01`

**What it does:** Significance threshold for calling a window as having a
statistically significant copy number difference from the expected ratio.

**Justification:** 0.01 is the **VarScan2 official default** (`--p-value [0.01]`
per the copynumber command documentation). It is more stringent than the somatic
calling p-value (0.05) because copy number analysis operates on thousands of
genomic windows genome-wide. A 0.05 threshold with ~50,000 analysed windows
would yield ~2,500 false positive windows by chance alone. At 0.01, this drops
to ~500 — still requiring segmentation to merge adjacent windows, but
substantially cleaner input for `copyCaller`.

Note: this is the **window-level** calling threshold. VarScan also uses a
hardcoded internal `p < 0.05` for segment change-point detection — these are
two separate p-values operating at different levels of the CNV algorithm.

---

### `MIN_SEGMENT_SIZE = 10`

**What it does:** Minimum number of consecutive windows that must support a CNV
call for a segment to be reported.

**Justification:** A single outlier window — from a mapping artefact, a
polymorphic deletion, or a read-depth spike near a repetitive element — would
otherwise generate spurious single-window CNV calls. Requiring at least 10
consecutive windows ensures that reported segments span a meaningful genomic
distance (typically 10–50 kb depending on window size), consistent with real CNV
events. Single-window calls almost invariably represent noise.

---

### `MAX_SEGMENT_SIZE = 100`

**What it does:** Maximum number of windows allowed in a single CNV segment.

**Justification:** Very large uniform segments (>100 windows) may represent
baseline normalisation artefacts — for example, a chromosome arm that is
uniformly offset due to aneuploidy misclassified as a local CNV. Capping segment
size forces such events to be broken into multiple segments, which are then more
accurately handled by downstream circular binary segmentation (CBS). In practice,
most focal amplifications and deletions of biological interest span 10–50
windows; whole-arm events are better detected at the segmentation level.

---

## 4. copyCaller Parameters

---

### `CNV_AMP_THRESHOLD = 0.25` and `CNV_DEL_THRESHOLD = 0.25`

**What they do:** Log2 ratio thresholds used by `copyCaller` to classify
segments as amplifications or deletions:
- A segment with log2 ratio ≥ `CNV_AMP_THRESHOLD` is called an amplification
- A segment with log2 ratio ≤ −`CNV_DEL_THRESHOLD` is called a deletion

**Justification:** 0.25 is the **VarScan2 official default** for both
parameters. In a fully diploid, 100% pure tumour, a single-copy gain (3 copies)
produces a log2 ratio of log2(3/2) ≈ +0.58, and a single-copy loss (1 copy)
produces log2(1/2) = −1.0. The 0.25 threshold is intentionally set well below
the single-copy gain signal to account for:

- Tumour impurity (stromal/immune cell dilution reduces the observed log2 ratio)
- Intra-tumour heterogeneity (subclonal CNVs produce attenuated signals)
- GC content and mappability biases that offset individual windows

For high-purity samples (>90%), raising both thresholds to 0.30–0.40 reduces
noise calls. For low-purity samples (<40%), lower to 0.15–0.20 to detect
diluted CNV signals.

---

### `CNV_RECENTER_UP = 0` and `CNV_RECENTER_DOWN = 0`

**What they do:** Shift the log2 ratio baseline before classification. If set,
`copyCaller` adds `CNV_RECENTER_UP` to or subtracts `CNV_RECENTER_DOWN` from
all log2 ratios before applying the amp/del thresholds.

**Justification:** VarScan `copynumber` assumes the diploid baseline is at
log2 = 0. This assumption breaks in tumours with widespread chromosomal
instability (e.g. high-grade serous ovarian cancer, osteosarcoma, small cell
lung cancer), where the majority of the genome is aneuploid. In these cases, the
median log2 ratio across all windows is systematically offset from 0.

**How to determine the correct recenter value:**

1. Run Stage 5 (`copynumber`) and Stage 6 (`copyCaller`) with the defaults
   (`RECENTER_UP=0`, `RECENTER_DOWN=0`)
2. Inspect the `*.copynumber` output — calculate the median log2 ratio
   across all autosomal windows with sufficient coverage
3. If the median is negative (e.g. −0.3), set `CNV_RECENTER_DOWN=0.3` and
   re-run Stage 6 only
4. If the median is positive (e.g. +0.2), set `CNV_RECENTER_UP=0.2` and
   re-run Stage 6 only

Failure to recenter in an unstable genome will cause the majority of gains and
losses to be misclassified — the most common source of systematic error in
VarScan CNV analysis.

---

## 5. Data-Ratio Calculation

---

### `data-ratio = normal_primary_mapped / tumour_primary_mapped` with `scale=6`

**What it does:** Normalises the VarScan copy number log2 ratios for global
depth differences between the normal and tumour libraries. Without this
correction, a tumour sequenced at 120× and a normal at 60× would show a
genome-wide apparent +1 copy gain that is purely a depth artefact.

**Primary mapped reads:** The ratio is calculated from `primary mapped` counts
(as reported in `samtools flagstat`) rather than total mapped counts. This
excludes secondary alignments and supplementary (split/chimeric) reads, which
are disproportionately elevated in chromosomally unstable CRC tumours due to
structural rearrangements. Using total mapped counts can offset the ratio by
2–8% in samples with high structural rearrangement burden, introducing a
systematic CNV baseline error.

**Justification for `scale=6`:** The data-ratio enters directly into the log2
ratio computation:

```
log2(tumour_depth / normal_depth / data-ratio)
```

Rounding to 2 decimal places introduces a quantisation error of up to 0.005 in
the ratio. This propagates to a log2 error of:

```
log2(1.005) ≈ 0.007 log2 units
```

While small per-window, this systematic offset shifts the baseline of the entire
copy number profile, causing CBS segmentation to misidentify the diploid
baseline. Six decimal places reduces this quantisation error to a negligible
~7 × 10⁻⁸ log2 units.

---

## 6. bam-readcount Parameters

---

### `BRC_MAP_QUAL = 10`

**What it does:** Minimum mapping quality for reads included in the bam-readcount
allele counts. Also applied as `-q BRC_MAP_QUAL` in `samtools mpileup` (Stage 2)
and as `--min-map-qual` in `VarScan somatic` and `VarScan copynumber`.

**Justification:** 10 is the threshold explicitly recommended in the **VarScan2
FAQ**:
> "A mapping quality of 0 indicates multiple equally probable mapping locations.
> Exclude these reads. **A minimum mapping quality of 10 is even better.**"

`MAPQ` scores 1–9 indicate the aligner has low but non-zero confidence in the
read's placement. These reads originate disproportionately from:
- Segmental duplications (e.g. pericentromeric regions, chromosome 22q11)
- Paralogous gene families (CYP2D6, SMN1, NOTCH2NL, PMS2) where somatic
  events are clinically significant
- Repeat-flanking boundaries where structural variants occur

Reads from these regions produce phantom allele counts that inflate VAF
estimates and generate false positive somatic calls. The same threshold is
applied consistently to mpileup, VarScan internal filtering, and bam-readcount
to ensure all stages operate on the same read set.

---

### `BRC_BASE_QUAL = 15`

**What it does:** Minimum base quality score at the variant position for a read
to contribute to the bam-readcount allele counts (`-b` parameter).

**Justification:** Base quality 15 corresponds to a per-base error probability
of ~3.2% (Phred: `10^(-15/10)`). The VarScan2 fpfilter documentation
explicitly recommends `-b 15` as the appropriate floor for this filter. This is
lower than the `-b 20` often used in variant calling (error probability 1%)
because:

- fpfilter uses the base quality distribution **pattern** (not just the presence
  of reads) to assess variant quality
- Being too stringent here (e.g. `-b 20`) discards reads near read ends where
  base qualities legitimately drop, reducing coverage at variant sites
  unnecessarily
- fpfilter then applies `--min-ref-basequal` and `--min-var-basequal` thresholds
  using the full per-base quality data from bam-readcount

Using `-b 20` is inconsistent with the fpfilter workflow and suppresses readcount
evidence at low-base-quality variant-containing positions.

---

## 7. fpfilter Parameters

---

### `--min-var-count = 3`

**What it does:** Minimum number of reads supporting the variant allele for a
call to pass fpfilter.

**Justification:** Even at a site with high depth, a variant supported by only
1–2 reads is likely a sequencing error or PCR artefact. Three reads is the
minimum evidence threshold recommended in the VarScan2 fpfilter documentation.
This is an absolute count floor independent of VAF — a variant at 15% VAF with
only 2 supporting reads at 13× depth should not be trusted.

---

### `--min-var-freq = MIN_VAR_FREQ (0.10)`

**What it does:** Re-applies the VAF filter at the fpfilter stage using the more
accurate bam-readcount allele counts rather than VarScan's internal pileup counts.

**Justification:** VarScan's internal VAF calculation and bam-readcount can
produce slightly different counts due to differences in read filtering.
Re-applying the same 10% VAF floor using the bam-readcount data ensures
consistency and catches rare cases where VarScan called a variant at slightly
above 10% VAF but the more carefully filtered bam-readcount data shows it below
threshold.

---

### `--min-ref-basequal` and `--min-var-basequal = BRC_BASE_QUAL (15)`

**What they do:** Minimum base quality for reference-supporting and
variant-supporting reads respectively when fpfilter evaluates the quality
distribution.

**Justification:** These mirror the `-b 15` bam-readcount threshold for internal
consistency. fpfilter assesses the quality score distribution of ref vs alt
supporting reads — if variant-supporting reads have systematically lower base
quality than reference-supporting reads at the same position, it is a hallmark
of a sequencing artefact. Setting both thresholds to 15 ensures that the quality
distributions being compared are computed on the same quality-filtered read sets.

---

## Summary Table

| Parameter | Value | Official Default | Primary Purpose | Risk if too low | Risk if too high |
|---|---|---|---|---|---|
| `MIN_COVERAGE` | 20 | 8 | Statistical reliability | FP from sampling noise | Missed calls in low-depth regions |
| `MIN_COVERAGE_NORMAL` | 10 | 8 | Germline assessment accuracy | Germline misclassified as somatic | Miss calls with uneven depth |
| `MIN_COVERAGE_TUMOR` | 20 | 6 | Somatic detection floor (WES: ≥2 reads at 10% VAF) | Unreliable VAF estimates; single-read FP | Miss subclonal variants in low-coverage exons |
| `MIN_BASE_QUAL` | 20 | 20 | Base quality floor (mpileup, somatic, CNV) | Low-quality base errors inflate allele counts | Coverage loss near read ends |
| `MIN_VAR_FREQ` | 0.10 | 0.10 | Sequencing error floor | Error-rate FP | Miss low-VAF subclonal mutations |
| `MIN_FREQ_FOR_HOM` | 0.75 | 0.75 | Zygosity accuracy | Het variants called hom | Miss true hom variants in impure tumours |
| `P_VALUE` | 0.99 | 0.99 | Detection gate (permissive; enables p-value calculation) | Miss true variants before classification | Setting to 1.0 disables p-value calculation entirely |
| `SOMATIC_P_VALUE` | 0.05 | 0.05 | Somatic classification | Germline called somatic | Miss somatic variants with borderline support |
| `STRAND_FILTER` | 1 | 1 | Strand-bias artefact removal | FFPE / single-strand artefacts pass | Reduced sensitivity in amplicon data |
| `PROCESS_P_VALUE` | 0.05 | — | High-confidence gate | More FP in `.hc.vcf` | Reduced sensitivity for weak somatic signals |
| `MIN_TUMOR_FREQ` | 0.10 | — | HC somatic threshold | Noise in `.hc.vcf` output | Miss subclonal high-confidence calls |
| `MAX_NORMAL_FREQ` | 0.05 | — | Germline exclusion | Germline variants in somatic output | Miss somatic calls in contaminated normals |
| `CNV_MIN_COVERAGE` | 20 | 20 | CNV window reliability | Noisy log2 ratios in low-coverage windows | Miss CNVs in under-sequenced regions |
| `CNV_P_VALUE` | 0.01 | 0.01 | CNV window significance | Noise windows called as CNV | Miss low-amplitude focal CNVs |
| `MIN_SEGMENT_SIZE` | 10 | 10 | Artefact rejection | Single-window artefacts as CNV calls | Miss very focal amplifications |
| `MAX_SEGMENT_SIZE` | 100 | 100 | Baseline artefact control | Arm-level events split excessively | Large artefacts masked as single segment |
| `CNV_AMP_THRESHOLD` | 0.25 | 0.25 | Amplification classification boundary | Noise segments called as amplifications | Miss low-purity or subclonal amplifications |
| `CNV_DEL_THRESHOLD` | 0.25 | 0.25 | Deletion classification boundary | Noise segments called as deletions | Miss low-purity or subclonal deletions |
| `CNV_RECENTER_UP` | 0 | 0 | Baseline correction (globally positive profile) | Gains misclassified as diploid | Over-correct; introduce false deletions |
| `CNV_RECENTER_DOWN` | 0 | 0 | Baseline correction (globally negative profile) | Losses misclassified as diploid | Over-correct; introduce false gains |
| `BRC_MAP_QUAL` | 10 | — | Ambiguous-read exclusion (FAQ: ≥10) | Mismapped reads inflate allele counts | Reduced coverage at low-complexity sites |
| `BRC_BASE_QUAL` | 15 | — | Low-quality base exclusion (fpfilter recommendation) | Error reads inflate readcounts | Coverage loss at read-end variant positions |

---

## Usage

Print the built-in help and quick-reference:

```bash
./target/release/varscan2_pipeline --help
```

Run the full pipeline:

```bash
# 1. Edit GENOMEIDX1 in varscan2_pipeline.rs, then: cargo build --release
# 2. Place BAM files in the working directory
# 3. Create sample_pairs.csv (tumour col 1, normal col 2)
# 4. Place VarScan.v2.3.9.jar in software/ and fpfilter.pl in scripts/
./target/release/varscan2_pipeline
```

Run a single stage:

```bash
./target/release/varscan2_pipeline --stage 5
```

Run a range of stages:

```bash
./target/release/varscan2_pipeline --from 3 --to 6
```

Resume after a failed or partial run (skip stages with matching SHA256 state):

```bash
./target/release/varscan2_pipeline --resume
./target/release/varscan2_pipeline --from 5 --resume
```

Preview what would run without executing anything:

```bash
./target/release/varscan2_pipeline --dry-run
./target/release/varscan2_pipeline --from 3 --to 8 --dry-run
```

Capture a full log while running in the background:

```bash
./target/release/varscan2_pipeline 2>&1 | tee varscan_run.log &
```

---

## Output Files

### Somatic variants — `somatic/`

| Pattern | Description |
|---------|-------------|
| `*.snp.vcf` | All raw SNP calls before classification |
| `*.indel.vcf` | All raw INDEL calls before classification |
| `*.snp.Somatic.hc.vcf` | High-confidence somatic SNPs |
| `*.snp.Germline.hc.vcf` | High-confidence germline SNPs |
| `*.snp.LOH.hc.vcf` | High-confidence loss-of-heterozygosity SNPs |
| `*.indel.Somatic.hc.vcf` | High-confidence somatic INDELs |
| `*.snp.Somatic.vcf` | All somatic SNPs (includes low-confidence) |

### Copy number — `copynumber/`

| Pattern | Description |
|---------|-------------|
| `*.copynumber` | Per-window raw log2 read-depth ratios |
| `*.copynumber.called` | Segmented CNV calls with amplification/deletion status |
| `*.copynumber.homdel` | Homozygous deletion regions |

### Filtered variants — `filtered/`

| Pattern | Description |
|---------|-------------|
| `*.fpfilter.vcf` | Final variants after false positive removal |

---

## Key Design Decisions

### Normal-first mpileup ordering (Stage 2)

VarScan `somatic` uses column order in the paired mpileup to assign variant
classification. Column 1 must be the normal sample; column 2 must be the tumour.
Reversing the order silently swaps all Somatic/Germline/LOH labels.

```bash
samtools mpileup ... "$normal_bam" "$tumor_bam" > out.mpileup
#                     ^-- col 1       ^-- col 2
```

### processSomatic must precede bam-readcount (Stage 4 before Stage 7)

`processSomatic` generates the `.hc.vcf` files that Stage 7 converts to `.var`
position lists. Without Stage 4 completing first, Stages 7, 8, and 9 have no
input and produce no output.

### BAM selection in bam-readcount (Stage 8)

| Variant class | BAM used | Reason |
|---|---|---|
| Somatic SNP / INDEL | tumour | Somatic mutations are absent from normal constitutional DNA |
| Germline SNP | normal | Germline variants are present in all cells; the normal BAM gives the cleanest signal |
| LOH SNP | tumour | Loss of heterozygosity is a tumour-specific allele loss event |

### Data-ratio precision (Stage 5)

`--data-ratio` normalises the CNV log2 baseline for coverage depth differences
between normal and tumour:

```
data-ratio = normal_mapped_reads / tumour_mapped_reads
```

`scale=6` in `bc` retains six decimal places. Rounding to two decimal places
introduces a systematic offset in the CNV log-ratio baseline.

### Consistent quality filters across mpileup, VarScan, and bam-readcount

`MIN_BASE_QUAL` (20) and `BRC_MAP_QUAL` (10) are applied uniformly via:
- `samtools mpileup -Q 20 -q 10` — excludes bases/reads at the pileup level
- `VarScan somatic --min-base-qual 20 --min-map-qual 10` — internal filter
- `VarScan copynumber --min-base-qual 20 --min-map-qual 10` — internal filter
- `bam-readcount -b 15 -q 10` — readcount filter (BRC_BASE_QUAL=15 per fpfilter docs)

This ensures every stage operates on a consistent read set. Without this
alignment, VarScan may report coverage figures that differ from IGV or other
tools, and bam-readcount may count reads that VarScan excluded or vice versa.

### Stage 4 output guard

After `processSomatic`, the pipeline calls `check_stage4_output()` which
counts `.hc.vcf` files in `somatic/` and aborts with a diagnostic message if
none are found. This prevents Stages 7–9 from silently producing empty outputs
when processSomatic filtering is too strict or Stage 3 produced no calls.

### copyCaller baseline recentering

In tumours with widespread chromosomal instability, the genome-wide log2 ratio
is offset from 0. Set `CNV_RECENTER_UP` or `CNV_RECENTER_DOWN` to the absolute
value of the median log2 ratio before re-running Stage 6. This is the most
common source of systematic error in VarScan CNV analysis and cannot be
corrected downstream.

### VCF sample labels via `VCF_SAMPLE_LIST`

Without `VCF_SAMPLE_LIST`, VarScan outputs VCF files with generic `NORMAL` and
`TUMOR` column headers. Setting this variable to a two-line file (normal sample
name, then tumour sample name) propagates correct identifiers into the VCF
`##SAMPLE` header and the genotype columns, preventing confusion when merging
multi-sample VCFs or importing into annotation tools.

### Flexible BAM / pairs-file naming via `PAIRS_SUFFIX` and `BAM_SUFFIX`

Sample names are derived by stripping `PAIRS_SUFFIX` from each pairs-file
entry. BAM files are then located as `$BAM_DIR/<sample>${BAM_SUFFIX}`. Setting
these two variables independently allows the same script to work with any
combination of pairs-file format and BAM naming convention.

### Duplicate read exclusion via `-F 0x400` (Stage 2)

`samtools mpileup` does not exclude duplicate-flagged reads by default. BAMs
produced by GATK/Picard `MarkDuplicates` retain flagged duplicates in the file
with the 0x400 bit set. Without explicit exclusion, these reads contribute to
depth counts and variant allele counts — inflating VAFs and undermining the
purpose of deduplication. The `-F 0x400` flag is applied to the mpileup command
to exclude all duplicate-flagged reads before any VarScan stage sees the data.

### WES target region restriction via `TARGET_BED` (Stage 2)

For whole exome sequencing, setting `TARGET_BED` to the capture kit BED file
passes it to `samtools mpileup` as `-l TARGET_BED`. This restricts the pileup
to captured exon regions only. Benefits:

- Mpileup file sizes drop from 50–200 GB (WGS) to ~500 MB (WES) per pair
- Off-target reads that generate noise in intronic and intergenic regions are excluded
- Alt contigs, unplaced scaffolds, and EBV/decoy sequences are naturally excluded
  without requiring a separate chromosome allowlist

For WGS, leave `TARGET_BED` empty and the full-genome pileup is generated.

### Primary mapped reads for CNV data-ratio (Stage 5)

The `data-ratio` for VarScan `copynumber` is computed from `primary mapped`
read counts (from `samtools flagstat`) rather than total mapped counts. Primary
mapped excludes secondary and supplementary alignments (split reads, chimeric
pairs). In CRC tumours with chromosomal instability, supplementary alignment
rates are elevated 2–8× above normal due to structural rearrangements. Using
total mapped counts biases the ratio and shifts the CNV log2 baseline
systematically across all windows.

### INDEL fpfilter scope limited to Somatic class (Stage 7)

Stage 7 converts `.hc.vcf` files to bam-readcount position lists. For INDELs,
only `*.indel.Somatic.hc.vcf` files are processed — Germline and LOH INDELs
are not sent through bam-readcount or fpfilter. Rationale:

- Somatic INDELs carry the highest false positive burden from local alignment
  artefacts near indel breakpoints; fpfilter is most impactful here
- Germline and LOH INDEL `.hc.vcf` files from Stage 4 are the final output for
  those classes and are suitable for direct downstream use

### BAM index validation before bam-readcount (Stage 8)

`bam-readcount` requires a `.bai` index alongside each BAM. Without one it
exits silently, producing an empty readcount file that cascades to empty
fpfilter output with no error message. The `check_bam_index()` function verifies
both `.bam.bai` and `.bai` naming conventions before launching any bam-readcount
job, aborting with an actionable error if an index is missing.

---

## Post-Pipeline Analysis

### Circular binary segmentation (CBS)

Apply CBS to `*.copynumber.called` files to produce smoothed, segment-level
copy number profiles:

```r
library(DNAcopy)
cn  <- read.table("sample.copynumber.called", header = TRUE, sep = "\t")
cna <- CNA(cn$log2ratio, cn$chrom, cn$position, data.type = "logratio")
seg <- segment(smooth.CNA(cna))
write.table(seg$output, "sample.segments.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```

### Variant annotation

```bash
vep --input_file  sample.snp.Somatic.hc.fpfilter.vcf \
    --output_file sample.snp.Somatic.hc.annotated.vcf \
    --vcf --cache --assembly GRCh38 --everything
```

---

## Troubleshooting

**All somatic/CNV stages skip with "Mpileup not found"**
Mpileup files are named `{normal}_{tumour}.mpileup`. Verify the column order
in your pairs file (normal first, tumour second) and confirm `PAIRS_SUFFIX`
correctly strips the extension to give the bare sample name.

**Pipeline aborts at Stage 4 guard: "processSomatic produced no .hc.vcf files"**
Three possible causes in order of likelihood:
1. Stage 3 produced no variant calls — check coverage and `MIN_COVERAGE` values
2. `processSomatic` filtered all variants — lower `MIN_TUMOR_FREQ` (try `0.05`)
   or raise `MAX_NORMAL_FREQ` (try `0.10`) for low-purity samples
3. Stage 3 was run without `--output-vcf 1` — processSomatic then outputs
   tab-separated `.hc` files instead of `.hc.vcf` files; re-run Stage 3

**processSomatic produces empty `.hc.vcf` files (no abort)**
Lower `MIN_TUMOR_FREQ` for low-purity tumours (e.g. `0.05`) or raise
`MAX_NORMAL_FREQ` if the normal sample has detectable contamination.

**CNV calls are entirely gains or entirely losses**
The log2 ratio baseline is offset — the tumour has widespread chromosomal
instability. Calculate the median log2 ratio from the `*.copynumber` output
and set `CNV_RECENTER_UP` or `CNV_RECENTER_DOWN` accordingly, then re-run
Stage 6 (copyCaller) only.

**bam-readcount produces no output**
Confirm that `.var` files in `snp-VAR/` are non-empty and that chromosome
names in the `.var` file match the reference naming convention (e.g. `chr1`
vs `1`). The pipeline validates BAI index presence before launching jobs
(Stage 8) — if this check was bypassed, verify that `.bai` files exist
alongside each BAM.

**Read counts in VarScan differ from IGV or samtools**
By default, IGV and samtools show all bases regardless of quality. VarScan
applies `MIN_BASE_QUAL=20` internally and mpileup applies `-Q 20`. To match
IGV counts, temporarily set `MIN_BASE_QUAL=0` and remove `-Q` from mpileup —
but do not use these settings for production analysis.

**Data-ratio defaults to 1.0 for all pairs**
The pipeline greps for `primary mapped` in flagstat output. This line was added
in samtools 1.13. For older versions, inspect a flagstats file manually and
change the grep pattern to match the available line (e.g. `mapped (`). Note
that older samtools `mapped (` includes secondary and supplementary alignments,
so the ratio will be a slight overestimate in structurally rearranged tumours.

**Mpileup files are very large / pipeline is slow**
For WES, set `TARGET_BED` to your capture kit BED file. Without it, mpileup
runs over the entire genome and produces files of 50–200 GB per pair. With
TARGET_BED set to a typical exome BED, this drops to ~500 MB per pair.

**bam-readcount aborts with "BAM index not found"**
Run `samtools index <sample>.bam` for each BAM file to generate the `.bai`
index. The pipeline checks for both `<file>.bam.bai` and `<file>.bai` naming
conventions. Ensure samtools is on your `$PATH` when indexing.

**VCF files show generic NORMAL/TUMOR column headers**
Set `VCF_SAMPLE_LIST` to a plain-text file containing the normal sample name
on line 1 and the tumour sample name on line 2. Re-run Stage 3 to regenerate
VCFs with correct column labels.

**fpfilter skips all readcount files**
The VCF filename is derived from the readcount filename by appending `.vcf`
and looking in `$SOMATICDIR`. Confirm the readcount base name (without
`.readcount`) exactly matches the VCF base name in `somatic/`.

---

## FAQ

**Why isn't there a `Cargo.toml` in older clones?**
The `Cargo.toml` was added after the initial commit. Pull the latest version
and run `cargo build --release`.

**Can I run this on macOS?**
Yes. The Rust binary uses the `sha2` crate for all SHA256 computation, so there
is no dependency on the Linux-only `sha256sum` utility. Install Rust via
`rustup`, then build normally.

**Do I need root / sudo?**
No. The pipeline only requires a writable working directory. Install samtools
and bam-readcount into a user-local prefix (e.g. `~/local/bin`) and add it to
`$PATH`.

---

## License

MIT — see [LICENSE](LICENSE).
