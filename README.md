# contamination_screen

Pairwise cross-sample contamination detector for tumour sequencing panel VCFs.

## Overview

This tool screens a cohort of annotated tumour VCF files for cross-sample
contamination by detecting shared rare somatic variants at consistent VAF ratios.
It is designed for targeted haematological oncology panel sequencing data
processed by the **Uranus clinical pipeline** (CUH Bioinformatics). The
pre-filtering steps exactly replicate the Uranus filtering logic so that the
variant set used for contamination screening is identical to the variant set
reported clinically.

## Expected input VCFs

Input VCFs must have been produced by the Uranus pipeline. Specifically, they
are expected to have already undergone the following steps, which are confirmed
by the presence of the corresponding command-line entries in the VCF header:

| Step | Command | Header tag |
|---|---|---|
| Variant calling | Sentieon TNhaplotyper2 (tumour-only) | `##SentieonCommandLine.TNhaplotyper2` |
| Soft-filtering | Sentieon TNfilter | `##SentieonCommandLine.TNfilter` |
| Normalisation | `bcftools norm -m -any --keep-sum AD` | `##bcftools_normCommand` |
| VEP annotation | Ensembl VEP with gnomAD, COSMIC, cohort prevalence, SpliceAI | `##VEP` |
| Cohort annotation | `bcftools annotate` (gnomAD, Prev_Count, RESCUE_LIST) | `##bcftools_annotateCommand` |

**The normalisation step (`bcftools norm -m -any`) must already have been
applied.** This tool does not re-normalise. Running it on un-normalised VCFs
will produce incorrect results because multiallelic records will not be split
and indel representations may differ between samples.

The VCF INFO field must contain a `CSQ` tag in standard Ensembl VEP format,
with the following subfields present at the indices expected by `bcftools
+split-vep`:

| Subfield | CSQ index (0-based) | Used for |
|---|---|---|
| `SYMBOL` | 1 | Gene label in output |
| `Consequence` | 2 | Synonymous variant filter |
| `gnomADe_AF` | 24 | Population frequency filter |
| `gnomADg_AF` | 27 | Population frequency filter |
| `Prev_Count_AC` | 31 | Cohort prevalence filter |

---

## Design logic

### The contamination signal

When one sample contaminates another (e.g. through index hopping, library
preparation cross-talk, or sample mix-up), variants from the source sample
appear in the recipient at a consistent fraction of the source VAF:

```
VAF_recipient ≈ VAF_source × contamination_fraction
```

If sample A (source) has a variant at 40% VAF and contaminates sample B at
50% level, that variant appears in B at ~20% VAF. If many variants from sample
A share this 2:1 ratio in sample B, that is almost certainly not random — it
is contamination.

This tool exploits that signature by computing `log2(VAF_B / VAF_A)` for every
shared variant in a pair and looking for a cluster of variants at a consistent
non-zero log2 ratio. The log2 scale is used because:

- It symmetrises the ratio (forward 2:1 = +1, reverse 1:2 = -1)
- Any contamination fraction produces a distinct sharp peak rather than a broad skew
- It is visually intuitive in histogram form

### Directionality

A peak at `log2_ratio < 0` means sample A's variants appear diluted in B —
**A is the contamination source**.  
A peak at `log2_ratio > 0` means B's variants appear diluted in A —
**B is the contamination source**.

The contamination fraction estimate is `2^|peak_log2_ratio|` (e.g. a peak at
-1 implies ~50% contamination of B by A).

### Central exclusion zone

Variants at similar VAF in both samples (log2 ratio ≈ 0) are excluded from
peak detection. These may represent variants shared for biological reasons (e.g.
same patient at different timepoints, clonal haematopoiesis present in both
samples) rather than contamination. The default exclusion window is ±0.3 log2
units (ratio range 0.81–1.23), configurable with `--central-excl`.

### Pre-filtering — Uranus clinical filter

Before comparison, each VCF is filtered to remove variants that add noise
rather than signal. The filtering exactly replicates the Uranus clinical
pipeline filter, applied in the same order and using the same commands:

| Excluded | Reason |
|---|---|
| gnomAD genome AF ≥ 0.002 | Common germline variants present in all samples |
| gnomAD exome AF ≥ 0.002 | As above |
| Cohort `Prev_Count_AC` > 853 | Highly recurrent in the Uranus cohort (artefacts or ubiquitous CH mutations) |
| Synonymous variants | Non-functional; unlikely to be specifically somatic |
| FORMAT/DP < 99 | Insufficient read depth for reliable VAF estimation |
| AF < 0.03 | Below minimum reportable VAF threshold |
| *Exception: GATA2 and TP53 synonymous variants are retained* | Clinically relevant markers in haematological malignancy |

The filtering pipeline is 6 steps, all piped with no intermediate files:

1. **`bcftools +split-vep --columns - -a CSQ -Ou -p 'CSQ_' -d`**  
   Exact Uranus split-vep command. Extracts all CSQ subfields into
   `CSQ_`-prefixed INFO tags (e.g. `CSQ_SYMBOL`, `CSQ_gnomADg_AF`),
   outputting one VCF record per transcript (`-d`/duplicate mode).  
   Using the `CSQ_` prefix avoids name conflicts with the existing
   `String`-typed standalone INFO tags (`gnomADg_AF`, `Prev_Count_AC`, etc.)
   that were added by the upstream annotation pipeline. The split-vep built-in
   `.*_AF` type rule automatically assigns `Float` to all gnomAD AF fields.

2. **`bcftools annotate -x INFO/CSQ`**  
   Exact Uranus command. Removes the now-redundant raw CSQ string, leaving
   only the expanded `CSQ_*` INFO tags.

3. **`bcftools annotate -x INFO/CSQ_Prev_Count_AC -h <hdr>`**  
   *Additional step not in Uranus — required for this tool only.*  
   split-vep assigns `String` to `CSQ_Prev_Count_AC` because the field name
   does not match any built-in Float/Integer type rule. A one-line header
   override recasts it to `Integer` so the `>853` arithmetic comparison in the
   next step works.

4. **`bcftools filter --soft-filter EXCLUDE -m + -e '(CSQ_Prev_Count_AC>853 || ...'`**  
   Exact Uranus filter expression. Soft-tags matching records by appending
   `EXCLUDE` to the FILTER column (`-m +` preserves existing FILTER values).

5. **`bcftools filter -e '(FORMAT/DP<99 || AF<0.03)'`**  
   Exact Uranus filter. Hard-removes low-depth and low-VAF variants. The
   `AF<0.03` threshold is consistent with `--min-af` (default 0.03).

6. **`bcftools view [-f PASS] -e 'FILTER~"EXCLUDE"'`**  
   Hard-removes all `EXCLUDE`-tagged records. In PASS-only mode (default),
   also restricts to variants that were originally PASS in the Sentieon
   TNfilter output. Pass `--include-non-pass` to retain soft-filtered variants.

Because split-vep `-d` creates one VCF record per transcript, the same variant
can appear multiple times. Most duplicates are removed by the downstream filters
(different transcripts have different `CSQ_Consequence` and `CSQ_SYMBOL`). Any
survivors are deduplicated when loading into memory, keeping the first
occurrence (VEP's primary/worst-consequence transcript).

Filtered VCFs are written to `results/filtered/` and reused on re-runs. Use
`--force-refilter` to regenerate them (e.g. after changing `--include-non-pass`).

### Pairwise comparison

For N samples, N×(N−1)/2 pairs are assessed. Variant DataFrames are loaded
into memory (typically < 1 MB per sample after filtering), so all comparisons
run from in-memory data with no temporary files. A multiprocessing pool is
used for cohorts with > 50 pairs; smaller cohorts run sequentially to avoid
subprocess pickle overhead.

### Flagging thresholds

A pair is flagged if **both**:

- `n_informative ≥ --min-shared` (default 10) — enough variants to be meaningful

**and at least one of**:

- `peak_count ≥ --peak-count` (default 8) — ≥ 8 variants cluster at the same ratio
- `peak_fraction ≥ --peak-fraction` (default 0.30) — ≥ 30% of informative shared variants fall in the same ratio bin

Thresholds should be tuned to cohort size and panel design. For small panels
(< 100 genes), `--min-shared` may need to be reduced.

---

## Requirements

- Python ≥ 3.9
- `bcftools` ≥ 1.14 with the `split-vep` plugin available
- `pandas` and `numpy` Python packages
- `matplotlib` (optional, for `--plots`)

Input VCFs must be bgzipped (`.vcf.gz`) and tabix-indexed (`.vcf.gz.tbi`).

## Installation

```bash
git clone <repo>
cd contamination
pip install pandas numpy matplotlib   # or: uv pip install pandas numpy matplotlib
```

## Usage

```
contamination_screen.py VCF_DIR [options]
```

### Required argument

| Argument | Description |
|---|---|
| `VCF_DIR` | Directory containing Uranus-annotated `*_annotated.vcf.gz` files |

### Key options

| Option | Default | Description |
|---|---|---|
| `--outdir / -o` | `results/` | Output directory |
| `--min-af` | `0.03` | VAF floor — variants below this in either sample are excluded from the informative set (consistent with Uranus AF < 0.03 filter) |
| `--min-shared` | `10` | Minimum informative shared variants required to assess a pair |
| `--peak-count` | `8` | Flag if dominant ratio bin contains ≥ this many variants |
| `--peak-fraction` | `0.30` | Flag if dominant ratio bin fraction ≥ this value |
| `--threads / -t` | `min(8, nCPU)` | Parallel threads (used when > 50 pairs) |
| `--include-non-pass` | off | Include soft-filtered (non-PASS) variants; by default only Sentieon PASS variants are used |
| `--plots` | off | Save log2-ratio histogram PNG for each flagged pair (requires matplotlib) |
| `--force-refilter` | off | Regenerate filtered VCFs even if they already exist |
| `--bin-width` | `0.2` | Histogram bin width in log2 units |
| `--central-excl` | `0.3` | Exclude ±this log2 window around 0 from peak detection |
| `--vcf-glob` | `*_annotated.vcf.gz` | Glob pattern to match VCF files in VCF_DIR |
| `--verbose / -v` | off | Debug logging |

### Examples

```bash
# Standard run: PASS-only variants, default thresholds
python contamination_screen.py /data/vcfs/ --outdir results/

# Include soft-filtered variants, generate histogram plots
python contamination_screen.py /data/vcfs/ \
    --include-non-pass \
    --plots \
    --outdir results_nonpass/

# More sensitive settings for small panels (< 100 genes)
python contamination_screen.py /data/vcfs/ \
    --min-shared 5 \
    --peak-count 4 \
    --peak-fraction 0.20 \
    --outdir results_sensitive/

# Force re-filter after changing --include-non-pass
python contamination_screen.py /data/vcfs/ \
    --force-refilter \
    --outdir results/
```

---

## Output files

```
results/
├── filtered/
│   ├── sample1_annotated.vcf.gz        Uranus-filtered VCF (reused on re-runs)
│   ├── sample1_annotated.vcf.gz.tbi
│   └── ...
├── summary.tsv                         One row per pair
├── matrix.tsv                          N×N directional contamination matrix
├── flagged_pairs/
│   └── SAMPLE_A__vs__SAMPLE_B.tsv      Shared variant detail for flagged pairs
└── plots/
    └── SAMPLE_A__vs__SAMPLE_B.png      Log2-ratio histogram (--plots only)
```

### `summary.tsv` columns

| Column | Description |
|---|---|
| `sample_a`, `sample_b` | Sample pair (in input file order) |
| `n_shared` | Total variants present in both samples (before VAF floor) |
| `n_informative` | Shared variants with VAF ≥ `--min-af` in both samples |
| `peak_log2_ratio` | Centre of the dominant ratio bin in log2 units |
| `peak_ratio` | 2^peak_log2_ratio — the implied VAF ratio between samples |
| `peak_count` | Number of variants in the dominant ratio bin |
| `peak_fraction` | peak_count / n_informative |
| `contamination_source` | Sample whose variants appear at higher VAF (the source) |
| `contamination_recipient` | Sample whose variants appear at lower VAF (the recipient) |
| `contamination_fraction` | Estimated fraction of recipient library derived from source |
| `flagged` | TRUE if pair meets the flagging thresholds |

### `matrix.tsv`

An N×N directional matrix. Rows = contamination recipient, columns =
contamination source. Cell value is `peak_count` for the source→recipient
direction. A clean cohort produces an all-zero matrix; a contaminated pair
produces a non-zero entry at `[recipient, source]`.

### Flagged pair detail TSVs

One file per flagged pair in `flagged_pairs/`. Contains all informative shared
variants with columns:

`chrom`, `pos`, `ref`, `alt`, `filter_a`, `filter_b`, `af_a`, `af_b`, `ratio`,
`log2_ratio`, `in_peak`, `gene`

Sorted by `log2_ratio` so the contamination peak cluster is visually
contiguous. `in_peak = True` marks the variants that fall within the dominant
ratio bin and drive the flagging.

---

## Interpretation notes

- A **sharp peak away from 0** in the log2-ratio histogram is the signature of
  contamination. The peak position gives the contamination fraction (2^|peak|)
  and its sign gives the direction.
- A **diffuse distribution centred on 0** is expected for unrelated samples
  with few shared rare variants — this is the null pattern.
- A peak **near 0** (excluded from detection by default) indicates variants
  shared at similar VAF in both samples. This is expected for serial samples
  from the same patient (clonal haematopoiesis or somatic mutations present at
  both timepoints) and does **not** indicate contamination.
- Very **low n_informative** (< `--min-shared`) pairs are reported in
  `summary.tsv` but not flagged; there are insufficient shared rare variants to
  draw conclusions.
- The `contamination_fraction` estimate assumes the source variant VAF
  represents the true allele fraction in the source library. For near-clonal
  or near-homozygous variants the estimate may be less accurate.
- This tool is designed for **unrelated samples** in a cohort. Pairs from the
  same patient at different timepoints will share many variants near ratio 1
  (which is not flagged) but may also show contamination signal if there was
  genuine cross-contamination.
