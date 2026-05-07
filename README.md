# contamination_screen

Pairwise cross-sample contamination detector for tumour sequencing panel VCFs.

## Overview

This tool screens a cohort of annotated tumour VCF files for cross-sample contamination by detecting shared rare somatic variants at consistent VAF ratios. It is designed for targeted panel sequencing data annotated with VEP/Ensembl CSQ fields, gnomAD population frequencies, and a cohort prevalence count.

## Design logic

### The contamination signal

When one sample contaminates another (e.g. through index hopping, library preparation cross-talk, or sample mix-up), variants from the source sample appear in the recipient at a consistent fraction of the source VAF:

```
VAF_recipient ≈ VAF_source × contamination_fraction
```

So if sample A (source) has a variant at 40% VAF and contaminates sample B at 50% level, that variant appears in B at ~20% VAF. If 30 variants from sample A share this 2:1 ratio in sample B, that is almost certainly not random — it is contamination.

This tool exploits that signature by computing `log2(VAF_B / VAF_A)` for every shared variant in a pair and looking for a cluster of variants at a consistent non-zero log2 ratio. The log2 scale is used because:
- It symmetrises the ratio (forward 2:1 = +1, reverse = -1)
- Any contamination fraction produces a distinct sharp peak rather than a broad skew
- Visually intuitive in histogram form

### Directionality

A peak at `log2_ratio < 0` means sample A's variants are diluted in B → **A is the source**.
A peak at `log2_ratio > 0` means B's variants are diluted in A → **B is the source**.

The contamination fraction estimate is `2^|peak_log2_ratio|` (e.g. a peak at −1 implies ~50% contamination of B by A).

### Central exclusion zone

Variants at similar VAF in both samples (log2 ratio ≈ 0) are excluded from peak detection. These may represent variants shared for biological reasons (same patient at different timepoints, true clonal overlap) rather than contamination. The default exclusion window is ±0.3 log2 units (ratio range 0.81–1.23), configurable with `--central-excl`.

### Pre-filtering

Before comparison, each VCF is filtered to remove variants that would add noise:

| Excluded | Reason |
|---|---|
| gnomAD genome AF ≥ 0.002 | Common germline variants appear in all samples |
| gnomAD exome AF ≥ 0.002 | As above |
| Cohort Prev_Count_AC > 853 | Recurrent artefacts or very common CH mutations |
| Synonymous variants | Non-functional, not specifically somatic |
| *Exception: GATA2 and TP53 synonymous variants are retained* | Clinically relevant markers |

This filtering requires a 4-step bcftools pipeline. `bcftools norm -m -any` has already been applied in the upstream clinical pipeline (confirmed via `bcftools_normCommand` in the VCF header), so no normalisation step is needed here. The pipeline:

1. `bcftools +split-vep --columns - -a CSQ -p CSQ_ -d` — the **exact clinical pipeline split-vep command**. Extracts all CSQ subfields into `CSQ_`-prefixed INFO tags, outputting one record per transcript (`-d`). The `CSQ_` prefix avoids conflicts with the existing String-typed INFO tags.
2. `bcftools annotate -x INFO/CSQ` — removes the now-redundant raw CSQ string.
3. `bcftools annotate -h` — recasts `CSQ_Prev_Count_AC` from `String` to `Integer` (split-vep's built-in type rules do not match `Prev_Count_AC`, so it defaults to String; arithmetic comparison `>853` requires Integer).
4. `bcftools filter --soft-filter EXCLUDE -m +` — the **exact expression from the clinical pipeline**, soft-tagging matching records with `FILTER=EXCLUDE`.
5. `bcftools filter -e '(FORMAT/DP<99 || AF<0.03)'` — hard-filters low-depth (DP < 99) and low-VAF (AF < 0.03) variants.
6. `bcftools view [-f PASS] -e 'FILTER~"EXCLUDE"'` — hard-filters by dropping `EXCLUDE`-tagged records; optionally restricts to originally-PASS records.

Because split-vep `-d` creates one VCF record per transcript, the same variant may appear multiple times. Most duplicates are removed by the downstream filters (different transcripts have different `CSQ_Consequence` and `CSQ_SYMBOL`). Any survivors are deduplicated in memory when loading, keeping the first occurrence (VEP's primary/worst-consequence transcript).

Filtered VCFs are written to `results/filtered/` and reused on re-runs. Use `--force-refilter` to regenerate them (e.g. after changing `--include-non-pass`).

### Pairwise comparison

For N samples, N×(N−1)/2 pairs are assessed. Variant DataFrames are loaded into memory (typically < 1 MB per sample after filtering), so all comparisons run from in-memory data with no temporary files. A multiprocessing pool is used for cohorts with > 50 pairs; smaller cohorts run sequentially to avoid pickle overhead.

### Flagging thresholds

A pair is flagged if **both**:
- `n_informative ≥ --min-shared` (default 10) — enough variants to be meaningful

**and at least one of**:
- `peak_count ≥ --peak-count` (default 8) — ≥ 8 variants cluster at the same ratio
- `peak_fraction ≥ --peak-fraction` (default 0.30) — ≥ 30% of informative shared variants at the same ratio

Thresholds should be tuned to cohort size and panel design. For small panels (< 100 genes), `--min-shared` may need to be reduced.

## Requirements

- Python ≥ 3.9
- `bcftools` ≥ 1.14 with `split-vep` plugin
- `pandas` and `numpy` Python packages
- `matplotlib` (optional, for `--plots`)

## Installation

```bash
git clone <repo>
cd contamination
pip install pandas numpy matplotlib   # or: uv pip install ...
```

## Usage

```
contamination_screen.py VCF_DIR [options]
```

### Required argument

| Argument | Description |
|---|---|
| `VCF_DIR` | Directory containing `*_annotated.vcf.gz` files |

### Key options

| Option | Default | Description |
|---|---|---|
| `--outdir / -o` | `results/` | Output directory |
| `--min-af` | `0.03` | VAF floor — variants below this in either sample are excluded from the informative set |
| `--min-shared` | `10` | Minimum informative shared variants to assess a pair |
| `--peak-count` | `8` | Flag if dominant ratio bin has ≥ this many variants |
| `--peak-fraction` | `0.30` | Flag if dominant ratio bin fraction ≥ this value |
| `--threads / -t` | `min(8, nCPU)` | Parallel threads (used when > 50 pairs) |
| `--include-non-pass` | off | Include soft-filtered (non-PASS) variants |
| `--plots` | off | Save log2-ratio histogram PNG for each flagged pair |
| `--force-refilter` | off | Regenerate filtered VCFs even if they exist |
| `--bin-width` | `0.2` | Histogram bin width in log2 units |
| `--central-excl` | `0.3` | Exclude ±this log2 window around 0 from peak detection |
| `--vcf-glob` | `*_annotated.vcf.gz` | Glob to match VCF files in VCF_DIR |
| `--verbose / -v` | off | Debug logging |

### Examples

```bash
# Standard run: PASS-only, default thresholds
python contamination_screen.py /data/vcfs/ --outdir results/

# Include soft-filtered variants, generate plots
python contamination_screen.py /data/vcfs/ \
    --include-non-pass \
    --plots \
    --outdir results_all/

# More sensitive for small panels
python contamination_screen.py /data/vcfs/ \
    --min-shared 5 \
    --peak-count 4 \
    --peak-fraction 0.20 \
    --outdir results_sensitive/

# Force re-filter (e.g. after switching --include-non-pass)
python contamination_screen.py /data/vcfs/ \
    --force-refilter \
    --outdir results/
```

## Output files

```
results/
├── filtered/
│   ├── sample1_annotated.vcf.gz        Pre-filtered VCF (reused on re-runs)
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
| `sample_a`, `sample_b` | Sample pair (alphabetical by input order) |
| `n_shared` | Total variants shared between the pair |
| `n_informative` | Shared variants with VAF ≥ `--min-af` in both |
| `peak_log2_ratio` | Centre of dominant ratio bin (log2 scale) |
| `peak_ratio` | 2^peak_log2_ratio — the implied VAF ratio |
| `peak_count` | Variants in the dominant ratio bin |
| `peak_fraction` | peak_count / n_informative |
| `contamination_source` | Sample contributing variants at higher VAF |
| `contamination_recipient` | Sample receiving contamination |
| `contamination_fraction` | Estimated fraction of recipient library from source |
| `flagged` | TRUE if pair meets flagging thresholds |

### `matrix.tsv`

An N×N matrix where rows are recipients and columns are sources. Cell value is the `peak_count` for the source→recipient direction. A clean cohort has all-zero cells; contaminated pairs show a non-zero entry in exactly one cell per pair (the source→recipient direction).

### Flagged pair detail TSVs

One file per flagged pair. Contains all informative shared variants with columns:
`chrom`, `pos`, `ref`, `alt`, `filter_a`, `filter_b`, `af_a`, `af_b`, `ratio`, `log2_ratio`, `in_peak`, `gene`

Sorted by `log2_ratio` so the peak cluster is visually contiguous. `in_peak = True` marks variants that fall within the dominant ratio bin.

## Interpretation notes

- A **sharp peak away from 0** in the log2-ratio histogram strongly suggests contamination. The peak position gives the contamination fraction and direction.
- A peak **near 0** (which is excluded from detection) may indicate biologically related samples (e.g. serial timepoints from the same patient) sharing somatic mutations at similar VAF — this is expected and is not contamination.
- Very **low n_informative** (< 10) pairs are reported but not flagged; insufficient shared rare variants to draw conclusions.
- The `contamination_fraction` estimate assumes the source variant VAF represents the true allele fraction in the source library. For very high-VAF or near-homozygous variants, the estimate may be less accurate.
