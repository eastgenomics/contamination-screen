# contamination_screen

Pairwise cross-sample contamination detector for tumour sequencing panel VCFs.

## Overview

This tool screens a cohort of annotated tumour VCF files for cross-sample
contamination by detecting shared variants at consistent VAF ratios. It is
designed for targeted haematological oncology panel sequencing data processed by
the **Uranus clinical pipeline** (CUH Bioinformatics), using VCFs output by the
**eggd_vep** stage.

## Expected input VCFs

Input VCFs must have been produced by the Uranus pipeline eggd_vep stage.
Specifically, they must have already undergone:

| Step | Command | Header tag |
|---|---|---|
| Variant calling | Sentieon TNhaplotyper2 (tumour-only) | `##SentieonCommandLine.TNhaplotyper2` |
| Soft-filtering | Sentieon TNfilter | `##SentieonCommandLine.TNfilter` |
| Normalisation | `bcftools norm -m -any --keep-sum AD` | `##bcftools_normCommand` |
| VEP annotation | Ensembl VEP (CSQ field) | `##INFO=<ID=CSQ,...>` |

**The normalisation step (`bcftools norm -m -any`) must already have been
applied.** This tool does not re-normalise.

The VCF INFO field must contain a `CSQ` tag in standard Ensembl VEP format,
with `SYMBOL` at index 1 (0-based) in the pipe-delimited format string.

---

## Design logic

### The contamination signal

When one sample contaminates another, variants from the source sample appear in
the recipient at a consistent fraction of the source VAF:

```
VAF_recipient ~ VAF_source x contamination_fraction
```

Cross-sample contamination primarily consists of the **source patient's germline
heterozygous SNPs** (~50% VAF in source) appearing at a diluted fraction in the
recipient. These variants span the full allele frequency spectrum (common to
rare) and include synonymous, intronic, and coding variants equally. They are
excellent contamination markers precisely because they are stable, high-depth,
and well-genotyped.

### Why no population frequency filter

Traditional somatic reporting filters (gnomAD AF, cohort prevalence, synonymous
exclusion) are **counterproductive** for contamination detection:

- gnomAD AF >= 0.002 removed 87% of contamination markers in testing
- Synonymous filter removed 47% of contamination markers
- These filters exist to identify somatic mutations for clinical reporting;
  contamination detection is a fundamentally different task

The tool applies only quality filters (depth and VAF floor) that ensure reliable
allele fraction estimates.

### Dual-peak detection

With no population filter, shared common germline variants between two unrelated
people will cluster at **ratio = 1** (both carry them at ~50% VAF). The
contamination signal clusters at a **different, consistent ratio** (e.g. 0.5 for
50% contamination).

The tool reports **two peaks** from the log2-ratio histogram:

1. **Overall peak** (tallest bin across full range) — typically ratio=1 from
   shared germline. Useful for detecting sample swaps (many variants at
   identical VAF).

2. **Non-unity peak** (tallest bin where |log2_ratio| > 0.3, i.e. ratio outside
   0.81-1.23) — this is the **contamination signal**. Flagging and
   directionality are based on this peak.

This separation ensures that shared germline variants (ratio=1) cannot mask the
contamination signal, regardless of how many common variants two people share.

### Directionality

From the non-unity peak:

- `peak_log2_ratio < 0`: sample A's variants appear diluted in B — **A is the
  source**
- `peak_log2_ratio > 0`: sample B's variants appear diluted in A — **B is the
  source**

The contamination fraction estimate is `2^|peak_log2_ratio|`.

### Pre-filtering

Only quality filters are applied — no biological/annotation filters:

| Filter | Reason |
|---|---|
| PASS only (default) | Removes Sentieon TNfilter soft-failures (strand bias, weak evidence, etc.) |
| FORMAT/DP >= 99 | Ensures reliable VAF estimation from adequate read depth |
| AF >= 0.03 | Below this VAF, allele fractions are noisy and unreliable |

The pipeline is just two bcftools commands:
```bash
bcftools view -f PASS input.vcf.gz -Ou \
| bcftools filter -e '(FORMAT/DP<99 || AF<0.03)' -Oz -o filtered.vcf.gz
```

### Flagging thresholds

A pair is flagged if **both**:

- `n_informative >= --min-shared` (default 10)

**and at least one of**:

- `peak_count >= --peak-count` (default 8) — the non-unity peak has >= 8 variants
- `peak_fraction >= --peak-fraction` (default 0.30) — >= 30% of informative
  shared variants fall in the non-unity peak

---

## Requirements

- Python >= 3.9
- `bcftools` >= 1.14
- `pandas` and `numpy`
- `matplotlib` (optional, for `--plots`)

Input VCFs must be bgzipped (`.vcf.gz`) and tabix-indexed (`.vcf.gz.tbi`).

## Installation

```bash
git clone <repo>
cd contamination
pip install pandas numpy matplotlib
```

## Usage

```
contamination_screen.py VCF_DIR [options]
```

### Key options

| Option | Default | Description |
|---|---|---|
| `--outdir / -o` | `results/` | Output directory |
| `--min-af` | `0.03` | VAF floor for informative variants |
| `--min-dp` | `99` | Minimum read depth |
| `--min-shared` | `10` | Minimum informative shared variants to assess |
| `--peak-count` | `8` | Flag if non-unity peak has >= N variants |
| `--peak-fraction` | `0.30` | Flag if non-unity peak fraction >= this |
| `--threads / -t` | `min(8, nCPU)` | Parallel threads |
| `--include-non-pass` | off | Include non-PASS variants |
| `--plots` | off | Generate histogram plots for flagged pairs |
| `--force-refilter` | off | Regenerate filtered VCFs |
| `--bin-width` | `0.2` | Histogram bin width (log2 units) |
| `--vcf-glob` | `*_annotated.vcf.gz` | File matching pattern |
| `--verbose / -v` | off | Debug logging |

### Examples

```bash
# Standard run
python contamination_screen.py /data/vcfs/ --outdir results/ --plots

# More sensitive (smaller panels)
python contamination_screen.py /data/vcfs/ \
    --min-shared 5 --peak-count 4 --outdir results_sensitive/
```

---

## Output files

```
results/
├── filtered/                          Quality-filtered VCFs (reused on re-runs)
├── summary.tsv                        One row per pair (dual-peak results)
├── matrix.tsv                         N x N directional contamination matrix
├── flagged_pairs/*.tsv                Variant-level detail for flagged pairs
└── plots/*.png                        Log2-ratio histograms (--plots only)
```

### `summary.tsv` columns

| Column | Description |
|---|---|
| `sample_a`, `sample_b` | Sample pair |
| `n_shared` | Total shared variants |
| `n_informative` | Shared variants with VAF >= min-af in both |
| `overall_log2`, `overall_count`, `overall_fraction` | Tallest peak across full range (usually ratio=1 from germline) |
| `peak_log2_ratio`, `peak_ratio` | Non-unity peak centre (the contamination signal) |
| `peak_count`, `peak_fraction` | Strength of non-unity peak |
| `contamination_source` | Source sample (higher VAF) |
| `contamination_recipient` | Recipient sample (diluted VAF) |
| `contamination_fraction` | Estimated fraction of recipient from source |
| `flagged` | TRUE if meets thresholds |

---

## Interpretation

- **Non-unity peak with high count:** cross-sample contamination at the implied
  fraction. Direction indicates which sample is the source.
- **Large overall peak at ratio=1:** many variants shared at similar VAF. In
  unrelated samples this suggests a sample swap or duplicate; in related samples
  (same patient) it reflects shared biology.
- **No significant peaks:** clean pair with minimal variant sharing.
- The tool is designed for **unrelated samples**. Same-patient pairs will have
  large overall peaks at ratio=1 (shared clonal mutations) which is expected and
  not flagged as contamination.
