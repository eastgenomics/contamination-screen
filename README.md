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

### Why no population frequency filter (below 0.40)

Traditional somatic reporting filters (gnomAD AF >= 0.002, cohort prevalence,
synonymous exclusion) are **counterproductive** for contamination detection.
Empirical testing on known contaminated pairs showed:

| Clinical filter | Contamination markers removed | Impact |
|---|---|---|
| gnomAD AF >= 0.002 | 35/40 (87.5%) | Devastating |
| Synonymous exclusion | 19/40 (47.5%) | Very harmful |
| Prev_Count_AC > 853 | 6/40 (15%) | Harmful |

These filters exist to identify somatic mutations for clinical reporting.
Contamination detection is a fundamentally different task — the contamination
signal IS germline variants, which these filters are designed to remove.

### The gnomAD >= 0.40 threshold

The one population-frequency filter retained is **gnomAD AF >= 0.40**.

**The problem it solves:** At high population frequencies, a substantial
proportion of individuals are homozygous-alt (AF ~ 1.0) while others are
heterozygous (AF ~ 0.5). When these two genotypes appear in different samples,
the ratio is:

```
AF_hom / AF_het = 1.0 / 0.5 = 2.0  (log2 = 1.0)
```

This exactly mimics 50% contamination and produces false positives.

**Empirical validation:**
- All 15 false-positive "contamination" markers between uncontaminated pairs had
  gnomAD AF > 0.40 (range: 0.43–0.81)
- Only 1 of 18 genuine contamination markers had gnomAD AF >= 0.40

**At gnomAD AF = 0.40:**
- P(homozygous-alt) = 0.40² = 16% of individuals
- P(heterozygous) = 2 × 0.40 × 0.60 = 48% of individuals
- Pairing a hom-alt with a het is common enough to produce many false 2:1 ratios

**Below gnomAD AF = 0.30:**
- P(homozygous-alt) = 0.30² = 9%
- Hom/het pairings still occur but at lower frequency; and at gnomAD AF 0.26–0.30,
  linked variants in the same gene (e.g. 5 PRPF8 SNPs on the same haplotype) can
  cluster at a single ratio. However, these linked clusters typically produce
  peak_count = 4–5 which is below the flagging threshold of 6.

The threshold is configurable with `--max-gnomad` (set to 1.0 to disable).

### Dual-peak detection

With no population filter below 0.40, shared common germline variants between
two unrelated people will cluster at **ratio = 1** (both carry them at ~50%
VAF). The contamination signal clusters at a **different, consistent ratio**
(e.g. 0.5 for 50% contamination).

The tool reports **two peaks** from the log2-ratio histogram:

1. **Overall peak** (tallest bin across full range) — typically ratio=1 from
   shared germline. Useful for detecting sample swaps (many variants at
   identical VAF).

2. **Non-unity peak** (tallest bin where |log2_ratio| > 0.3, i.e. ratio outside
   0.81–1.23) — this is the **contamination signal**. Flagging and
   directionality are based on this peak.

This separation ensures that shared germline variants (ratio=1) cannot mask the
contamination signal, regardless of how many common variants two people share.

### The unity zone boundary (|log2| <= 0.3)

The boundary between the "unity" region (shared germline) and "non-unity" region
(potential contamination) is set at **|log2_ratio| = 0.3**, corresponding to a
ratio range of 0.81–1.23.

**Derivation from measurement noise:**

For a germline heterozygous variant (AF ≈ 0.5) measured at read depth D, the
noise on log2(ratio) between two independent measurements is:

```
SE(log2_ratio) ≈ sqrt(2) × sqrt(0.25/D) / (0.5 × ln2)
```

| Depth | SE(log2_ratio) | 3σ range |
|---|---|---|
| 99 (minimum) | 0.205 | ±0.615 |
| 200 | 0.144 | ±0.433 |
| 500 | 0.091 | ±0.274 |
| 1000 (typical) | 0.065 | ±0.194 |
| 2000 | 0.046 | ±0.137 |

At the minimum depth (DP=99), the 0.3 threshold captures ±1.5σ (87% of germline
noise). At typical panel depth (DP=1000), it captures ±4.6σ (essentially all
noise).

**Empirical validation:** In test data, germline variants in the unity zone had
SD(log2_ratio) = 0.09, with all values within ±0.25 — well contained by the 0.3
boundary.

**Impact on detection range:**
- Detects contamination fractions from 0% to 81% (ratios 0–0.81 and 1.23–∞)
- Contamination above 81% (ratio 0.81–1.23) falls inside the unity zone but
  manifests as a near-complete sample swap and is captured by the **overall
  peak** — this is exactly what the dual-peak design handles.

### Directionality

From the non-unity peak:

- `peak_log2_ratio < 0`: sample A's variants appear diluted in B — **A is the
  source**
- `peak_log2_ratio > 0`: sample B's variants appear diluted in A — **B is the
  source**

The contamination fraction estimate is `2^|peak_log2_ratio|`.

### Pre-filtering

The pipeline applies quality filters and a single, targeted population-frequency
filter:

| Filter | Threshold | Rationale |
|---|---|---|
| PASS only (default) | Sentieon TNfilter | Removes soft-filtered variants (strand bias, weak evidence, etc.) |
| FORMAT/DP | >= 99 | Ensures reliable VAF estimation; at DP < 99, AF noise is too high |
| AF | >= 0.03 | Below this, allele fractions are unreliable for ratio calculation |
| gnomAD AF | >= 0.40 | Removes hom/het artefact that mimics 2:1 contamination ratio |

The pipeline uses 4 piped bcftools commands:
```bash
bcftools view -f PASS input.vcf.gz -Ou \
| bcftools filter -e '(FORMAT/DP<99 || AF<0.03)' -Ou \
| bcftools +split-vep --columns - -a CSQ -p CSQ_ -s worst -Ou \
| bcftools view -e 'CSQ_gnomADg_AF>=0.40 || CSQ_gnomADe_AF>=0.40' -Oz -o out.vcf.gz
```

The `+split-vep` step extracts CSQ subfields (including gnomAD AF and gene
SYMBOL) into properly-typed INFO tags, enabling the arithmetic comparison in
the final step. It uses `-s worst` (single worst-consequence transcript per
variant) to avoid duplicate records.

### Flagging threshold: `--peak-count` (default 6)

A pair is flagged when:

- `n_shared >= --min-shared` (default 10) — enough shared variants to assess
- **AND** `peak_count >= --peak-count` (default 6) — the non-unity peak has a
  statistically significant cluster

**Derivation from null model:**

Under the null hypothesis (no contamination), variants that fall outside the
unity zone are scattered across the non-unity region with no systematic
clustering. With a bin width of 0.2 log2 units and a total range of ±4 log2,
there are 37 non-unity bins. If k variants are distributed uniformly across
these 37 bins, the probability that the tallest bin reaches a given count
follows the maximum of a multinomial distribution.

Monte Carlo simulation (100,000 iterations, ~20 non-unity variants per pair):

| Threshold | P(exceed) per pair | Expected FP in 1128 pairs (48 samples) |
|---|---|---|
| >= 4 | 0.035 | 39.5 |
| >= 5 | 0.006 | 6.8 |
| >= 6 | 0.0004 | 0.4 |
| >= 7 | 0.00001 | 0.01 |
| >= 8 | ~0 | ~0 |

The default of **6** gives < 1 expected false positive for cohorts up to ~100
samples. For larger cohorts (96+ samples, 4560+ pairs), increase to 7.

**Linkage disequilibrium consideration:** Multiple variants in the same gene on
the same haplotype are not independent — they share the same genotype and
produce the same ratio. In testing, linked germline SNPs (e.g. 5 PRPF8 variants
on one haplotype) created false clusters of 4–5 variants at the same ratio.
The threshold of 6 absorbs this because linked clusters from a single gene
typically contribute ≤ 5 correlated variants, while genuine contamination
produces markers from many independent genomic loci.

---

## Requirements

- Python >= 3.9
- `bcftools` >= 1.14 with `split-vep` plugin
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
| `--max-gnomad` | `0.40` | Exclude variants with gnomAD AF >= this (removes hom/het artefact). Set to 1.0 to disable |
| `--min-af` | `0.03` | VAF floor applied during pre-filtering |
| `--min-dp` | `99` | Minimum read depth |
| `--min-shared` | `10` | Minimum shared variants to assess a pair |
| `--peak-count` | `6` | Flag if non-unity peak has >= N variants (statistically derived; see above) |
| `--threads / -t` | `min(8, nCPU)` | Parallel threads |
| `--include-non-pass` | off | Include non-PASS variants |
| `--plots` | off | Generate histogram plots for flagged pairs |
| `--max-output` | `10` | Maximum flagged pairs to write detail TSVs/plots for (ranked by peak_count). Set to 0 for no limit |
| `--force-refilter` | off | Regenerate filtered VCFs |
| `--bin-width` | `0.2` | Histogram bin width (log2 units) |
| `--vcf-glob` | `*_annotated.vcf.gz` | File matching pattern |
| `--verbose / -v` | off | Debug logging |

### Examples

```bash
# Standard run
python contamination_screen.py /data/vcfs/ --outdir results/ --plots

# More sensitive (smaller panels or lower contamination)
python contamination_screen.py /data/vcfs/ \
    --min-shared 5 --peak-count 5 --outdir results_sensitive/

# Larger cohorts (96+ samples)
python contamination_screen.py /data/vcfs/ --peak-count 7
```

---

## Output files

```
results/
├── filtered/                          Quality-filtered VCFs (reused on re-runs)
├── summary.tsv                        One row per pair (dual-peak results)
├── matrix.tsv                         N x N directional contamination matrix
├── flagged_pairs/*.tsv                Variant-level detail (top N by peak_count)
└── plots/*.png                        Log2-ratio histograms (top N, --plots only)
```

`summary.tsv` and `matrix.tsv` always report **all** pairs. Detail TSVs and
plots are only generated for the top `--max-output` (default 10) flagged pairs,
ranked by `peak_count` descending (most suspicious first). This keeps output
manageable for large cohorts. Set `--max-output 0` to output all flagged pairs.

If `--plate-layout` is provided, the matrix rows and columns are ordered by
plate position (column-major: A1, B1...H1, A2, B2...up to the last column),
making adjacent-well
contamination appear as non-zero entries near the diagonal. Without
`--plate-layout`, samples are ordered by input filename.

### `summary.tsv` columns

| Column | Description |
|---|---|
| `sample_a`, `sample_b` | Sample pair |
| `n_shared` | Total shared variants (after all quality + gnomAD filters) |
| `overall_log2`, `overall_count`, `overall_fraction` | Tallest peak across full range (usually ratio=1 from germline sharing) |
| `peak_log2_ratio`, `peak_ratio` | Non-unity peak centre (the contamination signal) |
| `peak_count`, `peak_fraction` | Strength of non-unity peak |
| `contamination_source` | Source sample (higher VAF) |
| `contamination_recipient` | Recipient sample (diluted VAF) |
| `contamination_fraction` | Estimated fraction of recipient library from source |
| `flagged` | TRUE if meets thresholds |

---

## Interpretation

- **Non-unity peak with peak_count >= 6:** cross-sample contamination at the
  implied fraction. Direction indicates which sample is the source. The peak_count
  reflects how many independent genomic loci confirm the signal.
- **Large overall peak at ratio=1:** many variants shared at similar VAF. In
  unrelated samples this suggests a sample swap or duplicate; in related samples
  (same patient) it reflects shared biology.
- **Non-unity peak_count of 4-5:** borderline — may represent linked germline
  variants (e.g. multiple SNPs in one gene on the same haplotype) rather than
  true contamination. Inspect the detail TSV to check whether variants span
  multiple chromosomes (contamination) or cluster in one gene (linkage).
- **No significant peaks:** clean pair with minimal variant sharing.
- The tool is designed for **unrelated samples**. Same-patient pairs will have
  large overall peaks at ratio=1 (shared clonal mutations) which is expected and
  not flagged as contamination.
- **Transitive contamination:** if sample Y is heavily contaminated by source Z,
  then Y will also flag against many other samples that independently carry Z's
  variants at germline het frequency. The diagnostic pattern is: one sample
  appears as recipient in many flagged pairs, all at the same fraction, and is
  never a source. The true source is the pair with the highest peak_count.
  Consider excluding confirmed-contaminated samples and re-running.

---

## Testing and validation

The tool was developed and validated using 3 VCFs from different patients
processed on the same Uranus sequencing run (26TULIP24), where cross-sample
contamination between adjacent samples was suspected.

| Pair | True status | n_shared | Non-unity peak_count | Flagged |
|---|---|---|---|---|
| A vs B | **Contaminated** (~50%) | 57 | 17 | **Yes** |
| B vs C | Clean | 31 | 5 (linked PRPF8) | No |
| C vs A | Clean | 25 | 4 (linked PRPF8) | No |

The contaminated pair shows 17 variants from multiple independent loci (PIK3CD,
CUX1, DNMT3A, ATM, IDH2, TP53, ASXL1, U2AF1, etc.) all at ratio 2:1,
correctly identifying ~50% contamination with the direction A→B. The clean pairs
show only 4-5 linked variants from single genes (PRPF8, SRCAP) which fall below
the flagging threshold.
