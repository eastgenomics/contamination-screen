# contamination_screen

Pairwise cross-sample contamination detector for tumour sequencing panel VCFs.

## Overview

This tool screens a cohort of annotated tumour VCF files for cross-sample
contamination by detecting shared variants at consistent VAF ratios. It is
designed for targeted haematological oncology panel sequencing data processed by
the **Uranus clinical pipeline** (CUH Bioinformatics), using VCFs output by the
**eggd_vep** stage.

Validated on a retrospective cohort of 100 consecutive haematological oncology
sequencing runs (97,310 pairwise comparisons, ~4,400 samples). With default
thresholds (`--peak-count 10`, `--n-shared-z 2.0`), 128 pairs were flagged
across 43 runs.

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
  peak_count = 4–5 which is below the flagging threshold of 10.

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

#### VEP configuration versions and gnomAD field availability

The gnomAD filter step uses whichever gnomAD AF fields are present in the VCF's
CSQ format string. Two configurations have been observed in Uranus pipeline
output:

| VEP config | gnomAD fields in CSQ | Runs (100-run cohort) | Approximate date boundary |
|---|---|---|---|
| Newer (eggd_vep ≥ v1.3) | `gnomADg_AF` + `gnomADe_AF` | RUN001–070 | Aug 2025 onwards |
| Older (eggd_vep < v1.3) | `gnomADe_AF` only | RUN071–100 | Before Aug 2025 |

The tool detects which fields are present at runtime using
`bcftools +split-vep -l` and builds the filter expression accordingly. If
neither field is present a warning is logged and the gnomAD step is skipped.

**Impact on `n_shared`:** runs annotated with the older VEP config produce
systematically fewer shared variants per pair (median `n_shared` ≈ 8) compared
to newer-config runs (median ≈ 18–20). This difference arises because
`gnomADe_AF` (exome-derived) assigns higher AF estimates to many coding variants
than `gnomADg_AF` (genome-derived), causing the older config to remove more
variants at the AF ≥ 0.40 threshold. The global SD for each group is internally
consistent (SD ≈ 4.9 for newer runs, SD ≈ 2.9 for older runs), but the two
populations cannot be pooled for z-scoring.

As a consequence, **`n_shared` z-scores are always computed within a single run**
(the default behaviour) rather than against a cross-run global distribution. A
global mean/SD blends the two VEP-config populations and produces z-scores that
are systematically biased in both directions.

### Dual flagging threshold: `peak_count` and `n_shared_z`

A pair is flagged when **all three** conditions are met:

- `n_shared >= --min-shared` (default 10) — enough shared variants to assess
- `peak_count >= --peak-count` (default 10) — the non-unity peak contains a
  statistically significant cluster of variants at a consistent non-unity ratio
- `n_shared_z >= --n-shared-z` (default 2.0) — the pair shares more variants
  than typical within this run (within-run z-score of `n_shared`)

The `n_shared_z` criterion adds a second, independent line of evidence:
contamination increases the number of detectable shared variants between a pair
(the source's variants appear in the recipient), so contaminated pairs have
elevated `n_shared` relative to clean pairs from the same run. Using both
criteria simultaneously substantially reduces false positives from linkage
disequilibrium and somatic LOH artefacts (see Known false positives below).

**`n_shared_z` is computed within each run** across all pairwise comparisons in
that cohort. It is stored in `summary.tsv` and can be inspected independently of
the flagging decision.

**Relationship between `peak_count` and `n_shared_z`:**

In the 100-run retrospective cohort (Pearson r within the 5–50% FREEMIX band):
- r(FREEMIX, peak_count) = 0.639
- r(FREEMIX, n_shared_z) = 0.569

Both metrics track the same underlying event, but `peak_count` is more
sensitive to the directionality of the signal (ratio clustering) while
`n_shared_z` reflects the magnitude of variant sharing.

**Linkage disequilibrium consideration:** Multiple variants in the same gene on
the same haplotype are not independent — they share the same genotype and
produce the same ratio. In testing, linked germline SNPs (e.g. 5 PRPF8 variants
on one haplotype) created false peak_count clusters of 4–5. These typically have
normal `n_shared_z` because the total number of shared variants is not elevated.
The dual threshold removes them: genuine contamination produces elevated signals
on both metrics simultaneously.

---

## Known false positive patterns

### 1. Somatic LOH / homozygous tumour variants at cancer gene loci

Tumour-only VCFs from high-purity haematological cancers can generate a
systematic false positive. If a patient has acquired homozygous somatic variants
(or loss of heterozygosity) at positions in cancer driver genes (BRAF, FLT3,
TET2, PRPF8, JAK1, etc.), these appear at AF ≈ 1.0 in their VCF. An unrelated
patient who is heterozygous germline at the same positions will have AF ≈ 0.5.
The ratio log2(0.5/1.0) = −1.0 appears for every such position, producing a
peak at exactly −1.0 that exactly mimics 50% contamination.

**Diagnostic features of this artefact:**
- `peak_log2_ratio` is exactly −1.0 (or +1.0) for virtually every pair involving
  the affected sample
- The same sample appears as "source" against many unrelated samples on the run,
  all at contamination_fraction = 0.5
- VerifyBamID FREEMIX is normal (~0.1–0.7%) for both samples — the BAM contains
  no foreign DNA
- In-peak variants cluster in haematological cancer driver genes
- In-peak variant AF in the "source" sample is ≈ 1.0 (homozygous), not ≈ 0.5
  (germline het)

The `n_shared_z` threshold provides partial protection: these artefact pairs
tend to have elevated `n_shared` because the homozygous tumour variants appear
in many paired samples, but the effect is run-wide rather than targeted, making
the z-score moderate rather than extreme. Definitive diagnosis requires inspecting
the variant-level detail TSV.

The gnomAD ≥ 0.40 filter removes the common-germline version of this artefact
(where both genotypes are population-common), but somatic variants at cancer
gene positions with population AF < 0.40 are not filtered and remain a source
of false positives.

### 2. Transitive contamination

If sample Y is heavily contaminated by source Z (~50% contamination), Y now
carries Z's variants at ~25% VAF. When Y is compared to other samples that
independently carry Z's variants at germline het frequency (~50%), the
contamination screen may flag Y as a "source" contaminating those samples, even
though Y is a recipient. This is the transitive artefact.

**Diagnostic pattern:** one sample appears as recipient in many high-confidence
flagged pairs at a consistent fraction, and also appears spuriously as source
in a small number of lower-confidence pairs. The true primary event is the pair
with the highest `peak_count` and `n_shared_z` involving the recipient sample.

### 3. Bone marrow chimerism (post-transplant samples)

Patients who have undergone allogeneic haematopoietic stem cell transplantation
(HSCT) have mixed donor/recipient DNA in their blood. This produces a signal
that is **mechanistically identical** to cross-sample contamination: the donor's
germline heterozygous variants appear at ~50% VAF alongside the patient's own
variants, giving a consistent non-unity log2-ratio peak at every pair comparison
between the chimeric sample and any other sample on the run.

The tool cannot distinguish bone marrow chimerism from true cross-sample
contamination using VCF data alone. Clinical context (transplant history) is
required.

**Diagnostic features:**
- The chimeric sample appears as recipient against many (often most) other
  samples on the run, all at contamination_fraction ≈ 0.5
- The signal is **reproducible across multiple sequencing runs** of the same
  patient — chimerism is a biological property of the sample, not a run-level
  artefact
- VerifyBamID FREEMIX will be elevated (typically 40–50%), reflecting the
  mixed DNA — this is the correct result, not a false positive by VerifyBamID
- The 'source' samples identified by the tool will be unrelated patients who
  happen to share common germline variants with the donor genome

**Known examples (retrospective 100-run cohort):**
- `25351K0007` — flags in RUN034 (26NGSHO2) and RUN037 (25NGSHO77); same
  patient sequenced twice, same signal reproduced at ~50% fraction both times
- `25238K0026` — same runs, same mechanism

These samples are listed in `retrospective/known_contaminated_samples.tsv`
alongside confirmed true contamination events.

---

## Correlation with VerifyBamID FREEMIX

VerifyBamID FREEMIX estimates cross-individual contamination from the BAM at
germline SNP positions. In the 100-run cohort (1,432 samples with both metrics):

| FREEMIX band | % also flagged by contamination screen |
|---|---|
| < 5% | ~14% |
| 5–25% | ~8% |
| > 25% | **89%** |

The two tools are **complementary, not redundant:**

- VerifyBamID detects recipients sensitively but cannot identify the source,
  quantify the fraction, or distinguish primary from transitive events
- The contamination screen identifies direction (source → recipient), estimates
  the fraction, and surfaces the primary event — but requires a clustering
  signal that may be absent at low contamination levels (< ~5% FREEMIX)
- **Sources** have near-normal FREEMIX (~0.3% median) because their own library
  is clean; the contamination signal is in their victim

Samples with FREEMIX > 5% but not flagged by the contamination screen are either
genuine low-level contamination events below the detection threshold, or reflect
background VerifyBamID noise in this panel type.

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
cd contamination-screen
python3 -m venv .venv
.venv/bin/pip install pandas numpy matplotlib
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
| `--peak-count` | `10` | Flag if non-unity peak has >= N variants |
| `--n-shared-z` | `2.0` | Flag only if within-run n_shared z-score exceeds this (combined with `--peak-count`) |
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
    --min-shared 5 --peak-count 8 --n-shared-z 1.5 --outdir results_sensitive/

# Larger cohorts (96+ samples)
python contamination_screen.py /data/vcfs/ --peak-count 10 --n-shared-z 2.0
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
making adjacent-well contamination appear as non-zero entries near the diagonal.
Without `--plate-layout`, samples are ordered by input filename.

### `summary.tsv` columns

| Column | Description |
|---|---|
| `sample_a`, `sample_b` | Sample pair |
| `n_shared` | Total shared variants (after all quality + gnomAD filters) |
| `n_shared_z` | Within-run z-score of `n_shared` across all pairs in this cohort |
| `overall_log2`, `overall_count`, `overall_fraction` | Tallest peak across full range (usually ratio=1 from germline sharing) |
| `peak_log2_ratio`, `peak_ratio` | Non-unity peak centre (the contamination signal) |
| `peak_count`, `peak_fraction` | Strength of non-unity peak |
| `contamination_source` | Source sample (higher VAF) |
| `contamination_recipient` | Recipient sample (diluted VAF) |
| `contamination_fraction` | Estimated fraction of recipient library from source |
| `flagged` | TRUE if meets all three thresholds (n_shared, peak_count, n_shared_z) |

---

## Interpretation

- **Flagged pair (peak_count ≥ 10, n_shared_z > 2.0):** strong evidence of
  cross-sample contamination at the implied fraction. Direction indicates which
  sample is the source. The peak_count reflects how many independent genomic
  loci confirm the signal; n_shared_z confirms the pair shares more variants
  than typical for this run.
- **Large overall peak at ratio=1:** many variants shared at similar VAF. In
  unrelated samples this suggests a sample swap or duplicate; in related samples
  (same patient) it reflects shared biology.
- **peak_count 7–9 with n_shared_z 1.5–2.0:** borderline — inspect the detail
  TSV. Check whether variants span multiple chromosomes (contamination) or
  cluster in one or two genes (linkage or LOH artefact).
- **High peak_count but normal FREEMIX, peak_log2_ratio exactly −1.0 against
  many samples:** likely somatic LOH artefact (see Known false positives).
  Inspect in-peak variant AFs — if source AF ≈ 1.0, this is not contamination.
- **No significant peaks:** clean pair with minimal variant sharing.
- The tool is designed for **unrelated samples**. Same-patient pairs will have
  large overall peaks at ratio=1 (shared clonal mutations) which is expected and
  not flagged as contamination.
- **Transitive contamination:** if sample Y is heavily contaminated by source Z,
  Y will also flag against many other samples. The diagnostic pattern is: one
  sample appears as recipient in many flagged pairs, all at the same fraction,
  and is never a source. The true source is the pair with the highest peak_count.
  Consider excluding confirmed-contaminated samples and re-running.

---

## Testing and validation

### Initial development (26TULIP24)

The tool was developed and validated using 3 VCFs from different patients
processed on the same Uranus sequencing run (26TULIP24), where cross-sample
contamination between adjacent plate wells was suspected. Confirmed by plate
adjacency (wells A4/B4), VerifyBamID FREEMIX (45.5% / 13.9%), and coverage.

| Pair | True status | n_shared | peak_count | n_shared_z | Flagged |
|---|---|---|---|---|---|
| A vs B | **Contaminated** (~50%) | 57 | 17 | 7.20 | **Yes** |
| B vs C | Clean | 31 | 5 (linked PRPF8) | 1.50 | No |
| C vs A | Clean | 25 | 4 (linked PRPF8) | 0.59 | No |

The contaminated pair shows 17 variants from multiple independent loci (PIK3CD,
CUX1, DNMT3A, ATM, IDH2, TP53, ASXL1, U2AF1, etc.) all at ratio 2:1,
correctly identifying ~50% contamination with the direction A→B. The clean pairs
show only 4–5 linked variants from single genes (PRPF8, SRCAP) which fall below
threshold on both peak_count and n_shared_z.

### Retrospective cohort (100 runs, May 2026)

The tool was applied retrospectively to 100 consecutive haematological oncology
sequencing runs (Mar 2025 – May 2026) using the `retrospective/` analysis
scripts. Key findings:

| Metric | Value |
|---|---|
| Runs screened | 100 |
| Total pairwise comparisons | 97,310 |
| Flagged pairs (peak_count ≥ 10, n_shared_z > 2.0) | 128 |
| Runs with ≥ 1 flagged pair | 43 / 100 |
| Strongest signal | peak_count=28, n_shared_z=13.4 (RUN049) |

**Correlation with VerifyBamID FREEMIX:** Pearson r = 0.19 across all samples;
r = 0.64 within the 5–50% FREEMIX band, confirming the tools detect the same
events. All primary contamination recipients (FREEMIX > 25%) were also flagged
by the screen; sources were correctly identified with near-normal FREEMIX.

**VEP config boundary:** Runs before approximately August 2025 (eggd_vep < v1.3,
`gnomADe_AF` only) have median n_shared ≈ 8 vs ≈ 18–20 for newer runs. This
arises from gnomADe assigning higher AF estimates to coding variants, causing
the AF ≥ 0.40 filter to remove more variants. Within-run z-scoring correctly
normalises for this difference; a global SD is not appropriate across the full
cohort.
