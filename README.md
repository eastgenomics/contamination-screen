# contamination_screen

Pairwise cross-sample contamination detector for tumour sequencing panel VCFs.

## Documentation

Supplementary HTML documents are available both in this repository and hosted online.

### Algorithm and code review

| Document | In repo | Online |
|---|---|---|
| **Algorithm briefing** — visual guide to the contamination signal, filtering, dual-peak detection, flagging thresholds, false positives, and validation | [Contamination_Screen_Code_Review_Briefing.html](Contamination_Screen_Code_Review_Briefing.html) | [view online](https://003-260510-docs-for-sharing.s3.eu-west-2.amazonaws.com/Contamination_Screen_Code_Review_Briefing.html) |

### Code explainer series

| Part | Topic | In repo | Online |
|---|---|---|---|
| 1 | High-level outline — 4 phases, data flow, function inventory | [Explainer 1](Contamination_Screen_Code_Explainer_1.html) | [view online](https://003-260510-docs-for-sharing.s3.eu-west-2.amazonaws.com/Contamination_Screen_Code_Explainer_1.html) |
| 2 | Module structure — imports, constants, naming conventions | [Explainer 2](Contamination_Screen_Code_Explainer_2.html) | [view online](https://003-260510-docs-for-sharing.s3.eu-west-2.amazonaws.com/Contamination_Screen_Code_Explainer_2.html) |
| 3 | Key concepts — log2 ratios, histograms, z-scores, subprocess pipes, threading | [Explainer 3](Contamination_Screen_Code_Explainer_3.html) | [view online](https://003-260510-docs-for-sharing.s3.eu-west-2.amazonaws.com/Contamination_Screen_Code_Explainer_3.html) |
| 4 | Function walkthrough — all 17 functions with purpose and logic | [Explainer 4](Contamination_Screen_Code_Explainer_4.html) | [view online](https://003-260510-docs-for-sharing.s3.eu-west-2.amazonaws.com/Contamination_Screen_Code_Explainer_4.html) |
| 5 | Line-by-line reference — dense annotations for the 4 hardest functions | [Explainer 5](Contamination_Screen_Code_Explainer_5.html) | [view online](https://003-260510-docs-for-sharing.s3.eu-west-2.amazonaws.com/Contamination_Screen_Code_Explainer_5.html) |

---

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

### How it works

For every pair of samples on a run, the tool merges their variant sets, computes
`log2(VAF_B / VAF_A)` for each shared variant, and builds a histogram. If one
sample's DNA is present in the other's library, the source's germline variants
appear at a consistent diluted fraction — creating a peak away from ratio=1.
The tool identifies this non-unity peak, estimates the contamination fraction,
and reports which sample is the source. See the **algorithm briefing** above for
a full visual explanation of the method, thresholds, and known limitations.

---

## Expected input VCFs

Input VCFs must have been produced by the Uranus pipeline eggd_vep stage and
must have already undergone:

| Step | Command | Header tag |
|---|---|---|
| Variant calling | Sentieon TNhaplotyper2 (tumour-only) | `##SentieonCommandLine.TNhaplotyper2` |
| Soft-filtering | Sentieon TNfilter | `##SentieonCommandLine.TNfilter` |
| Normalisation | `bcftools norm -m -any --keep-sum AD` | `##bcftools_normCommand` |
| VEP annotation | Ensembl VEP (CSQ field) | `##INFO=<ID=CSQ,...>` |

**The normalisation step (`bcftools norm -m -any`) must already have been
applied.** This tool does not re-normalise.

---

## Known false positive patterns

Four signal patterns can produce false or misleading results. See the algorithm
briefing for full diagnostic criteria.

- **Somatic LOH** — high-purity tumours with cancer driver genes at AF ≈ 1.0
  produce a false peak exactly at log2 = −1.0 against any het sample. Caught by
  the FREEMIX gate (FREEMIX will be normal, <1%).
- **Bone marrow chimerism** — post-HSCT mixed DNA is mechanically identical to
  contamination; requires clinical context to distinguish. FREEMIX will be
  elevated (18–47%).
- **Engineered cell line controls** — produce the highest signals in any cohort.
  Exclude with `--exclude "*Q*"` via `dx_fetch.py`.
- **Transitive contamination** — a heavily contaminated sample may spuriously
  appear as a source against other samples. Identified by inspecting pair
  topology in `summary.tsv`.

---

## Correlation with VerifyBamID FREEMIX

The two tools are complementary: VerifyBamID identifies recipients; this tool
identifies sources and estimates the fraction. Sources have near-normal FREEMIX
(~0.3%) because their own library is clean.

| FREEMIX band | % also flagged by contamination screen |
|---|---|
| < 5% | ~14% |
| 5–25% | ~8% |
| > 25% | **89%** |

---

## Validation

| Metric | Value |
|---|---|
| Runs screened (retrospective cohort) | 100 |
| Total pairwise comparisons | 97,310 |
| Flagged pairs at default thresholds | 128 |
| Runs with ≥ 1 flagged pair | 43 / 100 |
| Strongest signal | peak_count=28, n_shared_z=13.4 |
| Pearson r(FREEMIX, peak_count) within 5–50% FREEMIX band | 0.64 |

---

## Requirements

- Python >= 3.9
- `bcftools` >= 1.14 with `split-vep` plugin
- `pandas` and `numpy`
- `matplotlib` (optional, for `--plots`)
- `dxpy` (optional, for `dx_fetch.py` — fetching from DNAnexus)

Input VCFs must be bgzipped (`.vcf.gz`) and tabix-indexed (`.vcf.gz.tbi`).

## Installation

```bash
git clone <repo>
cd contamination-screen
python3 -m venv .venv
.venv/bin/pip install pandas numpy matplotlib dxpy
```

---

## Fetching data from DNAnexus

### dx_grab (dependency)

[dx_grab](https://github.com/eastgenomics/dx-grab) is a shared DNAnexus utility
library that `dx_fetch.py` depends on for all platform interactions — authentication,
project resolution, file discovery, archive handling, and download. It must be
present before running `dx_fetch.py`.

Clone it alongside this repo:

```bash
git clone https://github.com/eastgenomics/dx-grab.git ~/Documents/dx_grab
```

By default `dx_fetch.py` looks for dx_grab at `~/Documents/dx_grab`. If you
clone it elsewhere, set the `DXGRAB_DIR` environment variable:

```bash
export DXGRAB_DIR=/path/to/dx_grab
```

### dx_fetch.py

`dx_fetch.py` is a precursor helper that downloads everything needed for a
single MYE run from DNAnexus and prints the ready-to-run
`contamination_screen.py` command.

```bash
python dx_fetch.py PROJECT [--output DIR] [--exclude PATTERN]
                           [--yes] [--skip-archived] [--skip-existing]
                           [--dry-run]
```

`PROJECT` may be a project ID (`project-xxx...`), an exact project name, or a
glob pattern (e.g. `002_260423*MYE`) that must resolve to exactly one project.

Downloads into `<output>/` (default: `./<project_name>/`):

| Path | Contents |
|---|---|
| `vcfs/*.vcf.gz` + `*.tbi` | eggd_vep annotated VCFs, one per sample |
| `multiqc_general_stats.txt` | Plate well positions (for `--plate-layout`) |
| `multiqc_verifybamid.txt` | VerifyBamID FREEMIX values (for `--freemix-file`) |

Control samples (specimen IDs containing `Q`, e.g. `25357Q0020`) are excluded
by default. Pass `--no-exclude-controls` to override.

```bash
# List what would be downloaded (no download)
python dx_fetch.py '002_260423*MYE' --dry-run

# Download — output defaults to ./002_260423_A01303_0760_AHLTLCDRX7_MYE/
python dx_fetch.py '002_260423*MYE'

# Non-interactive: skip archived VCFs automatically
python dx_fetch.py '002_260423*MYE' --skip-archived

# Non-interactive: submit unarchive requests and wait
python dx_fetch.py '002_260423*MYE' --yes
```

On success, `dx_fetch.py` prints the `contamination_screen.py` command to run,
with `--plate-layout` and `--freemix-file` automatically filled in if those
files were found.

---

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
| `--n-shared-z` | `2.0` | Flag only if within-run n_shared z-score is at least this value (combined with `--peak-count`) |
| `--threads / -t` | `min(8, nCPU)` | Parallel threads |
| `--include-non-pass` | off | Include non-PASS variants |
| `--plots` | off | Generate histogram plots for flagged pairs |
| `--max-output` | `10` | Maximum flagged pairs to write detail TSVs/plots for (ranked by peak_count). Set to 0 for no limit |
| `--force-refilter` | off | Regenerate filtered VCFs even if cached |
| `--bin-width` | `0.2` | Histogram bin width (log2 units) |
| `--vcf-glob` | `*_annotated.vcf.gz` | File matching pattern |
| `--plate-layout` | _none_ | MultiQC general-stats file; orders matrix by plate well position (column-major) |
| `--freemix-file` | _none_ | MultiQC `multiqc_verifybamid.txt`; adds `recipient_freemix` column and gates flagging |
| `--freemix-threshold` | `0.15` | Minimum recipient FREEMIX fraction (0–1) required to flag when `--freemix-file` is provided |
| `--verbose / -v` | off | Debug logging |

### Examples

```bash
# Standard run
python contamination_screen.py /data/vcfs/ --outdir results/ --plots

# With FREEMIX gating (recommended)
python contamination_screen.py /data/vcfs/ \
    --plate-layout multiqc_general_stats.txt \
    --freemix-file multiqc_verifybamid.txt \
    --outdir results/ --plots

# More sensitive (smaller panels or lower contamination)
python contamination_screen.py /data/vcfs/ \
    --min-shared 5 --peak-count 8 --n-shared-z 1.5 --outdir results_sensitive/
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
ranked by `peak_count` descending. Set `--max-output 0` to output all flagged pairs.

If `--plate-layout` is provided, matrix rows and columns are ordered by plate
position (column-major: A1, B1...H1, A2...), making adjacent-well contamination
appear as non-zero entries near the diagonal.

### `summary.tsv` columns

| Column | Description |
|---|---|
| `sample_a`, `sample_b` | Sample pair |
| `n_shared` | Total shared variants after all filters |
| `n_shared_z` | Within-run z-score of `n_shared` across all pairs in this run |
| `overall_log2`, `overall_count`, `overall_fraction` | Tallest peak across full range (usually ratio=1 from germline sharing) |
| `peak_log2_ratio`, `peak_ratio` | Non-unity peak centre (the contamination signal) |
| `peak_count`, `peak_fraction` | Strength of non-unity peak |
| `contamination_source` | Source sample (higher VAF) |
| `contamination_recipient` | Recipient sample (diluted VAF) |
| `contamination_fraction` | Estimated fraction of recipient library from source |
| `recipient_freemix` | VerifyBamID FREEMIX for recipient (when `--freemix-file` supplied) |
| `flagged` | TRUE if all thresholds met |

---

## Interpretation

- **Flagged pair** (`peak_count ≥ 10`, `n_shared_z ≥ 2.0`): strong evidence of
  cross-sample contamination. `contamination_source` identifies the origin;
  `contamination_fraction` estimates the level.
- **Large overall peak at ratio=1**: many variants at identical VAF — suggests a
  sample swap or duplicate rather than contamination.
- **Borderline** (`peak_count 7–9`, `n_shared_z 1.5–2.0`): inspect the detail
  TSV. Variants spanning multiple chromosomes indicate contamination; variants
  clustering in one gene indicate linkage artefact.
- **High peak_count, normal FREEMIX, peak exactly −1.0 against many samples**:
  likely somatic LOH artefact — inspect in-peak source AFs (expect ~1.0, not 0.5).
- **No significant peaks**: clean pair.
