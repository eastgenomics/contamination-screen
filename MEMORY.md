# Session Memory — contamination_screen development
# 2026-05-07

## Project
`~/Documents/contamination/` — git repo
Tool: `contamination_screen.py` — pairwise cross-sample contamination detector for the Uranus haematological oncology panel sequencing pipeline (CUH Bioinformatics).

---

## Context that drove the work

Input VCFs are produced by the **`eggd_vep` stage** of the Uranus pipeline.

Three VCFs in `~/Downloads/` for patient TULIP24-5877:
- `26098K0076` suffix `-M-` (bone marrow, later timepoint) — 720 variants
- `26093K0005` suffix `-M-` (bone marrow, earlier timepoint) — 452 variants
- `26097K0043` suffix `-F-` (different sample type) — 1017 variants

Manual pairwise analysis (bcftools isec + VAF ratio) showed:
- A vs B: 307 shared variants; 40 PASS variants at ratio ~0.5 (AF_A ~50%, AF_B ~25%) — clonal expansion between timepoints, NOT contamination. Same-patient serial marrow samples. Genes: DNMT3A, ASXL1, CUX1, TP53, KRAS, ATM, IKZF1, NF1, PRPF8, SETBP1, etc.
- A vs C (F sample): different pattern — near-100% variants in A appearing at ~50% in C. Likely LOH/CN gain in marrow vs heterozygous in F sample. Most clonal expansion variants from A vs B are absent in C entirely.
- The tool was requested because the shared variants are thought to represent contamination between **unrelated** samples in a 48-sample cohort.

---

## Design decisions

### Contamination signal
Shared rare somatic variants at a consistent VAF ratio = contamination fingerprint.
`log2(VAF_B / VAF_A)` histogram; peak away from 0 = contamination.
- peak_log2 < 0: A is source (A's variants diluted in B)
- peak_log2 > 0: B is source
- No central exclusion zone: peak-finder scans full log2 range including ratio ~1 (sample swaps)
- contamination_fraction = 2^|peak_log2|

### Architecture
1. Pre-filter each VCF once → `results/filtered/` (reused on re-runs)
2. Load all filtered VCFs into memory as pandas DataFrames
3. N*(N-1)/2 pairwise comparisons (multiprocessing for >50 pairs, sequential for small cohorts)
4. Log2-ratio histogram peak finding per pair
5. Output: summary.tsv, matrix.tsv, flagged_pairs/*.tsv, plots/*.png (optional)

Parallel workers return stats dict only (no DataFrames). DataFrames for flagged pairs are recomputed in main process to avoid pickle overhead.

---

## The bcftools filter pipeline — key technical finding

The Uranus clinical pipeline uses:
```python
cmd = (
    f"bcftools +split-vep --columns - -a CSQ -Ou -p 'CSQ_' -d {vcf} | "
    f"bcftools annotate -x INFO/CSQ -o {output_vcf}"
)
```
followed by:
```
bcftools filter --soft-filter "EXCLUDE" -m + -e '(CSQ_Prev_Count_AC>853 || CSQ_gnomADg_AF >= 0.002 || CSQ_gnomADe_AF >= 0.002) || (CSQ_SYMBOL!="GATA2" || CSQ_Consequence!="synonymous_variant") && (CSQ_SYMBOL!="TP53" || CSQ_Consequence!="synonymous_variant") && (CSQ_Consequence=="synonymous_variant")'
bcftools filter -e '(FORMAT/DP<99 || AF<0.03)'
```

### Why `bcftools filter -e 'CSQ_...'` fails on the raw annotated VCFs
The gnomAD and Prev_Count fields are annotated as standalone INFO tags typed as **String** in the VCF header (e.g. `gnomADg_AF,Type=String`). bcftools cannot do arithmetic comparisons on String-typed fields.

### Why `bcftools +split-vep -c FIELDNAME:Type` fails
The field names (gnomADg_AF etc.) exist as String-typed INFO tags already. split-vep says "Existing header definitions will not be overwritten" and fails with "No such column" when trying to re-declare them. **Index-range notation works (`-c 1-2,24,27,31`) but name-based lookup does not for these fields.**

### Working solution
Using `-p 'CSQ_'` prefix produces `CSQ_gnomADg_AF` etc. — **different names** from the existing String INFO tags. No conflict. The .*_AF built-in type rule in split-vep automatically assigns Float to the gnomAD AF fields.

`CSQ_Prev_Count_AC` is still assigned String (no built-in pattern matches). Fix: `bcftools annotate -x INFO/CSQ_Prev_Count_AC -h prevcount.hdr` where the hdr file declares it as Integer.

### Full 6-step pipeline (in filter_vcf())
```
1. bcftools +split-vep --columns - -a CSQ -Ou -p 'CSQ_' -d    [exact Uranus]
2. bcftools annotate -x INFO/CSQ                               [exact Uranus]
3. bcftools annotate -x INFO/CSQ_Prev_Count_AC -h fix_hdr      [extra - type fix]
4. bcftools filter --soft-filter EXCLUDE -m + -e '...'         [exact Uranus]
5. bcftools filter -e '(FORMAT/DP<99 || AF<0.03)'              [exact Uranus]
6. bcftools view [-f PASS] -e 'FILTER~"EXCLUDE"'               [hard-filter]
```

`bcftools norm -m -any` is already in the VCF headers (`##bcftools_normCommand`) — confirmed present, not re-run.

### split-vep -d flag
Creates one VCF record per transcript (duplicate mode). Most duplicates are removed by downstream filters (different transcripts = different CSQ_Consequence/CSQ_SYMBOL). Survivors deduplicated in `load_vcf()` on (CHROM, POS, REF, ALT), keeping first record (VEP's worst-consequence transcript).

### Variant counts after full filter (PASS-only)
- `26093K0005`: 15 variants
- `26097K0043`: 8 variants
- `26098K0076`: 10 variants
(vs 452/1017/720 unfiltered)

---

## VCF structure notes (Uranus annotated VCFs)

- Sample name in `#CHROM` line column 10 (full pipeline ID e.g. `142820540-26098K0076-26TULIP24-5877-M-92197814`)
- CSQ field indices (0-based): SYMBOL=1, Consequence=2, gnomADe_AF=24, gnomADg_AF=27, Prev_Count_AC=31
- gnomADg_AF, gnomADe_AF, Prev_Count_AC all exist as standalone String-typed INFO tags AND in the CSQ string
- FORMAT fields include AF (Sentieon allele fraction), DP (depth), AD, GT, F1R2, F2R1, SB
- FILTER values include: PASS, orientation, weak_evidence, clustered_events, haplotype, strand_bias, base_qual, slippage, multiallelic, map_qual (soft-filter names from Sentieon TNfilter)

---

## Git log
```
b57e0d2  Update README: Uranus pipeline context, VCF preconditions, pipeline steps
d9ff92a  Match clinical split-vep command exactly: --columns -, -a CSQ, -d
f6cd6dc  Add DP/AF quality filter to match clinical pipeline
50344df  Fix filter pipeline: norm already applied upstream, drop that step
7d4793f  Initial implementation of pairwise contamination screen
```

---

## Known limitations / future work

- Variant counts after filtering are low (10-15 per sample PASS-only). For 48 unrelated samples, contamination signal should be stronger — many more shared rare variants expected between genuinely contaminated pairs vs. the same-patient samples used for testing here.
- The `--min-shared` default of 10 and `--peak-count` default of 8 may need tuning once run on the full 48-sample cohort.
- No central exclusion zone: a peak at ratio ~1 correctly flags same-patient sample pairs and swaps (previously this would have been silently suppressed). The analyst distinguishes same-patient sharing from true contamination by reviewing sample metadata.
- `load_vcf` reads `CSQ_SYMBOL` from the filtered VCF INFO field. If a variant has no CSQ annotation (intergenic), gene will be empty string.
- plots require matplotlib (optional dependency).

---

## Files
```
~/Documents/contamination/
├── contamination_screen.py    Main script (~550 lines)
├── README.md                  Usage, design, Uranus context
└── .git/
```
