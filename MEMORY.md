# Session Memory - contamination_screen development
# 2026-05-07

## Project
`~/Documents/contamination/` - git repo
Tool: `contamination_screen.py` - pairwise cross-sample contamination detector for the Uranus haematological oncology panel sequencing pipeline (CUH Bioinformatics).

---

## Context

Three VCFs in `~/Downloads/` from the same sequencing run (26TULIP24), **different patients**:
- `26098K0076` suffix `-M-` (142820540) - 720 variants unfiltered, 120 after quality filter
- `26093K0005` suffix `-M-` (142789636) - 452 variants unfiltered, 160 after quality filter
- `26097K0043` suffix `-F-` (142799546) - 1017 variants unfiltered, 122 after quality filter

Input VCFs are produced by the **`eggd_vep` stage** of the Uranus pipeline.

### Key finding: these are DIFFERENT patients, not serial samples
Initially misinterpreted as same-patient serial marrow samples (the 40 shared variants at ratio ~0.5 looked like clonal expansion). In fact they are three different patients on the same sequencing run. The shared variants represent genuine cross-sample contamination - one patient's germline variants leaking into another's library at ~50% contamination fraction.

### Sample ID format
`ACCESSION-LABID-RUN-PLATEPOS-TYPE-IDENTIFIER`
- 26TULIP24 = sequencing run name
- 5877 = plate position (shared across samples on same plate position??)
- M/F = sample type (not sex)
- 92197814 = a shared identifier (not patient-specific - possibly a panel or workflow ID)
- Patient identity is the lab ID: 26098K0076, 26093K0005, 26097K0043

---

## Design decisions

### Why no gnomAD/synonymous/Prev_Count filter
Cross-sample contamination consists primarily of the source patient's **germline heterozygous SNPs** (~50% VAF in source) appearing at a diluted fraction in the recipient. These variants:
- Have population-level gnomAD AFs (0.01-0.40 typically)
- Include synonymous, intronic, and coding variants equally
- May be recurrent in the Uranus haem-onc cohort (high Prev_Count_AC)

Testing against the 40 known contamination markers:
- gnomAD >= 0.002: removed **35/40** (87.5%) - DEVASTATING
- Synonymous filter: removed **19/40** (47.5%) - VERY HARMFUL
- Prev_Count_AC > 853: removed **6/40** (15%) - harmful
- DP < 99: removed **0/40** - safe to keep
- AF < 0.03: removed **0/40** - safe to keep (but limits low-contam sensitivity)

All annotation-based filters dropped. Only quality filters retained.

### Dual-peak detection
Without population filters, shared common germline variants between unrelated people cluster at ratio = 1 (both carry them at ~50%). This dominates the tallest histogram bin for every pair.

Solution: report TWO peaks:
1. **Overall peak** (tallest bin, typically ratio=1) - for swap/relatedness detection
2. **Non-unity peak** (tallest bin where |log2|>0.3) - the contamination signal

Flagging is based on the non-unity peak only. This ensures germline sharing cannot mask contamination.

### Architecture
1. Pre-filter: `bcftools view -f PASS | bcftools filter -e '(FORMAT/DP<99 || AF<0.03)'`
2. Load all filtered VCFs into memory as pandas DataFrames
3. N*(N-1)/2 pairwise comparisons (multiprocessing for >50 pairs)
4. Dual-peak log2-ratio histogram analysis per pair
5. Flag based on non-unity peak; report both peaks
6. Output: summary.tsv, matrix.tsv, flagged_pairs/*.tsv, plots/*.png

---

## Results on test data

With simplified pipeline (PASS + DP>=99 + AF>=0.03, no annotation filters):
- 120-160 variants per sample (vs 10-15 with old Uranus-mimic filters)
- **3/3 contaminated pairs flagged** (vs 0/3 previously)
- Non-unity peak: 13-18 variants at ratio 2:1 for each pair
- Direction correctly identified
- Contamination fraction correctly estimated at ~50%

---

## VCF structure notes (Uranus eggd_vep output)

- bcftools norm -m -any already applied (confirmed in header)
- CSQ field: SYMBOL at index 1 (0-based)
- FORMAT fields: AF (Sentieon allele fraction), DP, AD, GT
- FILTER values: PASS, orientation, weak_evidence, clustered_events, haplotype, strand_bias, base_qual, slippage, multiallelic, map_qual
- gnomADg_AF, gnomADe_AF, Prev_Count_AC exist as standalone String-typed INFO tags AND in CSQ

---

## Historical note: the split-vep journey

Earlier iterations tried to replicate the Uranus clinical filter pipeline exactly:
```
bcftools +split-vep --columns - -a CSQ -Ou -p 'CSQ_' -d | bcftools annotate -x INFO/CSQ
bcftools filter --soft-filter EXCLUDE -m + -e '(CSQ_Prev_Count_AC>853 || ...)'
bcftools filter -e '(FORMAT/DP<99 || AF<0.03)'
```

This required solving several bcftools issues:
- `CSQ_Prev_Count_AC` typed as String (split-vep has no built-in rule for it)
- Fix: `bcftools annotate -x INFO/CSQ_Prev_Count_AC -h fix_header.hdr` to recast as Integer
- split-vep `-c FIELDNAME` doesn't work (name collision with existing INFO tags)
- Fix: use index ranges (`-c 1-2,24,27,31`) or CSQ_ prefix (avoids conflicts)

All this became moot when we realised the clinical filters destroy the contamination signal.

---

## Git log
```
9b59b9f  Rewrite: drop gnomAD/synonymous/Prev_Count filters, add dual-peak detection
5f5eb86  Remove central exclusion zone from peak detection
e0bf0e8  Specify eggd_vep as the Uranus stage that produces the input VCFs
b57e0d2  Update README: Uranus pipeline context, VCF preconditions, pipeline steps
d9ff92a  Match clinical split-vep command exactly: --columns -, -a CSQ, -d
f6cd6dc  Add DP/AF quality filter to match clinical pipeline
50344df  Fix filter pipeline: norm already applied upstream, drop that step
7d4793f  Initial implementation of pairwise contamination screen
```

---

## Files
```
~/Documents/contamination/
├── contamination_screen.py    Main script (~500 lines)
├── README.md                  Usage, design, rationale
├── MEMORY.md                  This file
└── .git/
```
