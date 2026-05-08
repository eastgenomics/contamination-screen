# Session Memory — contamination_screen
# Last updated: 2026-05-07

## Project

`~/Documents/contamination/` — git repo, pushed to https://github.com/woook/contamination-screen (private)

Tool: `contamination_screen.py` — pairwise cross-sample contamination detector for the Uranus haematological oncology panel sequencing pipeline (CUH Bioinformatics).

---

## Key facts

### What the tool does
Detects cross-sample contamination in a cohort of tumour-only panel VCFs by finding shared variants at consistent VAF ratios. Uses dual-peak log2-ratio histogram analysis to separate contamination signal (non-unity ratio) from shared germline (ratio = 1).

### Input
VCFs from the **eggd_vep** stage of the Uranus pipeline. Must have had `bcftools norm -m -any` applied (confirmed in VCF header).

### Sample ID format (Uranus)
`ACCESSION-LABID-RUN-PLATEPOS-TYPE-IDENTIFIER`
- 26TULIP24 = sequencing run name
- 5877 = shared across all samples (panel/workflow ID, NOT patient-specific)
- M/F/U = sample type (Myeloid/FFPE/Unknown — NOT sex)
- Patient identity = lab ID (e.g. 26098K0076)
- 92197814 = shared identifier (NOT patient ID)

### The three test VCFs are DIFFERENT patients
Initially misinterpreted as same-patient serial samples. They are three unrelated patients on the same sequencing run (26TULIP24). The shared variants at ratio ~0.5 represent genuine cross-sample contamination, not clonal expansion.

---

## Filter pipeline (final)

```bash
bcftools view -f PASS input.vcf.gz -Ou \
| bcftools filter -e '(FORMAT/DP<99 || AF<0.03)' -Ou \
| bcftools +split-vep --columns - -a CSQ -p CSQ_ -s worst -Ou \
| bcftools view -e 'CSQ_gnomADg_AF>=0.40 || CSQ_gnomADe_AF>=0.40' -Oz -o out.vcf.gz
```

### Why NOT the clinical Uranus filters
The Uranus clinical filters (gnomAD >= 0.002, synonymous exclusion, Prev_Count_AC > 853) are designed for somatic variant reporting. They destroyed the contamination signal:
- gnomAD >= 0.002: removed 35/40 (87.5%) contamination markers
- Synonymous: removed 19/40 (47.5%)
- Prev_Count_AC > 853: removed 6/40 (15%)

Cross-sample contamination IS germline variants from the source patient.

### Why gnomAD >= 0.40 specifically
At gnomAD AF > 0.40, many individuals are homozygous-alt (AF~1.0) while others are het (AF~0.5). This creates a false 2:1 ratio mimicking contamination. Testing: ALL 15 false-positive markers had gnomAD > 0.40; only 1/18 real markers was above 0.40.

### Thresholds — all derived, not arbitrary
| Parameter | Value | Derivation |
|---|---|---|
| gnomAD AF | >= 0.40 | Hom/het artefact mechanism; empirically validated |
| Unity zone | \|log2\| <= 0.3 | Binomial noise at DP>=99: SE(log2)=0.20, captures ±1.5σ |
| peak_count | >= 6 | Multinomial null model: P(max>=6) < 0.001/pair, <1 FP in 1128 pairs |
| DP | >= 99 | Reliable VAF estimation |
| AF | >= 0.03 | Minimum reportable allele fraction |

### Linkage disequilibrium
Linked germline SNPs in the same gene (e.g. 5 PRPF8 variants on one haplotype) create false clusters of 4-5 at the same ratio. The peak_count >= 6 threshold absorbs this.

---

## 46-sample cohort results (run 26TULIP24)

### Confirmed contamination
**26098K0076 (B4) → 26093K0005 (A4)**: peak_count=17, ~50% contamination
- Plate-adjacent (same column)
- VerifyBamID FREEMIX: 45.5% (recipient), 13.9% (source) — bidirectional
- Recipient has lowest coverage in cohort (451x vs median 2468x)
- Consistent with physical mixing during library normalisation/pooling

### Transitive contamination problem
26093K0005 flagged as recipient in 13/19 pairs (all at fraction 0.50). These are NOT independent events — the contaminant's germline variants are also present as germline het in other patients, creating false 2:1 ratios. Diagnostic: one sample is recipient in many pairs at identical fraction.

### Other likely real events
- 26089K0042 (E1) → 26092K0050 (E2): plate-adjacent (same row), peak_count=8
- 26092K0050 (E2) → 26097K0043 (C2): same column, peak_count=10

### Other high-FREEMIX samples (from VerifyBamID)
- 26105Q0012 (A1): FREEMIX=20.7% — not strongly flagged by our tool
- 26092K0051 (H5): FREEMIX=7.0%

---

## Known limitations / future work

- **Transitive contamination**: a heavily contaminated sample flags against many others. Mitigation strategies discussed: iterative exclusion, variant uniqueness filter, one-recipient-many-sources pattern detection, graph-based minimum explanation. None implemented yet.
- **Bidirectional contamination**: tool reports dominant direction only. Both samples in the A4/B4 pair are contaminated by each other.
- **Contamination > 81%**: falls inside unity zone, detected by overall peak but not flagged as non-unity contamination. Near-complete swaps are captured differently.
- **--max-output** limits detail/plot output to top 10 flagged pairs by default.

---

## bcftools technical notes

- `bcftools +split-vep -c FIELDNAME` fails when INFO tags with same name exist (even if String-typed). Use `--columns -` with `-p CSQ_` prefix instead.
- `split-vep` `.*_AF` type rule assigns Float to gnomAD fields automatically.
- `CSQ_Prev_Count_AC` gets String (no matching rule) — needs manual reheader to Integer for arithmetic. Not needed in current pipeline (Prev_Count not used).
- Index-range notation (`-c 1-2,24,27,31`) works but name-based doesn't for these fields.

---

## Files
```
~/Documents/contamination/
├── contamination_screen.py    Main script (~480 lines)
├── README.md                  Comprehensive docs with threshold derivations
├── MEMORY.md                  This file
├── .gitignore                 Excludes data/results/pycache
└── .git/
```

GitHub: https://github.com/woook/contamination-screen (private)
