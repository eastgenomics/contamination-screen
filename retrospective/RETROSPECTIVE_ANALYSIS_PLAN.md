# Retrospective Contamination Screen: Last 100 Haemonc Runs

**Purpose:** Assess whether cross-sample contamination has occurred in the last
100 haematological oncology Uranus sequencing runs, using anonymised sample
identifiers to prevent incidental findings.

---

## Overview

**dx-grab** (`eastgenomics/dx-grab`) handles discovery and enumeration of VCFs
across DNAnexus projects. Its `--dry-run --json` mode gives a machine-readable
file manifest that we use to build the anonymisation map before any file
is downloaded. Downloads are then driven by a small wrapper that renames each
VCF to its anonymous ID on the way to disk — real sample names never appear in
any results file. TBI indices are always created locally after download.

`contamination_screen.py` runs per-run on each anonymised cohort. Results
stay anonymous throughout. De-anonymisation is only done offline, by hand, if
a specific flagged pair needs clinical follow-up.

---

## Why not use the `haem-vcf` dx-grab preset?

The `haem-vcf` preset targets **mutect2 pre-workbook VCFs** in
`*eggd_vcf_rescue*` folders. `contamination_screen.py` requires **eggd_vep
annotated VCFs** — the TNhaplotyper2-called, VEP-annotated output
(`*tnhaplotyper2_annotated.vcf.gz`) from the `eggd_vep` Uranus stage. These
contain the `CSQ` INFO field with gnomAD AF subfields (`CSQ_gnomADg_AF`,
`CSQ_gnomADe_AF`) that the filtering pipeline depends on. We use the same
project pattern (`002_2[56]*MYE`) but a different folder and filename.

---

## Prerequisites

```bash
# dx-grab
git clone https://github.com/eastgenomics/dx-grab.git
cd dx-grab
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
source .venv/bin/activate

# contamination screen dependencies (from this repo)
pip install pandas numpy matplotlib

# bcftools >= 1.14 with split-vep plugin (needed for filtering and indexing)
bcftools --version
bcftools plugin split-vep 2>/dev/null && echo "split-vep OK"

# DNAnexus auth
dx login
dx whoami
```

---

## Phase 1 — Discover All Available VCFs (Dry-Run)

Run projects follow the pattern `002_YYMMDD_*_MYE`. The year prefix means
`002_26*MYE` covers only 2026 runs. If fewer than 100 runs have occurred this
year, the search must also cover 2025 (`002_25*MYE`). The two searches are
combined and deduplicated before selecting the 100 most recent.

```bash
mkdir -p retrospective

# Search 2026 runs
python3 dx-grab.py \
    --project "002_26*MYE" \
    --folder "*eggd_vep*" \
    --name "*tnhaplotyper2_annotated.vcf.gz" \
    --exclude "*Q*" \
    --dry-run \
    --json > retrospective/vcfs_2026.json

# Count distinct runs found
python3 -c "
import json
files = json.load(open('retrospective/vcfs_2026.json'))
runs = {f['project_name'] for f in files}
print(f'2026 runs found: {len(runs)}')
"
```

If fewer than 100 runs are found in 2026, also search 2025:

```bash
python3 dx-grab.py \
    --project "002_25*MYE" \
    --folder "*eggd_vep*" \
    --name "*tnhaplotyper2_annotated.vcf.gz" \
    --exclude "*Q*" \
    --dry-run \
    --json > retrospective/vcfs_2025.json
```

Merge the two searches:

```bash
python3 -c "
import json
files_26 = json.load(open('retrospective/vcfs_2026.json'))
files_25 = json.load(open('retrospective/vcfs_2025.json')) if __import__('os').path.exists('retrospective/vcfs_2025.json') else []
all_files = files_26 + files_25
# Deduplicate by file_id
seen = set()
deduped = []
for f in all_files:
    if f['file_id'] not in seen:
        seen.add(f['file_id'])
        deduped.append(f)
json.dump(deduped, open('retrospective/all_vcfs.json', 'w'), indent=2)
runs = {f['project_name'] for f in deduped}
print(f'Total VCFs: {len(deduped)} across {len(runs)} runs')
controls = [f['name'] for f in deduped if 'Q' in f['name'].split('_tnhaplotyper2')[0]]
print(f'Control samples remaining (should be 0): {len(controls)}')
if controls: print(' ', controls)
"
```

The `--exclude "*Q*"` flag drops control samples (e.g.
`26Q98K0076_tnhaplotyper2_annotated.vcf.gz`) from both searches. Because
`all_vcfs.json` is the sole input to the mapping and download steps, controls
are excluded from everything that follows automatically.

---

## Phase 2 — Select the 100 Most Recent Runs

The project name pattern `002_YYMMDD_*` encodes the run date. Lexicographic
sort on project name = chronological order.

```python
# scripts/select_runs.py
import json
from collections import defaultdict
from pathlib import Path

files = json.load(open("retrospective/all_vcfs.json"))

by_project = defaultdict(list)
for f in files:
    by_project[f["project_name"]].append(f)

sorted_projects = sorted(by_project.keys(), reverse=True)  # newest first
if len(sorted_projects) < 100:
    print(f"WARNING: only {len(sorted_projects)} runs available (wanted 100)")
top_100 = sorted_projects[:100]

selected = []
for p in top_100:
    selected.extend(by_project[p])

json.dump(selected, open("retrospective/selected_vcfs.json", "w"), indent=2)
print(f"Selected {len(selected)} VCFs from {len(top_100)} runs")
for p in top_100:
    n_live     = sum(1 for f in by_project[p] if f["archival_state"] == "live")
    n_archived = sum(1 for f in by_project[p] if f["archival_state"] != "live")
    print(f"  {p}: {len(by_project[p])} VCFs  ({n_live} live, {n_archived} archived)")
```

```bash
python3 scripts/select_runs.py
```

Note the archived count — runs more than a few weeks old are likely fully
archived. See Phase 4 for how this is handled.

Save the run list (no sample names) as the run manifest:

```bash
python3 -c "
import json, csv
files = json.load(open('retrospective/selected_vcfs.json'))
projects = sorted({(f['project_name'], f['project_id']) for f in files})
with open('retrospective/run_manifest.tsv', 'w') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['run_index', 'project_name', 'project_id'])
    for i, (name, pid) in enumerate(projects, 1):
        w.writerow([f'{i:03d}', name, pid])
"
```

`run_manifest.tsv` contains no sample names and is safe to commit.

---

## Phase 3 — Build the Anonymisation Map

Assign each VCF a globally unique anonymous ID **before any download**.

```python
# scripts/build_mapping.py
import json, csv
from collections import defaultdict

files = json.load(open("retrospective/selected_vcfs.json"))

# Stable run order: newest run = RUN001
projects_sorted = sorted({f["project_name"] for f in files}, reverse=True)
proj_to_idx = {p: i + 1 for i, p in enumerate(projects_sorted)}

by_project = defaultdict(list)
for f in files:
    by_project[f["project_name"]].append(f)

rows = []
for proj_name, proj_files in sorted(by_project.items(), key=lambda x: proj_to_idx[x[0]]):
    run_idx = proj_to_idx[proj_name]
    for sample_idx, f in enumerate(sorted(proj_files, key=lambda x: x["name"]), start=1):
        real_sample = f["name"].replace("_tnhaplotyper2_annotated.vcf.gz", "")
        anon_id = f"RUN{run_idx:03d}_S{sample_idx:02d}"
        rows.append({
            "anon_id":        anon_id,
            "real_sample_id": real_sample,
            "file_id":        f["file_id"],
            "project_id":     f["project_id"],
            "project_name":   f["project_name"],
            "archival_state": f["archival_state"],
        })

import os
output_path = os.path.expanduser(
    "~/Downloads/contamination_screen/contamination_screen_mapping.tsv"
)
os.makedirs(os.path.dirname(output_path), exist_ok=True)
with open(output_path, "w") as fh:
    w = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"Mapping written to {output_path}")
print(f"Rows: {len(rows)}")
```

```bash
python3 scripts/build_mapping.py
```

**⚠️ Security:** `contamination_screen_mapping.tsv` links anonymous IDs to
real sample IDs (which are patient-linked). It is stored in
`~/Downloads/contamination_screen/` — outside this repository and not
committed to git. Do not copy, email, or share it alongside results.

---

## Phase 4 — Download with Anonymous Renaming

dx-grab is used for discovery only. Downloads are driven by `scripts/download_run.py`,
which writes each VCF under its anonymous name. TBI indices are never stored on
DNAnexus; they are always created locally with `bcftools index -t` after download.

### Phase 4a — Unarchive all archived files upfront (if needed)

Older runs are likely archived. `scripts/unarchive_all.py` submits unarchive
requests for **all** archived files across all selected runs in a single
batch (grouped by project, up to 1000 files per API call), then polls every 10
minutes in one loop until everything is live. Unarchiving typically takes
several hours.

```bash
# Only needed if select_runs.py reported any archived files
~/Documents/dx-grab/.venv/bin/python3 scripts/unarchive_all.py
```

Safe to re-run — already-live files are silently skipped. Ctrl+C aborts
polling but the unarchive requests remain active on DNAnexus; re-run to resume.

### Phase 4b — Download

`download_run.py` expects all files to be live. It skips any that are not with
a warning — run `unarchive_all.py` first if this happens.

```bash
# Example: download a single run manually
~/Documents/dx-grab/.venv/bin/python3 scripts/download_run.py \
    --run RUN001 \
    --mapping ~/Downloads/contamination_screen/contamination_screen_mapping.tsv \
    --outdir /tmp/RUN001 \
    --skip-existing
```

In normal operation `run_all.sh` calls `download_run.py` automatically for each run.

---

## Phase 5 — Per-Run Contamination Screen

```bash
#!/bin/bash
# scripts/run_all.sh
set -euo pipefail

MAPPING=~/Downloads/contamination_screen/contamination_screen_mapping.tsv
RESULTS="retrospective/results"
SCREEN="contamination_screen.py"

mkdir -p "$RESULTS"

for RUN in $(seq -w 001 100); do
    RUN_ID="RUN${RUN}"
    RUNDIR="$RESULTS/$RUN_ID"

    if [ -f "$RUNDIR/summary.tsv" ]; then
        echo "=== $RUN_ID already complete, skipping ==="
        continue
    fi

    echo "=== Processing $RUN_ID ==="

    TMPDIR=$(mktemp -d)
    trap "rm -rf '$TMPDIR'" EXIT

    python3 scripts/download_run.py \
        --run "$RUN_ID" \
        --mapping "$MAPPING" \
        --outdir "$TMPDIR" \
        --skip-existing

    N=$(ls "$TMPDIR"/*.vcf.gz 2>/dev/null | wc -l)
    if   [ "$N" -eq 0 ]; then
        echo "  No VCFs downloaded for $RUN_ID (all archived?), skipping"
        rm -rf "$TMPDIR"; trap - EXIT; continue
    elif [ "$N" -gt 96 ]; then PEAK_COUNT=8
    elif [ "$N" -gt 48 ]; then PEAK_COUNT=7
    else                        PEAK_COUNT=6
    fi

    echo "  $N samples -> --peak-count $PEAK_COUNT"

    python3 "$SCREEN" "$TMPDIR" \
        --outdir "$RUNDIR" \
        --vcf-glob "${RUN_ID}_S*_tnhaplotyper2_annotated.vcf.gz" \
        --peak-count "$PEAK_COUNT" \
        --threads 8 \
        --plots

    rm -rf "$TMPDIR"
    trap - EXIT

    echo "=== $RUN_ID complete ==="
done
```

```bash
bash scripts/run_all.sh 2>&1 | tee retrospective/run_all.log
```

The resume logic (checking for existing `summary.tsv`) makes the script safe
to re-run after interruption.

**Note on archived runs:** Run `unarchive_all.py` before starting `run_all.sh`
if any runs are archived. It submits all requests at once and waits; once it
completes, `run_all.sh` will find everything live.

---

## Phase 6 — Aggregate Results

```python
# scripts/aggregate.py
import pandas as pd
from pathlib import Path

frames = []
for path in sorted(Path("retrospective/results").glob("RUN*/summary.tsv")):
    run_id = path.parent.name
    df = pd.read_csv(path, sep="\t")
    df.insert(0, "run", run_id)
    frames.append(df)

combined = pd.concat(frames, ignore_index=True)
combined.to_csv("retrospective/results/all_runs_summary.tsv", sep="\t", index=False)

flagged = combined[combined["flagged"] == True].sort_values("peak_count", ascending=False)
flagged.to_csv("retrospective/results/flagged_pairs.tsv", sep="\t", index=False)

print(f"Runs completed:          {combined['run'].nunique()}")
print(f"Total pairs assessed:    {len(combined)}")
print(f"Flagged pairs:           {len(flagged)}")
print(f"Runs with ≥1 flag:       {flagged['run'].nunique()}")
print()
print("Top 20 flagged pairs:")
cols = ["run", "sample_a", "sample_b", "peak_count",
        "contamination_fraction", "contamination_source"]
print(flagged[cols].head(20).to_string(index=False))
```

All columns contain only anonymous IDs — safe to share within the team.

---

## Phase 7 — De-anonymise for Clinical Follow-Up (if needed)

Only perform this step if a specific flagged pair needs clinical or IG
follow-up. Treat any output as patient-linked data.

```python
import pandas as pd, os

mapping = pd.read_csv(
    os.path.expanduser("~/Downloads/contamination_screen/contamination_screen_mapping.tsv"),
    sep="\t"
)
anon_to_real = dict(zip(mapping["anon_id"], mapping["real_sample_id"]))
anon_to_run  = dict(zip(mapping["anon_id"], mapping["project_name"]))

flagged = pd.read_csv("retrospective/results/flagged_pairs.tsv", sep="\t")
flagged["real_source"]    = flagged["contamination_source"].map(anon_to_real)
flagged["real_recipient"] = flagged["contamination_recipient"].map(anon_to_real)
flagged["dnanexus_run"]   = flagged["contamination_source"].map(anon_to_run)

# ⚠ Patient-linked — handle under lab IG policy
flagged.to_csv("retrospective/results/flagged_pairs_IDENTIFIED.tsv", sep="\t", index=False)
```

---

## Interpretation Guide

| Signal | Interpretation |
|---|---|
| `flagged=True`, `peak_count ≥ 6` | Suspected contamination at the implied fraction — review detail TSV and histogram |
| One sample is recipient in many pairs, same fraction, never source | Transitive artefact: primary event is the pair with highest `peak_count` |
| Large `overall_count` at ratio ≈ 1 | Possible sample swap or duplicate |
| `peak_count` 4–5 | Borderline — check if variants span multiple chromosomes (contamination) or cluster in one gene (LD artefact, e.g. PRPF8) |

---

## Directory Structure

```text
contamination-screen/
├── contamination_screen.py
├── scripts/
│   ├── select_runs.py              # Phase 2
│   ├── build_mapping.py            # Phase 3
│   ├── download_run.py             # Phase 4
│   ├── run_all.sh                  # Phase 5
│   └── aggregate.py                # Phase 6
├── retrospective/
│   ├── vcfs_2026.json              # dx-grab output, 2026 runs
│   ├── vcfs_2025.json              # dx-grab output, 2025 runs (if needed)
│   ├── all_vcfs.json               # merged, deduplicated
│   ├── selected_vcfs.json          # top 100 runs
│   ├── run_manifest.tsv            # run index → project (no sample names) ✓ safe to commit
│   ├── run_all.log
│   └── results/
│       ├── RUN001/
│       │   ├── summary.tsv         # anonymous IDs only
│       │   ├── matrix.tsv
│       │   ├── flagged_pairs/
│       │   └── plots/
│       ├── RUN002/ ...
│       ├── all_runs_summary.tsv    # combined, anonymous IDs only ✓ safe to share
│       └── flagged_pairs.tsv       # anonymous IDs only ✓ safe to share

~/Downloads/contamination_screen/
└── contamination_screen_mapping.tsv   # real sample IDs — NOT in repo ⚠
```

---

## Estimated Cost

| Item | Estimate |
|---|---|
| dx-grab setup + dry-run | ~30 min |
| Phase 2–3 (selection + mapping) | ~30 min |
| Phase 4–5 (download + screen, 100 runs, live only) | ~4–8 hours wall-clock |
| Unarchiving delay (if older runs needed) | Several hours per batch |
| DNAnexus egress (~100 runs × 48 VCFs × 50 MB) | ~240 GB ≈ €24 |
