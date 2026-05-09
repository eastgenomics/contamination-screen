#!/usr/bin/env python3
"""Phase 6: Aggregate per-run summary.tsv files into combined results."""
import pandas as pd
from pathlib import Path

frames = []
for path in sorted(Path("retrospective/results").glob("RUN*/summary.tsv")):
    run_id = path.parent.name
    df = pd.read_csv(path, sep="\t")
    df.insert(0, "run", run_id)
    frames.append(df)

if not frames:
    print("No results found. Have any runs been screened?")
    raise SystemExit(1)

combined = pd.concat(frames, ignore_index=True)
combined.to_csv("retrospective/results/all_runs_summary.tsv", sep="\t", index=False)

flagged = combined[combined["flagged"] == True].sort_values("peak_count", ascending=False)
flagged.to_csv("retrospective/results/flagged_pairs.tsv", sep="\t", index=False)

print(f"Runs completed:          {combined['run'].nunique()}")
print(f"Total pairs assessed:    {len(combined)}")
print(f"Flagged pairs:           {len(flagged)}")
print(f"Runs with ≥1 flag:       {flagged['run'].nunique()}")

if not flagged.empty:
    print()
    print("Top 20 flagged pairs:")
    cols = ["run", "sample_a", "sample_b", "peak_count",
            "contamination_fraction", "contamination_source"]
    print(flagged[cols].head(20).to_string(index=False))
