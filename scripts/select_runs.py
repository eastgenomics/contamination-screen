#!/usr/bin/env python3
"""Phase 2: Select the 100 most recent runs from all_vcfs.json."""
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
