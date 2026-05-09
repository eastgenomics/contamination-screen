#!/usr/bin/env python3
"""Phase 3: Build anonymisation mapping from selected_vcfs.json."""
import json, csv, os
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

output_path = os.path.expanduser(
    "~/Downloads/contamination_screen/contamination_screen_mapping.tsv"
)
os.makedirs(os.path.dirname(output_path), exist_ok=True)
with open(output_path, "w") as fh:
    w = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"Mapping written to {output_path}")
print(f"Total samples: {len(rows)}")
