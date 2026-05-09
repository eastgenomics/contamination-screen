#!/usr/bin/env python3
"""
Phase 4: Download one run's VCFs with anonymous filenames and create TBI indices.

Run unarchive_all.py first if select_runs.py reported any archived files.
This script skips non-live files with a warning rather than waiting for them.
"""
import argparse, csv, os, subprocess, sys
import dxpy


def main():
    """Download and index VCFs for one retrospective run using anonymised IDs."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--run",           required=True,  help="e.g. RUN001")
    parser.add_argument("--mapping",       required=True,  help="Path to mapping TSV")
    parser.add_argument("--outdir",        required=True)
    parser.add_argument("--skip-existing", action="store_true")
    args = parser.parse_args()

    rows = []
    with open(args.mapping) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["anon_id"].startswith(args.run + "_"):
                rows.append(row)

    if not rows:
        print(f"No files found in mapping for {args.run}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    for row in rows:
        vcf_local = os.path.join(
            args.outdir, f"{row['anon_id']}_tnhaplotyper2_annotated.vcf.gz"
        )
        tbi_local = vcf_local + ".tbi"

        if args.skip_existing and os.path.exists(vcf_local) and os.path.exists(tbi_local):
            print(f"  Skipping (exists): {row['anon_id']}")
            continue

        # Re-check live state from DNAnexus rather than trusting cached mapping
        actual_state = dxpy.DXFile(row["file_id"], project=row["project_id"]).describe(
            fields={"archivalState": True}
        )["archivalState"]
        if actual_state != "live":
            print(f"  WARNING: {row['anon_id']} is not live "
                  f"(state={actual_state}) — run unarchive_all.py first",
                  file=sys.stderr)
            continue

        print(f"  {row['anon_id']} <- {row['real_sample_id']}")
        dxpy.download_dxfile(row["file_id"], vcf_local, project=row["project_id"])

        print(f"    Indexing...")
        subprocess.run(["bcftools", "index", "-t", vcf_local], check=True)

if __name__ == "__main__":
    main()
