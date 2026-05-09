#!/usr/bin/env python3
"""
Pre-flight unarchiving: submit unarchive requests for all archived files
across all selected runs in parallel, then poll until everything is live.

Run this before run_all.sh if select_runs.py reported any archived files.
Safe to re-run — already-live files are silently skipped.
"""
import json, sys, time
from collections import defaultdict
import dxpy
from dxpy.exceptions import DXAPIError

POLL_INTERVAL = 600  # 10 minutes
BATCH_SIZE    = 1000  # DNAnexus API limit per unarchive call


def main():
    """Submit unarchive requests for all archived VCFs and poll until live.

    Reads retrospective/selected_vcfs.json, submits unarchive requests
    batched by project, then polls every 10 minutes until all files are
    live. Safe to re-run — already-live files are silently skipped.
    """
    files = json.load(open("retrospective/selected_vcfs.json"))

    archived = [f for f in files if f["archival_state"] != "live"]
    if not archived:
        print("All files are live — nothing to unarchive.")
        return

    print(f"{len(archived)} archived file(s) across "
          f"{len({f['project_id'] for f in archived})} project(s).")

    # --- Submit all unarchive requests, batched by project ---
    by_project = defaultdict(list)
    for f in archived:
        by_project[f["project_id"]].append(f)

    submitted_ids = set()
    failed_ids    = set()
    for proj_id, proj_files in by_project.items():
        proj_name = proj_files[0]["project_name"]
        file_ids  = [f["file_id"] for f in proj_files]
        print(f"  {proj_name}: submitting {len(file_ids)} unarchive request(s)...")
        for i in range(0, len(file_ids), BATCH_SIZE):
            batch = file_ids[i : i + BATCH_SIZE]
            try:
                dxpy.api.project_unarchive(proj_id, {"files": batch})
                submitted_ids.update(batch)
            except DXAPIError as e:
                failed_ids.update(batch)
                print(f"  WARNING: unarchive failed for {proj_name}: {e}", file=sys.stderr)

    if failed_ids:
        print(
            f"\n{len(failed_ids)} file(s) failed to submit for unarchiving — "
            f"see warnings above. Re-run this script to retry.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"\nUnarchive requested for {len(submitted_ids)} file(s). "
          f"Unarchiving typically takes several hours.")
    print("Polling every 10 minutes — Ctrl+C to abort "
          "(requests remain active; re-run to resume polling).\n")

    # --- Poll all files in a single loop ---
    waiting = {
        f["file_id"]: f
        for f in archived
        if f["file_id"] in submitted_ids
    }

    try:
        while waiting:
            time.sleep(POLL_INTERVAL)

            newly_live = []
            still_waiting = {}
            for fid, f in waiting.items():
                try:
                    state = dxpy.DXFile(fid, project=f["project_id"]).describe(
                        fields={"archivalState": True}
                    )["archivalState"]
                except DXAPIError as e:
                    print(f"  WARNING: could not check {fid}: {e}", file=sys.stderr)
                    still_waiting[fid] = f
                    continue

                if state == "live":
                    newly_live.append(fid)
                else:
                    still_waiting[fid] = f

            ready = len(submitted_ids) - len(still_waiting)
            total = len(submitted_ids)
            ts    = time.strftime("%H:%M")
            print(f"[{ts}] {ready}/{total} files live "
                  f"({len(still_waiting)} still unarchiving)...")

            waiting = still_waiting

    except KeyboardInterrupt:
        print("\n\nAborted. Unarchiving continues on DNAnexus. "
              "Re-run this script to resume polling.")
        sys.exit(0)

    print(f"\nAll {len(submitted_ids)} file(s) are now live. "
          f"Ready to run run_all.sh.")


if __name__ == "__main__":
    main()
