#!/usr/bin/env python3
"""
dx_fetch.py — Fetch one MYE run's VCFs, plate layout, and FREEMIX data
from DNAnexus, ready for contamination_screen.py.

Downloads into a local output directory:
  <output>/vcfs/   *tnhaplotyper2_annotated.vcf.gz + .tbi for each sample
  <output>/multiqc_general_stats.txt   (plate layout)
  <output>/multiqc_verifybamid.txt     (FREEMIX / VerifyBamID)

On success, prints the ready-to-run contamination_screen.py command.

PROJECT may be:
  - A project ID:       project-BQbJpBj0bvyZqK7XG00yEv4Q
  - An exact name:      002_260423_A01303_0760_AHLTLCDRX7_MYE
  - A glob pattern:     002_260423*MYE   (must resolve to exactly one project)

Control samples (Q in specimen ID, e.g. 25357Q0020) are excluded by default.
Pass --no-exclude-controls to include them.

Exit codes:
  0  Success
  1  Error (auth failure, project not found, no VCFs, etc.)
  2  No VCF files found
"""

from __future__ import annotations

import argparse
import fnmatch
import os
import subprocess
import sys
from pathlib import Path

# -- dx-grab import -----------------------------------------------------------

_DXGRAB_DIR = Path(os.environ.get("DXGRAB_DIR", Path.home() / "Documents" / "dx_grab"))
if str(_DXGRAB_DIR) not in sys.path:
    sys.path.insert(0, str(_DXGRAB_DIR))

try:
    import dx_grab
except ImportError:
    print(
        f"ERROR: dx-grab not found at {_DXGRAB_DIR}.\n"
        "Clone it from https://github.com/eastgenomics/dx-grab",
        file=sys.stderr,
    )
    sys.exit(1)


# -- Error helper ------------------------------------------------------------

def _die(msg: str, code: int = 1) -> None:
    """Print an error message to stderr and exit."""
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


# -- Constants ----------------------------------------------------------------

_VCF_NAME_PAT      = "*tnhaplotyper2_annotated.vcf.gz"
_VCF_FOLDER_PAT    = "*/eggd_vep-*"
_PLATE_LAYOUT_NAME = "multiqc_general_stats.txt"
_FREEMIX_NAME      = "multiqc_verifybamid.txt"
_DEFAULT_EXCLUDE   = ["*Q*"]

# MultiQC folder preference — first match wins (compared case-insensitively)
_MULTIQC_FOLDER_PREFS = [
    "*eggd_multiqc*",   # eggd_MultiQC and its data subfolders
    "*multiqc_data*",   # legacy location at project root
]


# -- Argument parsing ---------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse and return command-line arguments."""
    p = argparse.ArgumentParser(
        prog="dx_fetch.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "project", metavar="PROJECT",
        help="DNAnexus project ID, exact name, or glob (e.g. '002_260423*MYE').",
    )
    p.add_argument(
        "--output", "-o", default=None, metavar="DIR",
        help="Local output directory. Default: ./<project_name>/",
    )
    p.add_argument(
        "--exclude", action="append", default=[], metavar="PATTERN",
        help="Exclude VCFs whose filename matches this glob. Repeatable. "
             "The default '*Q*' exclusion is applied separately via "
             "--no-exclude-controls and cannot be overridden with this flag.",
    )
    p.add_argument(
        "--no-exclude-controls", action="store_true",
        help="Disable the default '*Q*' exclusion for control samples.",
    )
    p.add_argument(
        "--yes", action="store_true",
        help="Automatically submit unarchive requests without prompting.",
    )
    p.add_argument(
        "--skip-archived", action="store_true",
        help="Skip archived VCFs instead of offering to unarchive them.",
    )
    p.add_argument(
        "--skip-existing", action="store_true",
        help="Skip files that already exist locally.",
    )
    p.add_argument(
        "--dry-run", action="store_true",
        help="List matched files without downloading anything.",
    )
    return p.parse_args()


# -- MultiQC file selection ---------------------------------------------------

def find_multiqc_file(dxpy, proj_dict: list[dict], filename: str) -> dict | None:
    """Find a MultiQC file in the project, preferring eggd_MultiQC folders.

    Uses dx_grab.find_files for discovery, then ranks candidates by folder
    preference: eggd_MultiQC (and nested data subfolders) > multiqc_data >
    anywhere else.
    """
    candidates = dx_grab.find_files(dxpy, proj_dict, filename, None)
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]

    def _rank(f: dict) -> int:
        """Return folder preference rank; lower is preferred."""
        folder = f["folder"].lower()
        for i, pat in enumerate(_MULTIQC_FOLDER_PREFS):
            if fnmatch.fnmatch(folder, pat):
                return i
        return len(_MULTIQC_FOLDER_PREFS)

    ranked = sorted(candidates, key=_rank)
    chosen, others = ranked[0], ranked[1:]
    print(f"  Multiple {filename!r} found; using: {chosen['folder']}/{chosen['name']}")
    for f in others:
        print(f"    (ignored) {f['folder']}/{f['name']}")
    return chosen


# -- TBI indexing -------------------------------------------------------------

def index_vcf(vcf_path: Path) -> bool:
    """Create a TBI index for a VCF.gz with bcftools index -t.

    Return True if the index already exists or was created successfully.
    Return False and print a warning if bcftools index fails.
    """
    tbi = Path(str(vcf_path) + ".tbi")
    if tbi.exists():
        return True
    result = subprocess.run(
        ["bcftools", "index", "-t", str(vcf_path)],
        capture_output=True,
    )
    if result.returncode != 0:
        print(
            f"  WARNING: bcftools index failed for {vcf_path.name}:\n"
            f"  {result.stderr.decode().strip()}",
            file=sys.stderr,
        )
        return False
    return True


# -- Output command -----------------------------------------------------------

def _print_command(outdir: Path, vcf_dir: Path,
                   plate_file: dict | None, freemix_file: dict | None,
                   dry_run: bool = False) -> None:
    """Print the ready-to-run contamination_screen.py command.

    Includes --plate-layout and --freemix-file only when the corresponding
    file dicts are non-None. Appends a NOTE for each omitted optional file.
    """
    prefix = "\n# Suggested command" + (
        " (dry run \u2014 adjust paths as needed):" if dry_run else ":"
    )
    print(prefix)
    parts = [
        f"python contamination_screen.py {vcf_dir}",
        f"    --outdir {outdir / 'results'}",
    ]
    if plate_file:
        parts.append(f"    --plate-layout {outdir / _PLATE_LAYOUT_NAME}")
    if freemix_file:
        parts.append(f"    --freemix-file {outdir / _FREEMIX_NAME}")
    print(" \\\n".join(parts))
    if not plate_file:
        print(f"# NOTE: {_PLATE_LAYOUT_NAME!r} not found \u2014 "
              "--plate-layout omitted (matrix ordered by filename).")
    if not freemix_file:
        print(f"# NOTE: {_FREEMIX_NAME!r} not found \u2014 "
              "--freemix-file omitted (FREEMIX gating disabled).")


# -- Main ---------------------------------------------------------------------

def main() -> None:
    """Authenticate, discover, download, index, and print the run command.

    Orchestrates the full fetch workflow: resolves the DNAnexus project,
    finds and filters VCF files, handles archival state, downloads VCFs and
    MultiQC files, indexes VCFs with bcftools, then prints the
    contamination_screen.py invocation.
    """
    args = parse_args()

    try:
        dxpy = dx_grab.check_auth()
    except (RuntimeError, ImportError) as e:
        _die(str(e))

    try:
        project_id, project_name = dx_grab.resolve_project(dxpy, args.project)
    except ValueError as e:  # raised by resolve_project() and find_projects() on bad project arg
        _die(str(e))

    print(f"Project: {project_name} ({project_id})")

    outdir  = Path(args.output) if args.output else Path(project_name)
    vcf_dir = outdir / "vcfs"

    # proj_dict: the single-element list dx_grab.find_files expects
    proj_dict = [{"id": project_id, "describe": {"name": project_name}}]

    # Build exclusion list
    exclude = list(args.exclude)
    if not args.no_exclude_controls:
        exclude += _DEFAULT_EXCLUDE
        print("Control samples excluded by default (--exclude '*Q*').")

    # --- Discover files ------------------------------------------------------

    print()
    vcf_files = dx_grab.find_files(dxpy, proj_dict, _VCF_NAME_PAT, _VCF_FOLDER_PAT)

    if exclude:
        before = len(vcf_files)
        vcf_files = [
            f for f in vcf_files
            if not any(fnmatch.fnmatch(f["name"].lower(), p.lower())
                       for p in exclude)
        ]
        n_excl = before - len(vcf_files)
        if n_excl:
            print(f"  Excluded {n_excl} file(s) matching: {', '.join(exclude)}")

    if not vcf_files:
        print(
            f"ERROR: No VCFs matching '{_VCF_NAME_PAT}' found in "
            f"'{_VCF_FOLDER_PAT}' in project '{project_name}'.",
            file=sys.stderr,
        )
        sys.exit(2)

    n_live    = sum(1 for f in vcf_files if f["archival_state"] == "live")
    n_nonlive = len(vcf_files) - n_live
    state_str = f"{n_live} live" + (f", {n_nonlive} archived" if n_nonlive else "")
    print(f"Found {len(vcf_files)} VCF(s)  ({state_str})")
    for f in sorted(vcf_files, key=lambda x: x["name"]):
        flag = "" if f["archival_state"] == "live" else f"  [{f['archival_state']}]"
        print(f"  {f['name']}  ({dx_grab.fmt_size(f['size'])}){flag}")

    print()
    plate_file  = find_multiqc_file(dxpy, proj_dict, _PLATE_LAYOUT_NAME)
    freemix_file = find_multiqc_file(dxpy, proj_dict, _FREEMIX_NAME)

    if plate_file:
        print(f"Plate layout: {plate_file['folder']}/{plate_file['name']}  "
              f"({dx_grab.fmt_size(plate_file['size'])})")
    else:
        print(f"Plate layout: NOT FOUND ({_PLATE_LAYOUT_NAME!r})")

    if freemix_file:
        print(f"FREEMIX data: {freemix_file['folder']}/{freemix_file['name']}  "
              f"({dx_grab.fmt_size(freemix_file['size'])})")
    else:
        print(f"FREEMIX data: NOT FOUND ({_FREEMIX_NAME!r})")

    if args.dry_run:
        print("\nDry run \u2014 nothing downloaded.")
        _print_command(outdir, vcf_dir, plate_file, freemix_file, dry_run=True)
        return

    # --- Download VCFs -------------------------------------------------------

    vcf_dir.mkdir(parents=True, exist_ok=True)
    for f in vcf_files:
        f["local_path"] = str(vcf_dir / f["name"])

    def _dl_vcfs(batch):
        """Download a batch of VCFs to vcf_dir."""
        dx_grab.download_files(dxpy, batch, str(vcf_dir),
                               skip_existing=args.skip_existing)

    print(f"\nDownloading VCFs to {vcf_dir}/")
    live_vcfs = [f for f in vcf_files if f["archival_state"] == "live"]
    if live_vcfs:
        _dl_vcfs(live_vcfs)

    vcf_files = dx_grab.handle_archives(
        dxpy, vcf_files,
        auto_yes=args.yes,
        skip_archived=args.skip_archived,
        on_live=_dl_vcfs,
    )

    # --- Index VCFs ----------------------------------------------------------

    print("\nIndexing VCFs with bcftools index -t ...")
    for f in sorted(vcf_files, key=lambda x: x["name"]):
        local = vcf_dir / f["name"]
        if local.exists():
            print(f"  {f['name']}")
            index_vcf(local)

    # --- Download MultiQC files ----------------------------------------------

    multiqc_files = [f for f in [plate_file, freemix_file] if f is not None]
    if multiqc_files:
        outdir.mkdir(parents=True, exist_ok=True)
        for f in multiqc_files:
            f["local_path"] = str(outdir / f["name"])

        def _dl_multiqc(batch):
            """Download a batch of MultiQC files to outdir."""
            dx_grab.download_files(dxpy, batch, str(outdir),
                                   skip_existing=args.skip_existing)

        live_mqc = [f for f in multiqc_files if f["archival_state"] == "live"]
        if live_mqc:
            _dl_multiqc(live_mqc)

        dx_grab.handle_archives(
            dxpy, multiqc_files,
            auto_yes=args.yes,
            skip_archived=args.skip_archived,
            on_live=_dl_multiqc,
        )

    # --- Done ----------------------------------------------------------------

    n_ready = sum(
        1 for f in vcf_files
        if (vcf_dir / f["name"]).exists()
        and Path(str(vcf_dir / f["name"]) + ".tbi").exists()
    )
    n_skipped = len(vcf_files) - n_ready
    print(f"\nDone. {n_ready} VCF(s) ready"
          + (f"  ({n_skipped} skipped/archived)" if n_skipped else "") + ".")

    # Only include MultiQC args in the command if the files are actually present
    # locally — they may have been skipped if archived and --skip-archived was set.
    plate_ready   = plate_file   if (outdir / _PLATE_LAYOUT_NAME).exists() else None
    freemix_ready = freemix_file if (outdir / _FREEMIX_NAME).exists()      else None
    _print_command(outdir, vcf_dir, plate_ready, freemix_ready)


if __name__ == "__main__":
    main()
