#!/usr/bin/env python3
"""
contamination_screen.py

Pairwise cross-sample contamination detector for tumour sequencing panel VCFs.

Screens a cohort of VCFs for cross-sample contamination by detecting shared
variants at consistent VAF ratios. Uses a dual-peak approach: reports both the
dominant ratio peak (typically ratio ~1 from shared germline) and the dominant
non-unity peak (the contamination signal).

Designed for Uranus pipeline (CUH Bioinformatics) haematological oncology panel
VCFs output by the eggd_vep stage. Input VCFs must already have been normalised
with bcftools norm -m -any.
"""

import argparse
import logging
import math
import os
import re
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# -- Constants ----------------------------------------------------------------

# Histogram parameters
_LOG2_RANGE = 4.0   # +/-4 log2 units spans ratios from 1:16 to 16:1
_BIN_WIDTH  = 0.2   # bin width in log2 units

# Default glob for VCF files
VCF_GLOB = "*_annotated.vcf.gz"

# Default gnomAD AF threshold: remove very common variants where hom/het
# differences between individuals create false 2:1 ratio signals
_MAX_GNOMAD_AF = 0.40

# Gene symbol is extracted by bcftools +split-vep as CSQ_SYMBOL (see filter_vcf)


# -- Argument parsing ---------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="contamination_screen.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "vcf_dir", type=Path,
        help="Directory containing annotated VCF.gz files",
    )
    p.add_argument(
        "--outdir", "-o", type=Path, default=Path("results"),
        help="Output directory (default: %(default)s)",
    )
    p.add_argument(
        "--min-af", type=float, default=0.03,
        help="Minimum VAF in both samples for a variant to be informative "
             "(default: %(default)s)",
    )
    p.add_argument(
        "--min-dp", type=int, default=99,
        help="Minimum read depth (FORMAT/DP) to retain a variant. The filter "
             "expression is DP < min_dp, so DP=99 is retained at default. "
             "This matches the clinical pipeline expression 'FORMAT/DP<99'. "
             "(default: %(default)s)",
    )
    p.add_argument(
        "--max-gnomad", type=float, default=_MAX_GNOMAD_AF,
        help="Exclude variants with gnomAD AF >= this value. Very common "
             "variants (AF > 0.40) produce false 2:1 ratio signals when one "
             "sample is hom-alt and the other is het. "
             "Set to 1.0 to disable. (default: %(default)s)",
    )
    p.add_argument(
        "--min-shared", type=int, default=10,
        help="Minimum number of informative shared variants required before "
             "assessing a pair (default: %(default)s)",
    )
    p.add_argument(
        "--peak-count", type=int, default=10,
        help="Flag a pair if the non-unity peak contains >= this many "
             "variants (default: %(default)s)",
    )
    p.add_argument(
        "--n-shared-z", type=float, default=2.0,
        help="Flag a pair only if the global z-score of n_shared exceeds this "
             "threshold (computed across all pairs in the run). Combined with "
             "--peak-count: both conditions must be met. (default: %(default)s)",
    )
    p.add_argument(
        "--threads", "-t", type=int, default=min(8, os.cpu_count() or 1),
        help="Parallel threads for pairwise comparisons (default: %(default)s)",
    )
    p.add_argument(
        "--include-non-pass", action="store_true",
        help="Include non-PASS variants (default: PASS only)",
    )
    p.add_argument(
        "--plots", action="store_true",
        help="Generate log2-ratio histogram plots for flagged pairs "
             "(requires matplotlib)",
    )
    p.add_argument(
        "--max-output", type=int, default=10,
        help="Maximum number of flagged pairs to write detail TSVs and plots "
             "for, ranked by peak_count (default: %(default)s). "
             "Set to 0 for no limit.",
    )
    p.add_argument(
        "--plate-layout", type=Path, default=None,
        help="MultiQC general stats file (multiqc_general_stats.txt) to order "
             "the matrix by plate position. Makes adjacent-well contamination "
             "patterns visible. If not provided, samples are ordered by filename.",
    )
    p.add_argument(
        "--vcf-glob", default=VCF_GLOB,
        help="Glob pattern for VCF files in vcf_dir (default: %(default)s)",
    )
    p.add_argument(
        "--bin-width", type=float, default=_BIN_WIDTH,
        help="Log2-ratio histogram bin width (default: %(default)s)",
    )
    p.add_argument(
        "--force-refilter", action="store_true",
        help="Re-run filtering even if filtered VCFs already exist",
    )
    p.add_argument(
        "--verbose", "-v", action="store_true",
        help="Enable debug logging",
    )
    p.add_argument(
        "--freemix-file", type=Path, default=None,
        metavar="FILE",
        help="Path to a MultiQC VerifyBamID output file (multiqc_verifybamid.txt) "
             "containing per-sample FREEMIX values. When provided, a pair is only "
             "flagged if the recipient's FREEMIX >= --freemix-threshold. "
             "The recipient_freemix column is always written to summary.tsv "
             "when this file is supplied, regardless of threshold.",
    )
    p.add_argument(
        "--freemix-threshold", type=float, default=0.15,
        metavar="FRAC",
        help="Minimum recipient FREEMIX fraction (0–1) required for flagging "
             "when --freemix-file is provided (default: %(default)s = 15%%).",
    )
    return p.parse_args()


# -- bcftools helpers ---------------------------------------------------------

def _run(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run a subprocess, logging the command at DEBUG level."""
    logging.debug("RUN: %s", " ".join(str(c) for c in cmd))
    return subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=check,
    )


def get_sample_name(vcf: Path) -> str:
    """Return the sample column name from the VCF #CHROM header line."""
    r = _run(["bcftools", "view", "-h", str(vcf)], check=False)
    if r.returncode != 0:
        logging.warning(
            "Could not read header from %s; using stem '%s'",
            vcf.name, vcf.stem,
        )
        return vcf.stem
    for line in r.stdout.decode().splitlines():
        if line.startswith("#CHROM"):
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 9:
                return parts[9]
    logging.warning(
        "Could not parse sample name from %s; using stem '%s'",
        vcf.name, vcf.stem,
    )
    return vcf.stem


def _drain_stderr(stream, collector: list) -> None:
    """Read all stderr from a subprocess pipe into collector list."""
    collector.append(stream.read())


def _csq_gnomad_fields(src: Path) -> List[str]:
    """
    Return which gnomAD AF subfields are present in the CSQ format string
    of this VCF. Runs bcftools +split-vep --list and filters for known names.
    """
    result = subprocess.run(
        ["bcftools", "+split-vep", "-l", str(src)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        return []
    available = {line.split("\t")[1] for line in result.stdout.strip().splitlines()
                 if "\t" in line}
    return [f for f in ("gnomADg_AF", "gnomADe_AF") if f in available]


def filter_vcf(src: Path, dst: Path, pass_only: bool, min_dp: int,
               min_af: float, max_gnomad: float) -> None:
    """
    Apply quality filters for contamination detection.

    Pipeline:

      1. bcftools view [-f PASS]
             Optionally restrict to PASS variants (Sentieon TNfilter).

      2. bcftools filter -e '(FORMAT/DP < min_dp || AF < min_af)'
             Remove low-depth and low-VAF variants where allele fractions are
             unreliable.

      3. bcftools +split-vep --columns - -a CSQ -p CSQ_ -s worst
             Extract CSQ subfields into CSQ_-prefixed INFO tags (worst
             transcript only). The .*_AF built-in type rule assigns Float to
             gnomAD AF fields, enabling arithmetic filtering in the next step.

      4. bcftools view -e 'CSQ_gnomADg_AF>=T || CSQ_gnomADe_AF>=T'
             Remove very common variants (gnomAD AF >= threshold). Only fields
             present in the VCF's CSQ format string are used; if neither
             gnomADg_AF nor gnomADe_AF is present this step is skipped with a
             warning (older VEP annotations may lack gnomAD genomes AF).
    """
    cmd1 = ["bcftools", "view", str(src)]
    if pass_only:
        cmd1 += ["-f", "PASS"]
    cmd1 += ["-Ou"]

    cmd2 = [
        "bcftools", "filter",
        "-e", f"(FORMAT/DP<{min_dp} || AF<{min_af})",
        "-Ou",
    ]

    cmd3 = [
        "bcftools", "+split-vep",
        "--columns", "-",
        "-a", "CSQ",
        "-p", "CSQ_",
        "-s", "worst",
        "-Ou",
    ]

    cmd4 = ["bcftools", "view"]
    if max_gnomad < 1.0:
        gnomad_fields = _csq_gnomad_fields(src)
        if gnomad_fields:
            clause = " || ".join(f"CSQ_{f}>={max_gnomad}" for f in gnomad_fields)
            cmd4 += ["-e", clause]
            logging.debug("gnomAD filter using: %s", clause)
        else:
            logging.warning(
                "%s: no gnomAD AF fields found in CSQ — skipping gnomAD filter",
                src.name,
            )
    cmd4 += ["-Oz", "-o", str(dst)]

    logging.debug("Filter pipeline for %s", src.name)
    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    p1.stdout.close()
    p3 = subprocess.Popen(cmd3, stdin=p2.stdout, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    p2.stdout.close()
    p4 = subprocess.Popen(cmd4, stdin=p3.stdout, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    p3.stdout.close()

    # Drain stderr from p1-p3 in background threads to prevent deadlock
    # if any process writes >64KB to stderr
    stderr_collectors: Dict[subprocess.Popen, list] = {p: [] for p in (p1, p2, p3)}
    drain_threads = [
        threading.Thread(target=_drain_stderr, args=(p.stderr, stderr_collectors[p]))
        for p in (p1, p2, p3)
    ]
    for t in drain_threads:
        t.start()

    _, p4_err = p4.communicate()
    for t in drain_threads:
        t.join()
    p1.wait(); p2.wait(); p3.wait()

    # Check returncodes BEFORE indexing
    for proc, name in [(p1, "view"), (p2, "filter/DP+AF"),
                       (p3, "split-vep"), (p4, "view/gnomAD")]:
        if proc.returncode != 0:
            stderr = b"".join(stderr_collectors.get(proc, []))
            raise RuntimeError(
                f"bcftools {name} failed for {src.name}:\n"
                f"{stderr.decode()}\n{p4_err.decode()}"
            )

    # Only index after all steps confirmed successful
    _run(["bcftools", "index", "-t", str(dst)])


# -- VCF loading --------------------------------------------------------------

def load_vcf(vcf: Path) -> pd.DataFrame:
    """
    Load a filtered VCF into a DataFrame.

    The filtered VCF has CSQ subfields expanded by split-vep into CSQ_-prefixed
    INFO tags. Gene symbol is read from CSQ_SYMBOL directly.

    Returns columns: chrom, pos, ref, alt, filter_col, af, gene
    """
    r = _run([
        "bcftools", "query",
        "-f", "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%AF]\t%INFO/CSQ_SYMBOL\n",
        str(vcf),
    ], check=False)

    if r.returncode != 0:
        stderr = r.stderr.decode()
        if "CSQ" in stderr or "not defined" in stderr:
            logging.error(
                "Failed to query %s: CSQ_SYMBOL not found. "
                "Re-run with --force-refilter to regenerate filtered VCFs.",
                vcf.name,
            )
        raise RuntimeError(f"bcftools query failed for {vcf.name}: {stderr}")

    rows = []
    for line in r.stdout.decode().splitlines():
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        chrom, pos_s, ref, alt, filt, af_s = parts[:6]
        gene = parts[6].strip() if len(parts) > 6 else ""
        if gene == ".":
            gene = ""

        try:
            af = float(af_s)
        except ValueError:
            continue

        rows.append({
            "chrom": chrom, "pos": int(pos_s),
            "ref": ref, "alt": alt,
            "filter_col": filt, "af": af, "gene": gene,
        })

    if not rows:
        return pd.DataFrame(columns=["chrom", "pos", "ref", "alt",
                                     "filter_col", "af", "gene"])
    df = pd.DataFrame(rows)
    df["pos"] = df["pos"].astype(int)
    # Deduplicate in case of multi-transcript records
    df = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt"], keep="first")
    return df


# -- Core comparison ----------------------------------------------------------

def _find_peaks(
    log2_ratios: np.ndarray,
    bin_width: float,
) -> dict:
    """
    Find peaks in the log2-ratio histogram using a dual-peak approach.

    Returns both the overall tallest peak AND the tallest peak in the non-unity
    region (|log2_ratio| > 0.3). This separates:
      - Shared germline variants (ratio ~1, log2 ~0) from
      - Contamination signal (consistent non-unity ratio)

    The non-unity peak is used for contamination flagging. The unity peak is
    reported for completeness (high count = potential swap or relatedness).

    Returns dict with keys:
        overall_log2, overall_count, overall_fraction
        nonunity_log2, nonunity_count, nonunity_fraction
    """
    result = {
        "overall_log2": float("nan"), "overall_count": 0, "overall_fraction": 0.0,
        "nonunity_log2": float("nan"), "nonunity_count": 0, "nonunity_fraction": 0.0,
    }

    if len(log2_ratios) == 0:
        return result

    edges = np.arange(
        -_LOG2_RANGE - bin_width / 2,
        _LOG2_RANGE + bin_width,
        bin_width,
    )
    counts, edges = np.histogram(log2_ratios, bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2.0
    n_total = len(log2_ratios)

    if counts.max() == 0:
        return result

    # Overall peak (tallest bin across full range)
    idx_overall = int(np.argmax(counts))
    result["overall_log2"] = float(centers[idx_overall])
    result["overall_count"] = int(counts[idx_overall])
    result["overall_fraction"] = float(counts[idx_overall]) / n_total

    # Non-unity peak (tallest bin where |log2| > 0.3, i.e. ratio outside 0.81-1.23)
    nonunity_mask = np.abs(centers) > 0.3
    nonunity_counts = counts.copy()
    nonunity_counts[~nonunity_mask] = 0

    if nonunity_counts.max() > 0:
        idx_nonunity = int(np.argmax(nonunity_counts))
        result["nonunity_log2"] = float(centers[idx_nonunity])
        result["nonunity_count"] = int(counts[idx_nonunity])
        result["nonunity_fraction"] = float(counts[idx_nonunity]) / n_total

    return result


def _compare(
    name_a: str, df_a: pd.DataFrame,
    name_b: str, df_b: pd.DataFrame,
    bin_width: float,
    include_shared_df: bool = False,
) -> Tuple[dict, Optional[pd.DataFrame]]:
    """
    Core pairwise comparison with dual-peak detection.

    Merges two variant DataFrames on (chrom, pos, ref, alt), computes
    log2(AF_B / AF_A), and finds both the overall peak and the non-unity
    peak (contamination signal).

    All variants have already passed the AF floor during pre-filtering,
    so no additional VAF threshold is applied here.

    Directionality (from non-unity peak):
        peak_log2_ratio < 0  =>  AF_A > AF_B  =>  A is the contamination source
        peak_log2_ratio > 0  =>  AF_B > AF_A  =>  B is the contamination source
    """
    merged = pd.merge(
        df_a[["chrom", "pos", "ref", "alt", "filter_col", "af", "gene"]].rename(
            columns={"filter_col": "filter_a", "af": "af_a", "gene": "gene"}),
        df_b[["chrom", "pos", "ref", "alt", "filter_col", "af"]].rename(
            columns={"filter_col": "filter_b", "af": "af_b"}),
        on=["chrom", "pos", "ref", "alt"],
        how="inner",
    )

    n_shared = len(merged)

    nan = float("nan")
    result: dict = {
        "sample_a": name_a, "sample_b": name_b,
        "n_shared": n_shared,
        # Overall peak (includes ratio=1 from shared germline)
        "overall_log2": nan, "overall_ratio": nan,
        "overall_count": 0, "overall_fraction": 0.0,
        # Non-unity peak (the contamination signal)
        "peak_log2_ratio": nan, "peak_ratio": nan,
        "peak_count": 0, "peak_fraction": 0.0,
        # Directionality (derived from non-unity peak)
        "contamination_source": "", "contamination_recipient": "",
        "contamination_fraction": nan,
        "flagged": False,
    }

    if n_shared < 2:
        return result, None

    # Filter out variants with AF=0 in either sample (would produce
    # infinite/undefined ratios and inflate n_shared without contributing signal)
    merged = merged[(merged["af_a"] > 0) & (merged["af_b"] > 0)].copy()
    n_valid = len(merged)
    if n_valid < 2:
        return result, None
    if n_valid < n_shared:
        logging.debug("  %d variant(s) with AF=0 excluded from ratio histogram",
                      n_shared - n_valid)
    result["n_shared"] = n_valid

    ratios = merged["af_b"] / merged["af_a"]
    merged["ratio"] = ratios.values
    merged["log2_ratio"] = np.log2(ratios.values)

    peaks = _find_peaks(merged["log2_ratio"].values, bin_width)

    # Update overall peak
    result.update({
        "overall_log2": peaks["overall_log2"],
        "overall_ratio": (
            2.0 ** peaks["overall_log2"]
            if not math.isnan(peaks["overall_log2"]) else nan
        ),
        "overall_count": peaks["overall_count"],
        "overall_fraction": peaks["overall_fraction"],
    })

    # Update non-unity peak (contamination signal)
    result.update({
        "peak_log2_ratio": peaks["nonunity_log2"],
        "peak_ratio": (
            2.0 ** peaks["nonunity_log2"]
            if not math.isnan(peaks["nonunity_log2"]) else nan
        ),
        "peak_count": peaks["nonunity_count"],
        "peak_fraction": peaks["nonunity_fraction"],
    })

    # Mark variants in non-unity peak bin
    plr = peaks["nonunity_log2"]
    if not math.isnan(plr) and peaks["nonunity_count"] > 0:
        half = bin_width / 2.0
        merged["in_peak"] = (
            (merged["log2_ratio"] >= plr - half) &
            (merged["log2_ratio"] <= plr + half)
        )
        # Directionality
        if plr < 0:
            result.update({
                "contamination_source": name_a,
                "contamination_recipient": name_b,
                "contamination_fraction": 2.0 ** plr,
            })
        else:
            result.update({
                "contamination_source": name_b,
                "contamination_recipient": name_a,
                "contamination_fraction": 2.0 ** (-plr),
            })
    else:
        merged["in_peak"] = False

    return result, (merged if include_shared_df else None)


def _worker(args_tuple: tuple) -> dict:
    """Thread worker. Returns result dict only."""
    name_a, name_b, df_cache, bin_width = args_tuple
    result, _ = _compare(name_a, df_cache[name_a], name_b, df_cache[name_b],
                         bin_width, include_shared_df=False)
    return result


# -- Output writers -----------------------------------------------------------

def _load_freemix(freemix_file: Path) -> dict:
    """Load a MultiQC multiqc_verifybamid.txt file.

    Returns a dict mapping sample name -> FREEMIX fraction (0–1).
    Handles both raw VerifyBamID output and MultiQC-aggregated format;
    expects columns 'Sample' and 'FREEMIX'.
    """
    df = pd.read_csv(freemix_file, sep="\t")
    if "Sample" not in df.columns or "FREEMIX" not in df.columns:
        raise ValueError(
            f"--freemix-file {freemix_file}: expected columns 'Sample' and 'FREEMIX', "
            f"got {list(df.columns)}"
        )
    df = df.drop_duplicates("Sample")
    logging.info("Loaded FREEMIX data for %d samples from %s",
                 len(df), freemix_file.name)
    return dict(zip(df["Sample"], df["FREEMIX"]))


def _apply_flags(
    results: List[dict],
    min_shared: int,
    peak_count_thresh: int,
    n_shared_z_thresh: float,
    freemix: Optional[dict] = None,
    freemix_threshold: float = 0.15,
) -> pd.DataFrame:
    """Convert results to DataFrame and flag pairs based on NON-UNITY peak.

    Flagging requires all of:
      - n_shared >= min_shared
      - peak_count >= peak_count_thresh
      - n_shared_z >= n_shared_z_thresh  (within-run z-score of n_shared)
      - recipient FREEMIX >= freemix_threshold  (only when freemix dict supplied)
    """
    df = pd.DataFrame(results)
    # Within-run z-score of n_shared across all pairs in this cohort
    n_mean = df["n_shared"].mean()
    n_std  = df["n_shared"].std()
    df["n_shared_z"] = (df["n_shared"] - n_mean) / n_std if n_std > 0 else 0.0

    base_flag = (
        (df["n_shared"] >= min_shared) &
        (df["peak_count"] >= peak_count_thresh) &
        (df["n_shared_z"] >= n_shared_z_thresh)
    )

    if freemix is not None:
        # Add recipient FREEMIX column (fraction 0–1) for all pairs
        df["recipient_freemix"] = df["contamination_recipient"].map(freemix)
        freemix_flag = df["recipient_freemix"] >= freemix_threshold
        n_missing = df.loc[base_flag, "recipient_freemix"].isna().sum()
        if n_missing:
            logging.warning(
                "%d pair(s) pass peak/z thresholds but have no FREEMIX data "
                "for their recipient — FREEMIX condition treated as not met.",
                n_missing,
            )
        df["flagged"] = base_flag & freemix_flag.fillna(False)
        logging.info(
            "FREEMIX filter applied (threshold >= %.0f%%): "
            "%d/%d base-flagged pairs retained",
            freemix_threshold * 100,
            df["flagged"].sum(), base_flag.sum(),
        )
    else:
        df["flagged"] = base_flag

    for col in ["overall_log2", "overall_ratio", "overall_fraction",
                "peak_log2_ratio", "peak_ratio",
                "peak_fraction", "contamination_fraction"]:
        if col in df.columns:
            df[col] = df[col].round(4)
    return df


def write_summary(df: pd.DataFrame, outdir: Path) -> None:
    out = outdir / "summary.tsv"
    df.to_csv(out, sep="\t", index=False, na_rep="NA")
    n = df["flagged"].sum()
    logging.info("Summary written: %s  (%d/%d pairs flagged)", out, n, len(df))


def _get_plate_order(plate_file: Path, sample_names: List[str]) -> List[str]:
    """
    Parse plate positions from a MultiQC general stats file and return
    sample names ordered by plate position (column-major: A1, B1...H1, A2...).

    Samples not found in the plate file are appended at the end in their
    original order.
    """
    if not plate_file.exists():
        logging.error("Plate layout file not found: %s", plate_file)
        sys.exit(1)

    gs = pd.read_csv(plate_file, sep="\t")
    col_row = "custom_content_samplesheet_wells-well_row"
    col_col = "custom_content_samplesheet_wells-well_column"

    if col_row not in gs.columns or col_col not in gs.columns:
        logging.warning("Plate layout columns not found in %s; using file order",
                        plate_file.name)
        return sample_names

    # Drop per-lane FASTQ rows (e.g. Sample_S1_L001) — these have no well
    # position and are a MultiQC-specific row type. We do not filter by
    # sample name format so the function works with any naming scheme.
    gs = gs[~gs["Sample"].str.contains(r"_S\d+_L\d+", na=False)]
    gs = gs.dropna(subset=[col_row, col_col])

    # Build position map: sample_name -> (column, row_letter)
    pos_map = {}
    for _, row in gs.iterrows():
        name = row["Sample"]
        if name in sample_names:
            pos_map[name] = (int(row[col_col]), row[col_row])

    # Sort by column then row
    ordered = sorted(
        [n for n in sample_names if n in pos_map],
        key=lambda n: (pos_map[n][0], pos_map[n][1]),
    )
    # Append any samples not found in plate file
    remaining = [n for n in sample_names if n not in pos_map]
    if remaining:
        logging.warning("%d sample(s) not found in plate layout; appended at end",
                        len(remaining))
    return ordered + remaining


def write_matrix(
    results: List[dict],
    sample_names: List[str],
    outdir: Path,
) -> None:
    """Write directional N x N contamination matrix (non-unity peak counts)."""
    idx = {name: i for i, name in enumerate(sample_names)}
    n = len(sample_names)
    matrix = np.zeros((n, n), dtype=int)

    for r in results:
        src = r.get("contamination_source", "")
        rec = r.get("contamination_recipient", "")
        if src and rec and src in idx and rec in idx:
            matrix[idx[rec], idx[src]] = r["peak_count"]

    df = pd.DataFrame(matrix, index=sample_names, columns=sample_names)
    df.index.name = "recipient \\ source"
    out = outdir / "matrix.tsv"
    df.to_csv(out, sep="\t")
    logging.info("Matrix written: %s", out)


def write_flagged_details(
    flagged: List[dict],
    shared_cache: Dict[tuple, pd.DataFrame],
    outdir: Path,
) -> None:
    """Write per-variant detail TSV for each flagged pair."""
    detail_dir = outdir / "flagged_pairs"
    detail_dir.mkdir(exist_ok=True)

    for r in flagged:
        na, nb = r["sample_a"], r["sample_b"]
        shared = shared_cache.get((na, nb))
        if shared is None or shared.empty:
            continue

        label_a = re.sub(r"[^A-Za-z0-9_\-]", "_", na)[:50]
        label_b = re.sub(r"[^A-Za-z0-9_\-]", "_", nb)[:50]
        out = detail_dir / f"{label_a}__vs__{label_b}.tsv"

        want = ["chrom", "pos", "ref", "alt",
                "filter_a", "filter_b",
                "af_a", "af_b", "ratio", "log2_ratio", "in_peak", "gene"]
        cols = [c for c in want if c in shared.columns]
        shared[cols].sort_values("log2_ratio").to_csv(out, sep="\t", index=False)
        logging.info("  Detail: %s", out.name)


def write_plots(
    flagged: List[dict],
    shared_cache: Dict[tuple, pd.DataFrame],
    outdir: Path,
    bin_width: float,
) -> None:
    """Generate log2-ratio histogram for each flagged pair."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logging.warning("matplotlib not available; skipping plots")
        return

    plot_dir = outdir / "plots"
    plot_dir.mkdir(exist_ok=True)

    bins = np.arange(-_LOG2_RANGE, _LOG2_RANGE + bin_width, bin_width)

    for r in flagged:
        na, nb = r["sample_a"], r["sample_b"]
        shared = shared_cache.get((na, nb))
        if shared is None or shared.empty or "log2_ratio" not in shared.columns:
            continue

        fig, ax = plt.subplots(figsize=(10, 4))
        ax.hist(shared["log2_ratio"], bins=bins,
                color="steelblue", edgecolor="white", linewidth=0.4, alpha=0.85)

        # Reference ratio lines (using axes fraction for y to avoid tight_layout issues)
        for lr, lbl in [(-1.0, "2:1"), (-2.0, "4:1"), (1.0, "1:2"), (2.0, "1:4")]:
            ax.axvline(lr, color="grey", linestyle=":", linewidth=0.8, alpha=0.6)
            ax.annotate(lbl, xy=(lr, 0.92),
                        xycoords=("data", "axes fraction"),
                        ha="center", fontsize=7, color="grey")

        # Non-unity peak line (contamination)
        plr = r.get("peak_log2_ratio", float("nan"))
        if not math.isnan(plr):
            src = r.get("contamination_source", "?")[:35]
            cf = r.get("contamination_fraction", float("nan"))
            ax.axvline(plr, color="crimson", linestyle="--", linewidth=1.8,
                       label=(
                           f"Contamination peak: log2 = {plr:.2f} "
                           f"(ratio = {2**plr:.3f})\n"
                           f"Source: {src}\n"
                           f"n = {r['peak_count']},  "
                           f"fraction = {r['peak_fraction']:.2f},  "
                           f"contam ~ {cf:.1%}"
                       ))

        # Overall peak line (likely germline sharing)
        olr = r.get("overall_log2", float("nan"))
        if not math.isnan(olr) and abs(olr - (plr if not math.isnan(plr) else 99)) > 0.1:
            ax.axvline(olr, color="orange", linestyle="-.", linewidth=1.2,
                       label=f"Overall peak: log2 = {olr:.2f} "
                             f"(n={r['overall_count']})")

        ax.set_xlabel("log2(VAF_B / VAF_A)")
        ax.set_ylabel("Variant count")
        ax.set_title(
            f"{na[:45]}  vs  {nb[:45]}\n"
            f"n_shared = {r['n_shared']}",
            fontsize=8,
        )
        ax.legend(fontsize=7, loc="upper left")
        plt.tight_layout()

        label_a = re.sub(r"[^A-Za-z0-9_\-]", "_", na)[:50]
        label_b = re.sub(r"[^A-Za-z0-9_\-]", "_", nb)[:50]
        out = plot_dir / f"{label_a}__vs__{label_b}.png"
        plt.savefig(out, dpi=150)
        plt.close(fig)
        logging.info("  Plot: %s", out.name)


# -- Main --------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)s  %(message)s",
        datefmt="%H:%M:%S",
    )

    if not args.vcf_dir.is_dir():
        logging.error("VCF directory not found: %s", args.vcf_dir)
        sys.exit(1)

    vcf_files = sorted(args.vcf_dir.glob(args.vcf_glob))
    if not vcf_files:
        logging.error('No VCF files matching "%s" in %s', args.vcf_glob, args.vcf_dir)
        sys.exit(1)

    n = len(vcf_files)
    n_pairs = n * (n - 1) // 2
    logging.info("Found %d VCF files -> %d pairwise comparisons", n, n_pairs)

    # -- Setup ----------------------------------------------------------------
    args.outdir.mkdir(parents=True, exist_ok=True)
    filtered_dir = args.outdir / "filtered"
    filtered_dir.mkdir(exist_ok=True)

    # -- Phase 1: Pre-filter each VCF ----------------------------------------
    logging.info("Phase 1: Pre-filtering VCFs (%s, DP>=%d, AF>=%.3f)...",
                 "PASS only" if not args.include_non_pass else "all variants",
                 args.min_dp, args.min_af)

    sample_names: List[str] = []
    filtered_vcfs: Dict[str, Path] = {}

    for vcf in vcf_files:
        name = get_sample_name(vcf)
        sample_names.append(name)
        dst = filtered_dir / vcf.name
        tbi = Path(str(dst) + ".tbi")
        cfg = Path(str(dst) + ".filter.json")

        # Build the current filter config for cache-invalidation comparison
        current_cfg = {
            "pass_only":   not args.include_non_pass,
            "min_dp":      args.min_dp,
            "min_af":      args.min_af,
            "max_gnomad":  args.max_gnomad,
        }

        def _cfg_matches() -> bool:
            if not cfg.exists():
                return False
            import json as _json
            try:
                return _json.loads(cfg.read_text()) == current_cfg
            except Exception:
                return False

        if not args.force_refilter and dst.exists() and tbi.exists() and _cfg_matches():
            logging.info("  Reusing: %s", dst.name)
        else:
            if dst.exists() and not _cfg_matches():
                logging.info(
                    "  Filter args changed — re-filtering: %s", vcf.name
                )
            else:
                logging.info("  Filtering: %s", vcf.name)
            try:
                filter_vcf(vcf, dst, pass_only=not args.include_non_pass,
                           min_dp=args.min_dp, min_af=args.min_af,
                           max_gnomad=args.max_gnomad)
                import json as _json
                cfg.write_text(_json.dumps(current_cfg))
            except RuntimeError as e:
                logging.error("Filter failed: %s", e)
                sys.exit(1)

        filtered_vcfs[name] = dst

    # -- Phase 2: Load filtered VCFs into DataFrames --------------------------
    logging.info("Phase 2: Loading filtered VCFs...")
    df_cache: Dict[str, pd.DataFrame] = {}

    for name in sample_names:
        df = load_vcf(filtered_vcfs[name])
        df_cache[name] = df
        logging.info("  %-55s %4d variants", name[:55], len(df))

    # -- Phase 3: Pairwise comparisons ----------------------------------------
    logging.info(
        "Phase 3: %d pairwise comparisons (%d thread%s)...",
        n_pairs, args.threads, "s" if args.threads != 1 else "",
    )

    # Use ThreadPoolExecutor to avoid DataFrame serialisation overhead.
    # pandas/numpy release the GIL for array operations, so threads
    # provide real parallelism here without the pickle cost of multiprocessing.
    tasks = [
        (na, nb, df_cache, args.bin_width)
        for na, nb in combinations(sample_names, 2)
    ]

    if args.threads > 1:
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            raw_results = list(executor.map(_worker, tasks))
    else:
        raw_results = [_worker(t) for t in tasks]

    # -- Phase 4: Flag and write outputs --------------------------------------
    logging.info("Phase 4: Writing outputs...")

    # Load FREEMIX data if supplied
    freemix_data = None
    if args.freemix_file is not None:
        freemix_data = _load_freemix(args.freemix_file)

    summary_df = _apply_flags(
        raw_results,
        min_shared=args.min_shared,
        peak_count_thresh=args.peak_count,
        n_shared_z_thresh=args.n_shared_z,
        freemix=freemix_data,
        freemix_threshold=args.freemix_threshold,
    )
    write_summary(summary_df, args.outdir)
    # Determine matrix ordering
    if args.plate_layout:
        matrix_order = _get_plate_order(args.plate_layout, sample_names)
        logging.info("Matrix ordered by plate position from %s", args.plate_layout.name)
    else:
        matrix_order = sample_names
    write_matrix(raw_results, matrix_order, args.outdir)

    flagged_keys = set(
        zip(
            summary_df.loc[summary_df["flagged"], "sample_a"],
            summary_df.loc[summary_df["flagged"], "sample_b"],
        )
    )
    for r in raw_results:
        r["flagged"] = (r["sample_a"], r["sample_b"]) in flagged_keys

    flagged = [r for r in raw_results if r["flagged"]]

    if flagged:
        # Sort by peak_count descending - most suspicious first
        flagged.sort(key=lambda r: r["peak_count"], reverse=True)

        # Limit detail/plot output to top N
        if args.max_output > 0:
            output_pairs = flagged[:args.max_output]
        else:
            output_pairs = flagged

        logging.info(
            "%d pair(s) flagged; writing details for top %d...",
            len(flagged), len(output_pairs),
        )

        # Compute shared DataFrames once for output pairs (avoids re-running
        # _compare in both write_flagged_details and write_plots)
        shared_cache: Dict[tuple, pd.DataFrame] = {}
        for r in output_pairs:
            na, nb = r["sample_a"], r["sample_b"]
            _, shared = _compare(na, df_cache[na], nb, df_cache[nb],
                                 args.bin_width, include_shared_df=True)
            shared_cache[(na, nb)] = shared

        write_flagged_details(output_pairs, shared_cache, args.outdir)
        if args.plots:
            write_plots(output_pairs, shared_cache, args.outdir, args.bin_width)
    else:
        logging.info("No pairs flagged.")

    logging.info("Done. Results in %s/", args.outdir)


if __name__ == "__main__":
    main()
