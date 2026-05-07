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
import tempfile
from itertools import combinations
from multiprocessing import Pool, cpu_count
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

# CSQ field index for SYMBOL (0-based, from standard Uranus VEP format string)
_CSQ_SYMBOL_IDX = 1


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
        help="Minimum read depth (FORMAT/DP) to retain a variant "
             "(default: %(default)s)",
    )
    p.add_argument(
        "--min-shared", type=int, default=10,
        help="Minimum number of informative shared variants required before "
             "assessing a pair (default: %(default)s)",
    )
    p.add_argument(
        "--peak-count", type=int, default=8,
        help="Flag a pair if the non-unity peak contains >= this many "
             "variants (default: %(default)s)",
    )
    p.add_argument(
        "--peak-fraction", type=float, default=0.30,
        help="Flag a pair if the non-unity peak fraction >= this value "
             "(default: %(default)s)",
    )
    p.add_argument(
        "--threads", "-t", type=int, default=min(8, cpu_count()),
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
    r = _run(["bcftools", "view", "-h", str(vcf)])
    for line in r.stdout.decode().splitlines():
        if line.startswith("#CHROM"):
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 9:
                return parts[9]
    return vcf.stem


def filter_vcf(src: Path, dst: Path, pass_only: bool, min_dp: int,
               min_af: float) -> None:
    """
    Apply quality filters for contamination detection.

    Minimal pipeline -- only filters that protect data quality without removing
    valid contamination markers:

      1. bcftools view [-f PASS]
             Optionally restrict to PASS variants (Sentieon TNfilter).

      2. bcftools filter -e '(FORMAT/DP < min_dp || AF < min_af)'
             Remove low-depth and low-VAF variants where allele fractions are
             unreliable.

    No gnomAD, Prev_Count_AC, or synonymous filters are applied. These would
    remove genuine contamination markers: cross-sample contamination consists
    primarily of the source patient's germline heterozygous SNPs (which have
    population-level gnomAD AFs) including synonymous variants.
    """
    cmd1 = ["bcftools", "view", str(src)]
    if pass_only:
        cmd1 += ["-f", "PASS"]
    cmd1 += ["-Ou"]

    cmd2 = [
        "bcftools", "filter",
        "-e", f"(FORMAT/DP<{min_dp} || AF<{min_af})",
        "-Oz", "-o", str(dst),
    ]

    logging.debug("Filter pipeline for %s", src.name)
    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    p1.stdout.close()

    _, p2_err = p2.communicate()
    p1.wait()

    for proc, name in [(p1, "view"), (p2, "filter")]:
        if proc.returncode != 0:
            stderr = proc.stderr.read().decode() if proc.stderr else ""
            raise RuntimeError(
                f"bcftools {name} failed for {src.name}:\n{stderr}\n{p2_err.decode()}"
            )

    _run(["bcftools", "index", "-t", str(dst)])


# -- VCF loading --------------------------------------------------------------

def load_vcf(vcf: Path) -> pd.DataFrame:
    """
    Load a filtered VCF into a DataFrame.

    Extracts CHROM, POS, REF, ALT, FILTER, AF, and the gene SYMBOL from the
    first transcript in the CSQ annotation.

    Returns columns: chrom, pos, ref, alt, filter_col, af, gene
    """
    r = _run([
        "bcftools", "query",
        "-f", "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%AF]\t%INFO/CSQ\n",
        str(vcf),
    ])

    rows = []
    for line in r.stdout.decode().splitlines():
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        chrom, pos_s, ref, alt, filt, af_s = parts[:6]
        csq_s = parts[6] if len(parts) > 6 else ""

        try:
            af = float(af_s)
        except ValueError:
            continue

        # Parse gene symbol from first CSQ transcript entry
        gene = ""
        if csq_s and csq_s != ".":
            first = csq_s.split(",")[0]
            csq_parts = first.split("|")
            if len(csq_parts) > _CSQ_SYMBOL_IDX:
                gene = csq_parts[_CSQ_SYMBOL_IDX]

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
    min_af: float,
    bin_width: float,
    include_shared_df: bool = False,
) -> Tuple[dict, Optional[pd.DataFrame]]:
    """
    Core pairwise comparison with dual-peak detection.

    Merges two variant DataFrames on (chrom, pos, ref, alt), applies the VAF
    floor, computes log2(AF_B / AF_A), and finds both the overall peak and the
    non-unity peak (contamination signal).

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
    inf = merged[(merged["af_a"] >= min_af) & (merged["af_b"] >= min_af)].copy()
    n_inf = len(inf)

    nan = float("nan")
    result: dict = {
        "sample_a": name_a, "sample_b": name_b,
        "n_shared": n_shared, "n_informative": n_inf,
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

    if n_inf < 2:
        return result, None

    ratios = (inf["af_b"] / inf["af_a"]).clip(1e-6, 1e6)
    inf["ratio"] = ratios.values
    inf["log2_ratio"] = np.log2(ratios.values)

    peaks = _find_peaks(inf["log2_ratio"].values, bin_width)

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
        inf["in_peak"] = (
            (inf["log2_ratio"] >= plr - half) &
            (inf["log2_ratio"] <= plr + half)
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
        inf["in_peak"] = False

    return result, (inf if include_shared_df else None)


def _worker(task: tuple) -> dict:
    """Multiprocessing worker. Returns result dict only."""
    name_a, df_a, name_b, df_b, min_af, bin_width = task
    result, _ = _compare(name_a, df_a, name_b, df_b,
                         min_af, bin_width,
                         include_shared_df=False)
    return result


# -- Output writers -----------------------------------------------------------

def _apply_flags(
    results: List[dict],
    min_shared: int,
    peak_count_thresh: int,
    peak_frac_thresh: float,
) -> pd.DataFrame:
    """Convert results to DataFrame and flag pairs based on NON-UNITY peak."""
    df = pd.DataFrame(results)
    df["flagged"] = (
        (df["n_informative"] >= min_shared) &
        (
            (df["peak_count"] >= peak_count_thresh) |
            (df["peak_fraction"] >= peak_frac_thresh)
        )
    )
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
    df_cache: Dict[str, pd.DataFrame],
    outdir: Path,
    min_af: float,
    bin_width: float,
) -> None:
    """Write per-variant detail TSV for each flagged pair."""
    detail_dir = outdir / "flagged_pairs"
    detail_dir.mkdir(exist_ok=True)

    for r in flagged:
        na, nb = r["sample_a"], r["sample_b"]
        _, shared = _compare(
            na, df_cache[na], nb, df_cache[nb],
            min_af, bin_width,
            include_shared_df=True,
        )
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
    df_cache: Dict[str, pd.DataFrame],
    outdir: Path,
    min_af: float,
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
        _, shared = _compare(
            na, df_cache[na], nb, df_cache[nb],
            min_af, bin_width,
            include_shared_df=True,
        )
        if shared is None or shared.empty or "log2_ratio" not in shared.columns:
            continue

        fig, ax = plt.subplots(figsize=(10, 4))
        ax.hist(shared["log2_ratio"], bins=bins,
                color="steelblue", edgecolor="white", linewidth=0.4, alpha=0.85)

        # Reference ratio lines
        for lr, lbl in [(-1.0, "2:1"), (-2.0, "4:1"), (1.0, "1:2"), (2.0, "1:4")]:
            ax.axvline(lr, color="grey", linestyle=":", linewidth=0.8, alpha=0.6)
            ax.text(lr, ax.get_ylim()[1] * 0.92, lbl,
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
            f"n_informative = {r['n_informative']}",
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

        if not args.force_refilter and dst.exists() and tbi.exists():
            logging.info("  Reusing: %s", dst.name)
        else:
            logging.info("  Filtering: %s", vcf.name)
            try:
                filter_vcf(vcf, dst, pass_only=not args.include_non_pass,
                           min_dp=args.min_dp, min_af=args.min_af)
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

    tasks = [
        (na, df_cache[na], nb, df_cache[nb],
         args.min_af, args.bin_width)
        for na, nb in combinations(sample_names, 2)
    ]

    if args.threads > 1 and n_pairs > 50:
        with Pool(args.threads) as pool:
            raw_results = pool.map(_worker, tasks)
    else:
        raw_results = [_worker(t) for t in tasks]

    # -- Phase 4: Flag and write outputs --------------------------------------
    logging.info("Phase 4: Writing outputs...")

    summary_df = _apply_flags(
        raw_results,
        min_shared=args.min_shared,
        peak_count_thresh=args.peak_count,
        peak_frac_thresh=args.peak_fraction,
    )
    write_summary(summary_df, args.outdir)
    write_matrix(raw_results, sample_names, args.outdir)

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
        logging.info("Writing details for %d flagged pair(s)...", len(flagged))
        write_flagged_details(
            flagged, df_cache, args.outdir,
            args.min_af, args.bin_width,
        )
        if args.plots:
            write_plots(
                flagged, df_cache, args.outdir,
                args.min_af, args.bin_width,
            )
    else:
        logging.info("No pairs flagged.")

    logging.info("Done. Results in %s/", args.outdir)


if __name__ == "__main__":
    main()
