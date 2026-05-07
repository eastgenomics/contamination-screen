#!/usr/bin/env python3
"""
contamination_screen.py

Pairwise cross-sample contamination detector for tumour sequencing panel VCFs.

Design rationale
----------------
Cross-sample contamination from an unrelated tumour leaves a characteristic
fingerprint in variant call data: rare somatic variants from the source sample
appear in the recipient at a consistent fraction of the source VAF. That fraction
is the contamination level (e.g. 50% contamination → shared variants appear at
~half the source VAF in the recipient).

This tool detects that signature by:

  1. Pre-filtering each VCF to remove high-frequency population variants and
     most synonymous variants (retaining GATA2 / TP53 synonymous as clinically
     relevant markers). This step uses a 4-stage bcftools pipeline to work
     around the fact that the gnomAD and Prev_Count INFO tags are typed as
     String in the VCF header: the conflicting tags are stripped, CSQ subfields
     are expanded into properly typed INFO tags via split-vep, and the resulting
     VCF is hard-filtered.

  2. Performing every pairwise comparison in the cohort (N*(N-1)/2 pairs).

  3. For each pair, computing log2(VAF_B / VAF_A) for all shared variants with
     VAF >= min_af in both samples, then scanning the log2-ratio histogram for
     a cluster of variants at a consistent non-unity ratio — the hallmark of
     contamination.

  4. Reporting the implied contamination direction and fraction for each flagged
     pair. A log2-ratio peak at a negative value means sample_a is the source
     (its variants appear at lower VAF in sample_b); a positive peak means
     sample_b is the source.

Output files
------------
  results/
  ├── filtered/            Pre-filtered VCFs (reused on re-runs)
  ├── summary.tsv          One row per pair; flagged pairs marked
  ├── matrix.tsv           N×N directional matrix [recipient, source] = peak_count
  ├── flagged_pairs/       Detailed variant table per flagged pair
  └── plots/               Log2-ratio histograms (--plots only)
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


# ── Constants ──────────────────────────────────────────────────────────────────

# The filter pipeline mimics the clinical pipeline exactly:
#
#   bcftools norm -m -any
#     Split multiallelic records and normalise indels.
#
#   bcftools +split-vep -c - -p CSQ_ -s worst
#     Extract ALL CSQ subfields for the worst-consequence transcript into
#     CSQ_-prefixed INFO tags (e.g. CSQ_SYMBOL, CSQ_gnomADg_AF, ...).
#     Using the CSQ_ prefix avoids name conflicts with the existing String-typed
#     INFO tags (gnomADg_AF, Prev_Count_AC, etc.) that were added by the upstream
#     annotation pipeline — those remain untouched as their original names.
#
#   bcftools annotate  (reheader CSQ_Prev_Count_AC String → Integer)
#     split-vep assigns String to Prev_Count_AC because the field name does not
#     match the built-in .*_AF Float pattern. A one-line header fix recasts it
#     so that the arithmetic comparison >853 works in the filter expression.
#
#   bcftools filter --soft-filter EXCLUDE -m +
#     The exact expression from the clinical pipeline. Soft-tags matching
#     records with FILTER=EXCLUDE (appended to existing FILTER value).
#
#   bcftools view [-f PASS] -e 'FILTER~"EXCLUDE"'
#     Hard-filter: drop EXCLUDE-tagged records. With --include-non-pass omitted
#     (default), also restrict to records that were originally PASS.

# All CSQ fields extracted by split-vep (- = all fields)
_SPLIT_VEP_COLS = "-"

# Prefix for split-vep INFO tags → CSQ_SYMBOL, CSQ_gnomADg_AF, …
_VEP_PREFIX = "CSQ_"

# Header line to recast CSQ_Prev_Count_AC from String to Integer
_PREVCOUNT_HDR = (
    '##INFO=<ID=CSQ_Prev_Count_AC,Number=1,Type=Integer,'
    'Description="Cohort prevalence AC (from CSQ)">'
)

# Soft-filter name written to the FILTER column by bcftools filter
_SOFT_FILTER_NAME = "EXCLUDE"

# Population / synonymous exclude expression — identical to the clinical pipeline.
# Applied with bcftools filter --soft-filter EXCLUDE -m + -e '...'
# Excludes:
#   - Cohort Prev_Count_AC > 853  (recurrent artefact / common CH)
#   - gnomAD genome AF >= 0.002
#   - gnomAD exome AF  >= 0.002
#   - Synonymous variants, UNLESS gene is GATA2 or TP53
_EXCLUDE_EXPR = (
    '(CSQ_Prev_Count_AC>853 || CSQ_gnomADg_AF >= 0.002 || CSQ_gnomADe_AF >= 0.002)'
    ' || '
    '(CSQ_SYMBOL!="GATA2" || CSQ_Consequence!="synonymous_variant")'
    ' && '
    '(CSQ_SYMBOL!="TP53" || CSQ_Consequence!="synonymous_variant")'
    ' && '
    '(CSQ_Consequence=="synonymous_variant")'
)

# Histogram parameters
_LOG2_RANGE   = 4.0   # ±4 log2 units spans ratios from 1:16 to 16:1
_BIN_WIDTH    = 0.2   # bin width in log2 units
_CENTRAL_EXCL = 0.3   # exclude bins within ±0.3 log2 of 0 (ratio 0.81–1.23)
                      # to avoid flagging variants shared at similar VAF
                      # (e.g. germline variants in related samples)

VCF_GLOB = "*_annotated.vcf.gz"


# ── Argument parsing ────────────────────────────────────────────────────────────

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
        "--min-shared", type=int, default=10,
        help="Minimum number of informative shared variants required before "
             "assessing a pair (default: %(default)s)",
    )
    p.add_argument(
        "--peak-count", type=int, default=8,
        help="Flag a pair if the dominant ratio bin contains >= this many "
             "variants (default: %(default)s)",
    )
    p.add_argument(
        "--peak-fraction", type=float, default=0.30,
        help="Flag a pair if the dominant ratio bin fraction >= this value "
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
        "--central-excl", type=float, default=_CENTRAL_EXCL,
        help="Exclude log2 ratios within ±this value of 0 from peak detection "
             "(default: %(default)s)",
    )
    p.add_argument(
        "--force-refilter", action="store_true",
        help="Re-run bcftools filtering even if filtered VCFs already exist",
    )
    p.add_argument(
        "--verbose", "-v", action="store_true",
        help="Enable debug logging",
    )
    return p.parse_args()


# ── bcftools helpers ────────────────────────────────────────────────────────────

def _run(cmd: List[str], stdin=None, check: bool = True) -> subprocess.CompletedProcess:
    """Run a subprocess, logging the command at DEBUG level."""
    logging.debug("RUN: %s", " ".join(str(c) for c in cmd))
    return subprocess.run(
        cmd, stdin=stdin,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=check,
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


def get_csq_format(vcf: Path) -> List[str]:
    """Parse the CSQ Format field order from the VCF INFO header."""
    r = _run(["bcftools", "view", "-h", str(vcf)])
    for line in r.stdout.decode().splitlines():
        if "##INFO=<ID=CSQ" in line:
            m = re.search(r"Format: ([^\"]+)\"", line)
            if m:
                return m.group(1).strip().split("|")
    return []


def filter_vcf(src: Path, dst: Path, pass_only: bool, tmpdir: Path) -> None:
    """
    Apply the population-frequency and synonymous-variant filter to src,
    writing a bgzipped, tabix-indexed VCF to dst.

    Pipeline (all steps piped, no intermediate files):

      Note: bcftools norm -m -any has already been applied in the upstream
      clinical pipeline (confirmed via bcftools_normCommand in the VCF header).
      Multiallelic records are therefore already split; no norm step is needed.

      1. bcftools +split-vep -c - -p CSQ_ -s worst
             Extract ALL CSQ subfields for the worst-consequence transcript into
             CSQ_-prefixed INFO tags (CSQ_SYMBOL, CSQ_Consequence,
             CSQ_gnomADe_AF, CSQ_gnomADg_AF, CSQ_Prev_Count_AC, …).
             The CSQ_ prefix avoids name conflicts with the existing String-typed
             INFO tags (gnomADg_AF, Prev_Count_AC, etc.) added by the upstream
             annotation pipeline — those remain untouched under their original
             names. The .*_AF pattern in split-vep's built-in type rules
             automatically assigns Float to CSQ_gnomADg_AF and CSQ_gnomADe_AF.

      2. bcftools annotate  (reheader CSQ_Prev_Count_AC String → Integer)
             CSQ_Prev_Count_AC is assigned String by split-vep (no built-in
             type rule matches the field name). A one-line header override
             recasts it to Integer so the >853 arithmetic comparison works.

      3. bcftools filter --soft-filter EXCLUDE -m +
             Exact clinical pipeline expression. Records matching the exclude
             criteria have EXCLUDE appended to their FILTER column.

      4. bcftools filter -e '(FORMAT/DP<99 || AF<0.03)'
             Hard-filter: remove low-depth (DP < 99) and low-VAF (AF < 0.03)
             variants. Matches the clinical pipeline quality threshold.
             Note: AF < 0.03 is consistent with --min-af (default 0.03).

      5. bcftools view [-f PASS] -e 'FILTER~"EXCLUDE"'
             Hard-filter: drop all EXCLUDE-tagged records. In PASS-only mode
             (-f PASS), also drop any records that were not originally PASS.
    """
    hdr_file = tmpdir / "prevcount.hdr"
    hdr_file.write_text(_PREVCOUNT_HDR + "\n")

    # ── Step 1: split-vep ────────────────────────────────────────────────────
    cmd1 = [
        "bcftools", "+split-vep", str(src),
        "-c", _SPLIT_VEP_COLS,
        "-s", "worst",
        "-p", _VEP_PREFIX,
        "-Ou",
    ]

    # ── Step 2: reheader CSQ_Prev_Count_AC String → Integer ──────────────────
    cmd2 = [
        "bcftools", "annotate",
        "-x", f"INFO/{_VEP_PREFIX}Prev_Count_AC",
        "-h", str(hdr_file),
        "-Ou",
    ]

    # ── Step 3: soft-filter (exact clinical expression) ──────────────────────
    cmd3 = [
        "bcftools", "filter",
        "--soft-filter", _SOFT_FILTER_NAME,
        "-m", "+",
        "-e", _EXCLUDE_EXPR,
        "-Ou",
    ]

    # ── Step 4: hard-filter low depth / low VAF ───────────────────────────────
    cmd4 = [
        "bcftools", "filter",
        "-e", "(FORMAT/DP<99 || AF<0.03)",
        "-Ou",
    ]

    # ── Step 5: hard-filter EXCLUDE tags (+ optional PASS restriction) ────────
    cmd5 = ["bcftools", "view", "-e", f'FILTER~"{_SOFT_FILTER_NAME}"']
    if pass_only:
        cmd5 += ["-f", "PASS"]
    cmd5 += ["-Oz", "-o", str(dst)]

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
    p5 = subprocess.Popen(cmd5, stdin=p4.stdout, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    p4.stdout.close()

    _, p5_err = p5.communicate()
    p1.wait(); p2.wait(); p3.wait(); p4.wait()

    for proc, name in [
        (p1, "split-vep"), (p2, "reheader"),
        (p3, "filter/EXCLUDE"), (p4, "filter/DP+AF"), (p5, "view/hard-filter"),
    ]:
        if proc.returncode != 0:
            stderr = proc.stderr.read().decode() if proc.stderr else ""
            raise RuntimeError(
                f"bcftools {name} failed for {src.name}:\n{stderr}\n{p5_err.decode()}"
            )

    _run(["bcftools", "index", "-t", str(dst)])


# ── VCF loading ─────────────────────────────────────────────────────────────────

def load_vcf(vcf: Path) -> pd.DataFrame:
    """
    Load a filtered VCF into a DataFrame.

    The filtered VCF has already had CSQ subfields extracted by split-vep,
    so CSQ_SYMBOL is available directly as an INFO tag — no raw CSQ parsing
    needed.

    Returns columns: chrom, pos, ref, alt, filter_col, af, gene
    """
    r = _run([
        "bcftools", "query",
        "-f", "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%AF]\t%INFO/CSQ_SYMBOL\n",
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
    return df


# ── Core comparison ─────────────────────────────────────────────────────────────

def _find_peak(
    log2_ratios: np.ndarray,
    bin_width: float,
    central_excl: float,
) -> dict:
    """
    Scan a log2-ratio distribution for the tallest bin outside the central
    exclusion zone.

    The central exclusion zone (|log2_ratio| < central_excl) covers variants
    at similar VAF in both samples, which could represent biologically shared
    mutations (e.g. same patient, different timepoints) rather than contamination.

    Returns:
        log2_ratio  : centre of the peak bin (nan if no peak found)
        count       : number of variants in that bin
        fraction    : count / total variants
    """
    if len(log2_ratios) == 0:
        return {"log2_ratio": float("nan"), "count": 0, "fraction": 0.0}

    edges = np.arange(
        -_LOG2_RANGE - bin_width / 2,
        _LOG2_RANGE + bin_width,
        bin_width,
    )
    counts, edges = np.histogram(log2_ratios, bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2.0

    # Zero out the central bins to avoid flagging shared biology
    masked = counts.astype(float).copy()
    masked[np.abs(centers) < central_excl] = 0.0

    if masked.max() == 0:
        return {"log2_ratio": float("nan"), "count": 0, "fraction": 0.0}

    idx = int(np.argmax(masked))
    return {
        "log2_ratio": float(centers[idx]),
        "count": int(counts[idx]),
        "fraction": float(counts[idx]) / len(log2_ratios),
    }


def _compare(
    name_a: str, df_a: pd.DataFrame,
    name_b: str, df_b: pd.DataFrame,
    min_af: float,
    bin_width: float,
    central_excl: float,
    include_shared_df: bool = False,
) -> Tuple[dict, Optional[pd.DataFrame]]:
    """
    Core pairwise comparison.

    Merges two variant DataFrames on (chrom, pos, ref, alt), applies the VAF
    floor, computes log2(AF_B / AF_A) for each shared variant, and identifies
    the dominant ratio cluster.

    Directionality:
        peak_log2_ratio < 0  →  AF_A > AF_B  →  A is the contamination source
        peak_log2_ratio > 0  →  AF_B > AF_A  →  B is the contamination source

    The contamination fraction is the implied fraction of the recipient's
    library that came from the source: 2^peak_log2_ratio (if source is A) or
    2^(-peak_log2_ratio) (if source is B).

    Returns:
        result dict  : statistics for this pair
        shared_df    : merged DataFrame with ratio columns (None if not requested
                       or if n_informative < 2)
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
        "peak_log2_ratio": nan, "peak_ratio": nan,
        "peak_count": 0, "peak_fraction": 0.0,
        "contamination_source": "", "contamination_recipient": "",
        "contamination_fraction": nan,
        "flagged": False,
    }

    if n_inf < 2:
        return result, None

    ratios = (inf["af_b"] / inf["af_a"]).clip(1e-6, 1e6)
    inf["ratio"] = ratios.values
    inf["log2_ratio"] = np.log2(ratios.values)

    peak = _find_peak(inf["log2_ratio"].values, bin_width, central_excl)

    result.update({
        "peak_log2_ratio": peak["log2_ratio"],
        "peak_ratio": (
            2.0 ** peak["log2_ratio"]
            if not math.isnan(peak["log2_ratio"]) else nan
        ),
        "peak_count": peak["count"],
        "peak_fraction": peak["fraction"],
    })

    if not math.isnan(peak["log2_ratio"]) and peak["count"] > 0:
        half = bin_width / 2.0
        inf["in_peak"] = (
            (inf["log2_ratio"] >= peak["log2_ratio"] - half) &
            (inf["log2_ratio"] <= peak["log2_ratio"] + half)
        )
        if peak["log2_ratio"] < 0:
            # AF_A > AF_B: A's variants appear diluted in B → A is source
            result.update({
                "contamination_source": name_a,
                "contamination_recipient": name_b,
                "contamination_fraction": 2.0 ** peak["log2_ratio"],
            })
        else:
            # AF_B > AF_A: B's variants appear diluted in A → B is source
            result.update({
                "contamination_source": name_b,
                "contamination_recipient": name_a,
                "contamination_fraction": 2.0 ** (-peak["log2_ratio"]),
            })
    else:
        inf["in_peak"] = False

    return result, (inf if include_shared_df else None)


def _worker(task: tuple) -> dict:
    """
    Multiprocessing worker. Unpacks task tuple and returns the result dict only
    (no DataFrame). DataFrames for flagged pairs are recomputed in the main
    process to avoid pickle overhead across all 1128 workers.
    """
    name_a, df_a, name_b, df_b, min_af, bin_width, central_excl = task
    result, _ = _compare(name_a, df_a, name_b, df_b,
                         min_af, bin_width, central_excl,
                         include_shared_df=False)
    return result


# ── Output writers ──────────────────────────────────────────────────────────────

def _apply_flags(
    results: List[dict],
    min_shared: int,
    peak_count_thresh: int,
    peak_frac_thresh: float,
) -> pd.DataFrame:
    """Convert results list to DataFrame and apply flagging logic."""
    df = pd.DataFrame(results)
    df["flagged"] = (
        (df["n_informative"] >= min_shared) &
        (
            (df["peak_count"] >= peak_count_thresh) |
            (df["peak_fraction"] >= peak_frac_thresh)
        )
    )
    for col in ["peak_log2_ratio", "peak_ratio",
                "peak_fraction", "contamination_fraction"]:
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
    """
    Write a directional N×N contamination matrix.

    Rows = contamination recipient, columns = contamination source.
    Cell value = peak_count for the (source→recipient) direction.
    A clean cohort will have a near-zero matrix; contaminated pairs will show
    non-zero cells in the expected source/recipient positions.
    """
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
    central_excl: float,
) -> None:
    """
    For each flagged pair, write a TSV of all shared informative variants with
    their AF_A, AF_B, ratio, log2_ratio, in_peak, gene, and FILTER status.
    Variants are sorted by log2_ratio so the peak cluster is visible.
    """
    detail_dir = outdir / "flagged_pairs"
    detail_dir.mkdir(exist_ok=True)

    for r in flagged:
        na, nb = r["sample_a"], r["sample_b"]
        _, shared = _compare(
            na, df_cache[na], nb, df_cache[nb],
            min_af, bin_width, central_excl,
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
    central_excl: float,
) -> None:
    """Generate a log2-ratio histogram for each flagged pair."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logging.warning("matplotlib not available; skipping plots (pip install matplotlib)")
        return

    plot_dir = outdir / "plots"
    plot_dir.mkdir(exist_ok=True)

    bins = np.arange(-_LOG2_RANGE, _LOG2_RANGE + bin_width, bin_width)

    for r in flagged:
        na, nb = r["sample_a"], r["sample_b"]
        _, shared = _compare(
            na, df_cache[na], nb, df_cache[nb],
            min_af, bin_width, central_excl,
            include_shared_df=True,
        )
        if shared is None or shared.empty or "log2_ratio" not in shared.columns:
            continue

        fig, ax = plt.subplots(figsize=(10, 4))
        ax.hist(shared["log2_ratio"], bins=bins,
                color="steelblue", edgecolor="white", linewidth=0.4, alpha=0.85)

        # Central exclusion zone
        ax.axvspan(-central_excl, central_excl, alpha=0.08, color="grey",
                   label=f"Central exclusion (±{central_excl} log₂)")

        # Reference ratio lines
        for lr, lbl in [(-1.0, "2:1"), (-2.0, "4:1"), (1.0, "1:2"), (2.0, "1:4")]:
            ax.axvline(lr, color="grey", linestyle=":", linewidth=0.8, alpha=0.6)
            ax.text(lr, ax.get_ylim()[1] * 0.92, lbl,
                    ha="center", fontsize=7, color="grey")

        # Peak line
        plr = r.get("peak_log2_ratio", float("nan"))
        if not math.isnan(plr):
            src = r.get("contamination_source", "?")[:35]
            rec = r.get("contamination_recipient", "?")[:35]
            cf = r.get("contamination_fraction", float("nan"))
            ax.axvline(plr, color="crimson", linestyle="--", linewidth=1.8,
                       label=(
                           f"Peak log₂(r) = {plr:.2f}  (ratio ≈ {2**plr:.3f})\n"
                           f"Source → {src}\n"
                           f"n = {r['peak_count']},  "
                           f"fraction = {r['peak_fraction']:.2f},  "
                           f"contam ≈ {cf:.1%}"
                       ))

        ax.set_xlabel("log₂(VAF_B / VAF_A)")
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


# ── Main ────────────────────────────────────────────────────────────────────────

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
    logging.info(
        "Found %d VCF files → %d pairwise comparisons", n, n_pairs
    )

    # ── Setup ────────────────────────────────────────────────────────────────────
    args.outdir.mkdir(parents=True, exist_ok=True)
    filtered_dir = args.outdir / "filtered"
    filtered_dir.mkdir(exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir_s:
        tmpdir = Path(tmpdir_s)

        # ── Phase 1: Pre-filter each VCF ─────────────────────────────────────────
        logging.info("Phase 1: Pre-filtering VCFs (%s)...",
                     "PASS only" if not args.include_non_pass else "all variants")

        sample_names: List[str] = []
        filtered_vcfs: Dict[str, Path] = {}
        csq_format: List[str] = []

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
                               tmpdir=tmpdir)
                except RuntimeError as e:
                    logging.error("Filter failed: %s", e)
                    sys.exit(1)

            filtered_vcfs[name] = dst

            if not csq_format:
                csq_format = get_csq_format(vcf)
                if csq_format:
                    logging.debug("CSQ format has %d fields", len(csq_format))

        # ── Phase 2: Load filtered VCFs into DataFrames ───────────────────────────
        logging.info("Phase 2: Loading filtered VCFs...")
        df_cache: Dict[str, pd.DataFrame] = {}

        for name in sample_names:
            df = load_vcf(filtered_vcfs[name])
            df_cache[name] = df
            logging.info("  %-55s %4d variants", name[:55], len(df))

        # ── Phase 3: Pairwise comparisons ─────────────────────────────────────────
        logging.info(
            "Phase 3: %d pairwise comparisons (%d thread%s)...",
            n_pairs, args.threads, "s" if args.threads != 1 else "",
        )

        tasks = [
            (na, df_cache[na], nb, df_cache[nb],
             args.min_af, args.bin_width, args.central_excl)
            for na, nb in combinations(sample_names, 2)
        ]

        if args.threads > 1 and n_pairs > 50:
            with Pool(args.threads) as pool:
                raw_results = pool.map(_worker, tasks)
        else:
            # Sequential for small cohorts — avoids pickle overhead
            raw_results = [_worker(t) for t in tasks]

        # ── Phase 4: Flag and write outputs ───────────────────────────────────────
        logging.info("Phase 4: Writing outputs...")

        summary_df = _apply_flags(
            raw_results,
            min_shared=args.min_shared,
            peak_count_thresh=args.peak_count,
            peak_frac_thresh=args.peak_fraction,
        )
        write_summary(summary_df, args.outdir)
        write_matrix(raw_results, sample_names, args.outdir)

        # Re-sync flagged status into raw_results for detail/plot writers
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
                args.min_af, args.bin_width, args.central_excl,
            )
            if args.plots:
                write_plots(
                    flagged, df_cache, args.outdir,
                    args.min_af, args.bin_width, args.central_excl,
                )
        else:
            logging.info("No pairs flagged.")

    logging.info("Done. Results in %s/", args.outdir)


if __name__ == "__main__":
    main()
