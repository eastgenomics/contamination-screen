#!/bin/bash
# Phase 5: Run contamination screen across all 100 runs.
# If cached filtered VCFs exist in results/RUN*/filtered/, they are used
# directly — no download required. Download only happens on first pass.
# If select_runs.py reported any archived files, run unarchive_all.py first.
set -euo pipefail

MAPPING=~/Downloads/contamination_screen/contamination_screen_mapping.tsv
RESULTS="retrospective/results"
SCREEN="contamination_screen.py"
DXPY=~/Documents/dx-grab/.venv/bin/python3

mkdir -p "$RESULTS"

for RUN in $(seq -w 001 100); do
    RUN_ID="RUN${RUN}"
    RUNDIR="$RESULTS/$RUN_ID"

    if [ -f "$RUNDIR/summary.tsv" ]; then
        echo "=== $RUN_ID already complete, skipping ==="
        continue
    fi

    echo "=== Processing $RUN_ID ==="

    # If filtered VCFs already exist, pass that directory directly —
    # contamination_screen.py will find the cached files and skip re-filtering.
    # This avoids re-downloading for re-screening runs.
    FILTERED_DIR="$RUNDIR/filtered"
    if [ -d "$FILTERED_DIR" ]; then
        N_CACHED=$(find "$FILTERED_DIR" -maxdepth 1 -type f -name '*.vcf.gz' | wc -l)
    else
        N_CACHED=0
    fi

    if [ "$N_CACHED" -gt 0 ]; then
        echo "  Using $N_CACHED cached filtered VCFs (no download needed)"
        VCF_SOURCE="$FILTERED_DIR"
    else
        TMPDIR=$(mktemp -d)
        trap "rm -rf '$TMPDIR'" EXIT

        "$DXPY" scripts/download_run.py \
            --run "$RUN_ID" \
            --mapping "$MAPPING" \
            --outdir "$TMPDIR" \
            --skip-existing

        N_CACHED=$(find "$TMPDIR" -maxdepth 1 -type f -name '*.vcf.gz' | wc -l)
        if [ "$N_CACHED" -eq 0 ]; then
            echo "  No VCFs for $RUN_ID (all archived?), skipping"
            rm -rf "$TMPDIR"; trap - EXIT; continue
        fi
        VCF_SOURCE="$TMPDIR"
    fi

    echo "  $N_CACHED samples"

    python3 "$SCREEN" "$VCF_SOURCE" \
        --outdir "$RUNDIR" \
        --vcf-glob "${RUN_ID}_S*_tnhaplotyper2_annotated.vcf.gz" \
        --threads 8 \
        --plots

    # Clean up temp download dir if used
    if [ -n "${TMPDIR:-}" ]; then
        rm -rf "$TMPDIR"
        trap - EXIT
    fi

    echo "=== $RUN_ID complete ==="
done
