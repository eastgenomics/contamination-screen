"""Unit tests for _compare directionality logic."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np
from contamination_screen import _compare


def _make_df(positions, af_value):
    """Create a minimal DataFrame with n variants at a fixed AF."""
    return pd.DataFrame({
        "chrom": ["1"] * len(positions),
        "pos": positions,
        "ref": ["A"] * len(positions),
        "alt": ["T"] * len(positions),
        "filter_col": ["PASS"] * len(positions),
        "af": [af_value] * len(positions),
        "gene": ["TEST"] * len(positions),
    })


def test_a_contaminates_b():
    """AF_A high (~0.5), AF_B low (~0.12) => A is source, peak_log2 < 0."""
    positions = list(range(1, 21))
    df_a = _make_df(positions, af_value=0.50)
    df_b = _make_df(positions, af_value=0.12)
    result, _ = _compare("sampleA", df_a, "sampleB", df_b, bin_width=0.2)
    assert result["contamination_source"] == "sampleA", f"Expected sampleA as source, got {result['contamination_source']}"
    assert result["contamination_recipient"] == "sampleB"
    assert result["non_unity_log2"] < 0, f"Expected negative log2, got {result['peak_log2_ratio']}"
    assert result["non_unity_count"] == 20
    print("PASS: test_a_contaminates_b")


def test_b_contaminates_a():
    """AF_B high (~0.5), AF_A low (~0.12) => B is source, peak_log2 > 0."""
    positions = list(range(1, 21))
    df_a = _make_df(positions, af_value=0.12)
    df_b = _make_df(positions, af_value=0.50)
    result, _ = _compare("sampleA", df_a, "sampleB", df_b, bin_width=0.2)
    assert result["contamination_source"] == "sampleB", f"Expected sampleB as source, got {result['contamination_source']}"
    assert result["contamination_recipient"] == "sampleA"
    assert result["non_unity_log2"] > 0, f"Expected positive log2, got {result['peak_log2_ratio']}"
    assert result["non_unity_count"] == 20
    print("PASS: test_b_contaminates_a")


def test_no_contamination_similar_af():
    """Similar AF in both => peak at unity zone, no non-unity peak."""
    positions = list(range(1, 21))
    np.random.seed(42)
    df_a = _make_df(positions, af_value=0.50)
    # Add small noise so they're not identical
    df_b = df_a.copy()
    df_b["af"] = 0.50 + np.random.normal(0, 0.02, 20)
    df_b["af"] = df_b["af"].clip(0.03, 1.0)
    result, _ = _compare("sampleA", df_a, "sampleB", df_b, bin_width=0.2)
    # Non-unity peak should be weak (most variants in unity zone)
    assert result["non_unity_count"] <= 5, f"Expected weak non-unity peak, got {result['peak_count']}"
    print("PASS: test_no_contamination_similar_af")


def test_contamination_fraction():
    """50% contamination should report fraction ~0.25 (ratio = 0.25)."""
    positions = list(range(1, 21))
    df_a = _make_df(positions, af_value=0.50)
    df_b = _make_df(positions, af_value=0.125)  # 25% of source AF
    result, _ = _compare("sampleA", df_a, "sampleB", df_b, bin_width=0.2)
    assert result["contamination_source"] == "sampleA"
    # Fraction should be ~0.25 (2^-2 = 0.25)
    assert 0.2 < result["contamination_fraction"] < 0.3, \
        f"Expected fraction ~0.25, got {result['contamination_fraction']}"
    print("PASS: test_contamination_fraction")


if __name__ == "__main__":
    test_a_contaminates_b()
    test_b_contaminates_a()
    test_no_contamination_similar_af()
    test_contamination_fraction()
    print("\nAll tests passed.")
