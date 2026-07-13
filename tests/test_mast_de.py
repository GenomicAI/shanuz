"""Tests for the MAST two-part hurdle branch of find_markers (test_use='mast')."""
import numpy as np
import pytest
import scipy.sparse as sp

from shanuz import create_shanuz_object, find_markers
from shanuz.preprocessing import normalize_data


@pytest.fixture
def hurdle_obj():
    """80 cells (A=40, B=40), 20 genes with planted signals of different kinds.

    g0  magnitude up in A          g2  detection up in A (present A / mostly-absent B)
    g1  magnitude up in B          g3  up in A among always-detected genes (continuous only)
    g8  null
    """
    rng = np.random.default_rng(0)
    n, G = 40, 20
    A = rng.poisson(1.0, size=(G, n)).astype(float)
    B = rng.poisson(1.0, size=(G, n)).astype(float)
    A[0] += 8                                    # g0 magnitude ↑ in A
    B[1] += 8                                    # g1 magnitude ↑ in B
    A[2], B[2] = rng.poisson(3.0, n), rng.poisson(0.1, n)      # g2 detection ↑ in A (non-perfect)
    A[3], B[3] = rng.poisson(5.0, n) + 5, rng.poisson(5.0, n)  # g3 always-detected, ↑ in A
    obj = create_shanuz_object(
        counts=sp.csc_matrix(np.hstack([A, B])),
        feature_names=[f"g{i}" for i in range(G)],
        cell_names=[f"c{i}" for i in range(2 * n)],
    )
    normalize_data(obj)
    obj.idents = ["A"] * n + ["B"] * n
    return obj


def test_mast_columns_and_marker_direction(hurdle_obj):
    res = find_markers(hurdle_obj, ident_1="A", test_use="mast",
                       min_pct=0.0, logfc_threshold=0.0)

    assert list(res.columns) == ["p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]
    # g0 up in A, g1 up in B (down in A) — both significant with the right sign.
    assert res.loc["g0", "p_val_adj"] < 0.05 and res.loc["g0", "avg_log2FC"] > 0
    assert res.loc["g1", "p_val_adj"] < 0.05 and res.loc["g1", "avg_log2FC"] < 0
    # A null gene survives Bonferroni as non-significant.
    assert res.loc["g8", "p_val_adj"] > 0.05
    # A planted marker ranks first; output sorted by p_val.
    assert res.index[0] in {"g0", "g1", "g2", "g3"}
    assert res["p_val"].is_monotonic_increasing


def test_mast_detects_detection_only_difference(hurdle_obj):
    """The hurdle's discrete component flags a gene that differs mainly in
    *detection rate* (present in A, mostly absent in B)."""
    res = find_markers(hurdle_obj, ident_1="A", test_use="mast",
                       min_pct=0.0, logfc_threshold=0.0)
    assert res.loc["g2", "pct.1"] > 0.8 and res.loc["g2", "pct.2"] < 0.3
    assert res.loc["g2", "p_val_adj"] < 0.05


def test_mast_detects_continuous_only_difference(hurdle_obj):
    """A gene detected in ~all cells but higher in A is caught by the continuous
    (Gaussian) component when detection carries no signal."""
    res = find_markers(hurdle_obj, ident_1="A", test_use="mast",
                       min_pct=0.0, logfc_threshold=0.0)
    assert res.loc["g3", "pct.1"] > 0.95 and res.loc["g3", "pct.2"] > 0.95
    assert res.loc["g3", "p_val_adj"] < 0.05 and res.loc["g3", "avg_log2FC"] > 0


def test_mast_accepts_latent_vars(hurdle_obj):
    """Cellular detection rate as a latent covariate (Seurat's CDR default)."""
    counts = hurdle_obj.assays["RNA"].layers["counts"]
    detected = (counts > 0)
    cdr = np.asarray(detected.mean(axis=0)).ravel()
    hurdle_obj.meta_data["cdr"] = cdr
    res = find_markers(hurdle_obj, ident_1="A", test_use="mast", latent_vars=["cdr"],
                       min_pct=0.0, logfc_threshold=0.0)
    assert res.loc["g0", "p_val_adj"] < 0.05


def test_mast_unknown_test_still_errors(hurdle_obj):
    with pytest.raises(ValueError, match="mast"):
        find_markers(hurdle_obj, ident_1="A", test_use="not_a_test")
