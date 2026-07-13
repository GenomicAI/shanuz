"""Tests for the bimod branch of find_markers (McDavid 2013 bimodal LRT)."""
import numpy as np
import pytest
import scipy.sparse as sp
from scipy.stats import norm

from shanuz import create_shanuz_object, find_markers
from shanuz.markers import _bimod_likelihood, _bimod_pvalue
from shanuz.preprocessing import normalize_data


def test_bimod_likelihood_matches_mcdavid_formula():
    """_bimod_likelihood equals the explicit point-mass + Gaussian log-likelihood."""
    x = np.array([0.0, 0.0, 1.0, 2.0, 3.0])
    x_pos = x[x > 0]
    alpha = len(x_pos) / len(x)  # 0.6
    expected = (
        (len(x) - len(x_pos)) * np.log(1 - alpha)
        + len(x_pos) * np.log(alpha)
        + norm.logpdf(x_pos, x_pos.mean(), np.std(x_pos, ddof=1)).sum()
    )
    assert np.isclose(_bimod_likelihood(x), expected)


def test_bimod_likelihood_edge_cases():
    # All zero -> only the detection term, with clamped alpha.
    assert np.isfinite(_bimod_likelihood(np.zeros(10)))
    # Single detected value (n_pos < 2) uses sd=1 rather than blowing up.
    assert np.isfinite(_bimod_likelihood(np.array([0.0, 0.0, 5.0])))
    # Identical detected values (sd == 0) is guarded, not NaN/inf.
    assert np.isfinite(_bimod_likelihood(np.array([0.0, 2.0, 2.0, 2.0])))


def test_bimod_pvalue_separates_distributions():
    rng = np.random.default_rng(1)
    same = _bimod_pvalue(rng.normal(2, 1, 60), rng.normal(2, 1, 60))
    diff = _bimod_pvalue(
        np.r_[np.zeros(30), rng.normal(3, 0.5, 20)], rng.normal(0.1, 0.3, 50)
    )
    assert 0.0 <= same <= 1.0 and 0.0 <= diff <= 1.0
    assert same > 0.1        # same distribution -> not significant
    assert diff < 1e-6       # strongly different -> highly significant


@pytest.fixture
def bimod_obj():
    """80 cells (A=40, B=40), 20 genes. g0 ↑A, g1 ↑B, g2 detection ↑A, g8 null."""
    rng = np.random.default_rng(1)
    n, G = 40, 20
    A = rng.poisson(1.0, size=(G, n)).astype(float)
    B = rng.poisson(1.0, size=(G, n)).astype(float)
    A[0] += 8
    B[1] += 8
    A[2], B[2] = rng.poisson(3.0, n), rng.poisson(0.1, n)  # detection ↑ in A
    obj = create_shanuz_object(
        counts=sp.csc_matrix(np.hstack([A, B])),
        feature_names=[f"g{i}" for i in range(G)],
        cell_names=[f"c{i}" for i in range(2 * n)],
    )
    normalize_data(obj)
    obj.idents = ["A"] * n + ["B"] * n
    return obj


def test_bimod_columns_and_marker_direction(bimod_obj):
    res = find_markers(bimod_obj, ident_1="A", test_use="bimod",
                       min_pct=0.0, logfc_threshold=0.0)

    assert list(res.columns) == ["p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]
    assert res.loc["g0", "p_val_adj"] < 0.05 and res.loc["g0", "avg_log2FC"] > 0
    assert res.loc["g1", "p_val_adj"] < 0.05 and res.loc["g1", "avg_log2FC"] < 0
    # bimod's Bernoulli component catches a detection-only difference.
    assert res.loc["g2", "pct.1"] > 0.8 and res.loc["g2", "pct.2"] < 0.4
    assert res.loc["g2", "p_val_adj"] < 0.05
    # Null gene stays non-significant after Bonferroni; output sorted by p_val.
    assert res.loc["g8", "p_val_adj"] > 0.05
    assert res.index[0] in {"g0", "g1", "g2"}
    assert res["p_val"].is_monotonic_increasing
