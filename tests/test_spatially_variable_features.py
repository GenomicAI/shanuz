"""Tests for find_spatially_variable_features (Moran's I)."""
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from shanuz import create_shanuz_object, find_spatially_variable_features
from shanuz.preprocessing import normalize_data
from shanuz.spatial.fov import create_fovs
from shanuz.spatial.variable_features import (
    _knn_weights,
    _morans_i,
    _morans_moments,
)


# ---------------------------------------------------------------------------
# Statistic correctness
# ---------------------------------------------------------------------------

def test_morans_i_matches_brute_force():
    """The vectorised (N/S0)·(zᵀWz)/(zᵀz) equals the explicit double sum."""
    rng = np.random.default_rng(0)
    n = 60
    xy = rng.uniform(0, 10, (n, 2))
    W = _knn_weights(xy, k=6)

    Z = rng.normal(size=(4, n))
    Z -= Z.mean(axis=1, keepdims=True)

    mine = _morans_i(Z, W)
    Wd = W.toarray()
    s0 = Wd.sum()
    brute = np.array([(n / s0) * (z @ Wd @ z) / (z @ z) for z in Z])
    np.testing.assert_allclose(mine, brute)


def test_morans_expected_value_is_minus_one_over_n_minus_1():
    rng = np.random.default_rng(1)
    xy = rng.uniform(0, 10, (50, 2))
    W = _knn_weights(xy, k=5)
    e_i, var = _morans_moments(W)
    assert np.isclose(e_i, -1.0 / (50 - 1))
    assert var > 0


def test_knn_weights_are_row_normalised():
    rng = np.random.default_rng(2)
    xy = rng.uniform(0, 10, (30, 2))
    W = _knn_weights(xy, k=5)
    rowsums = np.asarray(W.sum(axis=1)).ravel()
    np.testing.assert_allclose(rowsums, 1.0)
    assert W.diagonal().sum() == 0        # no self-weight


def test_knn_weights_reject_tiny_objects():
    with pytest.raises(ValueError, match="at least 3 cells"):
        _knn_weights(np.array([[0.0, 0.0], [1.0, 1.0]]), k=1)


# ---------------------------------------------------------------------------
# End-to-end
# ---------------------------------------------------------------------------

N_RAND = 22
RAND_GENES = [f"rand{i}" for i in range(N_RAND)]


@pytest.fixture
def spatial_obj():
    """120 cells; grad_x/grad_y/blob are spatially structured, rand* are not.

    The flat random genes deliberately dominate the library size. With only a
    couple of background genes, log-normalisation divides every gene by a
    spatially-structured total and leaks that structure into the random genes —
    a real compositional artefact, not a Moran's I error, but not what this test
    is measuring.
    """
    rng = np.random.default_rng(0)
    n = 120
    xy = rng.uniform(0, 10, (n, 2))
    genes = ["grad_x", "grad_y", "blob"] + RAND_GENES
    X = np.zeros((len(genes), n))
    X[0] = xy[:, 0]
    X[1] = xy[:, 1]
    X[2] = np.exp(-((xy[:, 0] - 5) ** 2 + (xy[:, 1] - 5) ** 2) / 4) * 10
    for i in range(N_RAND):
        X[3 + i] = rng.poisson(5, n)

    obj = create_shanuz_object(
        counts=sp.csc_matrix(np.abs(X) * 3),
        feature_names=genes,
        cell_names=[f"c{i}" for i in range(n)],
    )
    coords = pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1], "cell": obj.cell_names()})
    obj.images = create_fovs(coords, assay="RNA", default_name="rna")
    normalize_data(obj)
    return obj


def test_spatially_variable_ranks_structured_genes_first(spatial_obj):
    res = find_spatially_variable_features(spatial_obj, k=8)

    assert list(res.columns) == [
        "moransi", "moransi_pval", "moransi_padj", "moransi_rank",
    ]
    # The three spatially structured genes take the top three ranks.
    assert set(res.index[:3]) == {"grad_x", "grad_y", "blob"}
    assert (res.loc[["grad_x", "grad_y", "blob"], "moransi"] > 0.5).all()
    assert (res.loc[["grad_x", "grad_y", "blob"], "moransi_padj"] < 0.01).all()
    # Spatially random genes: I near zero and not significant.
    assert (res.loc[RAND_GENES, "moransi"].abs() < 0.2).all()
    assert (res.loc[RAND_GENES, "moransi_padj"] > 0.05).all()
    # Sorted by rank, 1 = most spatially variable.
    assert res["moransi_rank"].is_monotonic_increasing
    assert res["moransi_rank"].iloc[0] == 1


def test_results_written_to_feature_metadata(spatial_obj):
    res = find_spatially_variable_features(spatial_obj, k=8)
    meta = spatial_obj.assays["RNA"].meta_data
    for col in ("moransi", "moransi_pval", "moransi_padj", "moransi_rank"):
        assert col in meta.columns
    assert np.isclose(meta.loc["blob", "moransi"], res.loc["blob", "moransi"])


def test_features_argument_restricts_output(spatial_obj):
    res = find_spatially_variable_features(
        spatial_obj, features=["blob", "rand0"], k=8)
    assert set(res.index) == {"blob", "rand0"}
    assert res.index[0] == "blob"          # the structured one ranks first


def test_unknown_features_raise(spatial_obj):
    with pytest.raises(ValueError, match="None of the requested features"):
        find_spatially_variable_features(spatial_obj, features=["nope"], k=8)


def test_unsupported_method_raises(spatial_obj):
    with pytest.raises(NotImplementedError, match="markvariogram"):
        find_spatially_variable_features(spatial_obj, method="markvariogram")


def test_object_without_coordinates_raises():
    obj = create_shanuz_object(
        counts=sp.csc_matrix(np.ones((3, 5))),
        feature_names=["a", "b", "c"],
        cell_names=[f"c{i}" for i in range(5)],
    )
    with pytest.raises(ValueError, match="no spatial"):
        find_spatially_variable_features(obj)
