"""Regression guards for the spatial defects T-sp found against Seurat 5.5.1.

Each test here pins a number or a shape that shanuz got wrong and that the suite
did not notice: the whole 659-test run passed both before and after all three
fixes, because nothing asserted a spot radius, a ring's vertex count, or which
weight matrix Moran's I was built on.

Where a constant appears below with an "R:" comment, it was read off a live
Seurat 5.5.1 / SeuratObject 5.4.0 session on the same input, not derived from
shanuz's own output.
"""
import numpy as np
import pandas as pd
import pytest

from shanuz.plotting import _boundary_radius
from shanuz.spatial import create_centroids, create_fov, create_segmentation
from shanuz.spatial.variable_features import (
    _inverse_square_lag,
    _morans_i_from_lag,
    _moransi_table,
)


# ---------------------------------------------------------------------------
# Moran's I — the weight matrix
# ---------------------------------------------------------------------------

def _naive_morans_i(Z, xy):
    """Seurat's Moran's I written the obvious way: dense, unblocked, slow.

    RunMoransI builds `1 / pos.dist.mat ^ 2` whole, and Rfast2::moranI then
    row-standardises it and returns `cor(y, x) * sd(y) / sd(x)`. This is that,
    transcribed — the reference the blocked implementation has to reproduce.
    """
    from scipy.spatial.distance import squareform, pdist

    d = squareform(pdist(xy))
    with np.errstate(divide="ignore"):
        W = 1.0 / d**2
    np.fill_diagonal(W, 0.0)
    W = W / W.sum(axis=1, keepdims=True)
    out = []
    for z in Z:
        y = W @ z
        zc, yc = z - z.mean(), y - y.mean()
        out.append(float((yc @ zc) / (zc @ zc)))
    return np.array(out)


@pytest.fixture
def toy_field():
    """12x12 grid of cells carrying one smooth gradient and one noise gene."""
    rng = np.random.default_rng(7)
    gx, gy = np.meshgrid(np.arange(12.0), np.arange(12.0))
    xy = np.column_stack([gx.ravel(), gy.ravel()])
    smooth = xy[:, 0] + xy[:, 1]
    noise = rng.normal(size=len(xy))
    Z = np.vstack([smooth, noise, np.full(len(xy), 3.0)])
    return Z - Z.mean(axis=1, keepdims=True), xy


def test_inverse_square_reproduces_the_dense_reference(toy_field):
    Z, xy = toy_field
    Y, _, _, _ = _inverse_square_lag(Z, xy)
    got = _morans_i_from_lag(Z, Y)
    want = _naive_morans_i(Z, xy)
    # Row 2 is constant, so its I is undefined in both; compare the real genes.
    assert np.allclose(got[:2], want[:2], rtol=0, atol=1e-12)


def test_blocking_does_not_change_the_answer(toy_field, monkeypatch):
    """The blocked loop must be an optimisation, not an approximation."""
    import shanuz.spatial.variable_features as vf

    Z, xy = toy_field
    unblocked = _morans_i_from_lag(Z, _inverse_square_lag(Z, xy)[0])
    # One row per block is the most fragmented the loop can be.
    monkeypatch.setattr(vf, "_BLOCK_DOUBLES", 1)
    blocked = _morans_i_from_lag(Z, vf._inverse_square_lag(Z, xy)[0])
    assert np.allclose(unblocked[:2], blocked[:2], rtol=0, atol=1e-12)


def test_weight_moments_match_the_dense_matrix(toy_field):
    """S0/S1/S2 are accumulated blockwise; they must total what W really sums to."""
    from scipy.spatial.distance import squareform, pdist

    Z, xy = toy_field
    d = squareform(pdist(xy))
    with np.errstate(divide="ignore"):
        W = 1.0 / d**2
    np.fill_diagonal(W, 0.0)
    W = W / W.sum(axis=1, keepdims=True)
    sym = W + W.T
    _, s0, s1, s2 = _inverse_square_lag(Z, xy)
    assert s0 == pytest.approx(W.sum(), rel=1e-12)
    assert s1 == pytest.approx(0.5 * (sym**2).sum(), rel=1e-12)
    assert s2 == pytest.approx(((W.sum(1) + W.sum(0)) ** 2).sum(), rel=1e-12)


def test_knn_weights_give_a_different_statistic(toy_field):
    """Guards the default: kNN is close enough to look right and is not R's."""
    Z, xy = toy_field
    genes = ["smooth", "noise", "flat"]
    exact = _moransi_table(Z, xy, genes, k=10, weights="inverse_square")
    knn = _moransi_table(Z, xy, genes, k=10, weights="knn")
    assert exact.loc["smooth", "moransi"] != pytest.approx(
        knn.loc["smooth", "moransi"], rel=1e-6
    )
    # Both should still call the gradient spatial and the noise gene not.
    assert exact.loc["smooth", "moransi"] > exact.loc["noise", "moransi"]


def test_public_default_is_r_weights_not_knn():
    """Pins the *default*, not just the helper.

    An earlier version of this file only ever called ``_moransi_table`` with an
    explicit ``weights=``, so flipping the public default back to ``"knn"`` kept
    the whole file green. Go through the public function or guard nothing.
    """
    import scipy.sparse as sp_

    from shanuz import create_shanuz_object
    from shanuz.spatial import create_fovs, find_spatially_variable_features

    rng = np.random.default_rng(3)
    gx, gy = np.meshgrid(np.arange(10.0), np.arange(10.0))
    xy = np.column_stack([gx.ravel(), gy.ravel()])
    n = len(xy)
    gradient = xy[:, 0] + xy[:, 1]
    counts = np.vstack([gradient, rng.poisson(2.0, size=n), rng.poisson(2.0, size=n)])
    cells = [f"cell_{i}" for i in range(n)]
    feats = ["gradient", "noise1", "noise2"]
    obj = create_shanuz_object(
        sp_.csc_matrix(counts.astype(float)),
        assay="Xenium", feature_names=feats, cell_names=cells,
    )
    obj.images = create_fovs(pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1], "cell": cells}))

    got = find_spatially_variable_features(obj, layer="counts")["moransi"]
    Z = counts - counts.mean(axis=1, keepdims=True)
    want = pd.Series(_naive_morans_i(Z, xy), index=feats)
    knn = _moransi_table(Z, xy, feats, k=10, weights="knn")["moransi"]

    assert got.loc["gradient"] == pytest.approx(want.loc["gradient"], abs=1e-12)
    assert got.loc["gradient"] != pytest.approx(knn.loc["gradient"], rel=1e-6)


def test_coincident_cells_are_reported_not_silently_infinite(toy_field):
    Z, xy = toy_field
    xy = xy.copy()
    xy[1] = xy[0]
    with pytest.raises(ValueError, match="same coordinates"):
        _inverse_square_lag(Z, xy)


# ---------------------------------------------------------------------------
# Centroids — the auto radius
# ---------------------------------------------------------------------------

def _toy_coords():
    return pd.DataFrame(
        {"x": [1.0, 2, 3, 4], "y": [10.0, 20, 30, 40], "cell": list("abcd")}
    )


def test_centroids_auto_radius_matches_seuratobject():
    # R: CreateCentroids(coords) -> slot 'radius' == 0.165
    # .AutoRadius = 0.01 * mean(diff(range(x)), diff(range(y))) = 0.01*mean(3,30)
    assert create_centroids(_toy_coords()).radius() == pytest.approx(0.165)


def test_explicit_radius_still_wins():
    assert create_centroids(_toy_coords(), radius=7.0).radius() == 7.0


def test_fov_radius_is_none_but_the_boundary_carries_it():
    """Matches R, where Radius(FOV) is NULL and Radius(centroids) is not."""
    fov = create_fov(_toy_coords(), type_="centroids", assay="RNA")
    assert fov.radius() is None
    assert _boundary_radius(fov) == pytest.approx(0.165)


def test_the_plot_path_actually_reads_the_boundary_radius():
    """The wiring, not just the helper.

    Asserting ``_boundary_radius`` alone left the caller free to go back to
    ``fov.radius()`` — which is always None — with every test still passing. The
    spot renderer is the thing that was broken, so it is the thing to pin.
    """
    from shanuz.plotting import _spatial_panel, _spot_collection

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fov = create_fov(_toy_coords(), type_="centroids", assay="RNA")
    coords, radius, _ = _spatial_panel(fov, None)
    assert radius == pytest.approx(0.165)

    fig, ax = plt.subplots()
    try:
        coll = _spot_collection(ax, coords, radius, 1.0)
        # None is the silent fallback to a fixed-size scatter — the actual bug.
        assert coll is not None
        assert coll.get_offsets().shape == (4, 2)
    finally:
        plt.close(fig)


# ---------------------------------------------------------------------------
# Segmentation — closed rings
# ---------------------------------------------------------------------------

def test_segmentation_closes_each_ring():
    # R: CreateSegmentation on 4 vertices per cell -> GetTissueCoordinates has 5
    square = pd.DataFrame(
        {
            "x": [0.0, 1, 1, 0, 5, 6, 6, 5],
            "y": [0.0, 0, 1, 1, 5, 5, 6, 6],
            "cell": ["a"] * 4 + ["b"] * 4,
        }
    )
    tc = create_segmentation(square).get_tissue_coordinates()
    assert len(tc) == 10
    for cell in ("a", "b"):
        ring = tc.loc[cell]
        assert ring.iloc[0].tolist() == ring.iloc[-1].tolist()


def test_closing_a_ring_is_idempotent():
    closed = pd.DataFrame(
        {"x": [0.0, 1, 1, 0, 0], "y": [0.0, 0, 1, 1, 0], "cell": ["a"] * 5}
    )
    assert len(create_segmentation(closed).get_tissue_coordinates()) == 5


def test_concave_shapes_survive_closing():
    """An L must not come back as its convex hull — R keeps the notch too."""
    L = pd.DataFrame(
        {"x": [0.0, 2, 2, 1, 1, 0], "y": [0.0, 0, 1, 1, 2, 2], "cell": ["a"] * 6}
    )
    tc = create_segmentation(L).get_tissue_coordinates()
    assert len(tc) == 7  # R: 7
    assert [1.0, 1.0] in tc.values.tolist()


def test_degenerate_rings_are_left_alone():
    """Two vertices are a line, not a ring; closing it would invent an edge."""
    line = pd.DataFrame({"x": [0.0, 1], "y": [0.0, 1], "cell": ["a", "a"]})
    assert len(create_segmentation(line).get_tissue_coordinates()) == 2
