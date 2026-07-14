"""Tests for find_spatially_variable_features(method="markvariogram")."""
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from shanuz import create_shanuz_object, find_spatially_variable_features
from shanuz.preprocessing import normalize_data
from shanuz.spatial.fov import create_fovs
from shanuz.spatial.variable_features import (
    _band_weights,
    _mark_variogram,
    _nn_spacing,
)


def _grid(side: int) -> np.ndarray:
    """A side × side lattice of cells, one unit apart."""
    gx, gy = np.meshgrid(np.arange(side, dtype=float), np.arange(side, dtype=float))
    return np.column_stack([gx.ravel(), gy.ravel()])


def _centre(Z: np.ndarray) -> np.ndarray:
    return Z - Z.mean(axis=1, keepdims=True)


def _brute_gamma(z: np.ndarray, xy: np.ndarray, r: float, h: float) -> float:
    """γ(r) straight from the definition: an explicit loop over every cell pair."""
    num = den = 0.0
    for i in range(len(xy)):
        for j in range(i + 1, len(xy)):
            u = (float(np.hypot(*(xy[i] - xy[j]))) - r) / h
            if abs(u) < 1.0:
                w = 1.0 - u**2
                num += w * 0.5 * (z[i] - z[j]) ** 2
                den += w
    var = float(np.mean((z - z.mean()) ** 2))
    return (num / den) / var


# ---------------------------------------------------------------------------
# Statistic correctness
# ---------------------------------------------------------------------------

def test_mark_variogram_matches_brute_force():
    """The sparse identity equals the explicit pair-by-pair double sum."""
    rng = np.random.default_rng(0)
    xy = rng.uniform(0, 10, (60, 2))
    Z = _centre(rng.normal(size=(4, 60)))
    r, h = 3.0, 1.0

    mine = _mark_variogram(Z, _band_weights(xy, r=r, h=h))
    brute = np.array([_brute_gamma(z, xy, r, h) for z in Z])
    np.testing.assert_allclose(mine, brute)


def test_band_weights_only_hold_pairs_near_r():
    xy = _grid(10)
    r, h = 3.0, 1.0
    K = _band_weights(xy, r=r, h=h)

    assert (K != K.T).nnz == 0              # symmetric
    assert K.diagonal().sum() == 0          # no self-pairs

    D = np.linalg.norm(xy[:, None, :] - xy[None, :, :], axis=-1)
    Kd = K.toarray()
    in_band = np.abs(D - r) < h
    np.fill_diagonal(in_band, False)
    assert (Kd > 0).sum() == in_band.sum()  # exactly the pairs in the band
    # Epanechnikov: weight peaks at d = r and falls to zero at the band edge.
    np.testing.assert_allclose(Kd[in_band], 1.0 - ((D[in_band] - r) / h) ** 2)


def test_unstructured_marks_sit_at_one():
    """With no spatial structure, cells r apart differ as much as any two cells."""
    rng = np.random.default_rng(1)
    xy = _grid(20)
    Z = _centre(rng.normal(size=(20, len(xy))))

    gamma = _mark_variogram(Z, _band_weights(xy, r=3.0, h=1.0))
    assert abs(gamma.mean() - 1.0) < 0.05
    assert np.all(np.abs(gamma - 1.0) < 0.2)


def test_structured_marks_sit_below_one():
    """A smooth spatial gradient still has most of its variance intact at r."""
    xy = _grid(20)
    Z = _centre(np.vstack([xy[:, 0], xy[:, 1]]))

    gamma = _mark_variogram(Z, _band_weights(xy, r=3.0, h=1.0))
    assert np.all(gamma < 0.2)


def test_gamma_ignores_the_scale_of_expression():
    """Normalising by the gene's own variance makes γ shift- and scale-free."""
    rng = np.random.default_rng(2)
    xy = _grid(12)
    z = np.exp(-((xy[:, 0] - 6) ** 2 + (xy[:, 1] - 6) ** 2) / 8) + rng.normal(0, 0.1, len(xy))
    Z = _centre(np.vstack([z, 100.0 * z + 7.0]))

    gamma = _mark_variogram(Z, _band_weights(xy, r=3.0, h=1.0))
    assert np.isclose(gamma[0], gamma[1])


def test_gamma_ignores_the_scale_of_the_coordinates():
    """r is measured in nearest-neighbour spacings, so pixels and microns agree.

    This is where we part company with R, which takes r.metric in raw coordinate
    units and so answers differently on the same slide in different units.
    """
    rng = np.random.default_rng(3)
    xy = _grid(12)
    Z = _centre(rng.normal(size=(5, len(xy))))

    def gamma_at(coords, r_metric=3.0, bandwidth=1.0):
        spacing = _nn_spacing(coords)
        K = _band_weights(coords, r=r_metric * spacing, h=bandwidth * spacing)
        return _mark_variogram(Z, K)

    np.testing.assert_allclose(gamma_at(xy), gamma_at(xy * 1000.0))


def test_flat_gene_has_no_variogram():
    """A gene with no variance has nothing to decay, so γ is undefined."""
    xy = _grid(8)
    Z = _centre(np.vstack([np.zeros(len(xy)), np.arange(len(xy), dtype=float)]))

    gamma = _mark_variogram(Z, _band_weights(xy, r=2.0, h=1.0))
    assert np.isnan(gamma[0])
    assert np.isfinite(gamma[1])


def test_empty_distance_band_raises():
    """Reading the variogram further out than the slide is wide has no answer."""
    xy = _grid(6)
    Z = np.zeros((1, len(xy)))
    with pytest.raises(ValueError, match="No cell pairs"):
        _mark_variogram(Z, _band_weights(xy, r=500.0, h=1.0))


def test_nn_spacing_is_the_lattice_step():
    assert np.isclose(_nn_spacing(_grid(10)), 1.0)
    assert np.isclose(_nn_spacing(_grid(10) * 7.5), 7.5)


# ---------------------------------------------------------------------------
# End-to-end
# ---------------------------------------------------------------------------

N_RAND = 22
RAND_GENES = [f"rand{i}" for i in range(N_RAND)]
STRUCTURED = ["grad_x", "grad_y", "blob"]


@pytest.fixture
def spatial_obj():
    """A 24×24 lattice; grad_x/grad_y/blob are spatially structured, rand* are not.

    The slide has to be a good deal wider than the distance the variogram is read
    at, or ``r_metric=5`` spans a third of the tissue, few pairs land in the band
    and γ gets noisy. The blob is wide for the same reason: a bump only a cell or
    two across is *fully* decorrelated five cells out, and would score as
    unstructured — correctly, but that is not what this fixture is here to show.

    ``silent`` is never detected in any cell — a real thing on a panel, and the
    one gene that cannot be scored at all.
    """
    rng = np.random.default_rng(0)
    xy = _grid(24)
    n = len(xy)
    genes = STRUCTURED + RAND_GENES + ["silent"]

    X = np.zeros((len(genes), n))
    X[0] = xy[:, 0]
    X[1] = xy[:, 1]
    X[2] = np.exp(-((xy[:, 0] - 12) ** 2 + (xy[:, 1] - 12) ** 2) / 50) * 10
    for i in range(N_RAND):
        X[3 + i] = rng.poisson(5, n)
    # X[-1] ("silent") stays all zero.

    obj = create_shanuz_object(
        counts=sp.csc_matrix(np.abs(X) * 3),
        feature_names=genes,
        cell_names=[f"c{i}" for i in range(n)],
    )
    coords = pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1], "cell": obj.cell_names()})
    obj.images = create_fovs(coords, assay="RNA", default_name="rna")
    normalize_data(obj)
    return obj


def test_markvariogram_ranks_structured_genes_first(spatial_obj):
    res = find_spatially_variable_features(spatial_obj, method="markvariogram")

    assert list(res.columns) == ["markvariogram", "markvariogram_rank"]
    assert set(res.index[:3]) == set(STRUCTURED)
    # Structured genes keep most of their variance across five cell-widths;
    # unstructured ones have spent all of it and sit at γ ≈ 1.
    assert (res.loc[STRUCTURED, "markvariogram"] < 0.6).all()
    assert (res.loc[RAND_GENES, "markvariogram"] > 0.8).all()

    assert res["markvariogram_rank"].is_monotonic_increasing
    assert res["markvariogram_rank"].iloc[0] == 1


def test_undetected_gene_is_ranked_last(spatial_obj):
    res = find_spatially_variable_features(spatial_obj, method="markvariogram")
    assert np.isnan(res.loc["silent", "markvariogram"])
    assert res.index[-1] == "silent"
    assert res.loc["silent", "markvariogram_rank"] == len(res)


def test_moransi_also_survives_an_undetected_gene(spatial_obj):
    """The NaN score is ranked, not crashed on — both methods share that path."""
    res = find_spatially_variable_features(spatial_obj, method="moransi", k=8)
    assert np.isnan(res.loc["silent", "moransi"])
    assert res.index[-1] == "silent"


def test_both_methods_agree_on_the_top_gene(spatial_obj):
    """Two different statistics, one slide — they had better pick the same gene."""
    mv = find_spatially_variable_features(spatial_obj, method="markvariogram")
    mi = find_spatially_variable_features(spatial_obj, method="moransi", k=8)
    assert set(mv.index[:3]) == set(mi.index[:3]) == set(STRUCTURED)


def test_results_written_to_feature_metadata(spatial_obj):
    res = find_spatially_variable_features(spatial_obj, method="markvariogram")
    meta = spatial_obj.assays["RNA"].meta_data
    for col in ("markvariogram", "markvariogram_rank"):
        assert col in meta.columns
    assert np.isclose(meta.loc["blob", "markvariogram"], res.loc["blob", "markvariogram"])


def test_features_argument_restricts_output(spatial_obj):
    res = find_spatially_variable_features(
        spatial_obj, features=["blob", "rand0"], method="markvariogram")
    assert set(res.index) == {"blob", "rand0"}
    assert res.index[0] == "blob"


def test_r_metric_reads_the_variogram_further_out(spatial_obj):
    """The blob decorrelates with distance, so γ climbs as r grows."""
    near = find_spatially_variable_features(
        spatial_obj, method="markvariogram", r_metric=2.0)
    far = find_spatially_variable_features(
        spatial_obj, method="markvariogram", r_metric=6.0)
    assert near.loc["blob", "markvariogram"] < far.loc["blob", "markvariogram"]


def test_bandwidth_widens_the_band_without_moving_the_answer(spatial_obj):
    """A wider band averages over more pairs; the gene ranking is unmoved."""
    narrow = find_spatially_variable_features(
        spatial_obj, method="markvariogram", bandwidth=0.6)
    wide = find_spatially_variable_features(
        spatial_obj, method="markvariogram", bandwidth=2.0)
    assert not np.isclose(
        narrow.loc["blob", "markvariogram"], wide.loc["blob", "markvariogram"])
    assert set(narrow.index[:3]) == set(wide.index[:3])


def test_bad_parameters_raise(spatial_obj):
    with pytest.raises(ValueError, match="r_metric must be positive"):
        find_spatially_variable_features(
            spatial_obj, method="markvariogram", r_metric=0.0)
    with pytest.raises(ValueError, match="bandwidth must be positive"):
        find_spatially_variable_features(
            spatial_obj, method="markvariogram", bandwidth=-1.0)


def test_r_metric_beyond_the_slide_raises(spatial_obj):
    with pytest.raises(ValueError, match="No cell pairs"):
        find_spatially_variable_features(
            spatial_obj, method="markvariogram", r_metric=500.0)
