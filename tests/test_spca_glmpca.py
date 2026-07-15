"""Tests for the v0.5.0 reductions: run_spca (reduction.py) and glm_pca (glmpca.py)."""
import sys
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz.shanuz import create_shanuz_object  # noqa: E402
from shanuz.graph import Graph  # noqa: E402
from shanuz.preprocessing import (  # noqa: E402
    normalize_data,
    find_variable_features,
    scale_data,
)
from shanuz.reduction import run_pca, run_spca  # noqa: E402
from shanuz.glmpca import (  # noqa: E402
    glm_pca,
    _counts_for,
    _orthogonalize,
    _poisson_deviance,
    _nb_deviance,
    _estimate_theta,
)
from shanuz.neighbors import find_neighbors  # noqa: E402

pytest.importorskip("sklearn")


def _cosine(a, b):
    return abs(float(a @ b)) / (np.linalg.norm(a) * np.linalg.norm(b))


def _silhouette(embedding, labels):
    from sklearn.metrics import silhouette_score
    return float(silhouette_score(embedding, labels))


# ---------------------------------------------------------------------------
# run_spca
# ---------------------------------------------------------------------------

@pytest.fixture
def clustered():
    """90 cells in 3 clusters, each with its own block of elevated genes."""
    rng = np.random.default_rng(0)
    n_genes, per, k = 90, 30, 3
    n = per * k
    lam = np.full((n_genes, n), 3.0)
    for c in range(n):
        cluster = c // per
        lam[cluster * 20:(cluster + 1) * 20, c] += 8.0
    obj = create_shanuz_object(
        counts=sp.csc_matrix(rng.poisson(lam).astype(float)),
        assay="RNA",
        feature_names=[f"g{i}" for i in range(n_genes)],
        cell_names=[f"c{i}" for i in range(n)],
    )
    normalize_data(obj)
    find_variable_features(obj, selection_method="vst", nfeatures=60)
    scale_data(obj)
    run_pca(obj, n_pcs=10)
    labels = np.repeat(np.arange(k), per)
    return obj, labels


def _identity_graph(obj):
    cells = obj.cell_names()
    return Graph(sp.identity(len(cells), format="csc"), cell_names=cells,
                 assay_used="RNA")


def test_spca_with_an_identity_graph_is_pca(clustered):
    """G = I turns vᵀXᵀGXv back into vᵀXᵀXv, so the loadings must be PCA's."""
    obj, _ = clustered
    obj.graphs["identity"] = _identity_graph(obj)
    run_spca(obj, graph="identity", npcs=10)

    pca = obj.reductions["pca"].feature_loadings
    spca = obj.reductions["spca"].feature_loadings
    assert spca.shape == pca.shape
    # Eigenvector signs are arbitrary in both, so compare directions.
    for i in range(5):
        assert _cosine(spca[:, i], pca[:, i]) > 0.999


def test_spca_finds_what_the_graph_knows_and_pca_misses(clustered):
    """The supervision is the whole point: hand sPCA a grouping PCA can't see.

    Thirty coherent nuisance genes outvote ten genes carrying the real grouping,
    so PC1 tracks the nuisance. The graph knows the groups, and sPCA — which
    maximises vᵀXᵀGXv, i.e. the separation of the graph's blocks — puts them on
    SPC1.

    The nuisance axis runs through both groups equally — think cell-cycle phase
    spanning two cell types. That is what makes it a nuisance rather than a
    confounder, and it is the case sPCA is meant for: the strongest signal in the
    data is real, and is not the one you asked about.

    The grouping is deliberately *library-size neutral*: fifteen genes go up in
    one group and fifteen in the other, so the two groups differ in which genes
    are on rather than in how much RNA they hold. Let the grouping shift total
    counts instead and `normalize_data` divides that shift back out of every
    gene, smearing the grouping across the whole matrix — including the
    high-variance nuisance genes, which then drag their own variance onto the
    axis sPCA returns. It is the same compositional leak documented for Moran's I
    in `find_spatially_variable_features`, and it blunts the separation badly
    enough to hide the effect being tested.
    """
    rng = np.random.default_rng(1)
    per, n_genes = 50, 65
    n = 2 * per
    group = np.repeat([0, 1], per)
    nuisance = np.tile(np.linspace(0.0, 1.0, per), 2)   # identical within each group

    lam = np.full((n_genes, n), 3.0)
    lam[0:30] += 12.0 * nuisance                     # 30 genes track the nuisance
    lam[30:45] += 3.0 * group                        # 15 genes up in group 1
    lam[45:60] += 3.0 * (1 - group)                  # 15 genes up in group 0
    #  the remaining 5 genes are noise

    obj = create_shanuz_object(
        counts=sp.csc_matrix(rng.poisson(lam).astype(float)),
        assay="RNA",
        feature_names=[f"g{i}" for i in range(n_genes)],
        cell_names=[f"c{i}" for i in range(n)],
    )
    normalize_data(obj)
    find_variable_features(obj, selection_method="vst", nfeatures=n_genes)
    scale_data(obj)
    run_pca(obj, n_pcs=10)

    # A graph that knows the groups and nothing else: connected within, not across.
    block = (group[:, None] == group[None, :]).astype(float)
    obj.graphs["groups"] = Graph(sp.csc_matrix(block), cell_names=obj.cell_names(),
                                 assay_used="RNA")
    run_spca(obj, graph="groups", npcs=10)

    pc1 = obj.reductions["pca"].cell_embeddings[:, :1]
    spc1 = obj.reductions["spca"].cell_embeddings[:, :1]
    assert _silhouette(spc1, group) > 0.7       # sPCA puts the grouping first
    assert _silhouette(pc1, group) < 0.2        # PCA's leading axis is the nuisance


def test_spca_embeddings_feed_downstream_steps(clustered):
    obj, labels = clustered
    obj.graphs["identity"] = _identity_graph(obj)
    run_spca(obj, graph="identity", npcs=10)

    dr = obj.reductions["spca"]
    assert dr.cell_embeddings.shape == (len(labels), 10)
    assert dr.feature_loadings.shape[1] == 10
    assert np.isfinite(dr.cell_embeddings).all()
    assert dr.misc["spca_graph"] == "identity"
    assert (np.diff(dr.misc["eigenvalues"]) <= 1e-8).all()   # descending

    find_neighbors(obj, reduction="spca", graph_name="spca")
    assert "spca_snn" in obj.graphs


def test_spca_is_deterministic(clustered):
    obj, _ = clustered
    obj.graphs["identity"] = _identity_graph(obj)
    run_spca(obj, graph="identity", npcs=8, reduction_name="a")
    run_spca(obj, graph="identity", npcs=8, reduction_name="b")
    np.testing.assert_allclose(
        obj.reductions["a"].cell_embeddings, obj.reductions["b"].cell_embeddings)


def test_spca_without_a_graph_raises(clustered):
    obj, _ = clustered
    with pytest.raises(KeyError, match="not found"):
        run_spca(obj, graph="wsnn")


def test_spca_with_a_wrongly_sized_graph_raises(clustered):
    obj, _ = clustered
    obj.graphs["small"] = Graph(
        sp.identity(5, format="csc"),
        cell_names=obj.cell_names()[:5], assay_used="RNA")
    with pytest.raises(ValueError, match="cells"):
        run_spca(obj, graph="small")


# ---------------------------------------------------------------------------
# glm_pca
# ---------------------------------------------------------------------------

@pytest.fixture
def counts_obj():
    """60 cells in 3 clusters of raw Poisson counts — no normalisation at all.

    ``silent`` is detected in no cell: a real thing on a panel, and the gene most
    likely to produce a log(0) somewhere it shouldn't.
    """
    rng = np.random.default_rng(0)
    n_genes, per, k = 40, 20, 3
    n = per * k
    lam = np.full((n_genes, n), 2.0)
    for c in range(n):
        cluster = c // per
        lam[cluster * 10:(cluster + 1) * 10, c] += 10.0
    counts = rng.poisson(lam).astype(float)
    counts[-1, :] = 0.0                                  # the silent gene

    obj = create_shanuz_object(
        counts=sp.csc_matrix(counts),
        assay="RNA",
        feature_names=[f"g{i}" for i in range(n_genes - 1)] + ["silent"],
        cell_names=[f"c{i}" for i in range(n)],
    )
    return obj, np.repeat(np.arange(k), per)


def test_glmpca_recovers_clusters_from_raw_counts(counts_obj):
    obj, labels = counts_obj
    glm_pca(obj, n_components=3, seed=0)

    dr = obj.reductions["glmpca"]
    assert dr.cell_embeddings.shape == (len(labels), 3)
    assert np.isfinite(dr.cell_embeddings).all()
    assert _silhouette(dr.cell_embeddings, labels) > 0.5


def test_glmpca_deviance_only_ever_falls(counts_obj):
    """Backtracking rejects any step that doesn't improve the fit."""
    obj, _ = counts_obj
    glm_pca(obj, n_components=3, seed=0)

    deviance = obj.reductions["glmpca"].misc["deviance"]
    assert len(deviance) > 1
    assert (np.diff(deviance) <= 0).all()
    assert np.isfinite(deviance).all()


def test_glmpca_actually_fits_rather_than_stalling(counts_obj):
    """The saddle guard: U = V = 0 is stationary, so a fit started there goes
    nowhere while *looking* converged — the factors still orient themselves on the
    first step, so clusters show up in a plot even though the deviance has barely
    moved off its null value. Insist the fit does real work.
    """
    obj, _ = counts_obj
    glm_pca(obj, n_components=3, seed=0)
    dr = obj.reductions["glmpca"]

    deviance = dr.misc["deviance"]
    assert len(deviance) > 5                      # not "converged" on step one
    assert deviance[-1] < 0.9 * deviance[0]       # and the counts are better explained
    assert dr.misc["converged"]


def test_glmpca_survives_a_gene_detected_in_no_cell(counts_obj):
    """log(0) is the obvious way to blow this up; the intercept floors instead."""
    obj, _ = counts_obj
    glm_pca(obj, n_components=3, seed=0)

    dr = obj.reductions["glmpca"]
    assert np.isfinite(dr.feature_loadings).all()
    names = dr.features()
    silent = names.index("silent")
    others = [i for i in range(len(names)) if i != silent]
    # With no counts there is no information, so the ridge pulls it to zero: the
    # silent gene must not end up as a driver of any component.
    silent_weight = np.abs(dr.feature_loadings[silent]).max()
    assert silent_weight < np.abs(dr.feature_loadings[others]).max()
    assert silent_weight < 1e-3


def test_glmpca_is_deterministic(counts_obj):
    obj, _ = counts_obj
    glm_pca(obj, n_components=3, seed=7, reduction_name="a")
    glm_pca(obj, n_components=3, seed=7, reduction_name="b")
    np.testing.assert_allclose(
        obj.reductions["a"].cell_embeddings, obj.reductions["b"].cell_embeddings)


def test_glmpca_embeddings_feed_downstream_steps(counts_obj):
    obj, _ = counts_obj
    glm_pca(obj, n_components=3, seed=0)
    find_neighbors(obj, reduction="glmpca", graph_name="glmpca")
    assert "glmpca_snn" in obj.graphs


def test_glmpca_rejects_an_unknown_family(counts_obj):
    obj, _ = counts_obj
    with pytest.raises(NotImplementedError, match="not implemented"):
        glm_pca(obj, family="gaussian")


def test_glmpca_rejects_a_cell_with_no_counts():
    counts = np.ones((5, 4))
    counts[:, 2] = 0.0                                   # an empty droplet
    obj = create_shanuz_object(
        counts=sp.csc_matrix(counts), assay="RNA",
        feature_names=[f"g{i}" for i in range(5)],
        cell_names=[f"c{i}" for i in range(4)],
    )
    with pytest.raises(ValueError, match="zero total counts"):
        glm_pca(obj, n_components=2)


# ---------------------------------------------------------------------------
# glm_pca — negative binomial family
# ---------------------------------------------------------------------------

def _nb_counts(rng, mean, theta):
    """Draw NB counts with the given mean and dispersion (Var = μ + μ²/θ)."""
    p = theta / (theta + mean)
    return rng.negative_binomial(theta, p).astype(float)


@pytest.fixture
def overdispersed_obj():
    """60 cells in 3 clusters of *over-dispersed* counts — NB territory.

    Same cluster layout as ``counts_obj`` but drawn from NB(μ, θ=4) rather than
    Poisson, so the per-gene noise is well above what a Poisson fit expects. The
    dispersion is deliberately low enough to matter.
    """
    rng = np.random.default_rng(0)
    n_genes, per, k = 40, 20, 3
    n = per * k
    mean = np.full((n_genes, n), 2.0)
    for c in range(n):
        cluster = c // per
        mean[cluster * 10:(cluster + 1) * 10, c] += 10.0
    counts = _nb_counts(rng, mean, theta=4.0)

    obj = create_shanuz_object(
        counts=sp.csc_matrix(counts),
        assay="RNA",
        feature_names=[f"g{i}" for i in range(n_genes)],
        cell_names=[f"c{i}" for i in range(n)],
    )
    return obj, np.repeat(np.arange(k), per)


def test_glmpca_nb_recovers_clusters_from_overdispersed_counts(overdispersed_obj):
    obj, labels = overdispersed_obj
    glm_pca(obj, n_components=3, family="nb", seed=0)

    dr = obj.reductions["glmpca"]
    assert dr.cell_embeddings.shape == (len(labels), 3)
    assert np.isfinite(dr.cell_embeddings).all()
    assert _silhouette(dr.cell_embeddings, labels) > 0.5
    assert dr.misc["glmpca_family"] == "nb"


def test_glmpca_nb_reduces_to_poisson_at_large_fixed_theta(counts_obj):
    """θ → ∞ is the Poisson limit: the NB fit must land on the Poisson fit.

    Same counts, same seed; a negative binomial with a huge fixed dispersion has
    almost no ``μ²/θ`` term left, so every Fisher step matches its Poisson twin.
    """
    obj, _ = counts_obj
    glm_pca(obj, n_components=3, family="poisson", seed=0, reduction_name="pois")
    glm_pca(obj, n_components=3, family="nb", theta=1e8,
            optimize_theta=False, seed=0, reduction_name="nb")

    pois = obj.reductions["pois"].cell_embeddings
    nb = obj.reductions["nb"].cell_embeddings
    for i in range(3):
        assert _cosine(nb[:, i], pois[:, i]) > 0.999


def test_glmpca_nb_deviance_falls_with_theta_fixed(overdispersed_obj):
    """A moving θ re-scales the deviance, so monotonicity is only promised when
    θ is pinned. Hold it fixed and the backtracking guarantee is back."""
    obj, _ = overdispersed_obj
    glm_pca(obj, n_components=3, family="nb", theta=4.0,
            optimize_theta=False, seed=0)

    deviance = obj.reductions["glmpca"].misc["deviance"]
    assert len(deviance) > 1
    assert (np.diff(deviance) <= 1e-9).all()
    assert np.isfinite(deviance).all()
    assert deviance[-1] < 0.9 * deviance[0]


def test_glmpca_nb_learns_more_dispersion_on_noisier_data(counts_obj,
                                                          overdispersed_obj):
    """Estimated θ is a read-out of the noise: small when counts are over-dispersed,
    large (toward the Poisson limit) when they are merely Poisson."""
    pois_obj, _ = counts_obj
    over_obj, _ = overdispersed_obj
    glm_pca(pois_obj, n_components=3, family="nb", seed=0)
    glm_pca(over_obj, n_components=3, family="nb", seed=0)

    theta_pois = pois_obj.reductions["glmpca"].misc["theta"]
    theta_over = over_obj.reductions["glmpca"].misc["theta"]
    assert 0 < theta_over < theta_pois              # the noisier data pins θ lower


def test_glmpca_nb_stores_the_fitted_theta(overdispersed_obj, counts_obj):
    over_obj, _ = overdispersed_obj
    glm_pca(over_obj, n_components=3, family="nb", seed=0)
    theta = over_obj.reductions["glmpca"].misc["theta"]
    assert np.isfinite(theta) and theta > 0

    # Poisson is NB at θ = ∞; record that so misc["theta"] always means something.
    pois_obj, _ = counts_obj
    glm_pca(pois_obj, n_components=3, family="poisson", seed=0)
    assert pois_obj.reductions["glmpca"].misc["theta"] == np.inf


def test_glmpca_nb_is_deterministic(overdispersed_obj):
    obj, _ = overdispersed_obj
    glm_pca(obj, n_components=3, family="nb", seed=7, reduction_name="a")
    glm_pca(obj, n_components=3, family="nb", seed=7, reduction_name="b")
    np.testing.assert_allclose(
        obj.reductions["a"].cell_embeddings, obj.reductions["b"].cell_embeddings)


# --- the pieces, on their own ---------------------------------------------

def test_orthogonalize_reorders_without_changing_the_fit():
    """The rotation is exact: U·Vᵀ — all the model ever sees — is untouched."""
    rng = np.random.default_rng(3)
    U = rng.normal(size=(20, 4))
    V = rng.normal(size=(15, 4))

    loadings, factors = _orthogonalize(U, V)
    np.testing.assert_allclose(loadings @ factors.T, U @ V.T, atol=1e-10)

    # Loadings orthonormal; factors ordered biggest-first, as PCA readers expect.
    np.testing.assert_allclose(loadings.T @ loadings, np.eye(4), atol=1e-10)
    scale = np.linalg.norm(factors, axis=0)
    assert (np.diff(scale) <= 1e-8).all()


def test_poisson_deviance_is_zero_for_a_perfect_fit():
    Y = np.array([[0.0, 3.0], [7.0, 2.0]])
    assert _poisson_deviance(Y, Y) == pytest.approx(0.0)     # y = 0 must not blow up
    assert _poisson_deviance(Y, Y + 1.0) > 0.0


def test_nb_deviance_is_zero_for_a_perfect_fit():
    Y = np.array([[0.0, 3.0], [7.0, 2.0]])
    assert _nb_deviance(Y, Y, theta=5.0) == pytest.approx(0.0)   # y = 0 stays finite
    assert _nb_deviance(Y, Y + 1.0, theta=5.0) > 0.0


def test_nb_deviance_approaches_poisson_as_theta_grows():
    """The μ²/θ correction vanishes as θ → ∞, leaving the Poisson deviance."""
    rng = np.random.default_rng(4)
    Y = rng.poisson(3.0, size=(8, 6)).astype(float)
    mu = np.full_like(Y, 3.0)
    assert _nb_deviance(Y, mu, theta=1e10) == pytest.approx(
        _poisson_deviance(Y, mu), rel=1e-4)


def test_estimate_theta_recovers_a_known_dispersion():
    """Hand the estimator the true mean and a big sample; it should find θ."""
    rng = np.random.default_rng(5)
    theta_true, mu = 3.0, 5.0
    p = theta_true / (theta_true + mu)
    Y = rng.negative_binomial(theta_true, p, size=(200, 200)).astype(float)
    mu_mat = np.full_like(Y, mu)

    theta_hat = _estimate_theta(Y, mu_mat, theta=100.0)     # start far from truth
    assert theta_hat == pytest.approx(theta_true, rel=0.2)


def test_estimate_theta_stays_in_range_for_poisson_like_data():
    """No over-dispersion to find: θ climbs toward the ceiling, never negative."""
    rng = np.random.default_rng(6)
    Y = rng.poisson(4.0, size=(100, 100)).astype(float)
    mu = np.full_like(Y, 4.0)
    theta_hat = _estimate_theta(Y, mu, theta=10.0)
    assert theta_hat > 10.0                                 # pushed up, not down
    assert np.isfinite(theta_hat)


def test_counts_for_rejects_negative_values():
    """Handing GLM-PCA scaled data is the easy mistake; it should say so."""
    scaled = np.array([[-1.2, 0.4], [0.9, -0.3]])

    def getter(assay_obj, layer):
        return scaled, ["a", "b"]

    with pytest.raises(ValueError, match="negative values"):
        _counts_for(None, ["a", "b"], "scale.data", getter)
