"""``run_pca`` against ``RunPCA``'s irlba — solver accuracy, sdev, centring.

The suite was green across the solver change that took PC agreement with Seurat
from 15/30 to 30/30, so nothing here was pinning it. These tests do.

Every fixture below is deliberately **large enough to reach the branch under
test**: sklearn's ``PCA`` only switches to its randomized solver once
``max(shape) > 500`` and ``n_components < 0.8 * min(shape)``, so a small
fixture runs exact LAPACK either way and a "revert to sklearn" mutation
changes nothing at all. ``test_fixture_is_in_the_randomized_regime`` asserts
that precondition directly rather than leaving it to a comment.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from truecell.preprocessing import normalize_data, scale_data  # noqa: E402
from truecell.reduction import run_pca  # noqa: E402
from truecell.truecell import create_truecell_object  # noqa: E402

N_CELLS, N_GENES, N_PCS = 900, 700, 30


def _object(n_cells=N_CELLS, n_genes=N_GENES, seed=0):
    """Structured counts — a flat random matrix has a degenerate spectrum,
    which would make trailing-PC comparisons meaningless."""
    rng = np.random.default_rng(seed)
    loads = rng.normal(size=(n_genes, 12))
    scores = rng.normal(size=(12, n_cells)) * np.linspace(4, 1, 12)[:, None]
    lam = np.exp(loads @ scores * 0.35 + 1.5)
    counts = rng.poisson(lam).astype(float)
    genes = [f"g{i}" for i in range(n_genes)]
    cells = [f"c{i}" for i in range(n_cells)]
    obj = create_truecell_object(
        counts=counts, assay="RNA", feature_names=genes, cell_names=cells
    )
    normalize_data(obj)
    obj.get_assay().variable_features = genes
    scale_data(obj, features=genes)
    return obj, genes


def _scaled(obj, genes):
    mat = obj.get_assay().layer_data("scale.data", features=genes)
    return np.asarray(mat.todense()) if hasattr(mat, "todense") else np.asarray(mat)


def _exact(obj, genes, n_pcs=N_PCS):
    """``irlba(A = t(scale.data), nv = npcs)`` computed exactly."""
    u, d, vt = np.linalg.svd(_scaled(obj, genes).T, full_matrices=False)
    return u[:, :n_pcs] * d[:n_pcs], d[:n_pcs], vt[:n_pcs].T


def test_fixture_is_in_the_randomized_regime():
    """The guards below are only meaningful if sklearn would have gone
    randomized on this shape — otherwise reverting the solver is a no-op."""
    from sklearn.decomposition import PCA

    obj, genes = _object()
    data_t = _scaled(obj, genes).T
    p = PCA(n_components=N_PCS, random_state=42).fit(data_t)
    assert p._fit_svd_solver == "randomized", (
        f"fixture {data_t.shape} does not reach sklearn's randomized solver "
        f"(got {p._fit_svd_solver!r}); a solver mutation would be invisible"
    )


def test_every_pc_matches_an_exact_svd_not_just_the_leading_ones():
    """sklearn's randomized solver is accurate early and drifts late: on ifnb
    it reproduced 15 of 30 PCs above |r| = 0.99 against Seurat, with PC 28 at
    0.006. Trailing PCs are used — dims 1:30 is the default everywhere."""
    obj, genes = _object()
    run_pca(obj, n_pcs=N_PCS, features=genes, seed=42)
    got = np.asarray(obj.reductions["pca"].cell_embeddings)
    want, _d, _v = _exact(obj, genes)

    corr = np.array([
        abs(np.corrcoef(got[:, i], want[:, i])[0, 1]) for i in range(N_PCS)
    ])
    assert corr.min() > 1 - 1e-6, (
        f"{int((corr <= 1 - 1e-6).sum())} of {N_PCS} PCs drift from the exact "
        f"SVD; worst |r| = {corr.min():.6f} at PC {int(corr.argmin()) + 1}"
    )


def test_trailing_pc_variance_is_not_under_captured():
    """The randomized solver's signature is a *shrunken* trailing spectrum —
    a correlation check alone can miss a systematic scale error."""
    obj, genes = _object()
    run_pca(obj, n_pcs=N_PCS, features=genes, seed=42)
    got = np.asarray(obj.reductions["pca"].cell_embeddings)
    want, _d, _v = _exact(obj, genes)
    ratio = got[:, -1].std() / want[:, -1].std()
    assert abs(ratio - 1) < 1e-6, f"last PC's SD is {ratio:.5f}x the exact one"


def test_stdev_is_seurats_singular_value_formula():
    """``sdev <- d / sqrt(ncol(object) - 1)``, not the SD of the embedding.

    The two agree *exactly* on perfectly centred scale.data: if every gene has
    mean zero then ``1`` is orthogonal to the column space, so each left
    singular vector is already mean-zero. They part company as soon as
    centring is imperfect — which is the normal case, because ``ScaleData``
    clips at ``scale.max = 10`` *after* scaling. On ifnb they differ by 1.2e-4
    relative. So this fixture is deliberately un-centred; a centred one cannot
    tell the two formulas apart and the guard would be decorative.
    """
    obj, genes = _object()
    obj.get_assay().layers["scale.data"] = (
        _scaled(obj, genes) + np.linspace(1.0, 5.0, len(genes))[:, None]
    )
    run_pca(obj, n_pcs=N_PCS, features=genes, seed=42)

    got = np.asarray(obj.reductions["pca"].stdev)
    _e, d, _v = _exact(obj, genes)
    want = d / np.sqrt(len(obj.cell_names()) - 1)
    assert np.allclose(got, want, rtol=1e-9, atol=0)

    embedding_sd = np.sqrt(np.var(np.asarray(
        obj.reductions["pca"].cell_embeddings), axis=0, ddof=1))
    assert not np.allclose(got, embedding_sd, rtol=1e-9, atol=0), (
        "the two formulas agree on this fixture, so it cannot distinguish them"
    )


def test_run_pca_does_not_recentre_the_scaled_data():
    """``RunPCA.default`` hands ``t(scale.data)`` straight to irlba. sklearn's
    ``PCA`` subtracts the column means first — usually a near no-op on
    scale.data, but not when centring was skipped upstream."""
    obj, genes = _object()
    assay = obj.get_assay()
    shifted = _scaled(obj, genes) + np.linspace(1.0, 5.0, len(genes))[:, None]
    assay.layers["scale.data"] = shifted

    run_pca(obj, n_pcs=N_PCS, features=genes, seed=42)
    got = np.asarray(obj.reductions["pca"].cell_embeddings)
    u, d, _vt = np.linalg.svd(shifted.T, full_matrices=False)
    want = u[:, :N_PCS] * d[:N_PCS]
    corr = np.array([
        abs(np.corrcoef(got[:, i], want[:, i])[0, 1]) for i in range(N_PCS)
    ])
    assert corr.min() > 1 - 1e-6, (
        f"embedding does not match an un-centred SVD (worst |r| {corr.min():.6f})"
    )


def test_repeated_runs_are_bit_identical():
    """ARPACK is only deterministic because v0 is seeded."""
    obj, genes = _object()
    run_pca(obj, n_pcs=N_PCS, features=genes, seed=42, reduction_name="a")
    run_pca(obj, n_pcs=N_PCS, features=genes, seed=42, reduction_name="b")
    assert np.array_equal(
        np.asarray(obj.reductions["a"].cell_embeddings),
        np.asarray(obj.reductions["b"].cell_embeddings),
    )


def test_small_inputs_take_the_exact_path_and_still_agree():
    """ARPACK misbehaves when k approaches the rank, so small matrices route
    to a dense SVD. That branch needs its own check."""
    obj, genes = _object(n_cells=60, n_genes=40, seed=3)
    run_pca(obj, n_pcs=8, features=genes, seed=42)
    got = np.asarray(obj.reductions["pca"].cell_embeddings)
    want, _d, _v = _exact(obj, genes, n_pcs=8)
    corr = [abs(np.corrcoef(got[:, i], want[:, i])[0, 1]) for i in range(8)]
    assert min(corr) > 1 - 1e-9


@pytest.mark.parametrize("sparse", [False, True])
def test_sparse_and_dense_scale_data_give_the_same_embedding(sparse):
    """The old code ran ``TruncatedSVD`` for sparse and ``PCA`` for dense — two
    different algorithms behind one function."""
    import scipy.sparse as sp

    obj, genes = _object(n_cells=600, n_genes=520, seed=5)
    dense = _scaled(obj, genes)
    obj.get_assay().layers["scale.data"] = sp.csr_matrix(dense) if sparse else dense
    run_pca(obj, n_pcs=20, features=genes, seed=42)
    got = np.asarray(obj.reductions["pca"].cell_embeddings)
    u, d, _vt = np.linalg.svd(dense.T, full_matrices=False)
    want = u[:, :20] * d[:20]
    corr = [abs(np.corrcoef(got[:, i], want[:, i])[0, 1]) for i in range(20)]
    assert min(corr) > 1 - 1e-6
