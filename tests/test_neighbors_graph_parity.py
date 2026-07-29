"""``find_neighbors`` / ``find_clusters`` graphs against Seurat's own.

Every number pinned here was read off a live Seurat 5.5.1 run over the ifnb
RPCA embedding (13,999 cells, ``k.param = 20``, ``prune.SNN = 1/15``,
``nn.method = "rann"`` so the neighbour table is exact on both sides):

    nn   nnz = 279980 (= n*k)   row sums all 20   col sums [68, 24, 21, 23, 32]
    snn  nnz = 1120457   sum = 156062.938571   diagonal = 1 on all 13999 cells

The suite was green across all four fixes below, so none of it was pinning
these. The graphs are stored objects users read directly, and three of the four
defects were invisible to ``find_clusters`` — ``_sparse_to_igraph`` takes the
strict upper triangle, so it silently discarded the very diagonal that was
missing.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from truecell.clustering import _group_singletons  # noqa: E402
from truecell.neighbors import _build_snn, _knn_to_sparse  # noqa: E402


def _ranked(n=60, k=20, seed=0):
    """A neighbour table shaped like Seurat's: self first, k entries per row."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        others = rng.choice([j for j in range(n) if j != i], size=k - 1, replace=False)
        rows.append([i, *others.tolist()])
    return np.asarray(rows)


# ------------------------------------------------------------------ nn graph

def test_nn_graph_is_directed_with_exactly_k_entries_per_row():
    """Seurat's ``nn`` is ``sparseMatrix(i, j, x = 1)`` over the ranked table —
    no symmetrisation. ``mat + mat.T`` inflated nnz well past n*k."""
    n, k = 60, 20
    nn = _knn_to_sparse(_ranked(n, k), n).tocsr()
    assert nn.nnz == n * k, f"nnz {nn.nnz} != n*k {n * k} — graph was symmetrised"
    assert set(np.asarray(nn.sum(axis=1)).ravel().tolist()) == {float(k)}
    assert set(np.unique(nn.data).tolist()) == {1.0}


def test_nn_graph_keeps_the_in_degree_signal():
    """Symmetrising does not just change counts, it erases which cells are
    hubs: on ifnb, Seurat's column sums run 21–68 while every row sums to 20."""
    n, k = 60, 20
    nn = _knn_to_sparse(_ranked(n, k), n).tocsr()
    col = np.asarray(nn.sum(axis=0)).ravel()
    assert col.min() != col.max(), "column sums are constant — in-degree lost"
    assert not np.allclose(nn.toarray(), nn.toarray().T), "graph is symmetric"


# ----------------------------------------------------------------- snn graph

def test_snn_keeps_the_self_diagonal():
    """``ComputeSNN`` stores SNN[i,i] = k/(2k-k) = 1 for every cell; all 13999
    were present in Seurat's graph and none in truecell's."""
    n, k = 60, 20
    snn = _build_snn(_ranked(n, k), n, k, 1 / 15).tocsr()
    assert (snn.diagonal() != 0).sum() == n
    assert np.allclose(snn.diagonal(), 1.0)


def test_snn_jaccard_is_computed_in_double():
    """float32 put ~3e-08 on every weight. Seurat computes in double, and the
    weights are exact ratios of small integers, so the target is exactness."""
    n, k = 60, 20
    idx = _ranked(n, k)
    snn = _build_snn(idx, n, k, 0.0).tocsr()
    member = np.zeros((n, n))
    member[np.repeat(np.arange(n), k), idx.ravel()] = 1.0
    inter = member @ member.T
    exact = inter / (2 * k - inter)
    assert abs(snn.toarray() - exact).max() < 1e-15


def test_snn_prune_boundary_is_inclusive():
    """``ComputeSNN`` zeroes a weight only when ``it.value() < prune``, so an
    edge sitting exactly on the threshold survives.

    Note this is *not* a precision guard: float32 would pass it too, because a
    weak Python ``prune_snn`` is cast down to float32 for the comparison and
    both sides round identically. Which edges survive never depended on the
    dtype — only the stored weights did, which is what the test above pins.
    """
    k, n = 25, 40
    idx = np.zeros((n, k), dtype=int)
    idx[0] = [0, *range(2, 26)]                       # NN(0) = {0} ∪ {2..25}
    idx[1] = [1, *range(2, 23), 26, 27, 28]           # NN(1) = {1} ∪ {2..22} ∪ {26,27,28}
    for i in range(2, n):                             # |NN(0) ∩ NN(1)| = 21
        idx[i] = [i, *[(i + o) % n for o in range(1, k)]]

    snn = _build_snn(idx, n, k, 21 / 29).tocsr()
    assert snn[0, 1] == pytest.approx(21 / 29, abs=0, rel=1e-15)
    assert _build_snn(idx, n, k, 22 / 29).tocsr()[0, 1] == 0.0, "prune is not exclusive below"


# ------------------------------------------------------------ GroupSingletons

def _toy():
    ids = np.array(["0", "0", "0", "1", "1", "2"])   # "2" holds one cell
    snn = sp.csr_matrix(np.array([
        [1.0, .5, .5, .1, .1, .2],
        [.5, 1.0, .5, .1, .1, .2],
        [.5, .5, 1.0, .1, .1, .2],
        [.1, .1, .1, 1.0, .5, .9],
        [.1, .1, .1, .5, 1.0, .9],
        [.2, .2, .2, .9, .9, 1.0]]))
    return ids, snn


def test_group_singletons_absorbs_into_the_best_connected_cluster():
    """Seurat scores each candidate by mean SNN weight. Cell 5 averages 0.9 to
    cluster "1" and 0.2 to cluster "0"."""
    ids, snn = _toy()
    out = _group_singletons(ids, snn, True)
    assert out[5] == "1"
    assert (out == "1").sum() == 3 and "2" not in set(out.tolist())


def test_group_singletons_scores_by_mean_not_by_sum():
    """``sum(subSNN)/(nrow*ncol)`` — a bare sum would hand the singleton to
    whichever cluster is merely *larger*. The two rules must disagree here or
    the guard is decorative: cluster "0" wins on sum (1.2 > 1.0) and loses on
    mean (0.3 < 0.5), so only a mean-scorer answers "1"."""
    ids = np.array(["0", "0", "0", "0", "1", "1", "2"])
    w = np.zeros((7, 7))
    w[6, [0, 1, 2, 3]] = w[[0, 1, 2, 3], 6] = 0.30   # 4 cells: sum 1.2, mean 0.30
    w[6, [4, 5]] = w[[4, 5], 6] = 0.50               # 2 cells: sum 1.0, mean 0.50
    np.fill_diagonal(w, 1.0)
    assert _group_singletons(ids, sp.csr_matrix(w), True)[6] == "1"


def test_group_singletons_false_pools_them_without_truncation():
    """``ids`` is a width-1 unicode array here, so a naive assignment would
    store "s" instead of "singleton"."""
    ids, snn = _toy()
    out = _group_singletons(ids, snn, False)
    assert out[5] == "singleton", f"got {out[5]!r}"


def test_group_singletons_is_a_no_op_without_singletons():
    ids = np.array(["0", "0", "1", "1"])
    _, snn = _toy()
    assert np.array_equal(_group_singletons(ids, snn, True), ids)


def test_a_singleton_cannot_absorb_another_singleton():
    """Seurat fixes the candidate list *before* the loop, so one singleton
    never becomes another's target."""
    ids = np.array(["0", "0", "0", "1", "2"])
    w = np.full((5, 5), 0.05)
    w[4, 3] = w[3, 4] = 0.99        # the two singletons are each other's nearest
    snn = sp.csr_matrix(w)
    out = _group_singletons(ids, snn, True)
    assert set(out.tolist()) == {"0"}, f"singletons merged into each other: {out}"


# ----------------------------------------------------------- the UMAP consumer

def test_run_umap_zeroes_the_graph_diagonal(monkeypatch):
    """``RunUMAP.Graph`` opens with ``diag(x = object) <- 0``. Now that the SNN
    carries a diagonal of 1, failing to strip it feeds the layout n self-edges."""
    import umap.umap_ as uu

    from truecell.neighbors import find_neighbors
    from truecell.preprocessing import normalize_data, scale_data
    from truecell.reduction import run_pca
    from truecell.truecell import create_truecell_object
    from truecell.umap import run_umap

    rng = np.random.default_rng(0)
    genes = [f"g{i}" for i in range(40)]
    obj = create_truecell_object(
        counts=rng.poisson(4.0, size=(40, 80)).astype(float), assay="RNA",
        feature_names=genes, cell_names=[f"c{i}" for i in range(80)])
    normalize_data(obj)
    obj.get_assay().variable_features = genes
    scale_data(obj, features=genes)
    run_pca(obj, n_pcs=5, features=genes, seed=42)
    find_neighbors(obj, dims=range(5), k_param=10)
    assert (obj.graphs["RNA_snn"]._matrix.diagonal() != 0).sum() == 80

    seen = {}
    real = uu.simplicial_set_embedding

    def spy(*args, **kwargs):
        seen["graph"] = kwargs["graph"]
        return real(*args, **kwargs)

    monkeypatch.setattr(uu, "simplicial_set_embedding", spy)
    run_umap(obj, graph="RNA_snn", seed=42)
    assert (seen["graph"].tocsr().diagonal() != 0).sum() == 0, (
        "self-edges reached simplicial_set_embedding"
    )
