"""Weighted Nearest Neighbor (WNN) multimodal integration.

Port of Seurat's ``FindMultiModalNeighbors()`` (Hao et al., Cell 2021). Given
two or more modalities that already have their own reductions (e.g. RNA
``"pca"`` and ADT ``"apca"``), it learns a per-cell weight for each modality and
builds a joint weighted KNN (``wknn``) and SNN (``wsnn``) graph.

WNN is two coupled stages, and this module mirrors both:

1. ``FindModalityWeights`` (:func:`_modality_weights`) — how much does each cell
   trust each modality? Each cell's embedding is imputed from its *own*
   modality's neighbours and from the *other* modality's neighbours. Whichever
   imputation lands closer wins. The two distances go through an exponential
   kernel with a per-cell bandwidth, their ratio becomes a score, and a softmax
   over modalities turns the scores into weights that sum to 1.

2. ``MultiModalNN`` (:func:`_multi_modal_nn`) — the joint neighbour search. Each
   modality nominates candidate neighbours, the candidates are pooled per cell,
   and every candidate is scored by ``sum_r exp(-d_r / sigma_r) * weight_r``.
   The top ``k_nn`` become that cell's joint neighbours, and ``wknn``/``wsnn``
   are built from *that* ranking — not from a blend of per-modality graphs.

The weights are deliberately decisive rather than hovering near 0.5: the score
is clipped to 200 before the softmax, so ``exp(200)`` against ``exp(0)``
saturates and a cell with a clear preference commits to one modality. That
dynamic range is the point of WNN, and it is what makes the joint graph in
stage 2 behave differently from either modality alone.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import scipy.sparse as sp

from .graph import Graph
from .neighbors import _build_knn, _build_snn

# R's `cross.contant.list` (sic) — guards the score's denominator when the
# cross-modality kernel underflows to 0.
_CROSS_CONSTANT = 1e-4
# R clips the modality score to [0, 200] before the softmax (MinMax in
# FindModalityWeights). exp(200) is ~7e86 — large, but finite in float64.
_SCORE_MAX = 200.0


def find_multi_modal_neighbors(
    seurat,
    reduction_list: Sequence[str] = ("pca", "apca"),
    dims_list: Optional[Sequence[Sequence[int]]] = None,
    k_nn: int = 20,
    l2_norm: bool = True,
    knn_graph_name: str = "wknn",
    snn_graph_name: str = "wsnn",
    knn_range: int = 200,
    prune_snn: float = 1 / 15,
    sd_scale: float = 1.0,
    cross_constant: Optional[float] = None,
    smooth: bool = False,
    seed: int = 42,
) -> None:
    """Compute weighted-nearest-neighbour graphs across modalities.

    Mirrors R's ``FindMultiModalNeighbors(obj, reduction.list =
    list("pca","apca"), dims.list = list(1:30, 1:18))``. Stores:

    * ``seurat.graphs[knn_graph_name]`` — joint KNN graph
    * ``seurat.graphs[snn_graph_name]`` — joint SNN graph (feed to
      ``find_clusters(graph_name=...)`` or ``run_umap(graph=...)``)
    * ``seurat.meta_data["<modality>.weight"]`` — per-cell modality weights

    Parameters
    ----------
    reduction_list : reductions to combine, one per modality
    dims_list      : dims (0-indexed) to use from each reduction; default all
    k_nn           : number of joint neighbours to keep per cell
    l2_norm        : L2-normalise each embedding first (R's ``l2.norm``)
    knn_range      : candidate neighbours each modality nominates before the
                     joint re-ranking (R's ``knn.range``)
    prune_snn      : Jaccard prune threshold for the joint SNN graph
    sd_scale       : scaling on the per-cell kernel bandwidth (R's ``sd.scale``)
    cross_constant : denominator guard in the modality score; default 1e-4
    smooth         : average each cell's modality score over its neighbours

    Notes
    -----
    The weights are stored under each reduction's ``assay_used`` name, so an RNA
    ``"pca"`` and an ADT ``"apca"`` produce ``RNA.weight`` and ``ADT.weight`` —
    the same columns Seurat writes, and readable by the plotting functions as if
    they were features.
    """
    if len(reduction_list) < 2:
        raise ValueError("find_multi_modal_neighbors needs at least 2 reductions.")

    for r in reduction_list:
        if r not in seurat.reductions:
            raise KeyError(
                f"Reduction '{r}' not found. Compute it before WNN "
                "(e.g. run_pca for RNA, run_pca(reduction_name='apca') for ADT)."
            )

    cells = seurat.cell_names()
    n_cells = len(cells)
    cross_constant = _CROSS_CONSTANT if cross_constant is None else cross_constant

    if k_nn >= n_cells:
        raise ValueError(
            f"k_nn ({k_nn}) must be smaller than the number of cells ({n_cells})."
        )

    # Per-modality embeddings (optionally dim-subset, then L2-normalised).
    embs: list[np.ndarray] = []
    for m, r in enumerate(reduction_list):
        e = seurat.reductions[r].cell_embeddings
        if dims_list is not None and dims_list[m] is not None:
            e = e[:, list(dims_list[m])]
        e = np.asarray(e, dtype=float)
        embs.append(_l2_norm(e) if l2_norm else e)

    # Stage 1 — per-cell modality weights, plus the bandwidths and nearest
    # distances stage 2 reuses (R threads these through ModalityWeights@params).
    weights, sigmas, nearest_dist = _modality_weights(
        embs, k_nn, sd_scale, cross_constant, smooth, seed,
    )

    # Stage 2 — joint neighbour search, then graphs off that single ranking.
    select_nn = _multi_modal_nn(
        embs, weights, sigmas, nearest_dist, k_nn, knn_range, seed,
    )

    wknn = _knn_union_graph(select_nn, n_cells)
    wsnn = _build_snn(select_nn, n_cells, select_nn.shape[1], prune_snn).tocsc()

    assay_name = seurat.active_assay
    seurat.graphs[knn_graph_name] = Graph(matrix=wknn, cell_names=cells, assay_used=assay_name)
    seurat.graphs[snn_graph_name] = Graph(matrix=wsnn, cell_names=cells, assay_used=assay_name)

    # Store per-cell modality weights in meta_data, named by each modality's assay.
    for m, r in enumerate(reduction_list):
        assay_used = seurat.reductions[r].assay_used or r
        seurat.meta_data[f"{assay_used}.weight"] = weights[:, m]


# ---------------------------------------------------------------------------
# Stage 1 — modality weights (R: FindModalityWeights)
# ---------------------------------------------------------------------------

def _modality_weights(
    embs: list[np.ndarray],
    k_nn: int,
    sd_scale: float,
    cross_constant: float,
    smooth: bool,
    seed: int,
) -> tuple[np.ndarray, list[np.ndarray], list[np.ndarray]]:
    """Per-cell modality weights, bandwidths and nearest-neighbour distances.

    Returns ``(weights (n × n_mod), sigmas, nearest_dist)``.
    """
    n_mod = len(embs)
    n_cells = embs[0].shape[0]

    # R asks annoy for k.nn neighbours *including* the cell itself, so a cell is
    # imputed from k_nn - 1 others (PredictAssay drops the self column).
    nn_idx: list[np.ndarray] = []
    nearest_dist: list[np.ndarray] = []
    for e in embs:
        idx, dist = _build_knn(e, min(k_nn, n_cells), seed)
        nn_idx.append(idx)
        nearest_dist.append(dist[:, 1])   # distance to the first true neighbour

    # Kernel bandwidth per cell: how far away are the *least* similar cells this
    # cell still shares neighbours with? R builds this SNN with prune = 0 —
    # FindMultiModalNeighbors never passes its own prune.SNN down to it.
    sigmas: list[np.ndarray] = []
    for m in range(n_mod):
        snn_m = _build_snn(nn_idx[m], n_cells, nn_idx[m].shape[1], 0.0)
        sigma = _snn_bandwidth(snn_m, embs[m], k_nn, nearest_dist[m])
        sigmas.append(sigma * sd_scale)

    # Imputation distances: reconstruct each cell from its own modality's
    # neighbours, and from every other modality's neighbours.
    scores: list[np.ndarray] = []
    for m in range(n_mod):
        e = embs[m]
        within = _impute_dist(e, e[nn_idx[m][:, 1:]].mean(axis=1), nearest_dist[m])
        cross = np.column_stack([
            _impute_dist(e, e[nn_idx[mp][:, 1:]].mean(axis=1), nearest_dist[m])
            for mp in range(n_mod) if mp != m
        ])

        # Exponential kernel, then the ratio of own- to cross-modality affinity.
        # Small distance differences become large score differences here — this
        # is what gives the weights their range.
        sigma = np.clip(sigmas[m], 1e-12, None)
        within_k = np.exp(-within / sigma)
        cross_k = np.exp(-cross / sigma[:, None])
        score = within_k[:, None] / (cross_k + cross_constant)
        scores.append(np.clip(score, 0.0, _SCORE_MAX))

    if smooth:
        # R averages each cell's score over its own modality's neighbours. Its
        # indexing collapses the score matrix to a vector, which is exact for
        # the 2-modality case; the row-wise mean here generalises that.
        scores = [
            scores[m][nn_idx[m][:, 1:]].mean(axis=1)
            for m in range(n_mod)
        ]

    # Softmax across modalities.
    numerators = np.column_stack([np.exp(s).sum(axis=1) for s in scores])
    weights = numerators / numerators.sum(axis=1, keepdims=True)
    return weights, sigmas, nearest_dist


def _l2_norm(mat: np.ndarray) -> np.ndarray:
    """Scale each cell's embedding to unit length; zero out any non-finite row."""
    norms = np.linalg.norm(mat, axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = mat / norms
    out[~np.isfinite(out)] = 0.0
    return out


def _impute_dist(x: np.ndarray, y: np.ndarray, nearest_dist: np.ndarray) -> np.ndarray:
    """R's impute_dist: ||x - y|| minus the nearest-neighbour distance, ReLU'd.

    Subtracting d(cell, nearest neighbour) puts every cell on a common footing —
    the kernel then measures excess distance beyond an already-dense or already-
    sparse neighbourhood, rather than raw distance.
    """
    d = np.linalg.norm(x - y, axis=1) - nearest_dist
    return np.maximum(d, 0.0)


def _snn_bandwidth(
    snn: sp.spmatrix,
    emb: np.ndarray,
    k: int,
    nearest_dist: np.ndarray,
) -> np.ndarray:
    """Port of Seurat's SNN_SmallestNonzero_Dist (src/snn.cpp).

    For each cell, take the ``k`` SNN edges with the *smallest* Jaccard weights —
    the most distant cells it still shares any neighbours with — and average
    their embedding distances. Ties at the k-th weight are all included, and
    when ties push the set past ``k``, the ``k`` largest distances win.
    """
    snn = snn.tocsc()
    n = snn.shape[1]
    out = np.zeros(n, dtype=float)
    for i in range(n):
        lo, hi = snn.indptr[i], snn.indptr[i + 1]
        vals, idxs = snn.data[lo:hi], snn.indices[lo:hi]
        if vals.size == 0:
            continue
        order = np.argsort(vals, kind="stable")      # ascending, as in C++
        k_i = min(k, order.size)
        # Everything at or below the k-th smallest weight, ties included.
        threshold = vals[order[k_i - 1]]
        take = int(np.count_nonzero(vals <= threshold))
        d = np.linalg.norm(emb[idxs[order[:take]]] - emb[i], axis=1)
        if nearest_dist[i] > 0:
            d = np.maximum(d - nearest_dist[i], 0.0)
        if d.size > k_i:
            d = np.sort(d)[::-1][:k_i]
        out[i] = d.mean()
    return out


# ---------------------------------------------------------------------------
# Stage 2 — joint neighbour search (R: MultiModalNN)
# ---------------------------------------------------------------------------

def _multi_modal_nn(
    embs: list[np.ndarray],
    weights: np.ndarray,
    sigmas: list[np.ndarray],
    nearest_dist: list[np.ndarray],
    k_nn: int,
    knn_range: int,
    seed: int,
) -> np.ndarray:
    """Return each cell's ``k_nn`` joint neighbours as an (n × k_nn) index array.

    Every modality nominates ``knn_range`` candidates; the candidates are pooled
    per cell and re-scored with the *weighted* kernel affinity summed across
    modalities. A cell whose weight sits almost entirely on one modality ends up
    with that modality's neighbours — which is the intended behaviour, and why
    this cannot be reproduced by blending per-modality graphs after the fact.
    """
    n_mod = len(embs)
    n_cells = embs[0].shape[0]

    # Candidate pool per modality, self dropped (R: Indices(...)[, -1]).
    candidates: list[np.ndarray] = []
    for e in embs:
        idx, _ = _build_knn(e, min(knn_range, n_cells), seed)
        candidates.append(idx[:, 1:])

    select_nn = np.empty((n_cells, k_nn), dtype=int)
    for i in range(n_cells):
        pool = np.unique(np.concatenate([c[i] for c in candidates]))
        affinity = np.zeros(pool.size, dtype=float)
        for m in range(n_mod):
            d = np.linalg.norm(embs[m][pool] - embs[m][i], axis=1) - nearest_dist[m][i]
            np.maximum(d, 0.0, out=d)
            sigma = max(sigmas[m][i], 1e-12)
            affinity += np.exp(-d / sigma) * weights[i, m]
        # Highest joint affinity first; stable so ties fall to the lower index,
        # matching R's order(decreasing = TRUE).
        best = np.argsort(-affinity, kind="stable")[:k_nn]
        select_nn[i] = pool[best]
    return select_nn


def _knn_union_graph(select_nn: np.ndarray, n_cells: int) -> sp.csc_matrix:
    """Binary symmetric KNN graph from the joint neighbours, self included.

    R does ``A + t(A) - t(A) * A`` after setting the diagonal, i.e. a union
    rather than a sum: an edge exists if *either* cell listed the other.
    """
    n, k = select_nn.shape
    rows = np.repeat(np.arange(n), k)
    a = sp.csr_matrix(
        (np.ones(rows.size), (rows, select_nn.flatten())),
        shape=(n_cells, n_cells),
    )
    a = a + sp.eye(n_cells, format="csr")
    a.data[:] = 1.0                       # collapse any duplicate entries
    union = a + a.T - a.T.multiply(a)
    union.data[:] = 1.0
    return union.tocsc()
