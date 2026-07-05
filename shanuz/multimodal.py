"""Weighted Nearest Neighbor (WNN) multimodal integration.

Mirrors Seurat's FindMultiModalNeighbors() (Hao et al., Cell 2021). Given two
or more modalities that already have their own reductions (e.g. RNA ``"pca"``
and ADT ``"apca"``), it learns a per-cell weight for each modality and builds a
joint weighted KNN (``wknn``) and SNN (``wsnn``) graph.

The per-cell weights follow the roadmap's sanctioned approximation of the
Seurat kernel: a modality is trusted more for a cell when that cell's *own*
modality neighbours reconstruct its embedding better than the *other*
modality's neighbours do — a scale-invariant ratio within each modality's own
space, which sidesteps the incomparable-units problem across modalities.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import scipy.sparse as sp

from .graph import Graph
from .neighbors import _build_knn, _build_snn


def find_multi_modal_neighbors(
    seurat,
    reduction_list: Sequence[str] = ("pca", "apca"),
    dims_list: Optional[Sequence[Sequence[int]]] = None,
    k_nn: int = 20,
    knn_graph_name: str = "wknn",
    snn_graph_name: str = "wsnn",
    prune_snn: float = 1 / 15,
    seed: int = 42,
) -> None:
    """Compute weighted-nearest-neighbour graphs across modalities.

    Mirrors R's ``FindMultiModalNeighbors(obj, reduction.list =
    list("pca","apca"), dims.list = list(1:30, 1:18))``. Stores:

    * ``seurat.graphs[knn_graph_name]`` — weighted KNN graph
    * ``seurat.graphs[snn_graph_name]`` — weighted SNN graph (feed to
      ``find_clusters(graph_name=...)`` or ``run_umap(graph=...)``)
    * ``seurat.meta_data["<modality>.weight"]`` — per-cell modality weights

    Parameters
    ----------
    reduction_list : reductions to combine, one per modality
    dims_list      : dims (0-indexed) to use from each reduction; default all
    k_nn           : number of neighbours per modality
    prune_snn      : Jaccard prune threshold for the SNN graph
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
    n_mod = len(reduction_list)

    # Per-modality embeddings (optionally dim-subset).
    embs: list[np.ndarray] = []
    for m, r in enumerate(reduction_list):
        e = seurat.reductions[r].cell_embeddings
        if dims_list is not None and dims_list[m] is not None:
            e = e[:, list(dims_list[m])]
        embs.append(np.asarray(e, dtype=float))

    # Per-modality KNN. Request k_nn+1 so we can drop self (column 0) and keep
    # k_nn genuine neighbours for reconstruction; the +self set feeds the graph.
    nn_idx_self: list[np.ndarray] = []   # (n × k_nn+1), includes self at col 0
    nn_idx_other: list[np.ndarray] = []  # (n × k_nn),   excludes self
    for e in embs:
        idx, _ = _build_knn(e, min(k_nn + 1, n_cells), seed)
        nn_idx_self.append(idx)
        nn_idx_other.append(idx[:, 1:])

    # Reconstruction distances per modality.
    #   d_same  : cell reconstructed from its OWN modality's neighbours
    #   d_cross : cell reconstructed from the OTHER modality's neighbours
    #             (best / most competitive alternative when >2 modalities)
    eps = 1e-8
    thetas = np.zeros((n_cells, n_mod), dtype=float)
    for m in range(n_mod):
        e = embs[m]
        pred_same = e[nn_idx_other[m]].mean(axis=1)          # (n × dims)
        d_same = np.linalg.norm(e - pred_same, axis=1)

        d_cross = np.full(n_cells, np.inf)
        for mp in range(n_mod):
            if mp == m:
                continue
            pred_cross = e[nn_idx_other[mp]].mean(axis=1)
            d = np.linalg.norm(e - pred_cross, axis=1)
            d_cross = np.minimum(d_cross, d)

        # Modality affinity in (0,1): high when own neighbours reconstruct the
        # cell better (d_same small) than the other modality's neighbours.
        thetas[:, m] = d_cross / (d_same + d_cross + eps)

    # Per-cell weights across modalities (sum to 1).
    weights = thetas / np.clip(thetas.sum(axis=1, keepdims=True), eps, None)

    # Per-modality SNN graphs, weighted per cell and combined.
    combined = sp.csr_matrix((n_cells, n_cells), dtype=np.float64)
    knn_combined = sp.csr_matrix((n_cells, n_cells), dtype=np.float64)
    for m in range(n_mod):
        snn_m = _build_snn(nn_idx_self[m], n_cells, nn_idx_self[m].shape[1], prune_snn)
        knn_m = _knn_adjacency(nn_idx_self[m], n_cells)
        w = sp.diags(weights[:, m])          # row i scaled by this cell's weight
        combined = combined + (w @ snn_m.tocsr())
        knn_combined = knn_combined + (w @ knn_m.tocsr())

    # Symmetrize (row weighting breaks symmetry) → symmetric weighted graphs.
    wsnn = ((combined + combined.T) * 0.5).tocsc()
    wknn = ((knn_combined + knn_combined.T) * 0.5).tocsc()

    assay_name = seurat.active_assay
    seurat.graphs[knn_graph_name] = Graph(matrix=wknn, cell_names=cells, assay_used=assay_name)
    seurat.graphs[snn_graph_name] = Graph(matrix=wsnn, cell_names=cells, assay_used=assay_name)

    # Store per-cell modality weights in meta_data, named by each modality's assay.
    for m, r in enumerate(reduction_list):
        assay_used = seurat.reductions[r].assay_used or r
        seurat.meta_data[f"{assay_used}.weight"] = weights[:, m]


def _knn_adjacency(nn_idx: np.ndarray, n_cells: int) -> sp.csr_matrix:
    """Binary KNN adjacency (nn_idx includes self); not symmetrised here."""
    n, k = nn_idx.shape
    rows = np.repeat(np.arange(n), k)
    cols = nn_idx.flatten()
    data = np.ones(len(rows), dtype=np.float64)
    return sp.csr_matrix((data, (rows, cols)), shape=(n_cells, n_cells))
