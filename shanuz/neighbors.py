"""KNN and SNN graph construction.

Mirrors Seurat's FindNeighbors().
"""
from __future__ import annotations

from typing import Optional, Union

import numpy as np
import scipy.sparse as sp

from .graph import Graph


def find_neighbors(
    seurat,
    dims: Optional[Union[list[int], range]] = None,
    k_param: int = 20,
    assay: Optional[str] = None,
    reduction: str = "pca",
    graph_name: Optional[str] = None,
    nn_name: Optional[str] = None,
    prune_snn: float = 1 / 15,
    seed: int = 42,
) -> None:
    """Build KNN and SNN graphs from a low-dimensional embedding.

    Mirrors R's FindNeighbors(pbmc, dims = 1:10).
    Stores Graph objects in seurat.graphs[graph_name + '_nn'] and
    seurat.graphs[graph_name + '_snn'].

    Parameters
    ----------
    dims        : which PCs to use (0-indexed; default all available)
    k_param     : number of nearest neighbors
    reduction   : which reduction to use ('pca' by default)
    graph_name  : prefix for graph names (defaults to active assay name)
    prune_snn   : edges with Jaccard index below this are pruned (Seurat default 1/15)
    """
    assay_name = assay or seurat.active_assay

    # Get embeddings
    if reduction not in seurat.reductions:
        raise KeyError(f"Reduction '{reduction}' not found. Run run_pca() first.")
    dr = seurat.reductions[reduction]
    embeddings = dr.cell_embeddings  # (cells × dims)

    if dims is None:
        emb = embeddings
    else:
        dims_list = list(dims)
        emb = embeddings[:, dims_list]

    cells = seurat.cell_names()
    n_cells = len(cells)

    # Build KNN
    nn_idx, nn_dist = _build_knn(emb, k_param, seed)

    # Build KNN sparse graph (symmetric)
    knn_mat = _knn_to_sparse(nn_idx, n_cells)

    # Build SNN (shared nearest neighbor) sparse graph with Jaccard weights
    snn_mat = _build_snn(nn_idx, n_cells, k_param, prune_snn)

    prefix = graph_name or assay_name
    knn_name = f"{prefix}_nn"
    snn_name = f"{prefix}_snn"

    seurat.graphs[knn_name] = Graph(
        matrix=knn_mat, cell_names=cells, assay_used=assay_name
    )
    seurat.graphs[snn_name] = Graph(
        matrix=snn_mat, cell_names=cells, assay_used=assay_name
    )


# ------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------

def _build_knn(
    embeddings: np.ndarray,
    k: int,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (nn_idx, nn_dist) arrays of shape (n_cells, k)."""
    from sklearn.neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean", n_jobs=-1)
    nn.fit(embeddings)
    distances, indices = nn.kneighbors(embeddings)
    # Drop self (first column is always distance=0 to itself)
    return indices[:, 1:], distances[:, 1:]


def _knn_to_sparse(nn_idx: np.ndarray, n_cells: int) -> sp.csc_matrix:
    """Build a symmetric binary KNN adjacency matrix."""
    n, k = nn_idx.shape
    rows = np.repeat(np.arange(n), k)
    cols = nn_idx.flatten()
    data = np.ones(len(rows), dtype=np.float64)

    mat = sp.coo_matrix((data, (rows, cols)), shape=(n_cells, n_cells))
    # Symmetrize
    mat = mat + mat.T
    mat.data = np.ones_like(mat.data)
    return mat.tocsc()


def _build_snn(
    nn_idx: np.ndarray,
    n_cells: int,
    k: int,
    prune_snn: float,
) -> sp.csc_matrix:
    """Build sparse SNN graph with Jaccard similarity weights.

    SNN[i,j] = |NN(i) ∩ NN(j)| / |NN(i) ∪ NN(j)|
    Edges below prune_snn are dropped.
    """
    # Build set-based NN lookup for Jaccard computation
    # This is O(n*k) memory-efficient via sparse matrix multiplication

    n = nn_idx.shape[0]
    # Build boolean membership matrix: M[i,j] = 1 if j is a NN of i
    rows = np.repeat(np.arange(n), k)
    cols = nn_idx.flatten()
    vals = np.ones(len(rows), dtype=np.float32)
    # Include self
    self_rows = np.arange(n)
    self_cols = np.arange(n)
    self_vals = np.ones(n, dtype=np.float32)

    rows_all = np.concatenate([rows, self_rows])
    cols_all = np.concatenate([cols, self_cols])
    vals_all = np.concatenate([vals, self_vals])

    M = sp.csr_matrix(
        (vals_all, (rows_all, cols_all)), shape=(n_cells, n_cells), dtype=np.float32
    )

    # |NN(i) ∩ NN(j)| = (M @ M.T)[i,j]
    intersection = (M @ M.T).toarray().astype(float)

    # |NN(i)| = k+1 (including self)
    union = (k + 1) + (k + 1) - intersection

    # Jaccard
    jaccard = intersection / np.maximum(union, 1)

    # Prune and convert to sparse
    jaccard[jaccard < prune_snn] = 0
    np.fill_diagonal(jaccard, 0)

    return sp.csc_matrix(jaccard)
