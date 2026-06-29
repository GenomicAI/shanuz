"""Graph-based community detection (clustering).

Mirrors Seurat's FindClusters().
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp


def find_clusters(
    seurat,
    resolution: float = 0.5,
    algorithm: int = 1,
    graph_name: Optional[str] = None,
    random_seed: int = 0,
    n_iterations: int = -1,
) -> None:
    """Apply Louvain or Leiden clustering on the SNN graph.

    Mirrors R's FindClusters(pbmc, resolution = 0.5).
    Stores cluster assignments in seurat.meta_data['seurat_clusters']
    and updates seurat.idents.

    Parameters
    ----------
    resolution   : higher values give more / finer clusters
    algorithm    : 1 = Louvain (default), 2 = Louvain refined,
                   4 = Leiden
    graph_name   : SNN graph to use (defaults to '{assay}_snn')
    random_seed  : for reproducibility
    n_iterations : Leiden iterations (-1 = until stable)
    """
    assay_name = seurat.active_assay
    if graph_name is None:
        graph_name = f"{assay_name}_snn"
        if graph_name not in seurat.graphs:
            # Try knn graph
            graph_name = f"{assay_name}_nn"

    if graph_name not in seurat.graphs:
        raise KeyError(
            f"Graph '{graph_name}' not found. Run find_neighbors() first."
        )

    graph = seurat.graphs[graph_name]
    mat = graph._matrix  # scipy sparse (cells × cells)

    if algorithm == 4:
        labels = _leiden_clustering(mat, resolution, random_seed, n_iterations)
    else:
        labels = _louvain_clustering(mat, resolution, random_seed)

    cluster_series = pd.Categorical(
        [str(c) for c in labels],
        categories=[str(i) for i in sorted(set(labels))],
    )

    seurat.meta_data["seurat_clusters"] = cluster_series
    seurat.idents = cluster_series


# ------------------------------------------------------------------
# Louvain via igraph
# ------------------------------------------------------------------

def _louvain_clustering(
    mat: sp.spmatrix,
    resolution: float,
    seed: int,
) -> np.ndarray:
    """Louvain community detection using python-igraph."""
    import igraph as ig

    g = _sparse_to_igraph(mat)
    np.random.seed(seed)

    try:
        result = g.community_multilevel(
            weights="weight",
            resolution=resolution,
            return_levels=False,
        )
    except Exception:
        result = g.community_multilevel(resolution=resolution, return_levels=False)

    labels = np.array(result.membership)

    # Re-number clusters by size (largest = 0) — mirrors Seurat behavior
    labels = _renumber_by_size(labels)
    return labels


# ------------------------------------------------------------------
# Leiden via leidenalg
# ------------------------------------------------------------------

def _leiden_clustering(
    mat: sp.spmatrix,
    resolution: float,
    seed: int,
    n_iterations: int,
) -> np.ndarray:
    """Leiden community detection."""
    import igraph as ig
    import leidenalg

    g = _sparse_to_igraph(mat)

    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=resolution,
        weights="weight" if g.is_weighted() else None,
        n_iterations=n_iterations,
        seed=seed,
    )

    labels = np.array(partition.membership)
    labels = _renumber_by_size(labels)
    return labels


# ------------------------------------------------------------------
# Helper: sparse matrix → igraph Graph
# ------------------------------------------------------------------

def _sparse_to_igraph(mat: sp.spmatrix):
    """Convert scipy sparse adjacency matrix to an igraph Graph."""
    import igraph as ig

    mat = mat.tocoo()
    n = mat.shape[0]

    # Only upper triangle (undirected)
    mask = mat.row < mat.col
    rows = mat.row[mask].tolist()
    cols = mat.col[mask].tolist()
    weights = mat.data[mask].tolist()

    edges = list(zip(rows, cols))
    g = ig.Graph(n=n, edges=edges, directed=False)
    g.es["weight"] = weights
    return g


def _renumber_by_size(labels: np.ndarray) -> np.ndarray:
    """Renumber cluster labels so 0 is the largest cluster."""
    unique, counts = np.unique(labels, return_counts=True)
    order = unique[np.argsort(-counts)]
    mapping = {old: new for new, old in enumerate(order)}
    return np.array([mapping[l] for l in labels])
