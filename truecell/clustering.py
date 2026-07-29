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
    group_singletons: bool = True,
) -> None:
    """Apply Louvain or Leiden clustering on the SNN graph.

    Mirrors R's FindClusters(pbmc, resolution = 0.5).
    Stores cluster assignments in seurat.meta_data['seurat_clusters']
    and updates seurat.idents.

    Parameters
    ----------
    resolution   : higher values give more / finer clusters
    algorithm    : 1 = Louvain, 2 = Louvain (multilevel, igraph's default),
                   4 = Leiden. (3 = SLM is not implemented.)
    graph_name   : SNN graph to use (defaults to '{assay}_snn')
    random_seed  : for reproducibility
    n_iterations : Leiden iterations (-1 = until stable)
    group_singletons : absorb size-1 clusters into their best-connected
                   neighbour, as Seurat's ``GroupSingletons`` does. With
                   ``False`` they are all pooled into one ``"singleton"``
                   cluster instead — again matching Seurat.

    Notes
    -----
    Seurat runs its own modularity optimiser with ``n.start = 10`` restarts and
    keeps the highest-modularity partition; this runs a single pass of igraph's
    multilevel Louvain. On the same graph that makes truecell's partition land in
    a slightly shallower optimum — measurably so, but not necessarily a worse
    one. See the clustering section of ``tutorials/integration_vignette.md``.
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
    elif algorithm in (1, 2):
        # python-igraph's community_multilevel is the multilevel Louvain
        # algorithm (closest to Seurat's algorithm 1/2).
        labels = _louvain_clustering(mat, resolution, random_seed)
    elif algorithm == 3:
        raise NotImplementedError(
            "algorithm=3 (SLM) is not implemented. Use 1 or 2 (Louvain) or "
            "4 (Leiden)."
        )
    else:
        raise ValueError(
            f"Unknown algorithm {algorithm!r}. Use 1 or 2 (Louvain) or 4 (Leiden)."
        )

    str_labels = _group_singletons(
        np.asarray([str(c) for c in labels]), mat, group_singletons
    )
    present = sorted(set(str_labels), key=lambda s: (not s.isdigit(), s.isdigit() and int(s), s))
    cluster_series = pd.Categorical(str_labels, categories=present)

    seurat.meta_data["seurat_clusters"] = cluster_series
    seurat.idents = cluster_series


# ------------------------------------------------------------------
# Singleton absorption — Seurat's GroupSingletons
# ------------------------------------------------------------------

def _group_singletons(
    ids: np.ndarray,
    snn: sp.spmatrix,
    group_singletons: bool,
) -> np.ndarray:
    """Absorb size-1 clusters, mirroring Seurat's ``GroupSingletons``.

    Every cluster holding exactly one cell is reassigned to whichever
    *non*-singleton cluster it is most connected to, scored as the mean SNN
    weight from that cell to every cell of the candidate cluster (Seurat:
    ``sum(subSNN) / (nrow * ncol)``). With ``group_singletons = False`` they are
    instead pooled into a single ``"singleton"`` cluster.

    Seurat breaks ties by ``sample()`` under ``set.seed(1)``; an R RNG draw is
    not reproducible from Python, so the lowest-numbered cluster wins here.
    Ties need two candidates to share a mean connectivity exactly, which the
    Jaccard weights make rare.
    """
    ids = ids.copy()
    values, counts = np.unique(ids, return_counts=True)
    singletons = set(values[counts == 1].tolist())
    if not singletons:
        return ids
    if not group_singletons:
        # `ids` is a fixed-width unicode array sized to the labels it already
        # holds, so a plain assignment would truncate "singleton" to "si".
        ids = ids.astype(object)
        ids[np.isin(ids, list(singletons))] = "singleton"
        return ids.astype(str)

    # Candidate targets are fixed before the loop, as in Seurat: a cluster that
    # has just absorbed a singleton does not itself become a new target.
    targets = [v for v in values if v not in singletons]
    if not targets:
        return ids
    snn = snn.tocsr()
    members = {t: np.flatnonzero(ids == t) for t in targets}

    # Seurat iterates singletons in order of first appearance in `ids`.
    order = sorted(singletons, key=lambda s: int(np.flatnonzero(ids == s)[0]))
    for s in order:
        cell = int(np.flatnonzero(ids == s)[0])
        row = snn[cell].toarray().ravel()
        best = max(targets, key=lambda t: (row[members[t]].mean(), -_key(t)))
        ids[cell] = best
    return ids


def _key(label: str) -> float:
    """Sort key that keeps numeric cluster labels in numeric order."""
    try:
        return float(label)
    except ValueError:
        return float("inf")


# ------------------------------------------------------------------
# Louvain via igraph
# ------------------------------------------------------------------

def _seed_igraph(seed: int) -> None:
    """Seed igraph's own RNG (igraph does not use numpy's global RNG)."""
    import random as _random
    import igraph as ig

    _random.seed(seed)
    try:
        ig.set_random_number_generator(_random)
    except Exception:
        # Older igraph: fall back to seeding the stdlib RNG igraph reads from.
        pass


def _louvain_clustering(
    mat: sp.spmatrix,
    resolution: float,
    seed: int,
) -> np.ndarray:
    """Louvain community detection using python-igraph."""
    import warnings

    g = _sparse_to_igraph(mat)
    _seed_igraph(seed)

    try:
        result = g.community_multilevel(
            weights="weight",
            resolution=resolution,
            return_levels=False,
        )
    except TypeError:
        # Older igraph builds lack the `resolution` parameter. Keep the edge
        # weights — never silently downgrade to an unweighted clustering.
        warnings.warn(
            "This python-igraph version does not support the 'resolution' "
            "parameter for community_multilevel; clustering at the default "
            "resolution (1.0).",
            RuntimeWarning,
        )
        result = g.community_multilevel(weights="weight", return_levels=False)

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
    import leidenalg

    g = _sparse_to_igraph(mat)
    _seed_igraph(seed)

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
    return np.array([mapping[lab] for lab in labels])
