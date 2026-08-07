"""Graph-based community detection (clustering).

Mirrors Seurat's FindClusters().
"""
from __future__ import annotations

from collections.abc import Sequence
from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp


def _res_label(r: float) -> str:
    """Format a resolution as R's ``as.character()`` would, for the column name.

    Seurat names the column ``paste0(graph.name, "_res.", resolution)``, so the
    label is R's numeric-to-string conversion and not Python's. The two differ on
    the values users actually type: ``str(1.0)`` is ``"1.0"`` where R gives
    ``"1"``, which would put the partition in ``RNA_snn_res.1.0`` and leave a
    ported script reading ``RNA_snn_res.1`` with a ``KeyError``.

    R's rule (at the default ``scipen = 0``) is to render both fixed and
    scientific and keep whichever is shorter, ties going to fixed — so 0.001
    stays ``0.001`` while 0.0001 becomes ``1e-04``. Verified against R 4.x on 21
    values spanning 1e-7 to 1234567.
    """
    r = float(r)
    fixed = np.format_float_positional(r, trim="-")
    sci = np.format_float_scientific(r, trim="-", exp_digits=2)
    return sci if len(sci) < len(fixed) else fixed


def find_clusters(
    seurat,
    resolution: Union[float, Sequence[float]] = 0.5,
    algorithm: int = 1,
    graph_name: Optional[str] = None,
    random_seed: int = 0,
    n_iterations: int = -1,
    group_singletons: bool = True,
    cluster_name: Optional[Union[str, Sequence[str]]] = None,
) -> None:
    """Apply Louvain or Leiden clustering on the SNN graph.

    Mirrors R's ``FindClusters(pbmc, resolution = 0.5)``, including its
    vector form ``FindClusters(pbmc, resolution = c(0.4, 0.8, 1.2))``.

    Each resolution is written to its own metadata column, named
    ``{graph_name}_res.{resolution}`` as Seurat names it. ``seurat_clusters`` and
    the active identities are set from the **last** resolution in the sequence —
    last as given, not largest, which is what Seurat does.

    Parameters
    ----------
    resolution   : higher values give more / finer clusters. A sequence runs each
                   in turn, as R's ``resolution = c(...)`` does.
    algorithm    : 1 = Louvain, 2 = Louvain (multilevel, igraph's default),
                   4 = Leiden. (3 = SLM is not implemented.)
    graph_name   : SNN graph to use (defaults to '{assay}_snn')
    random_seed  : for reproducibility
    n_iterations : Leiden iterations (-1 = until stable)
    group_singletons : absorb size-1 clusters into their best-connected
                   neighbour, as Seurat's ``GroupSingletons`` does. With
                   ``False`` they are all pooled into one ``"singleton"``
                   cluster instead — again matching Seurat.
    cluster_name : override the generated column name(s); one name, or one per
                   resolution. Seurat's ``cluster.name``. ``seurat_clusters`` is
                   still written either way.

    Notes
    -----
    Seurat runs its own modularity optimiser with ``n.start = 10`` restarts and
    keeps the highest-modularity partition; this runs a single pass of igraph's
    multilevel Louvain. On the same graph that makes truecell's partition land in
    a slightly shallower optimum — measurably so, but not necessarily a worse
    one. See the clustering section of ``tutorials/integration_vignette.md``.

    Each resolution is clustered from the same seed rather than from a running
    RNG stream, so a partition does not depend on which resolutions preceded it
    or on the order they were given in. Verified against Seurat 5.5.1, where
    resolution 0.8 gives the same partition alone, in ``c(0.4, 0.8, 1.2)``, and
    in ``c(1.2, 0.8, 0.4)``.
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

    # Validated once rather than per iteration. This does *not* prevent a partial
    # write — the dispatch sits at the top of the loop body, so a bad `algorithm`
    # would raise on the first resolution either way, before anything is stored.
    # It is here so the check does not depend on the loop at all.
    #
    # A partial write is still reachable: if a *later* resolution fails inside
    # igraph, the earlier columns are already on the object. That is Seurat's
    # behaviour too, and is left alone.
    if algorithm == 3:
        raise NotImplementedError(
            "algorithm=3 (SLM) is not implemented. Use 1 or 2 (Louvain) or "
            "4 (Leiden)."
        )
    if algorithm not in (1, 2, 4):
        raise ValueError(
            f"Unknown algorithm {algorithm!r}. Use 1 or 2 (Louvain) or 4 (Leiden)."
        )

    # `np.number` is here for the numpy scalars a caller gets out of an array;
    # np.float64 subclasses float but np.float32 does not, and iterating one
    # raises rather than falling through to the sequence branch.
    if isinstance(resolution, (int, float, np.number)):
        resolutions = [float(resolution)]
    else:
        resolutions = [float(r) for r in resolution]
    if not resolutions:
        raise ValueError("`resolution` is empty; give at least one value.")

    if cluster_name is None:
        names = [f"{graph_name}_res.{_res_label(r)}" for r in resolutions]
    else:
        names = [cluster_name] if isinstance(cluster_name, str) else list(cluster_name)
        if len(names) != len(resolutions):
            raise ValueError(
                f"`cluster_name` has {len(names)} name(s) for "
                f"{len(resolutions)} resolution(s); give one per resolution."
            )

    cluster_series = None
    for res, name in zip(resolutions, names):
        if algorithm == 4:
            labels = _leiden_clustering(mat, res, random_seed, n_iterations)
        else:
            # python-igraph's community_multilevel is the multilevel Louvain
            # algorithm (closest to Seurat's algorithm 1/2).
            labels = _louvain_clustering(mat, res, random_seed)

        str_labels = _group_singletons(
            np.asarray([str(c) for c in labels]), mat, group_singletons
        )
        present = sorted(set(str_labels),
                         key=lambda s: (not s.isdigit(), s.isdigit() and int(s), s))
        cluster_series = pd.Categorical(str_labels, categories=present)
        seurat.meta_data[name] = cluster_series

    # The last resolution given, not the largest — Seurat takes the last column
    # of its results frame, so `resolution = c(1.2, 0.8, 0.4)` leaves the object
    # sitting on 0.4.
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
