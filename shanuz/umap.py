"""UMAP dimensionality reduction.

Mirrors Seurat's RunUMAP().
"""
from __future__ import annotations

from typing import Optional, Union

import numpy as np
import scipy.sparse as sp

from .dimreduc import DimReduc


def run_umap(
    seurat,
    dims: Optional[Union[list[int], range]] = None,
    reduction: str = "pca",
    graph: Optional[str] = None,
    n_components: int = 2,
    n_neighbors: int = 30,
    min_dist: float = 0.3,
    metric: str = "euclidean",
    reduction_name: str = "umap",
    reduction_key: str = "UMAP_",
    seed: int = 42,
    assay: Optional[str] = None,
) -> None:
    """Compute a UMAP embedding.

    Mirrors R's RunUMAP(pbmc, dims = 1:10) and RunUMAP(pbmc, graph = "wsnn").
    Stores a DimReduc in seurat.reductions[reduction_name].

    Two input modes (mutually exclusive):

    * ``reduction`` (default): embed the cells from a low-dimensional
      reduction (PCA/Harmony/ICA…). The fitted ``umap-learn`` model is stashed
      in ``dr.misc["umap_model"]`` for later transform-only projection.
    * ``graph``: embed a precomputed neighbour graph directly (e.g. the WNN
      ``"wsnn"`` graph from ``find_multi_modal_neighbors``), via UMAP's
      ``simplicial_set_embedding``. ``reduction``/``dims`` are ignored.

    Parameters
    ----------
    dims           : which dimensions of 'reduction' to use (0-indexed)
    reduction      : source reduction ('pca' by default)
    graph          : name of a precomputed graph in seurat.graphs to embed
                     (takes precedence over ``reduction`` when given)
    n_components   : output dimensions (2 for visualization)
    n_neighbors    : UMAP n_neighbors (Seurat default 30)
    min_dist       : UMAP min_dist (Seurat default 0.3)
    metric         : distance metric (reduction mode only)
    reduction_name : storage key in seurat.reductions
    seed           : random seed
    """
    assay_name = assay or seurat.active_assay
    cells = seurat.cell_names()
    dim_names = [f"{reduction_key}{i + 1}" for i in range(n_components)]

    if graph is not None:
        coords = _umap_from_graph(seurat, graph, n_components, min_dist, seed)
        seurat.reductions[reduction_name] = DimReduc(
            cell_embeddings=coords,
            assay_used=assay_name,
            key=reduction_key,
            cell_names=cells,
            feature_names=dim_names,
            misc={"umap_graph": graph},
        )
        return

    from umap import UMAP

    if reduction not in seurat.reductions:
        raise KeyError(f"Reduction '{reduction}' not found. Run run_pca() first.")

    embeddings = seurat.reductions[reduction].cell_embeddings  # (cells × n_dims)
    if dims is None:
        emb = embeddings
    else:
        emb = embeddings[:, list(dims)]

    reducer = UMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=seed,
    )
    umap_coords = reducer.fit_transform(emb)  # (cells × n_components)

    seurat.reductions[reduction_name] = DimReduc(
        cell_embeddings=umap_coords,
        assay_used=assay_name,
        key=reduction_key,
        cell_names=cells,
        feature_names=dim_names,
        misc={"umap_model": reducer},
    )


def _umap_from_graph(
    seurat,
    graph: str,
    n_components: int,
    min_dist: float,
    seed: int,
) -> np.ndarray:
    """Embed a precomputed affinity graph with UMAP's simplicial_set_embedding."""
    from sklearn.utils import check_random_state
    from umap.umap_ import find_ab_params, simplicial_set_embedding

    if graph not in seurat.graphs:
        raise KeyError(
            f"Graph '{graph}' not found. Run find_neighbors() or "
            "find_multi_modal_neighbors() first."
        )
    g = seurat.graphs[graph]
    mat = g.tocsr() if hasattr(g, "tocsr") else sp.csr_matrix(g)
    # UMAP expects a symmetric fuzzy-simplicial-set (affinity) graph.
    mat = mat.maximum(mat.T)

    a, b = find_ab_params(1.0, min_dist)
    coords, _ = simplicial_set_embedding(
        data=None,
        graph=mat.tocoo(),
        n_components=n_components,
        initial_alpha=1.0,
        a=a,
        b=b,
        gamma=1.0,
        negative_sample_rate=5,
        n_epochs=200,
        init="spectral",
        random_state=check_random_state(seed),
        metric="euclidean",
        metric_kwds={},
        densmap=False,
        densmap_kwds={},
        output_dens=False,
        output_metric="euclidean",
        output_metric_kwds={},
    )
    return np.asarray(coords)
