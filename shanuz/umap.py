"""UMAP dimensionality reduction.

Mirrors Seurat's RunUMAP().
"""
from __future__ import annotations

from typing import Optional, Union

import numpy as np

from .dimreduc import DimReduc


def run_umap(
    seurat,
    dims: Optional[Union[list[int], range]] = None,
    reduction: str = "pca",
    n_components: int = 2,
    n_neighbors: int = 30,
    min_dist: float = 0.3,
    metric: str = "euclidean",
    reduction_name: str = "umap",
    reduction_key: str = "UMAP_",
    seed: int = 42,
    assay: Optional[str] = None,
) -> None:
    """Compute UMAP embedding.

    Mirrors R's RunUMAP(pbmc, dims = 1:10).
    Stores a DimReduc in seurat.reductions[reduction_name].

    Parameters
    ----------
    dims           : which dimensions of 'reduction' to use (0-indexed)
    reduction      : source reduction ('pca' by default)
    n_components   : output dimensions (2 for visualization)
    n_neighbors    : UMAP n_neighbors (Seurat default 30)
    min_dist       : UMAP min_dist (Seurat default 0.3)
    metric         : distance metric
    reduction_name : storage key in seurat.reductions
    seed           : random seed
    """
    from umap import UMAP

    assay_name = assay or seurat.active_assay

    if reduction not in seurat.reductions:
        raise KeyError(f"Reduction '{reduction}' not found. Run run_pca() first.")

    dr = seurat.reductions[reduction]
    embeddings = dr.cell_embeddings  # (cells × n_dims)

    if dims is None:
        emb = embeddings
    else:
        dims_list = list(dims)
        emb = embeddings[:, dims_list]

    reducer = UMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=seed,
    )
    umap_coords = reducer.fit_transform(emb)  # (cells × n_components)

    cells = seurat.cell_names()
    dim_names = [f"{reduction_key}{i + 1}" for i in range(n_components)]

    dr_umap = DimReduc(
        cell_embeddings=umap_coords,
        assay_used=assay_name,
        key=reduction_key,
        cell_names=cells,
        feature_names=dim_names,
    )

    seurat.reductions[reduction_name] = dr_umap
