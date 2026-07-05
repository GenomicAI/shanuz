"""Dimensionality reduction: PCA and related helpers.

Mirrors Seurat's RunPCA() / RunSVD().
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import scipy.sparse as sp

from .dimreduc import DimReduc


def run_pca(
    seurat,
    n_pcs: int = 50,
    features: Optional[list[str]] = None,
    assay: Optional[str] = None,
    reduction_name: str = "pca",
    reduction_key: str = "PC_",
    seed: int = 42,
    layer: str = "scale.data",
) -> None:
    """Compute PCA on scaled data.

    Mirrors R's RunPCA(pbmc, features = VariableFeatures(object = pbmc)).
    Stores a DimReduc in seurat.reductions[reduction_name].

    Parameters
    ----------
    n_pcs          : number of principal components
    features       : genes to use (defaults to variable features)
    assay          : assay name (defaults to active assay)
    reduction_name : key for storage in seurat.reductions
    reduction_key  : prefix for dimension names (e.g. 'PC_')
    seed           : random seed for reproducibility
    layer          : which layer to take data from
    """
    from sklearn.decomposition import PCA, TruncatedSVD

    assay_name = assay or seurat.active_assay
    assay_obj = seurat.assays[assay_name]

    from .assay5 import Assay5
    from .assay import Assay

    # Determine feature set
    if features is None:
        if isinstance(assay_obj, Assay5):
            features = assay_obj.variable_features
            if not features:
                features = assay_obj._all_feature_names
        else:
            features = assay_obj.var_features or assay_obj._feature_names

    # Get scale.data for the selected features
    scaled = _get_scaled_data(assay_obj, features, layer)
    # scaled shape: (n_features_selected × n_cells)
    # PCA expects (n_samples × n_features) so we transpose: (n_cells × n_features)
    data_t = scaled.T  # (cells × features)

    n_pcs = min(n_pcs, min(data_t.shape) - 1)

    np.random.seed(seed)
    if sp.issparse(data_t):
        # TruncatedSVD for sparse (already centered is assumed)
        svd = TruncatedSVD(n_components=n_pcs, random_state=seed)
        embeddings = svd.fit_transform(data_t)
        loadings = svd.components_.T  # (features × n_pcs)
        explained = svd.explained_variance_ratio_
    else:
        pca = PCA(n_components=n_pcs, random_state=seed)
        embeddings = pca.fit_transform(data_t)  # (cells × n_pcs)
        loadings = pca.components_.T  # (features × n_pcs)
        explained = pca.explained_variance_ratio_

    # Per-PC standard deviation (sample SD, ddof=1) — matches Seurat's @stdev.
    stdev = np.sqrt(np.var(embeddings, axis=0, ddof=1))

    # Cell names
    cells = seurat.cell_names()

    dr = DimReduc(
        cell_embeddings=embeddings,
        feature_loadings=loadings,
        assay_used=assay_name,
        stdev=stdev,
        key=reduction_key,
        cell_names=cells,
        feature_names=list(features),
    )

    seurat.reductions[reduction_name] = dr


def run_ica(
    seurat,
    nics: int = 50,
    features: Optional[list[str]] = None,
    assay: Optional[str] = None,
    reduction_name: str = "ica",
    reduction_key: str = "ICA_",
    seed: int = 42,
    layer: str = "scale.data",
    max_iter: int = 200,
) -> None:
    """Independent Component Analysis on scaled data.

    Mirrors R's ``RunICA(obj, nics = 50)``. Stores a DimReduc (embeddings +
    loadings) under ``reduction_name``; ``find_neighbors`` / ``run_umap``
    already accept ``reduction="ica"``.
    """
    from sklearn.decomposition import FastICA

    assay_name = assay or seurat.active_assay
    assay_obj = seurat.assays[assay_name]

    from .assay5 import Assay5

    if features is None:
        if isinstance(assay_obj, Assay5):
            features = assay_obj.variable_features
            if not features:
                features = assay_obj._all_feature_names
        else:
            features = assay_obj.var_features or assay_obj._feature_names

    scaled = _get_scaled_data(assay_obj, features, layer)  # (features × cells)
    data_t = scaled.T  # (cells × features)
    if sp.issparse(data_t):
        data_t = data_t.toarray()

    nics = min(nics, min(data_t.shape))

    ica = FastICA(n_components=nics, random_state=seed, max_iter=max_iter)
    embeddings = ica.fit_transform(data_t)  # (cells × nics)
    loadings = ica.components_.T  # (features × nics)

    cells = seurat.cell_names()
    seurat.reductions[reduction_name] = DimReduc(
        cell_embeddings=embeddings,
        feature_loadings=loadings,
        assay_used=assay_name,
        key=reduction_key,
        cell_names=cells,
        feature_names=list(features),
    )


def run_tsne(
    seurat,
    dims: Optional[list[int]] = None,
    reduction: str = "pca",
    n_components: int = 2,
    perplexity: float = 30.0,
    reduction_name: str = "tsne",
    reduction_key: str = "tSNE_",
    seed: int = 42,
    assay: Optional[str] = None,
) -> None:
    """t-SNE embedding from an existing reduction.

    Mirrors R's ``RunTSNE(obj, dims = 1:10)``. Stores a DimReduc under
    ``reduction_name``.
    """
    from sklearn.manifold import TSNE

    assay_name = assay or seurat.active_assay

    if reduction not in seurat.reductions:
        raise KeyError(f"Reduction '{reduction}' not found. Run run_pca() first.")
    emb = seurat.reductions[reduction].cell_embeddings
    if dims is not None:
        emb = emb[:, list(dims)]

    tsne = TSNE(
        n_components=n_components,
        perplexity=perplexity,
        random_state=seed,
        init="pca",
    )
    coords = tsne.fit_transform(emb)

    cells = seurat.cell_names()
    dim_names = [f"{reduction_key}{i + 1}" for i in range(n_components)]
    seurat.reductions[reduction_name] = DimReduc(
        cell_embeddings=coords,
        assay_used=assay_name,
        key=reduction_key,
        cell_names=cells,
        feature_names=dim_names,
    )


def _get_scaled_data(assay_obj, features: list[str], layer: str) -> np.ndarray:
    """Extract scaled data for the given features as a (features × cells) array."""
    from .assay5 import Assay5
    from .assay import Assay
    from ._sparse import is_matrix_empty

    if isinstance(assay_obj, Assay5):
        if layer in assay_obj.layers:
            mat = assay_obj.layers[layer]
            # scale.data may be stored for a subset of features
            scaled_features = getattr(assay_obj, "_scaled_features", assay_obj._all_feature_names)
            feat_idx_map = {f: i for i, f in enumerate(scaled_features)}
            valid = [f for f in features if f in feat_idx_map]
            idx = [feat_idx_map[f] for f in valid]
            if sp.issparse(mat):
                return mat[idx, :].toarray().astype(float)
            return np.asarray(mat)[idx, :].astype(float)
        else:
            # Fall back to log-normalized data
            data = assay_obj.layers.get("data") or assay_obj.layers.get("counts")
            all_feats = assay_obj._all_feature_names
            feat_idx = [all_feats.index(f) for f in features if f in all_feats]
            if sp.issparse(data):
                sub = data[feat_idx, :].toarray().astype(float)
            else:
                sub = np.asarray(data)[feat_idx, :].astype(float)
            # Center in-place
            sub -= sub.mean(axis=1, keepdims=True)
            return sub
    else:
        if not is_matrix_empty(assay_obj.scale_data):
            sd = assay_obj.scale_data
            all_feats = assay_obj._feature_names
            feat_idx = [all_feats.index(f) for f in features if f in all_feats]
            if sp.issparse(sd):
                return sd[feat_idx, :].toarray().astype(float)
            return np.asarray(sd)[feat_idx, :].astype(float)
        else:
            # Fall back to data
            all_feats = assay_obj._feature_names
            feat_idx = [all_feats.index(f) for f in features if f in all_feats]
            mat = assay_obj.data
            if sp.issparse(mat):
                sub = mat[feat_idx, :].toarray().astype(float)
            else:
                sub = np.asarray(mat)[feat_idx, :].astype(float)
            sub -= sub.mean(axis=1, keepdims=True)
            return sub
