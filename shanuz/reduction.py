"""Dimensionality reduction: PCA, sPCA, ICA, t-SNE.

Mirrors Seurat's RunPCA() / RunSPCA() / RunICA() / RunTSNE().
"""
from __future__ import annotations

import warnings
from typing import Optional

import numpy as np
import scipy.sparse as sp

from .command import log_shanuz_command
from .dimreduc import DimReduc

# How many dropped features the warning spells out before summarising the rest.
# Seurat names every one; on a query missing a few thousand of a reference's
# variable features that is a wall of text nobody reads, and the count is the
# part that tells you whether to worry.
_MAX_NAMED_FEATURES = 20


def _truncated_svd(mat, n_components: int, seed: int):
    """A deterministic, accurate truncated SVD — R's ``irlba(A, nv = n)``.

    Seurat's ``RunPCA`` calls ``irlba``, which converges to the true leading
    singular triplets. ``sklearn.decomposition.PCA`` does not: once
    ``max(shape) > 500`` it silently switches to a *randomized* solver, which
    is accurate in the leading components and drifts badly in the trailing
    ones. On a 2,400-cell ifnb subsample the randomized solver reproduced only
    15 of 30 PCs above |r| = 0.99 against Seurat — PC 28 came back at 0.006 —
    while ARPACK matched all 30 to six decimals, six times faster than a dense
    ``np.linalg.svd``. Every default ``n_pcs=30`` neighbour graph, clustering
    and UMAP sits downstream of those trailing PCs, and Seurat v5's
    ``IntegrateEmbeddings`` corrects the embedding itself, so the drift lands
    directly in its output.

    ARPACK needs ``k`` comfortably below ``min(shape)``, so small inputs (test
    fixtures, tiny assays) take the exact dense path instead — cheap at that
    size and immune to the convergence failures ARPACK hits when ``k``
    approaches the matrix rank.
    """
    import scipy.sparse.linalg as sla

    if min(mat.shape) <= max(4 * n_components, 50):
        dense = mat.toarray() if sp.issparse(mat) else np.asarray(mat)
        u, d, vt = np.linalg.svd(dense, full_matrices=False)
        u, d, vt = u[:, :n_components], d[:n_components], vt[:n_components]
    else:
        # v0 fixes ARPACK's starting vector, which is what makes this
        # reproducible; svds returns ascending singular values.
        rng = np.random.default_rng(seed)
        u, d, vt = sla.svds(
            mat.asfptype() if sp.issparse(mat) else np.asarray(mat, dtype=float),
            k=n_components,
            v0=rng.standard_normal(min(mat.shape)),
        )
        order = np.argsort(-d)
        u, d, vt = u[:, order], d[order], vt[order]

    # Sign is arbitrary in an SVD, so pin it the way sklearn's svd_flip does:
    # make each component's largest-magnitude cell loading positive. Seurat
    # does not canonicalize, so signs may still differ from R's — compare
    # reductions with |correlation|, and do not "fix" a flipped PC.
    flip = np.sign(u[np.abs(u).argmax(axis=0), np.arange(u.shape[1])])
    flip[flip == 0] = 1.0
    return u * flip, d, vt * flip[:, None]


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
    assay_name = assay or seurat.active_assay
    assay_obj = seurat.assays[assay_name]

    features = _default_features(assay_obj, features)

    # Get scale.data for the selected features. `used` is what the layer
    # actually carried; it is what labels the loadings below, because a feature
    # the layer does not have has no row to be labelled.
    scaled, used = _prep_dr(assay_obj, features, layer)
    # scaled shape: (n_features_selected × n_cells)
    # irlba is handed t(scale.data), so transpose to (n_cells × n_features)
    data_t = scaled.T  # (cells × features)

    n_pcs = min(n_pcs, min(data_t.shape) - 1)

    np.random.seed(seed)
    u, d, vt = _truncated_svd(data_t, n_pcs, seed)
    embeddings = u * d                    # RunPCA.default: u %*% diag(d)
    loadings = vt.T                       # (features × n_pcs), pca.results$v
    # sdev <- pca.results$d / sqrt(max(1, ncol(object) - 1)), where `object` is
    # features × cells — so the denominator counts cells.
    stdev = d / np.sqrt(max(1, data_t.shape[0] - 1))

    # Cell names
    cells = seurat.cell_names()

    dr = DimReduc(
        cell_embeddings=embeddings,
        feature_loadings=loadings,
        assay_used=assay_name,
        stdev=stdev,
        key=reduction_key,
        cell_names=cells,
        feature_names=used,
    )

    seurat.reductions[reduction_name] = dr
    log_shanuz_command(
        seurat, "RunPCA", assay=assay_name,
        params={"n_pcs": n_pcs, "reduction_name": reduction_name, "seed": seed},
    )


def run_spca(
    seurat,
    graph: str,
    npcs: int = 50,
    features: Optional[list[str]] = None,
    assay: Optional[str] = None,
    reduction_name: str = "spca",
    reduction_key: str = "SPC_",
    seed: int = 42,
    layer: str = "scale.data",
) -> None:
    """Supervised PCA — the gene axes that best explain a cell-cell graph.

    Mirrors R's ``RunSPCA(obj, assay = "SCT", graph = "wsnn")``. Ordinary PCA
    picks the directions of greatest variance and knows nothing about which
    cells you consider neighbours. sPCA is handed a graph you already trust —
    typically the WNN graph from `find_multi_modal_neighbors`, which knows about
    protein as well as RNA — and finds the directions *in gene space* that best
    reproduce it. Where PCA maximises ``vᵀXᵀXv``, sPCA maximises ``vᵀXᵀGXv``:
    the same problem with the identity swapped for the graph. Set ``G = I`` and
    you get PCA back exactly.

    The point is the loadings. Because sPCA is still a linear map from genes to
    components, a query dataset can be pushed into a reference's graph-defined
    space with a single matrix multiply, which is what makes it the reduction
    Azimuth maps onto.

    Parameters
    ----------
    graph          : key in ``seurat.graphs`` — a cell × cell graph (e.g. "wsnn")
    npcs           : number of components
    features       : genes to use (defaults to variable features)
    assay          : assay name (defaults to active assay)
    reduction_name : key for storage in seurat.reductions
    reduction_key  : prefix for dimension names (e.g. 'SPC_')
    seed           : random seed (the eigensolver starts from a random vector)
    layer          : which layer to take data from

    Notes
    -----
    Seurat runs `irlba` on ``XᵀGX``, which is an SVD and so ranks components by
    ``|λ|``; we take the largest eigenvalues themselves, since ``vᵀXᵀGXv`` is
    what is being maximised and a graph can push some eigenvalues negative. With
    non-negative edge weights the leading eigenvalues are positive and the two
    orderings agree, so this only ever differs in the tail.
    """
    assay_name = assay or seurat.active_assay
    assay_obj = seurat.assays[assay_name]

    if graph not in seurat.graphs:
        raise KeyError(
            f"Graph '{graph}' not found. Run find_neighbors() or "
            f"find_multi_modal_neighbors() first. "
            f"Available graphs: {list(seurat.graphs)}"
        )

    features = _default_features(assay_obj, features)

    # `_prep_dr` densifies on the way out, so there is nothing left to convert
    # here — the `sp.issparse` branch that used to stand between these two lines
    # could not run.
    scaled, used = _prep_dr(assay_obj, features, layer)
    X = np.asarray(scaled.T, dtype=float)
    n_cells, n_features = X.shape

    G = seurat.graphs[graph].tocsr()
    if G.shape != (n_cells, n_cells):
        raise ValueError(
            f"Graph '{graph}' is {G.shape[0]}×{G.shape[1]} but the assay has "
            f"{n_cells} cells."
        )
    # A KNN graph need not be symmetric; the eigendecomposition needs it to be,
    # and (G + Gᵀ)/2 leaves the quadratic form vᵀXᵀGXv untouched anyway.
    G = (G + G.T) * 0.5

    npcs = min(npcs, n_features - 1)
    if npcs < 1:
        raise ValueError(
            f"Need at least 2 features to run sPCA; got {n_features}.")

    np.random.seed(seed)
    Z = X.T @ (G @ X)                                       # features × features
    Z = np.asarray(Z)
    Z = (Z + Z.T) * 0.5                                     # float asymmetry
    eigenvalues, loadings = _top_eigenvectors(Z, npcs, seed=seed)
    loadings = _flip_signs(loadings)                        # reproducible signs

    embeddings = X @ loadings                               # cells × npcs
    stdev = np.sqrt(np.var(embeddings, axis=0, ddof=1))

    seurat.reductions[reduction_name] = DimReduc(
        cell_embeddings=embeddings,
        feature_loadings=loadings,
        assay_used=assay_name,
        stdev=stdev,
        key=reduction_key,
        cell_names=seurat.cell_names(),
        feature_names=used,
        misc={"spca_graph": graph, "eigenvalues": eigenvalues},
    )


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

    features = _default_features(assay_obj, features)

    scaled, used = _prep_dr(assay_obj, features, layer)  # (features × cells)
    data_t = scaled.T  # (cells × features) — already dense, see run_pca

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
        feature_names=used,
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


def _default_features(assay_obj, features: Optional[list[str]]) -> list[str]:
    """The genes a reduction runs on: whatever was asked for, else the variable ones."""
    from .assay5 import Assay5

    if features is not None:
        return features
    if isinstance(assay_obj, Assay5):
        return assay_obj.variable_features or assay_obj._all_feature_names
    return assay_obj.var_features or assay_obj._feature_names


def _top_eigenvectors(Z: np.ndarray, k: int, seed: int = 42):
    """The k eigenvectors of symmetric ``Z`` with the largest eigenvalues."""
    if k >= Z.shape[0] - 1:
        eigenvalues, vectors = np.linalg.eigh(Z)
    else:
        from scipy.sparse.linalg import eigsh
        v0 = np.random.default_rng(seed).normal(size=Z.shape[0])
        eigenvalues, vectors = eigsh(Z, k=k, which="LA", v0=v0)
    order = np.argsort(eigenvalues)[::-1][:k]
    return eigenvalues[order], vectors[:, order]


def _flip_signs(vectors: np.ndarray) -> np.ndarray:
    """Pin each eigenvector's arbitrary sign so repeat runs agree.

    Convention (sklearn's): the largest-magnitude entry of every column is positive.
    """
    peak = np.argmax(np.abs(vectors), axis=0)
    signs = np.sign(vectors[peak, np.arange(vectors.shape[1])])
    signs[signs == 0] = 1.0
    return vectors * signs


def _layer_matrix_and_features(assay_obj, layer: str):
    """``(matrix, row labels, centre?)`` for ``layer``, or the fallback if absent.

    Each layer is labelled by its *own* feature list. Only ``scale.data`` holds a
    subset — ScaleData defaults to the variable features — so reading any other
    layer's rows through the scaled names, which is what a single shared name
    list amounts to, returns the wrong genes without raising.

    The third element says whether the caller must centre: ``scale.data`` is
    already centred, and the ``data``/``counts`` fallback is not.
    """
    from .assay5 import Assay5
    from ._sparse import is_matrix_empty

    if isinstance(assay_obj, Assay5):
        if layer in assay_obj.layers:
            names = assay_obj._layer_features.get(layer, assay_obj._all_feature_names)
            return assay_obj.layers[layer], list(names), False
        # Fall back to log-normalized data. Selected with an explicit `is None`
        # rather than `a or b`: `or` evaluates `bool(a)`, and scipy raises
        # "truth value of an array ... is ambiguous" for any matrix with more
        # than one entry — so this fallback raised on every call that reached
        # it, which is every `run_pca` on an Assay5 that had not been through
        # `scale_data()`.
        for name in ("data", "counts"):
            mat = assay_obj.layers.get(name)
            if mat is not None:
                names = assay_obj._layer_features.get(name, assay_obj._all_feature_names)
                return mat, list(names), True
        raise ValueError(
            f"assay has no {layer!r}, 'data' or 'counts' layer to draw from; "
            f"run normalize_data() and scale_data() first."
        )

    if not is_matrix_empty(assay_obj.scale_data):
        # `features("scale_data")` — not the assay's full feature list, which is
        # what stood here. scale_data holds only the scaled subset, so indexing
        # it by a gene's position in the *assay* reads a different gene's row,
        # or runs off the end of a matrix a fraction of the height.
        return assay_obj.scale_data, assay_obj.features("scale_data"), False
    return assay_obj.data, list(assay_obj._feature_names), True


def _layer_feature_names(assay_obj, layer: str) -> list[str]:
    """The features a reduction on ``layer`` could actually draw on."""
    return _layer_matrix_and_features(assay_obj, layer)[1]


def _select_features(
    names: list[str],
    features: list[str],
    layer: str,
    action: str = "reduction",
    stacklevel: int = 3,
) -> tuple[list[int], list[str]]:
    """``(row indices, labels)`` for ``features`` within ``names``.

    Mirrors Seurat's ``PrepDR5``: a feature the layer does not carry is dropped
    *with a warning that names it*, and none at all is an error. Both halves are
    the point — a caller that gets back a shorter matrix and no signal has no way
    to tell, and every row below the first drop then answers to the wrong gene.

    ``stacklevel`` counts out to the *user's* call. It is a parameter rather than
    a constant because the two callers sit at different depths, and a warning
    that points at a line inside this module tells the reader nothing.
    """
    pos = {f: i for i, f in enumerate(names)}
    kept = [f for f in features if f in pos]
    if not kept:
        raise ValueError(
            f"None of the {len(features)} requested features are present in "
            f"{layer!r}; nothing to {action}. Run scale_data() over them first."
        )
    if len(kept) < len(features):
        missing = [f for f in features if f not in pos]
        shown = ", ".join(missing[:_MAX_NAMED_FEATURES])
        if len(missing) > _MAX_NAMED_FEATURES:
            shown += f", ... ({len(missing) - _MAX_NAMED_FEATURES} more)"
        warnings.warn(
            f"The following {len(missing)} features are not present in "
            f"{layer!r}; running the {action} without them: {shown}",
            RuntimeWarning,
            stacklevel=stacklevel,
        )
    return [pos[f] for f in kept], kept


def _prep_dr(
    assay_obj, features: list[str], layer: str
) -> tuple[np.ndarray, list[str]]:
    """``(features × cells matrix, the features it actually holds)``.

    The caller is handed back the list that labels the rows it got. Returning
    only the matrix — which is what :func:`_get_scaled_data` did — leaves every
    caller to assume its own *request* labels the result, so a dropped feature
    shifts every row below it onto the wrong gene name in the stored loadings.
    """
    from ._sparse import as_dense

    mat, names, needs_centering = _layer_matrix_and_features(assay_obj, layer)
    idx, kept = _select_features(names, features, layer, stacklevel=4)

    sub = as_dense(mat[idx, :]).astype(float)
    if needs_centering:
        sub -= sub.mean(axis=1, keepdims=True)
    return sub, kept


def _get_scaled_data(assay_obj, features: list[str], layer: str) -> np.ndarray:
    """Extract scaled data for the given features as a (features × cells) array.

    Thin wrapper over :func:`_prep_dr` for callers that genuinely only want the
    numbers. Anything that goes on to *label* the rows must use `_prep_dr` and
    take the feature list with them.
    """
    return _prep_dr(assay_obj, features, layer)[0]
