"""Spatially variable feature detection.

Mirrors Seurat's ``FindSpatiallyVariableFeatures(obj, method = "moransi")``:
score each gene by its spatial autocorrelation (Moran's I) over a k-nearest-
neighbour graph built from the cells' tissue coordinates.

Pure NumPy/SciPy — the weight matrix comes from :func:`shanuz.spatial.spatial_knn`,
so any object with populated ``images`` works (Xenium / Visium / CosMx / MERSCOPE,
or ``from_anndata`` on a spatial ``obsm``).
"""
from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..composition import _bh
from .analysis import get_tissue_coordinates, spatial_knn


def _knn_weights(coords: np.ndarray, k: int, row_normalize: bool = True) -> sp.csr_matrix:
    """Sparse spatial weight matrix W from a k-nearest-neighbour graph.

    ``w_ij = 1`` when j is among i's k nearest neighbours (self excluded). With
    ``row_normalize`` each row sums to 1 — the usual row-standardised weights, so
    Moran's I is a plain neighbour-average correlation.
    """
    n = len(coords)
    if n < 3:
        raise ValueError("Need at least 3 cells to compute spatial autocorrelation.")
    k = min(k, n - 1)
    _, idx = spatial_knn(coords, k=k)
    rows = np.repeat(np.arange(n), idx.shape[1])
    cols = np.asarray(idx).ravel()
    W = sp.csr_matrix((np.ones(rows.size), (rows, cols)), shape=(n, n))
    W.setdiag(0.0)
    W.eliminate_zeros()
    if row_normalize:
        rs = np.asarray(W.sum(axis=1)).ravel()
        rs[rs == 0] = 1.0
        W = sp.diags(1.0 / rs) @ W
    return sp.csr_matrix(W)


def _morans_i(Z: np.ndarray, W: sp.csr_matrix) -> np.ndarray:
    """Moran's I for every row of ``Z`` (genes × cells, already mean-centred).

    ``I = (N / S0) · (zᵀ W z) / (zᵀ z)``, vectorised across genes.
    """
    n = W.shape[0]
    s0 = float(W.sum())
    # (W @ Z.T).T[g, i] = Σ_j W[i, j] · Z[g, j]  →  row-wise zᵀWz below.
    WZ = np.asarray((W @ Z.T).T)
    num = np.einsum("gi,gi->g", Z, WZ)
    den = np.einsum("gi,gi->g", Z, Z)
    out = np.full(Z.shape[0], np.nan)
    ok = den > 0
    out[ok] = (n / s0) * num[ok] / den[ok]
    return out


def _morans_moments(W: sp.csr_matrix) -> tuple[float, float]:
    """(expected I, variance of I) under the normality assumption.

    Depends only on the weight matrix, so it is computed once for all genes.
    """
    n = W.shape[0]
    s0 = float(W.sum())
    A = W + W.T
    s1 = 0.5 * float(A.multiply(A).sum())
    row = np.asarray(W.sum(axis=1)).ravel()
    col = np.asarray(W.sum(axis=0)).ravel()
    s2 = float(np.sum((row + col) ** 2))
    e_i = -1.0 / (n - 1)
    var = (n**2 * s1 - n * s2 + 3 * s0**2) / (s0**2 * (n**2 - 1)) - e_i**2
    return e_i, max(var, 0.0)


def find_spatially_variable_features(
    seurat,
    features: Optional[list[str]] = None,
    method: str = "moransi",
    k: int = 10,
    assay: Optional[str] = None,
    layer: Optional[str] = None,
    image: Optional[Union[str, Sequence[str]]] = None,
) -> pd.DataFrame:
    """Rank features by spatial autocorrelation (Moran's I).

    Mirrors R's ``FindSpatiallyVariableFeatures(obj, method = "moransi")``.

    Parameters
    ----------
    features : restrict to these genes (default: all in the assay). Passing the
               variable features keeps this fast on large panels.
    method   : only ``"moransi"`` is implemented (``"markvariogram"`` is planned).
    k        : neighbours per cell in the spatial weight graph.
    layer    : expression layer (default: the normalized ``data``).
    image    : image name(s) to draw coordinates from (default: all).

    Returns
    -------
    DataFrame indexed by gene with ``moransi`` (I, +ve = spatially clustered),
    ``moransi_pval`` (two-sided, normality assumption), ``moransi_padj``
    (Benjamini-Hochberg) and ``moransi_rank`` (1 = most spatially variable),
    sorted by rank. The same columns are also written into the assay's
    feature-level metadata, as ``find_variable_features`` does.

    Notes
    -----
    Run this on log-normalized ``data`` (the default). Be aware that when a few
    strongly spatial genes dominate a cell's library size, log-normalization
    divides every gene by a spatially-structured total and leaks that structure
    into otherwise-flat genes — inflating their I. That is a property of
    compositional normalization, not of the statistic; it is negligible on real
    panels but can bite on small synthetic ones.
    """
    from scipy.stats import norm

    from ..markers import _get_expression_matrix

    if method != "moransi":
        raise NotImplementedError(
            f"method={method!r} is not implemented; use 'moransi'."
        )

    coords = get_tissue_coordinates(seurat, image)
    if coords.empty:
        raise ValueError("Object has no spatial coordinates.")
    # A cell can appear once per image; keep its first placement.
    coords = coords.drop_duplicates(subset="cell")

    assay_name = assay or seurat.active_assay
    assay_obj = seurat.assays[assay_name]
    data, feature_names = _get_expression_matrix(assay_obj, layer)

    pos = {c: i for i, c in enumerate(seurat.cell_names())}
    keep = [c for c in coords["cell"] if c in pos]
    if len(keep) < 3:
        raise ValueError("Fewer than 3 cells have both coordinates and expression.")
    col_idx = [pos[c] for c in keep]
    xy = coords.set_index("cell").loc[keep, ["x", "y"]].to_numpy(dtype=float)

    if features is not None:
        want = set(features)
        rows = [i for i, f in enumerate(feature_names) if f in want]
        if not rows:
            raise ValueError("None of the requested features are in the assay.")
    else:
        rows = list(range(len(feature_names)))
    genes = [feature_names[i] for i in rows]

    sub = data[rows, :][:, col_idx]
    X = np.asarray(sub.toarray() if sp.issparse(sub) else sub, dtype=float)
    Z = X - X.mean(axis=1, keepdims=True)

    W = _knn_weights(xy, k=k)
    moran = _morans_i(Z, W)
    e_i, var = _morans_moments(W)

    if var > 0:
        z = (moran - e_i) / np.sqrt(var)
        pval = 2.0 * norm.sf(np.abs(z))
    else:
        pval = np.ones_like(moran)
    pval = np.where(np.isnan(moran), 1.0, pval)

    res = pd.DataFrame(
        {
            "moransi": moran,
            "moransi_pval": pval,
            "moransi_padj": _bh(pval),
        },
        index=genes,
    )
    # Rank 1 = most spatially variable (highest positive autocorrelation).
    res["moransi_rank"] = res["moransi"].rank(ascending=False, method="min").astype(int)
    res = res.sort_values("moransi_rank")

    _write_feature_meta(assay_obj, res)
    return res


def _write_feature_meta(assay_obj, res: pd.DataFrame) -> None:
    """Store the Moran's I columns in the assay's feature-level metadata."""
    meta = getattr(assay_obj, "meta_features", None)
    attr = "meta_features"
    if meta is None:
        meta = getattr(assay_obj, "meta_data", None)   # Assay5
        attr = "meta_data"
    if meta is None:
        return
    for col in ("moransi", "moransi_pval", "moransi_padj", "moransi_rank"):
        meta.loc[res.index, col] = res[col]
    setattr(assay_obj, attr, meta)
