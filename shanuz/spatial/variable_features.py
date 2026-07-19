"""Spatially variable feature detection.

Mirrors Seurat's ``FindSpatiallyVariableFeatures``: score every gene by how
strongly its expression is organised in space, then rank the genes by that
score. Both of R's methods are available through ``method=``:

``"moransi"``
    Moran's I — one number per gene saying whether neighbouring cells resemble
    each other. Comes with a p-value. By default the spatial weights are R's:
    inverse-squared distance between every pair of cells, row-standardised.
    ``weights="knn"`` swaps in a k-nearest-neighbour graph, which is an
    approximation but drops the cost from O(n²) to O(nk).

``"markvariogram"``
    The normalised mark variogram evaluated at a fixed distance — the average
    squared expression difference between cells that far apart, relative to the
    gene's own variance. No distributional assumption, and it asks a sharper
    question than Moran's I: not "is this gene autocorrelated at all" but "how
    much of its variance has decayed away by distance *r*".

Pure NumPy/SciPy — coordinates come from :func:`get_tissue_coordinates`, so any
object with populated ``images`` works (Xenium / Visium / CosMx / MERSCOPE, or
``from_anndata`` on a spatial ``obsm``).
"""
from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..composition import _bh
from .analysis import get_tissue_coordinates, spatial_knn

METHODS = ("moransi", "markvariogram")
WEIGHTS = ("inverse_square", "knn")

# Distances are computed in row blocks so that the O(n²) weight matrix is never
# materialised. The block is sized to keep one block of distances near this many
# doubles (~80 MB), which is what makes R's all-pairs weighting usable on a slide
# with tens of thousands of cells — R itself allocates the whole n × n matrix and
# runs out of memory well before that.
_BLOCK_DOUBLES = 10_000_000


# ---------------------------------------------------------------------------
# Moran's I
# ---------------------------------------------------------------------------

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


def _block_size(n: int) -> int:
    """How many rows of the weight matrix to hold in memory at once."""
    return int(np.clip(_BLOCK_DOUBLES // max(n, 1), 1, max(n, 1)))


def _inverse_square_base(xy: np.ndarray, start: int, stop: int) -> np.ndarray:
    """Rows ``start:stop`` of ``1 / d²``, with the diagonal already zeroed.

    ``cdist`` is used rather than the ``|a|² + |b|² - 2a·b`` expansion because the
    expansion loses roughly half the mantissa on coordinates far from the origin,
    and slide coordinates are exactly that — Xenium microns run into the tens of
    thousands. Matching R here is a matter of not throwing precision away.
    """
    from scipy.spatial.distance import cdist

    d2 = cdist(xy[start:stop], xy, metric="sqeuclidean")
    # w_ii = 0: R writes the zero in after inverting, which is the same thing.
    np.fill_diagonal(d2[:, start:stop], np.inf)
    if not d2.all():
        raise ValueError(
            "Two cells share the same coordinates, so an inverse-distance weight "
            "is infinite. Deduplicate the coordinates, or pass weights='knn'."
        )
    return 1.0 / d2


def _inverse_square_lag(Z: np.ndarray, xy: np.ndarray) -> tuple[np.ndarray, float, float, float]:
    """Row-standardised inverse-square spatial lag, and the moments of W.

    Returns ``(Y, s0, s1, s2)`` with ``Y[g, i] = Σ_j w_ij · Z[g, j]``. This is
    Seurat's weighting — ``weights <- 1 / pos.dist.mat ^ 2`` in ``RunMoransI``,
    row-standardised inside ``Rfast2::moranI`` — evaluated in row blocks so the
    n × n matrix is never held whole.

    Two passes are needed: w_ij = base_ij / rowsum_i, so every row's denominator
    must be known before any weight is final.
    """
    n = len(xy)
    block = _block_size(n)

    rowsum = np.empty(n)
    for s in range(0, n, block):
        e = min(s + block, n)
        rowsum[s:e] = _inverse_square_base(xy, s, e).sum(axis=1)
    rowsum[rowsum == 0] = 1.0

    Y = np.empty((Z.shape[0], n))
    colsum = np.zeros(n)
    s0 = 0.0
    s1 = 0.0
    for s in range(0, n, block):
        e = min(s + block, n)
        base = _inverse_square_base(xy, s, e)
        W = base / rowsum[s:e, None]
        Y[:, s:e] = Z @ W.T
        colsum += W.sum(axis=0)
        s0 += float(W.sum())
        # w_ij + w_ji = base_ij · (1/rowsum_i + 1/rowsum_j), the base being
        # symmetric — so the transpose never has to be formed.
        sym = base * (1.0 / rowsum[s:e, None] + 1.0 / rowsum[None, :])
        s1 += 0.5 * float(np.einsum("ij,ij->", sym, sym))
    # Every row of a row-standardised W sums to 1, so s2 = Σ_i (1 + colsum_i)².
    s2 = float(np.sum((1.0 + colsum) ** 2))
    return Y, s0, s1, s2


def _morans_i_from_lag(Z: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Moran's I per gene, from the expression and its spatial lag.

    ``I = Σ(y - ȳ)(z - z̄) / Σ(z - z̄)²`` — which is ``Rfast2``'s
    ``cor(y, x) · sd(y) / sd(x)`` with the standard deviations cancelled, the
    same number with one fewer place to lose precision. Valid because the lag was
    built from row-standardised weights, so S₀ = n and the usual ``n / S₀`` factor
    is 1.
    """
    Zc = Z - Z.mean(axis=1, keepdims=True)
    Yc = Y - Y.mean(axis=1, keepdims=True)
    num = np.einsum("gi,gi->g", Yc, Zc)
    den = np.einsum("gi,gi->g", Zc, Zc)
    out = np.full(Z.shape[0], np.nan)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    return out


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
    return _moments_from_sums(n, s0, s1, s2)


def _moments_from_sums(n: int, s0: float, s1: float, s2: float) -> tuple[float, float]:
    """(expected I, variance of I) under normality, from the three weight sums.

    Split out from :func:`_morans_moments` so the blocked inverse-square path can
    accumulate S₀/S₁/S₂ as it goes and share the closed form, rather than holding
    an n × n matrix just to total it up.
    """
    e_i = -1.0 / (n - 1)
    var = (n**2 * s1 - n * s2 + 3 * s0**2) / (s0**2 * (n**2 - 1)) - e_i**2
    return e_i, max(var, 0.0)


def _moransi_table(
    Z: np.ndarray,
    xy: np.ndarray,
    genes: list[str],
    k: int,
    weights: str = "inverse_square",
) -> pd.DataFrame:
    """Moran's I, its p-value and rank, for every row of ``Z``."""
    from scipy.stats import norm

    if weights == "inverse_square":
        Y, s0, s1, s2 = _inverse_square_lag(Z, xy)
        moran = _morans_i_from_lag(Z, Y)
        e_i, var = _moments_from_sums(len(xy), s0, s1, s2)
    else:
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
    res["moransi_rank"] = _rank(res["moransi"], ascending=False)
    return res.sort_values("moransi_rank")


# ---------------------------------------------------------------------------
# Mark variogram
# ---------------------------------------------------------------------------

def _nn_spacing(xy: np.ndarray) -> float:
    """Median nearest-neighbour distance — the natural length unit of a slide.

    Distances are quoted in these units so that ``r_metric`` means the same
    thing whether the coordinates arrive as Visium full-resolution pixels,
    Xenium microns, or anything else.
    """
    d, _ = spatial_knn(xy, k=1)
    spacing = float(np.median(d[:, 0]))
    if not np.isfinite(spacing) or spacing <= 0:
        raise ValueError(
            "Cannot measure a nearest-neighbour spacing: over half the cells sit "
            "on top of another cell."
        )
    return spacing


def _band_weights(xy: np.ndarray, r: float, h: float) -> sp.csr_matrix:
    """Symmetric sparse kernel weights on the cell pairs lying ≈ r apart.

    An Epanechnikov kernel over the distance error, ``w_ij = 1 − ((d_ij − r)/h)²``
    for ``|d_ij − r| < h`` and zero outside — a soft band of pairs around radius
    r, so that pairs land in it smoothly rather than by a hard cutoff.

    Only pairs inside the band are ever materialised. The full pairwise distance
    matrix is never formed, which is what keeps this affordable on a real slide.
    """
    from scipy.spatial import cKDTree

    n = len(xy)
    pairs = cKDTree(xy).query_pairs(r + h, output_type="ndarray")   # unique, i < j
    if pairs.size == 0:
        return sp.csr_matrix((n, n))

    d = np.linalg.norm(xy[pairs[:, 0]] - xy[pairs[:, 1]], axis=1)
    u = (d - r) / h
    inside = np.abs(u) < 1.0
    pairs, u = pairs[inside], u[inside]
    w = 1.0 - u**2

    i, j = pairs[:, 0], pairs[:, 1]
    K = sp.coo_matrix(
        (np.concatenate([w, w]), (np.concatenate([i, j]), np.concatenate([j, i]))),
        shape=(n, n),
    )
    return K.tocsr()


def _mark_variogram(Z: np.ndarray, K: sp.csr_matrix) -> np.ndarray:
    """Normalised mark variogram in one distance band, for every row of ``Z``.

    ``γ(r) = E[ ½·(m_i − m_j)² | d_ij ≈ r ] / Var(m)`` — a kernel-weighted
    (Nadaraya-Watson) average over the pairs in the band, divided by the gene's
    own variance so the scale of expression drops out. γ ≈ 1 means two cells r
    apart differ as much as two cells picked at random, i.e. no spatial
    structure at that distance; γ < 1 means they still resemble each other.

    ``Z`` must be mean-centred, so ``Var(m)`` is its row-wise mean square.

    The pairwise differences are never materialised — that array would be genes
    × pairs. Because K is symmetric with a zero diagonal,
    ``Σ_{i<j} K_ij·(m_i − m_j)² = Σ_i s_i·m_i² − mᵀKm`` with ``s = K·1``, which
    is two sparse products no matter how many pairs the band holds.
    """
    total = float(K.sum())                          # = 2 · Σ_{i<j} w_ij
    if total <= 0:
        raise ValueError(
            "No cell pairs lie in the distance band: lower r_metric or raise "
            "bandwidth."
        )
    s = np.asarray(K.sum(axis=1)).ravel()
    KZ = np.asarray((K @ Z.T).T)                    # genes × cells
    quad = np.einsum("gi,gi->g", Z, KZ)             # mᵀKm
    sq = (Z**2) @ s                                 # Σ_i s_i·m_i²
    gamma = (sq - quad) / total                     # the ½ cancels the 2 in `total`

    var = np.einsum("gi,gi->g", Z, Z) / Z.shape[1]
    out = np.full(Z.shape[0], np.nan)
    ok = var > 0                                    # a flat gene has no variogram
    out[ok] = gamma[ok] / var[ok]
    return out


def _markvariogram_table(
    Z: np.ndarray,
    xy: np.ndarray,
    genes: list[str],
    r_metric: float,
    bandwidth: float,
) -> pd.DataFrame:
    """The normalised mark variogram at ``r_metric``, and its rank."""
    if r_metric <= 0:
        raise ValueError("r_metric must be positive.")
    if bandwidth <= 0:
        raise ValueError("bandwidth must be positive.")

    spacing = _nn_spacing(xy)
    K = _band_weights(xy, r=r_metric * spacing, h=bandwidth * spacing)
    gamma = _mark_variogram(Z, K)

    res = pd.DataFrame({"markvariogram": gamma}, index=genes)
    # Rank 1 = most spatially variable: the *least* variance decayed by distance r.
    res["markvariogram_rank"] = _rank(res["markvariogram"], ascending=True)
    return res.sort_values("markvariogram_rank")


# ---------------------------------------------------------------------------
# Shared
# ---------------------------------------------------------------------------

def _rank(score: pd.Series, ascending: bool) -> pd.Series:
    """Competition rank over a score column, with unscorable genes sent last.

    A gene that is flat across the cells has no spatial score at all (NaN, from a
    zero denominator). It still needs a rank, or it could not be given one at
    all — ``na_option="bottom"`` parks those genes behind every scored gene.
    """
    return score.rank(
        ascending=ascending, method="min", na_option="bottom"
    ).astype(int)


def find_spatially_variable_features(
    seurat,
    features: Optional[list[str]] = None,
    method: str = "moransi",
    k: int = 10,
    weights: str = "inverse_square",
    assay: Optional[str] = None,
    layer: Optional[str] = None,
    image: Optional[Union[str, Sequence[str]]] = None,
    r_metric: float = 5.0,
    bandwidth: float = 1.0,
) -> pd.DataFrame:
    """Rank features by how spatially structured their expression is.

    Mirrors R's ``FindSpatiallyVariableFeatures(obj, method = ...)``.

    Parameters
    ----------
    features : restrict to these genes (default: all in the assay). Passing the
               variable features keeps this fast on large panels.
    method   : ``"moransi"`` (default) or ``"markvariogram"``.
    weights  : *moransi only* — ``"inverse_square"`` (default) reproduces R
               exactly: 1/d² between every pair of cells, row-standardised.
               ``"knn"`` uses a k-nearest-neighbour graph instead, an
               approximation that trades R's O(n²) cost for O(nk). Prefer it only
               when the slide is too large to wait on; see the note below on how
               far the two answers drift apart.
    k        : *moransi only, weights="knn" only* — neighbours per cell.
    r_metric : *markvariogram only* — the distance at which to read the
               variogram, **in units of the median nearest-neighbour spacing**.
               The default of 5 therefore means "five cells apart".
    bandwidth: *markvariogram only* — half-width of the distance band around
               ``r_metric``, in the same units. Widen it if the slide is sparse
               and the band catches too few pairs to average over.
    layer    : expression layer (default: the normalized ``data``).
    image    : image name(s) to draw coordinates from (default: all).

    Returns
    -------
    DataFrame indexed by gene, sorted so that rank 1 is the most spatially
    variable, with the columns for the chosen method:

    ``moransi``
        ``moransi`` (I; +ve = spatially clustered), ``moransi_pval`` (two-sided,
        normality assumption), ``moransi_padj`` (Benjamini-Hochberg) and
        ``moransi_rank``.
    ``markvariogram``
        ``markvariogram`` (γ at ``r_metric``; **lower** = more spatially
        structured, ≈ 1 = none) and ``markvariogram_rank``. There is no p-value —
        the variogram has no closed-form null, and R does not offer one either.

    The same columns are also written into the assay's feature-level metadata,
    as ``find_variable_features`` does.

    Notes
    -----
    Run this on log-normalized ``data`` (the default). Be aware that when a few
    strongly spatial genes dominate a cell's library size, log-normalization
    divides every gene by a spatially-structured total and leaks that structure
    into otherwise-flat genes — inflating their score. That is a property of
    compositional normalization, not of either statistic; it is negligible on
    real panels but can bite on small synthetic ones.

    The ``moransi`` **statistic** matches R to ~1e-15 with the default weights.
    The **p-value** deliberately does not. R runs a 999-permutation test through
    ``Rfast2::moranI``, which on a 248-gene panel returns just 14 distinct values
    and ties 233 genes at its 1/1025 floor — it cannot separate the most spatially
    variable gene from the two-hundredth. The normal-approximation p-value here is
    continuous, deterministic, and standard for Moran's I, so it is kept; matching
    R would cost information and gain nothing. ``weights="knn"`` changes the
    statistic itself, not just its cost. It is a decent approximation — on 2,000
    Xenium cells it tracks R at Pearson 0.986 and recovers 46 of R's top 50 genes
    — but it runs a median 1.23× high, differs by up to 0.14 in absolute I, and
    agrees on only 7 of R's top 10. Close enough to mislead, not close enough to
    call parity, which is why it is no longer the default.

    ``markvariogram`` differs from R's in two deliberate ways. R passes
    ``r.metric`` straight through to ``spatstat`` in raw coordinate units, so the
    same script gives different answers on pixel and micron coordinates; here r
    is measured in nearest-neighbour spacings and is scale-free. And γ is a
    kernel-weighted ratio estimator rather than ``spatstat``'s
    translation-corrected one, so the absolute γ values are close to but not
    identical with R's. The gene *ranking* — the thing the function is for — is
    what carries over.
    """
    from ..markers import _get_expression_matrix

    if method not in METHODS:
        raise NotImplementedError(
            f"method={method!r} is not implemented; use one of {list(METHODS)}."
        )
    if weights not in WEIGHTS:
        raise ValueError(
            f"weights={weights!r} is not recognised; use one of {list(WEIGHTS)}."
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

    if method == "moransi":
        res = _moransi_table(Z, xy, genes, k=k, weights=weights)
    else:
        res = _markvariogram_table(Z, xy, genes, r_metric=r_metric, bandwidth=bandwidth)

    _write_feature_meta(assay_obj, res)
    return res


def _write_feature_meta(assay_obj, res: pd.DataFrame) -> None:
    """Store the score columns in the assay's feature-level metadata."""
    meta = getattr(assay_obj, "meta_features", None)
    attr = "meta_features"
    if meta is None:
        meta = getattr(assay_obj, "meta_data", None)   # Assay5
        attr = "meta_data"
    if meta is None:
        return
    for col in res.columns:
        meta.loc[res.index, col] = res[col]
    setattr(assay_obj, attr, meta)
