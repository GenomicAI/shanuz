"""Preprocessing functions for single-cell RNA-seq analysis.

Mirrors Seurat's NormalizeData(), FindVariableFeatures(),
ScaleData(), and PercentageFeatureSet().
"""
from __future__ import annotations

import re
from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .command import log_shanuz_command
from .lazy import is_lazy


# Cells per block when streaming an on-disk layer. Peak memory for the
# streaming paths is a function of this, not of the dataset size.
_CELL_BLOCK = 10_000


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def _cell_blocks(matrix, block_size: int = _CELL_BLOCK):
    """Yield ``matrix`` in blocks of ``block_size`` cells, as CSC blocks.

    Accepts a :class:`~shanuz.lazy.LazyMatrix` or a scipy sparse matrix, so a
    reduction written once serves both the on-disk and the in-memory layer.
    That matters beyond tidiness: statistics computed by two implementations
    agree only to rounding, and where a *tie-break* consumes them, rounding
    decides which features come back.
    """
    if is_lazy(matrix):
        for _, _, block in matrix.col_blocks(block_size):
            yield block
        return
    csc = matrix.tocsc()
    for start in range(0, csc.shape[1], block_size):
        yield csc[:, start:min(start + block_size, csc.shape[1])]


def _get_assay(seurat, assay: Optional[str]):
    return seurat.assays[assay or seurat.active_assay]


def _get_layer(assay_obj, layer: Optional[str]):
    """Get a data matrix from an assay (Assay or Assay5)."""
    from .assay import Assay
    from .assay5 import Assay5

    if isinstance(assay_obj, Assay5):
        if layer is None:
            for candidate in ("data", "counts"):
                if candidate in assay_obj.layers:
                    return assay_obj.layers[candidate]
            return assay_obj.layers[assay_obj.default_layer]
        return assay_obj.layers[layer]
    else:
        if layer == "counts":
            return assay_obj.counts
        elif layer == "scale_data" or layer == "scale.data":
            return assay_obj.scale_data
        else:
            return assay_obj.data


def _set_layer(assay_obj, layer: str, value) -> None:
    """Store a data matrix in an assay."""
    from .assay import Assay
    from .assay5 import Assay5

    if isinstance(assay_obj, Assay5):
        assay_obj.set_layer_data(layer, value)
    else:
        if layer == "counts":
            assay_obj.counts = value
        elif layer in ("scale_data", "scale.data"):
            assay_obj.scale_data = np.asarray(value)
        else:
            assay_obj.data = value


# ------------------------------------------------------------------
# PercentageFeatureSet  (mirrors R PercentageFeatureSet())
# ------------------------------------------------------------------

def percentage_feature_set(
    seurat,
    pattern: str,
    col_name: Optional[str] = None,
    assay: Optional[str] = None,
    layer: str = "counts",
) -> None:
    """Add a metadata column with % of counts matching a gene name pattern.

    Mirrors R's PercentageFeatureSet(pbmc, pattern = "^MT-").
    Modifies seurat.meta_data in-place.
    """
    assay_obj = _get_assay(seurat, assay)
    mat = _get_layer(assay_obj, layer)

    from .assay import Assay
    from .assay5 import Assay5
    if isinstance(assay_obj, Assay5):
        feature_names = assay_obj._all_feature_names
    else:
        feature_names = assay_obj._feature_names

    # Find matching features
    rx = re.compile(pattern)
    match_mask = np.array([bool(rx.search(f)) for f in feature_names])

    if not match_mask.any():
        pct = np.zeros(mat.shape[1])
    else:
        # A lazy layer belongs on the sparse branch: indexing it returns a
        # sparse block, so it satisfies this branch's contract exactly, while
        # the dense one would `np.asarray` the whole store just to sum it.
        if sp.issparse(mat) or is_lazy(mat):
            total = np.array(mat.sum(axis=0)).flatten()
            matching = np.array(mat[match_mask, :].sum(axis=0)).flatten()
        else:
            total = mat.sum(axis=0)
            matching = mat[match_mask, :].sum(axis=0)
        total[total == 0] = 1
        pct = (matching / total) * 100.0

    if col_name is None:
        col_name = "percent.mt" if re.search(r"mt|mito", pattern, re.I) else "percent_feature"

    cells = seurat.cell_names()
    if isinstance(assay_obj, Assay5):
        cell_list = assay_obj._all_cell_names
    else:
        cell_list = assay_obj._cell_names

    seurat.meta_data[col_name] = pd.Series(pct, index=cell_list).reindex(cells).values


# ------------------------------------------------------------------
# NormalizeData  (mirrors R NormalizeData())
# ------------------------------------------------------------------

def normalize_data(
    seurat,
    normalization_method: str = "LogNormalize",
    scale_factor: float = 10000.0,
    assay: Optional[str] = None,
    margin: int = 1,
) -> None:
    """Log-normalize counts.

    Mirrors R's NormalizeData(pbmc, normalization.method = "LogNormalize",
    scale.factor = 10000). Modifies the assay in-place by adding / updating
    the 'data' layer.

    Parameters
    ----------
    normalization_method : 'LogNormalize', 'CLR', or 'RC'
    margin               : for CLR only, and matching Seurat's flag exactly —
                           normalize each feature across cells (1, Seurat's
                           default) or each cell across its features (2).
                           ADT/CITE-seq panels typically use margin=2.
    """
    assay_name = assay or seurat.active_assay
    assay_obj = seurat.assays[assay_name]
    counts = _get_layer(assay_obj, "counts")

    if normalization_method == "LogNormalize":
        normed = _log_normalize(counts, scale_factor)
    elif normalization_method == "CLR":
        normed = _clr_normalize(counts, margin=margin)
    elif normalization_method == "RC":
        normed = _rc_normalize(counts, scale_factor)
    else:
        raise ValueError(f"Unknown normalization_method: {normalization_method!r}")

    _set_layer(assay_obj, "data", normed)
    log_shanuz_command(
        seurat, "NormalizeData", assay=assay or seurat.active_assay,
        params={"normalization_method": normalization_method,
                "scale_factor": scale_factor, "margin": margin},
    )


def _log_normalize_lazy(counts, scale_factor: float, block_size: int = _CELL_BLOCK):
    """Stream LogNormalize over an on-disk layer, one cell-block at a time.

    Mirrors Seurat's ``LogNormalize.IterableMatrix``, which never materialises
    the matrix — ``t(t(data)/colSums(data))`` followed by ``log1p(data * scale)``
    stays a queued BPCells operation throughout.

    The whole transform is **value-only**: dividing a column by its total and
    taking ``log1p`` maps zero to zero, so which entries are non-zero never
    changes. That lets the sparsity pattern be carried straight across and only
    the values recomputed, block by block, so peak memory is one block plus the
    output rather than the dense ``features × cells`` array the generic path
    would build.
    """
    totals = counts.sum(axis=0)
    totals[totals == 0] = 1.0

    nnz = counts.nnz
    out_data = np.empty(nnz, dtype=float)
    out_indices = np.empty(nnz, dtype=np.int64)
    out_indptr = np.zeros(counts.ncol + 1, dtype=np.int64)

    offset = 0
    for start, stop, block in counts.col_blocks(block_size):
        end = offset + block.nnz
        # One scale factor per column, expanded to that column's non-zeros --
        # the same product `counts.dot(diags(scale / totals))` forms, in the
        # same order, so the result is bit-identical to the in-memory path.
        per_column = scale_factor / totals[start:stop]
        per_nonzero = np.repeat(per_column, np.diff(block.indptr))
        out_data[offset:end] = np.log1p(block.data.astype(float) * per_nonzero)
        out_indices[offset:end] = block.indices
        out_indptr[start + 1:stop + 1] = block.indptr[1:] + offset
        offset = end

    return sp.csc_matrix(
        (out_data, out_indices, out_indptr), shape=counts.shape
    )


def _log_normalize(counts, scale_factor: float = 10000.0):
    """Normalize counts per cell to scale_factor, then log1p."""
    if is_lazy(counts):
        return _log_normalize_lazy(counts, scale_factor)
    if sp.issparse(counts):
        cell_totals = np.array(counts.sum(axis=0)).flatten().astype(float)
        cell_totals[cell_totals == 0] = 1.0
        # Column-wise scaling via diagonal matrix
        scaler = sp.diags(scale_factor / cell_totals)
        normed = counts.dot(scaler).tocsc()
        normed = normed.astype(float)
        normed.data = np.log1p(normed.data)
        return normed
    else:
        counts = np.asarray(counts, dtype=float)
        cell_totals = counts.sum(axis=0)
        cell_totals[cell_totals == 0] = 1.0
        normed = counts / cell_totals[np.newaxis, :] * scale_factor
        return np.log1p(normed)


def _clr_normalize(counts, margin: int = 1):
    """Centered log-ratio normalization (for protein / ADT assays).

    Faithfully reproduces Seurat's CLR (``clr_function`` in NormalizeData):

        clr(x) = log1p( x / exp( sum(log1p(x[x > 0])) / length(x) ) )

    i.e. the geometric-mean denominator sums log1p over the *non-zero* entries
    but divides by the full length (zeros included). ``counts`` is
    features × cells; ``margin`` selects the direction, matching the axis R's
    ``CustomNormalize`` passes to ``apply(data, MARGIN = margin, clr_function)``:
      * 1 — normalize each feature across all cells (denominator over a row).
            R's ``apply`` MARGIN=1 is rows. Seurat's default.
      * 2 — normalize each cell across its features (denominator over a column).
            R's ``apply`` MARGIN=2 is columns. Recommended for small ADT panels.
    """
    if sp.issparse(counts):
        counts = counts.toarray().astype(float)
    else:
        counts = np.asarray(counts, dtype=float)

    log1p_pos = np.where(counts > 0, np.log1p(counts), 0.0)
    if margin == 2:
        # per-cell (column): length = number of features
        length = counts.shape[0]
        geo = np.exp(log1p_pos.sum(axis=0, keepdims=True) / length)
    else:
        # per-feature (row): length = number of cells
        length = counts.shape[1]
        geo = np.exp(log1p_pos.sum(axis=1, keepdims=True) / length)
    # geo >= 1 always (log1p >= 0), so the division is always well-defined.
    return np.log1p(counts / geo)


def _rc_normalize(counts, scale_factor: float = 10000.0):
    """Relative counts: counts / total * scale_factor (no log)."""
    if sp.issparse(counts):
        cell_totals = np.array(counts.sum(axis=0)).flatten().astype(float)
        cell_totals[cell_totals == 0] = 1.0
        scaler = sp.diags(scale_factor / cell_totals)
        return counts.dot(scaler).tocsc()
    else:
        counts = np.asarray(counts, dtype=float)
        cell_totals = counts.sum(axis=0)
        cell_totals[cell_totals == 0] = 1.0
        return counts / cell_totals[np.newaxis, :] * scale_factor


# ------------------------------------------------------------------
# FindVariableFeatures  (mirrors R FindVariableFeatures())
# ------------------------------------------------------------------

def find_variable_features(
    seurat,
    selection_method: str = "vst",
    nfeatures: int = 2000,
    assay: Optional[str] = None,
    layer: Optional[str] = None,
    mean_cutoff: tuple = (0.1, 8),
    dispersion_cutoff: tuple = (1, float("inf")),
) -> None:
    """Select highly variable features.

    Mirrors R's FindVariableFeatures(pbmc, selection.method = "vst",
    nfeatures = 2000). Modifies the assay in-place by setting
    var_features (Assay) or highly_variable in meta_data (Assay5).
    """
    assay_name = assay or seurat.active_assay
    assay_obj = seurat.assays[assay_name]

    if layer is not None:
        data = _get_layer(assay_obj, layer)
    else:
        from .assay5 import Assay5
        from ._sparse import is_matrix_empty
        if selection_method == "vst":
            # Seurat fits the vst mean-variance LOESS on RAW COUNTS, regardless of
            # whether NormalizeData() has run. Using normalized data here would
            # change which features are selected.
            if isinstance(assay_obj, Assay5):
                data = assay_obj.layers.get("counts")
                if data is None:
                    data = assay_obj.layers.get("data")
            else:
                data = assay_obj.counts
        else:
            # dispersion / mean.var.plot methods operate on log-normalized data
            if isinstance(assay_obj, Assay5):
                data = assay_obj.layers.get("data")
                if data is None:
                    data = assay_obj.layers.get("counts")
            else:
                data = assay_obj.data if not is_matrix_empty(assay_obj.data) else assay_obj.counts

    if selection_method == "vst":
        hvg_indices, means, variances, var_std = _vst_hvg(data, nfeatures)
    elif selection_method == "dispersion" or selection_method == "mvp":
        hvg_indices, means, variances, var_std = _dispersion_hvg(
            data, nfeatures, mean_cutoff, dispersion_cutoff
        )
    elif selection_method == "mean.var.plot":
        hvg_indices, means, variances, var_std = _dispersion_hvg(
            data, nfeatures, mean_cutoff, dispersion_cutoff
        )
    else:
        raise ValueError(f"Unknown selection_method: {selection_method!r}")

    from .assay5 import Assay5
    from .assay import Assay
    if isinstance(assay_obj, Assay5):
        feature_names = assay_obj._all_feature_names
    else:
        feature_names = assay_obj._feature_names

    hvg_names = [feature_names[i] for i in hvg_indices]

    # Store results
    if isinstance(assay_obj, Assay5):
        # Store HVF info in assay meta_data
        hvf_df = pd.DataFrame(
            {
                "means": means,
                "variances": variances,
                "variances.standardized": var_std,
                "highly_variable": np.zeros(len(feature_names), dtype=bool),
            },
            index=feature_names,
        )
        hvf_df.loc[hvg_names, "highly_variable"] = True
        for col in hvf_df.columns:
            assay_obj.meta_data[col] = hvf_df[col]
        assay_obj.variable_features = hvg_names
    else:
        assay_obj.var_features = hvg_names
    log_shanuz_command(
        seurat, "FindVariableFeatures", assay=assay or seurat.active_assay,
        params={"selection_method": selection_method, "nfeatures": nfeatures},
    )


def _vst_hvg(
    data,
    nfeatures: int = 2000,
    loess_span: float = 0.3,
    clip_max: Optional[float] = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Variance-stabilizing transformation HVG selection (Seurat's 'vst').

    Faithfully reproduces Seurat's algorithm on the raw counts:
      1. Per-gene mean and *sample* (N-1) variance of the counts.
      2. Fit a degree-2 LOESS of log10(variance) ~ log10(mean) to get the
         expected variance for each gene given its mean.
      3. Standardize each value with z = (x - mean) / expected_sd, clipping z
         to a maximum of ``clip_max`` (Seurat default sqrt(n_cells)). This clip
         is what stops a single high-count outlier cell from dominating, and is
         essential for counts-based VST to recover biologically variable genes.
      4. variance.standardized = sum(z_clipped^2) / (n_cells - 1).
      5. Rank genes by variance.standardized (descending).
    """
    n_genes, n_cells = data.shape
    if clip_max is None:
        clip_max = np.sqrt(n_cells)

    blockwise = is_lazy(data) or sp.issparse(data)
    if blockwise:
        # Seurat's VST.IterableMatrix gets these from BPCells::matrix_stats, a
        # streaming row reduction. Every non-zero carries its gene index in
        # `indices`, so both moments accumulate per gene across cell-blocks
        # without the matrix ever being whole in memory.
        #
        # Sparse layers take this path too, deliberately. `variance.standardized`
        # has exact ties, and which of two tied genes is selected is decided by
        # a tie-break -- so two code paths agreeing only to 1e-14 would rank
        # them differently, and an on-disk layer would silently return a
        # different feature set (and, through PCA, a different clustering) from
        # the in-memory one. One implementation cannot disagree with itself.
        row_sum = np.zeros(n_genes)
        row_sum_sq = np.zeros(n_genes)
        for block in _cell_blocks(data, _CELL_BLOCK):
            vals = block.data.astype(float)
            row_sum += np.bincount(block.indices, weights=vals, minlength=n_genes)
            row_sum_sq += np.bincount(block.indices, weights=vals * vals,
                                      minlength=n_genes)
        means = row_sum / n_cells
        mean_sq = row_sum_sq / n_cells
    else:
        d = np.asarray(data, dtype=float)
        means = d.mean(axis=1)
        mean_sq = (d ** 2).mean(axis=1)

    # Sample variance (denominator N-1), matching Seurat's SparseRowVar2.
    variances = (mean_sq - means ** 2) * (n_cells / (n_cells - 1))

    valid = variances > 0
    valid_idx = np.where(valid)[0]
    lm_valid = np.log10(means[valid_idx])
    lv_valid = np.log10(variances[valid_idx])

    # Seurat uses loess(degree = 2) for vst; use the local-quadratic fitter.
    fitted_valid = _loess2(lm_valid, lv_valid, frac=loess_span)

    expected_var = np.zeros(n_genes)
    expected_var[valid_idx] = 10.0 ** fitted_valid
    expected_sd = np.sqrt(expected_var)

    var_standardized = np.zeros(n_genes)
    if blockwise:
        # Accumulated per gene across cell-blocks rather than read off one
        # gene's contiguous row, for the same reason as the moments above.
        # Genes with no fitted SD are masked out at the end, matching the
        # `continue` in the dense loop.
        scorable = np.zeros(n_genes, dtype=bool)
        scorable[valid_idx] = expected_sd[valid_idx] > 0
        sd_safe = np.where(scorable, expected_sd, 1.0)

        clipped_sq = np.zeros(n_genes)
        nnz_per_gene = np.zeros(n_genes, dtype=np.int64)
        for block in _cell_blocks(data, _CELL_BLOCK):
            rows = block.indices
            nnz_per_gene += np.bincount(rows, minlength=n_genes)
            z = (block.data.astype(float) - means[rows]) / sd_safe[rows]
            np.minimum(z, clip_max, out=z)
            clipped_sq += np.bincount(rows, weights=z * z, minlength=n_genes)

        # Zero entries standardize to -mu/sd (negative → never clipped).
        zero_term = (n_cells - nnz_per_gene) * (means / sd_safe) ** 2
        var_standardized = np.where(
            scorable, (zero_term + clipped_sq) / (n_cells - 1), 0.0
        )
    else:
        for g in valid_idx:
            sd = expected_sd[g]
            if sd <= 0:
                continue
            z = (d[g, :] - means[g]) / sd
            np.minimum(z, clip_max, out=z)
            var_standardized[g] = np.dot(z, z) / (n_cells - 1)

    # `argsort(-v, stable)` rather than `argsort(v)[::-1]`: R's
    # `head(order(x, decreasing = TRUE), n)` breaks ties by ascending original
    # index, and reversing an ascending sort breaks them descending. Ties are
    # common enough here to straddle the cutoff and decide which genes are
    # selected at all.
    top_idx = np.argsort(-var_standardized, kind="stable")[:nfeatures]

    return top_idx, means, variances, var_standardized


def _loess_fit(x: np.ndarray, y: np.ndarray, frac: float = 0.3) -> np.ndarray:
    """Fit a LOESS curve using statsmodels (with robustness iterations).

    Returns predicted y values at the input x points (in original order).
    Uses bisquare robustness iterations to down-weight outliers, matching
    robust LOESS behavior for HVG detection.
    """
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess
        return lowess(y, x, frac=frac, it=3, delta=0.0, return_sorted=False)
    except ImportError:
        return _loess2(x, y, frac=frac)


def _loess2(x: np.ndarray, y: np.ndarray, frac: float = 0.3, batch_size: int = 200) -> np.ndarray:
    """Degree-2 (local quadratic) LOESS matching R's loess(degree=2, family='gaussian').

    For each query point x₀, fits:
        ŷ(x₀) = β₀  from  min Σᵢ w(xᵢ)·[yᵢ − β₀ − β₁(xᵢ−x₀) − β₂(xᵢ−x₀)²]²
    where w is the tricube kernel over the k nearest neighbours (k = ceil(frac·n)).

    The fit is evaluated once per **distinct** x and broadcast back, and each
    neighbourhood is chosen by distance rather than by position in the sorted
    array. Both are needed to hold R's invariant that a fitted value is a
    function of x alone: ties are common here (on pbmc3k 85 % of genes share a
    ``log10(mean)`` with another gene, the largest run being 627 genes), and
    windowing by position hands the members of a tied run different
    neighbourhoods, so their fits differ according to how an unstable sort
    happened to order them. Neighbourhoods still span observations, duplicates
    included, which is what ``span`` means to R.
    """
    n = len(x)
    k = max(int(np.ceil(frac * n)), 6)

    # Lexicographic, so that a tied run is ordered by y rather than by however
    # the caller happened to supply it. Without that the window sums accumulate
    # in a caller-dependent order and the fit moves in its last bits when the
    # rows are merely permuted.
    sort_idx = np.lexsort((y, x))
    xs = x[sort_idx]
    ys = y[sort_idx]

    # Query points are the distinct x values; `inverse` maps every observation
    # back to the fit computed for its own x.
    xu, inverse = np.unique(x, return_inverse=True)
    m = len(xu)
    fitted_u = np.empty(m)

    # Window start per query point, by distance. The optimal start is
    # non-decreasing in x₀, so one left-to-right sweep finds them all.
    starts = np.empty(m, dtype=np.intp)
    left = 0
    limit = max(n - k, 0)
    for q in range(m):
        x0 = xu[q]
        while left < limit and (x0 - xs[left]) > (xs[left + k] - x0):
            left += 1
        starts[q] = left

    for b0 in range(0, m, batch_size):
        b1 = min(b0 + batch_size, m)
        nb = b1 - b0

        lefts = starts[b0:b1]
        win_idx = lefts[:, np.newaxis] + np.arange(min(k, n))   # (nb, k)
        xw = xs[win_idx]   # (nb, k)
        yw = ys[win_idx]   # (nb, k)
        xi = xu[b0:b1, np.newaxis]   # (nb, 1)

        # Tricube weights
        d = np.abs(xw - xi)
        dmax = np.maximum(d.max(axis=1, keepdims=True), 1e-10)
        u = d / dmax
        w = np.maximum(1.0 - u ** 3, 0.0) ** 3   # (nb, k)

        # Local quadratic design matrix: columns = [1, (x-x₀), (x-x₀)²]
        dx = xw - xi   # (nb, k)
        A = np.stack([np.ones_like(dx), dx, dx ** 2], axis=2)   # (nb, k, 3)

        # Weighted normal equations: (AᵀWA)β = AᵀWy
        Aw = A * w[:, :, np.newaxis]   # (nb, k, 3)
        AtWA = np.einsum("nki,nkj->nij", Aw, A)   # (nb, 3, 3)
        AtWy = np.einsum("nki,nk->ni", Aw, yw)    # (nb, 3)

        try:
            # solve needs RHS as (..., M, 1), returns (..., M, 1)
            betas = np.linalg.solve(AtWA, AtWy[:, :, np.newaxis])[:, :, 0]  # (nb, 3)
            fitted_u[b0:b1] = betas[:, 0]
        except np.linalg.LinAlgError:
            for j in range(nb):
                try:
                    sol, *_ = np.linalg.lstsq(AtWA[j], AtWy[j], rcond=None)
                    fitted_u[b0 + j] = sol[0]
                except Exception:
                    fitted_u[b0 + j] = float(np.average(yw[j], weights=w[j]))

    return fitted_u[inverse]


def _dispersion_hvg(
    data,
    nfeatures: int = 2000,
    mean_cutoff: tuple = (0.1, 8),
    dispersion_cutoff: tuple = (1, float("inf")),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Mean-variance-plot (dispersion) based HVG selection.

    Bins genes by mean expression, normalizes dispersion within each bin.
    Mirrors Seurat v2's FindVariableGenes.
    """
    if sp.issparse(data):
        n_cells = data.shape[1]
        means = np.array(data.mean(axis=1)).flatten()
        mean_sq = np.array(data.power(2).mean(axis=1)).flatten()
    else:
        data = np.asarray(data, dtype=float)
        n_cells = data.shape[1]
        means = data.mean(axis=1)
        mean_sq = (data ** 2).mean(axis=1)

    # Sample variance (N-1 denominator), consistent with the vst path.
    variances = (mean_sq - means ** 2) * (n_cells / (n_cells - 1))
    eps = 1e-10
    dispersion = np.log(variances / (means + eps) + eps)

    # Bin by log mean
    log_means = np.log(means + eps)
    n_bins = 20
    bins = np.percentile(log_means, np.linspace(0, 100, n_bins + 1))
    bins[-1] += 1e-6
    bin_assign = np.digitize(log_means, bins) - 1
    bin_assign = np.clip(bin_assign, 0, n_bins - 1)

    dispersion_scaled = np.zeros(len(means))
    for b in range(n_bins):
        mask = bin_assign == b
        if mask.sum() > 1:
            d = dispersion[mask]
            std = d.std()
            if std > 0:
                dispersion_scaled[mask] = (d - d.mean()) / std
            else:
                dispersion_scaled[mask] = 0

    n = min(nfeatures, len(dispersion_scaled))
    # Stable, and on the negated values -- see the note in `_vst_hvg`.
    top_idx = np.argsort(-dispersion_scaled, kind="stable")[:n]
    return top_idx, means, variances, dispersion_scaled


# ------------------------------------------------------------------
# ScaleData  (mirrors R ScaleData())
# ------------------------------------------------------------------

def scale_data(
    seurat,
    features: Optional[list[str]] = None,
    vars_to_regress: Optional[list[str]] = None,
    assay: Optional[str] = None,
    do_scale: bool = True,
    do_center: bool = True,
    scale_max: float = 10.0,
    layer: str = "data",
) -> None:
    """Scale and optionally center expression data.

    Mirrors R's ScaleData(). Stores result in the 'scale.data' layer
    (Assay5) or scale_data slot (Assay v3).

    Parameters
    ----------
    features      : genes to scale (defaults to variable features)
    vars_to_regress : metadata columns to regress out before scaling
    do_scale      : standardize variance to 1
    do_center     : subtract mean
    scale_max     : clip scaled values at this magnitude
    """
    assay_name = assay or seurat.active_assay
    assay_obj = seurat.assays[assay_name]

    from .assay5 import Assay5
    from .assay import Assay

    if isinstance(assay_obj, Assay5):
        all_features = assay_obj._all_feature_names
        feat_idx_map = {f: i for i, f in enumerate(all_features)}
    else:
        all_features = assay_obj._feature_names
        feat_idx_map = {f: i for i, f in enumerate(all_features)}

    # Select features to scale
    if features is None:
        if isinstance(assay_obj, Assay5):
            features = assay_obj.variable_features or all_features
        else:
            features = assay_obj.var_features or all_features

    feat_idx = [feat_idx_map[f] for f in features if f in feat_idx_map]
    features_present = [all_features[i] for i in feat_idx]

    # Get log-normalized data for the selected features (features × cells)
    data = _get_layer(assay_obj, layer)
    if sp.issparse(data) or is_lazy(data):
        # Subset first, densify second -- on a lazy layer the reverse would
        # read the whole store off disk to keep a few thousand rows of it.
        sub = data[feat_idx, :].toarray().astype(float)
    else:
        sub = np.asarray(data)[feat_idx, :].astype(float)

    # Regress out covariates
    if vars_to_regress:
        sub = _regress_out(sub, seurat.meta_data, vars_to_regress, features_present)

    # Center
    if do_center:
        gene_means = sub.mean(axis=1, keepdims=True)
        sub = sub - gene_means

    # Scale (sample SD, ddof=1, matching Seurat's ScaleData)
    if do_scale:
        gene_stds = sub.std(axis=1, ddof=1, keepdims=True)
        gene_stds[gene_stds == 0] = 1.0
        sub = sub / gene_stds

    # Clip
    sub = np.clip(sub, -scale_max, scale_max)

    # Store scaled data
    scaled_sparse = sp.csc_matrix(sub)
    if isinstance(assay_obj, Assay5):
        # Add/replace scale.data layer; update internal cell/feature maps
        assay_obj.set_layer_data("scale.data", scaled_sparse, feature_names=features_present)
        # Store which features are scaled (needed for PCA)
        assay_obj._scaled_features = features_present
    else:
        assay_obj.scale_data = sub
    log_shanuz_command(
        seurat, "ScaleData", assay=assay or seurat.active_assay,
        params={"do_center": do_center, "do_scale": do_scale,
                "n_features": len(features_present)},
    )


def _regress_out(
    data: np.ndarray,
    meta_data: pd.DataFrame,
    vars_to_regress: list[str],
    feature_names: list[str],
) -> np.ndarray:
    """Regress covariates out of each gene's expression using OLS."""
    from sklearn.linear_model import LinearRegression

    covariates = meta_data[vars_to_regress].values.astype(float)
    # Add intercept
    X = np.column_stack([covariates, np.ones(covariates.shape[0])])

    residuals = np.zeros_like(data)
    for i in range(data.shape[0]):
        y = data[i, :]
        try:
            coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
            residuals[i, :] = y - X @ coef
        except Exception:
            residuals[i, :] = y
    return residuals
