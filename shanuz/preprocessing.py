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


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

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
        if sp.issparse(mat):
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


def _log_normalize(counts, scale_factor: float = 10000.0):
    """Normalize counts per cell to scale_factor, then log1p."""
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

    if sp.issparse(data):
        csr = data.tocsr()
        means = np.array(csr.mean(axis=1)).flatten()
        mean_sq = np.array(csr.power(2).mean(axis=1)).flatten()
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
    if sp.issparse(data):
        indptr = csr.indptr
        values = csr.data.astype(float)
        for g in valid_idx:
            sd = expected_sd[g]
            if sd <= 0:
                continue
            mu = means[g]
            row_vals = values[indptr[g]:indptr[g + 1]]
            n_zero = n_cells - row_vals.size
            # Zero entries standardize to -mu/sd (negative → never clipped).
            zero_term = n_zero * (mu / sd) ** 2
            z = (row_vals - mu) / sd
            np.minimum(z, clip_max, out=z)
            var_standardized[g] = (zero_term + np.dot(z, z)) / (n_cells - 1)
    else:
        for g in valid_idx:
            sd = expected_sd[g]
            if sd <= 0:
                continue
            z = (d[g, :] - means[g]) / sd
            np.minimum(z, clip_max, out=z)
            var_standardized[g] = np.dot(z, z) / (n_cells - 1)

    top_idx = np.argsort(var_standardized)[::-1][:nfeatures]

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
    """
    n = len(x)
    k = max(int(np.ceil(frac * n)), 6)
    half_k = k // 2

    sort_idx = np.argsort(x)
    xs = x[sort_idx]
    ys = y[sort_idx]
    fitted_s = np.empty(n)

    for b0 in range(0, n, batch_size):
        b1 = min(b0 + batch_size, n)
        nb = b1 - b0
        b_i = np.arange(b0, b1)

        # Sliding window: exactly k neighbours for each query point
        lefts = np.maximum(0, b_i - half_k)
        rights = np.minimum(n, lefts + k)
        lefts = np.maximum(0, rights - k)

        win_idx = lefts[:, np.newaxis] + np.arange(k)   # (nb, k)
        xw = xs[win_idx]   # (nb, k)
        yw = ys[win_idx]   # (nb, k)
        xi = xs[b_i, np.newaxis]   # (nb, 1)

        # Tricube weights
        d = np.abs(xw - xi)
        dmax = np.maximum(d.max(axis=1, keepdims=True), 1e-10)
        u = d / dmax
        w = np.maximum(1.0 - u ** 3, 0.0) ** 3   # (nb, k)

        # Local quadratic design matrix: columns = [1, (x-x₀), (x-x₀)²]
        dx = xw - xi   # (nb, k)
        A = np.stack([np.ones((nb, k)), dx, dx ** 2], axis=2)   # (nb, k, 3)

        # Weighted normal equations: (AᵀWA)β = AᵀWy
        Aw = A * w[:, :, np.newaxis]   # (nb, k, 3)
        AtWA = np.einsum("nki,nkj->nij", Aw, A)   # (nb, 3, 3)
        AtWy = np.einsum("nki,nk->ni", Aw, yw)    # (nb, 3)

        try:
            # solve needs RHS as (..., M, 1), returns (..., M, 1)
            betas = np.linalg.solve(AtWA, AtWy[:, :, np.newaxis])[:, :, 0]  # (nb, 3)
            fitted_s[b0:b1] = betas[:, 0]
        except np.linalg.LinAlgError:
            for j in range(nb):
                try:
                    sol, *_ = np.linalg.lstsq(AtWA[j], AtWy[j], rcond=None)
                    fitted_s[b0 + j] = sol[0]
                except Exception:
                    fitted_s[b0 + j] = float(np.average(yw[j], weights=w[j]))

    fitted = np.empty(n)
    fitted[sort_idx] = fitted_s
    return fitted


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
    top_idx = np.argsort(dispersion_scaled)[::-1][:n]
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
    if sp.issparse(data):
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
