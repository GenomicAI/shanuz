"""SCTransform — regularized negative-binomial normalization.

Mirrors Seurat's SCTransform() / the sctransform package (Hafemeister & Satija,
Genome Biology 2019). Counts are modelled per gene with a negative-binomial GLM
whose mean depends on cell sequencing depth; the per-gene parameters are
regularized across genes as smooth functions of average expression, and the
Pearson residuals of that model become the normalized (`scale.data`) values.

The result is stored as a new assay (default name "SCT") with:
  * ``scale.data`` — clipped Pearson residuals for the variable features,
  * ``data``       — log1p of the corrected counts (depth-adjusted),
  * ``counts``     — corrected counts,
plus the residual-variance-ranked variable features.

This is a faithful pure-Python/NumPy reimplementation of the algorithm. It is
not bit-identical to the C++/R sctransform (LOESS/theta.ml differ), but follows
the same model and produces equivalent residuals and variable-feature ranking.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .assay5 import Assay5
from .preprocessing import _loess_fit


# ---------------------------------------------------------------------------
# Per-gene model fitting
# ---------------------------------------------------------------------------

def _fit_poisson_theta(
    counts_sub: np.ndarray,
    log10_umi: np.ndarray,
    gene_chunk: int,
    n_iter: int = 25,
    tol: float = 1e-4,
):
    """Fit, for every gene, a Poisson GLM  log(mu) = b0 + b1·log10(umi)  and a
    moment estimate of the negative-binomial dispersion theta.

    ``counts_sub`` is genes × (subsampled) cells. The design matrix [1, log10_umi]
    is shared across genes, so IRLS is vectorised within gene chunks.
    """
    G, n = counts_sub.shape
    x = log10_umi
    x2 = x * x
    b0 = np.zeros(G)
    b1 = np.ones(G)
    theta = np.full(G, 10.0)

    for start in range(0, G, gene_chunk):
        end = min(start + gene_chunk, G)
        Y = counts_sub[start:end]                        # (g × n)
        g = Y.shape[0]
        bb0 = np.log(Y.mean(axis=1) + 1e-9)              # init intercept
        bb1 = np.ones(g)

        for _ in range(n_iter):
            eta = np.clip(bb0[:, None] + bb1[:, None] * x[None, :], -30, 30)
            mu = np.exp(eta)
            # Weighted normal equations with shared design columns [1, x].
            Wz = mu * eta + (Y - mu)                       # W·z, W = mu (Poisson)
            S00 = mu.sum(axis=1)
            S01 = (mu * x[None, :]).sum(axis=1)
            S11 = (mu * x2[None, :]).sum(axis=1)
            T0 = Wz.sum(axis=1)
            T1 = (Wz * x[None, :]).sum(axis=1)
            det = S00 * S11 - S01 * S01
            det = np.where(np.abs(det) < 1e-12, 1e-12, det)
            nb0 = (S11 * T0 - S01 * T1) / det
            nb1 = (S00 * T1 - S01 * T0) / det
            if np.max(np.abs(nb0 - bb0)) < tol and np.max(np.abs(nb1 - bb1)) < tol:
                bb0, bb1 = nb0, nb1
                break
            bb0, bb1 = nb0, nb1

        eta = np.clip(bb0[:, None] + bb1[:, None] * x[None, :], -30, 30)
        mu = np.exp(eta)
        # Moment estimate of theta:  E[(y-mu)^2] = mu + mu^2/theta.
        num = (mu * mu).sum(axis=1)
        den = ((Y - mu) ** 2 - mu).sum(axis=1)
        th = np.where(den > 0, num / den, 1e4)
        th = np.clip(th, 1e-2, 1e4)

        b0[start:end] = bb0
        b1[start:end] = bb1
        theta[start:end] = th

    return b0, b1, theta


def _regularize(log10_gmean: np.ndarray, values: np.ndarray, frac: float = 0.3) -> np.ndarray:
    """Smooth per-gene parameters as a function of log10 mean expression (LOESS)."""
    return _loess_fit(log10_gmean, values, frac=frac)


# ---------------------------------------------------------------------------
# SCTransform
# ---------------------------------------------------------------------------

def sctransform(
    seurat,
    assay: Optional[str] = None,
    new_assay_name: str = "SCT",
    n_cells: int = 5000,
    n_features: int = 3000,
    min_cells: int = 5,
    clip_range: Optional[tuple] = None,
    gene_chunk: int = 500,
    seed: int = 42,
    set_default: bool = True,
    verbose: bool = False,
) -> "object":
    """Run SCTransform and attach a normalized assay.

    Mirrors R's ``SCTransform(object)``. Fits the regularized NB model on the
    active assay's counts and stores the result as ``new_assay_name`` ("SCT").

    Parameters
    ----------
    n_cells     : cells subsampled for parameter estimation (Seurat default 5000).
    n_features  : number of variable features by residual variance (default 3000).
    min_cells   : drop genes detected in fewer than this many cells (default 5).
    clip_range  : residual clip (low, high); default (−√(N/30), √(N/30)).
    gene_chunk  : genes processed per vectorised batch (memory control).
    set_default : make the new assay the active assay.

    Returns ``seurat`` with the new assay added.
    """
    src = seurat.assays[assay or seurat.active_assay]
    if isinstance(src, Assay5):
        counts = src.layers.get("counts")
        all_genes = src._all_feature_names
    else:
        counts = src.counts
        all_genes = src._feature_names
    if counts is None:
        raise ValueError("SCTransform requires a counts layer.")

    counts = counts.tocsc() if sp.issparse(counts) else sp.csc_matrix(counts)
    cell_names = seurat.cell_names()
    G_all, N = counts.shape

    # Drop genes detected in too few cells.
    nnz_per_gene = np.diff(counts.tocsr().indptr)
    keep = np.where(nnz_per_gene >= min_cells)[0]
    counts = counts[keep, :]
    genes = [all_genes[i] for i in keep]
    G = len(genes)
    if verbose:
        print(f"  SCTransform: {G}/{G_all} genes kept (>= {min_cells} cells), {N} cells")

    total_umi = np.asarray(counts.sum(axis=0)).flatten().astype(float)
    total_umi[total_umi == 0] = 1.0
    log10_umi = np.log10(total_umi)
    gene_mean = np.asarray(counts.mean(axis=1)).flatten()
    log10_gmean = np.log10(gene_mean)

    # Subsample cells for the (expensive) parameter estimation.
    rng = np.random.default_rng(seed)
    if N > n_cells:
        sub = np.sort(rng.choice(N, n_cells, replace=False))
    else:
        sub = np.arange(N)
    counts_sub = counts[:, sub].toarray().astype(float)
    log10_umi_sub = log10_umi[sub]

    if verbose:
        print(f"  fitting NB model on {len(sub)} cells ...")
    b0, b1, theta = _fit_poisson_theta(counts_sub, log10_umi_sub, gene_chunk)
    del counts_sub

    # Regularize parameters across genes vs. mean expression.
    b0r = _regularize(log10_gmean, b0)
    b1r = _regularize(log10_gmean, b1)
    theta_r = 10.0 ** _regularize(log10_gmean, np.log10(theta))
    theta_r = np.clip(theta_r, 1e-2, 1e6)

    if clip_range is None:
        c = np.sqrt(N / 30.0)
        clip_lo, clip_hi = -c, c
    else:
        clip_lo, clip_hi = clip_range

    median_log10_umi = float(np.median(log10_umi))

    # Pass A: stream gene chunks over all cells → residual variance + corrected
    # counts. Residuals themselves are discarded to bound memory.
    res_var = np.zeros(G)
    corrected_blocks = []
    counts_csr = counts.tocsr()
    for start in range(0, G, gene_chunk):
        end = min(start + gene_chunk, G)
        y = counts_csr[start:end].toarray().astype(float)
        eta = np.clip(b0r[start:end, None] + b1r[start:end, None] * log10_umi[None, :], -30, 30)
        mu = np.exp(eta)
        sd = np.sqrt(mu + mu * mu / theta_r[start:end, None])
        z = np.clip((y - mu) / sd, clip_lo, clip_hi)
        res_var[start:end] = z.var(axis=1, ddof=1)

        mu_med = np.exp(b0r[start:end, None] + b1r[start:end, None] * median_log10_umi)
        sd_med = np.sqrt(mu_med + mu_med * mu_med / theta_r[start:end, None])
        corr = np.clip(mu_med + z * sd_med, 0.0, None)
        corrected_blocks.append(sp.csr_matrix(np.round(corr)))
    corrected = sp.vstack(corrected_blocks, format="csc")
    del corrected_blocks

    # Variable features by residual variance (descending).
    n_feat = min(n_features, G)
    top = np.argsort(res_var)[::-1][:n_feat]
    var_features = [genes[i] for i in top]

    # Pass B: residuals (scale.data) for the variable features only.
    top_sorted = np.sort(top)
    y = counts_csr[top_sorted].toarray().astype(float)
    eta = np.clip(b0r[top_sorted, None] + b1r[top_sorted, None] * log10_umi[None, :], -30, 30)
    mu = np.exp(eta)
    sd = np.sqrt(mu + mu * mu / theta_r[top_sorted, None])
    resid = np.clip((y - mu) / sd, clip_lo, clip_hi)
    scale_feats = [genes[i] for i in top_sorted]

    # Build the SCT assay: corrected counts + log1p(data); residual scale.data.
    sct = Assay5(
        layers={
            "counts": corrected,
            "data": _log1p_sparse(corrected),
        },
        feature_names=list(genes),
        cell_names=list(cell_names),
        key=f"{new_assay_name.lower()}_",
    )
    sct.set_layer_data("scale.data", sp.csc_matrix(resid), feature_names=scale_feats)
    sct._scaled_features = scale_feats
    sct.variable_features = var_features
    # Stash residual variance / model parameters in the assay meta.
    sct.meta_data["residual_variance"] = pd.Series(res_var, index=genes)
    sct.meta_data["theta"] = pd.Series(theta_r, index=genes)

    seurat.assays[new_assay_name] = sct
    if set_default:
        seurat.active_assay = new_assay_name
    if verbose:
        print(f"  SCT assay '{new_assay_name}' added: {G} genes, "
              f"{n_feat} variable features")
    return seurat


def _log1p_sparse(mat: sp.spmatrix) -> sp.csc_matrix:
    """log1p of a sparse matrix (zeros stay zero)."""
    out = mat.tocsc(copy=True).astype(float)
    out.data = np.log1p(out.data)
    return out
