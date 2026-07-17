"""SCTransform — regularized negative-binomial normalization.

Mirrors Seurat's SCTransform() / the sctransform package (Hafemeister & Satija,
Genome Biology 2019; Choudhary & Satija, Genome Biology 2022 for the "v2"
flavor). Counts are modelled per gene with a negative-binomial GLM whose mean
depends on cell sequencing depth; the per-gene parameters are regularized across
genes as smooth functions of average expression, and the Pearson residuals of
that model become the normalized (`scale.data`) values.

The result is stored as a new assay (default name "SCT") with:
  * ``scale.data`` — clipped Pearson residuals for the variable features,
  * ``data``       — log1p of the corrected counts (depth-adjusted),
  * ``counts``     — corrected counts,
plus the residual-variance-ranked variable features.

Flavors
-------
``vst_flavor="v2"`` (the default, matching Seurat 5's ``SCTransform``) fixes the
depth slope at ``log(10)``, excludes non-overdispersed genes from the
regularization and gives them a pure-Poisson model, and floors the residual
variance at ``(median non-zero UMI / 5)^2``. ``vst_flavor="v1"`` fits a free
slope and regularizes every gene.

Fidelity
--------
This is a pure-Python/NumPy reimplementation, verified against R's sctransform
0.4.3 on PBMC 3k. It is not bit-identical — R samples its step-1 genes at random
(so R does not reproduce itself across seeds either) — but the model parameters
and residuals track R closely: on that dataset the regularized intercept matches
at Spearman 1.0000, theta at 0.96, residual variance at 0.9997, and 99.7% of the
3,000 variable features agree. See ``tests/test_sctransform_r_fidelity.py``.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import brentq
from scipy.special import digamma, polygamma

from .assay5 import Assay5
from .preprocessing import _regress_out

_SQRT2PI = np.sqrt(2.0 * np.pi)
# R's bandwidths.c cuts the kernel off at sqrt(DELMAX) standard deviations.
_DELMAX = 1000.0
# R scales ksmooth's normal kernel so its quartiles sit at +/-0.25*bandwidth.
_KSMOOTH_SD = 0.3706506


# ---------------------------------------------------------------------------
# Gene summaries
# ---------------------------------------------------------------------------

def _row_gmean(counts_csr: sp.csr_matrix, eps: float = 1.0) -> np.ndarray:
    """Per-gene geometric mean: ``exp(mean(log(x + eps))) - eps``.

    Mirrors sctransform's ``row_gmean``. The regularization is smoothed against
    log10 of *this*, not the arithmetic mean — the two differ substantially on
    sparse counts and pick out different neighbourhoods of genes.
    """
    G, N = counts_csr.shape
    out = np.empty(G)
    log_eps = np.log(eps)
    indptr, data = counts_csr.indptr, counts_csr.data
    for i in range(G):
        row = data[indptr[i]:indptr[i + 1]]
        total = np.log(row + eps).sum() + (N - row.shape[0]) * log_eps
        out[i] = np.exp(total / N) - eps
    return out


def _row_var(counts_csr: sp.csr_matrix) -> np.ndarray:
    """Per-gene sample variance (ddof=1) of a sparse matrix, without densifying."""
    N = counts_csr.shape[1]
    mean = np.asarray(counts_csr.mean(axis=1)).ravel()
    sq = np.asarray(counts_csr.multiply(counts_csr).sum(axis=1)).ravel()
    return (sq - N * mean ** 2) / (N - 1)


# ---------------------------------------------------------------------------
# Bandwidth selection — port of R's bw.SJ / bw.nrd
# ---------------------------------------------------------------------------

def _bw_pair_cnts(x: np.ndarray, nb: int = 1000):
    """Bin pairwise differences, as R's ``band_den_bin``."""
    n = x.shape[0]
    rang = (x.max() - x.min()) * 1.01
    dd = rang / nb
    idx = (x / dd).astype(np.int64)  # C truncates toward zero
    cnt = np.zeros(nb, dtype=np.int64)
    for i in range(1, n):
        diff = np.abs(idx[i] - idx[:i])
        np.add.at(cnt, diff[diff < nb], 1)
    return dd, cnt


def _bw_phi(n: int, d: float, cnt: np.ndarray, h: float, order: int) -> float:
    """R's ``band_phi4_bin`` / ``band_phi6_bin`` (Scott 1992, eq. 6.67/6.69)."""
    delta = (np.arange(cnt.shape[0]) * d / h) ** 2
    keep = delta < _DELMAX
    delta = delta[keep]
    if order == 4:
        term = np.exp(-delta / 2) * (delta * delta - 6 * delta + 3)
        total = 2 * np.sum(term * cnt[keep]) + n * 3
        return total / (n * (n - 1) * h ** 5 * _SQRT2PI)
    term = np.exp(-delta / 2) * (delta ** 3 - 15 * delta ** 2 + 45 * delta - 15)
    total = 2 * np.sum(term * cnt[keep]) - 15 * n
    return total / (n * (n - 1) * h ** 7 * _SQRT2PI)


def _bw_sj(x: np.ndarray, nb: int = 1000) -> float:
    """Sheather-Jones "solve-the-equation" bandwidth — port of R's ``bw.SJ``.

    NumPy/SciPy have no equivalent, and the smoother's bandwidth sets how hard
    the per-gene parameters are pulled toward the trend, so approximating it
    would quietly change every residual. Verified against R to <1e-6 relative.
    """
    x = np.asarray(x, dtype=float)
    n = x.shape[0]
    d, cnt = _bw_pair_cnts(x, nb)

    q75, q25 = np.percentile(x, [75, 25])
    scale = min(np.std(x, ddof=1), (q75 - q25) / 1.349)
    a = 1.24 * scale * n ** (-1 / 7)
    b = 1.23 * scale * n ** (-1 / 9)
    c1 = 1.0 / (2 * np.sqrt(np.pi) * n)

    td = -_bw_phi(n, d, cnt, b, order=6)
    if not np.isfinite(td) or td <= 0:
        raise ValueError("sample is too sparse to find TD")

    alph2 = 1.357 * (_bw_phi(n, d, cnt, a, order=4) / td) ** (1 / 7)

    def f(h):
        return (c1 / _bw_phi(n, d, cnt, alph2 * h ** (5 / 7), order=4)) ** (1 / 5) - h

    hmax = 1.144 * scale * n ** (-1 / 5)
    lower, upper = 0.1 * hmax, hmax
    for _ in range(100):
        if f(lower) * f(upper) <= 0:
            break
        if f(lower) > 0:
            upper *= 1.2
        else:
            lower /= 1.2
    else:
        raise ValueError("no solution in the specified range of bandwidths")
    # R stops uniroot at tol=0.1*lower, so its answer is only defined to a few
    # tenths of a percent and depends on the path the bracket took. We solve
    # tightly instead: a reproducible root, well inside R's own tolerance.
    return float(brentq(f, lower, upper, xtol=1e-12, rtol=1e-12))


def _bw_nrd(x: np.ndarray) -> float:
    """R's ``bw.nrd`` (Scott's rule), used to weight the step-1 gene sample.

    Note the 1.34 divisor: ``bw.nrd`` and ``bw.nrd0`` scale the IQR by 1.34 while
    ``bw.SJ`` uses 1.349. The two only diverge when the IQR term wins the min, so
    a normal-ish sample will not show the difference.
    """
    n = x.shape[0]
    q75, q25 = np.percentile(x, [75, 25])
    lo = min(np.std(x, ddof=1), (q75 - q25) / 1.34)
    return 1.06 * lo * n ** (-1 / 5)


# ---------------------------------------------------------------------------
# Smoothing and outlier detection
# ---------------------------------------------------------------------------

def _ksmooth_normal(x: np.ndarray, y: np.ndarray, x_points: np.ndarray,
                    bandwidth: float) -> np.ndarray:
    """Nadaraya-Watson smoother with a normal kernel — R's ``ksmooth``."""
    bw = _KSMOOTH_SD * bandwidth
    cutoff = 4.0 * bw
    order = np.argsort(x)
    xs, ys = np.asarray(x, float)[order], np.asarray(y, float)[order]
    out = np.empty(len(x_points))
    for i, x0 in enumerate(np.asarray(x_points, float)):
        lo = np.searchsorted(xs, x0 - cutoff, "left")
        hi = np.searchsorted(xs, x0 + cutoff, "right")
        if hi <= lo:
            out[i] = np.nan
            continue
        w = np.exp(-0.5 * ((xs[lo:hi] - x0) / bw) ** 2)
        sw = w.sum()
        out[i] = np.dot(w, ys[lo:hi]) / sw if sw > 0 else np.nan
    return out


def _robust_scale(x: np.ndarray) -> np.ndarray:
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826
    return (x - med) / (mad + np.finfo(float).eps)


def _robust_scale_binned(y: np.ndarray, x: np.ndarray, breaks: np.ndarray) -> np.ndarray:
    bins = np.digitize(x, breaks)
    score = np.zeros(x.shape[0])
    for b in np.unique(bins):
        m = bins == b
        score[m] = _robust_scale(y[m])
    return score


def _is_outlier(y: np.ndarray, x: np.ndarray, th: float = 10.0) -> np.ndarray:
    """Genes whose fitted parameter is wild for their expression level.

    Port of ``sctransform:::is_outlier``: robust z-score within overlapping bins
    of gene mean, so a gene must look extreme under both binnings to count.
    These are dropped before smoothing so they cannot drag the trend.
    """
    finite = np.isfinite(y)
    if not finite.all():
        out = np.ones(y.shape[0], dtype=bool)
        if finite.sum() > 1:
            out[finite] = _is_outlier(y[finite], x[finite], th)
        return out
    bin_width = (x.max() - x.min()) * _bw_sj(x) / 2
    eps = np.finfo(float).eps * 10
    breaks1 = np.arange(x.min() - eps, x.max() + bin_width, bin_width)
    breaks2 = np.arange(x.min() - eps - bin_width / 2, x.max() + bin_width, bin_width)
    return np.minimum(np.abs(_robust_scale_binned(y, x, breaks1)),
                      np.abs(_robust_scale_binned(y, x, breaks2))) > th


# ---------------------------------------------------------------------------
# Per-gene model fitting
# ---------------------------------------------------------------------------

def _theta_ml(Y: np.ndarray, mu: np.ndarray, limit: int = 10,
              eps: float = np.finfo(float).eps ** 0.25) -> np.ndarray:
    """Maximum-likelihood NB dispersion given fitted means — ``MASS::theta.ml``.

    Vectorised over genes (rows). The moment estimator this replaced ranked
    genes almost independently of the MLE (Spearman 0.18 against R on PBMC 3k),
    which inverted the regularized theta and flattened every residual.
    """
    n = Y.shape[1]
    t0 = n / np.sum((Y / mu - 1) ** 2, axis=1)
    active = np.ones(t0.shape[0], dtype=bool)
    for _ in range(limit):
        t0 = np.abs(t0)
        tt = t0[:, None]
        score = np.sum(digamma(tt + Y) - digamma(tt) + np.log(tt) + 1
                       - np.log(tt + mu) - (Y + tt) / (mu + tt), axis=1)
        info = np.sum(-polygamma(1, tt + Y) + polygamma(1, tt) - 1 / tt
                      + 2 / (mu + tt) - (Y + tt) / (mu + tt) ** 2, axis=1)
        delta = np.divide(score, info, out=np.zeros_like(score), where=info != 0)
        t0 = np.where(active, t0 + delta, t0)
        active &= np.abs(delta) > eps
        if not active.any():
            break
    return np.maximum(t0, 0.0)


def _fit_poisson(counts_sub: np.ndarray, log10_umi: np.ndarray, gene_chunk: int,
                 n_iter: int = 100, tol: float = 1e-9):
    """v1 step 1: per-gene Poisson GLM ``log(mu) = b0 + b1*log10(umi)``, then theta.ml.

    ``counts_sub`` is genes x (subsampled) cells. The design matrix [1, log10_umi]
    is shared across genes, so IRLS is vectorised within gene chunks.
    """
    G = counts_sub.shape[0]
    x, x2 = log10_umi, log10_umi * log10_umi
    b0 = np.zeros(G)
    b1 = np.ones(G)
    theta = np.full(G, 10.0)

    for start in range(0, G, gene_chunk):
        end = min(start + gene_chunk, G)
        Y = counts_sub[start:end]
        bb0 = np.log(Y.mean(axis=1) + 1e-9)
        bb1 = np.ones(Y.shape[0])

        for _ in range(n_iter):
            eta = np.clip(bb0[:, None] + bb1[:, None] * x[None, :], -30, 30)
            mu = np.exp(eta)
            wz = mu * eta + (Y - mu)  # W*z with W = mu (Poisson)
            s00 = mu.sum(axis=1)
            s01 = (mu * x[None, :]).sum(axis=1)
            s11 = (mu * x2[None, :]).sum(axis=1)
            t0_ = wz.sum(axis=1)
            t1_ = (wz * x[None, :]).sum(axis=1)
            det = s00 * s11 - s01 * s01
            det = np.where(np.abs(det) < 1e-12, 1e-12, det)
            nb0 = (s11 * t0_ - s01 * t1_) / det
            nb1 = (s00 * t1_ - s01 * t0_) / det
            done = (np.max(np.abs(nb0 - bb0)) < tol and np.max(np.abs(nb1 - bb1)) < tol)
            bb0, bb1 = nb0, nb1
            if done:
                break

        mu = np.exp(np.clip(bb0[:, None] + bb1[:, None] * x[None, :], -30, 30))
        b0[start:end] = bb0
        b1[start:end] = bb1
        theta[start:end] = _theta_ml(Y, mu)

    return b0, b1, np.maximum(theta, 1e-10)


def _fit_nb_offset(counts_sub: np.ndarray, log10_umi: np.ndarray, gene_chunk: int,
                   n_outer: int = 10, n_inner: int = 25, tol: float = 1e-8):
    """v2 step 1: NB GLM with the depth slope fixed — ``glm.nb(y ~ 1 + offset(log_umi))``.

    Only the intercept and theta are free; the slope is pinned at ``log(10)`` on
    the log10-UMI scale, i.e. mu is proportional to sequencing depth. Alternates
    IRLS for the intercept with ``theta.ml``, as R's ``glm.nb`` does.
    """
    G = counts_sub.shape[0]
    offset = np.log(10.0) * log10_umi  # log(umi)
    b0 = np.zeros(G)
    theta = np.full(G, np.inf)

    for start in range(0, G, gene_chunk):
        end = min(start + gene_chunk, G)
        Y = counts_sub[start:end]

        # Poisson start: closed form for an intercept-only log-link model.
        bb = np.log(np.maximum(Y.sum(axis=1), 1e-9)) - np.log(np.exp(offset).sum())
        mu = np.exp(np.clip(bb[:, None] + offset[None, :], -30, 30))
        th = _theta_ml(Y, mu)

        for _ in range(n_outer):
            for _ in range(n_inner):
                eta = np.clip(bb[:, None] + offset[None, :], -30, 30)
                mu = np.exp(eta)
                w = mu / (1.0 + mu / np.maximum(th, 1e-10)[:, None])
                z = eta + (Y - mu) / mu
                nb = np.sum(w * (z - offset[None, :]), axis=1) / np.sum(w, axis=1)
                done = np.max(np.abs(nb - bb)) < tol
                bb = nb
                if done:
                    break
            mu = np.exp(np.clip(bb[:, None] + offset[None, :], -30, 30))
            th_new = _theta_ml(Y, mu)
            if np.max(np.abs(th_new - th)) < tol * np.maximum(1.0, np.max(th)):
                th = th_new
                break
            th = th_new

        # R caps theta at mean(y)/1e-4 rather than letting it run to infinity.
        th = np.minimum(th, Y.mean(axis=1) / 1e-4)
        b0[start:end] = bb
        theta[start:end] = th

    return b0, np.maximum(theta, 1e-10)


# ---------------------------------------------------------------------------
# Regularization
# ---------------------------------------------------------------------------

def _dispersion_par(log10_gmean: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """The quantity that gets smoothed — R's ``theta_regularization="od_factor"``.

    ``log10(1 + gmean/theta)``, the overdispersion factor, rather than
    ``log10(theta)`` (R's other option, and shanuz's old behaviour). The od
    factor is far closer to linear in expression, so the smoother tracks it
    instead of fighting it; smoothing log(theta) against the wrong x-axis is what
    left the regularized theta anti-correlated with R's.
    """
    return np.log10(1.0 + 10.0 ** log10_gmean / theta)


def _theta_from_dispersion_par(log10_gmean: np.ndarray, disp: np.ndarray) -> np.ndarray:
    """Invert :func:`_dispersion_par` once it has been smoothed."""
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        return 10.0 ** log10_gmean / (10.0 ** disp - 1.0)


def _reg_model_pars(model_pars: np.ndarray, log10_gmean_step1: np.ndarray,
                    log10_gmean: np.ndarray, bw_adjust: float, verbose: bool):
    """Smooth per-gene parameters against log10 geometric mean.

    ``model_pars`` is (step1 genes x params) with the dispersion parameter in
    column 0. Returns the smoothed parameters for every gene.
    """
    outliers = np.zeros(model_pars.shape[0], dtype=bool)
    for j in range(model_pars.shape[1]):
        outliers |= _is_outlier(model_pars[:, j], log10_gmean_step1)
    if outliers.any():
        if verbose:
            print(f"  {outliers.sum()} outliers ignored during regularization")
        model_pars = model_pars[~outliers]
        log10_gmean_step1 = log10_gmean_step1[~outliers]

    bw = _bw_sj(log10_gmean_step1) * bw_adjust
    x_points = np.clip(log10_gmean, log10_gmean_step1.min(), log10_gmean_step1.max())
    fit = np.empty((log10_gmean.shape[0], model_pars.shape[1]))
    for j in range(model_pars.shape[1]):
        fit[:, j] = _ksmooth_normal(log10_gmean_step1, model_pars[:, j], x_points, bw)
    return fit


# ---------------------------------------------------------------------------
# SCTransform
# ---------------------------------------------------------------------------

def sctransform(
    seurat,
    assay: Optional[str] = None,
    new_assay_name: str = "SCT",
    n_cells: int = 5000,
    n_genes: int = 2000,
    n_features: int = 3000,
    min_cells: int = 5,
    vars_to_regress: Optional[list[str]] = None,
    clip_range: Optional[tuple] = None,
    gene_chunk: int = 500,
    seed: int = 42,
    set_default: bool = True,
    vst_flavor: str = "v2",
    bw_adjust: float = 3.0,
    verbose: bool = False,
) -> "object":
    """Run SCTransform and attach a normalized assay.

    Mirrors R's ``SCTransform(object)``. Fits the regularized NB model on the
    active assay's counts and stores the result as ``new_assay_name`` ("SCT").

    Parameters
    ----------
    n_cells     : cells subsampled for parameter estimation (Seurat default 5000).
    n_genes     : genes used for step-1 estimation, sampled to spread evenly over
                  expression (Seurat default 2000). ``None`` uses every gene.
    n_features  : number of variable features by residual variance (default 3000).
    min_cells   : drop genes detected in fewer than this many cells (default 5).
    vars_to_regress : metadata columns (e.g. 'percent.mt') regressed out of the
                  Pearson residuals, mirroring SCTransform's vars.to.regress.
    clip_range  : residual clip for ``scale.data``; default (-√(N/30), √(N/30)).
                  Note this is *not* the clip used when ranking variable features
                  — see below.
    vst_flavor  : "v2" (default, as Seurat 5) or "v1". See the module docstring.
    bw_adjust   : multiplier on the Sheather-Jones smoothing bandwidth (R's 3).
    set_default : make the new assay the active assay.

    Returns ``seurat`` with the new assay added.
    """
    if vst_flavor not in ("v1", "v2"):
        raise ValueError(f"vst_flavor must be 'v1' or 'v2', got {vst_flavor!r}")

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
    counts_csr = counts.tocsr()
    if verbose:
        print(f"  SCTransform ({vst_flavor}): {G}/{G_all} genes kept "
              f"(>= {min_cells} cells), {N} cells")

    total_umi = np.asarray(counts.sum(axis=0)).ravel().astype(float)
    total_umi[total_umi == 0] = 1.0
    log10_umi = np.log10(total_umi)
    log10_gmean = np.log10(_row_gmean(counts_csr, eps=1.0))

    rng = np.random.default_rng(seed)

    # ---- step 1: which cells and genes estimate the model ----
    if N > n_cells:
        cells_step1 = np.sort(rng.choice(N, n_cells, replace=False))
        det = np.diff(counts_csr[:, cells_step1].tocsr().indptr)
        genes_step1 = np.where(det >= min_cells)[0]
    else:
        cells_step1 = np.arange(N)
        genes_step1 = np.arange(G)

    gene_amean = np.asarray(counts_csr.mean(axis=1)).ravel()
    gene_var = _row_var(counts_csr)
    overdispersion = gene_var - gene_amean
    poisson_genes = np.zeros(G, dtype=bool)

    if vst_flavor == "v2":
        # Genes whose variance does not exceed their mean carry no NB signal;
        # regularizing them drags the trend, so R models them as pure Poisson.
        poisson_genes = (overdispersion <= 0) | (gene_amean < 0.001)
        genes_step1 = genes_step1[overdispersion[genes_step1] > 0]
        if verbose:
            print(f"  {poisson_genes.sum()} poisson genes excluded from regularization")

    if n_genes is not None and n_genes < len(genes_step1):
        # Sample step-1 genes inversely to their density in log-gmean, so the
        # smoother sees the sparse tails, not just the crowded middle.
        x = log10_gmean[genes_step1]
        bw = _bw_nrd(x)
        dens = np.exp(-0.5 * ((x[:, None] - x[None, :]) / bw) ** 2).sum(axis=1)
        dens /= (len(x) * bw * _SQRT2PI)
        prob = 1.0 / (dens + np.finfo(float).eps)
        genes_step1 = rng.choice(genes_step1, size=n_genes, replace=False,
                                 p=prob / prob.sum())
        genes_step1 = np.sort(genes_step1)

    log10_gmean_step1 = log10_gmean[genes_step1]
    Y1 = counts_csr[genes_step1][:, cells_step1].toarray().astype(float)
    log10_umi_step1 = log10_umi[cells_step1]

    if verbose:
        print(f"  fitting NB model on {len(cells_step1)} cells x "
              f"{len(genes_step1)} genes ...")

    # ---- step 1: fit ----
    if vst_flavor == "v2":
        b0_s1, theta_s1 = _fit_nb_offset(Y1, log10_umi_step1, gene_chunk)
        b1_s1 = np.full(len(genes_step1), np.log(10.0))
    else:
        b0_s1, b1_s1, theta_s1 = _fit_poisson(Y1, log10_umi_step1, gene_chunk)
    del Y1

    # ---- step 2: regularize ----
    model_pars = np.column_stack([_dispersion_par(log10_gmean_step1, theta_s1),
                                  b0_s1, b1_s1])
    fit = _reg_model_pars(model_pars, log10_gmean_step1, log10_gmean, bw_adjust, verbose)

    theta_r = _theta_from_dispersion_par(log10_gmean, fit[:, 0])
    b0r, b1r = fit[:, 1], fit[:, 2]

    if vst_flavor == "v2":
        mean_cell_sum = total_umi.mean()
        with np.errstate(divide="ignore"):
            b0r = np.where(poisson_genes,
                           np.log(np.maximum(gene_amean, 1e-300)) - np.log(mean_cell_sum),
                           b0r)
        theta_r = np.where(poisson_genes, np.inf, theta_r)
        b1r = np.full(G, np.log(10.0))  # fix_slope

    bad = ~np.isfinite(theta_r) & ~poisson_genes
    theta_r = np.where(bad, 1e6, theta_r)
    theta_r = np.where(poisson_genes, np.inf, np.clip(theta_r, 1e-2, 1e6))

    # ---- step 3: residuals ----
    # Two different clips, as in R. Residual *variance* — which ranks the
    # variable features — comes from residuals clipped at sqrt(N); only the
    # stored scale.data is clipped to the much tighter sqrt(N/30). Using the
    # tight clip for both crushes exactly the marker genes that define rare
    # subsets, so the ranking loses them.
    res_clip = np.sqrt(N)
    if clip_range is None:
        c = np.sqrt(N / 30.0)
        clip_lo, clip_hi = -c, c
    else:
        clip_lo, clip_hi = clip_range

    min_var = -np.inf
    if vst_flavor == "v2":
        # R: (median non-zero UMI / 5)^2 — a floor on model variance that stops
        # near-zero fitted means from manufacturing enormous residuals.
        min_var = (np.median(counts_csr.data) / 5.0) ** 2
        if verbose:
            print(f"  min_variance = {min_var:.6g}")

    median_log10_umi = float(np.median(log10_umi))
    res_var = np.zeros(G)
    corrected_blocks = []

    for start in range(0, G, gene_chunk):
        end = min(start + gene_chunk, G)
        y = counts_csr[start:end].toarray().astype(float)
        mu = np.exp(np.clip(b0r[start:end, None] + b1r[start:end, None] * log10_umi[None, :],
                            -30, 30))
        var = mu + mu * mu / theta_r[start:end, None]
        var = np.maximum(var, min_var)
        z = (y - mu) / np.sqrt(var)
        res_var[start:end] = np.clip(z, -res_clip, res_clip).var(axis=1, ddof=1)

        # Corrected counts: the residual re-expressed at the median depth. R
        # uses the unclipped residual and no variance floor here.
        mu_med = np.exp(b0r[start:end, None] + b1r[start:end, None] * median_log10_umi)
        var_med = mu_med + mu_med * mu_med / theta_r[start:end, None]
        z_raw = (y - mu) / np.sqrt(mu + mu * mu / theta_r[start:end, None])
        corr = np.clip(np.round(mu_med + z_raw * np.sqrt(var_med)), 0.0, None)
        corrected_blocks.append(sp.csr_matrix(corr))

    corrected = sp.vstack(corrected_blocks, format="csc")
    del corrected_blocks

    # ---- variable features by residual variance ----
    n_feat = min(n_features, G)
    top = np.argsort(res_var)[::-1][:n_feat]
    var_features = [genes[i] for i in top]

    # ---- scale.data: residuals for the variable features only ----
    top_sorted = np.sort(top)
    y = counts_csr[top_sorted].toarray().astype(float)
    mu = np.exp(np.clip(b0r[top_sorted, None] + b1r[top_sorted, None] * log10_umi[None, :],
                        -30, 30))
    var = np.maximum(mu + mu * mu / theta_r[top_sorted, None], min_var)
    resid = np.clip((y - mu) / np.sqrt(var), clip_lo, clip_hi)
    scale_feats = [genes[i] for i in top_sorted]

    if vars_to_regress:
        resid = _regress_out(resid, seurat.meta_data, list(vars_to_regress), scale_feats)

    sct = Assay5(
        layers={"counts": corrected, "data": _log1p_sparse(corrected)},
        feature_names=list(genes),
        cell_names=list(cell_names),
        key=f"{new_assay_name.lower()}_",
    )
    sct.set_layer_data("scale.data", sp.csc_matrix(resid), feature_names=scale_feats)
    sct._scaled_features = scale_feats
    sct.variable_features = var_features
    sct.meta_data["residual_variance"] = pd.Series(res_var, index=genes)
    sct.meta_data["theta"] = pd.Series(theta_r, index=genes)
    sct.meta_data["gmean"] = pd.Series(10.0 ** log10_gmean, index=genes)

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
