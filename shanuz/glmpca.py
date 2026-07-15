"""GLM-PCA — dimensionality reduction on raw counts.

Mirrors R's ``RunGLMPCA`` (SeuratWrappers → the `glmpca` package), following
Townes et al. (2019), *Feature selection and dimension reduction for single-cell
RNA-seq based on a multinomial model*.

The usual pipeline log-normalises counts and then runs PCA, which quietly assumes
the transformed data is Gaussian with constant variance. Counts are not: a gene
averaging 0.1 UMIs and one averaging 100 have wildly different noise, the log
transform does not fix that, and the pseudocount you add to survive ``log(0)``
distorts exactly the low-expression genes where most of the zeros live. GLM-PCA
drops the transformation entirely and fits a low-rank model on the count scale::

    Y[g,c] ~ Poisson(μ[g,c])
    log μ[g,c] = a[g] + o[c] + Σ_l U[g,l]·V[c,l]

``a`` is a per-gene intercept, ``o`` a fixed per-cell offset (log library size, so
sequencing depth is a known quantity rather than a factor to be discovered), and
``U``/``V`` are the rank-L loadings and factors — the counterparts of PCA's
loadings and embeddings. There is no closed form, so the model is fitted by
Fisher scoring.

Two noise models are available. ``family="poisson"`` assumes the variance equals
the mean; real UMI counts are usually noisier than that (the same gene in two
copies of the same cell state still varies more than Poisson predicts), and a
handful of over-dispersed genes will otherwise dominate the fit. ``family="nb"``
swaps in a negative binomial::

    Y[g,c] ~ NB(μ[g,c], θ)      Var = μ + μ²/θ

with a single shared dispersion ``θ``: as ``θ → ∞`` the extra ``μ²/θ`` term
vanishes and NB collapses back onto Poisson, so NB never fits *worse* than
Poisson, only more forgivingly. ``θ`` can be pinned or estimated from the data by
maximum likelihood alongside the factors (see ``glm_pca``'s ``theta`` /
``optimize_theta``).
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import scipy.sparse as sp

from .dimreduc import DimReduc

FAMILIES = ("poisson", "nb")

_THETA_MIN, _THETA_MAX = 1e-2, 1e6      # keep the dispersion positive and finite

_ETA_CLIP = 30.0        # exp(30) ≈ 1e13; keeps a bad step from overflowing to inf
_EPS = 1e-10


def glm_pca(
    seurat,
    n_components: int = 10,
    features: Optional[list[str]] = None,
    assay: Optional[str] = None,
    reduction_name: str = "glmpca",
    reduction_key: str = "GLMPC_",
    family: str = "poisson",
    layer: str = "counts",
    max_iter: int = 100,
    tol: float = 1e-4,
    penalty: float = 1.0,
    learning_rate: float = 0.1,
    theta: float = 100.0,
    optimize_theta: bool = True,
    seed: int = 42,
) -> None:
    """Fit a Poisson GLM-PCA and store it as a DimReduc.

    Mirrors R's ``RunGLMPCA(obj, L = 10)``. Takes **raw counts**, not normalised
    or scaled data — the whole point is to model the counts as counts. Stores
    factors as ``cell_embeddings`` and loadings as ``feature_loadings``, so
    ``find_neighbors(obj, reduction="glmpca")`` and ``run_umap`` work downstream
    exactly as they do off PCA.

    Parameters
    ----------
    n_components  : rank of the fit (L) — the number of factors
    features      : genes to use (defaults to variable features)
    assay         : assay name (defaults to active assay)
    family        : noise model — ``"poisson"`` or ``"nb"`` (negative binomial)
    layer         : layer to read counts from
    max_iter      : maximum Fisher scoring iterations
    tol           : stop when the relative change in deviance falls below this
    penalty       : L2 ridge on U and V. U and V can trade scale freely
                    (``U·Vᵀ = (cU)·(V/c)ᵀ``); the ridge is what pins that down.
    learning_rate : initial Fisher step size. Halved on any step that fails to
                    lower the deviance, so this is an opening bid, not a
                    commitment.
    theta         : negative-binomial dispersion (``Var = μ + μ²/θ``). Ignored for
                    Poisson. With ``optimize_theta`` it is only the starting value;
                    otherwise it is held fixed at this number for the whole fit.
    optimize_theta: re-estimate ``θ`` by maximum likelihood between factor updates
                    (NB only). Turn off to fit at a dispersion you already trust —
                    that also restores strict monotone deviance, since a moving
                    ``θ`` re-scales the deviance under it.
    seed          : only used if the counts have no structure at all — the fit is
                    started deterministically from the data (see `_init_factors`)

    Notes
    -----
    Deviance falls monotonically by construction: a step that raises it (or
    overflows) is rejected outright and retried at half the step size. The full
    trace is kept in ``reduction.misc["deviance"]`` — if it is still dropping
    steeply at the end, raise ``max_iter``. ``misc["converged"]`` says whether the
    fit stopped because it was done or because it ran out of iterations.

    Fitting is dense in genes × cells. Pass a few thousand variable features
    rather than the whole transcriptome, as you would to `run_pca`.
    """
    if family not in FAMILIES:
        raise NotImplementedError(
            f"family={family!r} is not implemented; use one of {list(FAMILIES)}.")

    from .markers import _get_expression_matrix
    from .reduction import _default_features

    assay_name = assay or seurat.active_assay
    assay_obj = seurat.assays[assay_name]

    features = _default_features(assay_obj, features)
    Y = _counts_for(assay_obj, features, layer, _get_expression_matrix)

    n_genes, n_cells = Y.shape
    n_components = min(n_components, min(n_genes, n_cells) - 1)
    if n_components < 1:
        raise ValueError(
            f"Need at least 2 genes and 2 cells to fit GLM-PCA; got "
            f"{n_genes} × {n_cells}.")

    # Library size as a fixed offset: depth is known, not something to rediscover.
    totals = Y.sum(axis=0)
    if np.any(totals <= 0):
        raise ValueError(
            "Some cells have zero total counts, so they have no library size to "
            "offset by. Filter them out first.")
    offsets = np.log(totals / totals.mean())

    if family == "poisson":
        intercept, U, V, deviance, converged = _fit_poisson(
            Y, n_components, offsets,
            max_iter=max_iter, tol=tol, penalty=penalty,
            learning_rate=learning_rate, seed=seed,
        )
        fitted_theta = np.inf
    else:                                                   # family == "nb"
        intercept, U, V, deviance, converged, fitted_theta = _fit_nb(
            Y, n_components, offsets,
            max_iter=max_iter, tol=tol, penalty=penalty,
            learning_rate=learning_rate, theta=theta,
            optimize_theta=optimize_theta, seed=seed,
        )
    loadings, factors = _orthogonalize(U, V)

    seurat.reductions[reduction_name] = DimReduc(
        cell_embeddings=factors,
        feature_loadings=loadings,
        assay_used=assay_name,
        stdev=np.sqrt(np.var(factors, axis=0, ddof=1)),
        key=reduction_key,
        cell_names=seurat.cell_names(),
        feature_names=list(features),
        misc={
            "glmpca_family": family,
            "deviance": deviance,
            "converged": converged,
            "intercept": intercept,
            "theta": fitted_theta,
        },
    )


def _counts_for(assay_obj, features, layer, getter) -> np.ndarray:
    """The raw counts for ``features``, dense, genes × cells."""
    mat, all_features = getter(assay_obj, layer)
    index = {f: i for i, f in enumerate(all_features)}
    missing = [f for f in features if f not in index]
    if missing:
        raise ValueError(
            f"{len(missing)} requested feature(s) are not in the assay, "
            f"e.g. {missing[:3]}.")

    rows = [index[f] for f in features]
    Y = mat[rows, :]
    Y = Y.toarray() if sp.issparse(Y) else np.asarray(Y)
    Y = Y.astype(float)

    if not np.isfinite(Y).all():
        raise ValueError("Counts contain NaN or infinite values.")
    if Y.min() < 0:
        raise ValueError(
            f"GLM-PCA models counts with a Poisson likelihood, but layer "
            f"{layer!r} holds negative values (min {Y.min():.3g}) — this looks "
            f"like normalised or scaled data. Pass layer='counts'.")
    return Y


def _init_factors(Y, intercept, offsets, L, seed, target_spread: float = 0.5):
    """Start the factors from the SVD of the null model's residuals.

    ``U = V = 0`` is an exact saddle of the log-likelihood — the score for each
    block is a product with the *other* block, so at zero both vanish. Start the
    two near zero and the fit inches away from the saddle, the deviance barely
    moves on the first step, and any relative-improvement stopping rule declares
    victory on a model that has not fitted anything. (It still looks plausible:
    one step is enough to orient the factors, so cluster structure shows up in a
    plot while the deviance sits at its null value.)

    So do what Townes does instead: fit the intercept-only model, and let the
    leading singular vectors of what it *fails* to explain say where the
    structure is. Deterministic — ``seed`` is only reached if the residuals have
    no structure at all.
    """
    mu0 = np.exp(np.clip(intercept[:, None] + offsets[None, :],
                         -_ETA_CLIP, _ETA_CLIP))
    residual = (Y - mu0) / np.sqrt(mu0 + _EPS)          # Pearson scale
    residual = np.nan_to_num(residual, nan=0.0, posinf=0.0, neginf=0.0)

    k = min(L, min(residual.shape) - 1)
    if k < min(residual.shape) // 2:
        from scipy.sparse.linalg import svds
        P, S, Qt = svds(residual, k=k, random_state=seed)
        order = np.argsort(S)[::-1]                     # svds returns ascending
        P, S, Qt = P[:, order], S[order], Qt[order]
    else:
        P, S, Qt = np.linalg.svd(residual, full_matrices=False)

    root = np.sqrt(S[:L])
    U = P[:, :L] * root
    V = Qt[:L].T * root

    # Pearson residuals are z-scores. A rank-L slab of them is a wildly
    # overconfident opening bid for a log-fold-change, so damp it to something a
    # log link can survive and let Fisher scoring take it from there.
    spread = float(np.std(U @ V.T))
    if spread > 0:
        U *= np.sqrt(target_spread / spread)
        V *= np.sqrt(target_spread / spread)
    else:                                               # no structure to find
        rng = np.random.default_rng(seed)
        U = rng.normal(scale=1e-2, size=U.shape)
        V = rng.normal(scale=1e-2, size=V.shape)
    return U, V


def _poisson_deviance(Y: np.ndarray, mu: np.ndarray) -> float:
    """2·Σ[y·log(y/μ) − (y − μ)] — the Poisson goodness-of-fit being minimised."""
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(Y > 0, Y * np.log(Y / mu), 0.0)
    return float(2.0 * np.sum(term - (Y - mu)))


def _fit_poisson(
    Y: np.ndarray,
    L: int,
    offsets: np.ndarray,
    max_iter: int,
    tol: float,
    penalty: float,
    learning_rate: float,
    seed: int,
):
    """Fisher scoring for the Poisson log-bilinear model, with backtracking.

    Each block (intercept, loadings, factors) gets a diagonal Newton step: the
    score divided by the Fisher information. Under a log link and Poisson noise
    both are one matrix product each — the score is driven by the raw residual
    ``Y − μ``, and the information by ``μ`` itself.
    """
    # Intercept-only fit: exp(a_g + o_c) reproduces each gene's mean rate. A gene
    # detected in no cell floors at exp(-30) instead of log(0) = -inf.
    rate = Y.sum(axis=1) / np.exp(offsets).sum()
    intercept = np.log(np.maximum(rate, _EPS))

    U, V = _init_factors(Y, intercept, offsets, L, seed)

    def mean_of(intercept, U, V):
        eta = intercept[:, None] + offsets[None, :] + U @ V.T
        return np.exp(np.clip(eta, -_ETA_CLIP, _ETA_CLIP))

    mu = mean_of(intercept, U, V)
    deviance = [_poisson_deviance(Y, mu)]
    lr = learning_rate
    converged = False

    for _ in range(max_iter):
        keep = (intercept.copy(), U.copy(), V.copy(), mu)

        # Gene intercepts: score Σ_c(y − μ), information Σ_c μ.
        intercept = intercept + lr * (Y - mu).sum(axis=1) / np.maximum(
            mu.sum(axis=1), _EPS)
        mu = mean_of(intercept, U, V)

        # Loadings. The ridge enters both score and information, so a gene with
        # no signal left (μ ≈ 0, hence no information) is pulled toward zero
        # rather than left wherever it was initialised.
        R = Y - mu
        U = U + lr * (R @ V - penalty * U) / (mu @ V**2 + penalty)
        mu = mean_of(intercept, U, V)

        # Factors.
        R = Y - mu
        V = V + lr * (R.T @ U - penalty * V) / (mu.T @ U**2 + penalty)
        mu = mean_of(intercept, U, V)

        current = _poisson_deviance(Y, mu)
        if not np.isfinite(current) or current > deviance[-1]:
            intercept, U, V, mu = keep      # the step overshot — take it back
            lr *= 0.5
            if lr < 1e-6:
                break                       # no step size left; this is the fit
            continue

        improvement = abs(deviance[-1] - current) / (abs(deviance[-1]) + 0.1)
        deviance.append(current)
        if improvement < tol:
            converged = True
            break

    return intercept, U, V, np.array(deviance), converged


def _nb_deviance(Y: np.ndarray, mu: np.ndarray, theta: float) -> float:
    """2·Σ[y·log(y/μ) − (y+θ)·log((y+θ)/(μ+θ))] — the NB goodness-of-fit.

    The second term is the negative-binomial correction to the Poisson deviance;
    as ``θ → ∞`` it tends to ``y − μ`` and this reduces to `_poisson_deviance`. For
    ``y = 0`` the first term is 0 (as in Poisson) and the second stays finite, so a
    gene detected in no cell does not blow it up.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(Y > 0, Y * np.log(Y / mu), 0.0)
    correction = (Y + theta) * np.log((Y + theta) / (mu + theta))
    return float(2.0 * np.sum(term - correction))


def _moment_theta(Y: np.ndarray, mu: np.ndarray) -> float:
    """Method-of-moments dispersion: match ``Var = μ + μ²/θ`` in aggregate.

    ``(y − μ)² − μ`` estimates ``μ²/θ`` cell by cell, so ``θ ≈ Σμ² / Σ[(y−μ)² − μ]``.
    A non-positive denominator means the counts are *under*-dispersed relative to
    Poisson — there is no finite NB that fits tighter than Poisson, so hand back
    the ceiling (effectively Poisson). Cheap, always defined, and close enough to
    the optimum to drop Newton straight into the concave basin.
    """
    num = float(np.sum(mu ** 2))
    den = float(np.sum((Y - mu) ** 2 - mu))
    if den <= 0:
        return _THETA_MAX
    return float(np.clip(num / den, _THETA_MIN, _THETA_MAX))


def _estimate_theta(Y: np.ndarray, mu: np.ndarray, theta: float,
                    n_steps: int = 10, tol: float = 1e-3) -> float:
    """ML estimate of the shared NB dispersion given the current means.

    Newton–Raphson on the profile log-likelihood for ``θ`` with the factors held
    fixed — the same alternating scheme MASS's ``theta.ml`` uses, just shared
    across the whole matrix rather than one θ per group. The score and observed
    information are the digamma/trigamma sums below; a step that leaves the
    concave region (information ≤ 0) or the sane range is refused, so this only
    ever hands back a positive, finite θ.

    Newton on this likelihood is only reliable *near* the optimum — far out it is
    not concave and a raw step can fly off to nonsense. So ignore the incoming
    ``theta`` as a starting point and seed from the method-of-moments estimate,
    which is already in the right basin; ``theta`` still matters as the value the
    fit ran at, just not as where this search begins.
    """
    from scipy.special import digamma, polygamma

    th = _moment_theta(Y, mu)
    for _ in range(n_steps):
        score = np.sum(
            digamma(Y + th) - digamma(th) + np.log(th) + 1.0
            - np.log(mu + th) - (Y + th) / (mu + th))
        # Observed information: −∂²ℓ/∂θ². Positive near the optimum.
        info = np.sum(
            polygamma(1, th) - polygamma(1, Y + th) - 1.0 / th
            + 1.0 / (mu + th) - (Y - mu) / (mu + th) ** 2)
        if not np.isfinite(score) or not np.isfinite(info) or info <= 0:
            break                                       # not usefully concave here
        new = th + score / info
        if not np.isfinite(new):
            break
        new = float(np.clip(new, _THETA_MIN, _THETA_MAX))
        if abs(new - th) <= tol * th:
            th = new
            break
        th = new
    return th


def _fit_nb(
    Y: np.ndarray,
    L: int,
    offsets: np.ndarray,
    max_iter: int,
    tol: float,
    penalty: float,
    learning_rate: float,
    theta: float,
    optimize_theta: bool,
    seed: int,
):
    """Fisher scoring for the negative-binomial log-bilinear model.

    Mirrors `_fit_poisson` step for step; the only change is the noise model. A
    log link and NB(μ, θ) noise give an IRLS working weight ``w = μ / (1 + μ/θ)``
    and an effective residual ``(y − μ) / (1 + μ/θ)`` — the raw residual and ``μ``
    of the Poisson updates, each divided by the same ``1 + μ/θ``. As ``θ → ∞`` the
    divisor is 1 and every line below becomes its Poisson twin.

    When ``optimize_theta`` is on, ``θ`` is re-estimated by ML after each accepted
    factor step and the running deviance is re-based to that new ``θ`` (the accept
    test only ever compares two deviances measured at the *same* ``θ``). That
    re-basing is why the NB deviance trace is not promised to be monotone unless
    ``θ`` is held fixed.
    """
    rate = Y.sum(axis=1) / np.exp(offsets).sum()
    intercept = np.log(np.maximum(rate, _EPS))

    U, V = _init_factors(Y, intercept, offsets, L, seed)
    th = float(np.clip(theta, _THETA_MIN, _THETA_MAX))

    def mean_of(intercept, U, V):
        eta = intercept[:, None] + offsets[None, :] + U @ V.T
        return np.exp(np.clip(eta, -_ETA_CLIP, _ETA_CLIP))

    mu = mean_of(intercept, U, V)
    deviance = [_nb_deviance(Y, mu, th)]
    lr = learning_rate
    converged = False

    for _ in range(max_iter):
        keep = (intercept.copy(), U.copy(), V.copy(), mu)

        denom = 1.0 + mu / th                           # the NB down-weighting
        R = (Y - mu) / denom                            # effective residual
        W = mu / denom                                  # working weight

        intercept = intercept + lr * R.sum(axis=1) / np.maximum(
            W.sum(axis=1), _EPS)
        mu = mean_of(intercept, U, V)

        denom = 1.0 + mu / th
        R = (Y - mu) / denom
        W = mu / denom
        U = U + lr * (R @ V - penalty * U) / (W @ V**2 + penalty)
        mu = mean_of(intercept, U, V)

        denom = 1.0 + mu / th
        R = (Y - mu) / denom
        W = mu / denom
        V = V + lr * (R.T @ U - penalty * V) / (W.T @ U**2 + penalty)
        mu = mean_of(intercept, U, V)

        current = _nb_deviance(Y, mu, th)
        if not np.isfinite(current) or current > deviance[-1]:
            intercept, U, V, mu = keep      # the step overshot — take it back
            lr *= 0.5
            if lr < 1e-6:
                break
            continue

        improvement = abs(deviance[-1] - current) / (abs(deviance[-1]) + 0.1)
        deviance.append(current)

        if optimize_theta:
            th = _estimate_theta(Y, mu, th)
            # Re-base the baseline to the new θ so the next accept test is fair.
            deviance[-1] = _nb_deviance(Y, mu, th)

        if improvement < tol:
            converged = True
            break

    return intercept, U, V, np.array(deviance), converged, th


def _orthogonalize(U: np.ndarray, V: np.ndarray):
    """Rotate the fitted factors into PCA-like order, changing nothing that matters.

    The fit only ever sees ``U·Vᵀ``, and that product is untouched by
    ``U → U·A``, ``V → V·A⁻ᵀ``. So the raw U and V that fall out of Fisher
    scoring are one arbitrary member of a whole family, with no reason for
    component 1 to be the biggest or for the components to be uncorrelated. Pick
    the rotation that makes the loadings orthonormal and sorts the factors by the
    share of the linear predictor they carry — so "component 1" means what a
    reader coming from PCA expects it to mean.
    """
    Qu, Ru = np.linalg.qr(U)
    Qv, Rv = np.linalg.qr(V)
    W, S, Zt = np.linalg.svd(Ru @ Rv.T)     # S descending

    loadings = Qu @ W                       # orthonormal columns
    factors = (Qv @ Zt.T) * S               # ordered; these carry the scale

    # Sign is arbitrary too. Pin it, flipping both sides so U·Vᵀ is preserved.
    peak = np.argmax(np.abs(loadings), axis=0)
    signs = np.sign(loadings[peak, np.arange(loadings.shape[1])])
    signs[signs == 0] = 1.0
    return loadings * signs, factors * signs
