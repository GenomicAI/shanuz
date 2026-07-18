from __future__ import annotations

from typing import Optional

import numpy as np


class JackStrawData:
    """Stores JackStraw permutation test results for a DimReduc.

    Mirrors R's JackStraw / JackStrawData from jackstraw.R.
    """

    __slots__ = (
        "empirical_p_values",
        "fake_reduction_scores",
        "overall_p_values",
        "score",
        "method",
    )

    def __init__(
        self,
        empirical_p_values: Optional[np.ndarray] = None,
        fake_reduction_scores: Optional[np.ndarray] = None,
        overall_p_values: Optional[np.ndarray] = None,
        score: Optional[np.ndarray] = None,
        method: Optional[str] = None,
    ) -> None:
        self.empirical_p_values = empirical_p_values
        self.fake_reduction_scores = fake_reduction_scores
        self.overall_p_values = overall_p_values
        self.score = score
        self.method = method

    def is_empty(self) -> bool:
        return self.empirical_p_values is None

    def __repr__(self) -> str:
        if self.is_empty():
            return "JackStrawData(empty)"
        shape = self.empirical_p_values.shape if self.empirical_p_values is not None else "?"
        return f"JackStrawData(p_values shape={shape}, method={self.method!r})"


# ---------------------------------------------------------------------------
# JackStraw permutation test  (mirrors R JackStraw() / ScoreJackStraw())
# ---------------------------------------------------------------------------

def _scaled_matrix_for_reduction(seurat, dr, layer: str):
    """Return (features_used × cells) scaled matrix in the reduction's feature
    order, plus the feature-name list. Falls back to the assay's scale.data."""
    from .reduction import _get_scaled_data

    assay_obj = seurat.assays[dr.assay_used or seurat.active_assay]
    features = list(dr.features())
    if not features:
        # No stored loading features — use the assay's variable features.
        from .assay5 import Assay5
        features = (
            assay_obj.variable_features if isinstance(assay_obj, Assay5)
            else assay_obj.var_features
        ) or assay_obj.features()
    mat = _get_scaled_data(assay_obj, features, layer)
    # Own the buffer: jack_straw scrambles rows in place while building the null,
    # and must never disturb the layer it read them from.
    return np.array(mat, dtype=float, copy=True), features


def _refit_loadings(mat: np.ndarray, ndims: int, seed: int) -> np.ndarray:
    """Feature loadings from a fresh PCA of ``mat`` (features × cells).

    The null in R's ``JackRandom`` comes from re-running the whole PCA on the
    permuted matrix, not from projecting onto the existing basis — a refit basis
    can rotate to absorb some of the scrambled signal, which is what gives the
    null its true spread. Uses the same estimator as :func:`run_pca` so the null
    loadings are on the same scale as the observed ones.
    """
    from sklearn.decomposition import PCA

    pca = PCA(n_components=ndims, random_state=seed)
    pca.fit(mat.T)                      # PCA expects (cells × features)
    return pca.components_.T            # (features × ndims)


def _prop_test(x1: float, x2: float, n1: float, n2: float) -> float:
    """R's ``prop.test(x = c(x1, x2), n = c(n1, n2))`` p-value.

    A two-sample test for equality of proportions: Yates-corrected chi-square on
    one degree of freedom. Ported rather than approximated because
    ``ScoreJackStraw``'s output *is* this p-value — verified to reproduce R's
    values on pbmc3k to nine significant figures across the full range
    (1e-143 to 1.0).
    """
    from scipy.stats import chi2

    x = np.array([x1, x2], dtype=float)
    n = np.array([n1, n2], dtype=float)
    estimate = x / n
    delta = estimate[0] - estimate[1]
    yates = min(0.5, abs(delta) / np.sum(1.0 / n))
    p = x.sum() / n.sum()
    if p <= 0.0 or p >= 1.0:
        # Degenerate: every observation on one side. R returns NaN here; the
        # caller's "no features below threshold" guard handles the real case.
        return 1.0
    observed = np.column_stack([x, n - x])
    expected = np.column_stack([n * p, n * (1.0 - p)])
    stat = np.sum((np.abs(observed - expected) - yates) ** 2 / expected)
    return float(chi2.sf(stat, 1))


def jack_straw(
    seurat,
    reduction: str = "pca",
    dims: int = 20,
    num_replicate: int = 100,
    prop_freq: float = 0.01,
    layer: str = "scale.data",
    seed: int = 42,
) -> "JackStrawData":
    """Permutation test for the significance of PCA dimensions.

    Mirrors R's ``JackStraw()``: a small fraction (``prop_freq``) of features is
    permuted across cells, the **PCA is re-run on the permuted matrix**, and the
    permuted features' loadings in that refit basis form the null distribution
    per PC. Each observed loading is then assigned an empirical p-value. Results
    are stored on ``seurat.reductions[reduction].jackstraw`` and returned.

    The refit is the expensive part and it is not optional: an earlier version of
    this function built the null by projecting the permuted rows onto the
    *fixed* original embedding, which is far cheaper but produces a much tighter
    null — a fixed basis cannot rotate to absorb the scrambled signal, so the
    permuted loadings come out too small and ordinary noise features look
    extreme against them. On pbmc3k that inflated the count of "significant"
    features on the pure-noise PCs 14-20 from R's 0-5 to 109-203, and left
    :func:`score_jackstraw` unable to reject any PC at all.

    Cost scales as ``num_replicate`` full PCAs; ~1-2 minutes for the Seurat
    defaults on a 2000-feature, 2700-cell object. Lower ``num_replicate`` when
    iterating, but note it also sets the p-value resolution: the smallest
    non-zero empirical p is ``1 / (num_replicate * n_permuted)``.

    Parameters
    ----------
    reduction      : reduction to test (default 'pca')
    dims           : number of PCs to score
    num_replicate  : permutation replicates (Seurat default 100)
    prop_freq      : fraction of features permuted per replicate (default 0.01)
    layer          : scaled layer feeding the reduction
    seed           : RNG seed
    """
    if reduction not in seurat.reductions:
        raise KeyError(f"Reduction '{reduction}' not found. Run run_pca() first.")
    dr = seurat.reductions[reduction]

    X, features = _scaled_matrix_for_reduction(seurat, dr, layer)
    n_features, n_cells = X.shape
    ndims = int(min(dims, dr.cell_embeddings.shape[1]))

    # The observed statistic is the reduction's own feature loadings, exactly as
    # R takes Loadings(object[[reduction]], projected = FALSE).
    loadings = np.asarray(dr.feature_loadings, dtype=float)
    if loadings.size == 0:
        raise ValueError(
            f"Reduction '{reduction}' has no feature loadings; JackStraw needs "
            "them as the observed statistic. Re-run run_pca().")
    if loadings.shape[0] != n_features:
        raise ValueError(
            f"Reduction '{reduction}' has {loadings.shape[0]} loadings for "
            f"{n_features} scaled features — the reduction and the layer "
            "disagree about the feature set.")
    obs_stat = np.abs(loadings[:, :ndims])                         # features × ndims

    # R: sample(rownames, size = nrow * prop.use), floored, with a hard floor of 3.
    n_perm = max(3, int(n_features * prop_freq))
    if n_perm > n_features:
        raise ValueError(
            f"prop_freq={prop_freq} selects {n_perm} of only {n_features} features")

    rng = np.random.default_rng(seed)
    null_chunks = []
    for _ in range(num_replicate):
        idx = rng.choice(n_features, size=n_perm, replace=False)
        # Scramble in place and restore afterwards: copying the whole matrix per
        # replicate would dominate the runtime on a real object.
        saved = X[idx, :].copy()
        for r in idx:
            X[r, :] = rng.permutation(X[r, :])
        refit = _refit_loadings(X, ndims, seed)
        null_chunks.append(np.abs(refit[idx, :]))                  # n_perm × ndims
        X[idx, :] = saved
    null_all = np.vstack(null_chunks)                              # (R·n_perm) × ndims

    # R's EmpiricalP: the fraction of null loadings STRICTLY greater than the
    # observed one. searchsorted(side='right') counts null <= obs.
    empirical = np.empty((n_features, ndims))
    n_null = null_all.shape[0]
    for j in range(ndims):
        col = np.sort(null_all[:, j])
        ranks = np.searchsorted(col, obs_stat[:, j], side="right")
        empirical[:, j] = (n_null - ranks) / n_null

    js = JackStrawData(
        empirical_p_values=empirical,
        fake_reduction_scores=null_all,
        overall_p_values=None,
        score=obs_stat,
        method="jackstraw",
    )
    dr.jackstraw = js
    return js


def score_jackstraw(
    seurat,
    reduction: str = "pca",
    dims: Optional[int] = None,
    score_thresh: float = 1e-5,
) -> np.ndarray:
    """Aggregate per-feature JackStraw p-values into one p-value per PC.

    Mirrors R's ``ScoreJackStraw()``: count the features whose empirical p-value
    falls at or below ``score_thresh``, and test that count against the number
    expected under a uniform null (``floor(n_features * score_thresh)``) with a
    two-proportion test. A *small* returned value marks a significant PC. A PC
    with no feature below the threshold scores exactly 1, as in R.

    Not a distributional goodness-of-fit test. An earlier version used a
    one-sided KS test against Uniform(0, 1), which is enormously more sensitive:
    with thousands of features it returned p-values around 1e-112 or smaller for
    *every* PC on pbmc3k, including pure noise, so no PC ever failed and the
    function could not do the one job it exists for — telling you where to cut.
    """
    if reduction not in seurat.reductions:
        raise KeyError(f"Reduction '{reduction}' not found.")
    dr = seurat.reductions[reduction]
    js = dr.jackstraw
    if js is None or js.is_empty():
        raise ValueError("Run jack_straw() before score_jackstraw().")

    emp = js.empirical_p_values
    ndims = emp.shape[1] if dims is None else int(min(dims, emp.shape[1]))
    n_features = emp.shape[0]
    expected = float(np.floor(n_features * score_thresh))

    overall = np.ones(ndims)
    for j in range(ndims):
        observed = int((np.clip(emp[:, j], 0.0, 1.0) <= score_thresh).sum())
        if observed == 0:
            overall[j] = 1.0                 # R's explicit guard
        else:
            overall[j] = _prop_test(observed, expected, n_features, n_features)

    js.overall_p_values = overall
    return overall
