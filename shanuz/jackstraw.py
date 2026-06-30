from __future__ import annotations

from typing import Optional

import numpy as np
import scipy.sparse as sp


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
    return np.asarray(mat, dtype=float), features


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

    Mirrors R's JackStraw(): a small fraction (``prop_freq``) of features is
    permuted across cells and re-projected onto the (fixed) cell-embedding
    basis to build a null distribution of feature loadings per PC. Each
    observed loading is then assigned an empirical p-value. Results are stored
    on ``seurat.reductions[reduction].jackstraw`` and returned.

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

    emb = np.asarray(dr.cell_embeddings[:, :ndims], dtype=float)   # cells × ndims
    s2 = (emb ** 2).sum(axis=0)                                     # ndims
    s2[s2 == 0] = 1.0

    # Centre each feature across cells (PCA operates on centred data); the
    # projected loading of feature f on PC j is (x_f · score_j) / ||score_j||².
    Xc = X - X.mean(axis=1, keepdims=True)
    obs_stat = np.abs((Xc @ emb) / s2)                             # features × ndims

    rng = np.random.default_rng(seed)
    n_perm = max(1, int(np.ceil(prop_freq * n_features)))
    null_chunks = []
    for _ in range(num_replicate):
        idx = rng.choice(n_features, size=n_perm, replace=False)
        perm = Xc[idx, :].copy()
        for r in range(n_perm):
            perm[r, :] = rng.permutation(perm[r, :])
        null_chunks.append(np.abs((perm @ emb) / s2))             # n_perm × ndims
    null_all = np.vstack(null_chunks)                              # (R·n_perm) × ndims

    # Empirical p per (feature, PC) = fraction of null loadings ≥ observed.
    empirical = np.empty((n_features, ndims))
    n_null = null_all.shape[0]
    for j in range(ndims):
        col = np.sort(null_all[:, j])
        ranks = np.searchsorted(col, obs_stat[:, j], side="left")
        empirical[:, j] = (n_null - ranks) / n_null

    js = JackStrawData(
        empirical_p_values=empirical,
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

    Mirrors R's ScoreJackStraw(): a PC whose feature p-values are skewed toward
    zero (relative to a uniform null) is significant. We quantify that skew with
    a one-sided KS test against Uniform(0, 1); a *small* returned p-value marks
    a significant PC. Results are stored on the reduction's JackStrawData.
    """
    from scipy.stats import kstest

    if reduction not in seurat.reductions:
        raise KeyError(f"Reduction '{reduction}' not found.")
    dr = seurat.reductions[reduction]
    js = dr.jackstraw
    if js is None or js.is_empty():
        raise ValueError("Run jack_straw() before score_jackstraw().")

    emp = js.empirical_p_values
    ndims = emp.shape[1] if dims is None else int(min(dims, emp.shape[1]))

    overall = np.ones(ndims)
    for j in range(ndims):
        pj = np.clip(emp[:, j], 0.0, 1.0)
        # 'greater': sample CDF lies above uniform's → mass concentrated at low
        # p-values → significant PC. Small KS p-value ⇒ significant.
        overall[j] = kstest(pj, "uniform", alternative="greater").pvalue

    js.overall_p_values = overall
    return overall
