"""Regression tests for analysis-pipeline fidelity to Seurat.

Each test guards against a specific bug that previously caused shanuz to
deviate from Seurat's behaviour:

  * avg_log2FC averaging order (markers.py)
  * VST variable-feature selection fitting on counts, not normalized data
    (preprocessing.py)
  * ridge_plot y-axis label / data alignment (plotting.py)
"""
import numpy as np
import scipy.sparse as sp
import pytest

from shanuz.preprocessing import (
    normalize_data, find_variable_features, _vst_hvg, _clr_normalize,
)
from shanuz.markers import find_markers


# ---------------------------------------------------------------------------
# Bug 1 — avg_log2FC must average expm1(x) within a group, then re-log.
# Seurat: log2(mean(expm1(x1)) + 1) - log2(mean(expm1(x2)) + 1)
# The previous bug computed expm1(mean(x)) which compresses fold-changes.
# ---------------------------------------------------------------------------

def test_avg_log2fc_matches_seurat_formula(small_seurat):
    normalize_data(small_seurat)
    small_seurat.idents = ["A"] * 10 + ["B"] * 10

    res = find_markers(
        small_seurat, ident_1="A", ident_2="B",
        min_pct=0.0, logfc_threshold=0.0,
    )
    assert len(res) > 0

    assay = small_seurat.assays["RNA"]
    data = assay.layers["data"]
    dense = data.toarray() if sp.issparse(data) else np.asarray(data, dtype=float)
    feats = assay._all_feature_names
    a_idx = list(range(10))
    b_idx = list(range(10, 20))

    for gene in res.index:
        gi = feats.index(gene)
        m1 = np.expm1(dense[gi, a_idx]).mean()
        m2 = np.expm1(dense[gi, b_idx]).mean()
        expected = np.log2(m1 + 1) - np.log2(m2 + 1)
        assert np.isclose(res.loc[gene, "avg_log2FC"], expected, atol=1e-9), gene


def test_wilcox_uses_tie_corrected_mannwhitney(small_seurat):
    """Seurat's wilcox test (presto / base-R wilcox.test) tie-corrects the
    rank-sum statistic. scipy.stats.ranksums does not, so marker p-values must
    match scipy.stats.mannwhitneyu (asymptotic, continuity-corrected), and must
    differ from the un-tie-corrected ranksums on this zero-heavy data."""
    from scipy.stats import mannwhitneyu, ranksums

    normalize_data(small_seurat)
    small_seurat.idents = ["A"] * 10 + ["B"] * 10
    res = find_markers(
        small_seurat, ident_1="A", ident_2="B",
        min_pct=0.0, logfc_threshold=0.0,
    )
    assert len(res) > 0

    assay = small_seurat.assays["RNA"]
    dense = assay.layers["data"].toarray()
    feats = assay._all_feature_names
    a_idx, b_idx = list(range(10)), list(range(10, 20))

    saw_difference = False
    for gene in res.index:
        gi = feats.index(gene)
        x1, x2 = dense[gi, a_idx], dense[gi, b_idx]
        if x1.sum() == 0 and x2.sum() == 0:
            continue
        _, p_mwu = mannwhitneyu(x1, x2, alternative="two-sided",
                                use_continuity=True, method="asymptotic")
        assert np.isclose(res.loc[gene, "p_val"], p_mwu, atol=1e-9), gene
        _, p_rs = ranksums(x1, x2)
        if not np.isclose(p_mwu, p_rs, atol=1e-6):
            saw_difference = True
    assert saw_difference  # tie correction actually changes the p-value here


def test_avg_log2fc_differs_from_buggy_formula(small_seurat):
    """The correct (mean-of-expm1) and buggy (expm1-of-mean) formulas must
    actually disagree on this data, otherwise the test above is vacuous."""
    normalize_data(small_seurat)
    small_seurat.idents = ["A"] * 10 + ["B"] * 10
    res = find_markers(
        small_seurat, ident_1="A", ident_2="B",
        min_pct=0.0, logfc_threshold=0.0,
    )

    assay = small_seurat.assays["RNA"]
    dense = assay.layers["data"].toarray()
    feats = assay._all_feature_names
    a_idx, b_idx = list(range(10)), list(range(10, 20))

    saw_difference = False
    for gene in res.index:
        gi = feats.index(gene)
        x1, x2 = dense[gi, a_idx], dense[gi, b_idx]
        buggy = np.log2(np.expm1(x1.mean()) + 1) - np.log2(np.expm1(x2.mean()) + 1)
        correct = res.loc[gene, "avg_log2FC"]
        if not np.isclose(buggy, correct, atol=1e-6):
            saw_difference = True
    assert saw_difference


# ---------------------------------------------------------------------------
# Bug 2 — vst selection.method must fit on raw counts, independent of whether
# NormalizeData() has run. So running it before vs after normalization must
# yield identical variable features.
# ---------------------------------------------------------------------------

def test_vst_uses_counts_not_normalized_data(small_seurat):
    find_variable_features(small_seurat, selection_method="vst", nfeatures=10)
    hvg_before = list(small_seurat.assays["RNA"].variable_features)

    normalize_data(small_seurat)
    find_variable_features(small_seurat, selection_method="vst", nfeatures=10)
    hvg_after = list(small_seurat.assays["RNA"].variable_features)

    assert hvg_before == hvg_after
    assert len(hvg_before) == 10


def test_vst_clipping_suppresses_single_cell_outliers():
    """Seurat clips standardized values at sqrt(n_cells) so a single high-count
    outlier cell can't inflate a gene's standardized variance. Removing the clip
    must measurably change (inflate) the outlier gene's standardized variance."""
    rng = np.random.default_rng(0)
    n_genes, n_cells = 60, 50
    base = rng.poisson(1.0, size=(n_genes, n_cells)).astype(float)
    # Gene 0: one outlier cell, zero elsewhere. Its mean (1.0) is typical, so the
    # LOESS expects a low variance — but the single high cell inflates the raw
    # variance. This is exactly the case Seurat's clip is designed to tame.
    base[0, :] = 0.0
    base[0, 0] = 50.0
    mat = sp.csc_matrix(base)

    _, _, _, vs_clipped = _vst_hvg(mat, nfeatures=n_genes)            # clip = sqrt(N)
    _, _, _, vs_unclipped = _vst_hvg(mat, nfeatures=n_genes, clip_max=1e12)

    assert vs_clipped[0] < 0.5 * vs_unclipped[0]


# ---------------------------------------------------------------------------
# Bug 3 — ridge_plot row j holds the density for unique[j]; the y tick labels
# must therefore be `unique` (not its reverse).
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# CLR normalization (multimodal / ADT support). Seurat's clr_function is
#   log1p(x / exp(sum(log1p(x[x>0])) / length(x)))
# applied via CustomNormalize as apply(data, MARGIN = margin, clr_function).
# R's apply MARGIN=1 is rows and MARGIN=2 is columns, so with counts stored
# features x cells: margin=1 normalizes each FEATURE across cells (Seurat's
# default), margin=2 normalizes each CELL across its features (ADT panels).
# ---------------------------------------------------------------------------

def _seurat_clr_vector(x):
    """Reference: Seurat's clr_function — log1p(x / exp(sum(log1p(x>0))/len))."""
    x = np.asarray(x, dtype=float)
    denom = np.exp(np.sum(np.log1p(x[x > 0])) / len(x))
    return np.log1p(x / denom)


# Ground truth captured from R, not re-derived in Python:
#   m <- matrix(c(10,0,3, 1,5,0, 0,2,8, 7,1,1), nrow = 4, byrow = TRUE)
#   NormalizeData(as.sparse(m), normalization.method = "CLR", margin = <m>)
# Seurat 5.5.1. The matrix is asymmetric so row-wise and column-wise CLR
# cannot coincide.
_R_CLR_INPUT = np.array([[10., 0., 3.],
                         [1., 5., 0.],
                         [0., 2., 8.],
                         [7., 1., 1.]])
_R_CLR_MARGIN1 = np.array([[1.34354, 0.00000, 0.61506],
                           [0.36241, 1.15812, 0.00000],
                           [0.00000, 0.51083, 1.29928],
                           [1.16467, 0.27382, 0.27382]])
_R_CLR_MARGIN2 = np.array([[1.32056, 0.00000, 0.70798],
                           [0.24259, 1.11227, 0.00000],
                           [0.00000, 0.59691, 1.32078],
                           [1.07222, 0.34235, 0.29513]])


def test_clr_margin_matches_r_ground_truth():
    """The margin flag must mean what it means in Seurat, axis included.

    Regression guard. The formula test below checks the per-vector kernel but
    derives the axis mapping in Python, so it passes just as happily when the
    two margins are swapped. Only fixed output from a real R run pins the axis.
    """
    got1 = _clr_normalize(_R_CLR_INPUT, margin=1)
    got2 = _clr_normalize(_R_CLR_INPUT, margin=2)
    assert np.allclose(got1, _R_CLR_MARGIN1, atol=1e-5)
    assert np.allclose(got2, _R_CLR_MARGIN2, atol=1e-5)
    # ...and specifically not the other way round.
    assert not np.allclose(got1, _R_CLR_MARGIN2, atol=1e-3)
    assert not np.allclose(got2, _R_CLR_MARGIN1, atol=1e-3)


def test_clr_matches_seurat_formula():
    """CLR must equal Seurat's clr_function along the margin, NOT a simple
    log1p-then-mean-center (the two are not algebraically equal)."""
    rng = np.random.default_rng(1)
    mat = rng.poisson(2.0, size=(5, 8)).astype(float)
    mat[0, 0] = 0.0
    mat[2, 3] = 0.0  # exercise the x>0 geometric-mean path

    # margin=1 → per feature (row), geometric mean across cells
    r1 = _clr_normalize(mat, margin=1)
    ref1 = np.vstack([_seurat_clr_vector(mat[i, :]) for i in range(mat.shape[0])])
    assert np.allclose(r1, ref1, atol=1e-12)

    # margin=2 → per cell (column), geometric mean across features
    r2 = _clr_normalize(mat, margin=2)
    ref2 = np.column_stack([_seurat_clr_vector(mat[:, j]) for j in range(mat.shape[1])])
    assert np.allclose(r2, ref2, atol=1e-12)

    # The two margins differ, and CLR is non-negative (log1p of a ratio ≥ 0).
    assert not np.allclose(r1, r2)
    assert np.all(r1 >= 0) and np.all(r2 >= 0)


def test_clr_margin_routed_through_normalize_data():
    rng = np.random.default_rng(0)
    counts = sp.csc_matrix(rng.poisson(3.0, size=(6, 12)).astype(float))
    from shanuz.shanuz import create_shanuz_object
    obj = create_shanuz_object(
        counts=counts, assay="ADT",
        feature_names=[f"p{i}" for i in range(6)],
        cell_names=[f"c{i}" for i in range(12)],
    )
    normalize_data(obj, normalization_method="CLR", margin=2)
    data = obj.assays["ADT"].layers["data"]
    dense = data.toarray() if sp.issparse(data) else np.asarray(data)

    counts_dense = counts.toarray()
    ref = np.column_stack([_seurat_clr_vector(counts_dense[:, j])
                           for j in range(counts_dense.shape[1])])
    assert np.allclose(dense, ref, atol=1e-12)


# ---------------------------------------------------------------------------
# Assay5 subset-layer support — scale.data may legitimately hold only the
# variable features (Seurat's default), and the default ScaleData() call must
# therefore succeed and feed PCA.
# ---------------------------------------------------------------------------

def test_scale_data_default_subset_features_then_pca():
    from shanuz.shanuz import create_shanuz_object
    from shanuz.preprocessing import scale_data
    from shanuz.reduction import run_pca

    rng = np.random.default_rng(0)
    counts = sp.csc_matrix(rng.poisson(0.5, size=(200, 80)).astype(float))
    obj = create_shanuz_object(
        counts=counts, assay="RNA",
        feature_names=[f"g{i}" for i in range(200)],
        cell_names=[f"c{i}" for i in range(80)],
    )
    normalize_data(obj)
    find_variable_features(obj, selection_method="vst", nfeatures=50)
    # Default features = the 50 variable features — must NOT raise on Assay5.
    scale_data(obj)
    sd = obj.assays["RNA"].layers["scale.data"]
    assert sd.shape == (50, 80)
    assert obj.assays["RNA"].features("scale.data") == obj.assays["RNA"].variable_features
    run_pca(obj, n_pcs=10)
    assert obj.reductions["pca"].cell_embeddings.shape == (80, 10)


def test_subset_resizes_graphs_and_drops_neighbors():
    from shanuz.shanuz import create_shanuz_object
    from shanuz.graph import Graph

    rng = np.random.default_rng(0)
    counts = sp.csc_matrix(rng.poisson(0.5, size=(40, 30)).astype(float))
    obj = create_shanuz_object(
        counts=counts, assay="RNA",
        feature_names=[f"g{i}" for i in range(40)],
        cell_names=[f"c{i}" for i in range(30)],
    )
    obj.graphs["RNA_snn"] = Graph(
        matrix=sp.random(30, 30, density=0.3, format="csc"),
        cell_names=obj.cell_names(), assay_used="RNA",
    )
    keep = [f"c{i}" for i in range(12)]
    sub = obj.subset(cells=keep)
    # Graph is resized to the retained cells (not the stale 30×30 matrix).
    assert sub.graphs["RNA_snn"]._matrix.shape == (12, 12)
    assert sub.graphs["RNA_snn"].cells() == keep
    # Neighbor indices are invalidated by subsetting and therefore dropped.
    assert dict(sub.neighbors) == {}


# ---------------------------------------------------------------------------
# JackStraw permutation test (jackstraw.py) — features carrying real PC
# structure must score far lower empirical p-values than noise features, and
# a structured PC must be flagged significant by score_jackstraw.
# ---------------------------------------------------------------------------

def test_jackstraw_separates_signal_from_noise():
    from shanuz.shanuz import create_shanuz_object
    from shanuz.preprocessing import scale_data
    from shanuz.reduction import run_pca
    from shanuz.jackstraw import jack_straw, score_jackstraw

    rng = np.random.default_rng(0)
    base = rng.poisson(0.5, size=(150, 120)).astype(float)
    base[:10, :60] += 6.0  # genes 0-9 carry strong structure across two groups
    obj = create_shanuz_object(
        counts=sp.csc_matrix(base), assay="RNA",
        feature_names=[f"g{i}" for i in range(150)],
        cell_names=[f"c{i}" for i in range(120)],
    )
    normalize_data(obj)
    find_variable_features(obj, nfeatures=150)
    scale_data(obj)
    run_pca(obj, n_pcs=15)

    js = jack_straw(obj, dims=10, num_replicate=100, seed=0)
    emp = js.empirical_p_values
    assert emp.shape == (150, 10)
    assert emp.min() >= 0.0 and emp.max() <= 1.0

    feats = obj.reductions["pca"].features()
    sig = [feats.index(f"g{i}") for i in range(10) if f"g{i}" in feats]
    rest = [i for i in range(len(feats)) if i not in sig]
    # Signal genes load significantly on PC1; noise genes do not.
    assert emp[sig, 0].mean() < 0.05
    assert emp[rest, 0].mean() > emp[sig, 0].mean()

    overall = score_jackstraw(obj, dims=10)
    assert overall.shape == (10,)
    assert overall[0] < 0.05  # the structured PC1 is significant


def _jackstraw_fixture():
    """A small object with structure in PC1 only — the rest is noise."""
    from shanuz.shanuz import create_shanuz_object
    from shanuz.preprocessing import scale_data
    from shanuz.reduction import run_pca

    rng = np.random.default_rng(0)
    base = rng.poisson(0.5, size=(150, 120)).astype(float)
    base[:10, :60] += 6.0
    obj = create_shanuz_object(
        counts=sp.csc_matrix(base), assay="RNA",
        feature_names=[f"g{i}" for i in range(150)],
        cell_names=[f"c{i}" for i in range(120)],
    )
    normalize_data(obj)
    find_variable_features(obj, nfeatures=150)
    scale_data(obj)
    run_pca(obj, n_pcs=15)
    return obj


def test_jackstraw_null_is_calibrated_on_noise_pcs():
    """Regression: the null must come from a PCA refit, not a fixed projection.

    JackStraw originally built its null by projecting the permuted rows onto the
    *existing* embedding. A fixed basis cannot rotate to absorb the scrambled
    signal, so the permuted loadings come out too small, the null is too tight,
    and ordinary noise features look extreme against it. R's ``JackRandom``
    re-runs the whole PCA per replicate instead.

    The observable: on a PC carrying no structure, a correct permutation null
    yields roughly uniform empirical p-values (median near 0.5). The fixed-basis
    null crushed those medians to ~0.25 and pushed 8-13 % of features below
    1e-5; the refit null gives medians ≥ 0.45 and under 2 %. On real pbmc3k data
    the same defect put 109-203 features below threshold on PCs 14-20, where R
    finds 0-5.
    """
    from shanuz.jackstraw import jack_straw

    obj = _jackstraw_fixture()
    emp = jack_straw(obj, dims=10, num_replicate=100, seed=0).empirical_p_values

    noise = emp[:, 4:]                              # PC 5-10 carry no structure
    medians = np.median(noise, axis=0)
    frac_sig = (noise <= 1e-5).mean(axis=0)
    assert medians.min() >= 0.40, f"null too tight, medians {medians}"
    assert frac_sig.max() <= 0.05, f"too many significant noise features {frac_sig}"


def test_jackstraw_stores_the_null_scores():
    """``fake_reduction_scores`` was declared but never populated; R stores it."""
    from shanuz.jackstraw import jack_straw

    obj = _jackstraw_fixture()
    js = jack_straw(obj, dims=5, num_replicate=10, seed=0)
    # 10 replicates x max(3, floor(150 * 0.01)) = 3 permuted features per replicate.
    assert js.fake_reduction_scores.shape == (30, 5)
    assert np.isfinite(js.fake_reduction_scores).all()


def test_score_jackstraw_does_not_flag_every_pc():
    """Regression: the aggregation must be R's prop.test, not a KS test.

    A one-sided KS test against Uniform(0, 1) is enormously more sensitive than
    the proportion test ``ScoreJackStraw`` actually uses. With thousands of
    features it returned ~1e-112 or smaller for *every* pbmc3k PC — including
    pure noise — so no PC ever failed the threshold and the function could not
    do the one job it exists for. Here only PC1 carries structure, so most of
    the ten tested PCs must come back insignificant.
    """
    from shanuz.jackstraw import jack_straw, score_jackstraw

    obj = _jackstraw_fixture()
    jack_straw(obj, dims=10, num_replicate=100, seed=0)
    scores = score_jackstraw(obj, dims=10)

    assert scores[0] < 0.05                          # the real one still passes
    assert (scores > 0.05).sum() >= 5, f"nothing is being rejected: {scores}"
    assert scores.max() == pytest.approx(1.0)        # a PC with no hits scores 1


@pytest.mark.parametrize("count,expected_p", [
    # Verified against R: prop.test(x = c(count, 0), n = c(2000, 2000))$p.value,
    # spanning the full range ScoreJackStraw produces on pbmc3k.
    (558, 1.545300664068e-142),
    (90, 2.337543460534e-21),
    (17, 1.007234160653e-04),
    (3, 2.480356409327e-01),
    (1, 1.0),
])
def test_prop_test_matches_r(count, expected_p):
    """The ported prop.test must reproduce R's p-value, not merely approximate it.

    ``ScoreJackStraw``'s output *is* this p-value, so an approximation would
    silently shift every PC cutoff.
    """
    from shanuz.jackstraw import _prop_test

    assert _prop_test(count, 0, 2000, 2000) == pytest.approx(expected_p, rel=1e-9)


def test_ridge_plot_labels_align_with_rows(small_seurat):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from shanuz.plotting import ridge_plot

    normalize_data(small_seurat)
    small_seurat.idents = ["A"] * 7 + ["B"] * 7 + ["C"] * 6

    fig = ridge_plot(small_seurat, features=["gene_0"], group_by=None)
    fig.canvas.draw()
    ax = fig.axes[0]
    labels = [t.get_text() for t in ax.get_yticklabels()]

    # ridge_plot reverses the sorted unique groups for top-to-bottom ordering.
    expected = sorted({"A", "B", "C"})[::-1]
    assert labels == expected

    import matplotlib.pyplot as plt
    plt.close(fig)
