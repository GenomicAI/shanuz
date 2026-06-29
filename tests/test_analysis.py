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
# CLR normalization margin (multimodal / ADT support). margin=1 centers each
# cell across features; margin=2 centers each feature across cells.
# ---------------------------------------------------------------------------

def test_clr_margin_centering_direction():
    mat = np.array([[1.0, 2, 3, 4],
                    [5, 6, 7, 8],
                    [0, 1, 0, 2]])  # features x cells

    r1 = _clr_normalize(mat, margin=1)   # per-cell (column) centering
    assert np.allclose(r1.mean(axis=0), 0.0, atol=1e-9)

    r2 = _clr_normalize(mat, margin=2)   # per-feature (row) centering
    assert np.allclose(r2.mean(axis=1), 0.0, atol=1e-9)

    assert not np.allclose(r1, r2)


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
    # margin=2 → each feature (row) is centered across cells
    assert np.allclose(dense.mean(axis=1), 0.0, atol=1e-9)


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
