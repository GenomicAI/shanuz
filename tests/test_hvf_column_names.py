"""`find_variable_features` must write the columns `HVFInfo()` returns.

SeuratObject keeps per-feature HVF statistics in the assay's meta.data behind a
method-and-layer prefix (`vf_vst_counts_mean`) and strips it on the way out, so
the names a Seurat user actually sees are `HVFInfo()`'s. truecell has no prefix
layer — `assay.meta_data` *is* the user-facing table — which makes matching
`HVFInfo()` the only way a column can read the same in either language.

It did not. truecell wrote `means` / `variances` / `variances.standardized`
against Seurat's `mean` / `variance` / `variance.standardized`, dropped
`variance.expected` entirely, and wrote scaled *dispersions* from the mvp path
into a column named `variances.standardized`. Nothing failed: every consumer was
in-tree and used truecell's spelling, so the divergence only surfaced where the
two tools met — the pbmc3k handoff renamed all three columns on the way out to
make the join work, and carried a comment saying the divergence was real.

Measured against Seurat 5.5.1 / SeuratObject 5.4.0 on an Assay5:

    HVFInfo(o[["RNA"]])                 mean, variance, variance.expected,
                                        variance.standardized
    HVFInfo(o[["RNA"]], status = TRUE)  ... plus variable, rank
    HVFInfo(o2[["RNA"]], "mvp")         mvp.mean, mvp.dispersion,
                                        mvp.dispersion.scaled

`variable`/`rank` are deliberately not mirrored: truecell spells the first
`highly_variable` and keeps rank as the order of `assay.variable_features`.
"""
import numpy as np
import pytest
import scipy.sparse as sp

from truecell.preprocessing import find_variable_features, normalize_data

# Verbatim from `HVFInfo()` under Seurat 5.5.1 — see the module docstring.
VST_COLUMNS = ["mean", "variance", "variance.expected", "variance.standardized"]
MVP_COLUMNS = ["mvp.mean", "mvp.dispersion", "mvp.dispersion.scaled"]

# The spellings truecell used before. Listed so a partial revert fails loudly
# rather than leaving one method writing the old names beside the new ones.
RETIRED_COLUMNS = ["means", "variances", "variances.standardized"]


def _fitted(small_seurat, method):
    normalize_data(small_seurat)
    find_variable_features(small_seurat, selection_method=method, nfeatures=5)
    return small_seurat.assays["RNA"].meta_data


def test_vst_writes_the_hvfinfo_column_names(small_seurat):
    md = _fitted(small_seurat, "vst")
    missing = [c for c in VST_COLUMNS if c not in md.columns]
    assert not missing, (
        f"vst is missing {missing}; `HVFInfo()` returns exactly {VST_COLUMNS}"
    )


def test_mvp_writes_the_hvfinfo_column_names(small_seurat):
    """A different method reports different quantities, under different names.

    `HVFInfo(method = "mvp")` never returns a `variance.standardized`. Sharing
    one set of column names across both methods put scaled dispersions in a
    column claiming to hold standardized variances, which no reader downstream
    could detect.
    """
    md = _fitted(small_seurat, "mvp")
    missing = [c for c in MVP_COLUMNS if c not in md.columns]
    assert not missing, f"mvp is missing {missing}"
    leaked = [c for c in VST_COLUMNS if c in md.columns]
    assert not leaked, (
        f"mvp wrote vst column(s) {leaked}; those names promise a LOESS "
        f"variance fit that the dispersion path never computes"
    )


@pytest.mark.parametrize("method", ["vst", "mvp", "dispersion", "mean.var.plot"])
def test_no_method_writes_the_retired_column_names(small_seurat, method):
    md = _fitted(small_seurat, method)
    survivors = [c for c in RETIRED_COLUMNS if c in md.columns]
    assert not survivors, (
        f"{method} still writes {survivors}; these are truecell's old spellings "
        f"and have no counterpart in Seurat's HVFInfo()"
    )


def test_variance_expected_column_holds_the_loess_fit_not_a_placeholder():
    """Pins the column's *contents*, not just its presence.

    Where no value is clipped, the definitions collapse into an exact identity:

        variance.standardized = Σ((x-μ)/σ̂)² / (n-1) = variance / variance.expected

    so a `variance.expected` holding the observed variance, its square root, or
    a copy of some neighbouring column fails here even though a
    column-names-only check would pass. Genes that clip are excluded — for them
    the identity legitimately does not hold, and that is the point of the clip.

    Deliberately routed through `find_variable_features` and read back off
    `meta_data` rather than unpacked from `_vst_hvg`. Both mutations that swap
    the column's contents live at the *storage* site, and a test that calls the
    computation directly passes through every one of them: it would be checking
    that `_vst_hvg` returns the right five arrays while the wrong one gets
    written to the column.
    """
    from truecell import create_truecell_object

    rng = np.random.default_rng(0)
    n_genes, n_cells = 60, 40
    counts = sp.csc_matrix(rng.poisson(4.0, size=(n_genes, n_cells)).astype(float))
    obj = create_truecell_object(
        counts=counts,
        assay="RNA",
        feature_names=[f"gene{i}" for i in range(n_genes)],
        cell_names=[f"cell{j}" for j in range(n_cells)],
    )
    find_variable_features(obj, selection_method="vst", nfeatures=10)
    md = obj.assays["RNA"].meta_data

    mean = md["mean"].to_numpy()
    variance = md["variance"].to_numpy()
    expected = md["variance.expected"].to_numpy()
    var_std = md["variance.standardized"].to_numpy()

    dense = counts.toarray()
    sd = np.sqrt(np.where(expected > 0, expected, 1.0))
    z = (dense - mean[:, None]) / sd[:, None]
    unclipped = (z.max(axis=1) <= np.sqrt(n_cells)) & (expected > 0)
    assert unclipped.sum() > 10, "fixture clips too much to test the identity"

    np.testing.assert_allclose(
        var_std[unclipped], (variance / expected)[unclipped], rtol=1e-12,
        err_msg="variance.expected is not the LOESS fit the standardization used",
    )


def test_variable_feature_plot_reads_the_stored_stats(small_seurat):
    """The rename has to reach the consumers, or they degrade in silence.

    `variable_feature_plot` prefers the stored VST statistics and falls back to
    recomputing a raw dispersion from the data matrix when it cannot find them.
    A consumer left looking for `means` therefore keeps drawing a plot — a
    different, wrong plot, on a different y axis, with no error anywhere. The
    axis label is what distinguishes the two branches.
    """
    pytest.importorskip("matplotlib")
    from truecell.plotting import variable_feature_plot

    _fitted(small_seurat, "vst")
    fig = variable_feature_plot(small_seurat, label=False)
    try:
        assert fig.axes[0].get_ylabel() == "Standardized Variance", (
            "variable_feature_plot fell back to recomputing dispersions; it no "
            "longer finds the VST statistics find_variable_features stored"
        )
    finally:
        import matplotlib.pyplot as plt
        plt.close(fig)
