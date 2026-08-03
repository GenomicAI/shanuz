"""Guards for `find_markers`' pre-filter-then-densify path.

`find_markers` used to build a dense (all genes × cells) array per group before
it had decided which genes it was going to test. It now runs both pre-filters —
`min_pct` and `logfc_threshold` — on the sparse matrices and densifies only the
genes that survive, which on PBMC 3k is ~1.6k rows of 13.7k.

That split introduced two things worth pinning, neither of which the suite could
see when the change was made:

**A row remap.** The per-gene loops index the dense array by position in
`test_indices`, not by feature index. Getting that wrong attaches every
statistic to the wrong gene — a silent, plausible-looking table. Mutating the
loop index to a constant left the whole suite green for `roc`.

**A dense fallback.** Both `_row_pct_positive` and `_row_expm1_sum` branch on
sparsity, and nothing anywhere ran `find_markers` against a dense layer. Every
mutation of those dense branches survived the full suite.
"""
import copy

import numpy as np
import pytest
import scipy.sparse as sp

from truecell import create_truecell_object
from truecell.markers import _roc_auc, find_markers
from truecell.preprocessing import normalize_data

TESTS = ["wilcox", "t", "bimod", "LR", "mast", "negbinom", "roc"]


@pytest.fixture
def graded_object():
    """30 genes × 60 cells, with a graded split so the tested set is a subset.

    Genes 0-9 are strongly up in group A, 10-19 mildly, 20-29 not at all — so
    `min_pct` and `logfc_threshold` both bite and `test_indices` is a genuine
    subset of the features. A pre-filter that kept everything, or a remap that
    read the wrong row, would have nothing to disagree with otherwise.
    """
    rng = np.random.default_rng(4)
    n1 = n2 = 30
    counts = rng.poisson(2.0, size=(30, n1 + n2)).astype(float)
    counts[:10, :n1] += rng.poisson(20.0, size=(10, n1))
    counts[10:20, :n1] += rng.poisson(3.0, size=(10, n1))
    counts[25:, :] = 0.0  # dropped by min_pct in both groups
    obj = create_truecell_object(
        sp.csc_matrix(counts), assay="RNA",
        feature_names=[f"g{i}" for i in range(30)],
        cell_names=[f"c{i}" for i in range(n1 + n2)],
    )
    normalize_data(obj)
    obj.idents = ["A"] * n1 + ["B"] * n2
    return obj


def _densified(obj):
    """The same object with every layer held as a dense array."""
    dense = copy.deepcopy(obj)
    layers = dense.assays["RNA"].layers
    for name, mat in list(layers.items()):
        layers[name] = np.asarray(mat.toarray())
    return dense


# ---------------------------------------------------------------------------
# The dense fallback has to agree with the sparse path
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("test_use", TESTS)
def test_a_dense_layer_gives_the_same_table_as_a_sparse_one(graded_object, test_use):
    sparse_res = find_markers(graded_object, "A", "B", test_use=test_use)
    dense_obj = _densified(graded_object)
    assert not sp.issparse(dense_obj.assays["RNA"].layers["data"]), (
        "fixture is not exercising the dense branch"
    )
    dense_res = find_markers(dense_obj, "A", "B", test_use=test_use)

    assert list(dense_res.index) == list(sparse_res.index)
    for col in sparse_res.columns:
        np.testing.assert_allclose(
            dense_res[col].to_numpy(dtype=float),
            sparse_res[col].to_numpy(dtype=float),
            rtol=1e-12, atol=0,
            err_msg=f"{test_use}: column {col!r} differs between layers",
        )


def test_the_pre_filter_actually_drops_genes(graded_object):
    """If nothing is filtered, the remap tests below prove nothing."""
    res = find_markers(graded_object, "A", "B")
    n_features = len(graded_object.assays["RNA"]._all_feature_names)
    assert 0 < len(res) < n_features


# ---------------------------------------------------------------------------
# Every statistic has to belong to the gene it is labelled with
# ---------------------------------------------------------------------------

def _layer_rows(obj, gene, group):
    """One gene's expression in one group, straight off the layer."""
    features = list(obj.assays["RNA"]._all_feature_names)
    cells = list(obj.cell_names())
    row = features.index(gene)
    cols = [i for i, (c, g) in enumerate(zip(cells, obj.idents)) if g == group]
    return np.asarray(obj.assays["RNA"].layers["data"][row, cols].todense()).ravel()


def test_roc_statistics_belong_to_their_own_gene(graded_object):
    """`myAUC` is recomputed per row from the layer and must match.

    The ROC branch reads `mat1[i]` where `i` counts along `test_indices`. Nothing
    else in the suite would notice if that index were wrong, because a shifted
    table is still a well-formed table with plausible numbers in it.
    """
    res = find_markers(graded_object, "A", "B", test_use="roc")
    assert len(res) > 1, "need several rows for a mis-mapping to show"
    for gene in res.index:
        auc, power = _roc_auc(_layer_rows(graded_object, gene, "A"),
                              _layer_rows(graded_object, gene, "B"))
        assert res.loc[gene, "myAUC"] == pytest.approx(auc, rel=1e-12)
        assert res.loc[gene, "power"] == pytest.approx(power, rel=1e-12)


def test_restricting_features_keeps_each_gene_s_own_statistics(graded_object):
    """`features=` narrows `test_indices`, which is what the dense rows index by.

    Nothing in the suite passed `features=` to `find_markers`, so a restriction
    that shifted the rows relative to the labels would have gone unseen.
    """
    full = find_markers(graded_object, "A", "B")
    assert len(full) > 3
    wanted = list(full.index[::2])

    subset = find_markers(graded_object, "A", "B", features=wanted)
    assert list(subset.index) == wanted
    for col in full.columns:
        np.testing.assert_allclose(
            subset[col].to_numpy(dtype=float),
            full.loc[wanted, col].to_numpy(dtype=float),
            rtol=1e-12, atol=0, err_msg=f"column {col!r} shifted under features=",
        )


def test_wilcox_p_values_belong_to_their_own_gene(graded_object):
    """The same check for the default test, recomputed with scipy directly."""
    from scipy.stats import mannwhitneyu

    res = find_markers(graded_object, "A", "B", test_use="wilcox")
    assert len(res) > 1
    for gene in res.index:
        _, expected = mannwhitneyu(
            _layer_rows(graded_object, gene, "A"),
            _layer_rows(graded_object, gene, "B"),
            alternative="two-sided", use_continuity=True, method="asymptotic",
        )
        assert res.loc[gene, "p_val"] == pytest.approx(expected, rel=1e-12)
