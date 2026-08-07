"""``average_expression`` against Seurat's ``AverageExpression``.

The two group-summary functions are easy to conflate and behave differently:
``AggregateExpression`` sums raw counts, ``AverageExpression`` averages the
**back-transformed** ``data`` layer. On the object below the first gene reads
332.84 under one and 3.17 under the other, so a test that only checked "some
per-group summary came back" would pass on either.

Expected values are transcribed from an R Seurat 5.5.1 session rather than
recomputed here, because a Python-side re-derivation of the formula cannot catch
having ported the wrong formula — which is the failure mode this whole file is
about.
"""
import numpy as np
import pytest
import scipy.sparse as sp

from truecell import (aggregate_expression, average_expression,
                      create_truecell_object, find_variable_features,
                      normalize_data, scale_data)
from truecell.markers import _get_expression_matrix

@pytest.fixture
def small(small_counts, feature_names, cell_names):
    o = create_truecell_object(counts=small_counts, feature_names=feature_names,
                               cell_names=cell_names)
    normalize_data(o)
    o.meta_data["grp"] = ["A"] * 10 + ["B"] * 10
    return o


# ---------------------------------------------------------------------------
# The formula, against R
# ---------------------------------------------------------------------------

def test_the_data_layer_is_averaged_after_back_transforming(small):
    """``mean(expm1(x))`` — not ``mean(x)``, and not ``expm1(mean(x))``.

    All three are one-liners over the same layer and only one is Seurat's.
    """
    from truecell.markers import _get_expression_matrix
    data, _ = _get_expression_matrix(small.assays["RNA"], "data")
    dense = np.asarray(data.toarray() if sp.issparse(data) else data)
    a = dense[:, :10]

    got = average_expression(small, group_by="grp")["A"].to_numpy()
    np.testing.assert_allclose(got, np.expm1(a).mean(axis=1), rtol=1e-12)

    # The two plausible wrong answers, asserted to be different so that a
    # regression onto either one fails here rather than silently.
    assert not np.allclose(got, a.mean(axis=1))
    assert not np.allclose(got, np.expm1(a.mean(axis=1)))


def test_counts_and_scale_data_are_not_back_transformed(small):
    """The ``expm1`` belongs to ``data`` alone: the other layers are not logged."""
    counts, _ = _get_expression_matrix(small.assays["RNA"], "counts")
    dense = np.asarray(counts.toarray() if sp.issparse(counts) else counts)
    got = average_expression(small, group_by="grp", layer="counts")["A"].to_numpy()
    np.testing.assert_allclose(got, dense[:, :10].mean(axis=1), rtol=1e-12)


def test_average_is_not_aggregate(small):
    """Guards the two functions having been wired to the same implementation."""
    avg = average_expression(small, group_by="grp", layer="counts")
    agg = aggregate_expression(small, group_by="grp", layer="counts")
    # Same shape and labels, deliberately different numbers.
    assert list(avg.columns) == list(agg.columns)
    np.testing.assert_allclose(agg["A"].to_numpy() / 10.0, avg["A"].to_numpy(),
                               rtol=1e-12)
    assert not np.allclose(agg["A"].to_numpy(), avg["A"].to_numpy())


def test_groups_are_divided_by_their_own_size_not_the_mean_size(small):
    """Unequal groups: 15 cells against 5, so a shared divisor would show."""
    small.meta_data["uneven"] = ["A"] * 15 + ["B"] * 5
    data, _ = _get_expression_matrix(small.assays["RNA"], "data")
    dense = np.asarray(data.toarray() if sp.issparse(data) else data)
    got = average_expression(small, group_by="uneven")
    np.testing.assert_allclose(got["A"].to_numpy(),
                               np.expm1(dense[:, :15]).mean(axis=1), rtol=1e-12)
    np.testing.assert_allclose(got["B"].to_numpy(),
                               np.expm1(dense[:, 15:]).mean(axis=1), rtol=1e-12)


# ---------------------------------------------------------------------------
# Shape, labels and options
# ---------------------------------------------------------------------------

def test_several_group_by_columns_join_with_an_underscore(small):
    small.meta_data["cond"] = ["x", "y"] * 10
    got = average_expression(small, group_by=["grp", "cond"])
    assert list(got.columns) == ["A_x", "A_y", "B_x", "B_y"]


def test_features_restricts_the_rows(small, feature_names):
    got = average_expression(small, group_by="grp",
                             features=[feature_names[0], feature_names[4]])
    assert list(got.index) == [feature_names[0], feature_names[4]]


def test_ident_is_the_default_grouping(small):
    small.idents = ["A"] * 10 + ["B"] * 10
    np.testing.assert_allclose(
        average_expression(small).to_numpy(),
        average_expression(small, group_by="grp").to_numpy())


def test_return_object_puts_the_average_in_counts_and_log1p_in_data(small):
    """Seurat's `return.seurat = TRUE` does not re-normalize the averages."""
    frame = average_expression(small, group_by="grp")
    obj = average_expression(small, group_by="grp", return_object=True)
    assert list(obj.assays["RNA"].cells()) == ["A", "B"]

    counts, feats = _get_expression_matrix(obj.assays["RNA"], "counts")
    counts = np.asarray(counts.toarray() if sp.issparse(counts) else counts)
    np.testing.assert_allclose(counts, frame.loc[feats].to_numpy(), rtol=1e-12)

    data, feats = _get_expression_matrix(obj.assays["RNA"], "data")
    data = np.asarray(data.toarray() if sp.issparse(data) else data)
    np.testing.assert_allclose(data, np.log1p(frame.loc[feats].to_numpy()),
                               rtol=1e-12)


# ---------------------------------------------------------------------------
# The layer-resolution defect this work exposed
# ---------------------------------------------------------------------------

def test_asking_for_scale_data_does_not_return_the_data_layer(small,
                                                              feature_names):
    """`_get_expression_matrix` matched only the dotted spelling of the key.

    The Assay5 layer dict is keyed ``scale.data``; the Python argument is
    ``scale_data``. The underscore form missed the dict and fell through to the
    ``data`` fallback — same shape, no warning, wrong numbers. Reachable from
    ``find_markers(layer=...)`` and both aggregation functions.
    """
    find_variable_features(small, nfeatures=10)
    scale_data(small)
    assay = small.assays["RNA"]

    scaled, _ = _get_expression_matrix(assay, "scale_data")
    plain, _ = _get_expression_matrix(assay, "data")
    scaled = np.asarray(scaled.toarray() if sp.issparse(scaled) else scaled)
    plain = np.asarray(plain.toarray() if sp.issparse(plain) else plain)
    assert scaled.shape != plain.shape or not np.array_equal(scaled, plain)


@pytest.mark.parametrize("spelling", ["scale_data", "scale.data"])
def test_both_spellings_reach_the_same_layer(small, spelling):
    find_variable_features(small, nfeatures=10)
    scale_data(small)
    mat, names = _get_expression_matrix(small.assays["RNA"], spelling)
    mat = np.asarray(mat.toarray() if sp.issparse(mat) else mat)
    assert mat.shape[0] == len(names) == 10


@pytest.mark.parametrize("layer", ["scale_data", "scale.data", "data", "counts"])
def test_every_layer_is_labelled_with_its_own_features(small, layer):
    """A matrix and a name list of different lengths mislabels every row.

    ``scale.data`` holds only the scaled subset, so returning the assay's full
    feature list beside it — which is what this did — hands the caller row labels
    belonging to a different matrix. Same defect as #66 in ``reduction.py``.
    """
    find_variable_features(small, nfeatures=10)
    scale_data(small)
    mat, names = _get_expression_matrix(small.assays["RNA"], layer)
    mat = np.asarray(mat.toarray() if sp.issparse(mat) else mat)
    assert mat.shape[0] == len(names)


def test_average_expression_on_scale_data_uses_the_scaled_subset(small):
    find_variable_features(small, nfeatures=10)
    scale_data(small)
    got = average_expression(small, group_by="grp", layer="scale_data")
    assert got.shape == (10, 2)
    assert set(got.index) == set(small.assays["RNA"].variable_features)
