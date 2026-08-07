"""``aggregate_expression`` against Seurat's ``AggregateExpression``.

``tests/test_pseudobulk_conserved.py`` already exercises this function, but every
assertion there checks Python against a Python re-derivation of the same formula
(``agg["d1"] == dense[:, ::2].sum(axis=1)``). That shape of test cannot detect
having ported the wrong *convention*, which is exactly how the CLR defect
survived its own unit test — the kernel was right and the axis was not.

This file pins values transcribed from an R Seurat 5.5.1 session instead. It
found one defect on the first run: ``return_object=True`` left the raw sums in
the ``data`` layer, where Seurat normalizes them.

The counts matrix is hardcoded rather than drawn, so neither language's RNG is in
the loop and the two sides are provably looking at the same numbers.
"""
import numpy as np
import pytest
import scipy.sparse as sp

from truecell import (aggregate_expression, average_expression,
                      create_truecell_object, normalize_data)
from truecell.markers import _get_expression_matrix

# 8 genes x 9 cells, column-major to match R's matrix() fill order.
COUNTS = np.array(
    [5, 0, 3, 3, 7, 9, 3, 5,
     2, 4, 7, 6, 8, 8, 1, 6,
     7, 7, 8, 1, 5, 9, 8, 9,
     4, 3, 0, 3, 5, 0, 2, 3,
     8, 1, 3, 3, 3, 7, 0, 1,
     9, 9, 0, 4, 7, 3, 2, 7,
     2, 0, 0, 4, 5, 5, 6, 8,
     4, 1, 4, 9, 8, 1, 1, 7,
     9, 9, 3, 6, 7, 2, 0, 3], dtype=float).reshape((8, 9), order="F")

GRP = ["A", "A", "A", "B", "B", "B", "C", "C", "C"]
COND = ["x", "y", "x", "y", "x", "y", "x", "y", "x"]

# AggregateExpression(o, group.by = "grp")
R_BY_GRP = np.array([
    [14, 21, 15], [11, 13, 10], [18, 3, 7], [10, 10, 19],
    [20, 15, 20], [26, 10, 8], [12, 4, 7], [20, 11, 18]], dtype=float)

# AggregateExpression(o, group.by = c("grp", "cond"))
R_BY_GRP_COND = np.array([
    [12, 2, 8, 13, 11, 4], [7, 4, 1, 12, 9, 1], [11, 7, 3, 0, 3, 4],
    [4, 6, 3, 7, 10, 9], [12, 8, 3, 12, 12, 8], [18, 8, 7, 3, 7, 1],
    [11, 1, 0, 4, 6, 1], [14, 6, 1, 10, 11, 7]], dtype=float)

# GetAssayData(AggregateExpression(..., return.seurat = TRUE), layer = "data")
R_RETURN_SEURAT_DATA = np.array([
    [6.97513565517, 7.78936889097, 7.27469276703],
    [6.73422852209, 7.31005061772, 6.86957402540],
    [7.22624231975, 5.84594034510, 6.51334423359],
    [6.63903728447, 7.04788696809, 7.51093567067],
    [7.33153010791, 7.45306228629, 7.56220161124],
    [7.59374330606, 7.04788696809, 6.64669017025],
    [6.82114076979, 6.13289925255, 6.51334423359],
    [7.33153010791, 7.14311812261, 7.45689884166]])


@pytest.fixture
def obj():
    o = create_truecell_object(
        counts=sp.csc_matrix(COUNTS),
        feature_names=[f"g{i}" for i in range(1, 9)],
        cell_names=[f"c{i}" for i in range(1, 10)])
    normalize_data(o)
    o.meta_data["grp"] = GRP
    o.meta_data["cond"] = COND
    return o


def _layer(obj, layer):
    mat, feats = _get_expression_matrix(obj.assays["RNA"], layer)
    mat = np.asarray(mat.toarray() if sp.issparse(mat) else mat)
    return mat, feats


def test_sums_match_seurat(obj):
    got = aggregate_expression(obj, group_by="grp")
    assert list(got.columns) == ["A", "B", "C"]
    np.testing.assert_array_equal(got.to_numpy(), R_BY_GRP)


def test_multi_column_grouping_matches_seurat(obj):
    got = aggregate_expression(obj, group_by=["grp", "cond"])
    assert list(got.columns) == ["A_x", "A_y", "B_x", "B_y", "C_x", "C_y"]
    np.testing.assert_array_equal(got.to_numpy(), R_BY_GRP_COND)


def test_return_object_normalizes_the_data_layer(obj):
    """The defect this file found.

    Seurat's ``return.seurat = TRUE`` runs ``NormalizeData`` over the pseudobulk,
    so ``data`` is ``log1p(sums / colSums × 10000)``. Truecell left the raw sums
    there, which every downstream reader of that layer would have taken for
    normalized expression.
    """
    out = aggregate_expression(obj, group_by="grp", return_object=True)

    counts, feats = _layer(out, "counts")
    np.testing.assert_array_equal(counts, R_BY_GRP)

    data, feats = _layer(out, "data")
    np.testing.assert_allclose(data, R_RETURN_SEURAT_DATA, rtol=1e-10)

    # And specifically not the two plausible wrong answers.
    assert not np.allclose(data, R_BY_GRP)
    assert not np.allclose(data, np.log1p(R_BY_GRP))


def test_normalization_can_be_switched_off(obj):
    out = aggregate_expression(obj, group_by="grp", return_object=True,
                               normalization_method=None)
    data, _ = _layer(out, "data")
    np.testing.assert_array_equal(data, R_BY_GRP)


def test_scale_factor_is_honoured(obj):
    out = aggregate_expression(obj, group_by="grp", return_object=True,
                               scale_factor=1e6)
    data, _ = _layer(out, "data")
    expected = np.log1p(R_BY_GRP / R_BY_GRP.sum(axis=0) * 1e6)
    np.testing.assert_allclose(data, expected, rtol=1e-10)


def test_the_two_functions_normalize_their_objects_differently(obj):
    """The one place ``AggregateExpression`` and ``AverageExpression`` diverge.

    Aggregate runs a full library-size normalization; Average takes a plain
    ``log1p`` of the means. Wiring both to the same helper would pass every other
    test in this file.
    """
    agg = aggregate_expression(obj, group_by="grp", return_object=True)
    avg_frame = average_expression(obj, group_by="grp")
    avg = average_expression(obj, group_by="grp", return_object=True)

    agg_data, _ = _layer(agg, "data")
    avg_data, feats = _layer(avg, "data")

    np.testing.assert_allclose(agg_data, R_RETURN_SEURAT_DATA, rtol=1e-10)
    np.testing.assert_allclose(avg_data,
                               np.log1p(avg_frame.loc[feats].to_numpy()),
                               rtol=1e-10)
    assert not np.allclose(agg_data, avg_data)
