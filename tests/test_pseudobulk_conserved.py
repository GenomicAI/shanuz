"""Tests for AggregateExpression (pseudobulk) and FindConservedMarkers."""
import numpy as np
import pandas as pd
import scipy.sparse as sp
import pytest

from shanuz import aggregate_expression, find_conserved_markers, find_markers
from shanuz.preprocessing import normalize_data


# ---------------------------------------------------------------------------
# aggregate_expression
# ---------------------------------------------------------------------------

def test_aggregate_sums_counts_per_group(small_seurat, small_counts):
    small_seurat.meta_data["donor"] = ["d1", "d2"] * 10  # d1 = even, d2 = odd

    agg = aggregate_expression(small_seurat, group_by="donor")

    assert list(agg.columns) == ["d1", "d2"]
    assert agg.shape == (50, 2)

    dense = small_counts.toarray()
    even = dense[:, ::2].sum(axis=1)
    odd = dense[:, 1::2].sum(axis=1)
    np.testing.assert_allclose(agg["d1"].to_numpy(), even)
    np.testing.assert_allclose(agg["d2"].to_numpy(), odd)


def test_aggregate_by_ident(small_seurat):
    small_seurat.idents = ["A"] * 10 + ["B"] * 10
    agg = aggregate_expression(small_seurat, group_by="ident")
    assert set(agg.columns) == {"A", "B"}


def test_aggregate_multiple_group_by_joins_with_underscore(small_seurat):
    small_seurat.meta_data["donor"] = ["d1"] * 10 + ["d2"] * 10
    small_seurat.meta_data["cond"] = ["x", "y"] * 10
    agg = aggregate_expression(small_seurat, group_by=["donor", "cond"])
    assert set(agg.columns) == {"d1_x", "d1_y", "d2_x", "d2_y"}


def test_aggregate_return_object(small_seurat, small_counts):
    small_seurat.meta_data["donor"] = ["d1", "d2"] * 10
    obj = aggregate_expression(small_seurat, group_by="donor", return_object=True)

    assert len(obj) == 2  # one "cell" per group
    assert obj.cell_names() == ["d1", "d2"]
    assert "donor" in obj.meta_data.columns

    df = aggregate_expression(small_seurat, group_by="donor")
    pb = obj.assays["RNA"].layers["counts"]
    pb = pb.toarray() if sp.issparse(pb) else np.asarray(pb)
    np.testing.assert_allclose(pb, df.to_numpy())


def test_aggregate_features_subset(small_seurat):
    small_seurat.meta_data["donor"] = ["d1", "d2"] * 10
    genes = ["gene_0", "gene_5", "gene_9"]
    agg = aggregate_expression(small_seurat, group_by="donor", features=genes)
    assert list(agg.index) == genes


def test_aggregate_unknown_group_by_raises(small_seurat):
    with pytest.raises(KeyError):
        aggregate_expression(small_seurat, group_by="nope")


# ---------------------------------------------------------------------------
# find_conserved_markers
# ---------------------------------------------------------------------------

def test_conserved_markers_basic(small_seurat):
    normalize_data(small_seurat)
    small_seurat.idents = ["A"] * 10 + ["B"] * 10
    # Each batch holds both idents (5 A + 5 B).
    small_seurat.meta_data["batch"] = (["x", "y"] * 5) + (["x", "y"] * 5)

    res = find_conserved_markers(
        small_seurat, ident_1="A", grouping_var="batch",
        min_pct=0.0, logfc_threshold=0.0,
    )

    for pref in ("x_", "y_"):
        assert f"{pref}p_val" in res.columns
        assert f"{pref}avg_log2FC" in res.columns
    assert "max_pval" in res.columns and "combined_p_val" in res.columns

    # Genes must be markers in BOTH levels.
    x_only = find_markers(small_seurat.subset(
        cells=[c for c, b in zip(small_seurat.cell_names(),
                                 small_seurat.meta_data["batch"]) if b == "x"]),
        ident_1="A", min_pct=0.0, logfc_threshold=0.0)
    assert set(res.index).issubset(set(x_only.index))

    # combined_p_val is a valid probability and sorted ascending.
    assert ((res["combined_p_val"] >= 0) & (res["combined_p_val"] <= 1)).all()
    assert res["combined_p_val"].is_monotonic_increasing
    # max_pval is the worst of the per-level p-values.
    np.testing.assert_allclose(
        res["max_pval"].to_numpy(),
        res[["x_p_val", "y_p_val"]].max(axis=1).to_numpy(),
    )


def test_conserved_markers_skips_level_without_comparison_group(small_seurat):
    normalize_data(small_seurat)
    small_seurat.idents = ["A"] * 10 + ["B"] * 10
    # "lonely" contains only ident A -> no comparison group -> skipped.
    small_seurat.meta_data["cond"] = ["lonely"] * 5 + ["together"] * 15

    with pytest.warns(RuntimeWarning, match="lonely"):
        res = find_conserved_markers(
            small_seurat, ident_1="A", grouping_var="cond",
            min_pct=0.0, logfc_threshold=0.0,
        )

    # Only the "together" level survives; combined == that level's p_val.
    assert any(c.startswith("together_") for c in res.columns)
    assert not any(c.startswith("lonely_") for c in res.columns)
    np.testing.assert_allclose(
        res["combined_p_val"].to_numpy(), res["together_p_val"].to_numpy()
    )


def test_conserved_markers_unknown_grouping_var_raises(small_seurat):
    normalize_data(small_seurat)
    small_seurat.idents = ["A"] * 10 + ["B"] * 10
    with pytest.raises(KeyError):
        find_conserved_markers(small_seurat, ident_1="A", grouping_var="missing")
