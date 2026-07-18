"""Tests for the object-internals tutorial's helpers and pipeline.

Covers ``tutorials/pbmc3k_objects_tutorial.py``:

  * the anchor vocabulary (``digest``, ``name_anchor``, ``matrix_anchor``,
    ``r_safe_names``, ``ident_anchor``) on inputs whose right answer can be
    written down by hand;
  * ``compare_anchors``, which is the tutorial's actual instrument — if it
    reports a match where there is none, every number the tutorial prints is
    worthless, so it is tested against both agreement and each kind of
    disagreement;
  * the layer round trip and the pipeline wiring on a small synthetic object,
    so they are exercised without the 2,700-cell download.

The real-data numbers live in ``tests/test_tutorial_smoke.py`` and only run
when ``SHANUZ_TUTORIAL_SMOKE=1`` is set. Network-free.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz import create_assay5_object  # noqa: E402
from tutorials.pbmc3k_objects_tutorial import (  # noqa: E402
    _column_sum,
    assign_idents,
    batch_labels,
    compare_anchors,
    digest,
    ident_anchor,
    join_all_layers,
    layer_anchor,
    matrix_anchor,
    name_anchor,
    r_safe_names,
    split_join_roundtrip,
)


# ----------------------------------------------------------------------
# The anchor vocabulary
# ----------------------------------------------------------------------


def test_digest_is_order_sensitive():
    # Half of what the tutorial checks is that R and Python hold their cells in
    # the same order, so a set-like fingerprint would defeat the purpose.
    assert digest(["a", "b"]) != digest(["b", "a"])


def test_digest_is_stable_for_the_same_sequence():
    assert digest(["a", "b", "c"]) == digest(["a", "b", "c"])


def test_name_anchor_reports_length_head_and_tail():
    got = name_anchor([f"c{i}" for i in range(10)])
    assert got["n"] == 10
    assert got["head"] == ["c0", "c1", "c2"]
    assert got["tail"] == ["c7", "c8", "c9"]


def test_name_anchor_handles_sequences_shorter_than_the_window():
    got = name_anchor(["only"])
    assert got["n"] == 1
    assert got["head"] == ["only"]


def test_matrix_anchor_counts_true_nonzeros_not_stored_entries():
    # An explicitly stored zero must not count as an edge; this is what makes
    # nnz comparable against R, which stores a different set of them.
    mat = sp.csc_matrix(np.array([[1.0, 0.0], [0.0, 2.0]]))
    mat[0, 1] = 0.0  # stored explicitly
    got = matrix_anchor(mat)
    assert got["nnz"] == 2
    assert got["sum"] == pytest.approx(3.0)
    assert got["shape"] == [2, 2]


def test_r_safe_names_matches_read10x_mangling():
    # Read10X rewrites '_' to '-'; shanuz's loader does not, and the difference
    # belongs to the two loaders rather than to the object model.
    assert r_safe_names(["Y_RNA", "RP11-1_2", "CD3E"]) == ["Y-RNA", "RP11-1-2", "CD3E"]


def test_ident_anchor_counts_each_level():
    got = ident_anchor(["T", "B", "T", "Mono"])
    assert got["levels"] == ["B", "Mono", "T"]
    assert got["counts"] == {"B": 1, "Mono": 1, "T": 2}


def test_batch_labels_alternate_deterministically():
    # Deterministic on purpose: Seurat's vignette uses sample(), which would
    # put two different partitions on the two sides of the comparison.
    assert batch_labels(4) == ["batch1", "batch2", "batch1", "batch2"]


def test_column_sum_reports_a_non_numeric_column_instead_of_raising():
    # FetchData is under test, so the anchor must survive it returning junk —
    # this is the shape the original defect produced.
    frame = pd.DataFrame({"x": [sp.csc_matrix((1, 1))] * 2})
    assert "non-numeric" in _column_sum(frame, "x")


def test_column_sum_reports_a_missing_column():
    assert _column_sum(pd.DataFrame({"a": [1]}), "b") == "missing"


# ----------------------------------------------------------------------
# compare_anchors — the instrument itself
# ----------------------------------------------------------------------


def test_compare_anchors_matches_identical_trees():
    tree = {"a": 1, "b": {"c": "x", "d": [1, 2]}}
    result = compare_anchors(tree, dict(tree))
    assert result["match"].all()
    assert len(result) == 3  # leaves, not branches


def test_compare_anchors_flags_a_differing_leaf():
    result = compare_anchors({"a": 1, "b": 2}, {"a": 1, "b": 99}).set_index("field")
    assert not result.loc["b", "match"]
    assert result.loc["a", "match"]


def test_compare_anchors_flags_a_field_present_on_only_one_side():
    result = compare_anchors({"a": 1}, {"a": 1, "extra": 2})
    row = result.set_index("field").loc["extra"]
    assert not row["match"]


def test_compare_anchors_is_order_sensitive_by_default():
    # Layer names and cell orders are the whole point; a set comparison here
    # would silently accept a permuted object.
    result = compare_anchors({"layers": ["a", "b"]}, {"layers": ["b", "a"]})
    assert not result["match"].all()


def test_compare_anchors_honours_the_named_float_tolerance():
    # `reductions.pca_stdev_head` is one of two fields allowed to differ, and
    # only because the two PCAs use different solvers.
    py = {"reductions": {"pca_stdev_head": [1.0, 2.0]}}
    r = {"reductions": {"pca_stdev_head": [1.0001, 2.0002]}}
    assert compare_anchors(py, r)["match"].all()


def test_compare_anchors_still_rejects_a_gross_float_difference():
    py = {"reductions": {"pca_stdev_head": [1.0, 2.0]}}
    r = {"reductions": {"pca_stdev_head": [1.0, 5.0]}}
    assert not compare_anchors(py, r)["match"].all()


def test_compare_anchors_applies_no_tolerance_to_unlisted_floats():
    # The exception list is short on purpose — an across-the-board tolerance
    # would buy back all the sensitivity that makes this tutorial worth running.
    assert not compare_anchors({"x": 1.0}, {"x": 1.0001})["match"].all()


def test_compare_anchors_ignores_order_only_where_listed():
    py = {"join_all_layers": {"layers_after_join": ["data", "scale.data", "counts"]}}
    r = {"join_all_layers": {"layers_after_join": ["data", "counts", "scale.data"]}}
    assert compare_anchors(py, r)["match"].all()


def test_unordered_comparison_still_checks_membership():
    py = {"join_all_layers": {"layers_after_join": ["data", "counts"]}}
    r = {"join_all_layers": {"layers_after_join": ["data", "scale.data"]}}
    assert not compare_anchors(py, r)["match"].all()


# ----------------------------------------------------------------------
# The layer round trip
# ----------------------------------------------------------------------


def _assay(n_cells=6):
    mat = sp.csc_matrix(
        np.arange(1, 4 * n_cells + 1, dtype=float).reshape(4, n_cells))
    return create_assay5_object(
        counts=mat,
        feature_names=[f"g{i}" for i in range(4)],
        cell_names=[f"c{i}" for i in range(n_cells)],
        key="rna_",
    )


def test_split_join_roundtrip_reports_a_clean_round_trip():
    got = split_join_roundtrip(_assay(), batch_labels(6), layer="counts")
    assert got["layers_after_split"] == ["counts.batch1", "counts.batch2"]
    assert got["layers_after_join"] == ["counts"]
    assert got["layer_name_restored"]
    assert got["cell_order_restored"]
    assert got["matrix_restored"]


def test_split_join_roundtrip_notices_a_permuted_result(monkeypatch):
    # The helper must be able to say "no" — otherwise the three booleans it
    # reports are decoration. Force a join that returns the split order.
    from shanuz.assay5 import StdAssay

    original = StdAssay.join_layers

    def permuting_join(self, layers=None):
        joined = original(self, layers)
        name = joined.layers_list()[0]
        mat = joined.layers[name]
        joined.layers[name] = mat[:, ::-1]
        joined._layer_cells[name] = list(reversed(joined._layer_cells[name]))
        return joined

    monkeypatch.setattr(StdAssay, "join_layers", permuting_join)
    got = split_join_roundtrip(_assay(), batch_labels(6), layer="counts")
    assert not got["cell_order_restored"]
    assert not got["matrix_restored"]


def test_join_all_layers_survives_a_prepared_assay():
    # The no-argument call on an assay carrying layers of differing feature
    # counts — the call every real script makes, and the one that used to raise.
    assay = _assay()
    assay._add_layer("scale.data", np.ones((2, 6)),
                     feature_names=["g0", "g1"], cell_names=assay.cells())
    got = join_all_layers(assay, batch_labels(6), layer="counts")
    assert got["error"] is None
    assert set(got["layers_after_join"]) == {"counts", "scale.data"}


def test_join_all_layers_records_an_error_instead_of_raising(monkeypatch):
    # This is how the defect presented, and the helper has to survive it: a
    # tutorial that dies on the finding cannot report the finding.
    from shanuz.assay5 import StdAssay

    def exploding_join(self, layers=None):
        raise ValueError("incompatible dimensions for axis 0")

    monkeypatch.setattr(StdAssay, "join_layers", exploding_join)
    got = join_all_layers(_assay(), batch_labels(6), layer="counts")
    assert got["error"].startswith("ValueError")
    assert got["layers_after_join"] is None


# ----------------------------------------------------------------------
# Identity assignment
# ----------------------------------------------------------------------


def test_assign_idents_applies_the_first_matching_gate():
    from shanuz import create_shanuz_object, normalize_data

    # Cell 0 expresses CD3E and LYZ — "T" wins, being the earlier gate.
    # Cell 1 expresses only LYZ, cell 2 nothing.
    counts = sp.csc_matrix(np.array([
        [5.0, 0.0, 0.0],   # CD3E
        [0.0, 0.0, 0.0],   # MS4A1
        [7.0, 9.0, 0.0],   # LYZ
    ]))
    obj = create_shanuz_object(
        counts=counts, assay="RNA", min_cells=0, min_features=0,
        feature_names=["CD3E", "MS4A1", "LYZ"],
        cell_names=["c0", "c1", "c2"], project="P",
    )
    normalize_data(obj)
    assert list(assign_idents(obj)) == ["T", "Mono", "Other"]


def test_layer_anchor_describes_every_layer():
    assay = _assay()
    got = layer_anchor(assay)
    assert set(got) == {"counts"}
    assert got["counts"]["shape"] == [4, 6]
    assert got["counts"]["cells"]["n"] == 6
