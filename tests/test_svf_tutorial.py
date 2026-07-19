"""Tests for the spatial-statistics tutorial's helpers.

These cover the comparison instrument, not the science. A tutorial whose
comparator quietly agrees with everything prints a parity table that means
nothing, so `compare_anchors` and its tolerance lookup get tested the same way
any other code would.

The tutorial's own end-to-end run needs the Xenium download and lives behind
SHANUZ_TUTORIAL_SMOKE=1 in test_tutorial_smoke.py.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from tutorials.xenium_svf_tutorial import (  # noqa: E402
    FLOAT_TOLERANCES,
    _tolerance_for,
    compare_anchors,
    digest,
    name_anchor,
)


# ---------------------------------------------------------------------------
# digest / name_anchor
# ---------------------------------------------------------------------------

def test_digest_is_order_sensitive():
    """The whole point: a reordered list is a different object."""
    assert digest(["a", "b", "c"]) != digest(["a", "c", "b"])


def test_digest_is_stable_and_short():
    assert digest(["a", "b"]) == digest(["a", "b"])
    assert len(digest(["a", "b"])) == 12


def test_digest_separates_on_newlines_not_concatenation():
    """`ab, c` and `a, bc` must not collide."""
    assert digest(["ab", "c"]) != digest(["a", "bc"])


def test_name_anchor_records_ends_and_length():
    a = name_anchor([f"c{i}" for i in range(10)])
    assert a["n"] == 10
    assert a["head"] == ["c0", "c1", "c2"]
    assert a["tail"] == ["c7", "c8", "c9"]


def test_name_anchor_head_and_tail_overlap_when_short():
    a = name_anchor(["x", "y"])
    assert a["n"] == 2
    assert a["head"] == ["x", "y"]
    assert a["tail"] == ["x", "y"]


# ---------------------------------------------------------------------------
# compare_anchors
# ---------------------------------------------------------------------------

def test_identical_anchors_all_match():
    a = {"x": 1, "y": ["p", "q"], "nested": {"z": "s"}}
    table = compare_anchors(a, dict(a))
    assert table["match"].all()
    assert len(table) == 3


def test_a_difference_is_reported():
    table = compare_anchors({"x": 1}, {"x": 2})
    assert not table["match"].any()
    assert table.iloc[0]["python"] == 1
    assert table.iloc[0]["r"] == 2


def test_missing_on_either_side_is_a_mismatch_not_a_skip():
    """A field only one tool produces must not silently vanish from the table."""
    table = compare_anchors({"only_py": 1}, {"only_r": 2})
    assert len(table) == 2
    assert not table["match"].any()


def test_nested_fields_are_dotted():
    table = compare_anchors({"a": {"b": 1}}, {"a": {"b": 1}})
    assert table.iloc[0]["field"] == "a.b"


def test_list_order_matters():
    assert not compare_anchors({"x": ["a", "b"]}, {"x": ["b", "a"]})["match"].any()


def test_int_and_string_forms_of_a_name_compare_equal():
    """R writes cell ids as strings, JSON round-trips some as numbers."""
    assert compare_anchors({"x": ["1", "2"]}, {"x": [1, 2]})["match"].all()


# ---------------------------------------------------------------------------
# tolerances
# ---------------------------------------------------------------------------

def test_tolerance_is_looked_up_by_exact_field_name():
    assert _tolerance_for("toy.auto_radius") == FLOAT_TOLERANCES["toy.auto_radius"]
    # A prefix must not grant slack to a field nobody chose.
    assert _tolerance_for("toy.auto_radius_extra") is None
    assert _tolerance_for("toy") is None


def test_float_anchors_use_their_tolerance():
    table = compare_anchors({"toy": {"auto_radius": 0.165}},
                            {"toy": {"auto_radius": 0.165 + 1e-15}})
    assert table["match"].all()


def test_a_float_anchor_still_fails_outside_its_tolerance():
    table = compare_anchors({"toy": {"auto_radius": 0.165}},
                            {"toy": {"auto_radius": 0.166}})
    assert not table["match"].any()


def test_untoleranced_floats_are_compared_exactly():
    """Anything not named in FLOAT_TOLERANCES gets no slack at all."""
    table = compare_anchors({"container": {"n_cells": 1.0}},
                            {"container": {"n_cells": 1.0 + 1e-12}})
    assert not table["match"].any()


def test_every_declared_tolerance_names_a_real_anchor():
    """Guards against a tolerance outliving the anchor it was written for.

    A stale entry is worse than no entry: it looks like a considered decision.
    """
    from tutorials.xenium_svf_tutorial import toy_anchors

    produced = {f"toy.{k}" for k in toy_anchors()}
    declared_toy = {k for k in FLOAT_TOLERANCES if k.startswith("toy.")}
    assert declared_toy <= produced


def test_nan_never_matches():
    table = compare_anchors({"toy": {"auto_radius": float("nan")}},
                            {"toy": {"auto_radius": float("nan")}})
    assert not table["match"].any()
