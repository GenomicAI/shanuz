"""Which features a reduction ran on, and which names label its rows.

`_get_scaled_data` filtered the requested features down to the ones the layer
carried and returned only the matrix. Every caller then labelled that matrix with
its own *request*, so a dropped feature shifted every row below it onto the wrong
gene: `run_pca` stored 2,000 names against 1,980 loadings, `viz_dim_loadings`
drew the wrong gene on every bar, and `jack_straw` handed each p-value to its
neighbour. Nothing raised, and the numbers were all individually correct.

Seurat's `PrepDR5` does the same subsetting but warns naming the excluded
features, and its loadings take their rownames from the subset matrix, so the
labels cannot drift from the rows.

The v3 path had a second problem underneath: `Assay.scale_data` is a bare ndarray
holding only the scaled subset, with no record of which subset. Four places
guessed differently — `features("scale_data")` took the first n features,
`subset()` threw the layer away, `_scaled_feature_names` returned the full list,
and `_get_scaled_data` indexed it by position in the assay. The last one either
raised `IndexError` or read whichever gene happened to sit at that row.
"""
import warnings

import numpy as np
import pytest
import scipy.sparse as sp

from shanuz import (
    create_shanuz_object,
    find_variable_features,
    normalize_data,
    run_ica,
    run_pca,
    scale_data,
)
from shanuz.assay import Assay
from shanuz.jackstraw import _scaled_matrix_for_reduction

N_GENES, N_CELLS, N_VAR = 60, 40, 20


def _obj(use_v5=True, seed=0, tag="c"):
    rng = np.random.default_rng(seed)
    counts = sp.csc_matrix(rng.poisson(3.0, size=(N_GENES, N_CELLS)).astype(float))
    obj = create_shanuz_object(
        counts=counts,
        assay="RNA",
        feature_names=[f"g{i}" for i in range(N_GENES)],
        cell_names=[f"{tag}{j}" for j in range(N_CELLS)],
        use_v5=use_v5,
    )
    normalize_data(obj)
    find_variable_features(obj, nfeatures=N_VAR)
    scale_data(obj)
    return obj


def _interleaved(obj):
    """(requested, present) where the absent features sit *between* the present ones.

    Interleaving is the point. A run that appends the dropped features to the end
    of the list would still label row 0 correctly; scattering them means the first
    wrong label lands on the very first row.

    `present` is reversed out of the layer's own order so that "the features the
    caller asked for" and "the order the layer stores them in" are different
    lists. A filter that walks the layer instead of the request returns the same
    genes and the same count, and only the ordering gives it away.
    """
    assay = obj.assays["RNA"]
    scaled = list(assay.features("scale_data")) if isinstance(assay, Assay) else list(
        assay._scaled_features
    )
    absent = [f for f in assay.features() if f not in set(scaled)]
    present = scaled[:10][::-1]
    asked = [x for pair in zip(absent[:10], present) for x in pair]
    return asked, present


# ---------------------------------------------------------------------------
# 1. The drop is announced, and it names what it dropped
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("use_v5", [True, False])
def test_run_pca_warns_and_names_the_features_it_could_not_use(use_v5):
    obj = _obj(use_v5=use_v5)
    asked, _ = _interleaved(obj)

    with pytest.warns(RuntimeWarning) as record:
        run_pca(obj, n_pcs=3, features=asked)

    message = str(record[0].message)
    assert "10 features" in message, "the count is the part that says whether to worry"
    assert "scale.data" in message, "the layer names where to look"
    dropped = [f for f in asked if f not in set(obj.reductions["pca"].features())]
    assert all(f in message for f in dropped), (
        "a warning that does not name the genes leaves the user to diff two "
        "lists by hand, which is exactly the work the warning exists to save"
    )


def test_no_warning_when_every_requested_feature_is_present():
    """The guard must not cry wolf on the ordinary path."""
    obj = _obj()
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        run_pca(obj, n_pcs=3)


def test_requesting_nothing_the_layer_has_is_an_error_not_an_empty_matrix():
    obj = _obj()
    absent = [f for f in obj.assays["RNA"].features()
              if f not in set(obj.assays["RNA"]._scaled_features)]
    with pytest.raises(ValueError, match="None of the"):
        run_pca(obj, n_pcs=3, features=absent)


# ---------------------------------------------------------------------------
# 2. The labels match the rows — the defect that produced wrong answers
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("use_v5", [True, False])
def test_pca_loadings_are_labelled_with_the_features_actually_used(use_v5):
    """Row *i* of the loadings must be the gene `features()[i]` says it is.

    Checked two ways, because they fail to different mutations: the name list
    must equal the surviving features (a stored `list(features)` fails this),
    and the matrix must equal a run that asked only for those (a filter that
    keeps the right count in the wrong order fails this one).
    """
    obj = _obj(use_v5=use_v5)
    asked, present = _interleaved(obj)

    with pytest.warns(RuntimeWarning):
        run_pca(obj, n_pcs=3, features=asked)
    dr = obj.reductions["pca"]

    assert dr.features() == present
    assert dr.feature_loadings.shape[0] == len(dr.features()), (
        "names and rows must be the same length — DimReduc does not check this, "
        "which is how the mismatch survived"
    )

    clean = _obj(use_v5=use_v5)
    run_pca(clean, n_pcs=3, features=present)
    np.testing.assert_allclose(
        dr.feature_loadings,
        clean.reductions["pca"].feature_loadings,
        rtol=1e-12,
        atol=1e-12,
        err_msg="asking for extra absent features changed which rows were used",
    )


def test_ica_labels_its_loadings_too():
    """`run_ica` stored `list(features)` from the same copied line."""
    obj = _obj()
    asked, present = _interleaved(obj)
    with pytest.warns(RuntimeWarning):
        run_ica(obj, nics=3, features=asked)
    dr = obj.reductions["ica"]
    assert dr.features() == present
    assert dr.feature_loadings.shape[0] == len(dr.features())


def test_jackstraw_gets_back_the_features_the_matrix_holds():
    """The p-value for gene *i* is only about gene *i* if these agree.

    `_scaled_matrix_for_reduction` returned `dr.features()` unchanged next to a
    matrix `_get_scaled_data` had already filtered, so the two disagreed exactly
    when the reduction's feature list held something the layer did not.
    """
    obj = _obj()
    asked, present = _interleaved(obj)
    with pytest.warns(RuntimeWarning):
        run_pca(obj, n_pcs=3, features=asked)
    dr = obj.reductions["pca"]
    # Put the absent features back on the reduction, as a hand-built or
    # round-tripped DimReduc may carry them.
    dr._feature_names = list(asked)

    with pytest.warns(RuntimeWarning):
        mat, features = _scaled_matrix_for_reduction(obj, dr, "scale.data")
    assert features == present
    assert mat.shape[0] == len(features)


# ---------------------------------------------------------------------------
# 3. The v3 Assay's scale_data now knows its own rows
# ---------------------------------------------------------------------------

def test_v3_scale_data_reports_the_features_it_holds():
    """Not the assay's first n features, which is what stood here.

    ScaleData's subset is the *variable* features, scattered through the assay,
    so `_feature_names[:n]` was right only when the variable features happened to
    be the leading ones — never, on real data.
    """
    obj = _obj(use_v5=False)
    assay = obj.assays["RNA"]
    labels = assay.features("scale_data")

    assert len(labels) == assay.scale_data.shape[0]
    assert set(labels) == set(assay.var_features)
    assert labels != assay.features()[: len(labels)], (
        "fixture is degenerate: the variable features are the leading ones, so "
        "the old guess would pass"
    )


def test_v3_and_v5_read_identical_rows_out_of_scale_data():
    """The v3 path indexed a 20-row matrix by positions in a 60-gene list.

    Depending on where the variable features landed that either raised
    `IndexError` or returned another gene's row. Checked against the v5 path
    rather than a recorded matrix: the two architectures scale the same numbers,
    so agreement is the property, and v5 was the one that already worked.
    """
    from shanuz.reduction import _get_scaled_data

    v5, v3 = _obj(use_v5=True), _obj(use_v5=False)
    feats = list(v5.assays["RNA"]._scaled_features)

    m5 = _get_scaled_data(v5.assays["RNA"], feats, "scale.data")
    m3 = _get_scaled_data(v3.assays["RNA"], feats, "scale.data")
    np.testing.assert_allclose(m3, m5, rtol=1e-12, atol=1e-12)


def test_v3_subset_keeps_the_scaled_layer():
    """`subset` dropped scale_data unless it was full height — i.e. always.

    It subset every other layer by the assay's row indices, which address a
    different matrix than scale_data's, so throwing it away was the only safe
    thing to do without labels. With labels it can be subset by name.
    """
    obj = _obj(use_v5=False)
    assay = obj.assays["RNA"]
    scaled = list(assay.features("scale_data"))
    keep = scaled[:5] + [f for f in assay.features() if f not in set(scaled)][:5]

    sub = assay.subset(features=keep)
    assert sub.features("scale_data") == scaled[:5]
    assert sub.scale_data.shape[0] == 5
    rows = [scaled.index(f) for f in scaled[:5]]
    np.testing.assert_allclose(sub.scale_data, assay.scale_data[rows, :])


def test_v3_assay_refuses_to_invent_labels_for_a_subset_scale_data():
    """A subset scale_data with no names is unrecoverable — say so, don't guess."""
    names = [f"g{i}" for i in range(10)]
    cells = [f"c{j}" for j in range(4)]
    with pytest.raises(ValueError, match="holds a subset"):
        Assay(
            counts=sp.csc_matrix(np.ones((10, 4))),
            scale_data=np.zeros((3, 4)),
            feature_names=names,
            cell_names=cells,
        )

    with pytest.raises(ValueError, match="3 names for 2 rows"):
        Assay(
            counts=sp.csc_matrix(np.ones((10, 4))),
            scale_data=np.zeros((2, 4)),
            scaled_features=names[:3],
            feature_names=names,
            cell_names=cells,
        )


def test_v3_relabels_scale_data_when_the_layer_is_replaced():
    """Stale labels are worse than none: every reader trusts them."""
    obj = _obj(use_v5=False)
    assay = obj.assays["RNA"]
    assert len(assay.features("scale_data")) == N_VAR

    new_feats = assay.features()[:5]
    assay.set_assay_data("scale_data", np.zeros((5, N_CELLS)), new_feats)
    assert assay.features("scale_data") == new_feats


# ---------------------------------------------------------------------------
# 4. Anchors: a per-object drop would misalign the matrices, so it cannot happen
# ---------------------------------------------------------------------------

def test_integration_features_intersects_against_the_layer_not_the_assay():
    """The reference and query matrices are multiplied row-for-row.

    `_integration_features` filtered against `features()` — the assay's full
    list — while `_anchor_feature_matrix` filtered against the scaled layer, so a
    feature one object had never scaled survived the first filter and was dropped
    by the second, in that object alone. The two matrices then had different
    heights (a crash) or, with a compensating drop elsewhere, the same height and
    different genes.
    """
    from shanuz.anchors import _anchor_feature_matrix, _integration_features

    a, b = _obj(seed=1, tag="a"), _obj(seed=2, tag="b")
    # Force the objects to disagree: rescale `b` on a different feature set.
    b_feats = list(b.assays["RNA"]._scaled_features)[:12]
    scale_data(b, features=b_feats)

    shared = _integration_features([a, b], None, "scale.data")
    assert shared, "the two objects do share scaled features"
    a_scaled = set(a.assays["RNA"]._scaled_features)
    b_scaled = set(b.assays["RNA"]._scaled_features)
    assert set(shared) <= a_scaled & b_scaled

    # Nothing is dropped downstream, so the two matrices line up row-for-row.
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        ma = _anchor_feature_matrix(a, shared, "scale.data")
        mb = _anchor_feature_matrix(b, shared, "scale.data")
    assert ma.shape[0] == mb.shape[0] == len(shared)


def test_anchor_matrix_refuses_a_partial_feature_set():
    """Where a reduction warns, anchors must raise — the rows are paired."""
    from shanuz.anchors import _anchor_feature_matrix

    a = _obj(seed=1, tag="a")
    scaled = list(a.assays["RNA"]._scaled_features)
    absent = [f for f in a.assays["RNA"].features() if f not in set(scaled)][:3]
    with pytest.raises(ValueError, match="row-for-row"):
        _anchor_feature_matrix(a, scaled[:5] + absent, "scale.data")
