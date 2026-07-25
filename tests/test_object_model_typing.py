"""Behaviour behind the object model's type annotations.

Two of the fixes that cleared mypy from `shanuz.py` / `assay5.py` /
`reduction.py` / `spatial/fov.py` changed what the code does, not just what it
claims — a merge across assay classes, and `run_pca` on an unscaled Assay5. Both
are pinned here, because a type checker is advisory in this repo (`|| true` in
CI) and a count that drifts across matrix legs cannot be the only thing standing
between the object model and a regression.

The annotations themselves — `StdAssay`'s methods returning `Self` rather than
the abstract base — are mypy's business and are not asserted here. What *is*
asserted is the runtime property `Self` promises: that those methods build their
result from `self.__class__`, so a subclass gets its own class back. Hardcoding
`Assay5(...)` anywhere in them would keep every other test green and quietly make
the annotation a lie.
"""
import numpy as np
import pytest
import scipy.sparse as sp

from shanuz import (
    create_shanuz_object,
    find_variable_features,
    normalize_data,
    run_pca,
)
from shanuz.assay import Assay
from shanuz.assay5 import Assay5


def _obj(use_v5=True, tag="c", n_genes=80, n_cells=40, seed=0):
    rng = np.random.default_rng(seed)
    counts = sp.csc_matrix(rng.poisson(3.0, size=(n_genes, n_cells)).astype(float))
    return create_shanuz_object(
        counts=counts,
        assay="RNA",
        feature_names=[f"g{i}" for i in range(n_genes)],
        cell_names=[f"{tag}{j}" for j in range(n_cells)],
        use_v5=use_v5,
    )


# ---------------------------------------------------------------------------
# 1. Merging across assay classes
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("first_is_v5", [True, False])
def test_merging_across_assay_classes_names_the_mismatch(first_is_v5):
    """A v3 Assay and a v5 Assay5 cannot be merged, and now say so.

    They keep their cells in differently-named private slots, so the merge used
    to reach for one the other does not have and die partway through with
    `AttributeError: 'Assay' object has no attribute '_all_cell_names'` — a
    message about a private attribute, from inside a method the caller never
    named, after the cell names and metadata had already been concatenated.
    """
    a = _obj(use_v5=first_is_v5, tag="a")
    b = _obj(use_v5=not first_is_v5, tag="b")

    with pytest.raises(TypeError) as excinfo:
        a.merge(b)

    message = str(excinfo.value)
    assert "RNA" in message, "the message must name the assay that clashed"
    assert "Assay5" in message and "Assay" in message, (
        "the message must name both classes — knowing which two collided is "
        "the difference between a fix and a guess"
    )


def test_merging_matching_assay_classes_still_works():
    """The guard must not cost the ordinary case."""
    a, b = _obj(tag="a"), _obj(tag="b")
    merged = a.merge(b)
    assert len(merged.cell_names()) == 80
    assert type(merged.assays["RNA"]) is Assay5

    a3, b3 = _obj(use_v5=False, tag="a"), _obj(use_v5=False, tag="b")
    merged3 = a3.merge(b3)
    assert len(merged3.cell_names()) == 80
    assert type(merged3.assays["RNA"]) is Assay


# ---------------------------------------------------------------------------
# 2. run_pca's fallback when scale.data is absent
# ---------------------------------------------------------------------------

def test_run_pca_falls_back_to_data_on_an_unscaled_assay5():
    """`X or Y` calls `bool(X)`, and scipy refuses to answer for a matrix.

    `_get_scaled_data` picked its fallback layer with
    `layers.get("data") or layers.get("counts")`. Any sparse matrix with more
    than one entry raises `ValueError: The truth value of an array with more
    than one element is ambiguous` from `__bool__`, so the fallback raised on
    every call that reached it — which is every `run_pca` on an Assay5 that had
    not been through `scale_data()`.

    The v3 path never had the bug (it reads `assay_obj.data` directly), so the
    two architectures disagreed about whether an unscaled object could be
    reduced at all. They are checked against each other rather than against a
    recorded number: identical embeddings is the property, and it makes the
    v3 path the reference for what the v5 fallback was always meant to do.
    """
    v5, v3 = _obj(use_v5=True), _obj(use_v5=False)
    for obj in (v5, v3):
        normalize_data(obj)
        find_variable_features(obj, nfeatures=30)
        assert "scale.data" not in getattr(obj.assays["RNA"], "layers", {})
        run_pca(obj, n_pcs=5)

    np.testing.assert_array_equal(
        v5.reductions["pca"].cell_embeddings,
        v3.reductions["pca"].cell_embeddings,
        err_msg="the v5 fallback no longer agrees with the v3 one it mirrors",
    )


# ---------------------------------------------------------------------------
# 3. The runtime property `Self` promises
# ---------------------------------------------------------------------------

class _DerivedAssay(Assay5):
    """A stand-in for any future Assay5 specialization. Adds nothing."""


def _derived_from(assay):
    return _DerivedAssay(
        layers={k: v.copy() for k, v in assay.layers.items()},
        feature_names=list(assay._all_feature_names),
        cell_names=list(assay._all_cell_names),
        key=assay._key,
    )


@pytest.mark.parametrize("call", [
    pytest.param(lambda a: a._copy(), id="_copy"),
    pytest.param(lambda a: a.rename_cells([f"x{i}" for i in range(40)]), id="rename_cells"),
    pytest.param(lambda a: a.subset(features=[f"g{i}" for i in range(10)]), id="subset"),
    pytest.param(lambda a: a.cast_assay(to_sparse=False), id="cast_assay"),
    pytest.param(lambda a: a.join_layers(), id="join_layers"),
    pytest.param(lambda a: a.split_layers(["p"] * 20 + ["q"] * 20), id="split_layers"),
])
def test_stdassay_methods_return_the_callers_own_class(call):
    """These are annotated `Self`; that is only true if they build via `self.__class__`.

    They used to be annotated with the abstract `StdAssay`, which forced every
    caller downstream to widen — that is what put `Assay | StdAssay` where
    `Shanuz` wanted `Assay | Assay5` and produced four of the errors this change
    cleared. `Self` is the accurate type *because* of the property asserted
    here; hardcode `Assay5(...)` in any one of these and the annotation becomes
    a promise the code does not keep, with nothing else in the suite noticing.
    """
    derived = _derived_from(_obj().assays["RNA"])
    assert type(call(derived)) is _DerivedAssay


def test_merge_also_returns_the_callers_own_class():
    a, b = _obj(tag="a").assays["RNA"], _obj(tag="b").assays["RNA"]
    merged = _derived_from(a).merge(_derived_from(b))
    assert type(merged) is _DerivedAssay
