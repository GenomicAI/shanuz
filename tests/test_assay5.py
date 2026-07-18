import numpy as np
import pytest
import scipy.sparse as sp
from shanuz import Assay5, create_assay5_object


def test_create(small_assay5):
    assert len(small_assay5.cells()) == 20
    assert len(small_assay5.features()) == 50


def test_layers_list(small_assay5):
    names = small_assay5.layers_list()
    assert "counts" in names


def test_default_layer(small_assay5):
    assert small_assay5.default_layer == "counts"


def test_set_default_layer(small_assay5):
    new_data = sp.csc_matrix(np.ones((50, 20)))
    small_assay5.layers["data"] = new_data
    small_assay5.default_layer = "data"
    assert small_assay5.default_layer == "data"


def test_layer_data(small_assay5):
    mat = small_assay5.layer_data("counts")
    assert mat.shape == (50, 20)


def test_layer_data_subset(small_assay5):
    cells = [f"cell_{i}" for i in range(5)]
    feats = [f"gene_{i}" for i in range(10)]
    sub = small_assay5.layer_data("counts", cells=cells, features=feats)
    assert sub.shape == (10, 5)


def test_variable_features(small_assay5):
    feats = [f"gene_{i}" for i in range(5)]
    small_assay5.variable_features = feats
    assert small_assay5.variable_features == feats


def test_subset(small_assay5):
    cells = [f"cell_{i}" for i in range(5)]
    feats = [f"gene_{i}" for i in range(10)]
    sub = small_assay5.subset(cells=cells, features=feats)
    assert len(sub.cells()) == 5
    assert len(sub.features()) == 10


def test_cast_assay_to_dense(small_assay5):
    dense = small_assay5.cast_assay(to_sparse=False)
    for mat in dense.layers.values():
        assert isinstance(mat, np.ndarray)


def test_repr(small_assay5):
    assert "Assay5" in repr(small_assay5)


# ----------------------------------------------------------------------
# split / join layers — R fidelity
#
# Found by the object-internals tutorial (T-obj). Before it, `split_layers`
# and `join_layers` had no call sites and no tests anywhere in the package,
# and the round trip was not the identity in three separate ways.
# ----------------------------------------------------------------------


def _ordered_assay():
    """4 features x 6 cells whose values encode their position, so a permuted
    column is visible by inspection rather than only by a checksum."""
    from shanuz import create_assay5_object
    mat = sp.csc_matrix(np.arange(1, 25, dtype=float).reshape(4, 6))
    return create_assay5_object(
        counts=mat,
        feature_names=[f"g{i}" for i in range(4)],
        cell_names=[f"c{i}" for i in range(6)],
        key="rna_",
    )


def test_split_layers_uses_seurats_dot_separator():
    # Seurat names the parts `counts.batch1`; users match on that spelling.
    split = _ordered_assay().split_layers(["a", "b"] * 3, layer="counts")
    assert sorted(split.layers_list()) == ["counts.a", "counts.b"]


def test_join_layers_restores_the_original_layer_name():
    split = _ordered_assay().split_layers(["a", "b"] * 3, layer="counts")
    assert split.join_layers().layers_list() == ["counts"]


def test_join_layers_restores_the_original_cell_order():
    # The interleaved split is the point: joining in layer order would give
    # c0 c2 c4 c1 c3 c5 while the assay's own cell vector never moved, so the
    # matrix would be silently transposed against the metadata indexing it.
    #
    # Binds to whatever layer the join produced rather than naming "counts":
    # `cells()` returns the assay's full cell vector for a layer that does not
    # exist, so asking for "counts" after a join that produced "joined"
    # compares the object against itself and passes on the broken code.
    assay = _ordered_assay()
    split = assay.split_layers(["a", "b"] * 3, layer="counts")
    joined = split.join_layers()
    (name,) = joined.layers_list()
    assert joined.cells(name) == assay.cells()


def test_split_join_round_trip_is_the_identity():
    # Same reasoning: read back the layer that exists, so this fails on the
    # matrix's contents rather than on a KeyError for a differently-named one.
    assay = _ordered_assay()
    before = assay.layer_data("counts").toarray()
    joined = assay.split_layers(["a", "b"] * 3, layer="counts").join_layers()
    (name,) = joined.layers_list()
    assert np.array_equal(joined.layer_data(name).toarray(), before)


def test_join_layers_leaves_unsplit_layers_alone():
    # The no-argument call is R's idiom, and a prepared assay always carries
    # layers of differing feature counts. Blindly hstacking them raised
    # ValueError; `scale.data` must come through untouched — and must not be
    # mistaken for layer `scale` split into group `data`.
    assay = _ordered_assay()
    assay._add_layer("scale.data", np.ones((2, 6)),
                     feature_names=["g0", "g1"],
                     cell_names=assay.cells())
    joined = assay.split_layers(["a", "b"] * 3, layer="counts").join_layers()
    assert set(joined.layers_list()) == {"counts", "scale.data"}
    assert joined.layer_data("scale.data").shape == (2, 6)


def test_split_layers_generic_is_registered():
    # `split_layers` was declared in shanuz.generics but never registered for
    # any type, so the documented generic raised NotImplementedError while the
    # method it should dispatch to worked fine.
    from shanuz import generics as G
    split = G.split_layers(_ordered_assay(), ["a", "b"] * 3)
    assert sorted(split.layers_list()) == ["counts.a", "counts.b"]


def test_subset_keeps_the_split_provenance():
    # Without this, a subset of a split assay can never be rejoined:
    # `join_layers` finds no stems, returns the object unchanged, and reports
    # no error — the layers stay split forever.
    assay = _ordered_assay()
    split = assay.split_layers(["a", "b"] * 3, layer="counts")
    sub = split.subset(cells=["c0", "c2", "c4"])
    joined = sub.join_layers()
    assert joined.layers_list() == ["counts"]
    assert joined.cells("counts") == ["c0", "c2", "c4"]
