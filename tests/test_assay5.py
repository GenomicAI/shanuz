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
