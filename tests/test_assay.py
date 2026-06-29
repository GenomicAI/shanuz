import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from shanuz import Assay, create_assay_object


def test_create_from_counts(small_assay):
    assert len(small_assay.cells()) == 20
    assert len(small_assay.features()) == 50


def test_layer_data(small_assay):
    mat = small_assay.layer_data("counts")
    assert mat.shape == (50, 20)


def test_layer_data_subset(small_assay, feature_names, cell_names):
    sub = small_assay.layer_data("counts", cells=cell_names[:5], features=feature_names[:10])
    assert sub.shape == (10, 5)


def test_set_assay_data(small_assay):
    new_data = sp.csc_matrix(np.ones((50, 20)))
    small_assay.set_assay_data("data", new_data)
    mat = small_assay.layer_data("data")
    assert mat.shape == (50, 20)


def test_variable_features(small_assay, feature_names):
    small_assay.variable_features = feature_names[:5]
    assert small_assay.variable_features == feature_names[:5]


def test_calc_n(small_assay):
    df = small_assay.calc_n()
    assert "nCount" in df.columns
    assert "nFeature" in df.columns
    assert len(df) == 20


def test_rename_cells(small_assay, cell_names):
    new_names = [f"renamed_{i}" for i in range(20)]
    renamed = small_assay.rename_cells(new_names)
    assert renamed.cells() == new_names


def test_subset(small_assay, feature_names, cell_names):
    sub = small_assay.subset(cells=cell_names[:5], features=feature_names[:10])
    assert len(sub.cells()) == 5
    assert len(sub.features()) == 10


def test_merge(small_assay, feature_names):
    other = create_assay_object(
        counts=sp.csc_matrix(np.ones((50, 10))),
        feature_names=feature_names,
        cell_names=[f"other_{i}" for i in range(10)],
    )
    merged = small_assay.merge(other)
    assert len(merged.cells()) == 30


def test_key_validation():
    with pytest.raises(ValueError):
        create_assay_object(
            counts=sp.csc_matrix(np.ones((5, 5))),
            key="badkey",  # missing trailing _
        )


def test_repr(small_assay):
    r = repr(small_assay)
    assert "Assay" in r
