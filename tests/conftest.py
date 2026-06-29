"""Shared fixtures for shanuz tests."""
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


@pytest.fixture
def small_counts():
    """50 features × 20 cells, Poisson-distributed counts."""
    rng = np.random.default_rng(42)
    data = rng.poisson(1.0, size=(50, 20)).astype(float)
    return sp.csc_matrix(data)


@pytest.fixture
def feature_names():
    return [f"gene_{i}" for i in range(50)]


@pytest.fixture
def cell_names():
    return [f"cell_{i}" for i in range(20)]


@pytest.fixture
def small_assay(small_counts, feature_names, cell_names):
    from shanuz import create_assay_object
    return create_assay_object(
        counts=small_counts,
        feature_names=feature_names,
        cell_names=cell_names,
        key="rna_",
    )


@pytest.fixture
def small_assay5(small_counts, feature_names, cell_names):
    from shanuz import create_assay5_object
    return create_assay5_object(
        counts=small_counts,
        feature_names=feature_names,
        cell_names=cell_names,
        key="rna_",
    )


@pytest.fixture
def small_seurat(small_counts, feature_names, cell_names):
    from shanuz import create_shanuz_object
    return create_shanuz_object(
        counts=small_counts,
        assay="RNA",
        feature_names=feature_names,
        cell_names=cell_names,
        project="TestProject",
    )
