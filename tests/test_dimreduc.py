import numpy as np
import pytest
from shanuz import DimReduc


@pytest.fixture
def pca(cell_names):
    rng = np.random.default_rng(0)
    emb = rng.standard_normal((20, 10))
    return DimReduc(
        cell_embeddings=emb,
        cell_names=cell_names,
        assay_used="RNA",
        key="PC_",
    )


def test_cells(pca, cell_names):
    assert pca.cells() == cell_names


def test_embeddings_shape(pca):
    assert pca.embeddings().shape == (20, 10)


def test_loadings_empty(pca):
    assert pca.loadings().shape[0] == 0


def test_set_loadings(pca, feature_names):
    load = np.random.randn(50, 10)
    pca.set_loadings(load)
    assert pca.loadings().shape == (50, 10)


def test_is_global_default(pca):
    assert not pca.is_global()


def test_rename_cells(pca):
    new_names = [f"new_{i}" for i in range(20)]
    renamed = pca.rename_cells(new_names)
    assert renamed.cells() == new_names
    assert pca.cells()[0] != new_names[0]


def test_subset(pca, cell_names):
    sub = pca.subset(cells=cell_names[:5])
    assert sub.embeddings().shape == (5, 10)


def test_repr(pca):
    assert "DimReduc" in repr(pca)
