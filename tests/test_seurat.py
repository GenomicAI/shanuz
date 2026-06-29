import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from shanuz import DimReduc, Shanuz, create_shanuz_object


def test_create(small_seurat):
    assert len(small_seurat) == 20


def test_cell_names(small_seurat, cell_names):
    assert small_seurat.cell_names() == cell_names


def test_feature_names(small_seurat, feature_names):
    assert small_seurat.feature_names() == feature_names


def test_active_assay(small_seurat):
    assert small_seurat.active_assay == "RNA"


def test_meta_data_has_ncount(small_seurat):
    assert "nCount_RNA" in small_seurat.meta_data.columns


def test_add_meta_data(small_seurat):
    series = pd.Series(range(20), index=small_seurat.cell_names())
    small_seurat.add_meta_data(series, col_name="batch")
    assert "batch" in small_seurat.meta_data.columns


def test_idents_set_and_get(small_seurat, cell_names):
    small_seurat.idents = ["A"] * 10 + ["B"] * 10
    assert "A" in list(small_seurat.idents)
    assert "B" in list(small_seurat.idents)


def test_stash_ident(small_seurat):
    small_seurat.stash_ident("old_ident")
    assert "old_ident" in small_seurat.meta_data.columns


def test_rename_idents(small_seurat):
    small_seurat.idents = ["A"] * 20
    small_seurat.rename_idents({"A": "B"})
    assert all(x == "B" for x in small_seurat.idents)


def test_which_cells(small_seurat, cell_names):
    small_seurat.idents = ["A"] * 10 + ["B"] * 10
    a_cells = small_seurat.which_cells(ident="A")
    assert len(a_cells) == 10
    assert all(c in cell_names[:10] for c in a_cells)


def test_fetch_data(small_seurat, feature_names, cell_names):
    df = small_seurat.fetch_data([feature_names[0]], cells=cell_names[:5])
    assert feature_names[0] in df.columns
    assert len(df) == 5


def test_fetch_meta(small_seurat, cell_names):
    df = small_seurat.fetch_data(["nCount_RNA"], cells=cell_names[:3])
    assert "nCount_RNA" in df.columns


def test_subset_cells(small_seurat, cell_names):
    sub = small_seurat.subset(cells=cell_names[:5])
    assert len(sub) == 5


def test_subset_idents(small_seurat):
    small_seurat.idents = ["A"] * 10 + ["B"] * 10
    sub = small_seurat.subset(idents="A")
    assert len(sub) == 10


def test_rename_cells(small_seurat, cell_names):
    new_names = [f"new_{i}" for i in range(20)]
    renamed = small_seurat.rename_cells(new_names)
    assert renamed.cell_names() == new_names
    assert small_seurat.cell_names() == cell_names  # immutable


def test_merge(small_seurat, feature_names):
    other = create_shanuz_object(
        counts=sp.csc_matrix(np.ones((50, 10))),
        assay="RNA",
        feature_names=feature_names,
        cell_names=[f"other_{i}" for i in range(10)],
        project="OtherProject",
    )
    merged = small_seurat.merge(other, add_cell_ids=["orig", "other"])
    assert len(merged) == 30


def test_add_reduction(small_seurat, cell_names):
    rng = np.random.default_rng(3)
    emb = rng.standard_normal((20, 5))
    pca = DimReduc(cell_embeddings=emb, cell_names=cell_names, key="PC_", assay_used="RNA")
    small_seurat.reductions["pca"] = pca
    assert "pca" in small_seurat.reduction_names()


def test_embeddings(small_seurat, cell_names):
    rng = np.random.default_rng(3)
    emb = rng.standard_normal((20, 5))
    pca = DimReduc(cell_embeddings=emb, cell_names=cell_names, key="PC_", assay_used="RNA")
    small_seurat.reductions["pca"] = pca
    result = small_seurat.embeddings("pca")
    assert result.shape == (20, 5)


def test_default_assay_setter(small_seurat):
    with pytest.raises(KeyError):
        small_seurat.default_assay = "ATAC"


def test_repr(small_seurat):
    r = repr(small_seurat)
    assert "Shanuz" in r
    assert "TestProject" in r


def test_getattr_metadata(small_seurat):
    small_seurat.add_meta_data(
        pd.Series([1] * 20, index=small_seurat.cell_names()), col_name="batch"
    )
    col = small_seurat.batch
    assert len(col) == 20
