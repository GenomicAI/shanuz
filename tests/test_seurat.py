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


# ----------------------------------------------------------------------
# FetchData / metadata — R fidelity
#
# Found by the object-internals tutorial (T-obj). `test_fetch_data` above
# passed throughout: it asserts the column name and the row count, which a
# column of identical sparse-matrix objects satisfies perfectly.
# ----------------------------------------------------------------------


def _prepared(small_seurat):
    """The fixture with a `data` layer and a two-component reduction."""
    from shanuz import normalize_data
    normalize_data(small_seurat)
    emb = np.arange(40, dtype=float).reshape(20, 2)
    small_seurat.reductions["pca"] = DimReduc(
        cell_embeddings=emb,
        key="PC_",
        assay_used="RNA",
        cell_names=small_seurat.cell_names(),
    )
    return small_seurat


def test_fetch_data_returns_numbers_for_a_feature(small_seurat, feature_names):
    obj = _prepared(small_seurat)
    column = obj.fetch_data([feature_names[0]])[feature_names[0]]

    # `np.asarray(sparse)` yields a 0-d object array wrapping the matrix, so
    # the old `.flatten()` broadcast one csc_matrix down all 20 rows.
    assert column.dtype != object
    assert len(column) == 20
    expected = obj.get_assay().layer_data(
        layer="data", features=[feature_names[0]], cells=obj.cell_names())
    assert np.allclose(column.to_numpy(), expected.toarray().ravel())


def test_fetch_data_reads_the_data_layer_not_counts(small_seurat, feature_names):
    # R's FetchData defaults to `data`; deferring to the layered assay's own
    # default picked `counts`, so a fetched gene came back as raw integers.
    obj = _prepared(small_seurat)
    column = obj.fetch_data([feature_names[0]])[feature_names[0]].to_numpy()
    counts = obj.get_assay().layer_data(
        layer="counts", features=[feature_names[0]], cells=obj.cell_names())
    assert not np.allclose(column, counts.toarray().ravel())


def test_fetch_data_falls_back_to_counts_with_a_warning(small_seurat, feature_names):
    # Seurat: "data layer is not found and counts layer is used".
    with pytest.warns(UserWarning, match="counts layer is used"):
        small_seurat.fetch_data([feature_names[0]])


def test_fetch_data_addresses_an_embedding_by_its_key(small_seurat):
    # `FetchData(obj, "PC_1")` is how every vignette asks for a component;
    # shanuz accepted only the reduction name and raised KeyError on this.
    obj = _prepared(small_seurat)
    got = obj.fetch_data(["PC_1"])["PC_1"].to_numpy()
    assert np.array_equal(got, obj.reductions["pca"].cell_embeddings[:, 0])


def test_fetch_data_names_reduction_columns_by_key(small_seurat):
    # Not `pca_1`: R names embedding columns from the reduction's Key(), which
    # is the same name the branch above resolves.
    obj = _prepared(small_seurat)
    assert list(obj.fetch_data(["pca"]).columns) == ["PC_1", "PC_2"]


def test_object_has_orig_ident(small_seurat):
    # The first column of every Seurat object's metadata, and the default
    # identity class; scripts group and split on it.
    assert small_seurat.meta_data.columns[0] == "orig.ident"
    assert set(small_seurat.meta_data["orig.ident"]) == {"TestProject"}


def test_add_meta_data_accepts_a_plain_vector(small_seurat):
    # R's AddMetaData documents "a vector, list, or data.frame"; shanuz raised
    # TypeError on the vector, which is the form the vignettes use.
    small_seurat.add_meta_data(np.arange(20), "depth")
    assert small_seurat.meta_data["depth"].tolist() == list(range(20))


def test_add_meta_data_rejects_a_wrong_length_vector(small_seurat):
    with pytest.raises(ValueError, match="entries but the object"):
        small_seurat.add_meta_data(np.arange(5), "depth")
