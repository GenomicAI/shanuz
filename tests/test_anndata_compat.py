import numpy as np
import pytest

anndata = pytest.importorskip("anndata", reason="anndata not installed")

from shanuz.compat import as_anndata, from_anndata


def test_as_anndata(small_seurat):
    adata = as_anndata(small_seurat)
    assert adata.n_obs == 20
    assert adata.n_vars == 50


def test_as_anndata_obs(small_seurat):
    adata = as_anndata(small_seurat)
    assert "ident" in adata.obs.columns


def test_as_anndata_uns(small_seurat):
    adata = as_anndata(small_seurat)
    assert "project_name" in adata.uns


def test_from_anndata(small_seurat):
    adata = as_anndata(small_seurat)
    seurat2 = from_anndata(adata, assay="RNA")
    assert len(seurat2) == 20
    assert len(seurat2.feature_names()) == 50


def test_roundtrip_cells(small_seurat):
    adata = as_anndata(small_seurat)
    seurat2 = from_anndata(adata, assay="RNA")
    assert set(seurat2.cell_names()) == set(small_seurat.cell_names())


def test_roundtrip_features(small_seurat):
    adata = as_anndata(small_seurat)
    seurat2 = from_anndata(adata, assay="RNA")
    assert set(seurat2.feature_names()) == set(small_seurat.feature_names())
