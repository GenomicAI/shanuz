"""Tests for the additional reductions & UMAP graph path (v0.5.0 / v0.4.0).

  * run_tsne              (reduction.py)
  * run_ica               (reduction.py)
  * run_umap(graph=...)   (umap.py)  — embed a precomputed graph directly

Network-free; small synthetic object with three well-separated clusters.
"""
import sys
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz.shanuz import create_shanuz_object  # noqa: E402
from shanuz.preprocessing import (  # noqa: E402
    normalize_data,
    find_variable_features,
    scale_data,
)
from shanuz.reduction import run_pca, run_ica, run_tsne  # noqa: E402
from shanuz.neighbors import find_neighbors  # noqa: E402
from shanuz.umap import run_umap  # noqa: E402

pytest.importorskip("sklearn")
pytest.importorskip("umap")


def _clustered_object(seed=0, per=30):
    rng = np.random.default_rng(seed)
    G, K = 150, 3
    n = K * per
    mat = rng.gamma(0.3, size=(G, n)) + 0.05
    cells = []
    for ci in range(n):
        k = ci // per
        mat[k * 40:(k + 1) * 40, ci] += 6.0   # cluster-specific gene block
        cells.append(f"c{ci}")
    obj = create_shanuz_object(
        counts=sp.csc_matrix(rng.poisson(mat * 3000.0 / mat.sum(axis=0, keepdims=True)).astype(float)),
        assay="RNA", feature_names=[f"g{i}" for i in range(G)], cell_names=cells,
    )
    normalize_data(obj)
    find_variable_features(obj, selection_method="vst", nfeatures=100)
    scale_data(obj)
    run_pca(obj, n_pcs=15)
    return obj, n


def test_run_tsne():
    obj, n = _clustered_object()
    run_tsne(obj, dims=range(10), perplexity=15.0, seed=0)
    assert "tsne" in obj.reductions
    emb = obj.reductions["tsne"].cell_embeddings
    assert emb.shape == (n, 2)
    assert np.isfinite(emb).all()


def test_run_ica():
    obj, n = _clustered_object()
    run_ica(obj, nics=10, seed=0)
    assert "ica" in obj.reductions
    dr = obj.reductions["ica"]
    assert dr.cell_embeddings.shape == (n, 10)
    assert dr.feature_loadings.shape[1] == 10          # loadings stored
    assert np.isfinite(dr.cell_embeddings).all()
    # downstream neighbours accept reduction="ica"
    find_neighbors(obj, reduction="ica", graph_name="ica")
    assert "ica_snn" in obj.graphs


def test_run_umap_from_graph():
    obj, n = _clustered_object()
    find_neighbors(obj, dims=range(10), reduction="pca")   # builds RNA_snn
    run_umap(obj, graph="RNA_snn", seed=0)
    assert "umap" in obj.reductions
    dr = obj.reductions["umap"]
    assert dr.cell_embeddings.shape == (n, 2)
    assert np.isfinite(dr.cell_embeddings).all()
    assert dr.misc.get("umap_graph") == "RNA_snn"


def test_run_umap_missing_graph_raises():
    obj, _ = _clustered_object()
    with pytest.raises(KeyError):
        run_umap(obj, graph="does_not_exist")
