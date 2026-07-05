"""Tests for Weighted Nearest Neighbor multimodal integration (v0.4.0).

  * find_multi_modal_neighbors  (multimodal.py)

Builds a synthetic two-modality object with COMPLEMENTARY structure: RNA
separates cluster A from {B,C}; ADT separates C from {A,B}. Neither modality
alone recovers all three clusters, but WNN should. Network-free.
"""
import sys
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz.shanuz import create_shanuz_object  # noqa: E402
from shanuz.assay5 import create_assay5_object  # noqa: E402
from shanuz.preprocessing import (  # noqa: E402
    normalize_data,
    find_variable_features,
    scale_data,
)
from shanuz.reduction import run_pca  # noqa: E402
from shanuz.neighbors import find_neighbors  # noqa: E402
from shanuz.clustering import find_clusters  # noqa: E402
from shanuz.multimodal import find_multi_modal_neighbors  # noqa: E402

pytest.importorskip("sklearn")
from sklearn.metrics import adjusted_rand_score  # noqa: E402


def _complementary_object(seed=0, per=40):
    rng = np.random.default_rng(seed)
    groups = ("A", "B", "C")
    n = len(groups) * per

    # RNA: genes 0-29 high only in A  → RNA merges B and C.
    Grna = 120
    rna = rng.gamma(0.3, size=(Grna, n)) + 0.05
    # ADT: proteins 0-9 high only in C → ADT merges A and B.
    Padt = 30
    adt = rng.gamma(0.3, size=(Padt, n)) + 0.05

    labels, cells = [], []
    for ci in range(n):
        g = groups[ci // per]
        if g == "A":
            rna[0:30, ci] += 6.0
        if g == "C":
            adt[0:10, ci] += 6.0
        labels.append(g)
        cells.append(f"cell{ci}")

    rna_counts = rng.poisson(rna * 3000.0 / rna.sum(axis=0, keepdims=True))
    adt_counts = rng.poisson(adt * 1000.0 / adt.sum(axis=0, keepdims=True))

    obj = create_shanuz_object(
        counts=sp.csc_matrix(rna_counts.astype(float)), assay="RNA",
        feature_names=[f"gene{i}" for i in range(Grna)], cell_names=cells,
    )
    normalize_data(obj)
    find_variable_features(obj, selection_method="vst", nfeatures=100)
    scale_data(obj)
    run_pca(obj, n_pcs=15)

    obj.assays["ADT"] = create_assay5_object(
        counts=sp.csc_matrix(adt_counts.astype(float)),
        feature_names=[f"prot{i}" for i in range(Padt)],
        cell_names=cells, key="adt_",
    )
    normalize_data(obj, assay="ADT", normalization_method="CLR", margin=2)
    scale_data(obj, assay="ADT")
    run_pca(obj, assay="ADT", reduction_name="apca", n_pcs=10)

    return obj, np.array(labels)


def test_wnn_builds_graphs_and_weights():
    obj, labels = _complementary_object()
    find_multi_modal_neighbors(
        obj, reduction_list=["pca", "apca"],
        dims_list=[range(15), range(10)], k_nn=20,
    )

    n = len(labels)
    assert "wknn" in obj.graphs and "wsnn" in obj.graphs
    assert obj.graphs["wsnn"].shape == (n, n)

    # Per-cell modality weights present, in [0,1], and sum to 1 across modalities.
    assert "RNA.weight" in obj.meta_data.columns
    assert "ADT.weight" in obj.meta_data.columns
    w = obj.meta_data[["RNA.weight", "ADT.weight"]].to_numpy()
    assert np.all(w >= 0) and np.all(w <= 1)
    assert np.allclose(w.sum(axis=1), 1.0, atol=1e-6)


def test_wnn_recovers_structure_better_than_rna_alone():
    obj, labels = _complementary_object()

    # RNA-only clustering (merges B and C).
    find_neighbors(obj, dims=range(15), reduction="pca")
    find_clusters(obj, resolution=1.0, graph_name="RNA_snn")
    rna_clusters = obj.meta_data["seurat_clusters"].to_numpy()

    # WNN clustering on the joint wsnn graph.
    find_multi_modal_neighbors(
        obj, reduction_list=["pca", "apca"],
        dims_list=[range(15), range(10)], k_nn=20,
    )
    find_clusters(obj, resolution=1.0, graph_name="wsnn")
    wnn_clusters = obj.meta_data["seurat_clusters"].to_numpy()

    ari_rna = adjusted_rand_score(labels, rna_clusters)
    ari_wnn = adjusted_rand_score(labels, wnn_clusters)

    # WNN should recover the 3-cluster ground truth better than RNA alone.
    assert ari_wnn > ari_rna
    assert len(np.unique(wnn_clusters)) >= len(np.unique(rna_clusters))


def test_wnn_requires_two_reductions():
    obj, _ = _complementary_object()
    with pytest.raises(ValueError):
        find_multi_modal_neighbors(obj, reduction_list=["pca"])
