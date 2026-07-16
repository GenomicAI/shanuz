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
from shanuz.neighbors import _build_knn  # noqa: E402
from shanuz.multimodal import (  # noqa: E402
    find_multi_modal_neighbors,
    _impute_dist,
    _l2_norm,
    _modality_weights,
    _multi_modal_nn,
    _snn_bandwidth,
)

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


def test_wnn_rejects_k_nn_above_cell_count():
    obj, labels = _complementary_object()
    with pytest.raises(ValueError, match="smaller than the number of cells"):
        find_multi_modal_neighbors(
            obj, reduction_list=["pca", "apca"], k_nn=len(labels) + 1,
        )


# ---------------------------------------------------------------------------
# Modality weights (R: FindModalityWeights)
# ---------------------------------------------------------------------------

def test_weights_follow_the_informative_modality():
    """A is separable only in RNA, C only in ADT — the weights must say so.

    This is the behaviour the exponential kernel + softmax buys. A linear
    distance ratio gets the direction right but compresses every cell to ~0.5,
    which is indistinguishable from "no preference".
    """
    obj, labels = _complementary_object()
    find_multi_modal_neighbors(
        obj, reduction_list=["pca", "apca"],
        dims_list=[range(15), range(10)], k_nn=20,
    )
    rna_w = obj.meta_data["RNA.weight"].to_numpy()
    adt_w = obj.meta_data["ADT.weight"].to_numpy()

    assert rna_w[labels == "A"].mean() > adt_w[labels == "A"].mean()
    assert adt_w[labels == "C"].mean() > rna_w[labels == "C"].mean()


def test_weights_span_a_real_range():
    """Regression guard: weights must not collapse into a band around 0.5.

    The pre-port approximation returned d_cross / (d_same + d_cross), whose
    two-modality softmax sat inside ~0.46-0.53 on real data and hid which
    modality actually drove each cell.
    """
    obj, _ = _complementary_object()
    find_multi_modal_neighbors(
        obj, reduction_list=["pca", "apca"],
        dims_list=[range(15), range(10)], k_nn=20,
    )
    adt_w = obj.meta_data["ADT.weight"].to_numpy()
    assert adt_w.max() - adt_w.min() > 0.3
    assert np.isfinite(adt_w).all()


def test_l2_norm_unit_rows_and_zero_safety():
    mat = np.array([[3.0, 4.0], [0.0, 0.0], [1.0, 0.0]])
    out = _l2_norm(mat)
    assert np.allclose(out[0], [0.6, 0.8])
    assert np.allclose(out[1], [0.0, 0.0])      # zero row must not become NaN
    assert np.allclose(np.linalg.norm(out[[0, 2]], axis=1), 1.0)


def test_impute_dist_subtracts_nearest_and_floors_at_zero():
    x = np.array([[0.0, 0.0], [0.0, 0.0]])
    y = np.array([[3.0, 4.0], [0.6, 0.8]])      # raw distances 5.0 and 1.0
    out = _impute_dist(x, y, nearest_dist=np.array([2.0, 4.0]))
    assert np.allclose(out, [3.0, 0.0])         # second is ReLU'd, not negative


def test_snn_bandwidth_uses_the_smallest_weight_edges():
    """Bandwidth reads the *least* similar SNN partners, not the nearest ones."""
    # Cell 0 shares edges with 1 (strong/near) and 2 (weak/far).
    snn = sp.csc_matrix(np.array([
        [1.0, 0.9, 0.1],
        [0.9, 1.0, 0.0],
        [0.1, 0.0, 1.0],
    ]))
    emb = np.array([[0.0, 0.0], [1.0, 0.0], [10.0, 0.0]])
    nearest = np.zeros(3)

    # k=1 → only the single weakest edge of cell 0, which is cell 2 at d=10.
    out = _snn_bandwidth(snn, emb, k=1, nearest_dist=nearest)
    assert np.isclose(out[0], 10.0)

    # nearest_dist is subtracted off, ReLU'd.
    out2 = _snn_bandwidth(snn, emb, k=1, nearest_dist=np.array([4.0, 0.0, 0.0]))
    assert np.isclose(out2[0], 6.0)


def test_snn_bandwidth_includes_ties_then_keeps_largest():
    """Ties at the k-th weight all enter; the k largest distances then win."""
    # Cell 0 has three edges all of equal weight, at distances 1, 2 and 9.
    snn = sp.csc_matrix(np.array([
        [1.0, 0.5, 0.5, 0.5],
        [0.5, 1.0, 0.0, 0.0],
        [0.5, 0.0, 1.0, 0.0],
        [0.5, 0.0, 0.0, 1.0],
    ]))
    emb = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [9.0, 0.0]])
    # k=1, but the three-way tie pulls all three in; the largest (9) survives.
    out = _snn_bandwidth(snn, emb, k=1, nearest_dist=np.zeros(4))
    assert np.isclose(out[0], 9.0)


# ---------------------------------------------------------------------------
# Joint neighbour search (R: MultiModalNN)
# ---------------------------------------------------------------------------

def test_joint_graph_is_symmetric_and_binary():
    obj, labels = _complementary_object()
    find_multi_modal_neighbors(
        obj, reduction_list=["pca", "apca"],
        dims_list=[range(15), range(10)], k_nn=20,
    )
    wknn = obj.graphs["wknn"].tocsr()
    assert (abs(wknn - wknn.T) > 1e-9).nnz == 0     # union symmetrisation
    assert set(np.unique(wknn.data)) <= {1.0}       # binary adjacency
    assert np.allclose(wknn.diagonal(), 1.0)        # self included


def test_joint_neighbours_beat_either_modality_alone():
    """The joint ranking must not just reproduce one modality's neighbours."""
    obj, labels = _complementary_object()
    embs = [
        _l2_norm(obj.reductions["pca"].cell_embeddings[:, :15]),
        _l2_norm(obj.reductions["apca"].cell_embeddings[:, :10]),
    ]
    n = embs[0].shape[0]
    weights, sigmas, nearest = _modality_weights(
        embs, k_nn=20, sd_scale=1.0, cross_constant=1e-4, smooth=False, seed=42,
    )
    joint = _multi_modal_nn(embs, weights, sigmas, nearest, k_nn=20,
                            knn_range=50, seed=42)

    assert joint.shape == (n, 20)
    assert (joint != np.arange(n)[:, None]).all()   # self is never a neighbour

    # Joint neighbours should share a label more often than RNA's do, since RNA
    # alone cannot tell B from C.
    rna_only, _ = _build_knn(embs[0], 21, 42)
    same_joint = (labels[joint] == labels[:, None]).mean()
    same_rna = (labels[rna_only[:, 1:]] == labels[:, None]).mean()
    assert same_joint > same_rna
