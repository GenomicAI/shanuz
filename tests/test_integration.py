"""Tests for batch correction & integration (v0.2.0).

  * run_harmony        (integration.py)
  * integrate_layers   (integration.py)

Builds a small synthetic dataset with a clear cell-type structure plus an
injected batch effect, then checks that Harmony mixes the batches (lower
per-batch silhouette) while preserving cell-type separation. Network-free.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz.shanuz import create_shanuz_object  # noqa: E402
from shanuz.preprocessing import (  # noqa: E402
    normalize_data,
    find_variable_features,
    scale_data,
)
from shanuz.reduction import run_pca  # noqa: E402
from shanuz.integration import run_harmony, integrate_layers  # noqa: E402

pytest.importorskip("harmonypy")
from sklearn.metrics import silhouette_score  # noqa: E402


def _batched_object(seed=0):
    """Two cell types × two batches; genes 100-149 carry the batch effect."""
    rng = np.random.default_rng(seed)
    G, per = 200, 50
    mat = np.zeros((G, 4 * per))
    celltype, batch, cells = [], [], []
    c = 0
    for t in ("A", "B"):
        for b in ("1", "2"):
            for _ in range(per):
                base = rng.gamma(0.3, size=G) + 0.05
                if t == "A":
                    base[0:50] += 5.0
                else:
                    base[50:100] += 5.0
                if b == "2":
                    base[100:150] += 4.0          # batch-specific gene block
                mat[:, c] = rng.poisson(base * 3000.0 / base.sum())
                celltype.append(t)
                batch.append(b)
                cells.append(f"c{c}")
                c += 1

    meta = pd.DataFrame({"batch": batch, "celltype": celltype}, index=cells)
    obj = create_shanuz_object(
        counts=sp.csc_matrix(mat), assay="RNA",
        feature_names=[f"g{i}" for i in range(G)], cell_names=cells,
        meta_data=meta,
    )
    normalize_data(obj)
    find_variable_features(obj, selection_method="vst", nfeatures=150)
    scale_data(obj)
    run_pca(obj, n_pcs=20)
    return obj


def test_run_harmony_mixes_batches_and_preserves_celltype():
    obj = _batched_object()
    batch = obj.meta_data["batch"].to_numpy()
    celltype = obj.meta_data["celltype"].to_numpy()

    run_harmony(obj, group_by="batch", max_iter_harmony=20, seed=0)

    assert "harmony" in obj.reductions
    pca = obj.reductions["pca"].cell_embeddings
    harm = obj.reductions["harmony"].cell_embeddings
    assert harm.shape == pca.shape
    assert np.isfinite(harm).all()

    # Batches should be MORE mixed after Harmony (lower batch silhouette).
    sil_batch_pca = silhouette_score(pca, batch)
    sil_batch_harm = silhouette_score(harm, batch)
    assert sil_batch_harm < sil_batch_pca

    # Cell-type structure should survive (types still separable).
    assert silhouette_score(harm, celltype) > 0.0


def test_integrate_layers_harmony_matches_run_harmony():
    obj = _batched_object()
    integrate_layers(obj, method="harmony", group_by="batch",
                     new_reduction="int", max_iter_harmony=20, seed=0)
    assert "int" in obj.reductions
    assert obj.reductions["int"].cell_embeddings.shape == \
        obj.reductions["pca"].cell_embeddings.shape


def test_run_harmony_unknown_group_by_raises():
    obj = _batched_object()
    with pytest.raises(KeyError):
        run_harmony(obj, group_by="nonexistent_column")


# ----------------------------------------------------------------------
# CCA / RPCA anchor integration
# ----------------------------------------------------------------------

from shanuz.anchors import (  # noqa: E402
    find_integration_anchors,
    integrate_data,
    IntegrationAnchors,
)


def _one_batch_object(batch, seed=0, n_per=50, G=200):
    """A single dataset: two cell types; batch '2' carries a batch-effect block."""
    rng = np.random.default_rng(seed)
    mat = np.zeros((G, 2 * n_per))
    celltype, cells = [], []
    c = 0
    for t in ("A", "B"):
        for _ in range(n_per):
            base = rng.gamma(0.3, size=G) + 0.05
            if t == "A":
                base[0:50] += 5.0
            else:
                base[50:100] += 5.0
            if batch == "2":
                base[100:150] += 4.0          # batch-specific gene block
            mat[:, c] = rng.poisson(base * 3000.0 / base.sum())
            celltype.append(t)
            cells.append(f"b{batch}_c{c}")
            c += 1

    meta = pd.DataFrame({"celltype": celltype, "batch": batch}, index=cells)
    obj = create_shanuz_object(
        counts=sp.csc_matrix(mat), assay="RNA",
        feature_names=[f"g{i}" for i in range(G)], cell_names=cells,
        meta_data=meta,
    )
    normalize_data(obj)
    find_variable_features(obj, selection_method="vst", nfeatures=150)
    scale_data(obj)
    return obj, np.array(celltype)


def _object_pair(seed=0):
    ref, ct_ref = _one_batch_object("1", seed=seed)
    query, ct_query = _one_batch_object("2", seed=seed + 1)
    celltype = np.concatenate([ct_ref, ct_query])
    return [ref, query], celltype


def test_find_integration_anchors_returns_scored_anchors():
    objs, _ = _object_pair()
    anchors = find_integration_anchors(objs, reduction="cca", k_filter=0)
    assert isinstance(anchors, IntegrationAnchors)
    df = anchors.anchors
    assert list(df.columns) == ["dataset1", "cell1", "dataset2", "cell2", "score"]
    assert len(df) > 0
    assert (df["dataset1"] == 0).all()
    assert (df["dataset2"] == 1).all()
    assert df["score"].between(0.0, 1.0).all()


def test_find_integration_anchors_rejects_unknown_reduction():
    objs, _ = _object_pair()
    with pytest.raises(ValueError):
        find_integration_anchors(objs, reduction="mnn")


def test_find_integration_anchors_needs_two_objects():
    objs, _ = _object_pair()
    with pytest.raises(ValueError):
        find_integration_anchors(objs[:1])


def test_find_integration_anchors_is_deterministic():
    objs, _ = _object_pair()
    a1 = find_integration_anchors(objs, reduction="cca")
    a2 = find_integration_anchors(objs, reduction="cca")
    pd.testing.assert_frame_equal(a1.anchors, a2.anchors)


def _integrated_pca(merged, seed=0):
    scale_data(merged, assay="integrated")
    run_pca(merged, n_pcs=15, assay="integrated", seed=seed)
    return merged.reductions["pca"].cell_embeddings


def test_integrate_data_clusters_by_celltype_not_batch():
    objs, celltype = _object_pair()
    batch = np.array(["1"] * len(objs[0]) + ["2"] * len(objs[1]))

    anchors = find_integration_anchors(objs, reduction="cca")
    merged = integrate_data(anchors)

    assert "integrated" in merged.assays
    assert merged.active_assay == "integrated"

    emb = _integrated_pca(merged)
    # After anchor correction, cell type should dominate over batch.
    assert silhouette_score(emb, celltype) > silhouette_score(emb, batch)


def test_integrate_data_rpca_clusters_by_celltype():
    objs, celltype = _object_pair()
    batch = np.array(["1"] * len(objs[0]) + ["2"] * len(objs[1]))

    anchors = find_integration_anchors(objs, reduction="rpca")
    merged = integrate_data(anchors)
    emb = _integrated_pca(merged)
    assert silhouette_score(emb, celltype) > silhouette_score(emb, batch)


def test_integrate_data_leaves_the_reference_untouched():
    objs, _ = _object_pair()
    anchors = find_integration_anchors(objs, reduction="cca")
    features = anchors.anchor_features
    merged = integrate_data(anchors)

    ref = objs[0]
    ref_cells = ref.cell_names()
    got = merged.get_assay("integrated").layer_data(
        "data", features=features, cells=ref_cells
    )
    orig = ref.get_assay().layer_data("data", features=features, cells=ref_cells)
    orig = orig.toarray() if sp.issparse(orig) else np.asarray(orig)
    got = got.toarray() if sp.issparse(got) else np.asarray(got)
    assert np.allclose(got, orig)


def test_integrate_layers_cca_mixes_batches():
    obj = _batched_object()
    batch = obj.meta_data["batch"].to_numpy()
    celltype = obj.meta_data["celltype"].to_numpy()

    integrate_layers(obj, method="cca", group_by="batch", new_reduction="int_cca")

    assert "int_cca" in obj.reductions
    emb = obj.reductions["int_cca"].cell_embeddings
    assert emb.shape[0] == len(obj)
    assert np.isfinite(emb).all()
    # Integrated embedding should separate cell types more than batches.
    assert silhouette_score(emb, celltype) > silhouette_score(emb, batch)


def test_integrate_layers_rpca_produces_a_reduction():
    obj = _batched_object()
    integrate_layers(obj, method="rpca", group_by="batch", new_reduction="int_rpca")
    assert "int_rpca" in obj.reductions
    assert obj.reductions["int_rpca"].cell_embeddings.shape[0] == len(obj)


def test_integrate_layers_cca_requires_group_by():
    obj = _batched_object()
    with pytest.raises(ValueError):
        integrate_layers(obj, method="cca")
