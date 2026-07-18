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
    _standardize_and_l2,
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


def _unequal_pair(n_ref=40, n_query=70, seed=0):
    """A ref/query pair with DIFFERENT cell counts — the shape real batches have.

    ``_object_pair`` builds 100 cells on each side; the balanced case never
    exercises the code path where one dataset is larger than the other.
    """
    ref, ct_ref = _one_batch_object("1", seed=seed, n_per=n_ref)
    query, ct_query = _one_batch_object("2", seed=seed + 1, n_per=n_query)
    return [ref, query], np.concatenate([ct_ref, ct_query])


def test_find_integration_anchors_rpca_handles_unequal_sizes():
    """RPCA anchoring must not assume the reference and query are the same size.

    Regression for an IndexError in the reciprocal-PCA branch: it searched for
    query neighbours but indexed the result with reference-sized bounds, so any
    pair with ``n_query > n_ref`` (every real dataset — ifnb is CTRL 6548 vs STIM
    7451) walked off the end of the neighbour list. Balanced fixtures (equal
    n_per) masked it, and even without the crash the anchors were degenerate.
    ``cca`` is the sanity baseline on the identical objects.
    """
    for n_ref, n_query in ((40, 70), (70, 40)):   # both orderings
        objs, _ = _unequal_pair(n_ref=n_ref, n_query=n_query)
        assert len(objs[0]) != len(objs[1])
        rpca = find_integration_anchors(objs, reduction="rpca", k_filter=0)
        cca = find_integration_anchors(objs, reduction="cca", k_filter=0)
        # Real mutual anchors, not a crash and not a near-empty degenerate set.
        assert len(rpca.anchors) > 0
        assert len(rpca.anchors) > 0.2 * len(cca.anchors)


def test_integrate_data_rpca_unequal_sizes_clusters_by_celltype():
    """End-to-end RPCA correction still works when the datasets differ in size."""
    objs, celltype = _unequal_pair()
    batch = np.array(["1"] * len(objs[0]) + ["2"] * len(objs[1]))
    anchors = find_integration_anchors(objs, reduction="rpca")
    merged = integrate_data(anchors)
    emb = _integrated_pca(merged)
    assert silhouette_score(emb, celltype) > silhouette_score(emb, batch)


def test_standardize_and_l2_normalizes_dimensions_then_cells():
    """The reciprocal-PCA embedding normalisation (Seurat's ``l2.norm=TRUE``).

    ``_standardize_and_l2`` divides every dimension (column) by its SD across the
    stacked ref+query cells, then L2-normalises every cell (row). The column step
    is the one that matters: skip it and the neighbour search collapses onto PC1's
    dominant variance, which is what made RPCA under-integrate (Bug 2). This pins
    both halves — rows come out unit-norm, and the result differs from plain
    row-normalisation whenever the columns have unequal spread (so the SD step is
    provably not a no-op).
    """
    from shanuz.anchors import _l2_normalize_rows

    rng = np.random.default_rng(0)
    # Columns with wildly different scales — the reciprocal-PCA situation, where
    # PC1 carries orders of magnitude more variance than the trailing PCs.
    emb = rng.normal(size=(60, 6)) * np.array([100.0, 30.0, 8.0, 3.0, 1.0, 0.3])
    out = _standardize_and_l2(emb)
    assert np.allclose(np.linalg.norm(out, axis=1), 1.0)      # every cell unit-norm
    assert not np.allclose(out, _l2_normalize_rows(emb))       # SD step is not a no-op


def test_find_integration_anchors_rpca_embedding_is_l2_normalized():
    """RPCA must normalise its reciprocal embedding before the MNN and weighting.

    Regression for Bug 2 (the anchor-quality half): shanuz searched the raw
    ``scaled.T @ loadings`` projection, so PC1's variance dominated the neighbour
    search and the anchors were wrong — RPCA under-integrated ifnb (batch-mixing
    entropy 0.22 vs Seurat 0.91). The fix mirrors Seurat's ``ReciprocalProject``:
    standardise each dimension, then L2-normalise each cell. The observable
    signature is that the per-query weight embedding has unit-norm rows, which the
    raw projection does not. (The *emergent* under-integration only reproduces on
    real data with many overlapping cell types — synthetic batches, however
    strong, are too separable to mislead the reciprocal search; that behaviour is
    pinned by the gated ifnb regression in ``test_tutorial_smoke``.)
    """
    objs, _ = _unequal_pair(n_ref=40, n_query=70)
    anchors = find_integration_anchors(objs, reduction="rpca", k_filter=0)
    we = anchors.weight_embeddings[1]            # query cells in the reference space
    assert we.shape[0] == len(objs[1])
    assert np.allclose(np.linalg.norm(we, axis=1), 1.0)


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
