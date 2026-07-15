"""Tests for reference mapping, part two (v0.3.0).

  * project_umap  (mapping.py)
  * map_query     (mapping.py)

Builds an annotated reference with a fitted PCA + UMAP and an unlabelled query
(same two cell types, the query carrying a batch-effect gene block the reference
never sees), then checks that projecting the query lands each cell in the right
neighbourhood of the *reference's* UMAP and that ``map_query`` both annotates and
places the query in one call. Network-free.
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
from shanuz.umap import run_umap  # noqa: E402
from shanuz.dimreduc import DimReduc  # noqa: E402
from shanuz.transfer import find_transfer_anchors  # noqa: E402
from shanuz.mapping import map_query, project_umap  # noqa: E402


def _labelled_object(batch, seed=0, n_per=40, G=200):
    """Two cell types (A/B); batch '2' carries a batch-effect gene block."""
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


def _ref_query(seed=0):
    """Reference (with a fitted PCA + UMAP) and an unlabelled query.

    Both are scaled on a shared feature set so the query carries every feature
    the reference's PCA is built on.
    """
    reference, ct_ref = _labelled_object("1", seed=seed)
    query, ct_query = _labelled_object("2", seed=seed + 1)

    common = [
        f for f in reference.get_assay().variable_features
        if f in set(query.get_assay().variable_features)
    ]
    scale_data(reference, features=common)
    scale_data(query, features=common)
    run_pca(reference, n_pcs=20, features=common)
    run_umap(reference, n_neighbors=15, min_dist=0.3, seed=42)
    return reference, query, ct_ref, ct_query


def _umap_centroid_accuracy(query_coords, ref_coords, ct_ref, ct_query):
    """Classify each query cell by the nearest reference-type UMAP centroid."""
    classes = sorted(set(ct_ref))
    centroids = {c: ref_coords[ct_ref == c].mean(axis=0) for c in classes}
    preds = []
    for row in query_coords:
        d = {c: np.linalg.norm(row - centroids[c]) for c in classes}
        preds.append(min(d, key=d.get))
    return float(np.mean(np.asarray(preds) == np.asarray(ct_query)))


# ----------------------------------------------------------------------
# project_umap
# ----------------------------------------------------------------------


def test_project_umap_stores_reference_embedding():
    reference, query, _, _ = _ref_query()

    dr = project_umap(query, reference)

    assert isinstance(dr, DimReduc)
    assert "ref.umap" in query.reductions
    assert dr.cell_embeddings.shape == (len(query), 2)
    assert dr.cells() == query.cell_names()
    # The intermediate reference-PCA projection is kept for introspection.
    assert dr.misc["query_pca"].shape[0] == len(query)


def test_project_umap_places_query_near_matching_reference_type():
    reference, query, ct_ref, ct_query = _ref_query()
    dr = project_umap(query, reference)

    ref_umap = reference.reductions["umap"].cell_embeddings
    acc = _umap_centroid_accuracy(dr.cell_embeddings, ref_umap, ct_ref, ct_query)
    # Despite the query's batch block, projection lands cells with the right type.
    assert acc > 0.85


def test_project_umap_missing_pca_raises():
    reference, query, _, _ = _ref_query()
    del reference.reductions["pca"]
    with pytest.raises(KeyError):
        project_umap(query, reference)


def test_project_umap_requires_fitted_model():
    reference, query, _, _ = _ref_query()
    # A UMAP reduction with no stored model (e.g. embedded from a graph).
    reference.reductions["umap"] = DimReduc(
        cell_embeddings=reference.reductions["umap"].cell_embeddings,
        cell_names=reference.cell_names(),
        key="UMAP_",
        misc={"umap_graph": "wsnn"},
    )
    with pytest.raises(ValueError):
        project_umap(query, reference)


# ----------------------------------------------------------------------
# map_query
# ----------------------------------------------------------------------


def test_map_query_annotates_and_projects():
    reference, query, _, ct_query = _ref_query()
    anchors = find_transfer_anchors(reference, query, reduction="pcaproject")

    pred = map_query(anchors, refdata="celltype")

    # Labels transferred and written onto the query metadata.
    assert "predicted.id" in query.meta_data.columns
    assert "prediction.score.max" in query.meta_data.columns
    assert np.mean(query.meta_data["predicted.id"].to_numpy() == ct_query) > 0.85
    # Predictions are also returned.
    assert list(pred.index) == query.cell_names()
    # And the query is placed in the reference UMAP.
    assert "ref.umap" in query.reductions
    assert query.reductions["ref.umap"].cell_embeddings.shape == (len(query), 2)


def test_map_query_without_refdata_only_projects():
    reference, query, _, _ = _ref_query()
    anchors = find_transfer_anchors(reference, query)

    pred = map_query(anchors, refdata=None)

    assert pred is None
    assert "predicted.id" not in query.meta_data.columns
    assert "ref.umap" in query.reductions


def test_map_query_imputes_without_writing_metadata():
    reference, query, _, _ = _ref_query()
    anchors = find_transfer_anchors(reference, query)

    feats = ["g10", "g60"]
    ref_expr = reference.get_assay().layer_data("data", features=feats)
    ref_expr = ref_expr.toarray() if sp.issparse(ref_expr) else np.asarray(ref_expr)

    imputed = map_query(anchors, refdata=ref_expr, refdata_features=feats)

    # Imputation output is returned as a matrix, not written to metadata.
    assert list(imputed.index) == feats
    assert list(imputed.columns) == query.cell_names()
    assert "predicted.id" not in query.meta_data.columns
    # The UMAP projection still happens.
    assert "ref.umap" in query.reductions
