"""Tests for leverage-score sketching (v0.8.0 Scale).

  * leverage_score  (sketch.py)
  * sketch_data     (sketch.py)
  * project_data    (sketch.py)

Builds small clustered datasets (each cell type has its own elevated gene block)
and checks that leverage scores match the textbook definition, that sketching
draws a leverage-weighted subset that over-represents the rare states, and that
``project_data`` places every full-dataset cell back into the sketch's PCA/UMAP
and transfers the sketch's labels to the full data. Network-free.
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
from shanuz.reduction import run_pca, _default_features, _get_scaled_data  # noqa: E402
from shanuz.umap import run_umap  # noqa: E402
from shanuz.sketch import leverage_score, sketch_data, project_data  # noqa: E402


def _clustered_object(sizes=(60, 60, 60), seed=0, G=150, nfeatures=100):
    """Cells in K clusters; cluster k carries an elevated gene block."""
    rng = np.random.default_rng(seed)
    K = len(sizes)
    block = G // (K + 1)
    cols, celltype, cells = [], [], []
    c = 0
    for k, nk in enumerate(sizes):
        for _ in range(nk):
            base = rng.gamma(0.3, size=G) + 0.05
            base[k * block:(k + 1) * block] += 5.0
            cols.append(rng.poisson(base * 3000.0 / base.sum()))
            celltype.append(f"T{k}")
            cells.append(f"c{c}")
            c += 1

    mat = np.array(cols).T  # (G features × n cells)
    meta = pd.DataFrame({"celltype": celltype}, index=cells)
    obj = create_shanuz_object(
        counts=sp.csc_matrix(mat), assay="RNA",
        feature_names=[f"g{i}" for i in range(G)], cell_names=cells,
        meta_data=meta,
    )
    normalize_data(obj)
    find_variable_features(obj, selection_method="vst", nfeatures=nfeatures)
    scale_data(obj)
    return obj, np.array(celltype)


# ----------------------------------------------------------------------
# leverage_score
# ----------------------------------------------------------------------


def test_leverage_score_matches_textbook_definition():
    obj, _ = _clustered_object()

    # Build the exact matrix the function scores on, and its textbook leverage
    # (row norms of the left singular vectors).
    assay = obj.get_assay()
    feats = _default_features(assay, None)
    A = np.asarray(_get_scaled_data(assay, feats, "scale.data"), dtype=float).T
    U, s, _ = np.linalg.svd(A, full_matrices=False)
    tol = s[0] * max(A.shape) * 1e-8
    r = int((s > tol).sum())
    exact = (U[:, :r] ** 2).sum(axis=1)

    scores = leverage_score(obj, nsketch=10 ** 9)  # nsketch >= n → exact

    assert np.allclose(scores, exact, atol=1e-6)
    # Leverage scores are non-negative, each at most 1, and sum to the rank.
    assert scores.min() >= -1e-9
    assert scores.max() <= 1 + 1e-8
    assert scores.sum() == pytest.approx(r, abs=1e-6)
    # Also written back to metadata.
    assert "leverage.score" in obj.meta_data.columns


def test_leverage_score_flags_rare_cells():
    # 120 common cells, 6 rare cells in a distinct state.
    obj, ct = _clustered_object(sizes=(120, 6), G=120, nfeatures=80)
    scores = leverage_score(obj, nsketch=10 ** 9)

    assert scores[ct == "T1"].mean() > scores[ct == "T0"].mean()


def test_leverage_score_sketch_approximates_exact():
    # n = 810 cells, d = 20 features: nsketch=500 sits between d and n, so the
    # CountSketch path is exercised. A single-hash sketch is a subspace embedding
    # only once it has ~d² rows, so this uses nsketch ≈ d² — the regime where the
    # approximation tracks the exact scores (as it does at nsketch=5000 for the
    # million-cell case this is a proxy for).
    obj, _ = _clustered_object(sizes=(270, 270, 270), G=80, nfeatures=20)
    exact = leverage_score(obj, nsketch=10 ** 9, seed=0)
    approx = leverage_score(obj, nsketch=500, seed=0)

    assert np.corrcoef(exact, approx)[0, 1] > 0.9


# ----------------------------------------------------------------------
# sketch_data
# ----------------------------------------------------------------------


def test_sketch_data_returns_subset():
    obj, _ = _clustered_object()
    sk = sketch_data(obj, ncells=40, seed=0)

    assert len(sk) == 40
    assert set(sk.cell_names()) <= set(obj.cell_names())
    # Scores stashed on the source object; the sketch's assay is renamed.
    assert "leverage.score" in obj.meta_data.columns
    assert sk.active_assay == "sketch"
    assert "sketch" in sk.assays
    assert sk.misc["sketch"]["ncells"] == 40


def test_sketch_data_oversamples_rare():
    # 10% of cells are a distinct rare state; leverage sampling should enrich it.
    obj, ct = _clustered_object(sizes=(180, 20), G=120, nfeatures=80)
    sk = sketch_data(obj, ncells=60, seed=0)

    sketch_types = obj.meta_data.loc[sk.cell_names(), "celltype"].to_numpy()
    rare_full = float(np.mean(ct == "T1"))
    rare_sketch = float(np.mean(sketch_types == "T1"))
    assert rare_sketch > rare_full


def test_sketch_data_caps_at_available_cells():
    obj, _ = _clustered_object(sizes=(30, 30))
    sk = sketch_data(obj, ncells=10_000, seed=0)

    assert len(sk) == len(obj)


# ----------------------------------------------------------------------
# project_data
# ----------------------------------------------------------------------


def test_project_data_projects_full_into_sketch_pca():
    obj, _ = _clustered_object(sizes=(80, 80, 80))
    sk = sketch_data(obj, ncells=60, seed=0)
    run_pca(sk, n_pcs=20)

    project_data(obj, sk, project_umap=False)

    assert "pca.full" in obj.reductions
    full_emb = obj.reductions["pca.full"].cell_embeddings
    assert full_emb.shape == (len(obj), 20)

    # Projecting the sketch cells through the sketch's own loadings reproduces
    # their sketch-PCA coordinates (up to a per-PC offset → correlation ≈ 1).
    idx = [obj.cell_names().index(c) for c in sk.cell_names()]
    sk_emb = sk.reductions["pca"].cell_embeddings
    for d in range(3):
        r = abs(np.corrcoef(full_emb[idx, d], sk_emb[:, d])[0, 1])
        assert r > 0.99


def test_project_data_transfers_labels():
    obj, ct = _clustered_object(sizes=(80, 80, 80))
    sk = sketch_data(obj, ncells=90, seed=0)
    run_pca(sk, n_pcs=20)

    project_data(obj, sk, refdata={"predicted": "celltype"}, project_umap=False)

    assert "predicted" in obj.meta_data.columns
    assert "predicted.score" in obj.meta_data.columns
    acc = np.mean(obj.meta_data["predicted"].to_numpy() == ct)
    assert acc > 0.85


def test_project_data_projects_umap_when_model_present():
    obj, _ = _clustered_object(sizes=(80, 80, 80))
    sk = sketch_data(obj, ncells=90, seed=0)
    run_pca(sk, n_pcs=20)
    run_umap(sk, n_neighbors=15, min_dist=0.3, seed=42)

    project_data(obj, sk)  # project_umap=True by default

    assert "ref.umap" in obj.reductions
    assert obj.reductions["ref.umap"].cell_embeddings.shape == (len(obj), 2)


def test_project_data_missing_sketch_pca_raises():
    obj, _ = _clustered_object(sizes=(40, 40))
    sk = sketch_data(obj, ncells=40, seed=0)  # no run_pca on the sketch
    with pytest.raises(KeyError):
        project_data(obj, sk)
