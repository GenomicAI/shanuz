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


def test_integrate_layers_cca_not_implemented():
    obj = _batched_object()
    with pytest.raises(NotImplementedError):
        integrate_layers(obj, method="cca", group_by="batch")


def test_run_harmony_unknown_group_by_raises():
    obj = _batched_object()
    with pytest.raises(KeyError):
        run_harmony(obj, group_by="nonexistent_column")
