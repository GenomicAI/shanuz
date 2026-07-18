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
from shanuz.reduction import run_pca, _default_features  # noqa: E402
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


def test_leverage_score_matches_seurats_truncated_svd():
    """Seurat's definition, not the textbook one — see the regression test below."""
    obj, _ = _clustered_object()

    # R's small-data branch is rowSums(irlba(data, nv = 50)$v ^ 2): the leading 50
    # right singular vectors of the (features × cells) *data* layer, uncentered.
    assay = obj.get_assay()
    feats = _default_features(assay, None)
    A = assay.layer_data("data")
    idx = [assay.features().index(f) for f in feats]
    A = np.asarray(A[idx, :].todense() if sp.issparse(A) else np.asarray(A)[idx, :],
                   dtype=float)
    _, _, vt = np.linalg.svd(A, full_matrices=False)
    k = min(50, min(A.shape) - 1)
    expected = (vt[:k] ** 2).sum(axis=0)

    scores = leverage_score(obj)

    assert np.allclose(scores, expected, atol=1e-8)
    assert scores.min() >= -1e-9
    # Also written back to metadata.
    assert "leverage.score" in obj.meta_data.columns


def test_leverage_scores_sum_to_the_truncation_not_the_rank():
    """The regression test for the defect this tutorial found.

    Seurat truncates at 50 components, so the scores sum to 50. Shanuz used to
    whiten against the *full* rank, which sums to the rank instead — on a matrix
    with more features than 50 that is a completely different, and far flatter,
    distribution. Asserting the sum pins the truncation directly.
    """
    obj, _ = _clustered_object(sizes=(200, 200), G=600, nfeatures=400)

    scores = leverage_score(obj)

    assert scores.sum() == pytest.approx(50.0, abs=1e-6)
    assert scores.sum() != pytest.approx(400.0, abs=1.0)   # the old, wrong answer


def test_leverage_scores_are_not_flat():
    """The defect's *symptom*: near-uniform weights make sketching pointless.

    The old full-rank scores on a realistic matrix had a max/median of about 1.3
    against R's 6.5 — a spread so small that leverage sampling was uniform
    sampling with extra steps. Any implementation that saturates fails here.

    Thresholds are placed in the *gap*, not on either value. Measured on this
    fixture, whose 400 cells and 400 features make the old failure mode maximal
    (full rank saturates every score at exactly 1):

        ================  ==========  ============
                          max/median  CV
        ================  ==========  ============
        old (full rank)   1.000       0.0000
        fixed (rank 50)   1.525       0.1627
        ================  ==========  ============

    An earlier version of this test asserted ``> 1.5``, which sat directly on the
    fixed value and duly failed on Python 3.10 at 1.4956 while passing on 3.11
    and 3.12 — the SVD drifts ~2 % across versions. 1.25 leaves ~20 % either way.
    If this ever fails again, check that the *sum* test above still passes before
    touching the number: a genuine regression shows up there first.
    """
    obj, _ = _clustered_object(sizes=(200, 200), G=600, nfeatures=400)

    scores = leverage_score(obj)

    assert scores.max() / np.median(scores) > 1.25
    assert scores.std() / scores.mean() > 0.08


def test_leverage_score_picks_the_regime_the_way_seurat_does():
    """R switches to the sketched branch at ``ncells >= nsketch * 1.5``."""
    # Many more cells than features: sketching needs the matrix to be tall, and
    # the sketched branch additionally needs nsketch above 1.1x the feature count,
    # so a square fixture cannot reach that branch at all.
    obj, _ = _clustered_object(sizes=(500, 500), G=150, nfeatures=100)
    n = len(obj)                                     # 1000 cells, 100 features

    # ncells < nsketch * 1.5 → exact, so the scores sum to 50.
    exact = leverage_score(obj, nsketch=n, var_name=None)
    assert exact.sum() == pytest.approx(50.0, abs=1e-6)

    # ncells >= nsketch * 1.5 → the CountSketch/QR/JL path, which is on its own
    # scale (Seurat leaves the JL projection unscaled) and must *not* sum to 50.
    sketched = leverage_score(obj, nsketch=200, var_name=None)
    assert sketched.shape == exact.shape
    assert np.all(np.isfinite(sketched)) and sketched.min() >= 0
    assert sketched.sum() != pytest.approx(50.0, abs=1.0)


def test_leverage_score_bumps_nsketch_like_seurat():
    """Below 1.1x the feature count Seurat raises nsketch and warns, not errors."""
    obj, _ = _clustered_object(sizes=(500, 500), G=150, nfeatures=100)

    with pytest.warns(UserWarning, match="too close to the number of features"):
        scores = leverage_score(obj, nsketch=105, var_name=None)

    assert np.all(np.isfinite(scores)) and scores.min() >= 0


def test_leverage_score_refuses_a_square_matrix():
    obj, _ = _clustered_object(sizes=(60, 60), G=200, nfeatures=150)
    with pytest.raises(ValueError, match="too square"):
        leverage_score(obj, nsketch=50, var_name=None)


def test_leverage_score_flags_rare_cells():
    # 120 common cells, 6 rare cells in a distinct state.
    obj, ct = _clustered_object(sizes=(120, 6), G=120, nfeatures=80)
    scores = leverage_score(obj)

    assert scores[ct == "T1"].mean() > scores[ct == "T0"].mean()


def test_leverage_score_defaults_to_the_data_layer():
    """Seurat scores the log-normalized data, not scale.data.

    Worth pinning: the two give materially different geometries, and the default
    silently deciding which one you get is exactly how the original defect stayed
    invisible.
    """
    obj, _ = _clustered_object()

    default = leverage_score(obj, var_name=None)
    explicit = leverage_score(obj, layer="data", var_name=None)
    scaled = leverage_score(obj, layer="scale.data", var_name=None)

    assert np.allclose(default, explicit)
    assert not np.allclose(default, scaled)


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


def test_sketch_data_draws_towards_high_leverage_cells():
    """The mechanism, tested where it is actually deterministic.

    Deliberately *not* "the sketch over-represents rare cell type T1". Synthetic
    Poisson clusters do not reproduce that: real rare types (pDC, erythrocytes)
    are transcriptionally extreme, not merely scarce, and on ifnb they do get
    2-3x the average leverage — which is where the tutorial's smoke test checks
    it. What is true on any data is that sampling *proportional to leverage*
    favours high-leverage cells, and that is what this pins, averaged over seeds
    so a single unlucky draw cannot flip it.
    """
    obj, _ = _clustered_object(sizes=(300, 300), G=400, nfeatures=300)
    scores = leverage_score(obj, var_name=None)
    order = {c: i for i, c in enumerate(obj.cell_names())}
    n = len(obj)

    def mean_leverage(method):
        got = []
        for seed in range(20):
            sk = sketch_data(obj, ncells=n // 10, method=method, seed=seed)
            got.append(scores[[order[c] for c in sk.cell_names()]].mean())
        return float(np.mean(got))

    weighted = mean_leverage("LeverageScore")
    uniform = mean_leverage("Uniform")

    assert weighted > uniform
    assert weighted > scores.mean()
    # Uniform must land on the population mean — if it does not, the sampling
    # itself is biased and the comparison above proves nothing.
    assert uniform == pytest.approx(scores.mean(), rel=0.02)


def test_sketch_data_uniform_method_ignores_leverage():
    obj, _ = _clustered_object(sizes=(100, 100))
    sk = sketch_data(obj, ncells=50, method="Uniform", var_name="lev", seed=0)

    assert len(sk) == 50
    assert np.allclose(obj.meta_data["lev"].to_numpy(), 1.0)


def test_sketch_data_rejects_unknown_method():
    obj, _ = _clustered_object(sizes=(40, 40))
    with pytest.raises(ValueError, match="Uniform"):
        sketch_data(obj, ncells=10, method="Nonsense")


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


def test_project_data_votes_over_neighbours_not_anchors():
    """Seurat's ``ProjectData`` transfers by weighted kNN, not integration anchors.

    Pinned because the anchor route this replaced was *more* accurate on ifnb
    (0.936 vs R's 0.905) while being the wrong algorithm — and unusable at the
    scale sketching targets, since anchor finding against the full dataset is the
    cost sketching removes. Matching R took accuracy *down* to 0.903, at 98.1 %
    per-cell agreement. A regression would most likely look like an improvement,
    so the check is structural: no anchor object is built, and the scores are a
    normalized vote (bounded by 1) over exactly ``k_weight`` neighbours.
    """
    import shanuz.transfer as transfer_mod

    obj, ct = _clustered_object(sizes=(80, 80, 80))
    sk = sketch_data(obj, ncells=90, seed=0)
    run_pca(sk, n_pcs=20)

    calls = []
    original = transfer_mod.find_transfer_anchors
    transfer_mod.find_transfer_anchors = lambda *a, **k: calls.append(1)
    try:
        project_data(obj, sk, refdata={"predicted": "celltype"}, project_umap=False)
    finally:
        transfer_mod.find_transfer_anchors = original

    assert calls == [], "project_data must not build transfer anchors"
    scores = obj.meta_data["predicted.score"].to_numpy()
    assert scores.min() > 0.0 and scores.max() <= 1.0 + 1e-9
    acc = np.mean(obj.meta_data["predicted"].to_numpy() == ct)
    assert acc > 0.85


def test_project_data_accepts_a_bare_column_name():
    """R expands a bare string to {col: col}; a raw array is not accepted."""
    obj, _ = _clustered_object(sizes=(60, 60))
    sk = sketch_data(obj, ncells=60, seed=0)
    run_pca(sk, n_pcs=10)

    project_data(obj, sk, refdata="celltype", project_umap=False)

    assert "celltype.score" in obj.meta_data.columns
    with pytest.raises(KeyError, match="not found in the sketch"):
        project_data(obj, sk, refdata="no_such_column", project_umap=False)


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
