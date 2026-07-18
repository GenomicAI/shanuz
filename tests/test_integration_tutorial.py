"""Network-free tests for the ifnb integration tutorial (tutorials/ifnb_integration_tutorial.py).

Covers the pure metric helpers (silhouette mixing, cluster ARI, batch entropy,
scoreboard) directly, and drives the whole pipeline on a small synthetic
two-condition dataset with a planted batch effect — never touching the network
or the real ifnb download. The synthetic batches are deliberately *unequal*
(120 vs 150 cells) so the RPCA leg exercises the reciprocal-anchor path that a
balanced fixture would not.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import tutorials.ifnb_integration_tutorial as tut  # noqa: E402
from tutorials.ifnb_integration_tutorial import (  # noqa: E402
    mixing_metrics,
    cluster_ari,
    batch_entropy,
    build_scoreboard,
)


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

def test_mixing_metrics_separated_vs_mixed():
    rng = np.random.default_rng(0)
    n = 60
    batch = np.array(["A"] * n + ["B"] * n)
    # Batch-separated embedding: A near (0,0), B near (10,0).
    sep = np.vstack([rng.normal(0, 0.3, (n, 2)), rng.normal([10, 0], 0.3, (n, 2))])
    # Mixed embedding: both batches drawn from the same cloud.
    mixed = rng.normal(0, 1.0, (2 * n, 2))
    celltype = np.array((["X"] * (n // 2) + ["Y"] * (n // 2)) * 2)

    m_sep = mixing_metrics(sep, batch, celltype)
    m_mix = mixing_metrics(mixed, batch, celltype)
    # A separated batch has a much higher batch silhouette than a mixed one.
    assert m_sep["sil_batch"] > 0.5
    assert m_mix["sil_batch"] < m_sep["sil_batch"]
    assert set(m_sep) == {"sil_batch", "sil_celltype"}


def test_cluster_ari_identical_and_shuffled():
    a = np.array([0, 0, 1, 1, 2, 2])
    assert cluster_ari(a, a) == pytest.approx(1.0)
    # A relabeling is still a perfect partition match (ARI is label-invariant).
    assert cluster_ari(a, np.array([5, 5, 9, 9, 7, 7])) == pytest.approx(1.0)
    # An unrelated split scores far below 1.
    assert cluster_ari(a, np.array([0, 1, 0, 1, 0, 1])) < 0.5


def test_cluster_ari_rejects_bad_shapes():
    with pytest.raises(ValueError):
        cluster_ari(np.array([0, 1]), np.array([0, 1, 2]))
    with pytest.raises(ValueError):
        cluster_ari(np.array([]), np.array([]))


def test_batch_entropy_extremes():
    # Two clusters, each a single batch -> no mixing -> 0.
    clusters = np.array([0] * 10 + [1] * 10)
    batch = np.array(["A"] * 10 + ["B"] * 10)
    assert batch_entropy(clusters, batch) == pytest.approx(0.0)

    # Two clusters, each half A / half B -> perfect mixing -> 1.
    clusters2 = np.array([0, 1] * 10)
    batch2 = np.array(["A", "A", "B", "B"] * 5)
    assert batch_entropy(clusters2, batch2) == pytest.approx(1.0)

    # A single batch level is undefined (no mixing to measure).
    assert np.isnan(batch_entropy(clusters, np.array(["A"] * 20)))


def test_batch_entropy_partial_between_zero_and_one():
    # 70/30 split in every cluster -> partial mixing, strictly inside (0, 1).
    rng = np.random.default_rng(1)
    clusters = np.repeat(np.arange(4), 50)
    batch = np.where(rng.random(200) < 0.7, "A", "B")
    e = batch_entropy(clusters, batch)
    assert 0.0 < e < 1.0


def test_build_scoreboard_column_order():
    rows = [
        {"method": "uncorrected (PCA)", "sil_batch": 0.1, "sil_celltype": 0.1,
         "n_clusters": 8, "ari_celltype": 0.4, "batch_entropy": 0.2},
        {"method": "harmony", "sil_batch": 0.01, "sil_celltype": 0.2,
         "n_clusters": 7, "ari_celltype": 0.9, "batch_entropy": 0.98},
    ]
    board = build_scoreboard(rows)
    assert list(board["method"]) == ["uncorrected (PCA)", "harmony"]
    assert list(board.columns) == [
        "method", "sil_batch", "sil_celltype", "n_clusters",
        "ari_celltype", "batch_entropy"]


# ---------------------------------------------------------------------------
# Synthetic end-to-end pipeline (no network)
# ---------------------------------------------------------------------------

def _synthetic_ifnb(seed=0):
    """A tiny ifnb stand-in: 3 cell types x 2 UNEQUAL conditions + a batch block.

    Returns the loader's 4-tuple (counts, genes, cells, meta) so it can be
    substituted for shanuz.datasets.ifnb. CTRL has 120 cells, STIM 150 — unequal
    on purpose, to drive the RPCA reciprocal-anchor path.
    """
    rng = np.random.default_rng(seed)
    G = 300
    types = ["Mono", "T", "B"]
    sizes = {"CTRL": [45, 40, 35], "STIM": [55, 50, 45]}   # 120 vs 150
    cols, stim, anno, cells = [], [], [], []
    c = 0
    for b in ("CTRL", "STIM"):
        for ti, t in enumerate(types):
            for _ in range(sizes[b][ti]):
                base = rng.gamma(0.3, size=G) + 0.05
                base[ti * 50:(ti + 1) * 50] += 5.0          # cell-type marker block
                if b == "STIM":
                    base[200:250] += 4.0                     # batch (IFN-like) block
                cols.append(rng.poisson(base * 3000.0 / base.sum()))
                stim.append(b)
                anno.append(t)
                cells.append(f"cell{c}")
                c += 1
    counts = sp.csc_matrix(np.asarray(cols).T)               # G x n_cells
    genes = [f"g{i}" for i in range(G)]
    meta = pd.DataFrame({"stim": stim, "seurat_annotations": anno}, index=cells)
    return counts, genes, cells, meta


@pytest.fixture(scope="module")
def integrated():
    """Run the full pipeline once on the synthetic dataset, figures in a tmp dir."""
    pytest.importorskip("harmonypy")
    data = _synthetic_ifnb()
    mp = pytest.MonkeyPatch()
    mp.setattr(tut, "ifnb", lambda data_dir=None: data)
    tmp = Path(tempfile.mkdtemp())
    mp.setattr(tut, "FIGURES", tmp)
    obj, summary = tut.run_full(verbose=False, do_umap=False)
    yield tut, obj, summary, tmp
    mp.undo()


def test_load_builds_object_with_labels():
    data = _synthetic_ifnb()
    mp = pytest.MonkeyPatch()
    mp.setattr(tut, "ifnb", lambda data_dir=None: data)
    obj = tut.load_ifnb_object()
    mp.undo()
    assert tut.BATCH in obj.meta_data.columns
    assert tut.CELLTYPE in obj.meta_data.columns
    assert obj.meta_data[tut.BATCH].nunique() == 2
    assert (obj.meta_data[tut.BATCH] == "CTRL").sum() == 120
    assert (obj.meta_data[tut.BATCH] == "STIM").sum() == 150


def test_prep_writes_hvgs_and_pca(integrated):
    _tut, obj, _summary, tmp = integrated
    assert "pca" in obj.reductions
    # run_full writes the shared variable-feature list for the R reference.
    assert (tmp / "hvg_features.txt").exists()
    hvg = (tmp / "hvg_features.txt").read_text().split()
    assert len(hvg) > 0


def test_pipeline_stores_every_method_partition(integrated):
    _tut, obj, _summary, _tmp = integrated
    for key in ("pca", "harmony", "cca", "rpca"):
        assert f"clusters_{key}" in obj.meta_data.columns
        assert key in obj.reductions
        emb = obj.reductions[key].cell_embeddings
        assert emb.shape[0] == len(obj.cell_names())
        assert np.isfinite(emb).all()


def test_harmony_mixes_batches(integrated):
    _tut, obj, _summary, _tmp = integrated
    batch = np.asarray(obj.meta_data[tut.BATCH])
    pca_emb = obj.reductions["pca"].cell_embeddings
    harm_emb = obj.reductions["harmony"].cell_embeddings
    m_pca = mixing_metrics(pca_emb, batch, batch)
    m_harm = mixing_metrics(harm_emb, batch, batch)
    # Harmony lowers batch separation relative to the uncorrected PCA.
    assert m_harm["sil_batch"] < m_pca["sil_batch"]


def test_summarize_scoreboard_has_a_row_per_method(integrated):
    _tut, obj, summary, _tmp = integrated
    board = summary["scoreboard"]
    assert list(board["method"])[0] == "uncorrected (PCA)"
    assert set(tut.METHODS).issubset(set(board["method"]))
    assert summary["n_batches"] == 2
    assert summary["n_celltypes"] == 3


def test_report_concordance_reads_fabricated_r_calls(integrated):
    _tut, obj, _summary, tmp = integrated
    # Fabricate an R calls file that equals Python's own partitions -> ARI 1.0.
    cells = obj.cell_names()
    meta = obj.meta_data
    df = pd.DataFrame({
        "cell": cells,
        "stim": np.asarray(meta[tut.BATCH]),
        "seurat_annotations": np.asarray(meta[tut.CELLTYPE]),
        "R_pca": np.asarray(meta["clusters_pca"]),
        "R_harmony": np.asarray(meta["clusters_harmony"]),
        "R_cca": np.asarray(meta["clusters_cca"]),
        "R_rpca": np.asarray(meta["clusters_rpca"]),
    })
    df.to_csv(tmp / "r_calls.csv", index=False)
    out = tut.report_concordance(obj, verbose=False)
    assert out is not None
    assert out["harmony"]["ari_partitions"] == pytest.approx(1.0)
    # Python-vs-Python cell-type recovery equals R's on the same labels.
    assert out["harmony"]["py_ari_celltype"] == pytest.approx(
        out["harmony"]["r_ari_celltype"])


def test_report_concordance_absent_returns_none(integrated):
    _tut, obj, _summary, tmp = integrated
    missing = tmp / "does_not_exist.csv"
    assert tut.report_concordance(obj, r_calls_path=missing, verbose=False) is None
