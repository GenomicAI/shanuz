"""Tests for the v5 feature additions:

  * SCTransform              (sctransform.py)
  * AddModuleScore           (module_score.py)
  * CellCycleScoring         (module_score.py)
  * DotPlot                  (plotting.py)
  * extra DE tests LR/negbinom/roc (markers.py)

These use small synthetic data with injected structure — network-free.
"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz.shanuz import create_shanuz_object  # noqa: E402
from shanuz.preprocessing import normalize_data  # noqa: E402
from shanuz.markers import find_markers  # noqa: E402
from shanuz.module_score import add_module_score, cell_cycle_scoring, CC_GENES  # noqa: E402
from shanuz.sctransform import sctransform  # noqa: E402
from shanuz.reduction import run_pca  # noqa: E402


# ---------------------------------------------------------------------------
# SCTransform
# ---------------------------------------------------------------------------

def test_sctransform_creates_sct_assay_and_stabilizes_variance():
    rng = np.random.default_rng(0)
    G, N = 300, 400
    depth = rng.integers(500, 5000, size=N)
    base = rng.gamma(0.3, size=(G, 1))
    counts = rng.poisson(base * depth[None, :] / depth.mean()).astype(float)
    counts[:10, :200] += rng.poisson(8, size=(10, 200))  # structured genes
    obj = create_shanuz_object(
        counts=sp.csc_matrix(counts), assay="RNA",
        feature_names=[f"g{i}" for i in range(G)],
        cell_names=[f"c{i}" for i in range(N)],
    )
    sctransform(obj, n_cells=400, n_features=50, seed=0)

    assert "SCT" in obj.assay_names()
    assert obj.active_assay == "SCT"
    sct = obj.assays["SCT"]
    assert set(["counts", "data", "scale.data"]).issubset(sct.layers)

    # scale.data holds residuals for the variable features only (subset layer).
    assert sct.layers["scale.data"].shape == (50, N)
    assert len(sct.variable_features) == 50

    # Pearson residuals are variance-stabilized: mean ~0, std ~O(1).
    sd = sct.layers["scale.data"].toarray()
    assert abs(sd.mean()) < 0.2
    assert 0.5 < sd.std() < 3.0

    # Structured genes dominate the residual-variance ranking.
    vf = set(sct.variable_features)
    assert sum(f"g{i}" in vf for i in range(10)) >= 8

    # PCA runs on the SCT assay.
    run_pca(obj, n_pcs=10)
    assert obj.reductions["pca"].cell_embeddings.shape == (N, 10)


def test_sctransform_drops_low_detection_genes():
    rng = np.random.default_rng(1)
    counts = rng.poisson(0.5, size=(60, 80)).astype(float)
    counts[0, :] = 0.0          # never detected
    counts[1, :2] = 1.0         # detected in 2 cells only
    counts[1, 2:] = 0.0
    obj = create_shanuz_object(
        counts=sp.csc_matrix(counts), assay="RNA",
        feature_names=[f"g{i}" for i in range(60)],
        cell_names=[f"c{i}" for i in range(80)],
    )
    sctransform(obj, n_cells=80, n_features=20, min_cells=5, seed=0)
    sct_genes = set(obj.assays["SCT"]._all_feature_names)
    assert "g0" not in sct_genes      # 0 cells
    assert "g1" not in sct_genes      # < min_cells


# ---------------------------------------------------------------------------
# AddModuleScore
# ---------------------------------------------------------------------------

def test_add_module_score_higher_in_program_cells():
    rng = np.random.default_rng(0)
    G, N = 200, 100
    base = rng.poisson(0.5, size=(G, N)).astype(float)
    prog = [f"g{i}" for i in range(10)]
    for gi in range(10):                     # program: high only in cells 0-49
        base[gi, :50] += 12
    for gi in range(10, 80):                 # decoys: high uniformly (control bins)
        base[gi, :] += rng.poisson(6, size=N)
    obj = create_shanuz_object(
        counts=sp.csc_matrix(base), assay="RNA",
        feature_names=[f"g{i}" for i in range(G)],
        cell_names=[f"c{i}" for i in range(N)],
    )
    normalize_data(obj)
    add_module_score(obj, prog, name="Prog", seed=1)
    assert "Prog1" in obj.meta_data.columns
    score = obj.meta_data["Prog1"].values
    assert score[:50].mean() > score[50:].mean()


def test_add_module_score_dict_uses_keys_as_columns():
    rng = np.random.default_rng(2)
    obj = create_shanuz_object(
        counts=sp.csc_matrix(rng.poisson(1.0, size=(50, 40)).astype(float)),
        assay="RNA", feature_names=[f"g{i}" for i in range(50)],
        cell_names=[f"c{i}" for i in range(40)],
    )
    normalize_data(obj)
    add_module_score(obj, {"A": ["g0", "g1"], "B": ["g2", "g3"]}, seed=1)
    assert "A" in obj.meta_data.columns and "B" in obj.meta_data.columns


# ---------------------------------------------------------------------------
# CellCycleScoring
# ---------------------------------------------------------------------------

def test_cell_cycle_scoring_assigns_phases():
    rng = np.random.default_rng(0)
    s_genes = CC_GENES["s_genes"][:20]
    g2m_genes = CC_GENES["g2m_genes"][:20]
    genes = list(dict.fromkeys([f"g{i}" for i in range(60)] + s_genes + g2m_genes))
    gi = {g: i for i, g in enumerate(genes)}
    m = rng.poisson(0.4, size=(len(genes), 90)).astype(float)
    for g in s_genes:
        m[gi[g], :30] += 10        # S cells
    for g in g2m_genes:
        m[gi[g], 30:60] += 10      # G2M cells
    obj = create_shanuz_object(
        counts=sp.csc_matrix(m), assay="RNA",
        feature_names=genes, cell_names=[f"c{i}" for i in range(90)],
    )
    normalize_data(obj)
    cell_cycle_scoring(obj)

    assert {"S.Score", "G2M.Score", "Phase"}.issubset(obj.meta_data.columns)
    phase = obj.meta_data["Phase"].values
    # The injected S / G2M populations are recovered as the majority phase.
    assert (phase[:30] == "S").sum() >= 25
    assert (phase[30:60] == "G2M").sum() >= 25


# ---------------------------------------------------------------------------
# DotPlot
# ---------------------------------------------------------------------------

def test_dot_plot_returns_figure():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from shanuz.plotting import dot_plot

    rng = np.random.default_rng(0)
    genes = ["LYZ", "CD3D", "MS4A1"] + [f"g{i}" for i in range(10)]
    gi = {g: i for i, g in enumerate(genes)}
    base = rng.poisson(0.4, size=(len(genes), 90)).astype(float)
    for cells, g in zip([range(0, 30), range(30, 60), range(60, 90)],
                        ["LYZ", "CD3D", "MS4A1"]):
        for c in cells:
            base[gi[g], c] += 10
    obj = create_shanuz_object(
        counts=sp.csc_matrix(base), assay="RNA",
        feature_names=genes, cell_names=[f"c{i}" for i in range(90)],
    )
    normalize_data(obj)
    obj.idents = ["A"] * 30 + ["B"] * 30 + ["C"] * 30

    fig = dot_plot(obj, ["LYZ", "CD3D", "MS4A1"])
    assert fig is not None and len(fig.axes) >= 1
    import matplotlib.pyplot as plt
    plt.close(fig)


# ---------------------------------------------------------------------------
# Extra DE tests
# ---------------------------------------------------------------------------

def _two_group_object():
    rng = np.random.default_rng(0)
    a = rng.poisson(0.5, size=(40, 40))
    a[0, :] += 8                       # g0 marks group A
    b = rng.poisson(0.5, size=(40, 40))
    b[1, :] += 8                       # g1 marks group B
    obj = create_shanuz_object(
        counts=sp.csc_matrix(np.hstack([a, b]).astype(float)), assay="RNA",
        feature_names=[f"g{i}" for i in range(40)],
        cell_names=[f"c{i}" for i in range(80)],
    )
    normalize_data(obj)
    obj.idents = ["A"] * 40 + ["B"] * 40
    return obj


@pytest.mark.parametrize("test_use", ["LR", "negbinom"])
def test_de_pvalue_tests_rank_marker_top(test_use):
    obj = _two_group_object()
    res = find_markers(obj, ident_1="A", test_use=test_use,
                       min_pct=0.0, logfc_threshold=0.0)
    assert list(res.columns) == ["p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]
    assert res.index[0] == "g0"
    assert res.loc["g0", "p_val"] < 0.01


def test_roc_test_returns_auc_power_and_ranks_marker():
    obj = _two_group_object()
    res = find_markers(obj, ident_1="A", test_use="roc",
                       min_pct=0.0, logfc_threshold=0.0)
    assert set(["myAUC", "avg_diff", "power", "avg_log2FC", "pct.1", "pct.2"]) == set(res.columns)
    assert "p_val" not in res.columns
    # g0 perfectly separates the groups → top by power, AUC ~1.
    assert res.index[0] == "g0"
    assert res.loc["g0", "power"] > 0.9
    assert res.loc["g0", "myAUC"] > 0.95


def test_lr_test_accepts_latent_vars():
    obj = _two_group_object()
    rng = np.random.default_rng(3)
    obj.meta_data["cc"] = rng.standard_normal(len(obj))
    res = find_markers(obj, ident_1="A", test_use="LR", latent_vars=["cc"],
                       min_pct=0.0, logfc_threshold=0.0)
    assert res.index[0] == "g0"
