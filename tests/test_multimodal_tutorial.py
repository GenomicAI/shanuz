"""Tests for the multimodal CITE-seq tutorial's cell annotation and WNN section.

Builds a tiny synthetic two-assay (RNA + ADT) object with cleanly-separated
signal and checks the combined protein-priority / RNA-fallback gating in
annotate_cells(), the run_wnn() joint-clustering flow, and the figures the
walkthrough's Step 8 embeds. Network-free.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz.shanuz import create_shanuz_object  # noqa: E402
from shanuz.assay5 import create_assay5_object  # noqa: E402
from shanuz.preprocessing import (  # noqa: E402
    normalize_data, find_variable_features, scale_data,
)
from shanuz.reduction import run_pca  # noqa: E402
from shanuz.plotting import dim_plot, vln_plot  # noqa: E402
from tutorials.cbmc_citeseq_tutorial import annotate_cells, run_wnn  # noqa: E402
from tutorials.generate_multimodal_plots import _group_panel  # noqa: E402


def _two_assay_object(rna_levels, adt_levels, n=6):
    clusters = list(rna_levels)
    rna_genes = sorted({g for lv in rna_levels.values() for g in lv} | {"MALAT1"})
    adt_prots = sorted({p for lv in adt_levels.values() for p in lv})
    rg = {g: i for i, g in enumerate(rna_genes)}
    ag = {p: i for i, p in enumerate(adt_prots)}

    cells, idents = [], []
    rmat = np.zeros((len(rna_genes), len(clusters) * n))
    amat = np.zeros((len(adt_prots), len(clusters) * n))
    col = 0
    for cl in clusters:
        for _ in range(n):
            rmat[rg["MALAT1"], col] = 50
            for g, v in rna_levels[cl].items():
                rmat[rg[g], col] = v
            for p, v in adt_levels[cl].items():
                amat[ag[p], col] = v
            cells.append(f"cell_{col}")
            idents.append(cl)
            col += 1

    obj = create_shanuz_object(
        counts=sp.csc_matrix(rmat), assay="RNA",
        feature_names=rna_genes, cell_names=cells,
    )
    normalize_data(obj)
    obj.assays["ADT"] = create_assay5_object(
        counts=sp.csc_matrix(amat), feature_names=adt_prots,
        cell_names=cells, key="adt_",
    )
    normalize_data(obj, assay="ADT", normalization_method="CLR", margin=2)
    obj.idents = idents
    return obj


def test_annotate_cells_protein_then_rna_fallback():
    obj = _two_assay_object(
        rna_levels={
            "t": {}, "b": {}, "nk": {}, "mono": {},
            "plt": {"PPBP": 200, "PF4": 200},     # protein-panel-less → RNA
        },
        adt_levels={
            "t":   {"CD3": 200, "CD4": 200},
            "b":   {"CD19": 200},
            "nk":  {"CD16": 200, "CD56": 200},
            "mono": {"CD14": 200},
            "plt": {},
        },
    )
    anno = annotate_cells(obj)
    assert anno["t"] == "CD4 T"
    assert anno["b"] == "B"
    assert anno["nk"] == "NK"
    assert anno["mono"] == "CD14+ Mono"
    assert anno["plt"] == "Platelet"      # resolved by RNA, not protein


def test_annotate_cells_cd8_split():
    # A non-T cluster ("b") is needed so CD3 (margin=2 CLR) is actually
    # discriminative — otherwise a uniformly-high protein centers to ~0.
    obj = _two_assay_object(
        rna_levels={"cd4": {}, "cd8": {}, "b": {}},
        adt_levels={
            "cd4": {"CD3": 200, "CD4": 200, "CD8": 1},
            "cd8": {"CD3": 200, "CD8": 200, "CD4": 1},
            "b":   {"CD19": 200},
        },
    )
    anno = annotate_cells(obj)
    assert anno["cd4"] == "CD4 T"
    assert anno["cd8"] == "CD8 T"
    assert anno["b"] == "B"


# ---------------------------------------------------------------------------
# WNN section (run_wnn)
# ---------------------------------------------------------------------------

def _wnn_ready_object(seed=0, per=30):
    """RNA (workflow run through PCA) + CLR-normalised ADT, ready for run_wnn."""
    rng = np.random.default_rng(seed)
    groups = ("A", "B", "C")
    n = len(groups) * per
    Grna, Padt = 120, 30
    rna = rng.gamma(0.3, size=(Grna, n)) + 0.05
    adt = rng.gamma(0.3, size=(Padt, n)) + 0.05
    cells = []
    for ci in range(n):
        g = groups[ci // per]
        if g == "A":
            rna[0:30, ci] += 6.0      # RNA separates A from {B,C}
        if g == "C":
            adt[0:10, ci] += 6.0      # ADT separates C from {A,B}
        cells.append(f"cell{ci}")

    obj = create_shanuz_object(
        counts=sp.csc_matrix(rng.poisson(rna).astype(float)), assay="RNA",
        feature_names=[f"g{i}" for i in range(Grna)], cell_names=cells,
    )
    normalize_data(obj)
    find_variable_features(obj, selection_method="vst", nfeatures=100)
    scale_data(obj)
    run_pca(obj, n_pcs=15)

    obj.assays["ADT"] = create_assay5_object(
        counts=sp.csc_matrix(rng.poisson(adt).astype(float)),
        feature_names=[f"prot{i}" for i in range(Padt)], cell_names=cells, key="adt_",
    )
    normalize_data(obj, assay="ADT", normalization_method="CLR", margin=2)
    return obj, n


def test_run_wnn_builds_joint_graphs_umap_and_weights():
    obj, n = _wnn_ready_object()
    run_wnn(obj, rna_dims=range(10), resolution=1.0)

    # ADT reduction + joint graphs + joint UMAP produced.
    assert "apca" in obj.reductions
    assert "wknn" in obj.graphs and "wsnn" in obj.graphs
    assert obj.reductions["wnn_umap"].cell_embeddings.shape == (n, 2)
    assert np.isfinite(obj.reductions["wnn_umap"].cell_embeddings).all()

    # Joint clustering column + per-cell modality weights summing to 1.
    assert "wnn_clusters" in obj.meta_data.columns
    w = obj.meta_data[["RNA.weight", "ADT.weight"]].to_numpy()
    assert np.all(w >= 0) and np.all(w <= 1)
    assert np.allclose(w.sum(axis=1), 1.0, atol=1e-6)


# ---------------------------------------------------------------------------
# Step 8 figures (generate_multimodal_plots)
# ---------------------------------------------------------------------------

def test_wnn_umap_plots_on_the_joint_embedding():
    """Figure 08 — dim_plot has to accept the joint reduction, not just "umap"."""
    obj, _ = _wnn_ready_object()
    run_wnn(obj, rna_dims=range(10), resolution=1.0)
    fig = dim_plot(obj, reduction="wnn_umap", group_by="wnn_clusters", label=True)
    assert fig.axes and fig.axes[0].collections


def test_modality_weight_plots_as_a_metadata_feature():
    """Figure 10 — ADT.weight is a metadata column, not a gene.

    Seurat's VlnPlot resolves features against metadata; shanuz's vln_plot does
    the same via _get_expression, and that fallback is what the figure rides on.
    """
    obj, n = _wnn_ready_object()
    run_wnn(obj, rna_dims=range(10), resolution=1.0)
    obj.meta_data["ct"] = ["x" if i < n // 2 else "y" for i in range(n)]
    fig = vln_plot(obj, "ADT.weight", group_by="ct")
    assert fig.axes
    # The violin must span the real weight spread, not collapse to a point.
    lo, hi = fig.axes[0].get_ylim()
    assert hi > lo


def test_group_panel_colours_labels_consistently():
    """Figure 09 draws two embeddings; a cell type must keep its colour."""
    import matplotlib.pyplot as plt

    emb = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
    labels = np.array(["B", "A", "A", "B"])
    fig, axes = plt.subplots(1, 2)
    _group_panel(axes[0], emb, labels, "left")
    _group_panel(axes[1], emb[::-1], labels[::-1], "right", legend=True)

    # One collection per group, in sorted-label order, matching colours across
    # panels even though the second panel sees the cells in reverse.
    assert len(axes[0].collections) == len(axes[1].collections) == 2
    for left, right in zip(axes[0].collections, axes[1].collections):
        assert np.allclose(left.get_facecolor(), right.get_facecolor())
    assert axes[0].get_title() == "left"

    # The legend is the reliable key when a centroid label lands off-cluster.
    assert axes[0].get_legend() is None
    assert [t.get_text() for t in axes[1].get_legend().get_texts()] == ["A", "B"]
    plt.close(fig)
