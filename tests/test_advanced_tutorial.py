"""Tests for the advanced PBMC 8k subclustering tutorial's annotation helpers.

These exercise the cluster-annotation logic with small, cleanly-separated
synthetic data (no network / dataset download required). The magnitudes are
deliberately extreme so the assertions are robust to reasonable threshold
changes — they verify the *logic*, not exact cut-offs.
"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tutorials.pbmc8k_subclustering_tutorial import (  # noqa: E402
    annotate_clusters, annotate_tnk_subsets, BROAD_MARKERS, TNK_PANEL,
)


def _make_object(gene_levels, n_per_cluster=6, housekeeping=50):
    """Build a normalised Shanuz object from per-cluster gene expression.

    gene_levels: dict[cluster_label -> dict[gene -> raw_count]]. Every gene
    mentioned anywhere becomes a row; unmentioned genes are 0 in that cluster.
    A housekeeping gene ('MALAT1') is added to every cell so totals are sane.
    """
    from shanuz.shanuz import create_shanuz_object
    from shanuz.preprocessing import normalize_data

    clusters = list(gene_levels)
    genes = sorted({g for lv in gene_levels.values() for g in lv} | {"MALAT1"})
    gidx = {g: i for i, g in enumerate(genes)}

    cells, idents, cols = [], [], []
    mat = np.zeros((len(genes), len(clusters) * n_per_cluster))
    col = 0
    for cl in clusters:
        for _ in range(n_per_cluster):
            mat[gidx["MALAT1"], col] = housekeeping
            for g, v in gene_levels[cl].items():
                mat[gidx[g], col] = v
            cells.append(f"cell_{col}")
            idents.append(cl)
            col += 1

    obj = create_shanuz_object(
        counts=sp.csc_matrix(mat), assay="RNA",
        feature_names=genes, cell_names=cells,
    )
    normalize_data(obj)
    obj.idents = idents
    return obj


def test_annotate_clusters_assigns_by_relative_enrichment():
    obj = _make_object({
        "A": {"CD8A": 200, "CD8B": 200},          # CD8 T
        "B": {"MS4A1": 200, "CD79A": 200},        # B
        "C": {"LYZ": 200, "CD14": 200},           # CD14+ Mono
    })
    anno = annotate_clusters(obj, BROAD_MARKERS)
    assert anno["A"] == "CD8 T"
    assert anno["B"] == "B"
    assert anno["C"] == "CD14+ Mono"


def test_annotate_tnk_subsets_hierarchical_gating():
    obj = _make_object({
        # NK: CD3-negative, NKG7/GNLY high
        "nk":    {"NKG7": 200, "GNLY": 200, "KLRD1": 200},
        # CD8: CD3+ and CD8B high (definitive)
        "cd8":   {"CD3D": 30, "CD3E": 30, "CD8A": 200, "CD8B": 200, "GZMK": 60},
        # CD4 naive: CD3+, CCR7/SELL/LEF1 high, no CD8, no cytotoxic
        "naive": {"CD3D": 30, "CD3E": 30, "CCR7": 200, "SELL": 200, "LEF1": 200},
        # CD4 memory: CD3+, IL7R/S100A4 high, no naive markers, no CD8
        "mem":   {"CD3D": 30, "CD3E": 30, "IL7R": 200, "S100A4": 200},
    })
    anno = annotate_tnk_subsets(obj)
    assert anno["nk"] == "NK"
    assert anno["cd8"] == "CD8 T"
    assert anno["naive"] == "CD4 Naive"
    assert anno["mem"] == "CD4 Memory"


def test_annotate_tnk_cytotoxic_cd3_cell_is_cd8_not_cd4():
    """A CD3+ cytotoxic cluster lacking CD8B (MAIT/gamma-delta-like) must group
    with CD8/cytotoxic T, never CD4 — the bug the cyto gate fixes."""
    obj = _make_object({
        "cyto":  {"CD3D": 30, "CD3E": 30, "NKG7": 200, "GZMK": 200},  # CD8B absent
        "memref": {"CD3D": 30, "CD3E": 30, "IL7R": 200, "S100A4": 200},
    })
    anno = annotate_tnk_subsets(obj)
    assert anno["cyto"] == "CD8 T"


def test_tnk_panel_genes_used_by_annotator():
    # Guard: every gene the hierarchical gates read is declared in TNK_PANEL.
    for g in ["CD3D", "CD3E", "CD8A", "CD8B", "GNLY", "NKG7", "CCR7",
              "SELL", "LEF1", "IL7R", "S100A4", "GZMK"]:
        assert g in TNK_PANEL
