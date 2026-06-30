"""SCTransform Tutorial — PBMC 3k with Shanuz.

A Python port of Seurat's sctransform vignette
(https://satijalab.org/seurat/articles/sctransform_vignette) on the PBMC 3k
dataset. SCTransform replaces the NormalizeData -> FindVariableFeatures ->
ScaleData trio with a single regularized negative-binomial model, returning
Pearson residuals that more effectively remove technical (sequencing-depth and
percent-mito) effects. The vignette's headline result is that this sharper
normalization, run over more PCs (dims 1:30), resolves finer immune subsets
than the standard log-normalization workflow.

This script runs BOTH workflows on the same cells so their cluster resolution
can be compared directly:
  * SCT workflow  : sctransform(vars.to.regress="percent.mt") -> PCA -> dims 1:30
  * Std workflow  : LogNormalize -> VST -> ScaleData -> PCA -> dims 1:10

Usage
-----
    python tutorials/pbmc3k_sctransform_tutorial.py [--data-dir PATH]

The PBMC 3k dataset (~24 MB) downloads automatically to ~/.shanuz_data/pbmc3k.

References
----------
Hafemeister C, Satija R (2019) Genome Biology 20, 296.
Choudhary S, Satija R (2022) Genome Biology 23, 27. (sctransform v2)
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from shanuz.datasets import pbmc3k
from shanuz.shanuz import create_shanuz_object
from shanuz.preprocessing import (
    normalize_data, find_variable_features, scale_data, percentage_feature_set,
)
from shanuz.sctransform import sctransform
from shanuz.reduction import run_pca
from shanuz.neighbors import find_neighbors
from shanuz.clustering import find_clusters
from shanuz.umap import run_umap
from shanuz.markers import find_all_markers
from shanuz.plotting import _get_expression


# Markers shown in the vignette's FeaturePlots — these define the fine subsets
# SCTransform is meant to resolve.
VIGNETTE_MARKERS_1 = ["CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7"]
VIGNETTE_MARKERS_2 = ["CD3D", "ISG15", "TCL1A", "FCER2", "XCL1", "FCGR3A"]
VLN_MARKERS = ["CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"]

# Fine-grained lineage panel for relative-enrichment annotation. Several of
# these (CD8 effector vs naive, CD4 naive vs memory, NK bright vs dim) are the
# distinctions SCTransform is meant to sharpen.
FINE_MARKERS = {
    "Naive CD4 T":   ["IL7R", "CCR7", "LEF1", "SELL"],
    "Memory CD4 T":  ["IL7R", "S100A4", "IL32", "ANXA1"],
    "CD8 Naive/Mem": ["CD8A", "CD8B", "CCR7"],
    "CD8 Effector":  ["CD8A", "GZMK", "CCL5", "NKG7"],
    "B":             ["MS4A1", "CD79A", "TCL1A", "FCER2"],
    "CD14+ Mono":    ["CD14", "LYZ", "S100A8", "S100A9"],
    "FCGR3A+ Mono":  ["FCGR3A", "MS4A7"],
    "NK":            ["GNLY", "NKG7", "KLRD1", "XCL1"],
    "DC":            ["FCER1A", "CST3"],
    "pDC":           ["SERPINF1", "ITM2C"],
    "Platelet":      ["PPBP", "PF4"],
}

# Which fine lineages count as resolved T-cell subsets, for the comparison.
T_SUBSETS = {"Naive CD4 T", "Memory CD4 T", "CD8 Naive/Mem", "CD8 Effector"}


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def load_object(data_dir=None):
    counts, genes, cells = pbmc3k(data_dir=data_dir)
    pbmc = create_shanuz_object(
        counts=counts, assay="RNA", min_cells=3, min_features=200,
        project="pbmc3k", feature_names=genes, cell_names=cells,
    )
    percentage_feature_set(pbmc, pattern=r"^MT-", col_name="percent.mt")
    return pbmc


def run_sct_workflow(pbmc, dims=range(30), resolution=0.8, seed=0):
    """SCTransform workflow (mirrors the vignette: regress percent.mt, dims 1:30)."""
    sctransform(pbmc, vars_to_regress=["percent.mt"], n_features=3000, seed=42)
    run_pca(pbmc, n_pcs=50, features=pbmc.assays["SCT"].variable_features)
    find_neighbors(pbmc, dims=dims, k_param=20)
    find_clusters(pbmc, resolution=resolution, algorithm=1, random_seed=seed)
    run_umap(pbmc, dims=dims, seed=42)
    return pbmc


def run_std_workflow(pbmc, dims=range(10), resolution=0.8, seed=0):
    """Standard log-normalization workflow for comparison (dims 1:10)."""
    normalize_data(pbmc, normalization_method="LogNormalize", scale_factor=10000)
    find_variable_features(pbmc, selection_method="vst", nfeatures=2000)
    scale_data(pbmc, features=pbmc.assays["RNA"].variable_features)
    run_pca(pbmc, n_pcs=50, features=pbmc.assays["RNA"].variable_features)
    find_neighbors(pbmc, dims=dims, k_param=20)
    find_clusters(pbmc, resolution=resolution, algorithm=1, random_seed=seed)
    run_umap(pbmc, dims=dims, seed=42)
    return pbmc


def annotate_clusters(pbmc, marker_sets, assay=None):
    """Assign each cluster to the lineage whose markers are most *enriched*.

    Each marker's per-cluster mean is z-scored across clusters, so a cluster
    scores on a lineage by relative enrichment (matching the advanced tutorial).
    Reuse is allowed: several clusters may share a lineage.
    """
    idents = np.array([str(i) for i in pbmc.idents])
    clusters = sorted(set(idents), key=lambda x: int(x) if x.isdigit() else x)
    feats = set(pbmc.assays[assay or pbmc.active_assay]._all_feature_names)

    needed = {g for gs in marker_sets.values() for g in gs if g in feats}
    zmean = {}
    for g in needed:
        expr = _get_expression(pbmc, g, assay=assay)
        per_cluster = np.array([expr[idents == c].mean() for c in clusters])
        sd = per_cluster.std()
        zmean[g] = (per_cluster - per_cluster.mean()) / sd if sd > 1e-9 \
            else np.zeros(len(clusters))

    assignment = {}
    for ci, c in enumerate(clusters):
        best, best_score = "Unknown", -np.inf
        for lineage, genes in marker_sets.items():
            present = [g for g in genes if g in zmean]
            if not present:
                continue
            score = float(np.mean([zmean[g][ci] for g in present]))
            if score > best_score:
                best_score, best = score, lineage
        assignment[c] = best
    return assignment


# ---------------------------------------------------------------------------
# Full run
# ---------------------------------------------------------------------------

def section(title):
    print(f"\n{'=' * 64}\n  {title}\n{'=' * 64}")


def _n_t_subsets(anno):
    return len({lab for lab in anno.values() if lab in T_SUBSETS})


def run_full(data_dir=None, verbose=True):
    t0 = time.time()

    if verbose:
        section("1. Load PBMC 3k + percent.mt")
    base = load_object(data_dir)
    if verbose:
        print(f"  {len(base)} cells x {len(base.assays['RNA']._all_feature_names)} genes")

    # ---- SCTransform workflow ----
    if verbose:
        section("2. SCTransform workflow (regress percent.mt, PCA dims 1:30)")
    sct = load_object(data_dir)
    run_sct_workflow(sct)
    n_sct = sct.meta_data["seurat_clusters"].nunique()
    sct_assay = sct.assays["SCT"]
    if verbose:
        print(f"  SCT assay: {len(sct_assay._all_feature_names)} genes, "
              f"{len(sct_assay.variable_features)} variable features")
        print(f"  {n_sct} clusters at resolution 0.8 (dims 1:30)")
        print(f"  Top SCT variable features: {sct_assay.variable_features[:10]}")

    sct_anno = annotate_clusters(sct, FINE_MARKERS, assay="SCT")
    sct.meta_data["sct_clusters"] = [str(i) for i in sct.idents]
    sct.stash_ident("sct_clusters")
    sct.rename_idents(sct_anno)
    sct.meta_data["sct_celltype"] = [str(i) for i in sct.idents]

    # ---- Standard workflow ----
    if verbose:
        section("3. Standard LogNormalize workflow (PCA dims 1:10)")
    std = load_object(data_dir)
    run_std_workflow(std)
    n_std = std.meta_data["seurat_clusters"].nunique()
    if verbose:
        print(f"  {n_std} clusters at resolution 0.8 (dims 1:10)")
    std_anno = annotate_clusters(std, FINE_MARKERS, assay="RNA")
    std.meta_data["std_clusters"] = [str(i) for i in std.idents]
    std.stash_ident("std_clusters")
    std.rename_idents(std_anno)

    # ---- Comparison ----
    if verbose:
        section("4. SCTransform vs standard — cluster resolution")
        print(f"  {'workflow':<16}{'clusters':>10}{'T-subsets resolved':>22}")
        print(f"  {'-' * 46}")
        print(f"  {'SCTransform':<16}{n_sct:>10}{_n_t_subsets(sct_anno):>22}")
        print(f"  {'LogNormalize':<16}{n_std:>10}{_n_t_subsets(std_anno):>22}")
        print("\n  SCT cluster -> annotation:")
        for c, lab in sct_anno.items():
            print(f"    cluster {c:>2} -> {lab}")

    if verbose:
        section("5. SCT marker check (vignette FeaturePlot genes)")
    sct.idents = sct.meta_data["sct_clusters"].astype(str).tolist()
    sct_markers = find_all_markers(sct, only_pos=True, min_pct=0.25,
                                   logfc_threshold=0.25)
    sct.rename_idents(sct_anno)
    if verbose:
        for clid in sorted(sct_markers["cluster"].unique(), key=int):
            top = sct_markers[sct_markers["cluster"] == clid].nsmallest(4, "p_val")
            print(f"    cluster {clid}: " + ", ".join(top["gene"].tolist()))

    if verbose:
        section("Summary")
        print(f"  Total runtime: {time.time() - t0:.1f}s")
        print(f"\n  SCT:  {sct}")

    return sct, std, sct_anno, std_anno, sct_markers


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PBMC 3k SCTransform tutorial")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    run_full(data_dir=args.data_dir)
