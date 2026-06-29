"""Advanced PBMC 8k Tutorial — Clustering + Subclustering with Shanuz.

A more complex companion to the PBMC 3k tutorial. It reproduces the standard
Seurat guided-clustering workflow (Satija et al. 2015; Butler et al. 2018) on a
larger 10x Genomics dataset (~8,400 PBMCs, GRCh38) and then adds the advanced
*subclustering* step used throughout the Seurat reference papers to resolve
fine-grained immune subsets: the T/NK lymphoid compartment is isolated and
re-analysed from scratch (HVG -> PCA -> neighbours -> clusters -> UMAP) to
separate naive CD4, memory CD4, CD8, and NK populations that the global
clustering lumps together.

Usage
-----
    python tutorials/pbmc8k_subclustering_tutorial.py [--data-dir PATH]

The PBMC 8k dataset (~38 MB) downloads automatically to ~/.shanuz_data/pbmc8k.

References
----------
Satija R, Farrell JA, Gennert D, et al. (2015) Nature Biotechnology 33, 495-502.
Butler A, Hoffman P, Smibert P, et al. (2018) Nature Biotechnology 36, 411-420.
Hao Y, Hao S, Andersen-Nissen E, et al. (2021) Cell 184, 3573-3587.
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

from shanuz.datasets import pbmc8k
from shanuz.shanuz import create_shanuz_object
from shanuz.preprocessing import (
    normalize_data, find_variable_features, scale_data, percentage_feature_set,
)
from shanuz.reduction import run_pca
from shanuz.neighbors import find_neighbors
from shanuz.clustering import find_clusters
from shanuz.umap import run_umap
from shanuz.markers import find_all_markers
from shanuz.plotting import _get_expression


# ---------------------------------------------------------------------------
# Canonical marker panels
# ---------------------------------------------------------------------------

# Broad lineages — each cluster is assigned to the lineage whose markers it
# expresses most strongly (reuse allowed: several clusters can be "CD4 T").
BROAD_MARKERS = {
    "CD4 T":         ["IL7R", "CD3D", "CCR7"],
    "CD8 T":         ["CD8A", "CD8B", "GZMK"],
    "NK":            ["GNLY", "NKG7", "KLRD1"],
    "B":             ["MS4A1", "CD79A", "CD79B"],
    "CD14+ Mono":    ["CD14", "LYZ", "S100A8"],
    "FCGR3A+ Mono":  ["FCGR3A", "MS4A7"],
    "DC":            ["FCER1A", "CST3"],
    "Platelet":      ["PPBP", "PF4"],
}

# The lymphoid lineages we re-analyse (subcluster) together.
LYMPHOID_LINEAGES = {"CD4 T", "CD8 T", "NK"}

# Genes used by the hierarchical T/NK subset annotator.
TNK_PANEL = ["CD3D", "CD3E", "CD8A", "CD8B", "GNLY", "NKG7", "KLRD1",
             "CCR7", "SELL", "LEF1", "IL7R", "S100A4", "GZMK"]


# ---------------------------------------------------------------------------
# Core pipeline (shared with generate_advanced_plots.py)
# ---------------------------------------------------------------------------

def run_pipeline(pbmc, dims=range(10), resolution=0.5, n_pcs=50, k_param=20,
                 nfeatures=2000, normalize=True):
    """Run the standard workflow on a (possibly already-subset) Shanuz object."""
    if normalize:
        normalize_data(pbmc, normalization_method="LogNormalize", scale_factor=10000)
    find_variable_features(pbmc, selection_method="vst", nfeatures=nfeatures)
    scale_data(pbmc, features=pbmc.assays["RNA"]._all_feature_names)
    run_pca(pbmc, n_pcs=n_pcs, features=pbmc.assays["RNA"].variable_features,
            reduction_name="pca")
    find_neighbors(pbmc, dims=dims, k_param=k_param)
    find_clusters(pbmc, resolution=resolution, algorithm=1, random_seed=0)
    run_umap(pbmc, dims=dims, reduction_name="umap", seed=42)
    return pbmc


def load_object(data_dir=None):
    counts, genes, cells = pbmc8k(data_dir=data_dir)
    pbmc = create_shanuz_object(
        counts=counts, assay="RNA", min_cells=3, min_features=200,
        project="pbmc8k", feature_names=genes, cell_names=cells,
    )
    percentage_feature_set(pbmc, pattern=r"^MT-", col_name="percent.mt")
    return pbmc


def qc_filter(pbmc, max_features=2500, max_mt=5.0):
    md = pbmc.meta_data
    keep = (
        (md["nFeature_RNA"] > 200) &
        (md["nFeature_RNA"] < max_features) &
        (md["percent.mt"] < max_mt)
    )
    return pbmc.subset(cells=list(md.index[keep]))


def annotate_clusters(pbmc, marker_sets):
    """Assign each cluster to the lineage whose markers are most *enriched*.

    Each marker's per-cluster mean expression is z-scored across clusters, so a
    cluster scores on a lineage by how relatively enriched its markers are
    (CD8A is lower-magnitude than IL7R but still flags the CD8 cluster). Each
    lineage's score is the mean z-score of its present markers; argmax wins.

    Returns {cluster_label: lineage_name}. Reuse is allowed so several clusters
    can map to the same lineage (e.g. multiple T-cell clusters).
    """
    idents = np.array([str(i) for i in pbmc.idents])
    clusters = sorted(set(idents), key=lambda x: int(x) if x.isdigit() else x)
    feats = set(pbmc.assays["RNA"]._all_feature_names)

    needed = {g for gs in marker_sets.values() for g in gs if g in feats}
    # z-score of each marker's per-cluster mean across clusters
    zmean = {}
    for g in needed:
        expr = _get_expression(pbmc, g)
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


def annotate_tnk_subsets(sub):
    """Annotate T/NK subclusters with biologically-ordered gating.

    Flat argmax over marker panels fails here because the high-magnitude naive
    markers (CCR7/SELL) outscore the lower-magnitude but definitive CD8 markers.
    We instead gate in lineage-priority order on mean expression:

      1. NK            — CD3 low and NKG7/GNLY high (NK cells are CD3-negative)
      2. CD8 T         — CD8A/B detectable, or a CD3+ cytotoxic (NKG7/GZMK)
                         program (captures CD8 effector plus MAIT/gamma-delta T,
                         which are CD8-lineage cytotoxic cells, never CD4)
      3. CD4 Naive     — naive markers (CCR7/SELL/LEF1) high
      4. CD4 Memory    — otherwise (IL7R / S100A4)
    """
    idents = np.array([str(i) for i in sub.idents])
    clusters = sorted(set(idents), key=lambda x: int(x) if x.isdigit() else x)
    feats = set(sub.assays["RNA"]._all_feature_names)
    expr = {g: _get_expression(sub, g) for g in TNK_PANEL if g in feats}

    def m(genes, mask):
        vals = [expr[g][mask].mean() for g in genes if g in expr]
        return float(np.mean(vals)) if vals else 0.0

    assignment = {}
    for c in clusters:
        mask = idents == c
        cd3 = m(["CD3D", "CD3E"], mask)
        cd8 = m(["CD8A", "CD8B"], mask)
        nk = m(["GNLY", "NKG7"], mask)
        cyto = m(["NKG7", "GZMK"], mask)
        naive = m(["CCR7", "SELL", "LEF1"], mask)
        if cd3 < 0.75 and nk > 1.5:
            assignment[c] = "NK"
        elif cd8 > 0.6 or (cd3 >= 0.75 and cyto > 1.5):
            assignment[c] = "CD8 T"
        elif naive > 0.9:
            assignment[c] = "CD4 Naive"
        else:
            assignment[c] = "CD4 Memory"
    return assignment


# ---------------------------------------------------------------------------
# Pretty-printing helpers
# ---------------------------------------------------------------------------

def section(title):
    print(f"\n{'=' * 64}\n  {title}\n{'=' * 64}")


def top_markers_table(all_markers, n=3):
    out = []
    for cl in sorted(all_markers["cluster"].unique(), key=lambda x: int(x)):
        sub = all_markers[all_markers["cluster"] == cl].nsmallest(n, "p_val")
        out.append(f"    Cluster {cl}: " + ", ".join(sub["gene"].tolist()))
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Full run
# ---------------------------------------------------------------------------

def run_full(data_dir=None, verbose=True):
    """Execute the entire advanced workflow and return all artefacts."""
    t0 = time.time()

    if verbose:
        section("1. Load PBMC 8k + QC")
    pbmc = load_object(data_dir)
    n_raw = len(pbmc)
    pbmc = qc_filter(pbmc)
    if verbose:
        print(f"  {n_raw} cells -> {len(pbmc)} after QC "
              f"(nFeature 200-2500, percent.mt < 5)")

    if verbose:
        section("2. Standard workflow (normalize -> HVG -> PCA -> cluster -> UMAP)")
    run_pipeline(pbmc, dims=range(10), resolution=0.5)
    n_clusters = pbmc.meta_data["seurat_clusters"].nunique()
    if verbose:
        print(f"  {n_clusters} global clusters at resolution 0.5")
        print(f"  Top HVGs: {pbmc.assays['RNA'].variable_features[:8]}")

    if verbose:
        section("3. Marker genes per global cluster")
    all_markers = find_all_markers(pbmc, only_pos=True, min_pct=0.25,
                                   logfc_threshold=0.25)
    if verbose:
        print(top_markers_table(all_markers, n=3))

    if verbose:
        section("4. Broad lineage annotation")
    broad = annotate_clusters(pbmc, BROAD_MARKERS)
    pbmc.meta_data["broad_cluster"] = [str(i) for i in pbmc.idents]
    pbmc.stash_ident("global_clusters")
    pbmc.rename_idents(broad)
    pbmc.meta_data["broad_celltype"] = [str(i) for i in pbmc.idents]
    if verbose:
        for c, lin in broad.items():
            print(f"    cluster {c:>2} -> {lin}")
        dist = pd.Series(list(pbmc.idents)).value_counts()
        print("\n  Lineage sizes:")
        for ct, k in dist.items():
            print(f"    {ct}: {k}")

    # ----- Subclustering the T/NK lymphoid compartment -----
    if verbose:
        section("5. Subcluster the T/NK lymphoid compartment")
    lymphoid_clusters = [c for c, lin in broad.items() if lin in LYMPHOID_LINEAGES]
    pbmc.meta_data["global_clusters"]  # ensure present
    global_idents = pbmc.meta_data["global_clusters"].astype(str).values
    cells = pbmc.cell_names()
    lymphoid_cells = [c for c, g in zip(cells, global_idents) if g in set(lymphoid_clusters)]
    if verbose:
        print(f"  Global clusters {sorted(lymphoid_clusters, key=int)} "
              f"-> {len(lymphoid_cells)} T/NK cells")

    sub = pbmc.subset(cells=lymphoid_cells)
    # Re-analyse from counts; data layer is already normalised, so skip renorm.
    run_pipeline(sub, dims=range(10), resolution=0.6, n_pcs=30, normalize=False)
    n_sub = sub.meta_data["seurat_clusters"].nunique()
    if verbose:
        print(f"  {n_sub} subclusters at resolution 0.6")

    if verbose:
        section("6. Annotate T/NK subclusters")
    sub_markers = find_all_markers(sub, only_pos=True, min_pct=0.25,
                                   logfc_threshold=0.25)
    sub_anno = annotate_tnk_subsets(sub)
    sub.stash_ident("sub_clusters")
    sub.rename_idents(sub_anno)
    sub.meta_data["tnk_subset"] = [str(i) for i in sub.idents]
    if verbose:
        for c, lin in sub_anno.items():
            print(f"    subcluster {c:>2} -> {lin}")
        dist = pd.Series(list(sub.idents)).value_counts()
        print("\n  Subset sizes:")
        for ct, k in dist.items():
            print(f"    {ct}: {k}")
        print("\n  Subcluster top markers:")
        print(top_markers_table(sub_markers, n=4))

    if verbose:
        section("Summary")
        print(f"  Total runtime: {time.time() - t0:.1f}s")
        print(f"\n  Global:  {pbmc}")
        print(f"\n  T/NK:    {sub}")

    return pbmc, sub, all_markers, sub_markers, broad, sub_anno


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PBMC 8k subclustering tutorial")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    run_full(data_dir=args.data_dir)
