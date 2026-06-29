"""PBMC 3k Tutorial — Shanuz Python implementation.

Mirrors the official Seurat PBMC 3k guided clustering tutorial step-by-step:
  https://satijalab.org/seurat/articles/pbmc3k_tutorial

Each step prints summary statistics comparable to the R tutorial so you can
validate the Python results directly.

Usage
-----
    python tutorials/pbmc3k_tutorial.py [--data-dir PATH]

If --data-dir is not supplied the PBMC 3k dataset is downloaded automatically
to ~/.shanuz_data/pbmc3k (~24 MB, 10X Genomics).

References
----------
Hao et al. (2024) Nature Biotechnology — https://doi.org/10.1038/s41587-023-01767-y
Stuart et al. (2019) Cell — https://doi.org/10.1016/j.cell.2019.05.031
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Ensure the package root is on the path when running the script directly
# ---------------------------------------------------------------------------
_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from shanuz.datasets import pbmc3k
from shanuz.shanuz import create_shanuz_object
from shanuz.preprocessing import (
    normalize_data,
    find_variable_features,
    scale_data,
    percentage_feature_set,
)
from shanuz.reduction import run_pca
from shanuz.neighbors import find_neighbors
from shanuz.clustering import find_clusters
from shanuz.umap import run_umap
from shanuz.markers import find_markers, find_all_markers


# ---------------------------------------------------------------------------
# Known expected results from the R tutorial (used for validation)
# ---------------------------------------------------------------------------
EXPECTED = {
    # After CreateSeuratObject(min.cells=3, min.features=200)
    "n_features_raw": 13714,
    "n_cells_raw": 2700,
    # After QC filter (nFeature_RNA 200-2500, percent.mt < 5)
    "n_cells_filtered": 2638,
    # Top 10 HVGs from the R tutorial
    "top10_hvg": {
        "PPBP", "LYZ", "S100A9", "IGLL5", "GNLY",
        "FTL", "PF4", "FTH1", "GNG11", "S100A8",
    },
    # Canonical marker genes expected per cell type
    "canonical_markers": {
        "CD14+ Mono":    ["LYZ", "CD14", "S100A9"],
        "NK":            ["NKG7", "GNLY"],
        "B":             ["MS4A1", "CD79A"],
        "CD8 T":         ["CD8A"],
        "DC":            ["FCER1A"],
        "Platelet":      ["PPBP"],
    },
    # Number of clusters expected (resolution=0.5 → ~9)
    "n_clusters_expected": 9,
}


def validate(label: str, value, expected=None, atol: float = 0.05) -> None:
    """Print a validation line. Green check if matches, red ✗ otherwise."""
    if expected is None:
        print(f"  [INFO] {label}: {value}")
        return
    if isinstance(expected, set) and isinstance(value, (set, list)):
        overlap = set(value) & expected
        pct = len(overlap) / len(expected) * 100
        ok = pct >= 50
        mark = "OK" if ok else "FAIL"
        print(f"  [{mark}] {label}: overlap {len(overlap)}/{len(expected)} "
              f"({pct:.0f}%)  got={sorted(set(value))[:10]}")
    elif isinstance(expected, int):
        ok = abs(value - expected) <= max(1, int(expected * atol))
        mark = "OK" if ok else "FAIL"
        print(f"  [{mark}] {label}: {value}  (expected ~{expected})")
    else:
        print(f"  [INFO] {label}: {value}  (expected {expected})")


def section(title: str) -> None:
    print(f"\n{'=' * 60}")
    print(f"  {title}")
    print(f"{'=' * 60}")


# ---------------------------------------------------------------------------
# Main tutorial
# ---------------------------------------------------------------------------

def run_tutorial(data_dir: str | None = None) -> None:
    t0_total = time.time()

    # -----------------------------------------------------------------------
    section("1. Load Data")
    # -----------------------------------------------------------------------
    t0 = time.time()
    counts, genes, cells = pbmc3k(data_dir=data_dir)
    print(f"  Raw matrix: {counts.shape[0]} genes × {counts.shape[1]} cells  "
          f"({time.time() - t0:.1f}s)")

    # -----------------------------------------------------------------------
    section("2. Create Shanuz Object (min.cells=3, min.features=200)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    pbmc = create_shanuz_object(
        counts=counts,
        assay="RNA",
        min_cells=3,
        min_features=200,
        project="pbmc3k",
        feature_names=genes,
        cell_names=cells,
    )
    n_feat = len(pbmc.assays["RNA"]._all_feature_names)
    n_cells = len(pbmc.cell_names())
    print(f"  {n_feat} features × {n_cells} cells  ({time.time() - t0:.1f}s)")
    validate("n_features after min.cells=3", n_feat, EXPECTED["n_features_raw"])
    validate("n_cells after min.features=200", n_cells, EXPECTED["n_cells_raw"])

    # -----------------------------------------------------------------------
    section("3. QC Metrics")
    # -----------------------------------------------------------------------
    t0 = time.time()
    percentage_feature_set(pbmc, pattern=r"^MT-", col_name="percent.mt")
    md = pbmc.meta_data
    print(f"  nFeature_RNA: mean={md['nFeature_RNA'].mean():.0f}  "
          f"min={md['nFeature_RNA'].min():.0f}  "
          f"max={md['nFeature_RNA'].max():.0f}")
    print(f"  percent.mt:   mean={md['percent.mt'].mean():.2f}%  "
          f"max={md['percent.mt'].max():.2f}%  ({time.time() - t0:.1f}s)")

    # -----------------------------------------------------------------------
    section("4. Filter Cells (nFeature 200-2500, percent.mt < 5)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    keep = (
        (md["nFeature_RNA"] > 200) &
        (md["nFeature_RNA"] < 2500) &
        (md["percent.mt"] < 5)
    )
    keep_cells = list(md.index[keep])
    pbmc = pbmc.subset(cells=keep_cells)
    n_cells_filt = len(pbmc.cell_names())
    print(f"  {n_cells_filt} cells retained  ({time.time() - t0:.1f}s)")
    validate("n_cells after QC filter", n_cells_filt, EXPECTED["n_cells_filtered"])

    # -----------------------------------------------------------------------
    section("5. Normalize Data (LogNormalize, scale.factor=10000)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    normalize_data(pbmc, normalization_method="LogNormalize", scale_factor=10000)
    print(f"  Log-normalization complete  ({time.time() - t0:.1f}s)")

    # Quick sanity check: mean of normalized data should be ~1-3 log units
    rna = pbmc.assays["RNA"]
    norm = rna.layers["data"]
    mean_expr = float(np.array(norm.mean()))
    print(f"  Mean log-normalized expression: {mean_expr:.4f}")

    # -----------------------------------------------------------------------
    section("6. Find Variable Features (VST, nfeatures=2000)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    find_variable_features(pbmc, selection_method="vst", nfeatures=2000)
    hvg = pbmc.assays["RNA"].variable_features
    top10 = hvg[:10]
    print(f"  {len(hvg)} variable features selected  ({time.time() - t0:.1f}s)")
    print(f"  Top 10 HVGs: {top10}")
    validate("Top 10 HVG overlap with R tutorial", set(top10), EXPECTED["top10_hvg"])

    # -----------------------------------------------------------------------
    section("7. Scale Data")
    # -----------------------------------------------------------------------
    t0 = time.time()
    all_genes = pbmc.assays["RNA"]._all_feature_names
    scale_data(pbmc, features=all_genes)
    print(f"  Scaled {len(all_genes)} genes  ({time.time() - t0:.1f}s)")

    # -----------------------------------------------------------------------
    section("8. Run PCA (npc=50)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    run_pca(pbmc, n_pcs=50, features=hvg, reduction_name="pca")
    pca_emb = pbmc.reductions["pca"].cell_embeddings
    print(f"  PCA: {pca_emb.shape[0]} cells × {pca_emb.shape[1]} PCs  "
          f"({time.time() - t0:.1f}s)")
    stdev = pbmc.reductions["pca"].stdev
    print(f"  PC1 stdev={stdev[0]:.3f}  PC2 stdev={stdev[1]:.3f}  "
          f"PC10 stdev={stdev[9]:.3f}")
    # Top loadings of PC1 (should be ribosomal/mitochondrial or strong cell-type genes)
    loadings = pbmc.reductions["pca"].feature_loadings
    feat_names = pbmc.reductions["pca"]._feature_names
    top_pc1 = [feat_names[i] for i in np.argsort(np.abs(loadings[:, 0]))[::-1][:5]]
    print(f"  Top PC1 loading genes: {top_pc1}")

    # -----------------------------------------------------------------------
    section("9. Find Neighbors (dims=1:10, k=20)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    find_neighbors(pbmc, dims=range(10), k_param=20)
    print(f"  KNN+SNN graphs built  ({time.time() - t0:.1f}s)")
    print(f"  Graphs: {list(pbmc.graphs)}")

    # -----------------------------------------------------------------------
    section("10. Find Clusters (resolution=0.5)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    find_clusters(pbmc, resolution=0.5, algorithm=1, random_seed=0)
    cluster_counts = pbmc.meta_data["seurat_clusters"].value_counts().sort_index()
    n_clusters = len(cluster_counts)
    print(f"  {n_clusters} clusters found  ({time.time() - t0:.1f}s)")
    print("  Cells per cluster:")
    for c, n in cluster_counts.items():
        print(f"    Cluster {c}: {n} cells")
    validate("Number of clusters", n_clusters, EXPECTED["n_clusters_expected"])

    # -----------------------------------------------------------------------
    section("11. Run UMAP (dims=1:10)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    run_umap(pbmc, dims=range(10), reduction_name="umap", seed=42)
    umap_emb = pbmc.reductions["umap"].cell_embeddings
    print(f"  UMAP: {umap_emb.shape}  ({time.time() - t0:.1f}s)")
    print(f"  UMAP range: x=[{umap_emb[:,0].min():.2f},{umap_emb[:,0].max():.2f}]  "
          f"y=[{umap_emb[:,1].min():.2f},{umap_emb[:,1].max():.2f}]")

    # -----------------------------------------------------------------------
    section("12. Find Cluster 2 Markers")
    # -----------------------------------------------------------------------
    t0 = time.time()
    c2_markers = find_markers(pbmc, ident_1="2", only_pos=True)
    print(f"  Cluster 2 markers: {len(c2_markers)} genes  ({time.time() - t0:.1f}s)")
    if len(c2_markers) > 0:
        top5 = c2_markers.head(5)
        print("  Top 5 cluster 2 markers:")
        print(top5[["avg_log2FC", "pct.1", "pct.2", "p_val_adj"]].to_string())

    # -----------------------------------------------------------------------
    section("13. Find All Markers (only.pos=True, logfc.threshold=0.25)")
    # -----------------------------------------------------------------------
    t0 = time.time()
    all_markers = find_all_markers(
        pbmc, only_pos=True, min_pct=0.25, logfc_threshold=0.25
    )
    print(f"  Total marker genes: {len(all_markers)}  ({time.time() - t0:.1f}s)")
    top_per_cluster = (
        all_markers.groupby("cluster")
        .apply(lambda x: x.nsmallest(3, "p_val"))
        .reset_index(drop=True)
    )
    print("\n  Top 3 markers per cluster:")
    for cluster in sorted(all_markers["cluster"].unique(), key=lambda x: int(x)):
        sub = top_per_cluster[top_per_cluster["cluster"] == cluster]
        genes_str = ", ".join(sub["gene"].tolist())
        print(f"    Cluster {cluster}: {genes_str}")

    # -----------------------------------------------------------------------
    section("14. Validate Canonical Marker Expression")
    # -----------------------------------------------------------------------
    marker_genes_flat = [g for gs in EXPECTED["canonical_markers"].values() for g in gs]
    found_in_markers = set(all_markers["gene"].tolist())
    for cell_type, canon_genes in EXPECTED["canonical_markers"].items():
        hits = [g for g in canon_genes if g in found_in_markers]
        pct = len(hits) / len(canon_genes) * 100
        mark = "OK" if pct >= 50 else "WARN"
        print(f"  [{mark}] {cell_type}: {hits} ({pct:.0f}% canonical markers found)")

    # -----------------------------------------------------------------------
    section("15. Cell Type Annotation (R tutorial mapping)")
    # -----------------------------------------------------------------------
    # Map cluster labels to cell types using known marker patterns
    # This is R's RenameIdents() step
    cluster_to_celltype = _assign_cell_types(all_markers, pbmc)
    names = new_cluster_ids = {str(k): v for k, v in cluster_to_celltype.items()}
    pbmc.rename_idents(names)

    celltype_counts = pd.Series(list(pbmc.idents)).value_counts()
    print("\n  Cell type distribution:")
    for ct, n in celltype_counts.sort_values(ascending=False).items():
        print(f"    {ct}: {n} cells")

    # -----------------------------------------------------------------------
    section("Summary")
    # -----------------------------------------------------------------------
    total = time.time() - t0_total
    print(f"\n  Total runtime: {total:.1f}s")
    print(f"\n{pbmc}")

    return pbmc


def _assign_cell_types(
    all_markers: pd.DataFrame,
    pbmc,
) -> dict[int, str]:
    """Heuristically assign cell types to clusters based on top markers.

    Mirrors the manual annotation step in the R tutorial.
    """
    # Canonical markers for each cell type (ordered by specificity)
    markers_ref = {
        "Naive CD4 T":   ["IL7R", "CCR7"],
        "CD14+ Mono":    ["CD14", "LYZ"],
        "Memory CD4 T":  ["IL7R", "S100A4"],
        "B":             ["MS4A1"],
        "CD8 T":         ["CD8A"],
        "FCGR3A+ Mono":  ["FCGR3A", "MS4A7"],
        "NK":            ["GNLY", "NKG7"],
        "DC":            ["FCER1A", "CST3"],
        "Platelet":      ["PPBP"],
    }

    clusters = sorted(set(str(i) for i in pbmc.idents))
    cluster_top_genes: dict[str, set] = {}
    for cluster in clusters:
        sub = all_markers[all_markers["cluster"] == cluster].head(50)
        cluster_top_genes[cluster] = set(sub["gene"].tolist())

    assignment: dict[str, str] = {}
    used_types: set[str] = set()

    for cluster in clusters:
        top_genes = cluster_top_genes.get(cluster, set())
        best_type = "Unknown"
        best_score = 0
        for cell_type, canon in markers_ref.items():
            if cell_type in used_types:
                continue
            score = sum(1 for g in canon if g in top_genes)
            if score > best_score:
                best_score = score
                best_type = cell_type
        if best_score > 0:
            used_types.add(best_type)
        assignment[cluster] = best_type

    return assignment


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PBMC 3k tutorial")
    parser.add_argument(
        "--data-dir",
        default=None,
        help="Directory for PBMC3k data (default: ~/.shanuz_data/pbmc3k)",
    )
    args = parser.parse_args()
    run_tutorial(data_dir=args.data_dir)
