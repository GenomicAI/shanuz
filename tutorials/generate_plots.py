"""Generate all comparison plots for the PBMC 3k tutorial using shanuz.plotting.

Saves PNG figures to tutorials/figures/ for use in the tutorial README.

Usage
-----
    python tutorials/generate_plots.py [--data-dir PATH]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from shanuz.datasets import pbmc3k
from shanuz.shanuz import create_shanuz_object
from shanuz.preprocessing import (
    normalize_data, find_variable_features, scale_data, percentage_feature_set,
)
from shanuz.reduction import run_pca
from shanuz.neighbors import find_neighbors
from shanuz.clustering import find_clusters
from shanuz.umap import run_umap
from shanuz.markers import find_all_markers
from shanuz.plotting import (
    vln_plot, feature_plot, dim_plot, elbow_plot, feature_scatter,
    variable_feature_plot, viz_dim_loadings, dim_heatmap, do_heatmap, ridge_plot,
)

FIGURES = Path(__file__).parent / "figures"
FIGURES.mkdir(exist_ok=True)

CELL_TYPE_MAP = {
    "0": "Naive CD4 T",
    "1": "Memory CD4 T",
    "2": "CD14+ Mono",
    "3": "CD8 T",
    "4": "B",
    "5": "FCGR3A+ Mono",
    "6": "NK",
    "7": "DC",
    "8": "Platelet",
}


def _save(fig, name):
    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    import matplotlib.pyplot as plt
    plt.close(fig)
    print(f"  Saved {path.name}")


def run_pipeline(data_dir=None):
    print("Running pipeline...")
    counts, genes, cells = pbmc3k(data_dir=data_dir)
    pbmc = create_shanuz_object(
        counts=counts, assay="RNA", min_cells=3, min_features=200,
        project="pbmc3k", feature_names=genes, cell_names=cells,
    )
    percentage_feature_set(pbmc, pattern=r"^MT-", col_name="percent.mt")
    md = pbmc.meta_data
    keep = (md["nFeature_RNA"] > 200) & (md["nFeature_RNA"] < 2500) & (md["percent.mt"] < 5)
    pbmc = pbmc.subset(cells=list(md.index[keep]))

    normalize_data(pbmc, normalization_method="LogNormalize", scale_factor=10000)
    find_variable_features(pbmc, selection_method="vst", nfeatures=2000)
    hvg = pbmc.assays["RNA"].variable_features
    scale_data(pbmc, features=pbmc.assays["RNA"]._all_feature_names)
    run_pca(pbmc, n_pcs=50, features=hvg, reduction_name="pca")
    find_neighbors(pbmc, dims=range(10), k_param=20)
    find_clusters(pbmc, resolution=0.5, algorithm=1, random_seed=0)
    run_umap(pbmc, dims=range(10), reduction_name="umap", seed=42)
    all_markers = find_all_markers(pbmc, only_pos=True, min_pct=0.25, logfc_threshold=0.25)
    pbmc.rename_idents(CELL_TYPE_MAP)
    print("Pipeline complete.")
    return pbmc, hvg, all_markers


def main(data_dir=None):
    pbmc, hvg, all_markers = run_pipeline(data_dir)
    print("\nGenerating plots...")

    # 1. QC violin (before filtering — re-create pre-filter object just for QC)
    counts, genes, cells = pbmc3k(data_dir=data_dir)
    pbmc_raw = create_shanuz_object(
        counts=counts, assay="RNA", min_cells=3, min_features=200,
        project="pbmc3k_qc", feature_names=genes, cell_names=cells,
    )
    percentage_feature_set(pbmc_raw, pattern=r"^MT-", col_name="percent.mt")
    # pt_size=1.0 matches R's default of showing individual data points on violins
    _save(vln_plot(pbmc_raw, ["nFeature_RNA", "nCount_RNA", "percent.mt"], ncol=3,
                   figsize=(12, 4), pt_size=1.0), "01_qc_violin.png")

    # 2. QC scatter
    _save(feature_scatter(pbmc_raw, "nCount_RNA", "percent.mt",
                          group_by=None, figsize=(5.5, 4.5)), "02a_qc_scatter_mt.png")
    _save(feature_scatter(pbmc_raw, "nCount_RNA", "nFeature_RNA",
                          group_by=None, figsize=(5.5, 4.5)), "02b_qc_scatter_feat.png")

    # Combined QC scatter (side by side)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    import matplotlib
    fig_a = feature_scatter(pbmc_raw, "nCount_RNA", "percent.mt", figsize=(5, 4))
    fig_b = feature_scatter(pbmc_raw, "nCount_RNA", "nFeature_RNA", figsize=(5, 4))
    plt.close(fig_a); plt.close(fig_b)
    # Re-draw directly into subplots
    from shanuz.plotting import _get_expression, _get_groups, _palette
    for ax_idx, (f1, f2) in enumerate([("nCount_RNA", "percent.mt"),
                                        ("nCount_RNA", "nFeature_RNA")]):
        ax = axes[ax_idx]
        x = _get_expression(pbmc_raw, f1)
        y = _get_expression(pbmc_raw, f2)
        ax.scatter(x, y, s=4, alpha=0.5, color="#F8766D" if ax_idx == 0 else "#00BFC4",
                   linewidths=0)
        ax.set_xlabel(f1); ax.set_ylabel(f2)
        ax.set_title(f"{f1} vs {f2}")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    _save(fig, "02_qc_scatter.png")

    # 3. Variable features
    _save(variable_feature_plot(pbmc, label=True, n_label=10, figsize=(9, 5)),
          "03_variable_features.png")

    # 4. PCA loadings — viz_dim_loadings matches R's VizDimLoadings (bar charts)
    _save(viz_dim_loadings(pbmc, reduction="pca", dims=[1, 2], n_features=15,
                           figsize=(10, 6)), "04_pca_loadings.png")

    # 5. PCA dimplot — no labels to match R's DimPlot(pbmc, reduction="pca")
    _save(dim_plot(pbmc, reduction="pca", label=False,
                   title="PCA — coloured by cluster", figsize=(7, 6)),
          "05_pca_dimplot.png")

    # 6. Elbow plot
    _save(elbow_plot(pbmc, ndims=20, figsize=(7, 4)), "06_elbow_plot.png")

    # 7. UMAP clusters — no labels to match R's DimPlot(pbmc, reduction="umap")
    _save(dim_plot(pbmc, reduction="umap", label=False,
                   title="UMAP — coloured by cluster", figsize=(7, 6)),
          "07_umap_clusters.png")

    # 8. Feature plots (canonical markers)
    canon = ["MS4A1", "CD79A", "NKG7", "GNLY", "FCGR3A", "LYZ", "PPBP", "CD8A", "IL7R"]
    _save(feature_plot(pbmc, canon, reduction="umap", ncol=3,
                       figsize=(13, 11)), "08_feature_plots.png")

    # 9. Violin plots — MS4A1 + CD79A only, to match R's markerplots-1.png reference
    _save(vln_plot(pbmc, ["MS4A1", "CD79A"], figsize=(12, 4), ncol=2),
          "09_marker_violins.png")
    # Also save NKG7 + PF4 separately (matches R's markerplots-2.png)
    _save(vln_plot(pbmc, ["NKG7", "PF4"], layer="counts", figsize=(12, 4), ncol=2),
          "09b_marker_violins_counts.png")

    # 10. DoHeatmap — top 5 markers per cluster
    top_genes = (
        all_markers.groupby("cluster", group_keys=False)
        .apply(lambda x: x.nlargest(5, "avg_log2FC"))
        ["gene"].tolist()
    )
    top_genes = list(dict.fromkeys(top_genes))
    pbmc.rename_idents({v: k for k, v in CELL_TYPE_MAP.items()})  # restore cluster numbers
    _save(do_heatmap(pbmc, top_genes, figsize=(14, max(6, len(top_genes) * 0.32))),
          "10_marker_heatmap.png")
    pbmc.rename_idents(CELL_TYPE_MAP)  # restore cell type names

    # 11. Annotated UMAP (cell types)
    _save(dim_plot(pbmc, reduction="umap", label=True,
                   title="UMAP — Cell Type Annotations", figsize=(9, 7)),
          "11_umap_labeled.png")

    # 12. Ridge plot (bonus)
    _save(ridge_plot(pbmc, ["LYZ", "NKG7", "MS4A1", "CD8A"],
                     figsize=(12, 8), ncol=2), "12_ridge_plot.png")

    print(f"\nAll plots saved to {FIGURES}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate PBMC3k comparison plots")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    main(data_dir=args.data_dir)
