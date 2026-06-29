"""Generate all figures for the advanced PBMC 8k subclustering tutorial.

Runs the full workflow from pbmc8k_subclustering_tutorial.run_full() and renders
the shanuz-side figures to tutorials/figures_advanced/.

Usage
-----
    python tutorials/generate_advanced_plots.py [--data-dir PATH]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from tutorials.pbmc8k_subclustering_tutorial import run_full
from shanuz.plotting import (
    vln_plot, feature_plot, dim_plot, elbow_plot, do_heatmap,
)

FIGURES = Path(__file__).parent / "figures_advanced"
FIGURES.mkdir(exist_ok=True)


def _save(fig, name):
    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    import matplotlib.pyplot as plt
    plt.close(fig)
    print(f"  Saved {path.name}")


def _top_genes(markers, n):
    genes = (
        markers.groupby("cluster", group_keys=False)
        .apply(lambda x: x.nlargest(n, "avg_log2FC"))
        ["gene"].tolist()
    )
    return list(dict.fromkeys(genes))


def main(data_dir=None):
    pbmc, sub, all_markers, sub_markers, broad, sub_anno = run_full(
        data_dir=data_dir, verbose=False
    )
    print("\nGenerating advanced-tutorial figures...")

    # ---- QC (pre-filter snapshot) ----
    from tutorials.pbmc8k_subclustering_tutorial import load_object
    raw = load_object(data_dir)
    _save(vln_plot(raw, ["nFeature_RNA", "nCount_RNA", "percent.mt"], ncol=3,
                   figsize=(12, 4), pt_size=0.0), "01_qc_violin.png")

    # ---- Global clustering ----
    _save(elbow_plot(pbmc, ndims=30, figsize=(7, 4)), "02_elbow_plot.png")
    _save(dim_plot(pbmc, reduction="umap", group_by="global_clusters", label=True,
                   title="PBMC 8k — global clusters", figsize=(8, 6.5)),
          "03_umap_global_clusters.png")
    _save(dim_plot(pbmc, reduction="umap", group_by="broad_celltype", label=True,
                   title="PBMC 8k — broad cell types", figsize=(8.5, 6.5)),
          "04_umap_global_celltypes.png")

    # ---- Lineage marker feature plots ----
    lineage_markers = ["CD3D", "CD8A", "IL7R", "MS4A1", "LYZ", "FCGR3A",
                       "GNLY", "FCER1A", "PPBP"]
    _save(feature_plot(pbmc, lineage_markers, reduction="umap", ncol=3,
                       pt_size=1.2, figsize=(12, 10)), "05_lineage_featureplots.png")

    # ---- Global marker heatmap (top 5 per cluster) ----
    top_global = _top_genes(all_markers, 5)
    _save(do_heatmap(pbmc, top_global, group_by="global_clusters",
                     figsize=(14, max(6, len(top_global) * 0.22))),
          "06_global_markers_heatmap.png")

    # ---- T/NK subclustering highlight ----
    _save(dim_plot(sub, reduction="umap", group_by="sub_clusters", label=True,
                   title="T/NK compartment — subclusters", figsize=(8, 6.5)),
          "07_umap_tnk_subclusters.png")
    _save(dim_plot(sub, reduction="umap", group_by="tnk_subset", label=True,
                   title="T/NK compartment — annotated subsets", figsize=(8.5, 6.5)),
          "08_umap_tnk_subsets.png")

    # ---- T/NK subset markers ----
    tnk_markers = ["CCR7", "SELL", "IL7R", "S100A4", "CD8A", "GZMK", "GNLY", "NKG7"]
    _save(feature_plot(sub, tnk_markers, reduction="umap", ncol=4,
                       pt_size=1.5, figsize=(15, 7)), "09_tnk_subset_featureplots.png")
    _save(vln_plot(sub, ["CCR7", "S100A4", "CD8A", "GNLY"], group_by="tnk_subset",
                   ncol=2, figsize=(11, 7)), "10_tnk_subset_violins.png")

    top_sub = _top_genes(sub_markers, 6)
    _save(do_heatmap(sub, top_sub, group_by="sub_clusters",
                     figsize=(12, max(6, len(top_sub) * 0.22))),
          "11_tnk_markers_heatmap.png")

    print(f"\nAll advanced figures saved to {FIGURES}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate PBMC8k advanced figures")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    main(data_dir=args.data_dir)
