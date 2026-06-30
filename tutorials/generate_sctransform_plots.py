"""Generate all figures for the SCTransform (PBMC 3k) tutorial.

Runs pbmc3k_sctransform_tutorial.run_full() and renders the shanuz-side figures
to tutorials/figures_sctransform/.

Usage
-----
    python tutorials/generate_sctransform_plots.py [--data-dir PATH]
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

from tutorials.pbmc3k_sctransform_tutorial import (
    run_full, VIGNETTE_MARKERS_1, VIGNETTE_MARKERS_2, VLN_MARKERS,
)
from shanuz.plotting import dim_plot, feature_plot, vln_plot

FIGURES = Path(__file__).parent / "figures_sctransform"
FIGURES.mkdir(exist_ok=True)


def _save(fig, name):
    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    import matplotlib.pyplot as plt
    plt.close(fig)
    print(f"  Saved {path.name}")


def main(data_dir=None):
    sct, std, sct_anno, std_anno, sct_markers = run_full(data_dir=data_dir, verbose=False)
    print("\nGenerating SCTransform figures...")
    import matplotlib.pyplot as plt

    # 1-2. SCT UMAP — clusters and annotated cell types
    sct.idents = sct.meta_data["sct_clusters"].astype(str).tolist()
    _save(dim_plot(sct, reduction="umap", group_by="sct_clusters", label=True,
                   title="PBMC 3k — SCTransform clusters", figsize=(8, 6.5)),
          "01_sct_umap_clusters.png")
    _save(dim_plot(sct, reduction="umap", group_by="sct_celltype", label=True,
                   title="PBMC 3k — SCTransform cell types", figsize=(8.5, 6.5)),
          "02_sct_umap_celltypes.png")

    # 3-4. Vignette FeaturePlots (protein-marker panels on the SCT UMAP)
    _save(feature_plot(sct, VIGNETTE_MARKERS_1, reduction="umap", assay="SCT",
                       ncol=3, min_cutoff="q05", max_cutoff="q95", pt_size=1.2,
                       figsize=(13, 7)), "03_sct_featureplots_1.png")
    _save(feature_plot(sct, VIGNETTE_MARKERS_2, reduction="umap", assay="SCT",
                       ncol=3, min_cutoff="q05", max_cutoff="q95", pt_size=1.2,
                       figsize=(13, 7)), "04_sct_featureplots_2.png")

    # 5. Violin plots of the eight vignette markers (by cluster).
    #    pt_size>0 overlays jittered cells to match R's VlnPlot(pt.size = 0.2)
    #    speckled look (matplotlib sizes points by area, so the value differs).
    _save(vln_plot(sct, VLN_MARKERS, group_by="sct_clusters", assay="SCT",
                   ncol=4, pt_size=2.0, figsize=(15, 7)), "05_sct_violins.png")

    # 6. SCTransform vs standard log-normalization, side by side
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
    for ax, obj, key, title in (
        (axes[0], sct, "sct_clusters", f"SCTransform — {sct.meta_data['sct_clusters'].nunique()} clusters (dims 1:30)"),
        (axes[1], std, "std_clusters", f"LogNormalize — {std.meta_data['std_clusters'].nunique()} clusters (dims 1:10)"),
    ):
        sub = dim_plot(obj, reduction="umap", group_by=key, label=True, title=title)
        # Re-render onto the shared figure by copying the single-axes plot.
        src = sub.axes[0]
        for coll in src.collections:
            off = coll.get_offsets()
            ax.scatter(off[:, 0], off[:, 1], c=coll.get_facecolors(), s=6, linewidths=0)
        for txt in src.texts:
            ax.text(txt.get_position()[0], txt.get_position()[1], txt.get_text(),
                    fontsize=9, fontweight="bold", ha="center", va="center")
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xticks([])
        ax.set_yticks([])
        plt.close(sub)
    fig.tight_layout()
    _save(fig, "06_sct_vs_std_umap.png")

    print(f"\nAll SCTransform figures saved to {FIGURES}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate PBMC 3k SCTransform figures")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    main(data_dir=args.data_dir)
