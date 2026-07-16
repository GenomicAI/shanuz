"""Generate all figures for the multimodal CITE-seq tutorial.

Runs cbmc_citeseq_tutorial.run_full() and renders the shanuz-side figures to
tutorials/figures_multimodal/.

Usage
-----
    python tutorials/generate_multimodal_plots.py [--data-dir PATH]
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

from tutorials.cbmc_citeseq_tutorial import run_full
from shanuz.plotting import (
    dim_plot, feature_plot, ridge_plot, feature_scatter, vln_plot,
    _get_expression, _get_embedding, _palette,
)

FIGURES = Path(__file__).parent / "figures_multimodal"
FIGURES.mkdir(exist_ok=True)


def _save(fig, name):
    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    import matplotlib.pyplot as plt
    plt.close(fig)
    print(f"  Saved {path.name}")


def _panel(ax, obj, feature, assay, emb, title, continuous):
    """Draw one feature on the embedding (gray low cells, coloured high)."""
    expr = _get_expression(obj, feature, assay=assay)
    if continuous:                       # ADT CLR — use a robust percentile range
        vmin, vmax = np.percentile(expr, 5), np.percentile(expr, 99)
    else:                                # RNA — zeros stay gray
        vmin, vmax = 0.0, max(np.percentile(expr, 99), 1e-9)
    low = expr <= vmin
    ax.scatter(emb[low, 0], emb[low, 1], c="#D3D3D3", s=3, linewidths=0, rasterized=True)
    hi = ~low
    order = np.argsort(expr[hi])
    sc = ax.scatter(emb[hi][order, 0], emb[hi][order, 1], c=expr[hi][order], s=4,
                    cmap="YlOrRd", vmin=vmin, vmax=vmax, linewidths=0, rasterized=True)
    import matplotlib.pyplot as plt
    plt.colorbar(sc, ax=ax, shrink=0.6, pad=0.01, aspect=25)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def _group_panel(ax, emb, labels, title, legend=False):
    """Draw a categorical embedding — one colour per label, shared across panels.

    dim_plot() renders a whole Figure, so the RNA-vs-WNN comparison hand-rolls
    the panels to put both embeddings on one canvas. Colours are keyed off the
    sorted label set so a cell type keeps its colour in both panels.

    The centroid labels mirror DimPlot(label = TRUE) on the R side — including
    its habit of dropping a label into empty space when a type is split across
    the embedding, which is why the legend is the reliable key.
    """
    groups = sorted(set(labels))
    colours = dict(zip(groups, _palette(len(groups))))
    for g in groups:
        m = labels == g
        ax.scatter(emb[m, 0], emb[m, 1], c=colours[g], s=4, linewidths=0,
                   alpha=0.7, label=g, rasterized=True)
        ax.annotate(g, (emb[m, 0].mean(), emb[m, 1].mean()), fontsize=8,
                    fontweight="bold", ha="center", va="center")
    if legend:
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
                  markerscale=3, fontsize=9)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def main(data_dir=None):
    obj, all_markers, anno = run_full(data_dir=data_dir, verbose=False)
    print("\nGenerating multimodal figures...")
    import matplotlib.pyplot as plt

    # 1-2. RNA UMAP — clusters and protein-based cell types
    _save(dim_plot(obj, reduction="umap", group_by="rna_clusters", label=True,
                   title="CBMC — RNA clusters", figsize=(8, 6.5)),
          "01_rna_umap_clusters.png")
    _save(dim_plot(obj, reduction="umap", group_by="protein_celltype", label=True,
                   title="CBMC — cell types (protein + RNA)", figsize=(8.5, 6.5)),
          "02_rna_umap_celltypes.png")

    # 3. Surface-protein feature plots on the RNA UMAP
    proteins = ["CD3", "CD4", "CD8", "CD19", "CD14", "CD16", "CD56", "CD11c"]
    _save(feature_plot(obj, proteins, reduction="umap", assay="ADT", ncol=4,
                       min_cutoff="q05", max_cutoff="q95", pt_size=1.3,
                       figsize=(15, 7)), "03_adt_featureplots.png")

    # 4. Protein (ADT) vs RNA for the same marker, side by side
    pairs = [("CD19", "CD19"), ("CD3", "CD3E"), ("CD8", "CD8A"), ("CD14", "CD14")]
    emb = _get_embedding(obj, "umap")
    fig, axes = plt.subplots(len(pairs), 2, figsize=(9, 4 * len(pairs)))
    for i, (prot, gene) in enumerate(pairs):
        _panel(axes[i, 0], obj, prot, "ADT", emb, f"Protein: {prot}", continuous=True)
        _panel(axes[i, 1], obj, gene, "RNA", emb, f"RNA: {gene}", continuous=False)
    fig.tight_layout()
    _save(fig, "04_protein_vs_rna.png")

    # 5. ADT ridge plots by cell type
    _save(ridge_plot(obj, ["CD3", "CD19", "CD14", "CD56"], group_by="protein_celltype",
                     assay="ADT", ncol=2, figsize=(12, 9)), "05_adt_ridgeplots.png")

    # 6. ADT feature scatter — protein bivariates separate lineages
    _save(feature_scatter(obj, "CD4", "CD8", assay="ADT", group_by="protein_celltype",
                          figsize=(7, 5.5)), "06_adt_scatter_CD4_CD8.png")
    _save(feature_scatter(obj, "CD19", "CD3", assay="ADT", group_by="protein_celltype",
                          figsize=(7, 5.5)), "07_adt_scatter_CD19_CD3.png")

    # 8. WNN joint clusters on the joint embedding (run_full ran run_wnn already)
    _save(dim_plot(obj, reduction="wnn_umap", group_by="wnn_clusters", label=True,
                   title="CBMC — WNN joint clusters", figsize=(8, 6.5)),
          "08_wnn_umap_clusters.png")

    # 9. The point of WNN: same cells, same labels, RNA-only vs joint embedding
    labels = obj.meta_data["protein_celltype"].to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
    _group_panel(axes[0], _get_embedding(obj, "umap"), labels, "RNA-only UMAP")
    _group_panel(axes[1], _get_embedding(obj, "wnn_umap"), labels, "WNN joint UMAP",
                 legend=True)
    fig.tight_layout()
    _save(fig, "09_wnn_vs_rna_umap.png")

    # 10. The learned per-cell modality weights, read as a metadata "feature".
    # Seurat's WNN vignette plots RNA.weight; ADT.weight is 1 - RNA.weight, and
    # reads more directly as "how much this cell leans on protein".
    _save(vln_plot(obj, "ADT.weight", group_by="protein_celltype", figsize=(9, 5)),
          "10_adt_weight_by_celltype.png")

    print(f"\nAll multimodal figures saved to {FIGURES}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate CBMC CITE-seq figures")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    main(data_dir=args.data_dir)
