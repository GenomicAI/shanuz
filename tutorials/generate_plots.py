"""Generate all comparison plots for the PBMC 3k tutorial.

Saves PNG figures to tutorials/figures/ for use in the tutorial README.

Usage
-----
    python tutorials/generate_plots.py [--data-dir PATH]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

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

FIGURES = Path(__file__).parent / "figures"
FIGURES.mkdir(exist_ok=True)

# Cluster colour palette (matches Seurat default)
CLUSTER_COLORS = [
    "#F8766D", "#CD9600", "#7CAE00", "#00BE67",
    "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC", "#999999",
]

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

# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

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

    all_genes = pbmc.assays["RNA"]._all_feature_names
    scale_data(pbmc, features=all_genes)
    run_pca(pbmc, n_pcs=50, features=hvg, reduction_name="pca")
    find_neighbors(pbmc, dims=range(10), k_param=20)
    find_clusters(pbmc, resolution=0.5, algorithm=1, random_seed=0)
    run_umap(pbmc, dims=range(10), reduction_name="umap", seed=42)

    all_markers = find_all_markers(pbmc, only_pos=True, min_pct=0.25, logfc_threshold=0.25)
    pbmc.rename_idents(CELL_TYPE_MAP)
    print("Pipeline complete.")
    return pbmc, hvg, all_markers


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _save(fig, name):
    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path.name}")


def _cluster_colors(n):
    return CLUSTER_COLORS[:n]


# ---------------------------------------------------------------------------
# 1. QC Violin
# ---------------------------------------------------------------------------

def plot_qc_violin(pbmc):
    md = pbmc.meta_data
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    cols = ["nFeature_RNA", "nCount_RNA", "percent.mt"]
    titles = ["nFeature_RNA", "nCount_RNA", "percent.mt"]
    colors = ["#00BFC4", "#00BFC4", "#00BFC4"]
    for ax, col, title, color in zip(axes, cols, titles, colors):
        data = md[col].values
        parts = ax.violinplot(data, positions=[1], showmedians=True, showextrema=True)
        for pc in parts["bodies"]:
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
        parts["cmedians"].set_color("black")
        ax.set_title(title, fontsize=13)
        ax.set_xticks([])
        ax.set_ylabel(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("QC Metrics", fontsize=14, fontweight="bold")
    fig.tight_layout()
    _save(fig, "01_qc_violin.png")


# ---------------------------------------------------------------------------
# 2. QC Scatter
# ---------------------------------------------------------------------------

def plot_qc_scatter(pbmc):
    md = pbmc.meta_data
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    axes[0].scatter(md["nCount_RNA"], md["percent.mt"], s=4, alpha=0.5, color="#F8766D")
    axes[0].set_xlabel("nCount_RNA"); axes[0].set_ylabel("percent.mt")
    axes[0].set_title("nCount_RNA vs percent.mt")
    axes[0].spines["top"].set_visible(False); axes[0].spines["right"].set_visible(False)

    axes[1].scatter(md["nCount_RNA"], md["nFeature_RNA"], s=4, alpha=0.5, color="#00BFC4")
    axes[1].set_xlabel("nCount_RNA"); axes[1].set_ylabel("nFeature_RNA")
    axes[1].set_title("nCount_RNA vs nFeature_RNA")
    axes[1].spines["top"].set_visible(False); axes[1].spines["right"].set_visible(False)

    fig.tight_layout()
    _save(fig, "02_qc_scatter.png")


# ---------------------------------------------------------------------------
# 3. Variable Features Plot
# ---------------------------------------------------------------------------

def _get_layer(rna):
    d = rna.layers.get("data")
    return d if d is not None else rna.layers.get("counts")


def plot_hvg(pbmc):
    import scipy.sparse as sp
    rna = pbmc.assays["RNA"]
    data = _get_layer(rna)
    if sp.issparse(data):
        means = np.array(data.mean(axis=1)).flatten()
        mean_sq = np.array(data.power(2).mean(axis=1)).flatten()
    else:
        d = np.asarray(data, dtype=float)
        means = d.mean(axis=1); mean_sq = (d ** 2).mean(axis=1)
    variances = mean_sq - means ** 2

    hvg_set = set(rna.variable_features)
    feat_names = rna._all_feature_names
    is_hvg = np.array([f in hvg_set for f in feat_names])

    top10 = rna.variable_features[:10]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.scatter(means[~is_hvg], variances[~is_hvg], s=3, alpha=0.3, color="#AAAAAA",
               label="Other genes")
    ax.scatter(means[is_hvg], variances[is_hvg], s=4, alpha=0.6, color="#F8766D",
               label="Variable features")

    for gene in top10:
        if gene in feat_names:
            idx = feat_names.index(gene)
            ax.annotate(gene, (means[idx], variances[idx]), fontsize=7,
                        xytext=(4, 4), textcoords="offset points")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Average Expression (log scale)")
    ax.set_ylabel("Dispersion (log scale)")
    ax.set_title("Highly Variable Features  (top 2,000 shown in red)", fontsize=13)
    ax.legend(markerscale=3, fontsize=9)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    _save(fig, "03_variable_features.png")


# ---------------------------------------------------------------------------
# 4. PCA Loadings (PC1 + PC2)
# ---------------------------------------------------------------------------

def plot_pca_loadings(pbmc):
    dr = pbmc.reductions["pca"]
    loadings = dr.feature_loadings
    feat = dr._feature_names

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for i, ax in enumerate(axes):
        col = loadings[:, i]
        top_pos = np.argsort(col)[::-1][:10]
        top_neg = np.argsort(col)[:10]
        idx = np.concatenate([top_pos[::-1], top_neg])
        genes = [feat[j] for j in idx]
        vals = col[idx]
        colors = ["#F8766D" if v > 0 else "#00BFC4" for v in vals]
        y = np.arange(len(genes))
        ax.barh(y, vals, color=colors, edgecolor="none")
        ax.set_yticks(y); ax.set_yticklabels(genes, fontsize=8)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_title(f"PC{i+1} top loadings", fontsize=12)
        ax.set_xlabel("Loading score")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    _save(fig, "04_pca_loadings.png")


# ---------------------------------------------------------------------------
# 5. PCA DimPlot (PC1 vs PC2)
# ---------------------------------------------------------------------------

def plot_pca_dimplot(pbmc):
    emb = pbmc.reductions["pca"].cell_embeddings
    clusters = pbmc.meta_data["seurat_clusters"].astype(str).values
    unique = sorted(set(clusters), key=int)
    colors = _cluster_colors(len(unique))
    cmap = dict(zip(unique, colors))

    fig, ax = plt.subplots(figsize=(7, 6))
    for cl in unique:
        mask = clusters == cl
        ax.scatter(emb[mask, 0], emb[mask, 1], s=5, alpha=0.7,
                   color=cmap[cl], label=f"Cluster {cl}")
    ax.set_xlabel("PC_1"); ax.set_ylabel("PC_2")
    ax.set_title("PCA — coloured by cluster", fontsize=13)
    ax.legend(markerscale=3, fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    _save(fig, "05_pca_dimplot.png")


# ---------------------------------------------------------------------------
# 6. Elbow Plot
# ---------------------------------------------------------------------------

def plot_elbow(pbmc):
    stdev = pbmc.reductions["pca"].stdev
    pcs = np.arange(1, len(stdev) + 1)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(pcs, stdev, "o-", color="#F8766D", markersize=4, linewidth=1.5)
    ax.axvline(10, color="#666666", linestyle="--", linewidth=1, label="PC 10 (selected)")
    ax.set_xlabel("Principal Component"); ax.set_ylabel("Standard Deviation")
    ax.set_title("Elbow Plot — variance explained per PC", fontsize=13)
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    _save(fig, "06_elbow_plot.png")


# ---------------------------------------------------------------------------
# 7. UMAP coloured by cluster
# ---------------------------------------------------------------------------

def plot_umap_clusters(pbmc):
    emb = pbmc.reductions["umap"].cell_embeddings
    clusters = pbmc.meta_data["seurat_clusters"].astype(str).values
    unique = sorted(set(clusters), key=int)
    colors = _cluster_colors(len(unique))
    cmap = dict(zip(unique, colors))

    fig, ax = plt.subplots(figsize=(7, 6))
    for cl in unique:
        mask = clusters == cl
        ax.scatter(emb[mask, 0], emb[mask, 1], s=4, alpha=0.7,
                   color=cmap[cl], label=f"Cluster {cl}")
    ax.set_xlabel("UMAP_1"); ax.set_ylabel("UMAP_2")
    ax.set_title("UMAP — coloured by cluster", fontsize=13)
    ax.legend(markerscale=3, fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    _save(fig, "07_umap_clusters.png")


# ---------------------------------------------------------------------------
# 8. Feature plots — canonical marker genes
# ---------------------------------------------------------------------------

def plot_feature_plots(pbmc):
    canon = ["MS4A1", "CD79A", "NKG7", "GNLY", "FCGR3A", "LYZ", "PPBP", "CD8A", "IL7R"]
    rna = pbmc.assays["RNA"]
    feat_names = rna._all_feature_names
    import scipy.sparse as sp
    data = _get_layer(rna)

    emb = pbmc.reductions["umap"].cell_embeddings

    fig, axes = plt.subplots(3, 3, figsize=(13, 11))
    axes = axes.flatten()
    for i, gene in enumerate(canon):
        ax = axes[i]
        if gene in feat_names:
            idx = feat_names.index(gene)
            row = data[idx, :]
            expr = np.array(row.todense()).flatten() if sp.issparse(data) else np.asarray(row).flatten()
        else:
            expr = np.zeros(emb.shape[0])
        order = np.argsort(expr)
        sc = ax.scatter(emb[order, 0], emb[order, 1], c=expr[order], s=3,
                        cmap="YlOrRd", vmin=0, vmax=max(expr.max(), 0.01))
        plt.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
        ax.set_title(gene, fontsize=11, fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([])
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False); ax.spines["bottom"].set_visible(False)
    fig.suptitle("Canonical Marker Gene Expression on UMAP", fontsize=14, fontweight="bold")
    fig.tight_layout()
    _save(fig, "08_feature_plots.png")


# ---------------------------------------------------------------------------
# 9. Violin plots — key markers
# ---------------------------------------------------------------------------

def plot_marker_violins(pbmc):
    markers = ["MS4A1", "CD79A", "NKG7", "PF4"]
    rna = pbmc.assays["RNA"]
    feat_names = rna._all_feature_names
    data_layer = _get_layer(rna)
    import scipy.sparse as sp

    clusters = pbmc.meta_data["seurat_clusters"].astype(str).values
    unique = sorted(set(clusters), key=int)

    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    axes = axes.flatten()
    for i, gene in enumerate(markers):
        ax = axes[i]
        if gene in feat_names:
            idx = feat_names.index(gene)
            row = data_layer[idx, :]
            expr = np.array(row.todense()).flatten() if sp.issparse(data_layer) else np.asarray(row).flatten()
        else:
            expr = np.zeros(len(clusters))

        df = pd.DataFrame({"expr": expr, "cluster": clusters})
        grp = [df[df["cluster"] == cl]["expr"].values for cl in unique]
        colors = _cluster_colors(len(unique))
        parts = ax.violinplot(grp, positions=range(len(unique)), showmedians=True)
        for j, pc in enumerate(parts["bodies"]):
            pc.set_facecolor(colors[j])
            pc.set_alpha(0.75)
        parts["cmedians"].set_color("black")
        ax.set_xticks(range(len(unique)))
        ax.set_xticklabels(unique)
        ax.set_xlabel("Cluster"); ax.set_ylabel("Expression")
        ax.set_title(gene, fontsize=12, fontweight="bold")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.suptitle("Marker Gene Expression per Cluster", fontsize=14, fontweight="bold")
    fig.tight_layout()
    _save(fig, "09_marker_violins.png")


# ---------------------------------------------------------------------------
# 10. Top marker heatmap
# ---------------------------------------------------------------------------

def plot_marker_heatmap(pbmc, all_markers):
    import scipy.sparse as sp
    rna = pbmc.assays["RNA"]
    feat_names = rna._all_feature_names
    scaled = rna.layers.get("scale.data")

    top_genes = (
        all_markers.groupby("cluster", group_keys=False)
        .apply(lambda x: x.nlargest(5, "avg_log2FC"))
        ["gene"].tolist()
    )
    top_genes = list(dict.fromkeys(top_genes))  # deduplicate, preserve order

    clusters = pbmc.meta_data["seurat_clusters"].astype(str).values
    order = np.argsort(clusters.astype(int))
    sorted_clusters = clusters[order]

    def _extract_row(mat, idx):
        row = mat[idx, :]
        if sp.issparse(row):
            return np.asarray(row.todense()).flatten()
        arr = np.asarray(row).flatten()
        return arr

    rows = []
    for gene in top_genes:
        if gene in feat_names:
            idx = feat_names.index(gene)
            if scaled is not None:
                row = _extract_row(scaled, idx)
            else:
                d = _get_layer(rna)
                row = _extract_row(d, idx)
        else:
            row = np.zeros(len(clusters))
        rows.append(row[order])

    mat = np.array(rows)
    mat = np.clip(mat, -2.5, 2.5)

    fig, ax = plt.subplots(figsize=(14, max(6, len(top_genes) * 0.32)))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5,
                   interpolation="none")
    plt.colorbar(im, ax=ax, shrink=0.5, label="Scaled expression")

    ax.set_yticks(range(len(top_genes)))
    ax.set_yticklabels(top_genes, fontsize=7)

    # Cluster boundary lines + labels
    unique_cl = sorted(set(sorted_clusters), key=int)
    colors = _cluster_colors(len(unique_cl))
    prev = 0
    for cl in unique_cl:
        count = np.sum(sorted_clusters == cl)
        mid = prev + count // 2
        ax.axvline(prev - 0.5, color="white", linewidth=1)
        ax.text(mid, -1.5, cl, ha="center", va="top", fontsize=7, color=colors[int(cl)])
        prev += count

    ax.set_xlabel("Cells (sorted by cluster)")
    ax.set_title("Top 5 Marker Genes per Cluster (scaled expression)", fontsize=13)
    ax.set_xticks([])
    fig.tight_layout()
    _save(fig, "10_marker_heatmap.png")


# ---------------------------------------------------------------------------
# 11. UMAP with cell type labels
# ---------------------------------------------------------------------------

def plot_umap_labeled(pbmc):
    emb = pbmc.reductions["umap"].cell_embeddings
    idents = list(pbmc.idents)
    unique_types = list(dict.fromkeys(idents))
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_types)))
    cmap = dict(zip(unique_types, colors))

    fig, ax = plt.subplots(figsize=(9, 7))
    for ct in unique_types:
        mask = np.array([i == ct for i in idents])
        ax.scatter(emb[mask, 0], emb[mask, 1], s=5, alpha=0.7,
                   color=cmap[ct], label=ct)

    # Centroid labels
    for ct in unique_types:
        mask = np.array([i == ct for i in idents])
        cx, cy = emb[mask, 0].mean(), emb[mask, 1].mean()
        ax.text(cx, cy, ct, fontsize=8, fontweight="bold", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, ec="none"))

    ax.set_xlabel("UMAP_1"); ax.set_ylabel("UMAP_2")
    ax.set_title("UMAP — Cell Type Annotations", fontsize=13)
    ax.legend(markerscale=3, fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    _save(fig, "11_umap_labeled.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(data_dir=None):
    pbmc, hvg, all_markers = run_pipeline(data_dir)

    print("\nGenerating plots...")
    plot_qc_violin(pbmc)
    plot_qc_scatter(pbmc)
    plot_hvg(pbmc)
    plot_pca_loadings(pbmc)
    plot_pca_dimplot(pbmc)
    plot_elbow(pbmc)
    plot_umap_clusters(pbmc)
    plot_feature_plots(pbmc)
    plot_marker_violins(pbmc)
    plot_marker_heatmap(pbmc, all_markers)
    plot_umap_labeled(pbmc)

    print(f"\nAll plots saved to {FIGURES}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate PBMC3k comparison plots")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    main(data_dir=args.data_dir)
