"""Plotting module for Shanuz — Python equivalent of Seurat's plotting functions.

Mirrors the R Seurat plotting API:
  VlnPlot       → vln_plot
  FeaturePlot   → feature_plot
  DimPlot       → dim_plot
  ElbowPlot     → elbow_plot
  FeatureScatter → feature_scatter
  VariableFeaturePlot → variable_feature_plot
  DimHeatmap    → dim_heatmap
  DoHeatmap     → do_heatmap
  RidgePlot     → ridge_plot

All functions return a matplotlib Figure so the caller can save or display it:

    fig = dim_plot(pbmc, reduction="umap", label=True)
    fig.savefig("umap.png", dpi=150, bbox_inches="tight")
    # or in a Jupyter notebook just call the function — the figure displays inline
"""
from __future__ import annotations

from typing import Optional, Union

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Lazy matplotlib import — keeps the package importable without matplotlib
# ---------------------------------------------------------------------------

def _mpl():
    try:
        import matplotlib.pyplot as plt
        return plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for plotting. Install it with:\n"
            "  pip install matplotlib"
        ) from e


def _sns():
    try:
        import seaborn as sns
        return sns
    except ImportError:
        return None


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Seurat-like categorical colour palette (up to 36 groups)
_PALETTE_36 = [
    "#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4",
    "#00A9FF", "#C77CFF", "#FF61CC", "#ABA300", "#00C19A",
    "#E68613", "#0CB702", "#00B8E7", "#35A2FF", "#A3A500",
    "#F564E3", "#FF6C90", "#D39200", "#93AA00", "#00BA38",
    "#00C0AF", "#619CFF", "#DB72FB", "#FF65AC", "#B79F00",
    "#00BE6F", "#F0766E", "#E76BF3", "#00B0F6", "#A3A500",
    "#39B600", "#F8766D", "#00BFC4", "#C77CFF", "#FF61CC",
    "#00B8E7",
]


def _palette(n: int) -> list[str]:
    """Return n categorical colours."""
    if n <= len(_PALETTE_36):
        return _PALETTE_36[:n]
    import matplotlib.pyplot as plt
    cmap = plt.cm.get_cmap("tab20", n)
    return [f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"
            for r, g, b, *_ in (cmap(i) for i in range(n))]


def _get_assay_obj(obj, assay: Optional[str]):
    return obj.assays[assay or obj.active_assay]


def _get_data_matrix(assay_obj, layer: Optional[str] = None):
    """Return the gene-expression matrix (genes × cells) from an assay."""
    from .assay import Assay
    from .assay5 import Assay5
    import scipy.sparse as sp

    if isinstance(assay_obj, Assay5):
        if layer is not None and layer in assay_obj.layers:
            return assay_obj.layers[layer]
        for candidate in ("data", "counts"):
            if candidate in assay_obj.layers:
                return assay_obj.layers[candidate]
        return assay_obj.layers[assay_obj.default_layer]
    else:
        if layer in ("counts",):
            return assay_obj.counts
        elif layer in ("scale.data", "scale_data"):
            return assay_obj.scale_data
        else:
            d = assay_obj.data
            from ._sparse import is_matrix_empty
            return d if not is_matrix_empty(d) else assay_obj.counts


def _get_expression(obj, feature: str, assay: Optional[str] = None,
                    layer: Optional[str] = None) -> np.ndarray:
    """Return a 1-D expression vector for *feature* (gene or metadata col)."""
    import scipy.sparse as sp

    # Metadata columns (nFeature_RNA, percent.mt, nCount_RNA, …)
    if feature in obj.meta_data.columns:
        return obj.meta_data[feature].values.astype(float)

    assay_obj = _get_assay_obj(obj, assay)
    mat = _get_data_matrix(assay_obj, layer)
    feat_names = assay_obj._all_feature_names

    if feature not in feat_names:
        raise KeyError(f"Feature '{feature}' not found in assay or metadata.")

    idx = feat_names.index(feature)
    row = mat[idx, :]
    if sp.issparse(row):
        return np.asarray(row.todense()).flatten()
    return np.asarray(row).flatten()


def _get_embedding(obj, reduction: str) -> np.ndarray:
    """Return cell embeddings (n_cells × 2) for *reduction*."""
    if reduction not in obj.reductions:
        available = list(obj.reductions.keys())
        raise KeyError(
            f"Reduction '{reduction}' not found. Available: {available}"
        )
    return obj.reductions[reduction].cell_embeddings[:, :2]


def _get_groups(obj, group_by: Optional[str]) -> np.ndarray:
    """Return per-cell group labels as an array of strings."""
    if group_by is None:
        return np.array([str(i) for i in obj.idents])
    if group_by in obj.meta_data.columns:
        return obj.meta_data[group_by].astype(str).values
    raise KeyError(f"group_by column '{group_by}' not found in meta_data.")


def _strip_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _subplot_grid(n: int, ncol: Optional[int] = None):
    """Return (nrow, ncol) for n subplots."""
    if ncol is None:
        ncol = min(n, 4 if n > 4 else n)
    nrow = int(np.ceil(n / ncol))
    return nrow, ncol


# ---------------------------------------------------------------------------
# 1. vln_plot — VlnPlot
# ---------------------------------------------------------------------------

def vln_plot(
    obj,
    features: Union[str, list[str]],
    group_by: Optional[str] = None,
    assay: Optional[str] = None,
    layer: Optional[str] = None,
    pt_size: float = 0.0,
    ncol: Optional[int] = None,
    figsize: Optional[tuple] = None,
    palette: Optional[list] = None,
) -> "plt.Figure":
    """Violin plot of feature expression per cluster/identity.

    Mirrors R's ``VlnPlot(pbmc, features = c("LYZ", "CD3D"))``.

    Parameters
    ----------
    obj      : Shanuz object
    features : gene name(s) or metadata column(s)
    group_by : metadata column used for grouping (default: active idents)
    pt_size  : size of individual data points overlaid on violins (0 = none)
    ncol     : number of columns in subplot grid
    figsize  : figure size in inches; auto-computed if None
    palette  : list of colours per group
    """
    plt = _mpl()
    if isinstance(features, str):
        features = [features]

    groups = _get_groups(obj, group_by)
    unique = sorted(set(groups), key=lambda x: (int(x) if x.isdigit() else x))
    colors = palette or _palette(len(unique))
    cmap = dict(zip(unique, colors))

    nrow, nc = _subplot_grid(len(features), ncol)
    if figsize is None:
        figsize = (max(5, nc * 4), nrow * 3.5)

    fig, axes = plt.subplots(nrow, nc, figsize=figsize, squeeze=False)
    axes_flat = axes.flatten()

    for i, feat in enumerate(features):
        ax = axes_flat[i]
        expr = _get_expression(obj, feat, assay, layer)
        grp_data = [expr[groups == g] for g in unique]

        parts = ax.violinplot(grp_data, positions=range(len(unique)),
                              showmedians=True, showextrema=False)
        for j, pc in enumerate(parts["bodies"]):
            pc.set_facecolor(colors[j])
            pc.set_alpha(0.8)
        parts["cmedians"].set_color("black")
        parts["cmedians"].set_linewidth(1.5)

        if pt_size > 0:
            for j, g in enumerate(unique):
                jitter = np.random.uniform(-0.2, 0.2, size=(groups == g).sum())
                ax.scatter(j + jitter, expr[groups == g],
                           s=pt_size, alpha=0.4, color=colors[j], zorder=3)

        ax.set_xticks(range(len(unique)))
        ax.set_xticklabels(unique, rotation=45 if len(unique) > 6 else 0,
                           ha="right" if len(unique) > 6 else "center", fontsize=9)
        ax.set_ylabel("Expression")
        ax.set_title(feat, fontsize=11, fontweight="bold")
        _strip_axes(ax)

    for i in range(len(features), len(axes_flat)):
        axes_flat[i].set_visible(False)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 2. feature_plot — FeaturePlot
# ---------------------------------------------------------------------------

def feature_plot(
    obj,
    features: Union[str, list[str]],
    reduction: str = "umap",
    assay: Optional[str] = None,
    layer: Optional[str] = None,
    ncol: Optional[int] = None,
    order: bool = True,
    min_cutoff: Optional[float] = None,
    max_cutoff: Optional[float] = None,
    colormap: str = "YlOrRd",
    pt_size: float = 3.0,
    figsize: Optional[tuple] = None,
) -> "plt.Figure":
    """Visualise feature expression on a dimensionality reduction embedding.

    Mirrors R's ``FeaturePlot(pbmc, features = c("MS4A1", "LYZ"))``.

    Parameters
    ----------
    features    : gene name(s) or metadata column(s) to plot
    reduction   : which reduction to use ("umap", "pca", …)
    order       : plot cells with highest expression on top
    min_cutoff  : clip expression below this percentile (e.g. "q05")
    max_cutoff  : clip expression above this percentile
    colormap    : matplotlib colormap name for expression
    pt_size     : scatter point size
    """
    plt = _mpl()
    if isinstance(features, str):
        features = [features]

    emb = _get_embedding(obj, reduction)
    ax1_label = reduction.upper() + "_1"
    ax2_label = reduction.upper() + "_2"

    nrow, nc = _subplot_grid(len(features), ncol)
    if figsize is None:
        figsize = (nc * 4.5, nrow * 4)

    fig, axes = plt.subplots(nrow, nc, figsize=figsize, squeeze=False)
    axes_flat = axes.flatten()

    for i, feat in enumerate(features):
        ax = axes_flat[i]
        expr = _get_expression(obj, feat, assay, layer)

        vmin = np.percentile(expr, 5) if min_cutoff == "q05" else (min_cutoff or 0)
        vmax = np.percentile(expr, 95) if max_cutoff == "q95" else (max_cutoff or expr.max())
        vmax = max(vmax, vmin + 1e-9)

        # R Seurat default: non-expressing (zero) cells rendered in light gray,
        # expressing cells drawn on top sorted by expression level.
        zero_mask = expr <= vmin
        nonzero_mask = ~zero_mask
        ax.scatter(emb[zero_mask, 0], emb[zero_mask, 1],
                   c="#D3D3D3", s=pt_size, linewidths=0, rasterized=True)

        expr_nz = expr[nonzero_mask]
        emb_nz = emb[nonzero_mask]
        if order:
            sort_idx = np.argsort(expr_nz)
            expr_nz = expr_nz[sort_idx]
            emb_nz = emb_nz[sort_idx]

        sc = ax.scatter(emb_nz[:, 0], emb_nz[:, 1], c=expr_nz, s=pt_size,
                        cmap=colormap, vmin=vmin, vmax=vmax, linewidths=0,
                        rasterized=True)
        plt.colorbar(sc, ax=ax, shrink=0.6, pad=0.01, aspect=25)
        ax.set_xlabel(ax1_label, fontsize=8)
        ax.set_ylabel(ax2_label, fontsize=8)
        ax.set_title(feat, fontsize=11, fontweight="bold")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    for i in range(len(features), len(axes_flat)):
        axes_flat[i].set_visible(False)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 3. dim_plot — DimPlot
# ---------------------------------------------------------------------------

def dim_plot(
    obj,
    reduction: str = "umap",
    group_by: Optional[str] = None,
    label: bool = True,
    label_size: int = 9,
    pt_size: float = 4.0,
    alpha: float = 0.7,
    figsize: tuple = (7, 6),
    palette: Optional[list] = None,
    title: Optional[str] = None,
) -> "plt.Figure":
    """Plot cells in a reduced-dimension embedding coloured by identity.

    Mirrors R's ``DimPlot(pbmc, reduction = "umap", label = TRUE)``.

    Parameters
    ----------
    reduction : which reduction to use ("umap", "pca", …)
    group_by  : metadata column for colouring (default: active idents)
    label     : add centroid labels for each group
    label_size: font size for centroid labels
    pt_size   : scatter point size
    alpha     : point transparency
    """
    plt = _mpl()
    emb = _get_embedding(obj, reduction)
    groups = _get_groups(obj, group_by)
    unique = sorted(set(groups), key=lambda x: (int(x) if x.isdigit() else x))
    colors = palette or _palette(len(unique))
    cmap = dict(zip(unique, colors))

    ax1_label = reduction.upper() + "_1"
    ax2_label = reduction.upper() + "_2"

    fig, ax = plt.subplots(figsize=figsize)
    for g, color in zip(unique, colors):
        mask = groups == g
        ax.scatter(emb[mask, 0], emb[mask, 1], s=pt_size, alpha=alpha,
                   color=color, label=g, linewidths=0)

    if label:
        for g in unique:
            mask = groups == g
            cx, cy = emb[mask, 0].mean(), emb[mask, 1].mean()
            ax.text(cx, cy, g, fontsize=label_size, fontweight="bold",
                    ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.25", fc="white",
                              alpha=0.75, ec="none"))

    ax.set_xlabel(ax1_label); ax.set_ylabel(ax2_label)
    ax.set_title(title or f"{reduction.upper()} — coloured by "
                          f"{'ident' if group_by is None else group_by}",
                 fontsize=12)
    ax.legend(markerscale=2.5, fontsize=8,
              bbox_to_anchor=(1.01, 1), loc="upper left",
              frameon=False)
    _strip_axes(ax)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 4. elbow_plot — ElbowPlot
# ---------------------------------------------------------------------------

def elbow_plot(
    obj,
    reduction: str = "pca",
    ndims: int = 20,
    figsize: tuple = (7, 4),
) -> "plt.Figure":
    """Rank principal components by standard deviation.

    Mirrors R's ``ElbowPlot(pbmc)``.

    Parameters
    ----------
    reduction : name of the PCA-like reduction
    ndims     : number of PCs to show
    """
    plt = _mpl()
    if reduction not in obj.reductions:
        raise KeyError(f"Reduction '{reduction}' not found.")

    stdev = obj.reductions[reduction].stdev
    n = min(ndims, len(stdev))
    pcs = np.arange(1, n + 1)

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(pcs, stdev[:n], "o-", color="#F8766D", markersize=5, linewidth=1.8)
    ax.set_xlabel("Principal Component")
    ax.set_ylabel("Standard Deviation")
    ax.set_title("Elbow Plot", fontsize=12)
    ax.set_xticks(pcs[::max(1, n // 10)])
    _strip_axes(ax)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 5. feature_scatter — FeatureScatter
# ---------------------------------------------------------------------------

def feature_scatter(
    obj,
    feature1: str,
    feature2: str,
    group_by: Optional[str] = None,
    assay: Optional[str] = None,
    layer: Optional[str] = None,
    pt_size: float = 4.0,
    alpha: float = 0.6,
    figsize: tuple = (6, 5),
    palette: Optional[list] = None,
) -> "plt.Figure":
    """Scatter plot of two features coloured by identity.

    Mirrors R's ``FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")``.

    Parameters
    ----------
    feature1 / feature2 : gene names or metadata columns
    group_by : column for colouring; default: active idents
    """
    plt = _mpl()
    x = _get_expression(obj, feature1, assay, layer)
    y = _get_expression(obj, feature2, assay, layer)
    groups = _get_groups(obj, group_by)
    unique = sorted(set(groups), key=lambda v: (int(v) if v.isdigit() else v))
    colors = palette or _palette(len(unique))

    fig, ax = plt.subplots(figsize=figsize)
    for g, color in zip(unique, colors):
        mask = groups == g
        ax.scatter(x[mask], y[mask], s=pt_size, alpha=alpha,
                   color=color, label=g, linewidths=0)

    ax.set_xlabel(feature1); ax.set_ylabel(feature2)
    ax.set_title(f"{feature1} vs {feature2}", fontsize=12)
    ax.legend(markerscale=2.5, fontsize=8,
              bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)
    _strip_axes(ax)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 6. variable_feature_plot — VariableFeaturePlot
# ---------------------------------------------------------------------------

def variable_feature_plot(
    obj,
    assay: Optional[str] = None,
    log: bool = True,
    label: bool = True,
    n_label: int = 10,
    figsize: tuple = (9, 5),
) -> "plt.Figure":
    """Plot mean expression vs dispersion and highlight variable features.

    Mirrors R's ``VariableFeaturePlot(pbmc)``.

    Parameters
    ----------
    log     : use log10 axes
    label   : annotate the top *n_label* HVGs by name
    n_label : number of top HVGs to label
    """
    import scipy.sparse as sp
    plt = _mpl()

    assay_obj = _get_assay_obj(obj, assay)
    feat_names = assay_obj._all_feature_names
    hvg_set = set(assay_obj.variable_features)
    top_labeled = assay_obj.variable_features[:n_label]
    is_hvg = np.array([f in hvg_set for f in feat_names])

    # Prefer pre-computed VST stats stored by find_variable_features
    md = getattr(assay_obj, "meta_data", None)
    use_std_var = (md is not None
                   and "variances.standardized" in md.columns
                   and "means" in md.columns)

    if use_std_var:
        means = md["means"].values
        y_vals = md["variances.standardized"].values
        y_label = "Standardized Variance"
        # clip extreme standardized values for readability
        y_vals = np.clip(y_vals, 0, np.percentile(y_vals[y_vals > 0], 99.5))
        log = False  # standardized variance on linear scale matches R
    else:
        mat = _get_data_matrix(assay_obj)
        if sp.issparse(mat):
            means = np.array(mat.mean(axis=1)).flatten()
            sq_means = np.array(mat.power(2).mean(axis=1)).flatten()
        else:
            d = np.asarray(mat, dtype=float)
            means = d.mean(axis=1)
            sq_means = (d ** 2).mean(axis=1)
        y_vals = sq_means - means ** 2
        y_label = "Dispersion (log)" if log else "Dispersion"

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(means[~is_hvg], y_vals[~is_hvg], s=3, alpha=0.3,
               color="#AAAAAA", label="Other genes", linewidths=0)
    ax.scatter(means[is_hvg], y_vals[is_hvg], s=4, alpha=0.7,
               color="#F8766D", label="Variable features", linewidths=0)

    if label:
        for gene in top_labeled:
            if gene in feat_names:
                idx = feat_names.index(gene)
                ax.annotate(gene, (means[idx], y_vals[idx]),
                            fontsize=7, xytext=(3, 3),
                            textcoords="offset points", color="#333333")

    if log and not use_std_var:
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Average Expression (log)")
    else:
        ax.set_xlabel("Average Expression")

    ax.set_ylabel(y_label)
    ax.set_title(f"Highly Variable Features  ({len(hvg_set):,} selected)",
                 fontsize=12)
    ax.legend(markerscale=3, fontsize=9)
    _strip_axes(ax)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 7. viz_dim_loadings — VizDimLoadings
# ---------------------------------------------------------------------------

def viz_dim_loadings(
    obj,
    reduction: str = "pca",
    dims: Union[int, list[int]] = [1, 2],
    n_features: int = 15,
    ncol: Optional[int] = None,
    figsize: Optional[tuple] = None,
) -> "plt.Figure":
    """Horizontal bar charts of top positive and negative loading genes per PC.

    Mirrors R's ``VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")``.

    Parameters
    ----------
    dims       : PC indices (1-based) to visualise
    n_features : number of top genes to show per direction (positive + negative)
    """
    plt = _mpl()
    if isinstance(dims, int):
        dims = [dims]
    dims_0 = [d - 1 for d in dims]

    dr = obj.reductions[reduction]
    loadings = dr.feature_loadings
    feat_names = list(dr._feature_names) if hasattr(dr, "_feature_names") else []

    nrow, nc = _subplot_grid(len(dims), ncol)
    if figsize is None:
        figsize = (nc * 4.5, nrow * max(4, n_features * 0.35))

    fig, axes = plt.subplots(nrow, nc, figsize=figsize, squeeze=False)
    axes_flat = axes.flatten()

    for plot_i, (dim, d0) in enumerate(zip(dims, dims_0)):
        ax = axes_flat[plot_i]
        col = loadings[:, d0]

        top_pos_idx = np.argsort(col)[::-1][:n_features]
        top_neg_idx = np.argsort(col)[:n_features]
        # Combine: negatives at top (sorted most negative first), positives below
        idx = np.concatenate([top_pos_idx[::-1], top_neg_idx[::-1]])
        genes = [feat_names[i] if i < len(feat_names) else f"gene_{i}" for i in idx]
        vals  = col[idx]
        colors = ["#F8766D" if v > 0 else "#00BFC4" for v in vals]

        y_pos = np.arange(len(genes))
        ax.barh(y_pos, vals, color=colors, edgecolor="none", height=0.75)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(genes, fontsize=max(6, 9 - n_features // 10))
        ax.axvline(0, color="#333333", linewidth=0.8)
        ax.set_xlabel("Loading score")
        ax.set_title(f"{reduction.upper()} {dim}", fontsize=12, fontweight="bold")
        _strip_axes(ax)

        # Compact legend patches
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color="#F8766D", label="Positive"),
                            Patch(color="#00BFC4", label="Negative")],
                  fontsize=8, frameon=False, loc="lower right")

    for i in range(len(dims), len(axes_flat)):
        axes_flat[i].set_visible(False)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 8. dim_heatmap — DimHeatmap
# ---------------------------------------------------------------------------

def dim_heatmap(
    obj,
    reduction: str = "pca",
    dims: Union[int, list[int]] = 1,
    cells: int = 500,
    balanced: bool = True,
    ncol: Optional[int] = None,
    figsize: Optional[tuple] = None,
) -> "plt.Figure":
    """Heatmap of gene loadings for selected principal components.

    Shows the most extreme cells (highest / lowest scores) and the top
    loading genes for each PC — mirrors R's ``DimHeatmap(pbmc, dims = 1:6)``.

    Parameters
    ----------
    dims     : PC index (1-based int) or list of indices
    cells    : number of extreme cells to show per PC
    balanced : if True, take equal numbers from both extremes of the PC score
    """
    import scipy.sparse as sp
    plt = _mpl()

    if isinstance(dims, int):
        dims = [dims]
    dims_0 = [d - 1 for d in dims]  # convert to 0-based

    dr = obj.reductions[reduction]
    emb = dr.cell_embeddings
    loadings = dr.feature_loadings
    feat_names = dr._feature_names if hasattr(dr, "_feature_names") else []

    n_top = 10
    nrow, nc = _subplot_grid(len(dims), ncol)
    if figsize is None:
        figsize = (nc * 5, nrow * 5)

    fig, axes = plt.subplots(nrow, nc, figsize=figsize, squeeze=False)
    axes_flat = axes.flatten()

    for plot_i, (dim, d0) in enumerate(zip(dims, dims_0)):
        ax = axes_flat[plot_i]
        scores = emb[:, d0]
        col_loads = loadings[:, d0]

        # Select extreme cells
        n_half = cells // 2
        if balanced:
            top_cells = np.argsort(scores)[-n_half:]
            bot_cells = np.argsort(scores)[:n_half]
            sel_cells = np.concatenate([bot_cells, top_cells])
        else:
            sel_cells = np.argsort(np.abs(scores))[-cells:]
            sel_cells = sel_cells[np.argsort(scores[sel_cells])]

        # Select top-loading genes
        top_genes_idx = np.concatenate([
            np.argsort(col_loads)[::-1][:n_top],
            np.argsort(col_loads)[:n_top]
        ])
        top_genes_idx = np.unique(top_genes_idx)
        top_genes_idx = top_genes_idx[np.argsort(col_loads[top_genes_idx])[::-1]]

        # Build expression matrix: top genes × selected cells
        assay_obj = _get_assay_obj(obj, None)
        mat = _get_data_matrix(assay_obj, "scale.data")
        all_feats = assay_obj._all_feature_names

        rows = []
        gene_labels = []
        for gi in top_genes_idx:
            if gi < len(feat_names):
                gname = feat_names[gi]
            else:
                continue
            if gname in all_feats:
                aidx = all_feats.index(gname)
                row = mat[aidx, :]
                row = np.asarray(row.todense()).flatten() if sp.issparse(row) else np.asarray(row).flatten()
                rows.append(row[sel_cells])
                gene_labels.append(gname)

        if not rows:
            ax.set_visible(False)
            continue

        mat_sub = np.array(rows)
        mat_sub = np.clip(mat_sub, -2.5, 2.5)

        im = ax.imshow(mat_sub, aspect="auto", cmap="RdBu_r",
                       vmin=-2.5, vmax=2.5, interpolation="none")
        ax.set_yticks(range(len(gene_labels)))
        ax.set_yticklabels(gene_labels, fontsize=7)
        ax.set_xticks([])
        if balanced:
            ax.axvline(n_half - 0.5, color="white", linewidth=1.5)
        ax.set_title(f"PC {dim}", fontsize=11)
        plt.colorbar(im, ax=ax, shrink=0.6, label="Scaled expr.")

    for i in range(len(dims), len(axes_flat)):
        axes_flat[i].set_visible(False)

    fig.suptitle(f"DimHeatmap — {reduction.upper()}", fontsize=13,
                 fontweight="bold", y=1.01)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 8. do_heatmap — DoHeatmap
# ---------------------------------------------------------------------------

def do_heatmap(
    obj,
    features: list[str],
    group_by: Optional[str] = None,
    assay: Optional[str] = None,
    layer: str = "scale.data",
    label: bool = True,
    figsize: Optional[tuple] = None,
    palette: Optional[list] = None,
) -> "plt.Figure":
    """Expression heatmap of selected genes across cells, sorted by cluster.

    Mirrors R's ``DoHeatmap(pbmc, features = top10$gene)``.

    Parameters
    ----------
    features : list of gene names to show as rows
    group_by : column used to sort and colour cells (default: active idents)
    layer    : which data layer to use (default: "scale.data")
    label    : annotate cluster boundaries with group names
    """
    import scipy.sparse as sp
    plt = _mpl()

    assay_obj = _get_assay_obj(obj, assay)
    mat = _get_data_matrix(assay_obj, layer)
    all_feats = assay_obj._all_feature_names

    groups = _get_groups(obj, group_by)
    unique = sorted(set(groups), key=lambda x: (int(x) if x.isdigit() else x))
    colors = palette or _palette(len(unique))

    cell_order = np.argsort([unique.index(g) for g in groups])
    sorted_groups = groups[cell_order]

    def _extract(gene):
        if gene not in all_feats:
            return np.zeros(len(groups))
        idx = all_feats.index(gene)
        row = mat[idx, :]
        arr = np.asarray(row.todense()).flatten() if sp.issparse(row) else np.asarray(row).flatten()
        return arr

    rows = [_extract(g)[cell_order] for g in features]
    matrix = np.clip(np.array(rows), -2.5, 2.5)

    if figsize is None:
        figsize = (max(10, len(groups) / 200), min(14, max(4, len(features) * 0.28)))

    fig, axes = plt.subplots(
        2, 1, figsize=figsize,
        gridspec_kw={"height_ratios": [0.05, 1], "hspace": 0.02}
    )

    # Top colour bar showing cluster membership
    colour_row = np.array([unique.index(g) for g in sorted_groups])[np.newaxis, :]
    axes[0].imshow(colour_row, aspect="auto",
                   cmap=plt.matplotlib.colors.ListedColormap(colors),
                   vmin=0, vmax=len(unique) - 1, interpolation="none")
    axes[0].set_xticks([]); axes[0].set_yticks([])
    for spine in axes[0].spines.values():
        spine.set_visible(False)

    # Main heatmap
    im = axes[1].imshow(matrix, aspect="auto", cmap="RdBu_r",
                        vmin=-2.5, vmax=2.5, interpolation="none")
    axes[1].set_yticks(range(len(features)))
    axes[1].set_yticklabels(features, fontsize=max(5, 9 - len(features) // 20))
    axes[1].set_xticks([])

    # Cluster boundary lines + labels inside the colour bar row
    if label:
        prev = 0
        for gi, g in enumerate(unique):
            count = np.sum(sorted_groups == g)
            mid = prev + count / 2
            # Separator lines in both the colour bar and the heatmap
            axes[0].axvline(prev - 0.5, color="white", linewidth=1.2)
            axes[1].axvline(prev - 0.5, color="white", linewidth=0.8)
            # Labels sit inside axes[0] (the colour bar) — y=0.5 centres vertically
            axes[0].text(mid, 0.5, g, ha="center", va="center",
                         fontsize=7, color="white", fontweight="bold")
            prev += count

    plt.colorbar(im, ax=axes[1], shrink=0.4, pad=0.01, label="Scaled expression")
    axes[1].set_xlabel("Cells (sorted by cluster)")
    fig.suptitle("Expression Heatmap", fontsize=13, fontweight="bold", y=1.0)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 9. ridge_plot — RidgePlot
# ---------------------------------------------------------------------------

def ridge_plot(
    obj,
    features: Union[str, list[str]],
    group_by: Optional[str] = None,
    assay: Optional[str] = None,
    layer: Optional[str] = None,
    ncol: Optional[int] = None,
    figsize: Optional[tuple] = None,
    palette: Optional[list] = None,
) -> "plt.Figure":
    """Ridgeline (joy) plots of feature expression per group.

    Mirrors R's ``RidgePlot(pbmc, features = c("LYZ", "CD3D"))``.

    Requires ``scipy`` for KDE smoothing.
    """
    from scipy.stats import gaussian_kde
    plt = _mpl()

    if isinstance(features, str):
        features = [features]

    groups = _get_groups(obj, group_by)
    unique = sorted(set(groups), key=lambda x: (int(x) if x.isdigit() else x))[::-1]
    colors = palette or _palette(len(unique))

    nrow, nc = _subplot_grid(len(features), ncol)
    if figsize is None:
        figsize = (nc * 5, nrow * max(3, len(unique) * 0.5))

    fig, axes = plt.subplots(nrow, nc, figsize=figsize, squeeze=False)
    axes_flat = axes.flatten()

    for i, feat in enumerate(features):
        ax = axes_flat[i]
        expr = _get_expression(obj, feat, assay, layer)
        xmin, xmax = expr.min(), expr.max()
        xs = np.linspace(xmin - 0.1, xmax + 0.1, 300)

        for j, (g, color) in enumerate(zip(unique, colors[::-1])):
            vals = expr[groups == g]
            if vals.std() < 1e-9 or len(vals) < 3:
                continue
            kde = gaussian_kde(vals, bw_method="scott")
            ys = kde(xs)
            ys = ys / ys.max() * 0.9  # normalise height
            ax.fill_between(xs, j + ys, j, alpha=0.8, color=color)
            ax.plot(xs, j + ys, color="white", linewidth=0.5)

        ax.set_yticks(range(len(unique)))
        # Row j holds the density for unique[j] (see loop above), so the tick
        # labels must be `unique`, not its reverse.
        ax.set_yticklabels(unique, fontsize=9)
        ax.set_xlabel("Expression")
        ax.set_title(feat, fontsize=11, fontweight="bold")
        _strip_axes(ax)

    for i in range(len(features), len(axes_flat)):
        axes_flat[i].set_visible(False)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

__all__ = [
    "vln_plot",
    "feature_plot",
    "dim_plot",
    "elbow_plot",
    "feature_scatter",
    "variable_feature_plot",
    "viz_dim_loadings",
    "dim_heatmap",
    "do_heatmap",
    "ridge_plot",
]
