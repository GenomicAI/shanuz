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
  PlotPerturbScore  → plot_perturb_score
  MixscapeHeatmap   → mixscape_heatmap

All functions return a matplotlib Figure so the caller can save or display it:

    fig = dim_plot(pbmc, reduction="umap", label=True)
    fig.savefig("umap.png", dpi=150, bbox_inches="tight")
    # or in a Jupyter notebook just call the function — the figure displays inline
"""
from __future__ import annotations

import re
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


# Seurat's mixscape plots default to R colour names that matplotlib does not know.
# Keeping the R spellings in the signatures (so the call reads like the R one) means
# translating them here.
_R_COLOURS = {"orange2": "#EE9A00"}


def _r_colour(name):
    """Translate an R colour name to a matplotlib-readable one, else pass through."""
    if not isinstance(name, str):
        return name
    if name in _R_COLOURS:
        return _R_COLOURS[name]
    m = re.fullmatch(r"gr[ea]y(\d{1,3})", name)     # greyNN is NN% grey in R
    if m and int(m.group(1)) <= 100:
        level = round(int(m.group(1)) * 255 / 100)
        return f"#{level:02x}{level:02x}{level:02x}"
    return name


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


def _resolve_layer(assay_obj, layer: Optional[str] = None):
    """Return ``(matrix, feature_names)`` for *layer*, rows and names aligned.

    ``scale_data(obj, features = [...])`` scales only a subset, leaving a
    ``scale.data`` layer with fewer rows than the assay has features — so row *i*
    of the matrix is not feature *i* of the assay. Callers that index rows by gene
    must go through the layer's own feature list, which is what this returns.
    """
    from .assay5 import Assay5

    mat = _get_data_matrix(assay_obj, layer)
    if isinstance(assay_obj, Assay5):
        for name, candidate in assay_obj.layers.items():
            if candidate is mat:
                return mat, list(
                    assay_obj._layer_features.get(name, assay_obj._all_feature_names)
                )
        return mat, list(assay_obj._all_feature_names)
    if layer in ("scale.data", "scale_data"):
        return mat, list(assay_obj.features(layer="scale_data"))
    return mat, list(assay_obj.features())


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
    cells: Optional[list[str]] = None,
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
    cells    : restrict to these cells *and* show them in this exact order,
               rather than the default group-sorted order. Mirrors R's ``cells``
               argument; :func:`mixscape_heatmap` uses it to order cells by
               knockout probability.
    """
    import scipy.sparse as sp
    plt = _mpl()

    assay_obj = _get_assay_obj(obj, assay)
    mat, all_feats = _resolve_layer(assay_obj, layer)

    groups = _get_groups(obj, group_by)
    unique = sorted(set(groups), key=lambda x: (int(x) if x.isdigit() else x))
    colors = palette or _palette(len(unique))

    if cells is None:
        cell_order = np.argsort([unique.index(g) for g in groups])
    else:
        pos = {c: i for i, c in enumerate(obj.cell_names())}
        missing = [c for c in cells if c not in pos]
        if missing:
            raise KeyError(f"Cells not in object: {missing[:5]}")
        cell_order = np.array([pos[c] for c in cells], dtype=int)
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

    # Cluster boundary lines + labels inside the colour bar row.
    # Walk contiguous runs rather than whole groups: with the default sort each
    # group is one run, but an explicit `cells` order (e.g. by knockout
    # probability) interleaves them, and a run-based pass labels both correctly.
    if label:
        boundaries = [0] + [
            i for i in range(1, len(sorted_groups))
            if sorted_groups[i] != sorted_groups[i - 1]
        ] + [len(sorted_groups)]
        # Interleaved orders produce many short runs; label only the ones wide
        # enough to read, so the colour bar does not fill with overlapping text.
        min_run = 0 if cells is None else max(1, len(sorted_groups) // 50)
        for start, end in zip(boundaries[:-1], boundaries[1:]):
            # Separator lines in both the colour bar and the heatmap
            axes[0].axvline(start - 0.5, color="white", linewidth=1.2)
            axes[1].axvline(start - 0.5, color="white", linewidth=0.8)
            if end - start >= min_run:
                # Labels sit inside axes[0] (the colour bar) — y=0.5 centres vertically
                axes[0].text((start + end) / 2, 0.5, sorted_groups[start],
                             ha="center", va="center",
                             fontsize=7, color="white", fontweight="bold")

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
# dot_plot — DotPlot
# ---------------------------------------------------------------------------

def dot_plot(
    obj,
    features: Union[str, list[str]],
    group_by: Optional[str] = None,
    assay: Optional[str] = None,
    layer: str = "data",
    cols: tuple = ("lightgrey", "blue"),
    col_min: float = -2.5,
    col_max: float = 2.5,
    dot_min: float = 0.0,
    dot_scale: float = 6.0,
    scale: bool = True,
    figsize: Optional[tuple] = None,
) -> "plt.Figure":
    """Dot plot of feature expression across groups.

    Mirrors R's ``DotPlot(pbmc, features = c("LYZ", "CD3D"))``. For each
    (feature, group) the **dot size** encodes the fraction of cells in the group
    expressing the feature (counts > 0) and the **colour** encodes the average
    expression (z-scored across groups when ``scale=True``, matching Seurat).

    Parameters
    ----------
    features : gene name(s) to plot (x-axis).
    group_by : metadata column for grouping (default: active idents, y-axis).
    scale    : z-score each feature's average expression across groups.
    col_min/col_max : colour scale limits for the scaled average expression.
    dot_min  : minimum fraction-expressing to draw a dot.
    dot_scale: scaling factor for dot area.
    """
    plt = _mpl()
    import matplotlib as mpl

    if isinstance(features, str):
        features = [features]

    groups = _get_groups(obj, group_by)
    unique = sorted(set(groups), key=lambda x: (int(x) if x.isdigit() else x))

    # Per (feature, group): average expression and fraction expressing.
    avg = np.zeros((len(features), len(unique)))
    pct = np.zeros((len(features), len(unique)))
    for fi, feat in enumerate(features):
        expr = _get_expression(obj, feat, assay, layer)
        for gi, g in enumerate(unique):
            cell_vals = expr[groups == g]
            if cell_vals.size == 0:
                continue
            avg[fi, gi] = cell_vals.mean()
            pct[fi, gi] = float((cell_vals > 0).mean())

    # Scale average expression across groups (z-score per feature), like Seurat.
    color_vals = avg.copy()
    if scale:
        mu = avg.mean(axis=1, keepdims=True)
        sd = avg.std(axis=1, ddof=1, keepdims=True)
        sd[sd == 0] = 1.0
        color_vals = np.clip((avg - mu) / sd, col_min, col_max)
        vmin, vmax = col_min, col_max
    else:
        # log1p of average expression, as Seurat does for the unscaled case.
        color_vals = np.log1p(avg)
        vmin, vmax = float(color_vals.min()), float(color_vals.max())

    cmap = mpl.colors.LinearSegmentedColormap.from_list("dotplot", list(cols))
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    if figsize is None:
        figsize = (max(5, 0.7 * len(features) + 2), max(3, 0.5 * len(unique) + 1.5))
    fig, ax = plt.subplots(figsize=figsize)

    for fi in range(len(features)):
        for gi in range(len(unique)):
            frac = pct[fi, gi]
            if frac < dot_min or frac == 0:
                continue
            ax.scatter(
                fi, gi, s=(frac * dot_scale) ** 2 * np.pi,
                c=[cmap(norm(color_vals[fi, gi]))],
                edgecolors="black", linewidths=0.3, zorder=3,
            )

    ax.set_xticks(range(len(features)))
    ax.set_xticklabels(features, rotation=90, fontsize=9)
    ax.set_yticks(range(len(unique)))
    ax.set_yticklabels(unique, fontsize=9)
    ax.set_xlim(-0.5, len(features) - 0.5)
    ax.set_ylim(-0.5, len(unique) - 0.5)
    ax.set_axisbelow(True)
    ax.grid(True, color="0.9", linewidth=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # Colour bar (average expression) and a size legend (% expressed).
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.02, aspect=12)
    cbar.set_label("Average expression" + (" (scaled)" if scale else ""), fontsize=8)

    for frac in (0.25, 0.5, 0.75, 1.0):
        ax.scatter([], [], s=(frac * dot_scale) ** 2 * np.pi, c="grey",
                   edgecolors="black", linewidths=0.3, label=f"{int(frac * 100)}%")
    ax.legend(title="% expressed", bbox_to_anchor=(1.22, 0.0), loc="lower left",
              labelspacing=1.2, frameon=False, fontsize=8, title_fontsize=8)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Spatial plots — ImageDimPlot / ImageFeaturePlot
# ---------------------------------------------------------------------------

def image_dim_plot(
    obj,
    group_by: Optional[str] = None,
    image: Optional[Union[str, list]] = None,
    size: float = 1.0,
    cols: Optional[dict] = None,
    ncol: Optional[int] = None,
    flip_y: bool = True,
    figsize: Optional[tuple] = None,
) -> "plt.Figure":
    """Plot cell centroids in physical space, coloured by a grouping variable.

    Matplotlib equivalent of Seurat's ``ImageDimPlot`` — drawn directly from
    centroids (``obj.get_tissue_coordinates``), so it is immune to the
    ggplot2-4.x blank-render issue. One panel per image.

    Parameters
    ----------
    group_by : metadata column (default: active idents) for point colour.
    image    : image name(s) to draw (default: all).
    cols     : optional ``{group: colour}`` mapping.
    flip_y   : invert the y-axis so images match Seurat's orientation.
    """
    plt = _mpl()
    coords = obj.get_tissue_coordinates(image)
    if coords.empty:
        raise ValueError("Object has no spatial coordinates to plot.")

    labels = _get_groups(obj, group_by)
    lab_by_cell = dict(zip(obj.cell_names(), labels))
    coords = coords.assign(group=coords["cell"].map(lab_by_cell).astype(str))

    images = list(dict.fromkeys(coords["image"]))
    uniq = sorted(coords["group"].unique())
    if cols is None:
        cols = dict(zip(uniq, _palette(len(uniq))))

    nrow, ncol = _subplot_grid(len(images), ncol)
    if figsize is None:
        figsize = (5 * ncol, 4.5 * nrow)
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize, squeeze=False)
    axes_flat = axes.ravel()

    for ax, img in zip(axes_flat, images):
        d = coords[coords["image"] == img]
        for g, dd in d.groupby("group"):
            ax.scatter(dd["x"], dd["y"], s=size, c=[cols.get(g, "grey")],
                       label=g, linewidths=0)
        ax.set_title(str(img), fontsize=9, fontweight="bold")
        ax.set_aspect("equal"); ax.axis("off")
        if flip_y:
            ax.invert_yaxis()
    for ax in axes_flat[len(images):]:
        ax.axis("off")

    handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=6,
                          markerfacecolor=cols.get(g, "grey"), markeredgewidth=0)
               for g in uniq]
    fig.legend(handles, uniq, title=group_by or "ident", loc="center right",
               fontsize=8, title_fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0, 0.88, 1))
    return fig


def image_feature_plot(
    obj,
    feature: str,
    image: Optional[Union[str, list]] = None,
    size: float = 1.0,
    cmap: str = "viridis",
    assay: Optional[str] = None,
    layer: Optional[str] = None,
    ncol: Optional[int] = None,
    flip_y: bool = True,
    figsize: Optional[tuple] = None,
) -> "plt.Figure":
    """Plot cell centroids in physical space, coloured by a feature's expression.

    Matplotlib equivalent of Seurat's ``ImageFeaturePlot``. One panel per image.
    """
    plt = _mpl()
    coords = obj.get_tissue_coordinates(image)
    if coords.empty:
        raise ValueError("Object has no spatial coordinates to plot.")

    expr = _get_expression(obj, feature, assay, layer)
    expr_by_cell = dict(zip(obj.cell_names(), expr))
    coords = coords.assign(val=coords["cell"].map(expr_by_cell).astype(float))

    images = list(dict.fromkeys(coords["image"]))
    nrow, ncol = _subplot_grid(len(images), ncol)
    if figsize is None:
        figsize = (5 * ncol, 4.5 * nrow)
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize, squeeze=False)
    axes_flat = axes.ravel()
    vmax = float(np.nanmax(coords["val"])) if len(coords) else 1.0

    sc = None
    for ax, img in zip(axes_flat, images):
        d = coords[coords["image"] == img]
        sc = ax.scatter(d["x"], d["y"], s=size, c=d["val"], cmap=cmap,
                        vmin=0, vmax=vmax, linewidths=0)
        ax.set_title(str(img), fontsize=9, fontweight="bold")
        ax.set_aspect("equal"); ax.axis("off")
        if flip_y:
            ax.invert_yaxis()
    for ax in axes_flat[len(images):]:
        ax.axis("off")
    if sc is not None:
        fig.colorbar(sc, ax=axes_flat.tolist(), shrink=0.5, label=feature)
    return fig


# ---------------------------------------------------------------------------
# Visium tissue-image plots — SpatialDimPlot / SpatialFeaturePlot
# ---------------------------------------------------------------------------

def _resolve_fovs(obj, image: Optional[Union[str, list]] = None) -> dict:
    """The requested image slots, in order, as a ``{name: fov}`` dict."""
    images = getattr(obj, "images", None) or {}
    if not images:
        raise ValueError("Object has no spatial images to plot.")
    if image is None:
        names = list(images)
    elif isinstance(image, str):
        names = [image]
    else:
        names = list(image)
    missing = [n for n in names if n not in images]
    if missing:
        raise KeyError(
            f"No such image(s): {missing}. Available: {list(images)}."
        )
    return {n: images[n] for n in names}


def _spatial_panel(fov, resolution: Optional[str]):
    """Draw-space coordinates, spot radius and background image for one FOV.

    A ``VisiumV2`` carrying a tissue photo reports its coordinates in
    full-resolution pixels, so they are scaled into the image's own pixel space
    here. Any other FOV has no image; its coordinates are used as they stand.
    """
    img = fov.get_image()
    if img is None:
        return fov.get_tissue_coordinates(), fov.radius(), None
    return fov.scale_coordinates(resolution=resolution), fov.spot_radius(resolution), img


def _spot_collection(ax, coords, radius: Optional[float], pt_size_factor: float):
    """Spots as true-to-scale circles, or None when the spot size is unknown."""
    if radius is None or not np.isfinite(radius) or radius <= 0:
        return None
    from matplotlib.collections import EllipseCollection

    d = 2.0 * float(radius) * pt_size_factor
    offsets = np.column_stack([
        coords["x"].to_numpy(dtype=float),
        coords["y"].to_numpy(dtype=float),
    ])
    return EllipseCollection(
        widths=d, heights=d, angles=0.0, units="xy",
        offsets=offsets, offset_transform=ax.transData, linewidths=0.0,
    )


def _spatial_limits(ax, coords, radius: Optional[float], img, crop: bool) -> None:
    """Frame the panel, always with y increasing downward (image convention)."""
    if img is not None and not crop:
        h, w = img.shape[:2]
        ax.set_xlim(-0.5, w - 0.5)
        ax.set_ylim(h - 0.5, -0.5)
        return
    x = coords["x"].to_numpy(dtype=float)
    y = coords["y"].to_numpy(dtype=float)
    pad = (radius or 0.0) + 0.02 * max(float(np.ptp(x)), float(np.ptp(y)), 1.0)
    ax.set_xlim(x.min() - pad, x.max() + pad)
    ax.set_ylim(y.max() + pad, y.min() - pad)


def spatial_dim_plot(
    obj,
    group_by: Optional[str] = None,
    image: Optional[Union[str, list]] = None,
    cols: Optional[dict] = None,
    pt_size_factor: float = 1.6,
    size: float = 12.0,
    alpha: float = 1.0,
    image_alpha: float = 1.0,
    resolution: Optional[str] = None,
    crop: bool = True,
    ncol: Optional[int] = None,
    figsize: Optional[tuple] = None,
) -> "plt.Figure":
    """Plot spots on the tissue image, coloured by a grouping variable.

    Matplotlib equivalent of Seurat's ``SpatialDimPlot``. Each panel draws the
    H&E photo held by a :class:`~shanuz.spatial.visium.VisiumV2` image and
    overlays the spots at their true diameter. An image slot with no photo (a
    plain ``FOV``, or a Visium bundle loaded with ``image=False``) degrades to a
    bare scatter of the same spots — the plot still works, it just has no
    tissue underneath.

    Parameters
    ----------
    group_by : metadata column (default: active idents) for spot colour.
    image    : image name(s) to draw (default: all).
    cols     : optional ``{group: colour}`` mapping.
    pt_size_factor : scales the spots relative to their real diameter, as in
        Seurat. Ignored when the image has no ``scalefactors_json.json``, in
        which case ``size`` (a plain scatter point size) applies instead.
    image_alpha : opacity of the tissue photo — drop it to make spots pop.
    resolution : ``"hires"`` / ``"lowres"``; defaults to whatever was loaded.
    crop     : zoom to the spots rather than showing the whole slide.
    """
    plt = _mpl()
    fovs = _resolve_fovs(obj, image)

    labels = _get_groups(obj, group_by)
    lab_by_cell = dict(zip(obj.cell_names(), [str(v) for v in labels]))

    panels = {}
    for name, fov in fovs.items():
        coords, radius, img = _spatial_panel(fov, resolution)
        coords = coords.assign(group=[lab_by_cell.get(c) for c in coords.index])
        coords = coords[coords["group"].notna()]
        if not coords.empty:
            panels[name] = (coords, radius, img)
    if not panels:
        raise ValueError("No spots left to plot: the images share no cells with the object.")

    uniq = sorted({g for coords, _, _ in panels.values() for g in coords["group"]})
    if cols is None:
        cols = dict(zip(uniq, _palette(len(uniq))))

    nrow, ncol = _subplot_grid(len(panels), ncol)
    if figsize is None:
        figsize = (5 * ncol, 4.5 * nrow)
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize, squeeze=False)
    axes_flat = axes.ravel()

    for ax, (name, (coords, radius, img)) in zip(axes_flat, panels.items()):
        if img is not None:
            ax.imshow(img, alpha=image_alpha)
        face = [cols.get(g, "grey") for g in coords["group"]]
        coll = _spot_collection(ax, coords, radius, pt_size_factor)
        if coll is None:
            ax.scatter(coords["x"], coords["y"], s=size, c=face,
                       alpha=alpha, linewidths=0)
        else:
            coll.set_facecolor(face)
            coll.set_alpha(alpha)
            ax.add_collection(coll)
        _spatial_limits(ax, coords, radius, img, crop)
        ax.set_title(str(name), fontsize=9, fontweight="bold")
        ax.set_aspect("equal")
        ax.axis("off")
    for ax in axes_flat[len(panels):]:
        ax.axis("off")

    handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=6,
                          markerfacecolor=cols.get(g, "grey"), markeredgewidth=0)
               for g in uniq]
    fig.legend(handles, uniq, title=group_by or "ident", loc="center right",
               fontsize=8, title_fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0, 0.88, 1))
    return fig


def spatial_feature_plot(
    obj,
    feature: str,
    image: Optional[Union[str, list]] = None,
    cmap: str = "viridis",
    pt_size_factor: float = 1.6,
    size: float = 12.0,
    alpha: float = 1.0,
    image_alpha: float = 1.0,
    assay: Optional[str] = None,
    layer: Optional[str] = None,
    resolution: Optional[str] = None,
    crop: bool = True,
    ncol: Optional[int] = None,
    figsize: Optional[tuple] = None,
) -> "plt.Figure":
    """Plot spots on the tissue image, coloured by a feature's expression.

    Matplotlib equivalent of Seurat's ``SpatialFeaturePlot``. Behaves exactly
    like :func:`spatial_dim_plot` — same tissue background, same true-to-scale
    spots, same fallback when no image is stored — but colours the spots by a
    continuous value on a shared scale across panels.
    """
    plt = _mpl()
    fovs = _resolve_fovs(obj, image)

    expr = _get_expression(obj, feature, assay, layer)
    expr_by_cell = dict(zip(obj.cell_names(), expr))

    panels = {}
    for name, fov in fovs.items():
        coords, radius, img = _spatial_panel(fov, resolution)
        coords = coords.assign(val=[expr_by_cell.get(c, np.nan) for c in coords.index])
        coords = coords[coords["val"].notna()]
        if not coords.empty:
            panels[name] = (coords, radius, img)
    if not panels:
        raise ValueError("No spots left to plot: the images share no cells with the object.")

    vals = np.concatenate([c["val"].to_numpy(dtype=float) for c, _, _ in panels.values()])
    vmax = float(np.nanmax(vals)) if vals.size else 1.0
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0

    nrow, ncol = _subplot_grid(len(panels), ncol)
    if figsize is None:
        figsize = (5 * ncol, 4.5 * nrow)
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize, squeeze=False)
    axes_flat = axes.ravel()

    mappable = None
    for ax, (name, (coords, radius, img)) in zip(axes_flat, panels.items()):
        if img is not None:
            ax.imshow(img, alpha=image_alpha)
        v = coords["val"].to_numpy(dtype=float)
        coll = _spot_collection(ax, coords, radius, pt_size_factor)
        if coll is None:
            mappable = ax.scatter(coords["x"], coords["y"], s=size, c=v, cmap=cmap,
                                  vmin=0, vmax=vmax, alpha=alpha, linewidths=0)
        else:
            coll.set_array(v)
            coll.set_cmap(cmap)
            coll.set_clim(0, vmax)
            coll.set_alpha(alpha)
            ax.add_collection(coll)
            mappable = coll
        _spatial_limits(ax, coords, radius, img, crop)
        ax.set_title(str(name), fontsize=9, fontweight="bold")
        ax.set_aspect("equal")
        ax.axis("off")
    for ax in axes_flat[len(panels):]:
        ax.axis("off")

    if mappable is not None:
        fig.colorbar(mappable, ax=axes_flat.tolist(), shrink=0.5, label=feature)
    return fig


# ---------------------------------------------------------------------------
# Mixscape — PlotPerturbScore / MixscapeHeatmap
# ---------------------------------------------------------------------------

def _mixscape_scores(obj, target_gene_ident: str, assay: str) -> pd.DataFrame:
    """The stored perturbation-score frame for one target gene.

    R reads this from ``Tool(object, slot = "RunMixscape")``; :func:`run_mixscape`
    keeps the same frame under ``obj.misc["mixscape"][assay]["genes"][gene]``.
    """
    store = obj.misc.get("mixscape", {}).get(assay)
    if store is None:
        raise KeyError(
            f"No mixscape results for assay {assay!r}. Run run_mixscape(obj, "
            f"assay={assay!r}) first."
        )
    genes = store["genes"]
    if target_gene_ident not in genes:
        raise KeyError(
            f"{target_gene_ident!r} is not a target gene in this screen. "
            f"Available: {sorted(genes)}"
        )
    scores = genes[target_gene_ident].get("scores")
    if scores is None:
        raise ValueError(
            f"No perturbation score stored for {target_gene_ident!r} — its cells "
            "were called NP without a mixture fit (too few cells, or too few DE "
            "genes), so there is no score axis to plot."
        )
    return scores


def plot_perturb_score(
    obj,
    target_gene_ident: str,
    target_gene_class: str = "gene",
    mixscape_class: str = "mixscape_class",
    col: str = "orange2",
    split_by: Optional[str] = None,
    before_mixscape: bool = False,
    prtb_type: str = "KO",
    assay: str = "PRTB",
    seed: int = 0,
    figsize: Optional[tuple] = None,
) -> "plt.Figure":
    """Density of one target gene's perturbation scores (Seurat's ``PlotPerturbScore``).

    Mirrors ``PlotPerturbScore(object, target.gene.ident = "IFNGR2",
    target.gene.class = "gene", mixscape.class = "mixscape_class")``. This is the
    diagnostic that shows *why* mixscape split a guide the way it did: the score
    is each cell's projection onto the gene's perturbation axis, and the plot
    overlays the NT control density against the guide's own. A guide with a real
    effect is bimodal — one lobe sitting on the NT curve (the escapers) and one
    shifted away from it (the knockouts) — which is exactly the structure the
    mixture model is asked to find.

    Two views, per R's ``before.mixscape``:

    * ``before_mixscape=False`` (default) — colour by ``mixscape_class``, i.e.
      *after* the call: NT, ``"<gene> NP"``, and ``"<gene> KO"`` get their own
      curves, so you see where mixscape actually drew the line.
    * ``before_mixscape=True`` — colour by the raw guide label only, the view you
      would have had without mixscape: NT against the whole guide population.

    Cells are also drawn as a jittered strip — controls above the axis, the target
    gene below — so single-cell density is visible where the curves overlap.

    Parameters
    ----------
    target_gene_ident : the target gene to plot (must be one mixscape tested).
    target_gene_class : metadata column of per-cell guide class (default ``"gene"``).
    mixscape_class    : metadata column of mixscape classifications.
    col               : colour for the target gene / knockout class. Controls and
                        non-perturbed cells are fixed greys, as in R.
    split_by          : metadata column to facet on, for screens spanning more
                        than one cell type.
    before_mixscape   : colour by raw guide label instead of the mixscape call.
    prtb_type         : perturbation label used by :func:`run_mixscape` (``"KO"``).
    assay             : perturbation-signature assay the scores were stored under.
    seed              : random state for the jitter strip (determinism).

    Returns
    -------
    matplotlib Figure
    """
    from scipy.stats import gaussian_kde
    plt = _mpl()

    scores = _mixscape_scores(obj, target_gene_ident, assay)
    pvec = scores["pvec"].to_numpy(dtype=float)
    cells = list(scores.index)
    guide = scores[scores.columns[1]].astype(str).to_numpy()

    nt_labels = sorted(set(guide) - {target_gene_ident})
    nt_name = nt_labels[0] if nt_labels else "NT"

    col = _r_colour(col)
    if before_mixscape:
        classes = [nt_name, target_gene_ident]
        colours = {nt_name: _r_colour("grey49"), target_gene_ident: col}
        group = guide
    else:
        if mixscape_class not in obj.meta_data.columns:
            raise KeyError(
                f"Column {mixscape_class!r} not in meta_data — run run_mixscape "
                "first, or pass before_mixscape=True."
            )
        group = obj.meta_data[mixscape_class].reindex(cells).astype(str).to_numpy()
        ko_name = f"{target_gene_ident} {prtb_type}"
        np_name = f"{target_gene_ident} NP"
        classes = [nt_name, np_name, ko_name]
        colours = {
            nt_name: _r_colour("grey49"),
            np_name: _r_colour("grey79"),
            ko_name: col,
        }

    if split_by is None:
        facets = [(None, np.arange(len(cells)))]
    else:
        if split_by not in obj.meta_data.columns:
            raise KeyError(f"split_by column {split_by!r} not found in meta_data.")
        split_vals = obj.meta_data[split_by].reindex(cells).astype(str).to_numpy()
        facets = [
            (lvl, np.where(split_vals == lvl)[0])
            for lvl in sorted(set(split_vals))
        ]

    nrow, ncol = _subplot_grid(len(facets))
    if figsize is None:
        figsize = (5.5 * ncol, 4.0 * nrow)
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize, squeeze=False)
    axes_flat = axes.ravel()
    rng = np.random.default_rng(seed)

    for ax, (level, idx) in zip(axes_flat, facets):
        sub_p, sub_g = pvec[idx], group[idx]

        # Densities first — their peak sets the band the jitter strip occupies,
        # mirroring R's read of the built panel's y-range.
        top = 0.0
        curves = []
        for cls in classes:
            vals = sub_p[sub_g == cls]
            if vals.size < 2 or np.ptp(vals) <= 0:
                continue
            grid = np.linspace(sub_p.min(), sub_p.max(), 256)
            dens = gaussian_kde(vals)(grid)
            curves.append((cls, grid, dens))
            top = max(top, float(dens.max()))
        for cls, grid, dens in curves:
            ax.plot(grid, dens, color=colours[cls], linewidth=1.5, label=cls)

        if top <= 0:
            top = 1.0
        band = top / 10.0
        for cls in classes:
            vals = sub_p[sub_g == cls]
            if vals.size == 0:
                continue
            # Controls sit above the axis, the target gene's cells below.
            if cls == nt_name:
                y = rng.uniform(0.001, band, size=vals.size)
            else:
                y = rng.uniform(-band, 0.0, size=vals.size)
            ax.scatter(vals, y, color=colours[cls], s=1.0, linewidths=0)

        ax.axhline(0.0, color="black", linewidth=0.5)
        ax.set_xlabel("perturbation score")
        ax.set_ylabel("Cell density")
        if level is not None:
            ax.set_title(str(level), fontsize=11)
        _strip_axes(ax)
        ax.legend(frameon=False, fontsize=9)

    for ax in axes_flat[len(facets):]:
        ax.axis("off")

    fig.suptitle(f"Perturbation score — {target_gene_ident}",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    return fig


def mixscape_heatmap(
    obj,
    ident_1: str,
    ident_2: Optional[str] = None,
    balanced: bool = True,
    logfc_threshold: float = 0.25,
    assay: str = "RNA",
    max_genes: int = 100,
    test_use: str = "wilcox",
    max_cells_group: Optional[int] = None,
    order_by_prob: bool = True,
    mixscape_class: str = "mixscape_class",
    prtb_type: str = "KO",
    fc_name: str = "avg_log2FC",
    pval_cutoff: float = 0.05,
    seed: int = 0,
    **kwargs,
) -> "plt.Figure":
    """DE heatmap with cells ordered by knockout probability (``MixscapeHeatmap``).

    Mirrors ``MixscapeHeatmap(object, ident.1 = "NT", ident.2 = "IFNGR2 KO",
    balanced = TRUE, max.genes = 20)``. Where :func:`plot_perturb_score` shows the
    one-dimensional score mixscape split on, this shows the genes underneath it:
    the DE genes between the two classes, with every cell ordered left-to-right by
    its knockout posterior. Read together with the class colour bar, a clean
    screen shows the expression block turning on in step with the probability —
    the escapers at the low-probability end still looking like control.

    ``ident_1`` / ``ident_2`` are ``mixscape_class`` levels (e.g. ``"NT"``,
    ``"IFNGR2 KO"``, ``"IFNGR2 NP"``), which :func:`run_mixscape` also leaves as
    the active identity.

    Parameters
    ----------
    ident_1, ident_2 : the two classes to contrast (``ident_2=None`` → all others).
    balanced         : take up to ``max_genes`` genes from *each* direction of the
                       fold change; otherwise only up-regulated ones.
    max_genes        : cap on DE genes per direction.
    max_cells_group  : downsample each class to this many cells.
    order_by_prob    : order cells by ``<mixscape_class>_p_<type>``, highest first.
                       If False, cells are shuffled (as R does).
    mixscape_class   : metadata column of mixscape classifications; also names the
                       posterior column read for the ordering.
    prtb_type        : perturbation label used by :func:`run_mixscape` (``"KO"``).
    fc_name          : fold-change column in the DE table (``"avg_log2FC"``).
    seed             : random state for downsampling / shuffling (determinism).
    **kwargs         : forwarded to :func:`do_heatmap`.

    Returns
    -------
    matplotlib Figure
    """
    from .markers import find_markers
    from .preprocessing import scale_data

    if mixscape_class not in obj.meta_data.columns:
        raise KeyError(
            f"Column {mixscape_class!r} not in meta_data — run run_mixscape first."
        )

    # find_markers resolves ident_1/ident_2 against the active identity, so drive
    # the DE off the mixscape classes — restored afterwards.
    saved_ident = pd.Categorical(list(obj.idents))
    obj.idents = list(obj.meta_data[mixscape_class].astype(str))
    try:
        markers = find_markers(
            obj,
            ident_1=ident_1,
            ident_2=ident_2,
            assay=assay,
            test_use=test_use,
            only_pos=False,
            logfc_threshold=logfc_threshold,
        )
    finally:
        obj.idents = saved_ident

    if markers.empty or fc_name not in markers.columns:
        raise ValueError(
            f"No DE genes between {ident_1!r} and {ident_2 or 'the rest'!r}; "
            "nothing to plot."
        )

    sig = markers[markers["p_val"] < pval_cutoff]
    pos = list(sig.index[sig[fc_name] > logfc_threshold])[:max_genes]
    neg = list(sig.index[sig[fc_name] < -logfc_threshold])[:max_genes] if balanced else []
    marker_list = pos + neg
    if not marker_list:
        raise ValueError(
            f"No DE genes passed p_val < {pval_cutoff} and "
            f"|{fc_name}| > {logfc_threshold}; nothing to plot."
        )

    idents = [ident_1] + ([ident_2] if ident_2 is not None else [])
    klass = obj.meta_data[mixscape_class].astype(str)
    if ident_2 is None:
        keep = list(obj.meta_data.index[klass == ident_1])
        keep += list(obj.meta_data.index[klass != ident_1])
    else:
        keep = list(obj.meta_data.index[klass.isin(idents)])
    if not keep:
        raise ValueError(f"No cells in classes {idents}.")

    sub = obj.subset(cells=keep)
    rng = np.random.default_rng(seed)

    if max_cells_group is not None:
        sub_klass = sub.meta_data[mixscape_class].astype(str)
        picked: list[str] = []
        for lvl in sorted(set(sub_klass)):
            lvl_cells = list(sub.meta_data.index[sub_klass == lvl])
            if len(lvl_cells) > max_cells_group:
                lvl_cells = list(rng.choice(lvl_cells, max_cells_group, replace=False))
            picked += lvl_cells
        sub = sub.subset(cells=picked)

    scale_data(sub, features=marker_list, assay=assay)

    if order_by_prob:
        p_col = f"{mixscape_class}_p_{prtb_type.lower()}"
        if p_col not in sub.meta_data.columns:
            raise KeyError(
                f"Posterior column {p_col!r} not in meta_data — run run_mixscape "
                f"with prtb_type={prtb_type!r}, or pass order_by_prob=False."
            )
        # NT cells carry no posterior; R initialises the column to 0, so they sort
        # to the low-probability end rather than dropping out.
        probs = sub.meta_data[p_col].to_numpy(dtype=float)
        probs = np.nan_to_num(probs, nan=0.0)
        order = np.argsort(-probs, kind="stable")
    else:
        order = rng.permutation(len(sub.cell_names()))
    ordered_cells = [sub.cell_names()[i] for i in order]

    return do_heatmap(
        sub,
        features=marker_list,
        group_by=mixscape_class,
        assay=assay,
        cells=ordered_cells,
        **kwargs,
    )


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
    "dot_plot",
    "image_dim_plot",
    "image_feature_plot",
    "spatial_dim_plot",
    "spatial_feature_plot",
    "plot_perturb_score",
    "mixscape_heatmap",
]
