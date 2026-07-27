---
name: shanuz-plotting
description: Use when producing figures from a shanuz object — dim_plot, feature_plot, vln_plot, dot_plot, ridge_plot, elbow_plot, do_heatmap, dim_heatmap, feature_scatter, variable_feature_plot, viz_dim_loadings, the spatial image_* / spatial_* plots, and the Mixscape diagnostics. Covers what each maps to in Seurat, the shared arguments, and saving figures headlessly.
---

# Plotting in shanuz

Load the `shanuz` skill first.

**Every plotting function returns a `matplotlib.figure.Figure`** and draws
nothing on its own. You save or display it.

```python
import shanuz

fig = shanuz.dim_plot(obj, reduction="umap", label=True)
fig.savefig("umap.png", dpi=150, bbox_inches="tight")
```

Requires `pip install "shanuz[analysis]"` (matplotlib + seaborn). matplotlib is
an optional dependency, imported lazily — which is why the return annotation is
`Figure` under `if TYPE_CHECKING` and not a hard import.

In a script or CI, pick a non-interactive backend **before** importing pyplot:

```python
import matplotlib
matplotlib.use("Agg")
```

## The map from Seurat

| Seurat | shanuz | Shows |
|---|---|---|
| `DimPlot` | `dim_plot` | Cells on an embedding, coloured by identity |
| `FeaturePlot` | `feature_plot` | Expression of one or more features on an embedding |
| `VlnPlot` | `vln_plot` | Violin per group |
| `DotPlot` | `dot_plot` | Mean expression + fraction detecting, per group × feature |
| `RidgePlot` | `ridge_plot` | Ridgeline per group |
| `ElbowPlot` | `elbow_plot` | Standard deviation per PC |
| `FeatureScatter` | `feature_scatter` | Two features against each other |
| `VariableFeaturePlot` | `variable_feature_plot` | Mean–variance, HVGs highlighted |
| `VizDimLoadings` | `viz_dim_loadings` | Top loading genes per component |
| `DimHeatmap` | `dim_heatmap` | Cells × top loading genes, per component |
| `DoHeatmap` | `do_heatmap` | Expression heatmap, cells ordered by group |
| `ImageDimPlot` | `image_dim_plot` | Spatial cells coloured by identity |
| `ImageFeaturePlot` | `image_feature_plot` | Spatial cells coloured by expression |
| `SpatialDimPlot` | `spatial_dim_plot` | Visium spots over the H&E, by identity |
| `SpatialFeaturePlot` | `spatial_feature_plot` | Visium spots over the H&E, by expression |
| `PlotPerturbScore` | `plot_perturb_score` | Mixscape perturbation-score densities |
| `MixscapeHeatmap` | `mixscape_heatmap` | Mixscape DE genes, cells ordered by KO probability |

## Shared arguments

- `group_by=` — a metadata column. `None` uses the active identity.
- `assay=` / `layer=` — which assay and layer to read features from. This is how
  you plot protein instead of RNA: `feature_plot(obj, ["CD3"], assay="ADT")`.
  There is no `DefaultAssay(obj) <- "ADT"` step.
- `reduction=` — `"umap"`, `"pca"`, `"harmony"`, `"wnn_umap"`, …
- `ncol=` — panels per row when several features are given.
- `figsize=` — matplotlib inches.
- `palette=` / `cols=` — colours.
- `pt_size=`, `alpha=` — point size and opacity.

## The common calls

```python
# Embeddings
fig = shanuz.dim_plot(obj, reduction="umap", label=True, label_size=9, pt_size=4.0)
fig = shanuz.dim_plot(obj, reduction="umap", group_by="batch", label=False)

# Expression
fig = shanuz.feature_plot(obj, ["LYZ", "MS4A1", "NKG7"], reduction="umap", ncol=3,
                          min_cutoff="q05", max_cutoff="q95", colormap="YlOrRd")
fig = shanuz.vln_plot(obj, ["nFeature_RNA", "nCount_RNA", "percent.mt"], ncol=3)
fig = shanuz.ridge_plot(obj, ["CD3D", "LYZ"])
fig = shanuz.dot_plot(obj, canonical_markers, group_by="cell_type",
                      col_min=-2.5, col_max=2.5, dot_scale=6.0, scale=True)

# QC and dimensionality
fig = shanuz.feature_scatter(obj, "nCount_RNA", "percent.mt")
fig = shanuz.variable_feature_plot(obj, n_label=10)
fig = shanuz.elbow_plot(obj, ndims=50)
fig = shanuz.viz_dim_loadings(obj, dims=[1, 2], n_features=15)
fig = shanuz.dim_heatmap(obj, dims=list(range(9)), cells=500, balanced=True)

# Markers
top = markers.groupby("cluster").head(10)
fig = shanuz.do_heatmap(obj, list(top["gene"]), layer="scale.data")
```

`min_cutoff` / `max_cutoff` on `feature_plot` accept a number **or** a quantile
string (`"q05"`, `"q95"`) — the usual fix when one outlier cell flattens the
whole colour scale.

`do_heatmap` reads `scale.data` by default, so the genes you pass must have been
scaled. If your heatmap is empty, that is usually why: `scale_data` defaults to
the variable features, and your marker list may reach outside them.

```python
shanuz.scale_data(obj, features=list(set(shanuz.generics.variable_features(obj)) | set(genes)))
```

## Spatial and Mixscape

See `shanuz-spatial` and `shanuz-multimodal` for these in context.

```python
fig = shanuz.image_dim_plot(obj, group_by="cell_type", size=1.0, flip_y=True)
fig = shanuz.image_feature_plot(obj, feature="Slc17a7", cmap="viridis")
fig = shanuz.spatial_dim_plot(obj, group_by="seurat_clusters", pt_size_factor=1.6,
                              image_alpha=1.0, crop=True)
fig = shanuz.spatial_feature_plot(obj, feature="Hpca")

fig = shanuz.plot_perturb_score(obj, target_gene_ident="IFNGR2", assay="PRTB")
fig = shanuz.mixscape_heatmap(obj, ident_1="IFNGR2 KO", ident_2="NT", max_genes=100)
```

## Customising

The returned `Figure` is an ordinary matplotlib figure — reach into it.

```python
fig = shanuz.dim_plot(obj, reduction="umap")
ax = fig.axes[0]
ax.set_title("PBMC 3k — Louvain, res 0.5")
ax.set_xlabel("UMAP 1")
fig.savefig("umap.pdf", bbox_inches="tight")     # vector output for figures
```

Close figures in a loop or the process will hold every one of them:

```python
import matplotlib.pyplot as plt
for gene in genes:
    fig = shanuz.feature_plot(obj, gene)
    fig.savefig(f"{gene}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
```

## Traps

| Symptom | Cause |
|---|---|
| Nothing displays | The function returns a figure; it does not call `plt.show()`. |
| `ModuleNotFoundError: matplotlib` | Install `shanuz[analysis]`. |
| Hangs or errors in CI / headless | Set `matplotlib.use("Agg")` before importing pyplot. |
| `do_heatmap` blank | Genes not in `scale.data` — see above. |
| Colour scale washed out | One outlier cell; use `min_cutoff="q05"`, `max_cutoff="q95"`. |
| Plotting protein shows RNA | Pass `assay="ADT"`; there is no default-assay switch. |
| Spatial plot mirrored | Toggle `flip_y`. |
| Memory grows across a loop | Figures never closed. |
