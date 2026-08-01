# Plotting

Every function here returns a matplotlib `Figure`, so you save it, compose it, or
let a notebook display it:

```python
fig = truecell.dim_plot(pbmc, reduction="umap", label=True)
fig.savefig("umap.png", dpi=150, bbox_inches="tight")
```

matplotlib is **optional** — it lives in the `analysis` extra, and `import truecell`
works without it. That is why the `Figure` return annotation is an
`if TYPE_CHECKING:` import: it documents the return type without making the import
mandatory.

## Theme

Every plot on this page scales its text from one base size and takes its group
colours from one palette. Set them once instead of passing the same overrides to
each call:

```python
import truecell as tc

tc.set_theme(base_size=13, style="seurat")   # bigger text, cowplot's look
fig = tc.dim_plot(pbmc)

with tc.theme_context(base_size=8):          # scoped to the block
    panel = tc.feature_plot(pbmc, ["MS4A1", "LYZ"])
```

The default (`base_size=10`, no style) leaves matplotlib's rcParams untouched and
reproduces the sizes the module used before the theme existed.

::: truecell.plotting.set_theme

::: truecell.plotting.theme_context

::: truecell.plotting.get_theme

::: truecell.plotting.reset_theme

## Colour

::: truecell.plotting.hue_pal

## Plots

::: truecell.plotting.dim_plot

::: truecell.plotting.feature_plot

::: truecell.plotting.vln_plot

::: truecell.plotting.ridge_plot

::: truecell.plotting.dot_plot

::: truecell.plotting.feature_scatter

::: truecell.plotting.elbow_plot

::: truecell.plotting.variable_feature_plot

::: truecell.plotting.viz_dim_loadings

::: truecell.plotting.dim_heatmap

::: truecell.plotting.do_heatmap

## Spatial

::: truecell.plotting.image_dim_plot

::: truecell.plotting.image_feature_plot

::: truecell.plotting.spatial_dim_plot

::: truecell.plotting.spatial_feature_plot

## Mixscape diagnostics

::: truecell.plotting.plot_perturb_score

::: truecell.plotting.mixscape_heatmap
