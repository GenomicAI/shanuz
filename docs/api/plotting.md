# Plotting

Every function here returns a matplotlib `Figure`, so you save it, compose it, or
let a notebook display it:

```python
fig = truecell.dim_plot(pbmc, reduction="umap", label=True)
fig.savefig("umap.png", dpi=150, bbox_inches="tight")
```

matplotlib and seaborn are **optional** — they live in the `analysis` extra, and
`import truecell` works without them. That is why the `Figure` return annotation is
an `if TYPE_CHECKING:` import: it documents the return type without making the
import mandatory.

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
