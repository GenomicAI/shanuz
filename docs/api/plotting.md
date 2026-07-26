# Plotting

Every function here returns a matplotlib `Figure`, so you save it, compose it, or
let a notebook display it:

```python
fig = shanuz.dim_plot(pbmc, reduction="umap", label=True)
fig.savefig("umap.png", dpi=150, bbox_inches="tight")
```

matplotlib and seaborn are **optional** — they live in the `analysis` extra, and
`import shanuz` works without them. That is why the `Figure` return annotation is
an `if TYPE_CHECKING:` import: it documents the return type without making the
import mandatory.

::: shanuz.plotting.dim_plot

::: shanuz.plotting.feature_plot

::: shanuz.plotting.vln_plot

::: shanuz.plotting.ridge_plot

::: shanuz.plotting.dot_plot

::: shanuz.plotting.feature_scatter

::: shanuz.plotting.elbow_plot

::: shanuz.plotting.variable_feature_plot

::: shanuz.plotting.viz_dim_loadings

::: shanuz.plotting.dim_heatmap

::: shanuz.plotting.do_heatmap

## Spatial

::: shanuz.plotting.image_dim_plot

::: shanuz.plotting.image_feature_plot

::: shanuz.plotting.spatial_dim_plot

::: shanuz.plotting.spatial_feature_plot

## Mixscape diagnostics

::: shanuz.plotting.plot_perturb_score

::: shanuz.plotting.mixscape_heatmap
