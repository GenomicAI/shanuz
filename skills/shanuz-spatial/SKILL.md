---
name: shanuz-spatial
description: Use for spatial transcriptomics in shanuz — loading Xenium, Visium, CosMx or MERSCOPE data, the FOV/Centroids/Segmentation/VisiumV2 containers, tissue coordinates, nearest-neighbour and local-density statistics, niche assays, spatially variable features (Moran's I / mark variogram), composition tests, and the image_* / spatial_* plots.
---

# Spatial analysis in shanuz

Load the `shanuz` skill first. A spatial object is an ordinary `Shanuz` object
with one extra container: `obj.images`, holding an `FOV` (imaging-based) or a
`VisiumV2` (spot-based).

## Loading

```python
import shanuz

obj = shanuz.load_xenium(path, assay="Xenium", keep_controls=False)
obj = shanuz.load_visium(path, assay="Spatial", image=True,
                         image_resolution="lowres", filter_by_tissue=True)
obj = shanuz.load_cosmx(path, assay="Nanostring", fov_column="fov")
obj = shanuz.load_merscope(path, assay="Vizgen", keep_controls=False)
```

`keep_controls=False` drops the negative-control / blank probes — keep them only
when auditing the run's background.

Test data:

```python
from shanuz.datasets import xenium_mouse_brain, visium_mouse_brain
obj = shanuz.load_xenium(xenium_mouse_brain())     # returns a Path
obj = shanuz.load_visium(visium_mouse_brain())
```

Everything downstream — normalize, PCA, cluster — is the standard workflow. The
spatial functions are what you add on top.

## The containers

| Class | Holds |
|---|---|
| `FOV` | One field of view: boundaries + molecules for a set of cells |
| `Centroids` | Cell centre points, with `radius`, `nsides`, `theta` |
| `Segmentation` | Cell boundary polygons |
| `Molecules` | Individual transcript coordinates |
| `VisiumV2` | Spot coordinates **plus** the H&E image and its `ScaleFactors` |
| `ScaleFactors` | The `scalefactors_json.json` values |

```python
obj.image_names()
coords = shanuz.get_tissue_coordinates(obj, image=None)     # DataFrame: x, y, cell
img    = obj.images["slice1"]
shanuz.generics.radius(img)
shanuz.generics.get_image(img)                              # VisiumV2 only
```

Build them by hand when your coordinates come from somewhere else:

```python
fov = shanuz.create_fov(coords_df, type_="centroids", radius=10.0, assay="Xenium")
cen = shanuz.create_centroids(coords_df, nsides=0, radius=10.0)
seg = shanuz.create_segmentation(polygon_df)
fovs = shanuz.create_fovs(coords_df, fov="fov_column")      # one FOV per level
```

`coords` frames are `x`, `y` (plus `cell`), indexed or columned to match the
object's barcodes.

### The Visium radius, which is a real divergence

`scalefactors_json.json`'s `spot_diameter_fullres` is a **diameter**, and Seurat
stores it in a slot named `radius`. Reading it as a radius makes each spot
overlap its neighbour by 31 µm on a slide with a fixed 100 µm pitch — the
geometry settles it. **shanuz's `radius()` returns half of Seurat's, on purpose**,
and a test pins that divergence so it cannot be "fixed" back.

Also: R's `Radius()` returns `NULL` on a `VisiumV2` (there is no
`Radius.VisiumV2` method, only `VisiumV1`). shanuz's `VisiumV2.radius()` answers.
This is the one tutorial where Seurat is the one that is wrong.

## Spatial statistics

```python
# How far is each cell of type X from the nearest cell of type Y?
d = shanuz.nearest_neighbor_distance(obj, group_by="cell_type", reference="Astrocyte")

# How many of each type sit in each cell's k-neighbourhood?
n = shanuz.local_neighborhood(obj, group_by="cell_type", k=10)

# Is a type over- or under-represented between conditions?
c = shanuz.composition_test(obj, group_by="cell_type", split_by="condition")
```

All three return DataFrames and leave the object alone.

## Niches

```python
shanuz.build_niche_assay(obj, group_by="cell_type", k=20, niches=4,
                         assay_name="niche", cluster=True, seed=0)
fig = shanuz.image_dim_plot(obj, group_by="niches")
```

Each cell is described by the composition of its `k` nearest neighbours; those
profiles are clustered into `niches` tissue neighbourhoods. `k` sets the spatial
scale — a small `k` finds tight micro-environments, a large one finds regions.
Sweep it rather than accepting the default; `niches` likewise.

## Spatially variable features

```python
svf = shanuz.find_spatially_variable_features(
    obj, features=shanuz.generics.variable_features(obj),
    method="moransi", weights="inverse_square",
)
top = svf.head(20)          # rank 1 is the most spatially structured
```

**`weights="inverse_square"` is the default and reproduces R exactly** — 1/d²
between every pair of cells, row-standardised. It is O(n²). `weights="knn"` is an
approximation trading that for O(nk); reach for it only when the slide is too
large to wait on, and say which you used.

`method="markvariogram"` is the alternative. Its `r_metric=5.0` is **in units of
the median nearest-neighbour spacing** — the default means "five cells apart" —
and `bandwidth` is the half-width of the band around it. Widen `bandwidth` on a
sparse slide where the band catches too few pairs.

Pass `features=` (the variable features are a good default) — running on a full
panel is much slower and rarely tells you more.

Moran's I matches R to **1.6e-14**, and picks 10/10 of Seurat's top genes, on a
36,602-cell slide R cannot hold in memory.

## Plotting

```python
fig = shanuz.image_dim_plot(obj, group_by="cell_type", size=1.0, flip_y=True)
fig = shanuz.image_feature_plot(obj, feature="Slc17a7", cmap="viridis")

# Visium — spots drawn over the H&E at their true diameter
fig = shanuz.spatial_dim_plot(obj, group_by="seurat_clusters", pt_size_factor=1.6)
fig = shanuz.spatial_feature_plot(obj, feature="Hpca", image_alpha=1.0, crop=True)
```

`flip_y=True` is the default because image coordinates run top-down; flip it back
if your slide comes out mirrored against the source viewer.

Use `image_*` for imaging-based platforms (Xenium/CosMx/MERSCOPE) and
`spatial_*` for Visium, where there is a tissue image to draw under the spots.

## Traps

| Symptom | Cause |
|---|---|
| `obj.images` empty after loading | Wrong path level — point at the directory holding the matrix and `spatial/`, not one above. |
| Plot is mirrored vertically | `flip_y`; toggle it. |
| `find_spatially_variable_features` never returns | `weights="inverse_square"` is O(n²) over all cells. Restrict `features=`, or use `weights="knn"` and note the approximation. |
| Spots overlap on the H&E | Radius vs diameter — see above. shanuz is deliberately half of Seurat. |
| Control probes dominate the variable features | Loaded with `keep_controls=True`. |
| Niches look like cell types | `k` too small; the niche assay degenerates to the cell's own identity. |

## Reference

- [Xenium spatial](https://genomicai.github.io/shanuz/tutorials/xenium_spatial_tutorial/) — verified to 8 s.f. against R Seurat.
- [Spatial statistics & the container](https://genomicai.github.io/shanuz/tutorials/svf_vignette/) — 38 of 39 anchors exact.
- [Visium](https://genomicai.github.io/shanuz/tutorials/visium_vignette/) — 24 of 24 anchors, and the radius finding.
