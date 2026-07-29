# Spatial

Imaging-based and spot-based assays: Xenium, Visium, CosMx and MERSCOPE. The
container mirrors Seurat v5's — an `FOV` holding `Centroids`, `Segmentation`
and `Molecules` boundaries, or a `VisiumV2` holding the H&E image and its
`ScaleFactors`.

!!! warning "`radius` on a Visium FOV"
    Seurat stores `spot_diameter_fullres` in the FOV's `radius` slot — a
    diameter where a radius is named — and `Radius()` on its own `VisiumV2`
    returns `NULL`. `load_visium` stores a radius. The slide's fixed 100 µm spot
    pitch is what settles which reading is right; the working is in
    [the Visium vignette](../tutorials/visium_vignette.md).

`find_spatially_variable_features` computes Moran's I on R's inverse-square
distance weights, not on a kNN graph. Those give different answers, and the kNN
version was the bug.

## Loading

::: truecell.spatial.loaders.load_xenium

::: truecell.spatial.loaders.load_visium

::: truecell.spatial.loaders.load_cosmx

::: truecell.spatial.loaders.load_merscope

## Containers

::: truecell.spatial.fov.FOV

::: truecell.spatial.fov.create_fov

::: truecell.spatial.fov.create_fovs

::: truecell.spatial.centroids.Centroids

::: truecell.spatial.centroids.create_centroids

::: truecell.spatial.segmentation.Segmentation

::: truecell.spatial.segmentation.create_segmentation

::: truecell.spatial.molecules.Molecules

::: truecell.spatial.molecules.create_molecules

::: truecell.spatial.base.SpatialImage

::: truecell.spatial.visium.VisiumV2

::: truecell.spatial.visium.ScaleFactors

## Spatial analysis

::: truecell.spatial.analysis.get_tissue_coordinates

::: truecell.spatial.analysis.spatial_knn

::: truecell.spatial.analysis.nearest_neighbor_distance

::: truecell.spatial.analysis.local_neighborhood

::: truecell.spatial.analysis.build_niche_assay

::: truecell.spatial.variable_features.find_spatially_variable_features

::: truecell.composition.composition_test
