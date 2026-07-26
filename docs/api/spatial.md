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

::: shanuz.spatial.loaders.load_xenium

::: shanuz.spatial.loaders.load_visium

::: shanuz.spatial.loaders.load_cosmx

::: shanuz.spatial.loaders.load_merscope

## Containers

::: shanuz.spatial.fov.FOV

::: shanuz.spatial.fov.create_fov

::: shanuz.spatial.fov.create_fovs

::: shanuz.spatial.centroids.Centroids

::: shanuz.spatial.centroids.create_centroids

::: shanuz.spatial.segmentation.Segmentation

::: shanuz.spatial.segmentation.create_segmentation

::: shanuz.spatial.molecules.Molecules

::: shanuz.spatial.molecules.create_molecules

::: shanuz.spatial.base.SpatialImage

::: shanuz.spatial.visium.VisiumV2

::: shanuz.spatial.visium.ScaleFactors

## Spatial analysis

::: shanuz.spatial.analysis.get_tissue_coordinates

::: shanuz.spatial.analysis.spatial_knn

::: shanuz.spatial.analysis.nearest_neighbor_distance

::: shanuz.spatial.analysis.local_neighborhood

::: shanuz.spatial.analysis.build_niche_assay

::: shanuz.spatial.variable_features.find_spatially_variable_features

::: shanuz.composition.composition_test
