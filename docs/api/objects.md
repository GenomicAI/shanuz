# Objects

The container, before any analysis touches it. `Shanuz` holds one or more assays,
the reductions computed off them, the neighbour graphs, the per-cell metadata and
the command log — the same slots R's `Seurat` S4 class holds, as ordinary Python
classes with `__slots__`.

The one structural difference worth knowing up front: R's `Assay5` inherits from
`dgCMatrix`; here `Assay5` *wraps* a SciPy CSC matrix rather than subclassing it,
because subclassing `scipy.sparse` is a well-known trap. Everything you reach
through [the generics](generics.md) behaves the same either way.

The object model is checked against Seurat anchor by anchor in
[The Object Model Itself](../tutorials/objects_vignette.md) — 91 of 91 exact, no
tolerance.

## The top-level object

::: shanuz.shanuz.Shanuz

::: shanuz.shanuz.create_shanuz_object

## Assays

::: shanuz.assay5.Assay5

::: shanuz.assay5.create_assay5_object

::: shanuz.assay5.StdAssay

::: shanuz.assay.Assay

::: shanuz.assay.create_assay_object

## Reductions, graphs and neighbours

::: shanuz.dimreduc.DimReduc

::: shanuz.graph.Graph

::: shanuz.graph.as_graph

::: shanuz.neighbor.Neighbor

## Supporting structures

::: shanuz.jackstraw.JackStrawData

::: shanuz.logmap.LogMap

::: shanuz.mixins.key_mixin.KeyMixin

::: shanuz.command.ShanuzCommand

::: shanuz.command.log_shanuz_command
