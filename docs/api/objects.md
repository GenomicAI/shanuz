# Objects

The container, before any analysis touches it. `Truecell` holds one or more assays,
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

::: truecell.truecell.Truecell

::: truecell.truecell.create_truecell_object

## Assays

::: truecell.assay5.Assay5

::: truecell.assay5.create_assay5_object

::: truecell.assay5.StdAssay

::: truecell.assay.Assay

::: truecell.assay.create_assay_object

## Reductions, graphs and neighbours

::: truecell.dimreduc.DimReduc

::: truecell.graph.Graph

::: truecell.graph.as_graph

::: truecell.neighbor.Neighbor

## Supporting structures

::: truecell.jackstraw.JackStrawData

::: truecell.logmap.LogMap

::: truecell.mixins.key_mixin.KeyMixin

::: truecell.command.TruecellCommand

::: truecell.command.log_truecell_command
