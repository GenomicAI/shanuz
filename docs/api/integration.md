# Integration and reference mapping

Two different jobs that share machinery. **Integration** removes a batch effect
between datasets you intend to analyse together. **Reference mapping** leaves the
reference untouched and projects a query into it, carrying labels across.

`integrate_layers` is the v5 dispatcher — `method="harmony" | "cca" | "rpca"` —
and runs `integrate_embeddings`, the embedding-space algorithm, as Seurat v5
does. The v4 pair (`find_integration_anchors` + `integrate_data`, which corrects
expression rather than embeddings) is still available directly and is still what
you want if you are reproducing a v4 analysis. These were the same function once,
which was a bug: the v5 name ran the v4 algorithm.

Both anchor paths are compared against Seurat's own anchors, not just against the
clustering they produce, in [Anchor internals](../tutorials/anchors_vignette.md).

## Batch correction

::: truecell.integration.integrate_layers

::: truecell.integration.run_harmony

## Anchors, directly

::: truecell.anchors.find_integration_anchors

::: truecell.anchors.integrate_embeddings

::: truecell.anchors.integrate_data

::: truecell.anchors.IntegrationAnchors

## Reference mapping

::: truecell.transfer.find_transfer_anchors

::: truecell.transfer.transfer_data

::: truecell.transfer.TransferAnchors

::: truecell.mapping.map_query

::: truecell.mapping.project_umap
