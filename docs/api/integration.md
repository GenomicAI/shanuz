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

::: shanuz.integration.integrate_layers

::: shanuz.integration.run_harmony

## Anchors, directly

::: shanuz.anchors.find_integration_anchors

::: shanuz.anchors.integrate_embeddings

::: shanuz.anchors.integrate_data

::: shanuz.anchors.IntegrationAnchors

## Reference mapping

::: shanuz.transfer.find_transfer_anchors

::: shanuz.transfer.transfer_data

::: shanuz.transfer.TransferAnchors

::: shanuz.mapping.map_query

::: shanuz.mapping.project_umap
