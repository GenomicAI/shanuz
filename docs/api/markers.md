# Differential expression

`find_markers` implements all eight of Seurat's tests: `wilcox` (tie-corrected),
`t`, `bimod`, `LR`, `negbinom`, `roc`, `mast` and `deseq2`. Seven of them are
per-cell and reproduce Seurat's top 50 genes exactly on PBMC 3k; `deseq2` is
pseudobulk and deliberately does not, because it is answering a different
question — see [the DE vignette](../tutorials/de_vignette.md).

Two numbers to know before reading a result table:

- **`avg_log2FC` carries Seurat's pseudocount on the group *sum*, not the group
  mean.** Getting that backwards shifts every fold change and also changes which
  genes clear `logfc_threshold`, so it silently changes the returned gene set,
  not just a column.
- **`pct.1` and `pct.2` are rounded to three decimals**, by Seurat, inside
  `FindMarkers`. Anything comparing two runs gene-by-gene should not expect them
  closer than 5e-4.

## Per-cluster and per-pair tests

::: truecell.markers.find_markers

::: truecell.markers.find_all_markers

::: truecell.markers.find_conserved_markers

## Group summaries

`AggregateExpression` **sums raw counts** and is what pseudobulk differential
expression wants. `AverageExpression` **means the back-transformed values** and
is what a per-group expression summary wants. They are different functions, not
two scalings of one — see each docstring.

::: truecell.aggregate.aggregate_expression

::: truecell.aggregate.average_expression
