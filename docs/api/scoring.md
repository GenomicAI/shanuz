# Signature scoring

Score a gene programme per cell against a background of expression-matched
control genes. `cell_cycle_scoring` is `add_module_score` run twice, on the S and
G2/M lists, plus the discrete phase call.

The control genes are drawn at random from expression bins, so these scores carry
an RNG. Against R Seurat the per-cell phase call is 96.6 % concordant and the
continuous scores correlate at Pearson ≥ 0.998 — the residual is the control
draw, and nothing else. [Cell-cycle vignette](../tutorials/cellcycle_vignette.md).

::: truecell.module_score.add_module_score

::: truecell.module_score.cell_cycle_scoring

## The bundled gene lists

::: truecell.module_score.CC_GENES
