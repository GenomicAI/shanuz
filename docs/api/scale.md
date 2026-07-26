# Working at scale

Two independent answers to a dataset that will not fit: analyse a representative
subset, or keep the matrix on disk.

**Sketching** draws a leverage-weighted subset — rare states kept rather than
sampled away — analyses that, then extends the result back to every cell.
`leverage_score` gets the per-cell scores via a CountSketch, without a full SVD.

**`LazyMatrix`** is the on-disk path, memory-mapped compressed-sparse-column
arrays in BPCells' spirit but with no new dependency. A slice reads only the
cells it touches, `col_blocks` streams a million cells at bounded RAM, and it
drops straight into an `Assay5` layer. Against BPCells on PBMC 3k, shanuz's
on-disk and in-memory paths are bit-identical to each other; Seurat's differ by
1.0e-06. [The comparison](../tutorials/lazy_vignette.md).

## Sketching

::: shanuz.sketch.leverage_score

::: shanuz.sketch.sketch_data

::: shanuz.sketch.project_data

## Out-of-core matrices

::: shanuz.lazy.LazyMatrix

::: shanuz.lazy.write_lazy_matrix

::: shanuz.lazy.open_lazy_matrix

::: shanuz.lazy.is_lazy
