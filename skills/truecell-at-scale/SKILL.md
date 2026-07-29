---
name: truecell-at-scale
description: Use when a dataset is too large to analyse in RAM with truecell — leverage-score sketching (leverage_score, sketch_data, project_data) to analyse a representative subset and extend the result to every cell, and LazyMatrix (write_lazy_matrix, open_lazy_matrix, col_blocks) to keep the counts matrix memory-mapped on disk.
---

# Working at scale

Load the `truecell` skill first. Two independent tools; they compose.

- **Sketching** — reduce the number of *cells* you compute on, keeping the rare
  states a uniform sample would drop.
- **`LazyMatrix`** — keep the matrix *on disk*, so peak RAM is a function of the
  block size rather than the dataset size.

---

# 1 · Leverage-score sketching

## The pattern

```python
import truecell

# obj is normalized, with variable features selected
sketch = truecell.sketch_data(obj, ncells=50_000, method="LeverageScore")
# sketch is a NEW standalone Truecell object — rebind it

# Run the expensive analysis on the sketch only
truecell.scale_data(sketch)
truecell.run_pca(sketch, n_pcs=50)
truecell.find_neighbors(sketch, dims=range(30))
truecell.find_clusters(sketch, resolution=0.5)
truecell.run_umap(sketch, dims=range(30))

# Extend it to every cell in the full object
truecell.project_data(
    obj, sketch,
    reduction="pca", full_reduction="pca.full",
    umap_reduction="umap", full_umap_reduction="ref.umap",
    refdata={"projected.celltype": "seurat_clusters"},
)
```

## Why leverage weighting

Each cell is sampled without replacement with probability proportional to its
**statistical leverage** — how much of the data's variance it accounts for. Rare
populations are kept and indeed over-represented, where uniform sampling would
lose them. Leverage tracks rarity at **−0.929** in both truecell and Seurat.

`method="Uniform"` is the control, not a lesser option: run it to see what
leverage weighting is buying you on *your* data.

```python
scores = truecell.leverage_score(obj, features=hvg, nsketch=5000, seed=123)
```

`leverage_score` reads the **`data`** layer, not `scale.data` — matching Seurat.
Computed via a CountSketch, so no full SVD. Scores are written back to
`meta_data["leverage.score"]`.

## `sketch_data` differs from Seurat by design

Seurat stores the sketch as an extra assay on the same object. **truecell returns a
separate object**, matching its `subset` model. The sketch's active assay is
renamed to `sketched_assay` (default `"sketch"`) so the provenance is visible, and
`obj.misc["sketch"]` records how it was drawn.

## `project_data` does three things

1. **PCA** — every full-dataset cell is pushed through the sketch's loadings →
   `full.reductions[full_reduction]`.
2. **UMAP** — when `project_umap=True` and the sketch has a fitted UMAP model →
   `full.reductions[full_umap_reduction]`.
3. **Labels** — when `refdata=` is given: a weighted k-nearest-neighbour vote
   *inside the projected reduction*, with the sketch's own rows as reference
   (Seurat's `TransferSketchLabels`). Written to `full.meta_data`.

Step 3 is deliberately **not** the anchor path. Finding anchors between the
sketch and the full dataset costs exactly what sketching exists to avoid — on a
million-cell object the anchor route is unusable, not merely different. The two
agree per-cell 98.1 % of the time at matching accuracy.

`project_data` takes no `seed`: the label vote is deterministic.

## Sizing

- `ncells` — the sketch. Big enough that the rarest population you care about
  survives with enough cells to cluster. Sanity-check afterwards:
  `sketch.meta_data["known_type"].value_counts()`.
- `nsketch=5000` — the CountSketch dimension used to *compute* leverage. Separate
  knob; the default is fine unless leverage looks degenerate.

Measured against Seurat: exact-regime Spearman **1.000000**; projected label
accuracy 0.9050 on both sides.

**Cautionary note from this feature's own history.** `project_data` scored
*above* Seurat while it was broken. A divergence that flatters the port is a
reason to look harder, not to keep it.

---

# 2 · `LazyMatrix` — out-of-core matrices

BPCells-style, built on NumPy memory-mapping. **No new dependency.**

## The pattern

```python
from truecell import write_lazy_matrix, open_lazy_matrix, is_lazy

assay = obj.get_assay()

write_lazy_matrix(assay.layers["counts"], "counts.mat")   # persist to a directory
lazy = open_lazy_matrix("counts.mat")                     # maps, does not read
assay.set_layer_data("counts", lazy)                      # a LazyMatrix is a valid layer

# Slicing reads only the touched cells' non-zeros off disk, returning scipy sparse
block = assay.layer_data("counts", cells=obj.cell_names()[:1000])

# Reductions stream in one pass, without materialising
per_cell = lazy.sum(axis=0)
per_gene = lazy.mean(axis=1)
nnz_col  = lazy.nnz_per_col()

# The streaming primitive — a million cells at bounded RAM
for start, stop, chunk in lazy.col_blocks(block_size=50_000):
    ...        # chunk is a csc_matrix of those cells
```

Attributes: `shape`, `nrow`, `ncol`, `nnz`, `dtype`, `path`, `ndim`.
Methods: `sum`, `mean`, `nnz_per_col`, `nnz_per_row`, `col_blocks`, `to_scipy`,
`toarray`, `close`.

## Why compressed-sparse-**column**

truecell matrices are features × cells, and the operations that dominate at scale
— sketching, cell subsetting, per-cell normalisation — select *cells*, i.e.
columns. CSC stores each column contiguously, so reading an arbitrary set of
cells costs only their own non-zeros.

Selecting *features* (rows) still scans the touched columns, which is why
`m[np.ix_(genes, cells)]` narrows to cells first and applies the gene filter to a
small in-memory block.

## Rules

- **Build the object on the lazy layer, or swap it in — either works**, but check
  it stayed lazy: `is_lazy(obj.get_assay().layers["counts"])`. A constructor that
  densified on the way in was a real defect; the natural path is tested now, but
  verify rather than assume.
- **`np.asarray(lazy)` / `as_dense(lazy)` materialises the whole matrix.** That
  is the deliberate escape hatch — fine on a small dataset, fatal on the
  million-cell path. If a step is unexpectedly slow, check whether something
  densified.
- **Measure peak memory, don't infer it.** Going on disk once *raised* peak RAM
  4.6× because five functions densified the whole store. The bug was found only
  because someone measured and the number came out backwards.

truecell's on-disk and in-memory paths are **bit-identical**. (Seurat's differ by
1.0e-06 and pick a different variable feature — the comparison worth knowing if
you are reconciling the two.)

---

## Composing them

The intended million-cell shape: counts on disk as a `LazyMatrix`, a
leverage-weighted sketch pulled into memory for the expensive work, results
projected back to every cell.

```python
lazy = open_lazy_matrix("counts.mat")
obj.get_assay().set_layer_data("counts", lazy)
truecell.normalize_data(obj)
truecell.find_variable_features(obj)

sketch = truecell.sketch_data(obj, ncells=50_000)
# ... full analysis on `sketch` ...
truecell.project_data(obj, sketch, refdata={"celltype": "seurat_clusters"})
```

## Reference

- [Leverage-score sketching](https://genomicai.github.io/truecell/tutorials/sketch_vignette/) — both of Seurat's regimes, with the uniform control.
- [Out of core: `LazyMatrix` vs BPCells](https://genomicai.github.io/truecell/tutorials/lazy_vignette/) — 14 of 14 anchors, and the memory measurements.
