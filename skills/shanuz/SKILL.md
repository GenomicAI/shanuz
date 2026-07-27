---
name: shanuz
description: Use when writing, reading, debugging or reviewing Python single-cell analysis code that uses the shanuz package (a port of R Seurat) — creating Shanuz objects, the QC → normalize → HVG → scale → PCA → neighbours → clusters → UMAP → markers pipeline, or translating Seurat code to Python. Start here; it carries the API contracts that break code silently and routes to the task-specific shanuz skills.
---

# shanuz

`shanuz` is a Python port of [Seurat](https://satijalab.org/seurat/) v5 — the same
data structures and the same algorithms, checked against R Seurat 5.5.1 test by
test. Version **0.9.0**, Python **3.12+**, MIT.

- Docs: <https://genomicai.github.io/shanuz/> · Repo: <https://github.com/GenomicAI/shanuz>
- 105 public names, all exported from the package root **except the generics**
  (see contract 6 below).

```bash
pip install shanuz                 # core: objects, preprocessing, PCA, markers
pip install "shanuz[analysis]"     # + clustering, UMAP, plotting  ← the usual one
pip install "shanuz[anndata]"      # + AnnData interop
pip install "shanuz[integration]"  # + Harmony (harmonypy)
pip install "shanuz[deseq2]"       # + pseudobulk DESeq2
pip install "shanuz[all]"          # everything, incl. dev + docs tooling
```

## The six contracts

Almost every mistake made against this API is one of these. Check them before
writing anything.

**1. Analysis functions mutate in place and return `None`.**

```python
shanuz.normalize_data(pbmc)          # correct — call for effect
pbmc = shanuz.normalize_data(pbmc)   # WRONG — pbmc is now None
```

Applies to `normalize_data`, `find_variable_features`, `scale_data`,
`percentage_feature_set`, `run_pca` / `run_ica` / `run_spca` / `run_tsne` /
`glm_pca`, `find_neighbors`, `find_multi_modal_neighbors`, `find_clusters`,
`run_umap`, `run_harmony`, `integrate_layers`.

Rebind **only** for the functions that build a new object:
`subset`, `merge`, `sketch_data`, `integrate_data`.

A third group mutates in place *and* returns the same object it mutated
(`sctransform`, `add_module_score`, `cell_cycle_scoring`, `hto_demux`,
`multiseq_demux`, `calc_perturb_sig`, `run_mixscape`, `mixscape_lda`). Rebinding
is harmless there but is not the idiom — call for effect everywhere except the
four above.

Functions that compute a *result* rather than mutating return it: `find_markers`,
`find_all_markers`, `find_conserved_markers`, `aggregate_expression`,
`transfer_data`, `find_spatially_variable_features`, `composition_test`,
`leverage_score`, and every plotting function (→ `matplotlib.figure.Figure`).

**2. `dims` is 0-based.** `range(10)` here is `1:10` in R. This is the one
indexing difference in the API and it follows Python on purpose.

**3. Matrices are features × cells** — genes as rows, same as Seurat, transposed
relative to AnnData/scanpy. `create_shanuz_object(counts, ...)` expects
genes × cells.

**4. Expression lives in named layers**, not attributes: `counts` (raw),
`data` (log-normalized), `scale.data` (z-scored). Read them with
`obj.get_assay().layer_data("data")`.

**5. Names are snake_case ports of Seurat's**: `FindMarkers` → `find_markers`,
`RunPCA` → `run_pca`, `nn.method` → `nn_method`. Parameters that collide with
Python keywords get a trailing underscore (`lambda_`, `type_`).

**6. The generics are not top-level.** `cells`, `features`, `idents`,
`fetch_data`, `layer_data`, `layers`, `split_layers`, `join_layers`,
`embeddings`, `loadings`, `stdev`, `variable_features`, `which_cells`,
`rename_idents` and the rest live in `shanuz.generics`:

```python
import shanuz
shanuz.generics.features(pbmc)     # correct
shanuz.features(pbmc)              # AttributeError
```

Same for the loaders: `from shanuz.io import read_10x`,
`from shanuz.datasets import pbmc3k`, `from shanuz.compat.anndata import as_anndata`.
(The published API reference says everything is top-level; for the generics page
that is not true.)

## The canonical pipeline

```python
import shanuz
from shanuz.datasets import pbmc3k

counts, genes, cells = pbmc3k()            # caches to ~/.shanuz_data/ (~24 MB)
pbmc = shanuz.create_shanuz_object(
    counts=counts, feature_names=genes, cell_names=cells,
    project="pbmc3k", min_cells=3, min_features=200,
)

# QC
shanuz.percentage_feature_set(pbmc, pattern=r"^MT-", col_name="percent.mt")
md = pbmc.meta_data
keep = (md["nFeature_RNA"] > 200) & (md["nFeature_RNA"] < 2500) & (md["percent.mt"] < 5)
pbmc = pbmc.subset(cells=list(md.index[keep]))     # subset RETURNS a new object

# Normalize → select → scale
shanuz.normalize_data(pbmc, normalization_method="LogNormalize", scale_factor=10000)
shanuz.find_variable_features(pbmc, selection_method="vst", nfeatures=2000)
shanuz.scale_data(pbmc)                    # defaults to the variable features

# Reduce → graph → cluster → embed
shanuz.run_pca(pbmc, n_pcs=50)
shanuz.find_neighbors(pbmc, dims=range(10), k_param=20)   # writes RNA_nn, RNA_snn
shanuz.find_clusters(pbmc, resolution=0.5)                # writes seurat_clusters + idents
shanuz.run_umap(pbmc, dims=range(10), seed=42)

# Markers
markers = shanuz.find_all_markers(pbmc, only_pos=True, min_pct=0.25, logfc_threshold=0.25)

fig = shanuz.dim_plot(pbmc, reduction="umap", label=True)
fig.savefig("umap.png", dpi=150, bbox_inches="tight")
```

Where results land: `pbmc.meta_data` (per-cell columns), `pbmc.reductions`
(`"pca"`, `"umap"`), `pbmc.graphs` (`"RNA_nn"`, `"RNA_snn"`), `pbmc.idents`,
`pbmc.commands` (the audit log), `pbmc.misc` (stashed fit details).

## Which skill to load

| Task | Skill |
|---|---|
| Standard scRNA-seq run, QC thresholds, how many PCs, resolution choice | `shanuz-workflow` |
| Marker genes, the eight DE tests, pseudobulk, conserved markers | `shanuz-differential-expression` |
| Batch correction (Harmony/CCA/RPCA), label transfer, reference mapping | `shanuz-integration` |
| CITE-seq / WNN, cell hashing demultiplexing, pooled CRISPR (Mixscape) | `shanuz-multimodal` |
| Xenium / Visium / CosMx / MERSCOPE, niches, spatially variable features | `shanuz-spatial` |
| Datasets too big for RAM — leverage sketching, on-disk `LazyMatrix` | `shanuz-at-scale` |
| Any figure | `shanuz-plotting` |
| Porting existing R Seurat code, or comparing the two tools' output | `shanuz-from-seurat` |
| Contributing to shanuz itself — tests, lint, docs, fidelity method, release | `shanuz-dev` |

## Bundled reference

- [`reference/api-map.md`](reference/api-map.md) — every public function with its
  real signature and its Seurat equivalent. Read this before guessing a parameter name.
- [`reference/object-model.md`](reference/object-model.md) — the `Shanuz` /
  `Assay5` / `DimReduc` / `Graph` containers, the generics, subsetting, layers,
  AnnData interop.

## Known, deliberate differences from Seurat

Not bugs; do not "fix" them, and do not report them as regressions.

- **Louvain cluster counts drift by one.** Same algorithm, same resolution,
  different local optimum. PBMC 3k: 8 clusters to Seurat's 9 at ARI 0.938.
- **Variable-feature selection jitters at the boundary.** 1,998 of 2,000 genes
  shared on PBMC 3k; the two that swap sit at ranks ~1916–2016 where
  standardized variances agree to three decimals.
- **Anything with an RNG differs by its RNG and only by that.**
  `add_module_score` draws control genes at random (96.6 % phase concordance,
  Pearson ≥ 0.998 on the scores); `jack_straw` permutes.
- **shanuz's neighbour search is exact; Seurat's default `annoy` is approximate.**
  When comparing, pass `nn.method = "rann"` on the R side.

Full evidence, with the numbers: <https://genomicai.github.io/shanuz/fidelity/>.
