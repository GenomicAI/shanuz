---
name: truecell-workflow
description: Use when running or debugging a standard single-cell RNA-seq analysis in truecell — loading 10x data, QC filtering, LogNormalize vs SCTransform, variable-feature selection, choosing how many PCs, clustering resolution, UMAP, and cell-type annotation. Covers the decisions the pipeline forces, not just the call order.
---

# The standard truecell workflow

Load the `truecell` skill first for the API contracts (in-place mutation, 0-based
`dims`, features × cells). This skill covers the choices each step demands.

## The whole pipeline

```python
import truecell
from truecell.io import read_10x

counts, genes, cells = read_10x("filtered_gene_bc_matrices/hg19/")
obj = truecell.create_truecell_object(
    counts=counts, feature_names=genes, cell_names=cells,
    project="sample", min_cells=3, min_features=200,
)

truecell.percentage_feature_set(obj, pattern=r"^MT-", col_name="percent.mt")

md = obj.meta_data
keep = (md["nFeature_RNA"] > 200) & (md["nFeature_RNA"] < 2500) & (md["percent.mt"] < 5)
obj = obj.subset(cells=list(md.index[keep]))

truecell.normalize_data(obj)
truecell.find_variable_features(obj, selection_method="vst", nfeatures=2000)
truecell.scale_data(obj)

truecell.run_pca(obj, n_pcs=50)
truecell.find_neighbors(obj, dims=range(10), k_param=20)
truecell.find_clusters(obj, resolution=0.5)
truecell.run_umap(obj, dims=range(10), seed=42)

markers = truecell.find_all_markers(obj, only_pos=True, min_pct=0.25, logfc_threshold=0.25)
```

## Step 1 — Load and build

`create_truecell_object(counts, ...)` wants **features × cells**. `min_cells=3`
drops genes seen in fewer than 3 cells; `min_features=200` drops cells with
fewer than 200 detected genes. Both are applied at construction, before anything
else, exactly as `CreateSeuratObject` does.

Sources: `truecell.io.read_10x(dir)` for a 10x matrix directory,
`truecell.datasets.*` for the cached benchmark datasets,
`truecell.compat.anndata.from_anndata(adata)` for an existing `.h5ad`.

Construction writes `orig.ident`, `nCount_RNA`, `nFeature_RNA` into `meta_data`.

## Step 2 — QC

```python
truecell.percentage_feature_set(obj, pattern=r"^MT-",   col_name="percent.mt")   # human
truecell.percentage_feature_set(obj, pattern=r"^mt-",   col_name="percent.mt")   # mouse
truecell.percentage_feature_set(obj, pattern=r"^RP[SL]", col_name="percent.rb")
```

The pattern is a **Python regex** matched against feature names, so
`^MT-` — not R's `"^MT-"` with `grepl` semantics, though in practice they agree.
Check it caught something before filtering on it:

```python
assert obj.meta_data["percent.mt"].max() > 0, "pattern matched no genes — wrong case/species?"
```

There is no `subset(subset = expr)` string form. Build a mask over `meta_data`
and pass barcodes:

```python
md = obj.meta_data
obj = obj.subset(cells=list(md.index[(md["nFeature_RNA"] > 200) & (md["percent.mt"] < 5)]))
```

**Pick thresholds from the data, not from the tutorial.** 200 / 2500 / 5 % are
the PBMC 3k numbers. Plot first:

```python
fig = truecell.vln_plot(obj, ["nFeature_RNA", "nCount_RNA", "percent.mt"], ncol=3)
fig = truecell.feature_scatter(obj, "nCount_RNA", "percent.mt")
```

## Step 3 — Normalize: LogNormalize or SCTransform

Two arms. Do not mix them in one object's downstream steps.

**LogNormalize** (the default arm):

```python
truecell.normalize_data(obj, normalization_method="LogNormalize", scale_factor=10000)
truecell.find_variable_features(obj, selection_method="vst", nfeatures=2000)
truecell.scale_data(obj)                   # variable features only, as Seurat does
```

**SCTransform** (regularized NB Pearson residuals — one call replaces all three):

```python
truecell.sctransform(obj, vst_flavor="v2")     # writes an "SCT" assay, sets it default
truecell.run_pca(obj, n_pcs=50)                # then continue with more PCs, e.g. dims=range(30)
```

`vst_flavor="v2"` is Seurat 5's model and the default; `"v1"` is the 2019 one.
SCTransform normally supports **more PCs** downstream (30 rather than 10) because
the residuals carry more usable structure.

Regressing covariates out:

```python
truecell.scale_data(obj, vars_to_regress=["percent.mt", "nCount_RNA"])   # LogNormalize arm
truecell.sctransform(obj, vars_to_regress=["percent.mt"])                # SCT arm
```

For protein (ADT) or hashtag assays use CLR, not LogNormalize — see
`truecell-multimodal` for the `margin` argument, which is the thing people get
wrong.

## Step 4 — Variable features

```python
truecell.find_variable_features(obj, selection_method="vst", nfeatures=2000)
hvg = truecell.generics.variable_features(obj)
fig = truecell.variable_feature_plot(obj)
```

`selection_method`:
- `"vst"` — default. Honours `nfeatures`.
- `"mvp"` / `"dispersion"` — honours `mean_cutoff` and `dispersion_cutoff`
  instead, and **ignores `nfeatures`**. That split is Seurat's, not a quirk here.

Expect ~1,998 of 2,000 genes to match R on the same data; the pair that swaps
sits at the selection boundary.

## Step 5 — How many PCs

Three tools, in increasing cost:

```python
fig = truecell.elbow_plot(obj, ndims=50)                    # eyeball the knee
fig = truecell.dim_heatmap(obj, dims=list(range(9)), cells=500)
js = truecell.jack_straw(obj, dims=20, num_replicate=100)   # permutation test
scores = truecell.score_jackstraw(obj)
```

On PBMC 3k both tools keep **13 PCs**. JackStraw's cutoff moves with the seed —
truecell lands on 12/13/14/15 across 60 seeds with mode 13 — so treat |Δ| ≤ 2 as
agreement, not as a discrepancy to chase.

Being generous with PCs costs little; being stingy loses rare populations.
10 PCs is the PBMC 3k tutorial's answer, not a universal one.

## Step 6 — Graph and clusters

```python
truecell.find_neighbors(obj, dims=range(10), k_param=20)    # → graphs["RNA_nn"], ["RNA_snn"]
truecell.find_clusters(obj, resolution=0.5, algorithm=1, random_seed=0)
```

- `dims=range(10)` is R's `1:10`.
- `algorithm`: **1** Louvain (default) · **2** Louvain multilevel · **4** Leiden.
  **3 (SLM) is not implemented** — asking for it is an error, not a silent fallback.
- `resolution` up → more clusters. Sweep it rather than defending one value:

```python
for r in (0.2, 0.4, 0.6, 0.8, 1.0):
    truecell.find_clusters(obj, resolution=r)
    obj.meta_data[f"clusters_res{r}"] = obj.meta_data["seurat_clusters"]
```

Results land in `meta_data["seurat_clusters"]` **and** the active identity —
each `find_clusters` call overwrites both, which is why the sweep stashes a copy.

Expect a cluster count within one of Seurat's. Both tools optimise the same
modularity and land in different local optima; Seurat runs 10 restarts, truecell a
single multilevel pass. On ifnb RPCA the coarser truecell partition scored
**ARI 0.92 against the annotations to Seurat's 0.74** — a different count is not
a worse answer.

## Step 7 — UMAP

```python
truecell.run_umap(obj, dims=range(10), seed=42)          # from a reduction
truecell.run_umap(obj, graph="RNA_snn", seed=42)         # from a precomputed graph
```

`n_neighbors=30` and `min_dist=0.3` are the Seurat defaults. UMAP is for looking
at, never for deciding: cluster on the graph, not on the 2-D coordinates.

## Step 8 — Markers and annotation

```python
markers = truecell.find_all_markers(obj, only_pos=True, min_pct=0.25, logfc_threshold=0.25)
top = markers.groupby("cluster").head(10)

fig = truecell.feature_plot(obj, ["MS4A1", "CD3D", "LYZ", "NKG7"], reduction="umap", ncol=2)
fig = truecell.do_heatmap(obj, list(top["gene"]))
fig = truecell.dot_plot(obj, canonical_markers)

obj = obj.rename_idents({
    "0": "Naive CD4 T", "1": "CD14+ Mono", "2": "Memory CD4 T", "3": "B",
})
obj.meta_data["cell_type"] = obj.idents
```

`rename_idents` returns the object — rebind it. Stash the numeric labels first
with `obj.stash_ident("seurat_clusters_orig")` if you may want them back.

Full DE detail — all eight tests, pseudobulk, conserved markers — in
`truecell-differential-expression`.

## Common failures

| Symptom | Cause |
|---|---|
| `AttributeError: 'NoneType' object has no attribute …` | Rebound an in-place function: `obj = truecell.run_pca(obj)`. |
| QC filter drops every cell | Filtered on a `percent.mt` column of all zeros — wrong pattern case for the species. |
| `KeyError: 'RNA_snn'` in `find_clusters` | `find_neighbors` not run, or run on a different assay; pass `graph_name=`. |
| `run_pca` fails or gives noise | `scale.data` missing (no `scale_data` call), or scaled on too few features. |
| Off-by-one vs R | `dims=range(1, 11)` instead of `range(10)`. |
| `find_clusters(algorithm=3)` | SLM is not implemented; use 1, 2 or 4. |
| Clusters change every run | `random_seed` / `seed` left at different values between runs. |
| Only ~2,000 genes in `scale.data` | Expected — `scale_data` defaults to variable features. Pass `features=truecell.generics.features(obj)` for all genes. |

## Reference runs

- [PBMC 3k guided clustering](https://genomicai.github.io/truecell/tutorials/pbmc3k_tutorial/) — this workflow end to end, with the R comparison at each step.
- [PBMC 8k subclustering](https://genomicai.github.io/truecell/tutorials/advanced_pbmc8k_subclustering/) — plus hierarchical gating and a second clustering pass on the T/NK compartment.
- [SCTransform](https://genomicai.github.io/truecell/tutorials/sctransform_vignette/) — the alternative arm, with the fitted model compared per gene.
