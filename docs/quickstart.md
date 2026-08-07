# Quickstart

The PBMC 3k guided-clustering workflow, end to end, next to the Seurat code it
mirrors. Every Python line below was run to produce the output shown.

```bash
pip install "truecell[analysis]"
```

## Load and build the object

The dataset downloads on first call into `~/.truecell_data/` (~24 MB). Loaders
return the raw pieces — a counts matrix, gene names, cell names — and
`create_truecell_object` assembles them, applying the same `min_cells` /
`min_features` prefilter `CreateSeuratObject` does.

=== "Python"

    ```python
    import truecell
    from truecell.datasets import pbmc3k

    counts, genes, cells = pbmc3k()
    pbmc = truecell.create_truecell_object(
        counts=counts, feature_names=genes, cell_names=cells,
        project="pbmc3k", min_cells=3, min_features=200,
    )
    truecell.percentage_feature_set(pbmc, pattern=r"^MT-", col_name="percent.mt")
    ```

=== "R"

    ```r
    library(Seurat)

    pbmc <- CreateSeuratObject(Read10X(data.dir), project = "pbmc3k",
                               min.cells = 3, min.features = 200)
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    ```

## Quality control

=== "Python"

    ```python
    md = pbmc.meta_data
    keep = (md["nFeature_RNA"] > 200) & (md["nFeature_RNA"] < 2500) & (md["percent.mt"] < 5)
    pbmc = pbmc.subset(cells=list(md.index[keep]))
    ```

=== "R"

    ```r
    pbmc <- subset(pbmc, subset = nFeature_RNA > 200 &
                                  nFeature_RNA < 2500 & percent.mt < 5)
    ```

Both keep **2,638** of the 2,700 cells — the same 2,638 barcodes, not just the
same count.

## Normalize, select features, scale

=== "Python"

    ```python
    truecell.normalize_data(pbmc, normalization_method="LogNormalize", scale_factor=10000)
    truecell.find_variable_features(pbmc, selection_method="vst", nfeatures=2000)
    truecell.scale_data(pbmc, features=truecell.generics.features(pbmc))
    ```

=== "R"

    ```r
    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize",
                          scale.factor = 10000)
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    pbmc <- ScaleData(pbmc, features = rownames(pbmc))
    ```

!!! note "`scale_data` mutates in place"
    Seurat's functions return a modified object; truecell's write into the one you
    pass and return `None`. `pbmc = truecell.normalize_data(pbmc)` will leave you
    holding `None` — a difference worth internalising early. `subset` is the
    exception: it returns a new object, as it must.

`find_variable_features` picks **1,998 of the same 2,000 genes** Seurat picks.
The two that differ sit at the selection boundary, where the standardized
variances agree to three decimals; see [Fidelity](fidelity.md#what-actually-differs).

## Reduce, cluster, embed

=== "Python"

    ```python
    truecell.run_pca(pbmc, n_pcs=50)
    truecell.find_neighbors(pbmc, dims=range(10), k_param=20)
    truecell.find_clusters(pbmc, resolution=0.5, algorithm=1, random_seed=0)
    truecell.run_umap(pbmc, dims=range(10), seed=42)
    ```

=== "R"

    ```r
    pbmc <- RunPCA(pbmc, npcs = 50)
    pbmc <- FindNeighbors(pbmc, dims = 1:10, k.param = 20, nn.method = "rann")
    pbmc <- FindClusters(pbmc, resolution = 0.5, algorithm = 1, random.seed = 0)
    pbmc <- RunUMAP(pbmc, dims = 1:10, seed.use = 42)
    ```

!!! warning "`dims` is 0-based here, 1-based in R"
    `range(10)` and `1:10` are the same ten PCs. This is the one indexing
    difference in the API, and it follows Python rather than R on purpose.

!!! tip "Use `nn.method = \"rann\"` when comparing against R"
    Seurat's default neighbour search is `annoy`, which is approximate; truecell's
    is exact. Leaving the default in place compares two different neighbour
    tables and reports a difference that belongs to `annoy` rather than to
    either implementation. That trap cost one of the verify scripts a false
    negative of 182 SNN edges.

```
clusters: [692, 515, 458, 344, 301, 159, 155, 14]
```

Eight clusters, against Seurat's nine on the same data, at **ARI 0.899** — the
extra one is a 32-cell dendritic-cell population Seurat's deeper modularity
search separates.

## Markers

=== "Python"

    ```python
    markers = truecell.find_all_markers(
        pbmc, only_pos=True, min_pct=0.25, logfc_threshold=0.25,
    )
    ```

=== "R"

    ```r
    markers <- FindAllMarkers(pbmc, only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25)
    ```

```
cluster   gene  avg_log2FC
      1   FCN1    4.070133
      1 S100A8    6.607992
      3  CD79A    6.911221
      3  MS4A1    5.718520
      7 TMEM40   11.633066
      7 ITGA2B   12.070139
```

`find_markers` runs all eight of Seurat's tests via `test_use=`. On a shared
cell assignment, seven of them reproduce Seurat's top 50 genes exactly and
`avg_log2FC` agrees to 7.1e-15 — [the DE vignette](tutorials/de_vignette.md)
has the full table.

## Plot

Every plotting function returns a matplotlib `Figure`.

```python
fig = truecell.dim_plot(pbmc, reduction="umap", label=True)
fig.savefig("umap.png", dpi=150, bbox_inches="tight")
```

## Next

- [The full PBMC 3k tutorial](tutorials/pbmc3k_tutorial.md) — the same workflow
  with QC plots, elbow plot, feature plots, heatmaps and cell-type annotation,
  and the R comparison at each step.
- [All eighteen tutorials](tutorials/README.md).
- [API reference](api/index.md).
