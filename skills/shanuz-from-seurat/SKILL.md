---
name: shanuz-from-seurat
description: Use when porting R Seurat code to Python shanuz, translating a Seurat function or argument name, or comparing the two tools' numeric output. Carries the full name map, the translation rules for R idioms with no direct Python equivalent, and the traps that make a port look wrong when it is right.
---

# Porting Seurat to shanuz

Load the `shanuz` skill first for the API contracts.

## The five translation rules

**1. `FunctionName` → `function_name`, `arg.name` → `arg_name`.**
`FindVariableFeatures(selection.method=)` → `find_variable_features(selection_method=)`.
Names colliding with Python keywords take a trailing underscore: `lambda` →
`lambda_`, `type` → `type_`.

**2. Seurat returns a modified object; shanuz mutates in place and returns `None`.**

```r
pbmc <- NormalizeData(pbmc)          # R
```
```python
shanuz.normalize_data(pbmc)          # Python — do NOT reassign
```

Rebind only for `subset`, `merge`, `sketch_data`, `integrate_data`, and the
`rename_idents` / `stash_ident` family.

**3. `dims = 1:10` → `dims=range(10)`.** 0-based, and this is the only indexing
difference in the API.

**4. `rownames(obj)` → `shanuz.generics.features(obj)`; `colnames(obj)` →
`obj.cell_names()`.** Matrices are features × cells in both tools.

**5. `subset(obj, subset = expr)` has no string-expression form.** Build a mask
over `meta_data` and pass barcodes:

```r
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 5)
```
```python
md = pbmc.meta_data
pbmc = pbmc.subset(cells=list(md.index[(md["nFeature_RNA"] > 200) & (md["percent.mt"] < 5)]))
```

## Slot and accessor map

| R | Python |
|---|---|
| `pbmc@meta.data` | `pbmc.meta_data` |
| `pbmc[["RNA"]]` | `pbmc.assays["RNA"]` |
| `pbmc[["pca"]]` | `pbmc.reductions["pca"]` |
| `pbmc@graphs$RNA_snn` | `pbmc.graphs["RNA_snn"]` |
| `pbmc@images$slice1` | `pbmc.images["slice1"]` |
| `pbmc@commands` | `pbmc.commands` |
| `pbmc@misc` | `pbmc.misc` |
| `Idents(pbmc)` | `pbmc.idents` |
| `Idents(pbmc) <- "col"` | `pbmc.set_ident(cells, ident)` / set from a metadata column |
| `Cells(pbmc)` | `pbmc.cell_names()` |
| `Features(pbmc)` / `rownames(pbmc)` | `shanuz.generics.features(pbmc)` |
| `LayerData(obj, "data")` | `obj.get_assay().layer_data("data")` |
| `Layers(obj)` | `obj.get_assay().layers_list()` |
| `JoinLayers(obj)` | `obj.get_assay().join_layers()` |
| `split(obj, f = obj$batch)` | `obj.get_assay().split_layers(f)` |
| `Embeddings(obj, "pca")` | `obj.embeddings("pca")` |
| `Loadings(obj[["pca"]])` | `shanuz.generics.loadings(obj.reductions["pca"])` |
| `Stdev(obj, "pca")` | `shanuz.generics.stdev(obj.reductions["pca"])` |
| `VariableFeatures(obj)` | `shanuz.generics.variable_features(obj)` |
| `FetchData(obj, vars)` | `obj.fetch_data(vars)` |
| `AddMetaData(obj, df)` | `obj.add_meta_data(df)` |
| `WhichCells(obj, idents=)` | `obj.which_cells(ident=)` |
| `RenameIdents(obj, ...)` | `obj.rename_idents({...})` — rebind |
| `merge(a, b)` | `a.merge(b)` |
| `DefaultAssay(obj) <- "ADT"` | pass `assay="ADT"` to the call |
| `Read10X(dir)` | `shanuz.io.read_10x(dir)` |
| `Radius(img)` | `shanuz.generics.radius(img)` |
| `GetTissueCoordinates(obj)` | `shanuz.get_tissue_coordinates(obj)` |

## Function map

| Task | R (Seurat) | Python (shanuz) |
|---|---|---|
| Create object | `CreateSeuratObject(counts, min.cells, min.features)` | `create_shanuz_object(counts, min_cells=, min_features=)` |
| % mito | `PercentageFeatureSet(obj, pattern="^MT-")` | `percentage_feature_set(obj, pattern=r"^MT-")` |
| Normalize | `NormalizeData(obj, normalization.method, scale.factor)` | `normalize_data(obj, normalization_method=, scale_factor=)` |
| CLR (ADT) | `NormalizeData(obj, method="CLR", margin=2)` | `normalize_data(obj, normalization_method="CLR", margin=2)` |
| SCTransform | `SCTransform(obj, vars.to.regress="percent.mt")` | `sctransform(obj, vars_to_regress=["percent.mt"])` |
| HVGs | `FindVariableFeatures(obj, selection.method, nfeatures)` | `find_variable_features(obj, selection_method=, nfeatures=)` |
| Scale | `ScaleData(obj, features)` | `scale_data(obj, features=)` |
| PCA | `RunPCA(obj, npcs)` | `run_pca(obj, n_pcs=)` |
| ICA / t-SNE | `RunICA` / `RunTSNE` | `run_ica` / `run_tsne` |
| Supervised PCA | `RunSPCA(obj, graph="wsnn")` | `run_spca(obj, graph="wsnn")` |
| GLM-PCA | `RunGLMPCA(obj, L=10)` | `glm_pca(obj, n_components=10, family="poisson"\|"nb")` |
| JackStraw | `JackStraw` / `ScoreJackStraw` | `jack_straw` / `score_jackstraw` |
| Neighbors | `FindNeighbors(obj, dims=1:10)` | `find_neighbors(obj, dims=range(10))` |
| Cluster | `FindClusters(obj, resolution)` | `find_clusters(obj, resolution=, algorithm=)` |
| UMAP | `RunUMAP(obj, dims=1:10)` | `run_umap(obj, dims=range(10))` |
| WNN | `FindMultiModalNeighbors(reduction.list, dims.list)` | `find_multi_modal_neighbors(obj, reduction_list=, dims_list=)` |
| Harmony | `RunHarmony(obj, "batch")` | `run_harmony(obj, "batch")` or `integrate_layers(obj, method="harmony", group_by="batch")` |
| v5 integration | `IntegrateLayers(obj, method=CCAIntegration)` | `integrate_layers(obj, method="cca", group_by="batch")` |
| v4 anchors | `FindIntegrationAnchors(list)` → `IntegrateData(anchors)` | `find_integration_anchors(objs)` → `integrate_data(anchors)` |
| Transfer anchors | `FindTransferAnchors(reference, query)` | `find_transfer_anchors(reference, query)` |
| Transfer labels | `TransferData(anchors, refdata=ref$celltype)` | `transfer_data(anchors, refdata="celltype")` |
| Map query | `MapQuery(anchors, query, reference, refdata=list(...))` | `map_query(anchors, refdata="celltype")` |
| Project UMAP | `ProjectUMAP(query, reference, reduction.model="umap")` | `project_umap(query, reference)` |
| Leverage / sketch | `LeverageScore` / `SketchData` / `ProjectData` | `leverage_score` / `sketch_data` / `project_data` |
| On-disk matrix | `write_matrix_dir` / `open_matrix_dir` (BPCells) | `write_lazy_matrix` / `open_lazy_matrix` |
| Hashing | `HTODemux` / `MULTIseqDemux` | `hto_demux` / `multiseq_demux` |
| Mixscape | `CalcPerturbSig` / `RunMixscape` / `MixscapeLDA` | `calc_perturb_sig` / `run_mixscape` / `mixscape_lda` |
| Markers | `FindMarkers(obj, ident.1)` | `find_markers(obj, ident_1=)` |
| All markers | `FindAllMarkers(obj, only.pos, logfc.threshold)` | `find_all_markers(obj, only_pos=, logfc_threshold=)` |
| Conserved markers | `FindConservedMarkers(obj, ident.1, grouping.var)` | `find_conserved_markers(obj, ident_1=, grouping_var=)` |
| Pseudobulk | `AggregateExpression(obj, group.by)` | `aggregate_expression(obj, group_by=)` |
| DESeq2 / MAST / bimod | `FindMarkers(test.use="DESeq2"\|"MAST"\|"bimod")` | `find_markers(test_use="deseq2"\|"mast"\|"bimod")` — `deseq2` also needs `sample_col=` |
| Module / cell cycle | `AddModuleScore` / `CellCycleScoring` | `add_module_score` / `cell_cycle_scoring` |
| Spatial loaders | `LoadXenium` / `Load10X_Spatial` / `LoadNanostring` / `LoadVizgen` | `load_xenium` / `load_visium` / `load_cosmx` / `load_merscope` |
| Niches / SVF | `BuildNicheAssay` / `FindSpatiallyVariableFeatures` | `build_niche_assay` / `find_spatially_variable_features` |
| Plots | `DimPlot` / `FeaturePlot` / `VlnPlot` / `DoHeatmap` … | `dim_plot` / `feature_plot` / `vln_plot` / `do_heatmap` … |

## Idioms with no direct translation

**`DefaultAssay(obj) <- "ADT"`** — there is no stateful default-assay switch for
analysis calls. Pass `assay="ADT"` to whichever function needs it.

**`obj$col <- value`** — `obj.meta_data["col"] = value`, or
`obj.add_meta_data(series, col_name="col")`.

**Adding an assay:**
```r
pbmc[["ADT"]] <- CreateAssayObject(counts = adt)
```
```python
from shanuz.assay5 import create_assay5_object
obj.assays["ADT"] = create_assay5_object(counts=adt, feature_names=proteins,
                                         cell_names=obj.cell_names(), key="adt_")
```

**`SplitObject(obj, split.by = "batch")`** — v5 splits *layers*, not objects:
`obj.get_assay().split_layers(f)`, and `join_layers()` to put them back.

**`FindClusters(algorithm = 3)`** — SLM is not implemented. Use 1 (Louvain),
2 (Louvain multilevel), or 4 (Leiden).

## Comparing the two tools honestly

If you are checking a port against R, these differences are real and expected.
Chasing them wastes days.

**Pass `nn.method = "rann"` on the R side.** Seurat's default neighbour search is
`annoy`, which is approximate; shanuz's is exact. Comparing the defaults compares
two different neighbour tables and reports a difference belonging to `annoy`.
Skipping this once cost a verify script a false negative of 182 SNN edges.

**Hand both sides the same cell assignment** when testing anything downstream of
clustering. Otherwise a clustering difference surfaces as a DE difference and
gets blamed on the wrong function.

**Expect these, and do not report them as regressions:**

| Difference | Size | Why |
|---|---|---|
| Louvain cluster count | ±1 | Same algorithm, different local optimum. PBMC 3k: 8 vs 9 at ARI 0.938. Seurat runs 10 restarts; shanuz a single multilevel pass. |
| Variable features | ~2 of 2,000 | Boundary jitter — the swapped genes agree on standardized variance to three decimals. |
| `add_module_score`, `jack_straw` | RNG-sized | Random control genes / permutations. 96.6 % phase concordance, Pearson ≥ 0.998. Prove differences distribution-against-distribution over matched seeds, never from a single pair. |
| JackStraw PC cutoff | \|Δ\| ≤ 2 | R seeds each replicate from its loop index and is deterministic at 13; shanuz seeds from `seed` and has mode 13 over 60 seeds. |
| `deseq2` top-50 overlap | 15–32 genes | A divergence measurement, not a parity target. Reaching 50 would mean `sample_col` had stopped being honoured. |
| Visium spot radius | shanuz = half of Seurat's | `spot_diameter_fullres` is a diameter; Seurat stores it in a `radius` slot. Here Seurat is the one that is wrong. |
| R's `clara` | arm64 vs x86_64 | R's own function is architecture-dependent. shanuz targets IEEE/x86_64 semantics on purpose. |

**And these should match, so investigate if they don't:** `avg_log2FC` (7.1e-15),
CLR (4.2e-15), Moran's I (1.6e-14), the object-model accessors (91 of 91 anchors
exact), seven of the eight DE tests reproducing Seurat's top 50 exactly.

## A worked port

```r
library(Seurat)
pbmc <- CreateSeuratObject(Read10X(dir), project = "pbmc3k",
                           min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 &
                              nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, npcs = 50)
pbmc <- FindNeighbors(pbmc, dims = 1:10, nn.method = "rann")
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25)
```

```python
import shanuz
from shanuz.io import read_10x

counts, genes, cells = read_10x(dir)
pbmc = shanuz.create_shanuz_object(counts=counts, feature_names=genes, cell_names=cells,
                                   project="pbmc3k", min_cells=3, min_features=200)
shanuz.percentage_feature_set(pbmc, pattern=r"^MT-", col_name="percent.mt")

md = pbmc.meta_data
pbmc = pbmc.subset(cells=list(md.index[
    (md["nFeature_RNA"] > 200) & (md["nFeature_RNA"] < 2500) & (md["percent.mt"] < 5)]))

shanuz.normalize_data(pbmc)
shanuz.find_variable_features(pbmc, selection_method="vst", nfeatures=2000)
shanuz.scale_data(pbmc, features=shanuz.generics.features(pbmc))
shanuz.run_pca(pbmc, n_pcs=50)
shanuz.find_neighbors(pbmc, dims=range(10))
shanuz.find_clusters(pbmc, resolution=0.5)
shanuz.run_umap(pbmc, dims=range(10))
markers = shanuz.find_all_markers(pbmc, only_pos=True, min_pct=0.25)
```

Both keep the same 2,638 barcodes — the same barcodes, not just the same count.

Full evidence: <https://genomicai.github.io/shanuz/fidelity/>.
