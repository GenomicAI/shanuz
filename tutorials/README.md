# Shanuz Tutorials

Three end-to-end tutorials covering increasingly complex single-cell analysis workflows,
each pairing **R Seurat** code side-by-side with the equivalent **Python Shanuz** code.

---

## Tutorial Overview

| # | Tutorial | Dataset | Key Concepts | Complexity |
|---|----------|---------|--------------|-----------|
| 1 | [PBMC 3k — Guided Clustering](pbmc3k_tutorial.md) | 3,000 PBMCs · 10x Genomics (2016) | QC · Normalization · HVG/VST · PCA · Louvain · UMAP · Markers | Beginner |
| 2 | [PBMC 8k — Advanced Subclustering](advanced_pbmc8k_subclustering.md) | 8,400 PBMCs · GRCh38 · 10x Genomics | All of Tutorial 1 + subclustering, hierarchical cell-type gating, T/NK annotation | Intermediate |
| 3 | [CBMC CITE-seq — Multimodal](multimodal_citeseq.md) | 8,600 CBMCs · RNA + 13 surface proteins | Multi-assay objects · CLR normalization · Protein feature plots · RNA-protein comparison | Advanced |
| 4 | [PBMC 3k — SCTransform](sctransform_vignette.md) | 3,000 PBMCs · 10x Genomics (2016) | Regularized NB normalization · Pearson residuals · `vars.to.regress` · 30-PC workflow · SCT-vs-LogNormalize | Advanced |
| 5 | [Xenium — Spatial (R vs Python)](xenium_spatial_tutorial.md) | 36,602 cells · 10x Xenium mouse brain (CTX+HP) | `load_xenium` · `ImageDimPlot`/`ImageFeaturePlot` · nearest-neighbour distance · local density · `BuildNicheAssay` · `composition_test` — verified to 8 s.f. vs R Seurat | Spatial |

---

## Quick Start

Clone the repo and install dependencies once, then pick any tutorial script:

```bash
git clone https://github.com/GenomicAI/shanuz.git
cd shanuz
uv venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
uv pip install -e ".[analysis]"
```

Each tutorial has a **Python script** that runs the analysis and prints validation output,
and a **figure-generation script** that writes plots to a `figures_*/` subfolder.

---

## Tutorial 1 — PBMC 3k Guided Clustering

> **Walkthrough:** [`pbmc3k_tutorial.md`](pbmc3k_tutorial.md)

A step-by-step port of the official
[Seurat PBMC 3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial).
Covers every step from raw counts to annotated cell types — the entry point for new users.

```bash
python tutorials/pbmc3k_tutorial.py    # runs analysis + validation
python tutorials/generate_plots.py     # writes tutorials/figures/
```

**What you'll learn:**
- Load 10x Genomics data and create a Shanuz object
- Compute QC metrics (percent mitochondrial reads) and filter low-quality cells
- Normalize counts (`LogNormalize`), select highly variable genes (VST), and scale data
- Run PCA, build the KNN/SNN neighbor graph, and cluster with Louvain
- Embed with UMAP and visualize clusters
- Find cluster markers with Wilcoxon rank-sum and annotate 9 cell types

**Key output figures** (in `tutorials/figures/`):

| Figure | Description |
|--------|-------------|
| `01_qc_violin.png` | QC metrics violin plot |
| `03_variable_features.png` | HVG mean–variance plot |
| `07_umap_clusters.png` | UMAP coloured by cluster |
| `08_feature_plots.png` | Canonical marker feature plots |
| `10_marker_heatmap.png` | Top-10 markers per cluster heatmap |
| `11_umap_labeled.png` | UMAP with cell-type labels |

**Validation checkpoints:**

| Step | Expected | Status |
|------|----------|--------|
| Features after filtering | 13,714 | ✅ |
| Cells after QC | 2,638 | ✅ |
| HVG top-10 overlap (≥50%) | PPBP, LYZ, S100A9 … | ✅ 9/10 |
| Clusters at resolution 0.5 | 9 | ✅ |
| All 6 canonical cell types recovered | CD4 T, CD8 T, B, NK, Mono, DC | ✅ |

---

## Tutorial 2 — PBMC 8k Advanced Subclustering

> **Walkthrough:** [`advanced_pbmc8k_subclustering.md`](advanced_pbmc8k_subclustering.md)

A larger dataset (~8,400 cells) that demonstrates the standard Seurat **subclustering**
workflow. After global clustering and broad annotation, the T/NK lymphoid compartment is
isolated and re-analysed to resolve four functionally distinct subsets that global
clustering merges into one cluster.

```bash
python tutorials/pbmc8k_subclustering_tutorial.py   # analysis + validation
python tutorials/generate_advanced_plots.py         # writes tutorials/figures_advanced/
```

**What you'll learn:**
- Run the full pipeline on a larger dataset and use z-scored relative enrichment for
  automated cluster annotation (avoids bias from high-magnitude housekeeping genes)
- Subset a cell lineage and re-run HVG → PCA → neighbors → clusters → UMAP within it
- Apply hierarchical gating to annotate T/NK subsets:
  NK cells (CD3⁻ NKG7/GNLY⁺) → CD8 T (CD8A/B⁺ or cytotoxic CD3⁺, incl. γδ T and MAIT) →
  CD4 Naive (CCR7/SELL⁺) → CD4 Memory

**Key output figures** (in `tutorials/figures_advanced/`):

| Figure | Description |
|--------|-------------|
| `03_umap_global_clusters.png` | UMAP — global clusters |
| `04_umap_global_celltypes.png` | UMAP — broad cell types |
| `07_umap_tnk_subclusters.png` | T/NK subcluster UMAP |
| `08_umap_tnk_subsets.png` | T/NK subsets annotated |
| `09_tnk_subset_featureplots.png` | CD3, CD4, CD8A, CCR7, NKG7 feature plots |
| `11_tnk_markers_heatmap.png` | T/NK subset marker heatmap |

---

## Tutorial 3 — CBMC CITE-seq Multimodal

> **Walkthrough:** [`multimodal_citeseq.md`](multimodal_citeseq.md)

A Python port of Seurat's
[multimodal vignette](https://satijalab.org/seurat/articles/multimodal_vignette)
on the **CBMC CITE-seq** dataset (GSE100866, Stoeckius et al. 2017): ~8,600 cord-blood
mononuclear cells profiled simultaneously for the transcriptome and **13 surface proteins**.

```bash
python tutorials/cbmc_citeseq_tutorial.py     # analysis + validation
python tutorials/generate_multimodal_plots.py # writes tutorials/figures_multimodal/
```

**What you'll learn:**
- Load and align RNA + ADT matrices (handles the mouse spike-in gene prefix)
- Attach a second `ADT` assay to the Shanuz object
- CLR-normalize the surface protein counts (`margin=2` centers per-protein across cells)
- Overlay protein expression on the RNA UMAP — protein gives smooth, high-SNR signal
  where the encoding mRNA is sparse
- Annotate cells with a protein-priority gating strategy + RNA fallback for populations
  the 13-protein panel cannot resolve (platelets, erythroid, pDC, cycling cells)

**Key output figures** (in `tutorials/figures_multimodal/`):

| Figure | Description |
|--------|-------------|
| `01_rna_umap_clusters.png` | UMAP — RNA clusters |
| `02_rna_umap_celltypes.png` | UMAP — protein + RNA cell types |
| `03_adt_featureplots.png` | 8 surface protein feature plots on RNA UMAP |
| `04_protein_vs_rna.png` | Protein vs RNA side-by-side (CD19, CD3, CD8, CD14) |
| `05_adt_ridgeplots.png` | ADT ridge plots by cell type |
| `06_adt_scatter_CD4_CD8.png` | CD4 vs CD8 protein scatter |
| `07_adt_scatter_CD19_CD3.png` | CD19 vs CD3 protein scatter |

---

## Tutorial 4 — PBMC 3k SCTransform

> **Walkthrough:** [`sctransform_vignette.md`](sctransform_vignette.md)

A Python port of Seurat's
[sctransform vignette](https://satijalab.org/seurat/articles/sctransform_vignette).
`SCTransform` replaces the `NormalizeData → FindVariableFeatures → ScaleData`
trio with a single regularized negative-binomial model whose **Pearson
residuals** are the normalized values, then clusters over 30 PCs. The tutorial
runs SCTransform *and* the standard log-normalization workflow on the same cells
to compare their resolution.

```bash
python tutorials/pbmc3k_sctransform_tutorial.py   # analysis + SCT-vs-std comparison
python tutorials/generate_sctransform_plots.py    # writes tutorials/figures_sctransform/
```

**What you'll learn:**
- Run `sctransform(vars_to_regress=["percent.mt"])` — one call replacing three,
  producing a new `SCT` assay with corrected counts and residual `scale.data`
- Cluster over 30 PCs (the vignette's deeper embedding) and read the resulting
  cytotoxicity gradient: Naive CD4 → Memory CD4 → CD8 Effector → NK
- Compare SCTransform against log-normalization on identical cells

**Key output figures** (in `tutorials/figures_sctransform/`):

| Figure | Description |
|--------|-------------|
| `01_sct_umap_clusters.png` | UMAP — SCT clusters |
| `02_sct_umap_celltypes.png` | UMAP — annotated cell types |
| `03_sct_featureplots_1.png` | CD8A, GZMK, CCL5, S100A4, ANXA1, CCR7 |
| `04_sct_featureplots_2.png` | CD3D, ISG15, TCL1A, FCER2, XCL1, FCGR3A |
| `05_sct_violins.png` | Vignette marker violins per cluster |
| `06_sct_vs_std_umap.png` | SCTransform vs LogNormalize, side by side |

**Accuracy vs R:** exact match on the model, the **3,000 variable features**, the
30-PC embedding, and the recovered biology (CD8-effector/CD4/NK split, marker
patterns). Exact cluster count differs — `sctransform` is a faithful but not
bit-identical reimplementation (LOESS/`theta.ml` differ) and clustering/UMAP use
different libraries, the same caveat as the other tutorials.

---

## Tutorial 5 — Xenium Spatial (R vs Python)

> **Walkthrough:** [`xenium_spatial_tutorial.md`](xenium_spatial_tutorial.md)

The spatial counterpart to the others: a 10x **Xenium** section (mouse brain
coronal CTX+HP subset — the dataset in Seurat's spatial vignette, **36,602 cells
× 248 genes**) analysed side by side in R Seurat and Shanuz. It reproduces the
style of an internal Xenium mast-cell/neighbourhood workflow on fully public
data, exercising the spatial Seurat-parity layer end to end.

```bash
python tutorials/generate_spatial_plots.py     # auto-downloads ~20 MB → figures_spatial/
Rscript tutorials/xenium_spatial_verify.R      # R reference figures + r_reference.json
python tutorials/compare_xenium_anchors.py     # prints the R-vs-Python parity table
```

**What you'll learn:**
- `load_xenium` — build an object with expression **and** centroids, keeping only
  `Gene Expression` features (like `LoadXenium`)
- Deterministic marker-panel cell typing (the `KIT+ TPSAB1+`-style rule)
- `image_dim_plot` / `image_feature_plot` — plot cells in tissue space (immune to
  the `ggplot2` 4.x `ImageDimPlot` blank-render bug)
- `nearest_neighbor_distance` / `local_neighborhood` — `FNN::get.knn` idioms
- `build_niche_assay` — neighbourhood-composition niches (`BuildNicheAssay`)
- `composition_test` — Fisher + BH enrichment across a spatial split

**Key output figures** (in `tutorials/figures_spatial/`, `r_*` = R Seurat):

| Figure | Description |
|--------|-------------|
| `03_image_celltype.png` | Marker cell types in tissue space |
| `05_image_feature_Slc17a7.png` | Excitatory-neuron marker in space |
| `06_image_niches.png` | Neighbourhood niches |
| `07_image_focal.png` | Focal (Vascular) cells highlighted |

**Accuracy vs R:** every deterministic anchor matches **to 8 significant
figures** — cell counts, all cell-type counts, nearest-neighbour distances, local
density, and the composition test (log2 ratios, Fisher p, BH padj, χ² p).
Clustering / UMAP / niche layout are stochastic and agree in structure only, the
same caveat as the other tutorials.

---

## R Seurat → Shanuz API Quick Reference

| Task | R (Seurat) | Python (Shanuz) |
|------|-----------|-----------------|
| Create object | `CreateSeuratObject(counts, min.cells, min.features)` | `create_shanuz_object(counts, min_cells, min_features)` |
| % mito genes | `PercentageFeatureSet(pbmc, pattern="^MT-")` | `percentage_feature_set(pbmc, pattern=r"^MT-")` |
| Normalize | `NormalizeData(pbmc, method, scale.factor)` | `normalize_data(pbmc, normalization_method, scale_factor)` |
| CLR (ADT) | `NormalizeData(pbmc, method="CLR", margin=2)` | `normalize_data(pbmc, normalization_method="CLR", margin=2)` |
| SCTransform | `SCTransform(pbmc, vars.to.regress="percent.mt")` | `sctransform(pbmc, vars_to_regress=["percent.mt"])` |
| HVGs | `FindVariableFeatures(pbmc, selection.method, nfeatures)` | `find_variable_features(pbmc, selection_method, nfeatures)` |
| Scale | `ScaleData(pbmc, features)` | `scale_data(pbmc, features)` |
| PCA | `RunPCA(pbmc, features, npcs)` | `run_pca(pbmc, features, n_pcs)` |
| Neighbors | `FindNeighbors(pbmc, dims)` | `find_neighbors(pbmc, dims, k_param)` |
| Cluster | `FindClusters(pbmc, resolution)` | `find_clusters(pbmc, resolution, algorithm)` |
| UMAP | `RunUMAP(pbmc, dims)` | `run_umap(pbmc, dims)` |
| Markers | `FindMarkers(pbmc, ident.1)` | `find_markers(pbmc, ident_1)` |
| All markers | `FindAllMarkers(pbmc, only.pos, logfc.threshold)` | `find_all_markers(pbmc, only_pos, logfc_threshold)` |
| Conserved markers | `FindConservedMarkers(pbmc, ident.1, grouping.var)` | `find_conserved_markers(pbmc, ident_1, grouping_var)` |
| Pseudobulk | `AggregateExpression(pbmc, group.by)` | `aggregate_expression(pbmc, group_by)` |
| Pseudobulk DESeq2 | `FindMarkers(pbmc, test.use="DESeq2")` | `find_markers(pbmc, ident_1, test_use="deseq2", sample_col=...)` |
| MAST hurdle DE | `FindMarkers(pbmc, test.use="MAST")` | `find_markers(pbmc, ident_1, test_use="mast")` |
| Bimodal LRT DE | `FindMarkers(pbmc, test.use="bimod")` | `find_markers(pbmc, ident_1, test_use="bimod")` |
| Rename idents | `RenameIdents(pbmc, new.ids)` | `pbmc.rename_idents(mapping_dict)` |
| Subset cells | `subset(pbmc, subset = condition)` | `pbmc.subset(cells=keep_list)` |
| Add assay | `pbmc[["ADT"]] <- CreateAssayObject(counts)` | `obj.assays["ADT"] = create_assay5_object(counts, key="adt_")` |
| Switch assay | `DefaultAssay(cbmc) <- "ADT"` | `feature_plot(..., assay="ADT")` |
| Access metadata | `pbmc@meta.data` | `pbmc.meta_data` |
| Access assay | `pbmc[["RNA"]]` | `pbmc.assays["RNA"]` |
| Active idents | `Idents(pbmc)` | `pbmc.idents` |

### Pseudobulk & conserved markers

Two multi-sample DE helpers (mirroring Seurat's `AggregateExpression` and
`FindConservedMarkers`):

```python
from shanuz import aggregate_expression, find_conserved_markers

# Pseudobulk: sum raw counts per (cell type × donor) → features × groups DataFrame.
# Pass return_object=True to get a Shanuz object with one "cell" per group instead
# (the standard input for pyDESeq2-style sample-level testing).
pb = aggregate_expression(obj, group_by=["cell_type", "donor"])

# Conserved markers: genes up in cluster "B" in *every* condition. Runs FindMarkers
# per level of grouping_var, keeps genes significant in all, and combines their
# p-values with Fisher's method (the `combined_p_val` column; `max_pval` is the
# worst single-condition p-value).
cons = find_conserved_markers(obj, ident_1="B", grouping_var="condition",
                              only_pos=True)
cons.head()   # per-condition stats + max_pval + combined_p_val, sorted by combined_p_val
```

Pseudobulk DESeq2 (`find_markers(test_use="deseq2")`) tests **between conditions**
rather than between clusters: set `obj.idents` to the two conditions, aggregate to
one profile per replicate (`sample_col`), and run DESeq2 on those samples. Needs
`pip install shanuz[deseq2]`:

```python
obj.idents = obj.meta_data["condition"]              # e.g. "stim" vs "ctrl"
de = find_markers(obj, ident_1="stim", ident_2="ctrl",
                  test_use="deseq2", sample_col="donor")
de.head()   # p_val / avg_log2FC (DESeq2 log2FoldChange, +ve = up in stim) / pct.1 / pct.2 / p_val_adj
```

### Plotting

| R (Seurat) | Python (Shanuz) |
|-----------|-----------------|
| `VlnPlot(pbmc, features, slot)` | `vln_plot(pbmc, features, layer)` |
| `FeaturePlot(pbmc, features)` | `feature_plot(pbmc, features, assay)` |
| `DimPlot(pbmc, reduction, label, pt.size)` | `dim_plot(pbmc, reduction, label, pt_size)` |
| `ElbowPlot(pbmc)` | `elbow_plot(pbmc)` |
| `FeatureScatter(pbmc, feature1, feature2)` | `feature_scatter(pbmc, feature1, feature2)` |
| `VariableFeaturePlot(pbmc)` | `variable_feature_plot(pbmc)` |
| `VizDimLoadings(pbmc, dims, reduction)` | `viz_dim_loadings(pbmc, dims, reduction)` |
| `DimHeatmap(pbmc, dims, cells)` | `dim_heatmap(pbmc, dims, cells)` |
| `DoHeatmap(pbmc, features)` | `do_heatmap(pbmc, features)` |
| `RidgePlot(pbmc, features, ncol)` | `ridge_plot(pbmc, features, ncol)` |
| `ImageDimPlot(obj, group.by)` | `image_dim_plot(obj, group_by)` |
| `ImageFeaturePlot(obj, features)` | `image_feature_plot(obj, feature)` |

### Spatial

| R (Seurat) | Python (Shanuz) |
|-----------|-----------------|
| `LoadXenium(dir)` / `Load10X_Spatial` / `LoadNanostring` | `load_xenium(dir)` / `load_visium(dir)` / `load_cosmx(dir)` |
| `LoadVizgen(dir)` (MERSCOPE) | `load_merscope(dir)` — drops `Blank-*` controls by default |
| `GetTissueCoordinates(obj)` | `get_tissue_coordinates(obj)` |
| `FNN::get.knn(coords, k)` / `get.knnx` | `spatial_knn(coords, k, query)` |
| `FNN::get.knn` (nearest same-type) | `nearest_neighbor_distance(obj, group_by, reference)` |
| *(hand-rolled neighbourhood counts)* | `local_neighborhood(obj, group_by, reference, k)` |
| `BuildNicheAssay(obj, fov, group.by, niches.k)` | `build_niche_assay(obj, group_by, k, niches)` |
| *(hand-rolled Fisher + `p.adjust`)* | `composition_test(obj, group_by, split_by)` |

> **Plot output:** R renders to the graphics device automatically. Shanuz functions return a
> `matplotlib.Figure` — call `fig.savefig("out.png")` to save or display inline in Jupyter.
