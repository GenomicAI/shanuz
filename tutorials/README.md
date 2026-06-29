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

## R Seurat → Shanuz API Quick Reference

| Task | R (Seurat) | Python (Shanuz) |
|------|-----------|-----------------|
| Create object | `CreateSeuratObject(counts, min.cells, min.features)` | `create_shanuz_object(counts, min_cells, min_features)` |
| % mito genes | `PercentageFeatureSet(pbmc, pattern="^MT-")` | `percentage_feature_set(pbmc, pattern=r"^MT-")` |
| Normalize | `NormalizeData(pbmc, method, scale.factor)` | `normalize_data(pbmc, normalization_method, scale_factor)` |
| CLR (ADT) | `NormalizeData(pbmc, method="CLR", margin=2)` | `normalize_data(pbmc, normalization_method="CLR", margin=2)` |
| HVGs | `FindVariableFeatures(pbmc, selection.method, nfeatures)` | `find_variable_features(pbmc, selection_method, nfeatures)` |
| Scale | `ScaleData(pbmc, features)` | `scale_data(pbmc, features)` |
| PCA | `RunPCA(pbmc, features, npcs)` | `run_pca(pbmc, features, n_pcs)` |
| Neighbors | `FindNeighbors(pbmc, dims)` | `find_neighbors(pbmc, dims, k_param)` |
| Cluster | `FindClusters(pbmc, resolution)` | `find_clusters(pbmc, resolution, algorithm)` |
| UMAP | `RunUMAP(pbmc, dims)` | `run_umap(pbmc, dims)` |
| Markers | `FindMarkers(pbmc, ident.1)` | `find_markers(pbmc, ident_1)` |
| All markers | `FindAllMarkers(pbmc, only.pos, logfc.threshold)` | `find_all_markers(pbmc, only_pos, logfc_threshold)` |
| Rename idents | `RenameIdents(pbmc, new.ids)` | `pbmc.rename_idents(mapping_dict)` |
| Subset cells | `subset(pbmc, subset = condition)` | `pbmc.subset(cells=keep_list)` |
| Add assay | `pbmc[["ADT"]] <- CreateAssayObject(counts)` | `obj.assays["ADT"] = create_assay5_object(counts, key="adt_")` |
| Switch assay | `DefaultAssay(cbmc) <- "ADT"` | `feature_plot(..., assay="ADT")` |
| Access metadata | `pbmc@meta.data` | `pbmc.meta_data` |
| Access assay | `pbmc[["RNA"]]` | `pbmc.assays["RNA"]` |
| Active idents | `Idents(pbmc)` | `pbmc.idents` |

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

> **Plot output:** R renders to the graphics device automatically. Shanuz functions return a
> `matplotlib.Figure` — call `fig.savefig("out.png")` to save or display inline in Jupyter.
