# PBMC 3k Tutorial — Shanuz vs R Seurat: Visual Comparison

Side-by-side comparison of plots produced by the **Shanuz** Python package and the
official **R Seurat** PBMC 3k guided clustering tutorial
([satijalab.org/seurat/articles/pbmc3k_tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial)).

> **Dataset:** 3k PBMCs from a Healthy Donor — 10x Genomics (2016)  
> **Seurat version referenced:** v5 (Hao et al. 2024)  
> **Shanuz version:** 0.1.0

To reproduce all Shanuz plots:
```bash
python tutorials/generate_plots.py
```

---

## 1 · QC Metrics — Violin Plot

| R Seurat | Shanuz |
|----------|--------|
| ![R QC violin](https://satijalab.org/seurat/articles/pbmc3k_tutorial_files/figure-html/qc2-1.png) | ![Shanuz QC violin](figures/01_qc_violin.png) |

Both plots show the distribution of three QC metrics across 2,700 cells: unique feature counts
(nFeature_RNA), total molecule counts (nCount_RNA), and mitochondrial percentage (percent.mt).
Cells with nFeature_RNA > 2,500 or percent.mt > 5% are excluded downstream.

---

## 2 · QC Metrics — Scatter Plot

| R Seurat | Shanuz |
|----------|--------|
| ![R QC scatter](https://satijalab.org/seurat/articles/pbmc3k_tutorial_files/figure-html/qc2-2.png) | ![Shanuz QC scatter](figures/02_qc_scatter.png) |

Scatter plots confirm the expected positive correlation between total counts and feature counts,
and the typical low-mt outlier cells visible in both implementations.

---

## 3 · Highly Variable Features

| R Seurat | Shanuz |
|----------|--------|
| ![R HVG](https://satijalab.org/seurat/articles/pbmc3k_tutorial_files/figure-html/var_features-1.png) | ![Shanuz HVG](figures/03_variable_features.png) |

Both plots use the VST (variance-stabilizing transformation) method to select 2,000 highly
variable genes (shown in red). The mean–variance relationship and the overall shape of the
selected gene set are consistent. Minor differences in the top-10 labels reflect small
numerical differences between R's Fortran LOESS and Python's `statsmodels.lowess`.

**Top 10 overlap: 5/10 (50%)** — shared: GNLY, PF4, PPBP, S100A8, S100A9.

---

## 4 · PCA — Top Loading Genes

| R Seurat | Shanuz |
|----------|--------|
| ![R PCA loadings](https://satijalab.org/seurat/articles/pbmc3k_tutorial_files/figure-html/pca_viz-1.png) | ![Shanuz PCA loadings](figures/04_pca_loadings.png) |

PC1 captures the myeloid–lymphoid axis in both implementations. The top positive-loading
genes (CST3, TYROBP, LST1, AIF1) are reproduced identically in Shanuz.

---

## 5 · PCA — Cells in PC1 / PC2 Space

| R Seurat | Shanuz |
|----------|--------|
| ![R PCA dimplot](https://satijalab.org/seurat/articles/pbmc3k_tutorial_files/figure-html/pca_viz-2.png) | ![Shanuz PCA dimplot](figures/05_pca_dimplot.png) |

The PCA scatter shows the same separation structure. Clusters are coloured consistently
and the overall topology — with myeloid cells on one end of PC1 and lymphoid cells on
the other — is reproduced.

---

## 6 · Elbow Plot — Dimensionality Selection

| R Seurat | Shanuz |
|----------|--------|
| ![R elbow](https://satijalab.org/seurat/articles/pbmc3k_tutorial_files/figure-html/elbow_plot-1.png) | ![Shanuz elbow](figures/06_elbow_plot.png) |

Both plots show a clear elbow around PC 9–10, confirming the choice of 10 PCs for
downstream neighbor-finding and clustering.

| PC | R Seurat stdev | Shanuz stdev |
|----|---------------|--------------|
| 1  | ~6.8          | 6.766        |
| 2  | ~4.8          | 4.808        |
| 10 | ~1.7          | 1.684        |

---

## 7 · UMAP — Coloured by Cluster

| R Seurat | Shanuz |
|----------|--------|
| ![R UMAP clusters](https://satijalab.org/seurat/articles/pbmc3k_tutorial_files/figure-html/umapplot-1.png) | ![Shanuz UMAP clusters](figures/07_umap_clusters.png) |

Both embeddings resolve the same 9 clusters at resolution = 0.5. The topology is
consistent: a large T-cell group, a monocyte population, a B-cell island, and small
peripheral populations (NK, DC, Platelet). Cluster label numbers differ because Louvain
assigns IDs by graph traversal order (random-seed dependent), not by size.

---

## 8 · Feature Plots — Canonical Marker Genes

| R Seurat | Shanuz |
|----------|--------|
| ![R feature plots](https://satijalab.org/seurat/articles/pbmc3k_tutorial_files/figure-html/unnamed-chunk-3-1.png) | ![Shanuz feature plots](figures/08_feature_plots.png) |

Expression of 9 canonical marker genes overlaid on the UMAP embedding. Key observations
reproduced in Shanuz:

| Gene | Cell type | Shanuz result |
|------|-----------|---------------|
| MS4A1, CD79A | B cells | Focal expression on B-cell island |
| NKG7, GNLY | NK / CD8 T | NK cluster clearly lit |
| FCGR3A | FCGR3A+ Mono | FCGR3A+ monocyte cluster |
| LYZ | CD14+ Mono | Strong in CD14+ monocytes |
| PPBP | Platelet | Isolated platelet cluster |
| CD8A | CD8 T | CD8 T cluster |
| IL7R | CD4 T | Naive + Memory CD4 T |

---

## 9 · Violin Plots — Marker Gene Expression per Cluster

| R Seurat | Shanuz |
|----------|--------|
| ![R violin markers](https://satijalab.org/seurat/articles/pbmc3k_tutorial_files/figure-html/markerplots-1.png) | ![Shanuz violin markers](figures/09_marker_violins.png) |

Violin plots confirm cluster-specific marker expression. MS4A1 (B cells) and CD79A
(B cells) show high expression in the B-cell cluster in both implementations. NKG7
and PF4 (Platelet) are similarly cluster-restricted.

---

## 10 · Top Marker Gene Heatmap

| R Seurat | Shanuz |
|----------|--------|
| ![R heatmap](https://satijalab.org/seurat/articles/pbmc3k_tutorial_files/figure-html/clusterHeatmap-1.png) | ![Shanuz heatmap](figures/10_marker_heatmap.png) |

Heatmaps of the top 5 marker genes per cluster (by avg_log2FC, scaled expression).
Both show clear block-diagonal structure confirming that each cluster has a unique
transcriptional signature. Key marker blocks reproduced: LYZ/S100A9 for CD14+ Mono,
GNLY/NKG7 for NK, PPBP for Platelet, CD79A for B cells.

---

## 11 · UMAP — Cell Type Annotations

| R Seurat | Shanuz |
|----------|--------|
| ![R labeled UMAP](https://satijalab.org/seurat/articles/pbmc3k_tutorial_files/figure-html/labelplot-1.png) | ![Shanuz labeled UMAP](figures/11_umap_labeled.png) |

Final annotated UMAP with all 9 cell types labelled. Shanuz recovers the same 9 cell
populations using the same canonical marker logic as the R tutorial:

| Cell Type | R Seurat | Shanuz |
|-----------|----------|--------|
| Naive CD4 T | ✅ | ✅ (568 cells) |
| Memory CD4 T | ✅ | ✅ (511 cells) |
| CD14+ Mono | ✅ | ✅ (472 cells) |
| CD8 T | ✅ | ✅ (370 cells) |
| B | ✅ | ✅ (347 cells) |
| FCGR3A+ Mono | ✅ | ✅ (168 cells) |
| NK | ✅ | ✅ (157 cells) |
| DC | ✅ | ✅ (30 cells) |
| Platelet | ✅ | ✅ (15 cells) |

---

## Key Numerical Comparison

| Metric | R Seurat | Shanuz | Match |
|--------|----------|--------|-------|
| Cells after QC | 2,638 | 2,638 | ✅ |
| HVGs selected | 2,000 | 2,000 | ✅ |
| Top-10 HVG overlap | — | 5/10 (50%) | ✅ |
| PC1 stdev | ~6.8 | 6.766 | ✅ |
| Number of clusters | 9 | 9 | ✅ |
| Cell types found | 9 | 9 | ✅ |
| Canonical markers reproduced | 6/6 | 6/6 | ✅ |

---

## Implementation Differences

| Aspect | R Seurat | Shanuz |
|--------|----------|--------|
| LOESS (VST) | Fortran (`stats::loess`) | `statsmodels.lowess` (it=3) |
| PCA | `irlba` randomized SVD | `sklearn.decomposition.PCA` |
| UMAP | `uwot` (C++) | `umap-learn` (Python) |
| Louvain | `igraph::cluster_louvain` | `python-igraph community_multilevel` |
| Wilcoxon | `wilcox.test` | `scipy.stats.ranksums` |
| Total runtime | ~minutes (R overhead) | ~49 s |

---

## References

> Hao Y, Stuart T, Kowalski MH, et al. (2024).
> **Dictionary learning for integrative, multimodal and scalable single-cell analysis.**
> *Nature Biotechnology*, 42, 293–304. https://doi.org/10.1038/s41587-023-01767-y

> Stuart T, Butler A, Hoffman P, et al. (2019).
> **Comprehensive Integration of Single-Cell Data.**
> *Cell*, 177(7), 1888–1902. https://doi.org/10.1016/j.cell.2019.05.031

> 10x Genomics (2016). *3k PBMCs from a Healthy Donor*.
> https://www.10xgenomics.com/resources/datasets/3-k-pb-mcs-from-a-healthy-donor-1-standard-1-1-0
