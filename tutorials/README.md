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
| Supervised PCA | `RunSPCA(pbmc, graph="wsnn")` | `run_spca(pbmc, graph="wsnn")` |
| GLM-PCA | `RunGLMPCA(pbmc, L=10)` | `glm_pca(pbmc, n_components=10)` — Poisson; `family="nb"` for negative binomial |
| Neighbors | `FindNeighbors(pbmc, dims)` | `find_neighbors(pbmc, dims, k_param)` |
| Cluster | `FindClusters(pbmc, resolution)` | `find_clusters(pbmc, resolution, algorithm)` |
| UMAP | `RunUMAP(pbmc, dims)` | `run_umap(pbmc, dims)` |
| Harmony integration | `RunHarmony(pbmc, "batch")` | `run_harmony(pbmc, "batch")` / `integrate_layers(pbmc, method="harmony", group_by="batch")` |
| CCA/RPCA anchors | `FindIntegrationAnchors(list, reduction="cca")` → `IntegrateData(anchors)` | `find_integration_anchors(objs, reduction="cca")` → `integrate_data(anchors)` |
| CCA/RPCA layers | `IntegrateLayers(obj, method=CCAIntegration)` | `integrate_layers(obj, method="cca", group_by="batch")` (or `"rpca"`) |
| Transfer anchors | `FindTransferAnchors(reference, query, reduction="pcaproject")` | `find_transfer_anchors(reference, query, reduction="pcaproject")` (or `"cca"`) |
| Transfer labels | `TransferData(anchors, refdata=reference$celltype)` | `transfer_data(anchors, refdata="celltype")` |
| Project into reference UMAP | `ProjectUMAP(query, reference, reduction.model="umap")` | `project_umap(query, reference)` |
| Map query (annotate + place) | `MapQuery(anchors, query, reference, refdata=list(id="celltype"))` | `map_query(anchors, refdata="celltype")` |
| Leverage scores | `LeverageScore(obj)` | `leverage_score(obj)` |
| Sketch a large dataset | `SketchData(obj, ncells=5000, method="LeverageScore")` | `sketch_data(obj, ncells=5000)` |
| Project sketch → full data | `ProjectData(obj, sketched.assay="sketch", reduction="pca")` | `project_data(full, sketch, refdata={"cluster_full": "seurat_clusters"})` |
| Matrix to disk (out-of-core) | `write_matrix_dir(mat, "counts.mat")` (BPCells) | `write_lazy_matrix(mat, "counts.mat")` |
| Open on-disk matrix | `open_matrix_dir("counts.mat")` (BPCells) | `open_lazy_matrix("counts.mat")` |
| Demultiplex hashtags | `HTODemux(obj, assay="HTO")` | `hto_demux(obj, assay="HTO")` |
| Demultiplex (MULTI-seq) | `MULTIseqDemux(obj, assay="HTO")` | `multiseq_demux(obj, assay="HTO")` |
| Perturbation signature | `CalcPerturbSig(obj, gd.class="gene", nt.cell.class="NT")` | `calc_perturb_sig(obj, labels="gene", nt_class="NT")` |
| Mixscape (CRISPR KO calls) | `RunMixscape(obj, labels="gene", nt.class.name="NT")` | `run_mixscape(obj, labels="gene", nt_class="NT")` |
| Mixscape LDA (guide separation) | `MixscapeLDA(obj, labels="gene", nt.label="NT")` | `mixscape_lda(obj, labels="gene", nt_class="NT")` |
| Mixscape perturbation score | `PlotPerturbScore(obj, target.gene.ident="IFNGR2")` | `plot_perturb_score(obj, target_gene_ident="IFNGR2")` |
| Mixscape DE heatmap | `MixscapeHeatmap(obj, ident.1="NT", ident.2="IFNGR2 KO")` | `mixscape_heatmap(obj, ident_1="NT", ident_2="IFNGR2 KO")` |
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

### Beyond PCA — supervised PCA and GLM-PCA

Two reductions that answer questions PCA cannot. Both store a `DimReduc` exactly
as `run_pca` does, so `find_neighbors`, `find_clusters` and `run_umap` take them
with a `reduction=` argument and nothing else changes.

```python
from shanuz import run_spca, glm_pca, find_neighbors, run_umap

# Supervised PCA: the gene axes that best explain a graph you already trust.
find_multi_modal_neighbors(cbmc, reduction_list=["pca", "apca"])   # builds "wsnn"
run_spca(cbmc, graph="wsnn", npcs=50)
run_umap(cbmc, reduction="spca", dims=range(30))

# GLM-PCA: a Poisson fit straight on the raw counts — no log, no pseudocount.
glm_pca(pbmc, n_components=10)
find_neighbors(pbmc, reduction="glmpca", dims=range(10))
print(pbmc.reductions["glmpca"].misc["converged"])

# Negative binomial for over-dispersed counts; θ is estimated by ML.
glm_pca(pbmc, n_components=10, family="nb")
print(pbmc.reductions["glmpca"].misc["theta"])          # the fitted dispersion
```

**`run_spca`** maximises `vᵀXᵀGXv` where PCA maximises `vᵀXᵀXv` — the same problem
with the identity swapped for a cell-cell graph `G`. So it finds the gene
directions that reproduce a neighbourhood structure you have already decided is
the right one (typically the WNN graph, which knows about protein as well as RNA).
Hand it `G = I` and you get PCA back exactly. Its value is that the result is a
plain linear map from genes to components, so a query dataset can be projected
into a reference's space with one matrix multiply — which is why Azimuth maps onto
sPCA and not PCA.

**`glm_pca`** models counts as counts. The standard pipeline log-normalises and
then runs PCA, which quietly assumes constant-variance Gaussian data; a gene
averaging 0.1 UMIs and one averaging 100 do not have remotely comparable noise,
and the pseudocount you add to survive `log(0)` distorts exactly the
low-expression genes where most of the zeros are. GLM-PCA drops the transform and
fits `log μ = a[g] + o[c] + U·Vᵀ` with a Poisson likelihood, holding the log
library size as a fixed offset so sequencing depth is a known quantity rather than
something the factors have to spend themselves discovering.

Real UMI counts are usually noisier than Poisson allows — the same gene in two
copies of one cell state still varies more than the mean — and a handful of such
genes will dominate a Poisson fit. `family="nb"` swaps in a negative binomial,
`Var = μ + μ²/θ`, with a single shared dispersion `θ` estimated by maximum
likelihood alongside the factors (pass `optimize_theta=False` and a `theta=` to
pin it). As `θ → ∞` it collapses back onto Poisson, so NB never fits worse — only
more forgivingly. The fitted `θ` lands in `misc["theta"]`.

> Two things to know. `glm_pca` fits densely in genes × cells, so pass a few
> thousand variable features, as you would to `run_pca`. And check
> `misc["converged"]`; if it is `False`, or `misc["deviance"]` is still falling
> steeply at the end, raise `max_iter`. (With `family="nb"` and `θ` being
> estimated, the deviance is re-scaled as `θ` moves, so it is only strictly
> monotone when you pin `θ` with `optimize_theta=False`.)

### Integrating datasets — Harmony, CCA, RPCA

Two batches of the same tissue rarely line up: a shared cell type sits in a
different place in each dataset's PCA, so cells cluster by batch before they
cluster by biology. Shanuz offers two remedies.

**Harmony** (`run_harmony`) corrects an embedding you already have — it takes the
joint PCA and iteratively pulls the batches together while keeping cell types
apart. It is fast and needs only a `group_by` column.

**Anchor integration** (Seurat's CCA/RPCA) works from scratch, without assuming
the datasets even share a coordinate system. It builds a *shared* space for a
pair of datasets — by canonical correlation (`reduction="cca"`, the SVD of the
cross-covariance) or reciprocal PCA (`reduction="rpca"`) — then keeps only the
*mutual* nearest neighbours across datasets as **anchors**: cell *i* here and
cell *j* there, each among the other's closest matches, are almost certainly the
same state seen twice. Anchors are scored for neighbourhood consistency, filtered
against the raw expression, and then used to pull every query dataset onto the
reference.

```python
from shanuz import (
    find_integration_anchors, integrate_data, integrate_layers,
    run_harmony, scale_data, run_pca, find_clusters, run_umap,
)

# --- Harmony: correct an existing joint PCA in place -------------------------
run_harmony(pbmc, group_by="batch")             # stores reductions["harmony"]
run_umap(pbmc, reduction="harmony", dims=range(30))

# --- CCA/RPCA anchors: a list of per-batch objects → a corrected assay -------
anchors = find_integration_anchors([ref, query], reduction="cca")   # or "rpca"
merged = integrate_data(anchors)                # active assay is now "integrated"
scale_data(merged)
run_pca(merged)                                 # clusters by cell type, not batch
find_clusters(merged, resolution=0.5)

# --- One-call Seurat v5 path: split one object by batch, integrate, embed ----
integrate_layers(pbmc, method="cca", group_by="batch", new_reduction="integrated")
run_umap(pbmc, reduction="integrated", dims=range(30))
```

**`find_integration_anchors`** is *reference-based*: `objects[reference]` (index
0 by default) is the anchor every other dataset is corrected onto. **CCA** shines
when the datasets share structure but differ globally (cross-species, cross-
technology); **RPCA** is stricter and faster, a better fit when the batches are
already similar or very large. Both return the same `IntegrationAnchors`, which
is also what v0.3.0's reference mapping is built to consume.

> The correction is applied to the log-normalised `data` of the shared anchor
> features and stored as an `"integrated"` assay — so `scale_data` + `run_pca`
> on it is the natural next step. The reference dataset is left untouched.

### Reference mapping — annotating a query from an atlas

Integration mixes several datasets into one shared space. Reference mapping is
the asymmetric cousin: you keep an annotated atlas fixed and *borrow* its labels
for a new, unlabelled dataset. The anchor machinery is the same — build a shared
space, find mutual nearest neighbours, score and filter them — but the reference
is never moved and the anchors carry information *reference → query*.

```python
from shanuz import find_transfer_anchors, transfer_data

# reference: annotated atlas (has a "celltype" column). query: new, unlabelled.
# Both normalized + find_variable_features + scale_data, as usual.
anchors = find_transfer_anchors(reference, query, reduction="pcaproject")

# Classification: predict the query's cell types from the reference labels.
pred = transfer_data(anchors, refdata="celltype")
query.add_meta_data(pred["predicted.id"], col_name="predicted.celltype")
query.add_meta_data(pred["prediction.score.max"], col_name="prediction.score")
# pred also has one prediction.score.<class> column per reference class.

# Imputation: carry reference expression (features × ref-cells) onto the query.
ref_expr = reference.get_assay().layer_data("data", features=["CD3D", "MS4A1"])
imputed = transfer_data(anchors, refdata=ref_expr, refdata_features=["CD3D", "MS4A1"])
```

**`reduction="pcaproject"`** (the default) projects the query through the
*reference's* PCA loadings — the principal axes are learned once on the reference
and the query is pushed through the same map. Because those axes never saw the
query, batch-specific structure the reference lacks simply lands nowhere, which
is what makes projection robust for annotation. **`reduction="cca"`** learns a
joint space instead, for the harder cross-modality / cross-species cases.

> `transfer_data` weights each query cell's anchors with the same
> distance-weighted, score-scaled Gaussian kernel `integrate_data` uses, so a
> query cell surrounded by confident, consistent anchors of one type gets a
> high `prediction.score.max`; an ambiguous one gets a low score you can filter
> on.

#### Placing the query in the reference UMAP

Annotating the query is half the job; the other half is *seeing* it on the atlas
you already know how to read. `project_umap` runs the reference's **fitted** UMAP
model in transform-only mode, so the query cells land in the reference's existing
embedding rather than in a fresh, unrelated one. `map_query` composes the whole
workflow — transfer the labels *and* project the UMAP — in a single call.

```python
from shanuz import run_pca, run_umap, find_transfer_anchors, map_query, project_umap

# The reference needs a fitted PCA + UMAP; run_umap stashes the umap-learn model
# in reference.reductions["umap"].misc["umap_model"] for transform-only projection.
run_pca(reference)
run_umap(reference)                       # embeds from "pca", keeps the model

anchors = find_transfer_anchors(reference, query, reduction="pcaproject")

# One call: transfer_data writes predicted.id / prediction.score.* onto the
# query's metadata, and project_umap places it in the reference UMAP.
pred = map_query(anchors, refdata="celltype")
query.reductions["ref.umap"]              # the query, on the reference's UMAP

# Or just the projection, without label transfer:
project_umap(query, reference)            # -> query.reductions["ref.umap"]
```

`project_umap` is itself a two-step map: it projects the query through the
*reference's* PCA loadings into the reference's PC space (the same "project into a
space the query never helped define" logic as `pcaproject` above), then runs the
reference's UMAP model's `.transform`. A query cell that resembles reference
T cells is pinned near the reference's T-cell island and optimised to sit there —
so the query overlays the atlas, batch block and all.

### Sketching — analysing a million cells through a small subset

A huge atlas is mostly redundant: the common states are thousands of near-identical
cells stacked on top of each other, while the rare states — usually the interesting
ones — are a handful of points each. Sketching picks a small, information-dense
subset, does the expensive clustering / UMAP on *that*, and projects the answers
back onto every cell. The subset is drawn by **leverage**, not uniformly: a cell in
a dense redundant cloud has low leverage (drop it and nothing changes), a cell in a
sparse distinctive corner has high leverage (it is the only evidence that corner
exists) — so leverage sampling keeps the rare states a uniform sample would lose.

```python
from shanuz import sketch_data, project_data, run_pca, run_umap
from shanuz.neighbors import find_neighbors
from shanuz.clustering import find_clusters

# full: a large, normalized + scaled object. Draw a leverage-weighted sketch.
sketch = sketch_data(full, ncells=5000)   # a standalone object; assay -> "sketch"
# full.meta_data["leverage.score"] now holds the per-cell scores.

# Do the heavy analysis on the small sketch.
run_pca(sketch, n_pcs=30)
find_neighbors(sketch, dims=range(30))
find_clusters(sketch, resolution=0.8)     # -> sketch.meta_data["seurat_clusters"]
run_umap(sketch, dims=range(30))          # keeps the umap-learn model

# Extend the sketch's PCA, UMAP and cluster labels back to *every* cell.
project_data(full, sketch, refdata={"cluster_full": "seurat_clusters"})
full.reductions["pca.full"]               # every cell in the sketch's PC space
full.reductions["ref.umap"]               # every cell on the sketch's UMAP
full.meta_data["cluster_full"]            # every cell's transferred cluster
```

`leverage_score` computes the scores with a sparse **CountSketch** rather than a
full SVD of the data — the whole reason sketching scales — so it stays cheap on
millions of cells. `project_data` is the mirror image of reference mapping: the
sketch plays the role of the reference, and each full cell is pushed through the
sketch's PCA loadings (and its fitted UMAP model) exactly as `project_umap` places
a query on an atlas, with `find_transfer_anchors` + `transfer_data` carrying the
cluster labels across.

> `leverage_score` and `sketch_data` default to `nsketch=5000`; the CountSketch is
> a faithful embedding once it has more rows than roughly the squared feature
> count, which the default comfortably clears for the few-thousand variable
> features a sketch runs on.

### Lazy on-disk matrices — keeping a million cells out of RAM

Sketching shrinks *how many cells you analyse*; a **lazy matrix** shrinks *how much
of the matrix is in memory at once*. A dense million-by-twenty-thousand `float64`
matrix is 160 GB; even the sparse counts leave no room for the copies each step
makes. Seurat's answer is the `BPCells` package, which keeps the matrix on disk and
streams over it. `LazyMatrix` is the shanuz analogue, built on NumPy's
memory-mapping — **no new dependency**.

A matrix is written to a directory as the three memory-mapped arrays of a
compressed-sparse-column matrix (scipy's `csc_matrix` layout). Opening it maps
those arrays without reading them; a slice pulls only the touched cells off disk
and hands back an ordinary `scipy.sparse` block, so it drops straight into an assay
layer:

```python
from shanuz import write_lazy_matrix, open_lazy_matrix

assay = obj.get_assay()

# Persist the counts layer out-of-core, then map it back in and swap it in place.
write_lazy_matrix(assay.layers["counts"], "counts.mat")
lazy = open_lazy_matrix("counts.mat")
assay.set_layer_data("counts", lazy)          # a LazyMatrix is a valid layer

# Slicing reads only the selected cells' non-zeros off disk...
block = assay.layer_data("counts", cells=obj.cell_names()[:1000])

# ...and reductions stream in a single pass without materialising the matrix.
per_cell = lazy.sum(axis=0)                    # nCount, one pass over the store
per_gene = lazy.mean(axis=1)                   # per-feature mean

# col_blocks is the streaming primitive: walk a million cells at bounded RAM.
for start, stop, chunk in lazy.col_blocks(block_size=50_000):
    ...                                        # chunk is a csc_matrix of those cells
```

`LazyMatrix` stores **columns** (cells) contiguously, because the operations that
dominate at scale — sketching, cell subsetting, per-cell normalisation — select
cells, and CSC makes reading an arbitrary set of columns cost only their own
non-zeros. `as_dense(lazy)` / `np.asarray(lazy)` still materialise the whole thing
when you genuinely need it — the escape hatch you keep for the small datasets and
avoid on the million-cell path.

### Cell hashing — demultiplexing pooled samples

Cell Hashing (Stoeckius et al.) tags each *sample* with a distinct
antibody-oligo **hashtag** before pooling the samples on one lane. Every droplet
then carries a little vector of hashtag counts saying which sample it came from —
and droplets that caught two cells carry two tags. `hto_demux` (Seurat's
`HTODemux`) turns that hashtag matrix back into per-cell calls: **singlet** (one
tag), **doublet** (two or more), or **negative** (none).

The hashtags live in their own assay alongside the RNA. `hto_demux` learns each
tag's positive cutoff from the data — it CLR-normalizes the counts, k-means
clusters the cells (`k = n_hashtags + 1`), fits a negative binomial to each tag's
*background* (the cluster where it is least expressed), and thresholds at the
0.99 quantile:

```python
from shanuz import hto_demux

# `obj` has an "HTO" assay of hashtag counts (features = hashtags, columns = cells).
hto_demux(obj, assay="HTO")                    # positive_quantile=0.99 by default

# Per-cell results land in meta_data under the Seurat column names:
obj.meta_data["HTO_classification.global"]     # "Singlet" / "Doublet" / "Negative"
obj.meta_data["HTO_maxID"]                      # the top hashtag per cell
obj.meta_data["HTO_classification"]            # tag name (singlet) or "HTOa_HTOb" (doublet)
obj.meta_data["hash.ID"]                        # tag name / "Doublet" / "Negative"

# hash.ID is also set as the active identity, so keeping only the clean singlets
# of one sample is a one-liner:
singlets = obj.meta_data["HTO_classification.global"] == "Singlet"
sample3 = obj.subset(cells=obj.meta_data.index[
    singlets & (obj.meta_data["hash.ID"] == "HTO-3")
].tolist())

# The learned per-hashtag cutoffs are kept for inspection:
obj.misc["hto_demux"]["HTO"]["cutoffs"]         # {"HTO-1": 6.0, "HTO-2": 5.0, ...}
```

By default `hto_demux` CLR-normalizes internally; if you already ran
`normalize_data(obj, normalization_method="CLR", margin=2, assay="HTO")`, pass
`normalize=False` to reuse that `data` layer. The negative binomial — not a fixed
threshold — is what makes the call robust to each antibody's own staining
background and each run's own depth.

**MULTI-seq — the other demultiplexer.** MULTI-seq (McGinnis et al.) uses
lipid-anchored barcodes instead of antibodies, and Seurat's `MULTIseqDemux` calls
it a different way: rather than a background fit, it reads each barcode's cutoff
straight off the *shape* of its distribution. `multiseq_demux` is a drop-in
alternative that often disagrees with `hto_demux` at the margins, so it is handy
as a second opinion:

```python
from shanuz import multiseq_demux

# Same "HTO" assay of barcode counts. For each barcode a Gaussian KDE exposes its
# background and positive modes; the cutoff sits a fraction `quantile` between them.
multiseq_demux(obj, assay="HTO", quantile=0.7)     # 0.7 is the Seurat default

obj.meta_data["MULTI_ID"]                          # barcode name / "Doublet" / "Negative"
obj.meta_data["MULTI_classification"]              # a character copy of MULTI_ID
obj.misc["multiseq_demux"]["HTO"]["thresholds"]    # learned per-barcode cutoffs

# Don't want to guess `quantile`? Let it sweep for the value that maximises the
# singlet rate, peeling off negatives and re-thresholding until it settles:
multiseq_demux(obj, assay="HTO", autothresh=True)  # ignores `quantile`
```

`MULTI_ID` is set as the active identity, so subsetting one sample's singlets works
exactly as with `hash.ID` above. As with `hto_demux`, pass `normalize=False` to
reuse an existing CLR `data` layer instead of recomputing it.

### Pooled CRISPR screens — Mixscape

In a pooled CRISPR screen every cell carries a guide RNA, but carrying a guide is
not the same as being perturbed: some cells escape the knockout and look just like
controls. Mixscape (Papalexi, Mimitou et al.) separates the true knockouts (KO)
from those non-perturbed escapers (NP) so downstream analysis runs on genuinely
perturbed cells. It is a two-step workflow — build a per-cell **perturbation
signature**, then classify against it — mirroring Seurat's `CalcPerturbSig` +
`RunMixscape`. It expects a guide-assignment column (each cell's target gene, with
the controls labelled `"NT"`) and a computed reduction (`pca`):

```python
from shanuz import calc_perturb_sig, run_mixscape

# 1. Local perturbation signature: subtract each cell's 20 nearest NT controls
#    (in PCA space) to cancel cell-cycle / depth / batch variation. Stored as a
#    new "PRTB" assay. `split_by` keeps neighbours within a replicate.
calc_perturb_sig(obj, assay="RNA", labels="gene", nt_class="NT",
                 reduction="pca", ndims=15, num_neighbors=20)   # → obj.assays["PRTB"]

# 2. Per guide: DE vs NT picks the response genes, then an iterative 2-component
#    Gaussian mixture over the perturbation score splits KO from NP.
run_mixscape(obj, assay="PRTB", labels="gene", nt_class="NT", de_assay="RNA")

obj.meta_data["mixscape_class"]           # "IFNGR2 KO" / "IFNGR2 NP" / "NT"
obj.meta_data["mixscape_class.global"]    # "KO" / "NP" / "NT"
obj.meta_data["mixscape_class_p_ko"]      # KO posterior per guide cell (NaN for NT)
obj.misc["mixscape"]["PRTB"]["genes"]     # per-gene DE-gene / iteration / KO counts
```

`mixscape_class` is set as the active identity, so `obj.subset(idents="IFNGR2 KO")`
pulls just the confirmed knockouts. A guide with too few cells, or too few DE genes
to show a phenotype (`min_de_genes`, default 5), has all its cells called NP. For a
knock-down rather than a knockout screen, pass `prtb_type="KD"` — the class labels
and the `mixscape_class_p_kd` posterior column follow the name.

#### Separating the guide populations — `mixscape_lda`

`run_mixscape` asks which *cells* are perturbed. The complementary question is how
the guide *populations* differ from each other and from control, and `mixscape_lda`
(Seurat's `MixscapeLDA`) answers it with a single supervised map on which each guide
class forms its own cloud. It needs only the `PRTB` signature — it groups cells by
their raw guide label, so the KO/NP calls are not used and it can follow
`calc_perturb_sig` directly:

```python
from shanuz import mixscape_lda

# Per guide: DE vs NT picks its response genes, a PCA is fit on that guide's cells
# plus the NT cells, and every cell is projected onto the guide's npcs-dim subspace.
# The blocks are concatenated and one LDA is fit with the guide as the class.
mixscape_lda(obj, assay="PRTB", labels="gene", nt_class="NT", npcs=10)

obj.reductions["lda"]                      # n_classes - 1 discriminant dimensions
obj.reductions["lda"].misc["genes_used"]   # guides that contributed a block
obj.meta_data["lda_assignments"]           # predicted guide class per cell
obj.meta_data["LDAP_IFNGR2"]               # posterior for that class, per cell
```

Plot it like any other reduction — `dim_plot(obj, reduction="lda", group_by="gene")`.
A guide needs at least `npcs + 1` DE genes to contribute a block (it cannot support
`npcs` components otherwise) and is skipped if it falls short; if no guide clears the
bar, `mixscape_lda` raises and you should lower `npcs`. Because the grouping is by
guide rather than by mixscape class, NP escapers stay in their guide's group — and,
being control-like, they generally land on top of the NT cloud, which is itself a
useful read on how much of a guide escaped.

#### Checking the calls — `plot_perturb_score` and `mixscape_heatmap`

Mixscape's KO/NP split is a threshold on one number per cell: the perturbation
score, each cell's projection onto its guide's response axis. `plot_perturb_score`
(Seurat's `PlotPerturbScore`) draws that axis, overlaying the NT control density
against the guide's own — the single most useful check that a guide worked.

```python
from shanuz import plot_perturb_score, mixscape_heatmap

fig = plot_perturb_score(obj, target_gene_ident="IFNGR2")
fig.savefig("perturb_score.png", dpi=150, bbox_inches="tight")
```

A guide with a real effect is **bimodal**: one lobe sitting on the NT curve (the
escapers) and one shifted away from it (the knockouts) — exactly the structure the
mixture model is asked to find. A guide that simply did not work is one curve on
top of the controls, and no threshold will rescue it. By default the curves are
coloured by `mixscape_class`, so you see where mixscape actually drew the line;
pass `before_mixscape=True` for the raw guide label instead, which is the view you
would have had without mixscape at all. For a screen spanning several cell types,
`split_by="celltype"` facets it.

`mixscape_heatmap` (Seurat's `MixscapeHeatmap`) then shows the genes underneath
that score — the DE genes between two classes, with every cell ordered by its
knockout probability:

```python
fig = mixscape_heatmap(obj, ident_1="NT", ident_2="IFNGR2 KO",
                       max_genes=20, balanced=True)
```

Read with the class colour bar along the top, a clean screen shows the expression
block turning on in step with the probability, the low-probability escapers at one
end still looking like control. `ident_1` / `ident_2` are `mixscape_class` levels
(`"NT"`, `"IFNGR2 KO"`, `"IFNGR2 NP"`), which `run_mixscape` also leaves as the
active identity. `balanced=True` takes up to `max_genes` from each direction of
the fold change rather than only the up-regulated ones, and `max_cells_group`
downsamples each class for a large screen.

Both plots read the perturbation score that `run_mixscape` stores in
`obj.misc["mixscape"]["PRTB"]["genes"]`, so they need `run_mixscape` to have run —
unlike `mixscape_lda`, which needs only `calc_perturb_sig`. A guide whose cells
were called NP without a mixture fit (too few cells, or too few DE genes) has no
score axis, and `plot_perturb_score` says so rather than drawing an empty panel.

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
| `SpatialDimPlot(obj, group.by)` | `spatial_dim_plot(obj, group_by)` — spots over the H&E image |
| `SpatialFeaturePlot(obj, features)` | `spatial_feature_plot(obj, feature)` |

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
| `FindSpatiallyVariableFeatures(obj, method="moransi")` | `find_spatially_variable_features(obj, k=10)` — Moran's I |
| `FindSpatiallyVariableFeatures(obj, method="markvariogram", r.metric=5)` | `find_spatially_variable_features(obj, method="markvariogram", r_metric=5)` — `r_metric` is in cell spacings, not pixels |
| *(hand-rolled Fisher + `p.adjust`)* | `composition_test(obj, group_by, split_by)` |
| `GetImage(obj[["slice1"]])` | `obj.images["spatial"].get_image()` — the Visium H&E image |
| `ScaleFactors(obj[["slice1"]])` | `obj.images["spatial"].scale_factors` |
| `SpatialDimPlot(obj)` | `spatial_dim_plot(obj)` |
| `SpatialFeaturePlot(obj, features)` | `spatial_feature_plot(obj, feature)` |

> **Plot output:** R renders to the graphics device automatically. Shanuz functions return a
> `matplotlib.Figure` — call `fig.savefig("out.png")` to save or display inline in Jupyter.

#### Visium tissue images

`load_visium` reads the H&E tissue image and `scalefactors_json.json` by default,
giving each image slot a `VisiumV2` (Seurat v5's class) instead of a bare `FOV`:

```python
from shanuz import load_visium

obj = load_visium("visium_out/")            # image=True by default
fov = obj.images["spatial"]

fov.get_image().shape                       # (H, W, 3) — tissue_hires_image.png
fov.scale_factors.spot                      # spot diameter, full-resolution pixels
fov.radius()                                # spot radius, full-resolution pixels

# Spot coordinates stay in FULL-RESOLUTION pixels, so spatial_knn /
# nearest_neighbor_distance / Moran's I keep measuring real distances.
# Convert to the stored image's pixel space only when you draw:
xy = fov.scale_coordinates()                # x, y scaled onto the PNG
r  = fov.spot_radius()                      # matching spot radius, in image pixels
```

Pass `image_resolution="lowres"` for the smaller PNG, `image=False` to skip it
entirely, or `filter_by_tissue=True` to keep only spots with `in_tissue == 1`. A
bundle with no PNG still loads — you get a plain `FOV`, exactly as before.

#### Plotting spots on the tissue

`spatial_dim_plot` / `spatial_feature_plot` (`SpatialDimPlot` / `SpatialFeaturePlot`)
draw the H&E photo and overlay the spots on top of it — the scaling above happens
for you:

```python
from shanuz import spatial_dim_plot, spatial_feature_plot

fig = spatial_dim_plot(obj, group_by="seurat_clusters")
fig = spatial_feature_plot(obj, "Gad1", image_alpha=0.4)   # fade the tissue
fig.savefig("visium.png", dpi=150, bbox_inches="tight")
```

Spots are drawn at their **true diameter** (from `spot_diameter_fullres`), not as
fixed-size points, so they stay registered against the tissue at any zoom.
`pt_size_factor=` (default 1.6, as in Seurat) scales them; `crop=False` shows the
whole slide instead of zooming to the spots; `resolution="lowres"` draws the
smaller PNG.

Neither function needs an image to work. Plot an object loaded with `image=False`
and you get a plain scatter of the same spots — useful for Xenium/CosMx, or for a
Visium bundle whose PNG is missing.

#### Finding spatially variable genes

`find_spatially_variable_features` (`FindSpatiallyVariableFeatures`) ranks genes by
how strongly their expression is organised in space. Both of R's methods are here:

```python
from shanuz import find_spatially_variable_features

# Moran's I — is this gene autocorrelated at all? Comes with a p-value.
svf = find_spatially_variable_features(obj, k=10)
svf.head()          # moransi, moransi_pval, moransi_padj, moransi_rank

# Mark variogram — how much of the gene's variance has decayed by distance r?
mv = find_spatially_variable_features(obj, method="markvariogram", r_metric=5)
mv.head()           # markvariogram, markvariogram_rank

fig = spatial_feature_plot(obj, mv.index[0])   # the most spatially variable gene
```

Both return a table sorted so that **rank 1 is the most spatially variable**, and
both also write their columns into the assay's feature metadata, the way
`find_variable_features` does. On a large panel, pass `features=` to score only the
variable genes — the mark variogram is the heavier of the two.

The two statistics ask different questions. Moran's I gives one number for the
whole slide and a p-value with it. The mark variogram is read *at a distance*:
`markvariogram` is the average squared expression difference between cells about
`r_metric` apart, divided by the gene's own variance. So **≈ 1 means no spatial
structure** (two cells that far apart differ as much as any two cells picked at
random) and **below 1 means they still resemble each other**. There is no p-value —
the variogram has no closed-form null, and R does not offer one either.

> **`r_metric` is not in R's units.** R passes `r.metric` straight through to
> `spatstat` in raw coordinate units, so the same script answers differently on a
> slide measured in pixels and one in microns. Here `r_metric` is measured in
> nearest-neighbour spacings, so the default of 5 means *"five cells apart"* on any
> slide. Widen `bandwidth=` (also in spacings) if the tissue is sparse and too few
> cell pairs land near `r_metric` to average over.
