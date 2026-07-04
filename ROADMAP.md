# Shanuz Roadmap

This document tracks features planned for future releases, organized by milestone.
Each item includes the R Seurat equivalent, implementation notes, and dependencies
so any item can be picked up and scoped independently.

**What v0.1.0 already covers** (not listed below):  
LogNormalize · CLR · VST · ScaleData (+ covariate regression) · SCTransform ·
AddModuleScore · CellCycleScoring · PCA · UMAP · KNN/SNN · Louvain/Leiden ·
FindMarkers/FindAllMarkers (wilcox/t/LR/negbinom/roc) · JackStraw ·
DimPlot · FeaturePlot · VlnPlot · DotPlot · ElbowPlot · DoHeatmap · DimHeatmap ·
FeatureScatter · VariableFeaturePlot · RidgePlot ·
CITE-seq storage (Assay5 multi-layer) · AnnData interop ·
Spatial data structures (FOV/Centroids/Segmentation/Molecules)

---

## v0.2.0 — Batch Correction & Integration

> **Why first:** batch effects are unavoidable in real datasets; Harmony is the
> most widely-used, has a pip-installable Python package, and has a well-defined
> scope. CCA/RPCA follow naturally.

### Harmony integration
- **R:** `RunHarmony(obj, group.by.vars = "batch")` (via `harmony` package)
- **Python dep:** `harmonypy` (pip)
- **Plan:**
  1. `run_harmony(seurat, group_by, theta, lambda_, sigma, nclust, max_iter, random_seed)` → stores `"harmony"` in `obj.reductions`
  2. Input: PCA embeddings from `obj.reductions["pca"].cell_embeddings`
  3. Output: `DimReduc("harmony", embeddings)` — same shape as PCA
  4. Downstream: pass `reduction="harmony"` to `find_neighbors` / `run_umap`
- **Tests:** corrected embeddings have lower silhouette separation by batch than raw PCA

### CCA / RPCA integration (`IntegrateData` v4 API)
- **R:** `FindIntegrationAnchors(list, reduction="cca")` → `IntegrateData(anchors)`
- **Plan:**
  1. `find_integration_anchors(objects, dims, reduction, k_anchor, k_filter)` → `IntegrationAnchors`
  2. `integrate_data(anchors, dims)` → corrected `"integrated"` assay on merged object
  3. CCA: `scipy.linalg.svd` on cross-covariance; RPCA: project each dataset into
     shared PCA space, find mutual nearest neighbours (MNN) via `sklearn.neighbors`
- **Tests:** integrated embedding clusters by cell type, not by batch

### `IntegrateLayers` (Seurat v5 API)
- **R:** `IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca")`
- **Plan:** thin dispatch wrapper over the individual integration functions above;
  accepts `method` kwarg (`"harmony"`, `"cca"`, `"rpca"`)
- **Dep:** Harmony and CCA/RPCA functions above

---

## v0.3.0 — Reference Mapping & Label Transfer

> **Why:** label transfer from a curated atlas to a query dataset is a standard
> first-annotation step; all machinery (KNN, PCA, anchors) already exists.

### `FindTransferAnchors`
- **R:** `FindTransferAnchors(reference, query, dims = 1:30)`
- **Plan:**
  1. Project query into reference PCA space (`reference.reductions["pca"].loadings`)
  2. Find MNN between reference and projected query cells (sklearn `NearestNeighbors`)
  3. Score and filter anchors by consistency (same logic as CCA anchors above)
  4. Returns `TransferAnchors` object (stores anchor pairs + weights)
- **Dep:** CCA/RPCA from v0.2.0

### `TransferData`
- **R:** `TransferData(anchors, refdata = reference$celltype)`
- **Plan:**
  1. For each query cell, take a weighted average of reference labels/embeddings
     across its anchors
  2. Categorical labels → predicted label + prediction score per class
  3. Continuous (e.g. gene expression imputation) → weighted mean
- **Tests:** majority of cells get correct label when query == reference (self-transfer)

### `MapQuery` + `ProjectUMAP`
- **R:** `MapQuery(query, reference, refmodel)` then `ProjectUMAP(...)`
- **Plan:** compose `FindTransferAnchors` + `TransferData` + UMAP projection
  (transform-only mode of `umap-learn`, using reference's fitted UMAP model)
- **Note:** `umap-learn` exposes `UMAP.transform(query_embeddings)` — store the
  fitted model in `obj.reductions["umap"].misc["umap_model"]` after `run_umap`

---

## v0.4.0 — Weighted Nearest Neighbor (WNN)

> **Why:** shanuz already stores RNA + ADT assays; WNN is the natural joint
> analysis step for CITE-seq data and is well-scoped.

### `FindMultiModalNeighbors`
- **R:** `FindMultiModalNeighbors(obj, reduction.list = list("pca","apca"), dims.list = list(1:30, 1:18))`
- **Plan:**
  1. For each modality: compute KNN graph and within-modality prediction error
     (how well each cell's neighbours can predict its own embedding)
  2. Cell-specific modality weights: `w_RNA[i] = 1 - err_RNA[i] / (err_RNA[i] + err_ADT[i])`
     (Seurat's formula; approximated with local k-NN reconstruction error)
  3. Build weighted SNN: `SNN_wnn = w * SNN_RNA + (1-w) * SNN_ADT` per cell row
  4. Store as `"wknn"` and `"wsnn"` graphs; store per-cell weights in `meta_data`
- **Tests:** WNN clusters CBMC data closer to protein-defined ground truth than RNA alone

### WNN UMAP + clustering
- After `find_multi_modal_neighbors`, `find_clusters(reduction="wknn")` and
  `run_umap(graph="wsnn")` should work automatically since they accept any graph.
- **Plan:** verify `find_clusters` and `run_umap` accept `graph` kwarg and route
  correctly; add `reduction="wknn"` path to `run_umap` if not present.
- **Tutorial:** extend the CBMC CITE-seq tutorial (Tutorial 3) with WNN section.

---

## v0.5.0 — Additional Dimensionality Reductions

> Small self-contained additions; each is a single function wrapping a scipy/sklearn call.

### `run_tsne`
- **R:** `RunTSNE(obj, dims = 1:10)`
- **Python dep:** `scikit-learn` (`TSNE`) — already a dep
- **Plan:** mirrors `run_umap`; stores `DimReduc("tsne", embeddings)`

### `run_ica`
- **R:** `RunICA(obj, nics = 30)`
- **Python dep:** `sklearn.decomposition.FastICA`
- **Plan:** stores `DimReduc("ica", embeddings + loadings)`; `find_neighbors` and
  `run_umap` already accept `reduction="ica"`

### `run_spca` (supervised PCA)
- **R:** `RunSPCA(obj, assay, graph)`
- **Plan:** weighted PCA where gene loadings are regularised toward graph-defined
  gene modules; uses `scipy.linalg.svd` on `W @ X` where `W` is the gene-graph
  Laplacian smoothing matrix

### `glm_pca` (GLM-PCA)
- **R:** `RunGLMPCA` (via SeuratWrappers + glmpca)
- **Python dep:** `glmpca-py` (pip) or implement Poisson log-bilinear model directly
- **Plan:** wrap `glmpca.glmpca(Y, L=k)` → `DimReduc("glmpca", ...)`

---

## v0.6.0 — Pseudobulk DE & Advanced Marker Methods

### `AggregateExpression` (pseudobulk)
- **R:** `AggregateExpression(obj, group.by = c("celltype","donor"))`
- **Plan:**
  1. Sum raw counts per (`celltype`, `donor`) combination → pseudobulk count matrix
  2. Return a new `Shanuz` object (one "cell" per group) or a plain `pd.DataFrame`
  3. Intended input for `DESeq2`-style testing (see below)

### DESeq2-style pseudobulk DE
- **R:** `FindMarkers(obj, test.use = "DESeq2")` (via `DESeq2` R package)
- **Python dep:** `pydeseq2` (pip)
- **Plan:** add `"deseq2"` branch in `find_markers` dispatch; aggregates counts
  with `AggregateExpression` then runs `pydeseq2.DeseqDataSet`

### MAST
- **R:** `FindMarkers(obj, test.use = "MAST")`
- **Python dep:** `rpy2` (call R's MAST) or pure-Python hurdle model
- **Plan:** implement two-part hurdle model (logistic for detection + Gaussian for
  magnitude given detection) using `statsmodels`; no R dep required
- **Note:** `statsmodels` is already a dep (used by LR/negbinom tests)

### `FindConservedMarkers`
- **R:** `FindConservedMarkers(obj, ident.1, grouping.var)`
- **Plan:** run `find_markers` independently per group, then combine p-values
  (Fisher's method: `scipy.stats.combine_pvalues`) and report genes significant
  in all groups; returns combined + per-group stats columns

### `bimod` test (likelihood-ratio on bimodal model)
- **R:** `FindMarkers(obj, test.use = "bimod")`
- **Plan:** McDavid 2013 bimodal likelihood-ratio test; fits a mixture of a
  point mass at zero and a Gaussian; add `"bimod"` branch in `find_markers`

---

## v0.7.0 — Spatial Transcriptomics

> Data structures already exist (`FOV`, `Centroids`, `Segmentation`, `Molecules`).
> This milestone adds the loaders and spatial-aware analysis functions.
>
> **Partially delivered** (branch `feature/spatial-seurat-parity`): `load_xenium`,
> `load_visium`, `load_cosmx`; spatial-aware `from_anndata` (rebuilds `images`
> from `obsm['spatial']`); `get_tissue_coordinates`, `spatial_knn`,
> `nearest_neighbor_distance`, `local_neighborhood`, `build_niche_assay`;
> `image_dim_plot` / `image_feature_plot`; `composition_test`;
> `add_module_score(search=)`. Still open: `FindSpatiallyVariableFeatures`
> (Moran's I) and `load_merscope`.

### Loaders
| Function | Technology | File format |
|---|---|---|
| `load_10x_visium` | 10x Visium | `tissue_positions.csv`, `filtered_feature_bc_matrix/` |
| `load_xenium` | 10x Xenium | `transcripts.csv.gz`, `cell_feature_matrix/` |
| `load_cosmx` | NanoString CosMx | `*_exprMat_file.csv`, `*_fov_positions_file.csv` |
| `load_merscope` | Vizgen MERSCOPE | `cell_by_gene.csv`, `cell_metadata.csv` |

Each loader returns a `Shanuz` object with coordinates stored in an `FOV` slot.

### `FindSpatiallyVariableFeatures`
- **R:** `FindSpatiallyVariableFeatures(obj, method = "moransi")`
- **Plan:**
  1. **Moran's I** (primary): spatial autocorrelation statistic using the cell
     adjacency/distance weight matrix; compute with `esda.Moran` (`pysal` dep) or
     implement directly: `I = (N/W) * (z @ W @ z) / (z @ z)`
  2. **Markvariogram** (alternative): from the `Trendsceek` method
  3. Store results in `assay.var` (feature metadata), analogous to `FindVariableFeatures`
- **Python dep:** `libpysal` or pure NumPy implementation

### Spatial plots
| Function | R equivalent | Notes |
|---|---|---|
| `spatial_dim_plot` | `SpatialDimPlot` | cells on tissue image coloured by cluster/ident |
| `spatial_feature_plot` | `SpatialFeaturePlot` | gene expression on tissue |
| `image_dim_plot` | `ImageDimPlot` | for sub-cellular (Xenium/CosMx) data |
| `image_feature_plot` | `ImageFeaturePlot` | feature expression on segmentation image |

Implementation: scatter on tissue PNG background (`matplotlib.imshow` + `scatter`).

---

## v0.8.0 — Scale & Performance

### BPCells-style lazy/on-disk matrices
- **R:** `BPCells` package enables out-of-core analysis on millions of cells
- **Python dep:** `bpcells-python` (if available) or `zarr` / `h5py` lazy arrays
- **Plan:** `Assay5` layers can already hold any array-like object; extend
  `_get_data_matrix` to handle `zarr.Array` and `h5py.Dataset` transparently
  (they support numpy-style slicing). The key change is avoiding `.toarray()` /
  `.todense()` calls in hot paths — audit and gate these behind shape checks.

### `SketchData` (leverage-score sub-sampling)
- **R:** `SketchData(obj, ncells = 5000, method = "LeverageScore")`
- **Plan:**
  1. Compute per-cell leverage scores from PCA: `lev[i] = ||U[i,:]||² / k` where
     `U` is the left singular matrix of the scaled data
  2. Sample `ncells` cells with probability proportional to leverage scores
  3. Return subset `Shanuz` object; use `ProjectData` to extend results back to
     the full dataset
- **Tests:** sketched subset preserves cluster structure of the full dataset

### `ProjectData`
- **R:** `ProjectData(obj, reference = sketch, reduction = "pca")`
- **Plan:** project the full dataset into sketch's PCA/UMAP space using
  `sklearn.neighbors` or `umap-learn`'s `transform` method

---

## v0.9.0 — Specialized Assay Methods

### `HTODemux` / `MULTIseqDemux` (cell hashing)
- **R:** `HTODemux(obj)`, `MULTIseqDemux(obj)`
- **Plan:**
  - `hto_demux`: k-means (`k = n_hashtags + 1`) on CLR-normalised HTO counts;
    per-hashtag negative-binomial threshold to call positive cells; classify as
    singlet / doublet / negative
  - `multiseq_demux`: quantile-based bar-code classification
  - Both store `HTO_classification`, `HTO_maxID`, `nCount_HTO` in `meta_data`

### Mixscape (pooled CRISPR screen analysis)
- **R:** `RunMixscape(obj, target.gene.ident, nt.class.name)`
- **Plan:**
  1. `calc_perturbation_score`: per-cell log-fold-change vs negative-control cells
     projected onto a perturbation-response PCA (PRTB PCA)
  2. `run_mixscape`: Gaussian mixture model (2-component) per guide RNA;
     classify each cell as "perturbed" or "escaped"; uses `sklearn.mixture.GaussianMixture`
  3. `plot_mixscape`: violin of perturbation scores split by classification

---

## v0.10.0 — Package Infrastructure

### PyPI publication
- **Steps:**
  1. Add `build` and `twine` to `[dev]` extras in `pyproject.toml`
  2. Set up `__version__` import from `importlib.metadata` (replace the hard-coded string)
  3. Publish to TestPyPI first: `python -m build && twine upload --repository testpypi dist/*`
  4. Then PyPI: `twine upload dist/*`
- **Goal:** `pip install shanuz` works from PyPI

### GitHub Actions CI
- **File:** `.github/workflows/ci.yml`
- **Matrix:** Python 3.10, 3.11, 3.12 on ubuntu-latest
- **Jobs:**
  1. `lint` — `ruff check shanuz/`
  2. `test` — `pytest tests/ -v --cov=shanuz --cov-report=xml`
  3. `build` — `python -m build` (verify wheel builds cleanly)
- **Triggers:** push to `main`, all PRs

### Type annotations
- Add `from __future__ import annotations` to all modules (already done on some)
- Annotate all public function signatures (`mypy --strict` clean)
- Add `mypy` to the CI lint job

### Documentation site
- **Tool:** MkDocs + mkdocstrings (Material theme)
- **Structure:**
  ```
  docs/
    index.md          # Overview
    installation.md
    api/              # Auto-generated from docstrings
    tutorials/        # Symlinks to tutorials/*.md
    changelog.md
  ```
- **Deploy:** GitHub Actions → GitHub Pages on every push to `main`

### Changelog
- **File:** `CHANGELOG.md` at repo root
- Follow [Keep a Changelog](https://keepachangelog.com) format
- Populate retroactively for v0.1.0 from git log

---

## Dependency budget

Each milestone's new `pip` deps:

| Milestone | New deps |
|-----------|----------|
| v0.2.0 | `harmonypy` |
| v0.3.0 | *(none — uses sklearn already present)* |
| v0.4.0 | *(none)* |
| v0.5.0 | `glmpca-py` (optional) |
| v0.6.0 | `pydeseq2` (optional) |
| v0.7.0 | `libpysal` (optional, for Moran's I) |
| v0.8.0 | `zarr` (optional) |
| v0.9.0 | *(none — sklearn already present)* |
| v0.10.0 | `build`, `twine`, `mkdocs`, `mkdocstrings`, `ruff`, `mypy` (all dev-only) |

Optional deps go in a new `[spatial]`, `[integration]`, or `[all]` extra in
`pyproject.toml` so the base install stays lightweight.

---

## Priority order

If milestones are too large, these are the highest-value individual items:

1. **Harmony** (`v0.2.0`) — single function, well-scoped, huge real-world need
2. **WNN** (`v0.4.0`) — directly extends the existing CITE-seq tutorial
3. **GitHub Actions CI** (`v0.10.0`) — protects quality before any new feature lands
4. **`FindTransferAnchors` / `TransferData`** (`v0.3.0`) — enables atlas-based annotation
5. **Visium loader + `SpatialFeaturePlot`** (`v0.7.0`) — opens a new data modality
6. **`AggregateExpression` + DESeq2** (`v0.6.0`) — unlocks multi-sample DE
7. **`SketchData`** (`v0.8.0`) — enables million-cell datasets
