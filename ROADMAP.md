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

### Harmony integration — ✅ delivered
- Implemented in `shanuz/integration.py` as `run_harmony(...)`; stores a
  `DimReduc("harmony")` and is verified to lower per-batch silhouette while
  preserving cell-type separation (`tests/test_integration.py`). Enable with
  `pip install shanuz[integration]` (adds the `harmonypy` dep).
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

### `IntegrateLayers` (Seurat v5 API) — ✅ harmony method delivered
- Implemented as `integrate_layers(obj, method="harmony", group_by=...)` in
  `shanuz/integration.py`. The `"cca"` / `"rpca"` methods raise
  `NotImplementedError` pending the CCA/RPCA anchor work above.
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

### `FindMultiModalNeighbors` — ✅ delivered
- Implemented as `find_multi_modal_neighbors(...)` in `shanuz/multimodal.py`.
  Per-cell modality weights use the sanctioned scale-invariant approximation
  (own- vs cross-modality reconstruction distance); stores `wknn`/`wsnn` graphs
  and `<assay>.weight` columns. Verified on synthetic complementary-modality
  data to recover structure RNA alone cannot (`tests/test_multimodal_wnn.py`).
- **R:** `FindMultiModalNeighbors(obj, reduction.list = list("pca","apca"), dims.list = list(1:30, 1:18))`
- **Plan:**
  1. For each modality: compute KNN graph and within-modality prediction error
     (how well each cell's neighbours can predict its own embedding)
  2. Cell-specific modality weights: `w_RNA[i] = 1 - err_RNA[i] / (err_RNA[i] + err_ADT[i])`
     (Seurat's formula; approximated with local k-NN reconstruction error)
  3. Build weighted SNN: `SNN_wnn = w * SNN_RNA + (1-w) * SNN_ADT` per cell row
  4. Store as `"wknn"` and `"wsnn"` graphs; store per-cell weights in `meta_data`
- **Tests:** WNN clusters CBMC data closer to protein-defined ground truth than RNA alone

### WNN UMAP + clustering — ✅ delivered (tutorial pending)
- `find_clusters(graph_name="wsnn")` already routed correctly; `run_umap` now
  accepts a `graph=` kwarg that embeds a precomputed graph via UMAP's
  `simplicial_set_embedding` (`tests/test_reductions_extra.py`).
- **Still open:** extend the CBMC CITE-seq tutorial (Tutorial 3) with a WNN section.

---

## v0.5.0 — Additional Dimensionality Reductions

> Small self-contained additions; each is a single function wrapping a scipy/sklearn call.

### `run_tsne` — ✅ delivered
- Implemented in `shanuz/reduction.py` (`run_tsne`), mirrors `run_umap`; stores
  `DimReduc("tsne")`. Tested in `tests/test_reductions_extra.py`.
- **R:** `RunTSNE(obj, dims = 1:10)`
- **Python dep:** `scikit-learn` (`TSNE`) — already a dep

### `run_ica` — ✅ delivered
- Implemented in `shanuz/reduction.py` (`run_ica`); stores `DimReduc("ica",
  embeddings + loadings)`; `find_neighbors`/`run_umap` accept `reduction="ica"`.
  Tested in `tests/test_reductions_extra.py`.
- **R:** `RunICA(obj, nics = 30)`
- **Python dep:** `sklearn.decomposition.FastICA`

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

### `AggregateExpression` (pseudobulk) — ✅ delivered
- Implemented as `aggregate_expression(...)` in `shanuz/aggregate.py`. Sums raw
  counts per group via a single sparse `counts @ indicator` matmul; `group_by`
  accepts one or more metadata columns (joined with `"_"`, as Seurat does) or
  `"ident"`. Returns a features×groups `pd.DataFrame` (a `dict` for multiple
  assays), or a `Shanuz` object with one "cell" per group when
  `return_object=True` (`tests/test_pseudobulk_conserved.py`).
- **R:** `AggregateExpression(obj, group.by = c("celltype","donor"))`
- Intended input for `DESeq2`-style testing (see below).

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

### `FindConservedMarkers` — ✅ delivered
- Implemented as `find_conserved_markers(...)` in `shanuz/markers.py`. Runs
  `find_markers` independently within each level of `grouping_var`, keeps genes
  that are markers in *every* level, and combines their per-level p-values with
  Fisher's method (`scipy.stats.combine_pvalues`). Output has per-level prefixed
  stats plus `max_pval` and `combined_p_val` (sorted by the latter); levels
  lacking a comparison group are skipped with a warning
  (`tests/test_pseudobulk_conserved.py`).
- **R:** `FindConservedMarkers(obj, ident.1, grouping.var)`

### `bimod` test (likelihood-ratio on bimodal model)
- **R:** `FindMarkers(obj, test.use = "bimod")`
- **Plan:** McDavid 2013 bimodal likelihood-ratio test; fits a mixture of a
  point mass at zero and a Gaussian; add `"bimod"` branch in `find_markers`

---

## v0.7.0 — Spatial Transcriptomics

> **Largely delivered** on branch `feature/spatial-seurat-parity`. The data
> structures (`FOV`, `Centroids`, `Segmentation`, `Molecules`) plus these
> loaders and analysis functions are done and validated end-to-end against R
> Seurat in [Tutorial 5](tutorials/xenium_spatial_tutorial.md) (deterministic
> anchors match to 8 significant figures):
>
> - **Loaders:** `load_xenium`, `load_visium`, `load_cosmx` (each returns a
>   `Shanuz` object with coordinates in an `FOV` slot); spatial-aware
>   `from_anndata` (rebuilds `images` from `obsm['spatial']`)
> - **Analysis:** `get_tissue_coordinates`, `spatial_knn`,
>   `nearest_neighbor_distance`, `local_neighborhood`, `build_niche_assay`
>   (`BuildNicheAssay`), `composition_test`, `add_module_score(search=)`
> - **Plots:** `image_dim_plot` (`ImageDimPlot`), `image_feature_plot`
>   (`ImageFeaturePlot`) — sub-cellular Xenium/CosMx centroids
>
> Remaining open items:

### `load_merscope`
- **R:** `LoadVizgen(data.dir)` (Vizgen MERSCOPE)
- **File format:** `cell_by_gene.csv`, `cell_metadata.csv`
- **Plan:** mirror `load_cosmx` — read the cell×gene matrix + metadata
  (`center_x` / `center_y`) into a `Shanuz` object with populated `images`.

### `FindSpatiallyVariableFeatures`
- **R:** `FindSpatiallyVariableFeatures(obj, method = "moransi")`
- **Plan:**
  1. **Moran's I** (primary): spatial autocorrelation statistic using the cell
     adjacency/distance weight matrix (weights from `spatial_knn`); compute with
     `esda.Moran` (`pysal` dep) or directly: `I = (N/W) * (z @ W @ z) / (z @ z)`
  2. **Markvariogram** (alternative): from the `Trendsceek` method
  3. Store results in `assay.var` (feature metadata), analogous to `FindVariableFeatures`
- **Python dep:** `libpysal` or pure NumPy implementation

### Tissue-image spatial plots (Visium H&E)
| Function | R equivalent | Notes |
|---|---|---|
| `spatial_dim_plot` | `SpatialDimPlot` | clusters/ident over the H&E tissue image |
| `spatial_feature_plot` | `SpatialFeaturePlot` | gene expression over the tissue image |

The sub-cellular centroid variants (`image_dim_plot` / `image_feature_plot`) are
done; these add the Visium tissue-PNG background (`matplotlib.imshow` + `scatter`).

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

### PyPI publication — ✅ delivered
- `pip install shanuz` works: published as [`shanuz` v0.1.1](https://pypi.org/project/shanuz/)
  (`build` + `twine` added to `[dev]` extras; published to TestPyPI then PyPI;
  verified with a clean-venv install + import + mini-pipeline smoke test)
- Still open: replace the hard-coded `__version__` string in `shanuz/__init__.py`
  with an `importlib.metadata.version("shanuz")` lookup, so the version is
  defined in exactly one place (`pyproject.toml`)

### GitHub Actions CI — ✅ delivered
- **File:** `.github/workflows/ci.yml`
- **Matrix:** Python 3.10, 3.11, 3.12 on ubuntu-latest, via `astral-sh/setup-uv`
- **Jobs:** `ruff check shanuz/` (advisory — pre-existing lint debt not yet
  cleared, so it doesn't gate the build) and `pytest tests/ -q`
- **Triggers:** push to `main`, all PRs
- Still open: a dedicated `build` job (`python -m build`, verifies the wheel
  itself builds cleanly) and coverage reporting (`--cov=shanuz --cov-report=xml`)

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

1. ~~**Harmony** (`v0.2.0`)~~ — ✅ delivered (`run_harmony` / `integrate_layers`)
2. ~~**WNN** (`v0.4.0`)~~ — ✅ delivered (`find_multi_modal_neighbors` + `run_umap(graph=)`); CBMC tutorial section still open
3. ~~**GitHub Actions CI** (`v0.10.0`)~~ — ✅ delivered
4. **`FindTransferAnchors` / `TransferData`** (`v0.3.0`) — enables atlas-based annotation (next-cycle candidate; needs CCA/RPCA first)
5. **`FindSpatiallyVariableFeatures` (Moran's I) + `SpatialFeaturePlot`** (`v0.7.0`) — the remaining spatial gaps (loaders + niche/neighbourhood analysis already delivered)
6. ~~**`AggregateExpression`**~~ ✅ + **DESeq2** (`v0.6.0`) — `aggregate_expression`
   and `find_conserved_markers` delivered; DESeq2/MAST/bimod unlock multi-sample DE
7. **`SketchData`** (`v0.8.0`) — enables million-cell datasets
