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

## v0.5.0 — Additional Dimensionality Reductions — ✅ complete

> t-SNE and ICA are single functions wrapping a scikit-learn call. sPCA and
> GLM-PCA are not — both are implemented directly against NumPy/SciPy, and
> neither added a dependency.

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

### `run_spca` (supervised PCA) — ✅ delivered
- Implemented in `shanuz/reduction.py` as `run_spca(obj, graph="wsnn")`
  (`tests/test_spca_glmpca.py`).
- **The plan previously written here described the wrong algorithm** — a gene-graph
  Laplacian smoothing a gene × gene matrix. Seurat's `RunSPCA` takes a **cell × cell**
  graph (the documented call is `RunSPCA(reference, assay = "SCT", graph = "wsnn")`)
  and eigendecomposes `XᵀGX`, which is features × features. Corrected here.
- **What it does:** ordinary PCA maximises `vᵀXᵀXv` and knows nothing about which
  cells you consider neighbours. sPCA swaps the identity for a graph you already
  trust and maximises `vᵀXᵀGXv` — the gene axes that best reproduce that graph.
  Pass `G = I` and PCA falls back out exactly, which is how the implementation is
  tested (loadings match `run_pca` to a cosine > 0.999).
- **Why it matters:** the output is a *linear map from genes to components*, so a
  query dataset can be pushed into a reference's graph-defined space with one
  matrix multiply. That is why Azimuth maps onto sPCA rather than PCA, and it is
  the reduction v0.3.0's reference mapping will want.
- **One departure from R:** Seurat runs `irlba` (an SVD) on `XᵀGX`, ranking
  components by `|λ|`; we take the largest eigenvalues themselves, since `vᵀXᵀGXv`
  is the quantity being maximised and a graph can push eigenvalues negative. With
  non-negative edge weights the leading eigenvalues are positive, so the two
  orderings differ only in the tail.
- **R:** `RunSPCA(obj, assay, graph)`

### `glm_pca` (GLM-PCA) — ✅ delivered (Poisson)
- Implemented in `shanuz/glmpca.py` as `glm_pca(obj, n_components=10)`, following
  Townes et al. (2019). Pure NumPy/SciPy — **no `glmpca-py` dependency**, in
  keeping with how MAST and bimod were done (`tests/test_spca_glmpca.py`).
- **What it does:** log-normalise-then-PCA assumes the transformed counts are
  Gaussian with constant variance. They are not, and the pseudocount needed to
  survive `log(0)` distorts exactly the low-expression genes where the zeros live.
  GLM-PCA drops the transform and fits a low-rank model on the count scale:
  `Y[g,c] ~ Poisson(μ)`, `log μ = a[g] + o[c] + Σ_l U[g,l]·V[c,l]`, with the log
  library size as a *fixed* offset `o` so sequencing depth is a known quantity
  rather than a factor to be rediscovered. Factors land in `cell_embeddings`,
  loadings in `feature_loadings`, so `find_neighbors(reduction="glmpca")` and
  `run_umap` work downstream unchanged.
- **Fitting:** Fisher scoring, alternating over intercept → loadings → factors.
  Each block gets a diagonal Newton step (score ÷ Fisher information); under a log
  link and Poisson noise both are one matrix product. Any step that fails to lower
  the deviance is rejected and retried at half the step size, so the deviance falls
  monotonically by construction; the trace is kept in `misc["deviance"]`.
- **Initialisation is load-bearing, not a detail.** `U = V = 0` is an *exact saddle*
  of the log-likelihood — each block's score is a product with the other, so both
  vanish there. Starting near zero (the obvious choice, and `glmpca`'s own default)
  leaves the fit inching away from the saddle, and any relative-improvement stopping
  rule then declares convergence on a model that has fitted nothing. It fails
  *convincingly*: one step is enough to orient the factors, so clusters separate
  cleanly in a plot while the deviance sits at its null value. So the factors are
  seeded from the SVD of the intercept-only model's residuals instead, as Townes
  recommends. `test_glmpca_actually_fits_rather_than_stalling` guards it.
- **Still open:** `family="nb"` (negative binomial) raises `NotImplementedError`.
  Poisson understates the overdispersion in most scRNA-seq; NB needs a dispersion
  parameter estimated alongside the factors.
- **Scale:** the fit is dense in genes × cells. Pass a few thousand variable
  features, as you would to `run_pca`.
- **R:** `RunGLMPCA(obj, L = 10)` (via SeuratWrappers + glmpca)

---

## v0.6.0 — Pseudobulk DE & Advanced Marker Methods — ✅ complete

### `AggregateExpression` (pseudobulk) — ✅ delivered
- Implemented as `aggregate_expression(...)` in `shanuz/aggregate.py`. Sums raw
  counts per group via a single sparse `counts @ indicator` matmul; `group_by`
  accepts one or more metadata columns (joined with `"_"`, as Seurat does) or
  `"ident"`. Returns a features×groups `pd.DataFrame` (a `dict` for multiple
  assays), or a `Shanuz` object with one "cell" per group when
  `return_object=True` (`tests/test_pseudobulk_conserved.py`).
- **R:** `AggregateExpression(obj, group.by = c("celltype","donor"))`
- Intended input for `DESeq2`-style testing (see below).

### DESeq2-style pseudobulk DE — ✅ delivered
- Implemented as the `test_use="deseq2"` branch of `find_markers` (`_deseq2_pseudobulk`
  in `shanuz/markers.py`). Sums counts to one pseudobulk profile per (group ×
  `sample_col`) — the `AggregateExpression` operation — then fits
  `pydeseq2.DeseqDataSet(design="~condition")` and contrasts group 1 vs group 2.
  Returns Seurat-shaped columns (`p_val`/`avg_log2FC`/`pct.1`/`pct.2`/`p_val_adj`);
  warns below 2 replicates per group. Enable with `pip install shanuz[deseq2]`
  (`tests/test_deseq2_pseudobulk.py`).
- **R:** `FindMarkers(obj, test.use = "DESeq2")` (via `DESeq2` R package)
- **Python dep:** `pydeseq2` (pip, optional `[deseq2]` extra)

### MAST — ✅ delivered
- Implemented as the `test_use="mast"` branch of `find_markers` (`_mast_pvalue` in
  `shanuz/markers.py`): a pure-Python two-part hurdle LRT — a logistic model of
  detection (`expr > 0`) plus a Gaussian model of magnitude among detected cells,
  each `~ group (+ latent)`. The combined statistic is the sum of the two
  components' LR statistics on the sum of their df (components with no signal
  drop out). No R dep — `statsmodels` is already present. Pass the cellular
  detection rate via `latent_vars` to match Seurat's CDR covariate
  (`tests/test_mast_de.py`).
- **R:** `FindMarkers(obj, test.use = "MAST")`

### `FindConservedMarkers` — ✅ delivered
- Implemented as `find_conserved_markers(...)` in `shanuz/markers.py`. Runs
  `find_markers` independently within each level of `grouping_var`, keeps genes
  that are markers in *every* level, and combines their per-level p-values with
  Fisher's method (`scipy.stats.combine_pvalues`). Output has per-level prefixed
  stats plus `max_pval` and `combined_p_val` (sorted by the latter); levels
  lacking a comparison group are skipped with a warning
  (`tests/test_pseudobulk_conserved.py`).
- **R:** `FindConservedMarkers(obj, ident.1, grouping.var)`

### `bimod` test (likelihood-ratio on bimodal model) — ✅ delivered
- Implemented as the `test_use="bimod"` branch of `find_markers` (`_bimod_pvalue`
  / `_bimod_likelihood` in `shanuz/markers.py`), a faithful port of Seurat's
  `DifferentialLRT`/`bimodLikData`: each group's expression is modelled as a
  point mass at zero (Bernoulli detection rate) plus a Gaussian on the detected
  values, and `2·(logLik₁ + logLik₂ − logLik_pooled)` is tested as χ²(df=3). Pure
  Python, no new dep (`tests/test_bimod_de.py`).
- **R:** `FindMarkers(obj, test.use = "bimod")`

---

## v0.7.0 — Spatial Transcriptomics

> **Delivered.** The data structures (`FOV`, `Centroids`, `Segmentation`,
> `Molecules`) plus these loaders and analysis functions are done and validated
> end-to-end against R Seurat in
> [Tutorial 5](tutorials/xenium_spatial_tutorial.md) (deterministic anchors match
> to 8 significant figures):
>
> - **Loaders:** `load_xenium`, `load_visium`, `load_cosmx` (each returns a
>   `Shanuz` object with coordinates in an `FOV` slot); spatial-aware
>   `from_anndata` (rebuilds `images` from `obsm['spatial']`)
> - **Analysis:** `get_tissue_coordinates`, `spatial_knn`,
>   `nearest_neighbor_distance`, `local_neighborhood`, `build_niche_assay`
>   (`BuildNicheAssay`), `composition_test`, `add_module_score(search=)`,
>   `find_spatially_variable_features` (`FindSpatiallyVariableFeatures`) with both
>   the **moransi** and **markvariogram** methods
> - **Plots:** `image_dim_plot` (`ImageDimPlot`), `image_feature_plot`
>   (`ImageFeaturePlot`) — sub-cellular Xenium/CosMx centroids; `spatial_dim_plot`
>   (`SpatialDimPlot`), `spatial_feature_plot` (`SpatialFeaturePlot`) — Visium
>   spots over the H&E tissue image
>
> This milestone is complete.

### `load_merscope` — ✅ delivered
- Implemented in `shanuz/spatial/loaders.py` as `load_merscope(...)`, mirroring
  `load_cosmx`: reads `cell_by_gene.csv` + `cell_metadata.csv` (`center_x` /
  `center_y`) into a `Shanuz` object with populated `images` (one per `fov`).
  Drops `Blank-*` control barcodes by default (as `LoadVizgen` does; override
  with `keep_controls=True`) and tolerates both the named and unnamed cell-id
  column layouts Vizgen emits. Verified to flow through the whole spatial stack
  (`spatial_knn` → `nearest_neighbor_distance` → `build_niche_assay`)
  (`tests/test_merscope_loader.py`).
- **R:** `LoadVizgen(data.dir)` (Vizgen MERSCOPE)
- **File format:** `cell_by_gene.csv`, `cell_metadata.csv`

### `FindSpatiallyVariableFeatures` — ✅ delivered (both methods)
- Implemented as `find_spatially_variable_features(...)` in
  `shanuz/spatial/variable_features.py`, dispatching on `method=`. Pure
  NumPy/SciPy — no `libpysal` and no `spatstat` equivalent needed.

**`method="moransi"`** (default) builds a row-standardised sparse KNN weight
matrix from `spatial_knn`, then computes `I = (N/S0)·(zᵀWz)/(zᵀz)` vectorised
across all genes in one sparse matmul. Significance uses the closed-form
`E[I] = −1/(N−1)` and normality-assumption variance (computed once, since it
depends only on W) → z-score → two-sided p, plus BH adjustment. Writes
`moransi` / `moransi_pval` / `moransi_padj` / `moransi_rank` into the assay's
feature metadata (as `find_variable_features` does). Validated against a
brute-force double sum.

**`method="markvariogram"`** computes the normalised mark variogram
`γ(r) = E[½·(m_i − m_j)² | d_ij ≈ r] / Var(m)` — the expression difference
between cells about `r` apart, relative to the gene's own variance. `γ ≈ 1`
means two cells `r` apart differ as much as two picked at random (no structure);
`γ < 1` means they still resemble each other. Writes `markvariogram` /
`markvariogram_rank`; rank 1 = lowest γ. No p-value — the variogram has no
closed-form null, and R does not offer one either. Also validated against a
brute-force loop over every cell pair (`tests/test_markvariogram.py`).

- **Two deliberate departures from R,** both documented in the docstring:
  - **`r_metric` is in nearest-neighbour spacings, not raw coordinate units.** R
    passes `r.metric` straight through to `spatstat`, so the same script answers
    differently on a slide in pixels and in microns, and the default of 5 is only
    meaningful if you know your coordinate scale. Here `r_metric=5` means "five
    cells apart" on any slide.
  - **γ is a kernel-weighted (Nadaraya-Watson) ratio estimator**, not
    `spatstat`'s translation-corrected one, so absolute γ values are close to but
    not identical with R's. The gene *ranking* — what the function is for —
    carries over.
- **Performance:** the pairwise differences are never materialised (that array
  would be genes × pairs). Because the kernel matrix K is symmetric with a zero
  diagonal, `Σ_{i<j} K_ij·(m_i − m_j)² = Σ_i s_i·m_i² − mᵀKm` with `s = K·1`,
  which is two sparse products regardless of how many pairs land in the band. A
  `cKDTree` range query means only pairs near `r` are ever built.
- **Caveat documented in the docstring:** when a few strongly spatial genes
  dominate library size, log-normalisation leaks their structure into flat genes
  and inflates their score — a property of compositional normalisation, not of
  either statistic.
- **R:** `FindSpatiallyVariableFeatures(obj, method = "moransi" | "markvariogram")`

### Visium tissue image (`VisiumV2`) — ✅ delivered (data layer)
- `shanuz/spatial/visium.py` adds `VisiumV2` (an `FOV` subclass mirroring Seurat
  v5's) carrying the H&E image, the `ScaleFactors` from `scalefactors_json.json`,
  and which resolution is stored. `load_visium(..., image=True)` (the default) now
  reads `spatial/tissue_{hires,lowres}_image.png` + `spatial/scalefactors_json.json`
  and returns `VisiumV2` images; bundles with neither still load as a plain `FOV`,
  so this is backwards compatible. Also adds `filter_by_tissue=` (drops
  `in_tissue == 0` spots from the matrix as well as the coordinates), and registers
  the previously unimplemented `get_image` generic on `SpatialImage`
  (`tests/test_visium_image.py`).
- **Coordinate convention:** spot coordinates stay in **full-resolution pixels**
  (the space `tissue_positions.csv` uses), so `spatial_knn` /
  `nearest_neighbor_distance` / Moran's I keep seeing real, image-independent
  distances. `VisiumV2.scale_coordinates()` and `.spot_radius()` convert to the
  stored image's pixel space on demand — one multiply, at draw time.
- **Deps:** none added. The PNG is read with matplotlib (or Pillow) lazily; with
  neither installed the loader warns and skips the image rather than failing.
- **R:** `Load10X_Spatial(data.dir)`, `GetImage(obj[["slice1"]])`, `ScaleFactors()`

### Tissue-image spatial plots (Visium H&E) — ✅ delivered
| Function | R equivalent | Notes |
|---|---|---|
| `spatial_dim_plot` | `SpatialDimPlot` | clusters/ident over the H&E tissue image |
| `spatial_feature_plot` | `SpatialFeaturePlot` | gene expression over the tissue image |

- Both in `shanuz/plotting.py`. Each panel `imshow`s the photo held by the
  `VisiumV2` image, then overlays the spots on top of it
  (`tests/test_spatial_plots.py`).
- **Spots are drawn at their true diameter**, not as scatter points: they are an
  `EllipseCollection` in `units="xy"`, so a spot's size is expressed in *data*
  units and stays registered against the tissue when the axes are zoomed or the
  resolution changes. A `scatter(s=…)` size is in points² and would drift.
  `pt_size_factor=` (default 1.6, as in Seurat) scales relative to the real spot.
- **The image anchors the coordinate space.** Coordinates come from
  `VisiumV2.scale_coordinates()` (fullres → image pixels) and the diameter from
  `.spot_radius()`; `imshow` then supplies the frame. `crop=` (default `True`)
  zooms to the spots rather than the whole slide, `image_alpha=` fades the
  tissue, and `resolution=` overrides which PNG is drawn.
- **Two independent fallbacks, because the bundle has two independent holes.** No
  PNG (a plain `FOV`, or `load_visium(image=False)`) → a bare scatter of the same
  spots, y-axis still pointing down so it looks the same. No
  `scalefactors_json.json` → the photo still draws, but with no
  `spot_diameter_fullres` there is nothing to size the spots *to*, so they degrade
  to fixed-size scatter points.
- **R:** `SpatialDimPlot(obj)`, `SpatialFeaturePlot(obj, features = "Gad1")`

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
| v0.5.0 | *(none — sPCA and GLM-PCA are pure NumPy/SciPy; `glmpca-py` proved unnecessary)* |
| v0.6.0 | `pydeseq2` (optional) |
| v0.7.0 | *(none — Moran's I is pure NumPy/SciPy; the Visium PNG uses matplotlib, already in `[analysis]`)* |
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
4. **`FindTransferAnchors` / `TransferData`** (`v0.3.0`) — enables atlas-based annotation (next-cycle candidate; needs CCA/RPCA first). `run_spca` is now in place, which is the reduction Azimuth maps onto.
5. ~~**`FindSpatiallyVariableFeatures`** + **`SpatialFeaturePlot`**~~ ✅ (`v0.7.0`) — all four loaders, niche/neighbourhood analysis, both spatially-variable-feature methods (Moran's I and markvariogram), the `VisiumV2` tissue-image data layer and the `spatial_*` H&E plots delivered; **v0.7.0 is complete**
6. ~~**`AggregateExpression` + DESeq2**~~ ✅ (`v0.6.0`) — `aggregate_expression`,
   `find_conserved_markers`, and pseudobulk DESeq2 (`test_use="deseq2"`) delivered;
   MAST (`test_use="mast"`) and bimod (`test_use="bimod"`) too — **v0.6.0 complete**
7. **`SketchData`** (`v0.8.0`) — enables million-cell datasets
8. ~~**`run_spca` + `glm_pca`**~~ ✅ (`v0.5.0`) — **v0.5.0 is complete**; the only
   gap left in it is GLM-PCA's negative-binomial family
