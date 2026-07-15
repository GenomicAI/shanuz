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

## v0.2.0 — Batch Correction & Integration — ✅ complete

> **Why first:** batch effects are unavoidable in real datasets; Harmony is the
> most widely-used, has a pip-installable Python package, and has a well-defined
> scope. CCA/RPCA follow naturally.
>
> **Status:** all three integration paths delivered — Harmony (`run_harmony`),
> CCA/RPCA anchors (`find_integration_anchors` / `integrate_data`), and the
> `integrate_layers` dispatcher (`method="harmony"|"cca"|"rpca"`).

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

### CCA / RPCA integration (`IntegrateData` v4 API) — ✅ delivered
- Implemented in `shanuz/anchors.py` as `find_integration_anchors(...)` +
  `integrate_data(...)`, with the `IntegrationAnchors` result container. Anchors
  are mutual nearest neighbours in a shared CCA (SVD of the cross-covariance
  `AᵀB`) or reciprocal-PCA space, scored by neighbourhood consistency and
  filtered against the raw expression space; `integrate_data` corrects each
  query onto the reference with a Gaussian-weighted sum of anchor correction
  vectors and returns a merged object carrying an active `"integrated"` assay.
  Verified (`tests/test_integration.py`) to cluster by cell type, not batch.
- **R:** `FindIntegrationAnchors(list, reduction="cca")` → `IntegrateData(anchors)`
- **Implemented:**
  1. `find_integration_anchors(objects, reduction, dims, k_anchor, k_filter, k_score, reference)` → `IntegrationAnchors`
  2. `integrate_data(anchors, new_assay, k_weight, sd_weight)` → corrected `"integrated"` assay on merged object
  3. CCA: `numpy.linalg.svd` on the cross-covariance; RPCA: reciprocal PCA
     projections, MNN via `sklearn.neighbors`
- **Note:** *reference-based* — anchors link each dataset to `reference=0`, and
  every other dataset is corrected onto it (one of Seurat's supported modes).
  A full guide-tree over all pairwise anchors is a later refinement.
- **Reuse:** the same `IntegrationAnchors` object is what v0.3.0's reference
  mapping (`FindTransferAnchors` / `TransferData`) is built to consume.

### `IntegrateLayers` (Seurat v5 API) — ✅ delivered (harmony + cca + rpca)
- Implemented as `integrate_layers(obj, method=..., group_by=...)` in
  `shanuz/integration.py`. `method="harmony"` corrects an existing reduction;
  `method="cca"` / `"rpca"` split the object by `group_by`, run the anchor
  pipeline above, and store the batch-corrected embedding as a new reduction.
- **R:** `IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca")`
- **Dep:** Harmony and CCA/RPCA functions above

---

## v0.3.0 — Reference Mapping & Label Transfer — ✅ complete

> **Why:** label transfer from a curated atlas to a query dataset is a standard
> first-annotation step; all machinery (KNN, PCA, anchors) already exists.
>
> **Status:** all three pieces delivered — `find_transfer_anchors`
> (pcaproject / cca) and `transfer_data` (classification + imputation) in
> `shanuz/transfer.py`, plus `project_umap` / `map_query` in `shanuz/mapping.py`
> (placing the query in the reference UMAP). Reuses the `anchors.py` machinery
> end-to-end and adds no dependency.

### `FindTransferAnchors` — ✅ delivered
- Implemented as `find_transfer_anchors(reference, query, reduction="pcaproject")`
  in `shanuz/transfer.py`, returning a `TransferAnchors` container (anchor pairs
  `cell1/cell2/score` + the query embedding used for weighting). Default
  **pcaproject** projects the query through the reference's PCA loadings (computed
  on the shared anchor features, so the reference need not already carry a `pca`
  reduction); **cca** builds a jointly-learned space. MNN → score → filter reuse
  the same helpers as `find_integration_anchors` (`tests/test_transfer.py`).
- **R:** `FindTransferAnchors(reference, query, dims = 1:30)`
- **Dep:** CCA/RPCA anchors from v0.2.0 (`shanuz/anchors.py`)

### `TransferData` — ✅ delivered
- Implemented as `transfer_data(anchors, refdata)` in `shanuz/transfer.py`. Builds
  a per-query-cell weight over the anchors (the same distance-weighted,
  score-scaled Gaussian kernel `integrate_data` uses), then either **classifies**
  (a metadata column or 1-D label array → `predicted.id` +
  `prediction.score.<class>` per class, rows summing to 1, + `prediction.score.max`)
  or **imputes** (a 2-D `features × reference-cells` matrix → predicted query
  expression). Verified to recover a query's true cell types across an injected
  batch block (`tests/test_transfer.py`).
- **R:** `TransferData(anchors, refdata = reference$celltype)`
- **Tests:** query cells get the correct transferred label at >85% accuracy;
  per-class scores form a distribution; imputation recovers the marker split.

### `MapQuery` + `ProjectUMAP` — ✅ delivered
- Implemented in `shanuz/mapping.py` as `project_umap(query, reference)` and the
  `map_query(anchors, refdata="celltype")` convenience. `project_umap` is a
  two-step map: project the query through the reference's PCA loadings into the
  reference's PC space, then run the reference's *fitted* UMAP model in
  transform-only mode (`umap-learn`'s `UMAP.transform`), so the query lands in the
  reference's existing embedding. It aligns the loadings to only the reference-PCA
  features the query actually carries scaled, so a query missing a few variable
  features still projects. Result stored as `query.reductions["ref.umap"]`.
- `map_query` composes the whole workflow from a `TransferAnchors`: `transfer_data`
  the labels onto `query.meta_data` (`predicted.id` / `prediction.score.*`), then
  `project_umap` the query into the reference UMAP — one call from anchors to an
  annotated, atlas-placed query. Verified (`tests/test_mapping.py`) that projected
  query cells land nearest the matching reference type's UMAP centroid (>85%)
  despite an injected batch block, and that `map_query` annotates + places in one
  call.
- **R:** `MapQuery(anchorset, query, reference, refdata = ...)` then `ProjectUMAP(...)`
- **Note:** the fitted model lives in `reference.reductions["umap"].misc["umap_model"]`,
  stored by `run_umap` when embedding from a reduction (not a graph).

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

### `glm_pca` (GLM-PCA) — ✅ delivered (Poisson + negative binomial)
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
- **Negative binomial:** `family="nb"` fits `Y ~ NB(μ, θ)`, `Var = μ + μ²/θ`.
  Poisson understates the overdispersion in most scRNA-seq, and a few noisy genes
  then dominate a Poisson fit; NB down-weights them. The Fisher scoring loop is the
  Poisson one with a single divisor `1 + μ/θ` on both the residual and the working
  weight — as `θ → ∞` it collapses back onto Poisson exactly. The shared dispersion
  `θ` is estimated by maximum likelihood between factor updates (`optimize_theta`,
  MASS `theta.ml`-style Newton seeded from a method-of-moments estimate) or pinned
  at a value you pass. Stored in `misc["theta"]` (`inf` for a Poisson fit). A moving
  `θ` re-scales the deviance, so the monotone-deviance guarantee holds only when `θ`
  is held fixed (`optimize_theta=False`).
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

> **Status:** ✅ **complete.** Leverage-score sketching (`sketch_data` /
> `project_data` + `leverage_score`, `shanuz/sketch.py`) draws an information-dense
> subset of a huge dataset and extends the sketch's analysis back to every cell;
> BPCells-style lazy on-disk matrices (`LazyMatrix`, `shanuz/lazy.py`) keep a matrix
> out-of-core and stream over it. Both add **no new dependency** — sketching reuses
> the `mapping.py` / `transfer.py` machinery, and `LazyMatrix` is built on NumPy's
> memory-mapping alone.

### BPCells-style lazy/on-disk matrices — ✅ delivered
- Implemented as `LazyMatrix` in `shanuz/lazy.py`, with `write_lazy_matrix(matrix,
  path)` to persist a matrix out-of-core and `open_lazy_matrix(path)` to map it
  back (`is_lazy(x)` to test). A matrix is stored as a directory of three
  **memory-mapped** `.npy` arrays — the `data` / `indices` / `indptr` triple of a
  compressed-sparse-**column** matrix, exactly scipy's `csc_matrix` layout — plus a
  JSON header. Opening maps the arrays without reading them; peak memory is the
  slices you touch, not the whole matrix (`tests/test_lazy.py`).
- **A faithful drop-in for a sparse layer.** Indexing (`m[:, cells]`,
  `m[np.ix_(genes, cells)]`, slices, boolean masks) reads only the touched columns
  off disk and returns an ordinary `scipy.sparse.csc_matrix`, so every hot path
  that already accepts a sparse layer keeps working. `Assay5` layers hold any
  array-like, so `assay.set_layer_data("counts", lazy)` just works, and the
  full-densify helpers (`as_dense` / `np.asarray`) route through `__array__` — the
  deliberate "materialise everything" escape hatch you avoid on the huge path.
- **Streaming reductions.** `sum` / `mean` over axis 0 (per-cell) or 1 (per-feature)
  run in a single pass over the on-disk arrays; `nnz_per_col()` (the `nFeature`
  count) comes straight from `indptr`; `col_blocks(block_size)` walks the matrix in
  cell-blocks — the primitive for processing a million cells at bounded RAM.
- **R:** `BPCells` package enables out-of-core analysis on millions of cells.
- **Design note (CSC, not CSR):** shanuz matrices are `features × cells` and the
  scale-out ops (sketching, cell subsetting, per-cell normalisation) select
  **columns**, which CSC makes cheap. Row-only (feature) subsetting still scans the
  touched columns; keeping both orientations on disk (as BPCells does) to make that
  cheap too, and a full `.toarray()` audit of the hot paths to gate densification
  behind shape checks, are the natural follow-ons.

### `SketchData` (leverage-score sub-sampling) — ✅ delivered
- Implemented as `sketch_data(obj, ncells=5000, method="LeverageScore")` in
  `shanuz/sketch.py`, with the leverage computation exposed on its own as
  `leverage_score(obj)`. Each cell is sampled **without replacement** with
  probability proportional to its statistical leverage, so the rare states a
  uniform sample would drop are kept (indeed over-represented, which the tests
  check directly). Returns a standalone subset `Shanuz` object whose active assay
  is renamed to `"sketch"`; the scores are written back onto the source object's
  `leverage.score` metadata (`tests/test_sketch.py`).
- **The leverage computation is the point, not a detail.** The exact leverage of
  cell *i* is `ℓ_i = ‖U_i‖²` (the row norm of the left singular vectors), an
  *O(n·d²)* SVD — the very cost sketching exists to avoid. So (as Seurat's
  `LeverageScore` does) a sparse **CountSketch** `S` embeds the *n* rows into a
  small `nsketch × d` matrix `B` with `BᵀB ≈ AᵀA` in a single pass over the
  non-zeros; whitening `A` with `B`'s right singular vectors gives an approximately
  orthonormal `Z` and `ℓ_i ≈ ‖Z_i‖²`. When `nsketch ≥ n` no sketch is taken and
  the scores are exact — validated against the textbook `‖U_i‖²` to `1e-6`.
- **One caveat, documented in the docstring and a test:** a single-hash CountSketch
  is a faithful subspace embedding only once it has ~`d²` rows, so `nsketch` should
  comfortably exceed the feature count (the default 5000 does, for the few-thousand
  variable features a real sketch runs on).
- **R:** `SketchData(obj, ncells = 5000, method = "LeverageScore")`
- **Divergence from R:** Seurat stores the sketch as an extra assay on the same
  object; here it is a separate object (as the plan called for), which fits
  shanuz's `subset` model and keeps the full/sketch data cleanly apart.

### `ProjectData` — ✅ delivered
- Implemented as `project_data(full, sketch, ...)` in `shanuz/sketch.py`, the
  inverse of `sketch_data`. Projects every full-dataset cell through the *sketch's*
  PCA loadings (stored as `full.reductions["pca.full"]`) — the same transform-only
  linear map `project_umap` uses — and, when the sketch carries a fitted UMAP,
  through that model too (`full.reductions["ref.umap"]`, via `project_umap`).
  Optional `refdata` (`{new_col: sketch_col}` or a column name) carries the
  sketch's labels to the full data via `find_transfer_anchors` + `transfer_data`.
  Verified that projected sketch cells reproduce their sketch-PCA coordinates and
  that transferred labels recover the full data's cell types at >85%
  (`tests/test_sketch.py`).
- **R:** `ProjectData(obj, reference = sketch, reduction = "pca")`
- **Reuse:** no new machinery — the PCA/UMAP projection is `mapping.py`'s and the
  label transfer is `transfer.py`'s; `project_data` is the composition.

---

## v0.9.0 — Specialized Assay Methods

### `HTODemux` / `MULTIseqDemux` (cell hashing)
- **R:** `HTODemux(obj)`, `MULTIseqDemux(obj)`
- **`hto_demux` — ✅ delivered** (`shanuz/hto.py`): k-means (`k = n_hashtags + 1`)
  on CLR-normalised HTO counts, then a per-hashtag **negative-binomial** fit
  (maximum likelihood) to the tag's *lowest-expressing* cluster — its background —
  thresholded at `positive_quantile` (0.99) to call positive cells. Cells positive
  for zero / one / many hashtags are classified singlet / doublet / negative.
  Writes the Seurat columns `HTO_maxID`, `HTO_secondID`, `HTO_margin`,
  `HTO_classification`, `HTO_classification.global` plus a convenient `hash.ID`
  (also set as the active identity); learned cutoffs are stashed in
  `obj.misc["hto_demux"]`. `normalize=False` reuses a prior
  `normalize_data(method="CLR")` layer. Built on the existing CLR helper +
  sklearn k-means — no new dependency.
- **Still open:** `multiseq_demux` (McGinnis et al. — quantile-sweep bar-code
  classification); it shares `hto_demux`'s normalize-then-threshold skeleton.
  Seurat's `"clara"` k-medoids clustering option is not ported (`kfunc="kmeans"`
  only).

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
| v0.8.0 | *(none — sketching reuses existing machinery; `LazyMatrix` is built on NumPy memory-mapping alone)* |
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
4. ~~**`FindTransferAnchors` / `TransferData` / `MapQuery` / `ProjectUMAP`**~~ ✅
   (`v0.3.0`) — `shanuz/transfer.py` (`find_transfer_anchors` pcaproject/cca +
   `transfer_data` classification/imputation) and `shanuz/mapping.py`
   (`project_umap` + `map_query`) deliver atlas-based annotation and place the
   query in the reference UMAP. Built on the v0.2.0 anchor machinery — **v0.3.0 is
   complete**.
5. ~~**`FindSpatiallyVariableFeatures`** + **`SpatialFeaturePlot`**~~ ✅ (`v0.7.0`) — all four loaders, niche/neighbourhood analysis, both spatially-variable-feature methods (Moran's I and markvariogram), the `VisiumV2` tissue-image data layer and the `spatial_*` H&E plots delivered; **v0.7.0 is complete**
6. ~~**`AggregateExpression` + DESeq2**~~ ✅ (`v0.6.0`) — `aggregate_expression`,
   `find_conserved_markers`, and pseudobulk DESeq2 (`test_use="deseq2"`) delivered;
   MAST (`test_use="mast"`) and bimod (`test_use="bimod"`) too — **v0.6.0 complete**
7. ~~**`SketchData`** + **BPCells-style lazy matrices** (`v0.8.0`)~~ — ✅ delivered.
   `sketch_data` / `project_data` + `leverage_score` (`shanuz/sketch.py`):
   leverage-weighted subsampling for million-cell datasets, and projection of the
   sketch's PCA/UMAP/labels back to the full data. `LazyMatrix` (`shanuz/lazy.py`):
   out-of-core, memory-mapped compressed-sparse-column storage with lazy column
   reads, streaming block reductions, and a clean `Assay5`-layer drop-in — no new
   dependency. **v0.8.0 is complete.**
8. ~~**`run_spca` + `glm_pca`** (Poisson + negative binomial)~~ ✅ (`v0.5.0`) —
   **v0.5.0 is complete**; GLM-PCA now fits both `family="poisson"` and
   `family="nb"` (dispersion estimated by ML), closing the last gap in it
