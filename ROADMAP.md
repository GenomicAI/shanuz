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
- **RPCA (found in T6, `integration_vignette.md`):** the `reduction="rpca"` path had
  two real defects the synthetic tests missed, **both now fixed** — a crash on
  unequal batch sizes (fixed in #41) and a 4× under-integration (batch-mix 0.222 →
  **0.867** vs Seurat's 0.914; fixed by per-object scaling + Seurat's
  reciprocal-embedding SD/L2 normalization + disabling the RPCA anchor filter, all
  with regression tests). `reduction="cca"` and `run_harmony` match Seurat to three
  decimals; RPCA now integrates well, ~5% off Seurat on the residual (exact-vs-annoy
  NN, sklearn-vs-irlba PCA).
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

### `FindMultiModalNeighbors` — ✅ delivered (full port)
- Implemented as `find_multi_modal_neighbors(...)` in `shanuz/multimodal.py`.
  Both WNN stages are ported from the R/C++ source; stores `wknn`/`wsnn` graphs
  and `<assay>.weight` columns. Verified on synthetic complementary-modality
  data to recover structure RNA alone cannot (`tests/test_multimodal_wnn.py`).
- **R:** `FindMultiModalNeighbors(obj, reduction.list = list("pca","apca"), dims.list = list(1:30, 1:18))`
- **Stage 1 — `FindModalityWeights`** (`_modality_weights`): L2-normalise each
  embedding; impute each cell from its own modality's neighbours and from the
  other's; `d = ||x - x_hat|| - d(nearest neighbour)`, ReLU'd; per-cell kernel
  bandwidth from the SNN graph; `exp(-d / sigma)`; score
  `within / (cross + 1e-4)` **clipped to [0, 200]**; softmax across modalities.
- **Stage 2 — `MultiModalNN`** (`_multi_modal_nn`): each modality nominates
  `knn_range = 200` candidates, they are **unioned** per cell, and each is
  scored by `sum_r exp(-d_r / sigma_r) * weight_r`. The top `k_nn` become the
  cell's joint neighbours; `wknn`/`wsnn` are built from *that* ranking.
- **Tests:** WNN clusters CBMC data closer to protein-defined ground truth than RNA alone

#### Superseded: the earlier weighting approximation
The first implementation approximated stage 1 with a linear distance ratio,
`theta_m = d_cross / (d_same + d_cross)`, normalised across modalities — and
skipped stage 2 entirely, blending per-modality SNN graphs as
`SNN_wnn = w * SNN_RNA + (1-w) * SNN_ADT` instead. It was monotone in the right
quantity, so the *direction* was right, but the linear ratio has no dynamic
range: on CBMC every cell landed in 0.46–0.53, versus R's 0.21–0.65. A weight
pinned near 0.5 cannot say "this cell is decided by protein", which is the one
thing WNN exists to say. The exponential kernel plus the clipped softmax is what
supplies the range, and the joint neighbour search is what turns it into a
graph. Caught by putting the shanuz and Seurat weight violins side by side in
Tutorial 3 — the tutorial figure was the test.

Four details in the R source that the natural reading gets wrong, all load-bearing:
- `FindMultiModalNeighbors` never passes `prune.SNN` down to `FindModalityWeights`,
  so the **bandwidth's SNN graph uses `prune = 0`**, not `1/15`.
- `SNN_SmallestNonzero_Dist` (`src/snn.cpp`) averages distances to the
  **least**-similar SNN partners, and on ties at the k-th weight keeps the
  **k largest** distances.
- `PredictAssay` silently drops the self column, so imputation averages
  `k_nn - 1` neighbours.
- R requests `k.nn` neighbours **including** self, not `k.nn + 1`.

- **Known scaling limit:** `_multi_modal_nn` loops per cell over the unioned
  candidate pool. Fine at CBMC's ~8.6k cells (Tutorial 3 runs end to end in
  ~58s), but it is the bottleneck if WNN is pointed at a much larger dataset.
  Vectorising it is the obvious follow-on if that comes up.

#### Fixed: CLR's `margin` flag was inverted against Seurat
Not a WNN bug, but found through WNN and worth recording next to it. Seurat's
`CustomNormalize` applies the CLR kernel via
`apply(data, MARGIN = margin, clr_function)`, and R's `apply` treats `MARGIN = 1`
as rows and `MARGIN = 2` as columns. With counts stored features × cells that
makes Seurat's `margin=1` **per-feature across cells** (its default) and
`margin=2` **per-cell across features** (what ADT panels want).
`_clr_normalize` had the two swapped, so `normalize_data(..., margin=2)` computed
what Seurat's `margin=1` computes. The per-vector kernel was always exact,
including Seurat's quirk of summing `log1p` over non-zero entries while dividing
by the full length — only the axis was wrong.

Consequences, all now resolved: the ADT matrix feeding `apca` was the wrong
transform, so every `ADT.weight` was off (erythroid 0.62 vs Seurat's 0.40,
platelet 0.53 vs 0.30, and the whole table biased high); and the annotation
thresholds in `cbmc_citeseq_verify.R` had been retuned around the discrepancy and
documented as a CLR "scale" difference between the languages, which it was not.
Both scripts now share one set of thresholds.

`hto_demux` and `multiseq_demux` were **accidentally correct** — they passed
`margin=2` and, under the inverted code, got per-hashtag-across-cells, which is
numerically what Seurat's hashing vignette does at its default `margin=1`. Their
defaults moved 2 → 1 in the same change so their behaviour is unchanged; that
coupling is the trap in this fix.

It survived because `test_clr_matches_seurat_formula` verified the kernel but
derived the axis mapping in Python, restating the same inverted assumption it
was meant to check — right formula, assumed axis, exactly the failure mode of the
superseded weighting above. `test_clr_margin_matches_r_ground_truth` replaces
that with fixed output captured from a real Seurat 5.5.1 run; the old
implementation fails it on both margins.

**Breaking:** any caller passing `margin` explicitly to `normalize_data`,
`hto_demux` or `multiseq_demux` gets different output than before. Callers who
relied on the defaults are unaffected.

### WNN UMAP + clustering — ✅ delivered
- `find_clusters(graph_name="wsnn")` already routed correctly; `run_umap` now
  accepts a `graph=` kwarg that embeds a precomputed graph via UMAP's
  `simplicial_set_embedding` (`tests/test_reductions_extra.py`).
- CBMC CITE-seq tutorial (Tutorial 3) Step 8 covers WNN end to end, with R-side
  figures for the joint embedding and the modality weights.
- **Agreement with R:** both sides find 16 RNA clusters and 21 WNN clusters at
  `resolution = 0.6` and resolve the same nine lineages with the same
  multiplicities. Eight of nine per-cell-type `ADT.weight` means match Seurat to
  0.02 or better; progenitor is the exception at 0.06, on a 146-cell population
  the two sides do not cut identically. shanuz's neighbour search is exact where
  R's is approximate (annoy) and the Louvain implementations differ, so cluster
  boundaries still move slightly — small populations feel it most.

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
- **`hto_demux` — ✅ delivered** (`shanuz/hto.py`): clustering into
  `k = n_hashtags + 1` groups (`kfunc="clara"` by default, as in Seurat, or
  `"kmeans"` — see below) on CLR-normalised HTO counts, then a per-hashtag **negative-binomial** fit
  (maximum likelihood) to the tag's *lowest-expressing* cluster — its background —
  thresholded at `positive_quantile` (0.99) to call positive cells. Cells positive
  for zero / one / many hashtags are classified singlet / doublet / negative.
  Writes the Seurat columns `HTO_maxID`, `HTO_secondID`, `HTO_margin`,
  `HTO_classification`, `HTO_classification.global` plus a convenient `hash.ID`
  (also set as the active identity); learned cutoffs are stashed in
  `obj.misc["hto_demux"]`. `normalize=False` reuses a prior
  `normalize_data(method="CLR")` layer. Built on the existing CLR helper, with
  `clara` in-tree and sklearn k-means optional — no new dependency.
- **`multiseq_demux` — ✅ delivered** (`shanuz/multiseq.py`): the MULTI-seq
  chemistry (McGinnis et al.). Reuses `hto_demux`'s CLR extraction, then thresholds
  each barcode straight off its distribution shape rather than a background fit — a
  Gaussian KDE over a 100-point grid exposes the background and positive modes as
  the two tallest local maxima and the cutoff sits a fraction `quantile` (0.7)
  between them. Cells clearing zero / one / many barcodes are singlet / doublet /
  negative. `autothresh=True` runs McGinnis's iterative `deMULTIplex` sweep — pick
  the `q` maximising the singlet rate, peel off negatives, re-threshold the
  remainder up to `maxiter` rounds. Writes `MULTI_ID` (also the active identity)
  and `MULTI_classification`; thresholds stashed in `obj.misc["multiseq_demux"]`.
  KDE via `scipy.stats.gaussian_kde` — no new dependency.
- **`kfunc="clara"` — ✅ delivered** (`shanuz/_clara.py`): Seurat's default
  clustering for `HTODemux`, a port of the `clara` C routine in R's **cluster**
  package (2.1.8.2). CLARA is k-medoids for data too big for PAM: it draws
  `nsamples` (100, Seurat's default) sub-samples of `min(n, 40 + 2k)` cells, runs
  PAM on each, assigns every cell to the resulting medoids, and keeps the
  sub-sample with the lowest total dissimilarity. Ported rather than taken from a
  library because the details that decide the answer are all non-standard: clara
  draws from its **own** 16-bit LCG (`rngR = FALSE`), so `set.seed` cannot reach it
  and `clara` correctly takes no `seed` argument; at `pamLike = FALSE` the swap
  rule is the pre-2011 one the C source itself calls "a bit illogical", *not*
  `pam()`'s; ties break **last**-wins in BUILD but **first**-wins in SWAP; and
  cluster numbering follows first appearance. **`kfunc="clara"` is the default**,
  matching Seurat; `"kmeans"` remains available. `hto_demux` first shipped
  defaulting to `"kmeans"`, so this is a behaviour change for callers who never
  passed `kfunc` — one for the CHANGELOG when v0.10.0 adds it. The switch moves
  ~1% of calls on synthetic panels (rising with tag count — ~3.5% at 12 tags,
  where clara is also the *more* accurate of the two) because only each tag's
  least-expressing cluster feeds the background fit. Both scale linearly in cells;
  clara costs a roughly constant 4× (~1.3 s vs ~0.3 s at 100k cells). No new
  dependency; `clara` needs no sklearn.

#### Note: `clara`'s answer depends on the CPU it was compiled for
clara accepts a swap on *any* improvement below zero — R genuinely takes swaps
worth `-2.2e-16` — so a one-ulp difference in a single distance flips a swap, then
the winning sub-sample, then the whole clustering. That ulp is a compiler's
choice: `clara.c` built for **arm64** contracts `clk += d*d` into a fused
multiply-add (one rounding), while the same source built for baseline **x86_64**,
whose ISA has no FMA, rounds twice. Compiling the real `clara.c` for both targets
and running them against each other confirms it — same source, same input,
materially different clusterings on ~2% of random inputs and ~7% of realistic
hashtag panels. **R's clara is therefore not reproducible across architectures**,
and "bit-exact to R" is not a well-defined target.

shanuz follows plain IEEE double arithmetic, which is what numpy gives on every
platform and what `clara.c` gives on x86_64. Against that reference the port is
exact — 200/200 pathological and 40/40 realistic HTO cases — and the only inputs
where it differs from an arm64 R are exactly the ones where `clara.c` disagrees
with itself across architectures. Emulating FMA to chase the arm64 answer would
simply break fidelity the other way; there is no choice that satisfies both.
Three details in `_clara.py` exist solely to hold that reference — `np.cumsum`
rather than `np.sum` in `_selec`, the per-`j` accumulation of `dz` in `_bswap2`,
and `dz`'s h-major layout — and provably have **no** observable effect on any
architecture-stable input (no fixture among 1490 can distinguish them). They are
not redundant: without them the port drifts from the IEEE reference in the chaotic
regime. Don't "simplify" them.

### Mixscape (pooled CRISPR screen analysis) — ✅ delivered
- **R:** `CalcPerturbSig(obj, ...)` → `RunMixscape(obj, labels, nt.class.name)` →
  `MixscapeLDA(obj, labels, nt.label)`; plots `PlotPerturbScore(obj,
  target.gene.ident)` and `MixscapeHeatmap(obj, ident.1, ident.2)`
- **`calc_perturb_sig` — ✅ delivered** (`shanuz/mixscape.py`): Seurat's
  `CalcPerturbSig`. For every cell, the mean expression of its `num_neighbors`
  (20) nearest **non-targeting (NT)** control cells — in the first `ndims` of a
  reduction (default `pca`), optionally within each `split_by` batch — is
  subtracted from its own expression. The residual local perturbation signature
  (shared technical variation cancelled) is stored as a new assay (default
  `"PRTB"`). Neighbours via `sklearn.neighbors.NearestNeighbors` — no new dep.
- **`run_mixscape` — ✅ delivered** (`shanuz/mixscape.py`): Seurat's `RunMixscape`.
  Per target gene: (1) DE vs NT on `de_assay` (reuses `find_markers`) picks the
  perturbation-response genes; a gene with fewer than `min_de_genes` (5) is all
  NP. (2) On the signature over those genes, a *perturbation vector*
  (mean KO − mean NT) projects every cell to a single perturbation score, and an
  **iterative 2-component `sklearn.mixture.GaussianMixture`** splits the knockout
  (KO) mode from the non-perturbed (NP) mode — the vector is recomputed from each
  new KO set and re-fit until it stabilises (up to `iter_num` rounds). Writes
  `mixscape_class` (`"<gene> KO"` / `"<gene> NP"` / `"NT"`, also the active
  identity), `mixscape_class.global`, and `mixscape_class_p_<type>`; per-gene
  bookkeeping stashed in `obj.misc["mixscape"]`.
- **Two documented departures from R:** the mixture is fit on the pooled NT +
  gene-cell scores (so the NP mode is anchored by the controls), and the
  signature is read from the assay's `data` layer directly — Seurat's pre-`ScaleData`
  centring is a provable no-op for the KO/NP calls (it leaves the perturbation
  vector unchanged and only globally shifts the scores the re-fit mixture absorbs).
- **`mixscape_lda` — ✅ delivered** (`shanuz/mixscape.py`): Seurat's `MixscapeLDA`
  (its `PrepLDA` + `RunLDA`). The complementary question to `run_mixscape` — not
  *which cells are perturbed* but *how do the guide populations differ from each
  other and from control* — answered with one supervised map. Per guide: DE vs NT
  picks its response genes (a guide needs `npcs + 1`, default 10, to contribute);
  on the signature assay over those genes a PCA is fit on that guide's cells **plus**
  the NT cells, and **every** cell is projected onto that guide's `npcs`-dim
  subspace. The per-guide blocks are concatenated and a single
  `sklearn.discriminant_analysis.LinearDiscriminantAnalysis` is fit with the guide
  label as the class. Stores an `lda` reduction (key `LDA_`, `n_classes − 1` dims,
  loadings + assignments/posterior/`genes_used` in `misc`) and the
  `lda_assignments` / `LDAP_<class>` metadata columns.
- **Prerequisite is only `calc_perturb_sig`:** faithful to R, the LDA groups cells
  by their **raw guide label**, not `mixscape_class` — the KO/NP calls are not
  consumed and no cells are dropped, so NP escapers stay in their guide's group
  (and, being control-like, generally land on NT).
- **Two more documented departures from R:** Seurat's `ScaleData` → `RunPCA` →
  `ProjectCellEmbeddings` chain is composed directly (centre/scale each response
  gene against the guide-plus-NT reference, project through the reference PCA
  loadings — the same map, `scale_max=10` clip included); and MASS's leave-one-out
  CV posterior (`lda(..., CV = TRUE)`, stashed in `misc` by R and read by nothing)
  is not computed, only the resubstitution assignment and posterior.
- **`plot_perturb_score` — ✅ delivered** (`shanuz/plotting.py`): Seurat's
  `PlotPerturbScore`. The diagnostic for *why* mixscape split a guide the way it
  did — the perturbation score is the one axis the mixture is fit on, and the plot
  overlays the NT control density against the guide's own. A guide with a real
  effect is bimodal: one lobe on the NT curve (the escapers), one shifted away
  (the knockouts). `before_mixscape=False` (default) colours by `mixscape_class`
  (NT / `"<gene> NP"` / `"<gene> KO"`, i.e. where mixscape drew the line);
  `before_mixscape=True` colours by the raw guide label, the pre-mixscape view.
  Cells are also drawn as a jittered strip (controls above the axis, target gene
  below), and `split_by` facets a multi-cell-type screen. Densities via
  `scipy.stats.gaussian_kde` (R uses `geom_density`) — no new dep.
- **The perturbation score is now persisted** (`shanuz/mixscape.py`):
  `run_mixscape` previously computed the score and discarded it. R keeps it in the
  `Tool(object, "RunMixscape")` slot precisely so `PlotPerturbScore` can read it
  back, so it is now stored per gene under `obj.misc["mixscape"][assay]["genes"]`
  as a `scores` frame (a `pvec` column plus the guide-label column, indexed by
  cell — R's `gv` data.frame). Two details are faithful to R's source rather than
  to the obvious reading: the score is kept from the **first** iteration
  (R's `if (n.iter == 0)`), before any KO/NP split feeds back into the
  perturbation vector — a later round would shift the axis; and it is normalised
  by `vec · vec` (R's `ProjectVec` divides by `v2 %*% v2`), which the classifier
  never needed (a positive constant cannot move the split) but a plotted axis does.
  Genes that short-circuit to NP without a mixture fit store `scores = None`.
- **`mixscape_heatmap` — ✅ delivered** (`shanuz/plotting.py`): Seurat's
  `MixscapeHeatmap`. The genes underneath the score: DE genes between two
  `mixscape_class` levels (e.g. `"NT"` vs `"IFNGR2 KO"`), with every cell ordered
  left-to-right by its knockout posterior, so the expression block is seen turning
  on in step with the probability. `balanced=True` takes up to `max_genes` from
  each fold-change direction; `max_cells_group` downsamples per class;
  `order_by_prob=False` shuffles instead (as R does). Delegates to `do_heatmap`
  exactly as R delegates to `DoHeatmap`.
- **`do_heatmap` gained R's `cells` argument** to support that ordering (restrict
  to those cells *and* plot them in that exact order). Two latent bugs surfaced
  and were fixed with it: a `scale.data` layer scaled over a feature subset
  (`scale_data(obj, features = [...])`) has fewer rows than the assay has
  features, so indexing rows by the assay's full feature list read the wrong
  genes — layer-aware resolution now goes through `_resolve_layer`; and the
  group-label pass assumed each group was one contiguous block, which an explicit
  `cells` order interleaves, so it now walks runs.
- **Departure from R:** the R colour names in the signatures (`col = "orange2"`,
  and the fixed `grey49` / `grey79` for controls and escapers) are kept so the
  call reads like the R one, and translated to hex for matplotlib.
- **Still open:** `DEenrichRPlot` (mixscape's third plot) is not ported — it calls
  the enrichR web service for pathway enrichment, which would put a live network
  dependency in the test path.

---

## v0.10.0 — Package Infrastructure

### PyPI publication — ✅ delivered
- `pip install shanuz` works: published as [`shanuz`](https://pypi.org/project/shanuz/),
  currently 0.2.0 (`build` + `twine` added to `[dev]` extras; published to
  TestPyPI then PyPI; verified with a clean-venv install + import + mini-pipeline
  smoke test)
- ~~Still open: replace the hard-coded `__version__` string~~ — ✅ delivered:
  `shanuz/__init__.py` reads `importlib.metadata.version("shanuz")`, with a
  PEP 440-valid `0.0.0+unknown` fallback for a source tree that has no installed
  distribution. `pyproject.toml` is now the only place the version is written.
  **Tradeoff worth knowing:** that lookup resolves through the *install-time*
  `dist-info` snapshot, so an editable install keeps reporting the old version
  after a bump until it is reinstalled — a way to be wrong that the hard-coded
  string did not have. `tests/test_packaging.py::test_version_matches_pyproject`
  exists to make it loud rather than silent.
- **Still open — decide what the `>=` dependency floors mean.** `pyproject`
  declares only lower bounds (`pandas>=2.0`, `numpy>=1.24`, `anndata>=0.10`,
  `umap-learn>=0.5`, `scikit-learn>=1.3`, …), so a fresh install resolves to
  whatever shipped this week and no two installs need agree. This has now caused
  four distinct failures, none of them theoretical:
  1. tutorial figures drifting against the committed R panels;
  2. `mypy` aborting on `anndata` 0.13.2's PEP 695 syntax and silently checking
     nothing (fixed by dropping `python_version`);
  3. the mypy baseline spanning 81–85 across the CI matrix, because each leg
     resolves different versions;
  4. **pandas 3 crashing `pbmc3k_tutorial.py` on every fresh install** while a
     developer venv pinned at pandas 2 stayed green (#36).
  Note the asymmetry that makes this expensive: the *maintainer* never sees it.
  Existing venvs hold old versions; only new users get the break. Pinning upper
  bounds is not obviously right either — it locks users out of new releases and
  ages badly. The decision to make is which of {upper bounds, a tested lockfile
  for CI, a scheduled "resolve latest" CI leg} this project wants; the status quo
  is "find out from a user".

- **Still open — cut a release.** The tags stop at 0.2.0 (2026-07-05) while
  milestones v0.3.0–v0.9.0 have all landed on `main`, so `pip install shanuz`
  currently ships almost none of the README's feature list: of 22 advertised
  entry points sampled, 19 are absent from the published wheel (reference
  mapping, sketching, `LazyMatrix`, hashing, Mixscape, `run_spca`/`glm_pca`,
  pseudobulk DE, MERSCOPE/Visium). The README now says so out loud, but the real
  fix is a release. This is the highest-value remaining infra item.

### GitHub Actions CI — ✅ delivered
- **File:** `.github/workflows/ci.yml`
- **Matrix:** Python 3.12, 3.13 on ubuntu-latest, via `astral-sh/setup-uv`
  (was 3.10–3.12; moved to track [SPEC 0](https://scientific-python.org/specs/spec-0000/),
  which had already retired 3.10 in Oct 2024 and 3.11 in Oct 2025).
- **Open — add the 3.14 leg when `harmonypy` ships a cp314 wheel.** 3.14 was in
  this change and came out again: it is the *only* package in the set without
  one (manylinux cp39–cp313 only), and both ways around it are worse than
  waiting. Building from source needs BLAS plus a CMake-fetched armadillo that
  `ubuntu-latest` lacks — that is the failure we hit. Forcing a wheels-only
  resolve instead backtracks to `harmonypy` 0.2.0, which depends on torch and
  pulls triton and 24 `nvidia-*` packages. Everything else resolves clean: 95
  packages, wheels only, on `x86_64-manylinux_2_28`/3.14. Recheck with
  `uv pip compile pyproject.toml --extra all --python-platform
  x86_64-manylinux_2_28 --python-version 3.14 --only-binary :all:`.
- **Verify wheels on the *target* platform, not the dev machine.** The 3.14 leg
  was signed off locally on macOS arm64, where the harmonypy source build
  succeeds because Accelerate supplies BLAS — so a clean local install proved
  nothing about Linux, and CI found it in 23 seconds. `--python-platform` +
  `--only-binary :all:` is the check that would have caught it; the hand-picked
  PyPI wheel survey that preceded it missed harmonypy entirely.
- **Open — the matrix no longer varies dependencies.** Dropping 3.10 removed
  something nobody had written down: because the floors are `>=`, that leg was
  resolving numpy 2.2.6 / scipy 1.15.3 / pandas 2.3.3 where the others got
  2.4.6 / 1.18.0 / 3.0.3, so it was the *only* coverage of an older scientific
  stack. All three legs now resolve one identical set — three interpreters, one
  dependency version. This surfaced through `test_leverage_scores_are_not_flat`,
  whose docstring had attributed its 3.10 failure to CPython when the real
  variable was the resolved numpy. The fix is a dedicated oldest-supported-deps
  leg (install against the declared floors, not the latest), not keeping an EOL
  Python around as an accidental proxy for it.
- **Triggers:** push to `main`, all PRs
- **`test` job:** `ruff check shanuz` and `mypy`, both advisory (pre-existing debt
  not yet cleared, so neither gates the build), then `pytest tests/ -q` with
  `--cov=shanuz` (80%); the coverage XML is uploaded as an artifact from the 3.12
  leg rather than sent to a third-party service
- ~~Still open: a dedicated `build` job and coverage reporting~~ — ✅ delivered:
  the **`build` job** runs `uv build` (sdist + wheel), `twine check` on both, then
  installs the wheel into a clean venv and imports it **from outside the source
  tree**, asserting `__version__` matches `pyproject.toml`. That last step covers
  two things the `test` job structurally cannot: the test job installs `-e .[all]`,
  so it never proves the *wheel* is complete or that a plain `pip install shanuz`
  works without the optional scientific stack; and `__version__` now comes from
  installed metadata, which only a real non-editable install exercises.

### Type annotations
- ~~Add `mypy` to the CI lint job~~ — ✅ delivered: advisory, on the same footing
  as ruff. `[tool.mypy]` in `pyproject.toml` pins `files = ["shanuz"]` so a bare
  `mypy` checks exactly what CI checks (ruff takes its scope from the command
  line, where `ruff check shanuz` = 75 and `ruff check .` = 163 are both correct
  and easy to confuse), and `ignore_missing_imports` silences the stub-less
  scientific stack — without it the run is 222 errors, 139 of them purely
  "this third party has no stubs".
- **Baseline to work down:** ~80 errors in 15 modules on default settings, ~540
  in 47 under `--strict`. Neither is a fixed target: both drift with the
  interpreter and with the dependency versions each resolves against the `>=`
  floors in `pyproject.toml` (the CI matrix spans 81–85 on default settings), so
  read them as a scale rather than a score. The package already ships `py.typed`,
  so these annotations are what a downstream `mypy` trusts.
- **Start here — 21 of the 81 are `[name-defined]`,** and they are real rather
  than pedantic: string annotations naming symbols that are never in module scope
  (`-> "plt.Figure"` across 17 plotting functions where `plt` is only imported
  inside function bodies; `graph.py`/`neighbor.py` annotating each other's
  classes to dodge a circular import). Runtime is unaffected — string annotations
  are never evaluated, which is why 444 tests pass — but `typing.get_type_hints()`
  raises `NameError` on them, which **blocks the documentation site below**, since
  mkdocstrings resolves annotations. Fix with `if TYPE_CHECKING:` imports; keep
  the lazy runtime imports exactly as they are (they keep matplotlib optional and
  the import graph acyclic).
- Add `from __future__ import annotations` to all modules (already done on some)
- Annotate all public function signatures (`mypy --strict` clean)

### Tutorial coverage — the R-fidelity net for everything after PR #10
- **Why this matters most.** The two real defects ever found in this port — the
  CLR margin inversion (#32) and the SCTransform model (#37) — were both caught by
  a tutorial with an R side-by-side, *not* by the test suite, which was green
  through both. Meanwhile 24 feature PRs landed after #10 (integration, reference
  mapping, sketching, lazy matrices, HTO/MULTI-seq, Mixscape, spatial) and **none
  of them has ever been compared to real Seurat** — their tests assert
  self-consistency on synthetic `default_rng` fixtures. At the start of the
  initiative 36 of 103 public exports appeared in a runnable tutorial; 67 did
  not. As of T-obj it is **81 of 104**. Closing the rest is the highest-leverage
  correctness work left.
- **Plan:** ~12 new side-by-side tutorials in four waves, one tutorial per PR, in
  the existing shape (`<name>_tutorial.py` + `<name>_verify.R` + `<name>.md` +
  `figures_<name>/`). Wave 1 = integration (ifnb), cell hashing (GSE108313),
  reference mapping (panc8), Mixscape (GSE153056). Wave 2 = cell-cycle / module
  scores (THP-1), dim-reduction extras (pbmc3k), leverage sketching (ifnb) and
  object internals (pbmc3k). Remaining topics: the DE-test suite (Wave 3 — needs
  the Bioconductor trio) and spatial/scale. **The visualization gallery is no
  longer worth its own tutorial**: after T-obj, `dot_plot` is the only plotting
  export not exercised somewhere, so it should be folded into an existing one.
- **Wave 0 — ✅ delivered (#38).** The data plumbing every side-by-side needs:
  `shanuz.datasets` loaders for the raw-source datasets, `tutorials/export_seuratdata.R`
  for the two SeuratData-only ones (`ifnb`/`panc8`, verified to round-trip R's
  counts exactly), and the R deps (`SeuratData` + `harmony`).
- **Wave 1 — ✅ complete.**
  - **T7 cell hashing — ✅ delivered (#39).** `hto_demux` / `multiseq_demux`
    on GSE108313 (`hashing_vignette.md`). Result: `HTODemux` is **99.81 %**
    call-concordant with Seurat on identical input — the first real-data
    confirmation of the CLR fix (#32) and `clara` default (#34). `MULTIseqDemux`
    lands at 94.67 %; the gap is a real KDE-implementation difference (scipy
    `gaussian_kde` vs R `density()` — bandwidth *and* grid), logged not papered
    over. No defect found — the demuxers hold up.
  - **T9 Mixscape — ✅ delivered (#40).** `calc_perturb_sig` / `run_mixscape`
    / `mixscape_lda` on GSE153056 (`mixscape_vignette.md`). On a shared
    variable-feature basis, per-cell class concordance is **97.45 %** (KO/NP/NT and
    the full `<gene> KO`/`NP` class) — all NT cells agree, the same 14 guides read
    zero-effect on both sides, strong IFN-γ hits ≥97 %. Divergence is isolated to
    the weak boundary guides (MYC/SPI1/BRD4/CUL3) where the EM mixture is
    init-sensitive — a real method-level residual on a far more stochastic pipeline
    than the demuxers, not a bug. No defect found.
  - **T6 integration — ✅ delivered (#41). First defects of the initiative.**
    `run_harmony` / `integrate_layers` on ifnb (`integration_vignette.md`). Harmony
    and CCA reproduce Seurat's batch mixing and cell-type recovery to **three
    decimals** (batch-mixing entropy py/R 0.991 & 0.990/0.991). But investigating
    RPCA surfaced **two real bugs**, both now fixed: (1) a crash on unequal batch
    sizes — the reciprocal-PCA MNN args were mis-ordered, `IndexError` whenever
    n_query > n_ref, masked by the balanced synthetic test fixtures — **fixed in
    #41**; (2) a 4× under-integration (batch-mix 0.222 vs 0.914) — **fixed in its
    own follow-up PR** by matching Seurat's reciprocal-PCA construction: per-object
    scaling, the `ReciprocalProject` SD/L2 embedding normalization, and disabling
    the RPCA anchor filter, lifting batch-mix to **0.867** (py ARI→celltype 0.444 →
    0.677). The old "quality, not count" read was wrong — it was mostly count
    (global scaling under-found anchors) plus the missing embedding normalization.
  - **T8 reference mapping — ✅ delivered (#43). Wave 1's last tutorial.**
    `find_transfer_anchors` / `transfer_data` / `map_query` / `project_umap` on
    panc8 (`refmap_vignette.md`). Transferring `celltype` from a CEL-seq2 reference
    to a SMART-seq2 query on a shared HVG basis: per-cell label concordance with
    Seurat is **98.71 %** (2,363 of 2,394 cells), and each tool is ~98.5 % accurate
    against the held-out truth (shanuz 0.9845, Seurat 0.9879). Every abundant type
    is ≥98 %; the whole error budget is the rare types (<10 reference cells) where
    both tools stumble the *same way* — a small single-tech reference's honest
    limit, not a divergence. **No defect found** — the transfer stack ports
    faithfully on its first real-data benchmark, the initiative's other valid
    outcome.
- **Wave 2 — in progress.**
  - **T-cc cell-cycle & module scoring — ✅ delivered (#44). Wave 2's first.**
    `add_module_score` / `cell_cycle_scoring` on THP-1 (`cellcycle_vignette.md`),
    a proliferating line with real S/G2/M populations (unlike resting PBMCs). On
    identical counts + shared resolved gene lists, per-cell **Phase concordance is
    96.62 %** (20,028/20,729 cells) and the S/G2M/module scores correlate at
    Pearson ≥ 0.998. Both functions sample control genes at random and NumPy's RNG
    is not R's, so the scores are not bit-identical by construction — the residual
    is purely that RNG (the discrete Phase is robust to it), the same "don't chase
    the RNG" story as `clara` (hashing) and the MULTI-seq KDE. No defect found.
  - **T-dr dim-reduction extras — ✅ delivered (#45). Two more defects.**
    `jack_straw` / `score_jackstraw` / `run_ica` / `run_tsne` on pbmc3k
    (`dimreduc_vignette.md`). ICA reproduces R's subspace (Hungarian-matched
    \|r\| **0.982**) and t-SNE preserves its input's neighbourhoods as well as R's
    does (30-NN retention 0.470 vs 0.477) — both clean. **JackStraw was not.**
    Two independent bugs, both fixed here: (1) the permutation null was built by
    projecting the scrambled rows onto the *fixed* embedding instead of re-running
    the PCA as R's `JackRandom` does, making the null far too tight — on the
    pure-noise PCs 14-20 that put **109-203** of 2000 features below p ≤ 1e-5
    where R finds **0-5**; (2) `score_jackstraw` aggregated with a one-sided KS
    test rather than R's `prop.test`, so its largest score across all 20 PCs was
    **8.1e-112** and no PC ever failed. Together they made shanuz keep **all 20**
    PCs where Seurat keeps 13. After the fix both keep **13**, and the residual is
    permutation scatter (13-15 across seeds; R fixes its per-replicate seeds and
    is deterministic). `JackStrawData.fake_reduction_scores`, declared but never
    populated, is filled too.
  - **T-sk leverage-score sketching — ✅ delivered. Two more defects.**
    `leverage_score` / `sketch_data` / `project_data` on ifnb
    (`sketch_vignette.md`), exercising **both** of Seurat's regimes by moving
    `nsketch` rather than the dataset. After the fixes the exact regime is a
    per-cell match (Spearman **1.000000**, max abs diff 3.4e-6 — below R's own
    unseeded-`irlba` noise), leverage tracks cell-type rarity at Spearman
    **−0.929** in both tools, and `project_data` reaches Seurat's accuracy
    exactly (**0.9050** each). **Two bugs, both fixed here:** (1) leverage was
    whitened against the *full rank* rather than Seurat's rank-50 truncation, so
    on 2000 HVGs the scores were crushed to a max/median of **1.34** against R's
    **6.48** — uniform sampling scores 1.00, so leverage sampling had become an
    expensive way to sample uniformly; (2) `project_data` transferred labels
    through the integration anchors, where `ProjectData` uses a weighted k-NN
    vote in the projected reduction — the anchor route scored *better* (0.936 vs
    0.905), which is why it survived, but it costs exactly what sketching exists
    to remove and is unusable at the scale the API targets. Fixing it moved the
    headline number **down** and per-cell agreement **up** to 98.1 %.
  - **T-obj object internals — ✅ delivered. Eleven more defects, the largest
    haul of the initiative.** `Cells`/`Features`, the layered assay, `Key`,
    `Embeddings`/`Loadings`/`Stdev`, `Graphs`, `FetchData`, `Idents`/`WhichCells`/
    `RenameIdents`/`subset`, and the command log, on pbmc3k
    (`objects_vignette.md`). The container rather than an algorithm, which makes
    it the sharpest net in the series: nothing here is stochastic, so **89 of the
    91 anchors are compared with no tolerance at all** — orders, names,
    dimensions, keys and non-zero counts either match or they do not. Coverage
    went from 36/103 exports at the start of the initiative to **81/104**.
    **Eleven bugs, all fixed here.** Five in the layered assay: `split`/`JoinLayers`
    was not a round trip, returning a layer named `joined` whose columns were in
    *split* order while the assay's own cell vector never moved — the right
    numbers in the wrong columns, silently misaligned against the metadata that
    indexes them, with every shape and checksum intact. The no-argument
    `join_layers()` additionally raised `ValueError` on any prepared assay, and
    `generics.split_layers` was declared but never registered for any type.
    Three in `FetchData`: `np.asarray` on a sparse matrix wraps it rather than
    densifying, so every one of 2,700 rows held a copy of the whole matrix — on
    the most-called accessor in Seurat, on the default assay class, guarded by a
    test that asserted only the column name and row count. Plus `PC_1` was
    unaddressable and an unqualified fetch read `counts` where R reads `data`.
    Three in the bookkeeping: `log_shanuz_command` had **zero call sites** so
    `obj.commands` was always empty against Seurat's five entries; no
    `orig.ident`; and `add_meta_data` rejected the plain vector R documents.
    **`join_layers`/`split_layers` had zero call sites and zero tests before
    this** — the defining feature of the v5 object model, never once run.
    **Open, deliberately not fixed here:** `find_neighbors` symmetrizes the kNN
    graph where Seurat's is directed (nnz 75,740 against exactly 2,700 × 20 =
    54,000; degree 20-83, mean 28.1), and `_build_snn` drops the self-edge Seurat
    keeps (~2,700 of the 4,044-edge gap). Both change what clustering consumes
    and need their own comparison.
- **Expect bugs, and read a mismatch as a bug report.** Wave 1 went T7, T9 and T8
  clean, while **T6 found the first two defects**, **T-dr the next two**,
  **T-sk two more** and **T-obj eleven** —
  exactly the point: a green synthetic suite (balanced batches, self-consistent
  fixtures) hid a crash, a 4× under-integration, a mis-specified permutation null,
  the wrong significance test, a flattened sampling weight and a label transfer
  through the wrong machinery — all of which one real Seurat comparison exposed
  immediately. JackStraw and leverage sketching are the sharpest cases: JackStraw's
  tests asserted only that signal genes score lower than noise genes, which stayed
  true while the function was recommending every PC; `leverage_score`'s asserted
  its own full-rank definition to six decimals, and its sampling fixture was too
  small to be in the algorithm's regime at all. The
  known-good tolerance is narrow — deterministic values match exactly, Louvain
  cluster counts drift ±1 — so anything outside that band gets investigated,
  not written up as an expected difference. This repo has twice let a real defect
  hide behind a documented "language difference" caveat.
- **A fix can make the headline number worse, and still be the fix.** T-sk's
  `project_data` bug scored *above* Seurat before it was corrected. Fidelity to
  the reference implementation is the goal; when a divergence flatters us, that is
  a reason to look harder, not to keep it. Where a residual really is RNG, prove
  it distribution-against-distribution over matched seeds — single-run pairs were
  actively misleading on T-sk's sketch composition.

### Documentation site
- **Do the type annotations first.** mkdocstrings resolves annotations, and
  `typing.get_type_hints()` currently raises `NameError` on the plotting and
  `graph`/`neighbor` signatures — the site would hit that on day one.
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

### Changelog — ✅ delivered
- **File:** [`CHANGELOG.md`](CHANGELOG.md) at repo root, in
  [Keep a Changelog](https://keepachangelog.com) format
- Populated retroactively for 0.1.0, 0.1.1 and 0.2.0 — but deliberately **not**
  straight from the git log as this item originally said: there is a
  `chore: bump version to 0.1.2` commit, yet 0.1.2 was never tagged and never
  reached PyPI, so it is not a release and gets no entry. Taking the log at face
  value would have documented a version that never existed. (0.1.0 was tagged but
  never published; 0.1.1 was the first release on PyPI.)
- Carries a standing note that the milestones on this page are **not** releases,
  because the two sequences have diverged far enough to mislead: the tags stop at
  0.2.0 while the milestones run to v0.9.0, and a milestone can straddle releases
  (v0.7.0's spatial loaders shipped in 0.1.1; the rest of v0.7.0 has not shipped).
- Everything since 0.2.0 sits under `[Unreleased]`, including the three entries
  queued by earlier PRs: #32's `BREAKING` CLR margin fix, #33's clara
  cross-architecture caveat, and #34's `hto_demux` default change.

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
2. ~~**WNN** (`v0.4.0`)~~ — ✅ delivered (`find_multi_modal_neighbors` + `run_umap(graph=)`); full two-stage port, CBMC tutorial section complete
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
