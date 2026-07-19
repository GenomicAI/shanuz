# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

> **Milestones are not releases.** [`ROADMAP.md`](ROADMAP.md) tracks progress in
> `v0.N.0` *milestones*, named after the slice of Seurat they port. Those are
> planning labels and never become version numbers on their own — the milestone
> sequence even runs `v0.9.0` → `v0.10.0` with nothing in between. The versions
> below are the ones actually tagged and released, and the two have drifted a
> long way apart: the tags stop at 0.2.0 while the milestones have run through
> v0.9.0. Everything in between sits under [Unreleased]. A milestone can also
> span releases — v0.7.0's spatial loaders shipped in 0.1.1 while the rest of it
> is still unreleased — so the two sequences do not line up item for item.

## [Unreleased]

Work from six milestones — reference mapping, extra reductions, pseudobulk DE,
spatial, scale, and the specialized assays — plus one breaking fix. All of it is
on `main`; none of it is on PyPI.

### Fixed

*Spatial statistics and the spatial container, audited against Seurat 5.5.1*

- **`find_spatially_variable_features(method="moransi")` used the wrong spatial
  weights.** Seurat builds `1/d²` between every pair of cells and
  `Rfast2::moranI` row-standardises it; shanuz used a k-nearest-neighbour graph.
  It was a *good* approximation — Pearson 0.986 against R, 46 of R's top 50
  genes — which is exactly why nothing caught it, but it ran a median 1.23× high
  and recovered only **7 of R's top 10**, the part of the ranking anyone reads.
  R's weighting is now the default and matches to **1.6e-14 with 10/10**,
  evaluated in row blocks so the n × n matrix is never materialised: the full
  36,602-cell slide runs in 5.3 s at 0.95 GB where `RunMoransI` needs a 10.7 GB
  allocation. `weights="knn"` keeps the old path, documented as an approximation.
  The p-value deliberately stays a normal approximation — R's 999-permutation
  test returns 14 distinct values and ties 233 of 248 genes at its floor.
- **`Centroids` never carried a radius.** `SeuratObject` always computes one
  (`.AutoRadius`, 1% of the mean bounding-box dimension — 42.83 on the Xenium
  mouse brain). shanuz left it `None`, and because `_spot_collection` returns
  `None` for a `None` radius, every true-to-scale spot renderer silently fell
  back to a fixed-size scatter on every FOV not built from a Visium
  `scalefactors_json.json`. `_spatial_panel` also read the radius off the FOV,
  where R keeps `NULL`; it now reads the default boundary, as R does.
- **`Segmentation` stored polygons open.** R closes each ring by repeating the
  first vertex — a square is five rows, not four. Now closed, idempotently, with
  concave shapes preserved.
- **The pbmc3k figure generator mislabelled every cell type.** Its hardcoded
  cluster→cell-type map had `1↔2` and `3↔4` transposed, putting monocyte names
  on the T-cell compartment in every labelled figure, including the annotated
  UMAP that heads the tutorial. The R code printed beside it in
  `pbmc3k_tutorial.md` had the correct order throughout. Corrected in both,
  figures regenerated, and now guarded by a test that checks each label against
  its own discriminative marker instead of trusting the map.

*The object model, audited against Seurat 5.5.1 (#48)*

Eleven fidelity defects, found by the first tutorial to compare the **container**
rather than an algorithm. `join_layers` / `split_layers` had zero call sites and
zero tests before this — the defining feature of the v5 object model, never once
run.

- **`split` / `JoinLayers` was not a round trip.** The join returned a layer
  named `joined` (Seurat restores the original name) whose columns were in the
  *split's* order rather than the assay's. The assay's own cell vector never
  moves during a split, so the matrix came back **silently misaligned against
  the metadata that indexes it** — ask for cell `c1`'s column, get `c2`'s. Every
  shape, sum and checksum was intact, because every value was still present.
  `join_layers()` with no arguments — the only call a real script makes — also
  raised `ValueError` on any prepared assay, hstacking `counts`, `data` and a
  variable-features-only `scale.data` together regardless of feature count. The
  fix groups layers by the stem they were split from and records that provenance
  at split time, because the name cannot be parsed back: Seurat's own
  `scale.data` contains the separator. Split parts are now named `counts.batch1`,
  Seurat's spelling, not `counts_batch1`.
- **`shanuz.generics.split_layers` was declared but never registered** for any
  type, so the documented generic raised `NotImplementedError` while the method
  it should have dispatched to worked fine.
- **`fetch_data` returned objects instead of numbers.** `np.asarray` on a sparse
  matrix yields a 0-d *object* array wrapping it, not its contents, so
  `.flatten()` broadcast one `csc_matrix` down every row — 2,700 copies of the
  whole matrix in place of 2,700 expression values, on the most-called accessor
  in Seurat and on the default assay class. Its test asserted the column name
  and the row count, both of which that satisfies.
- **`fetch_data` could not address an embedding column by its key.** `PC_1`
  raised `KeyError`; only the reduction name worked, and it emitted `pca_1`
  rather than the `Key()`-derived `PC_1` R uses. Both now work.
- **`fetch_data` read `counts` where R reads `data`**, returning raw integers
  where every vignette shows normalized expression. It now defaults to `data`
  and, when there is no `data` layer, falls back to `counts` *with a warning*,
  as Seurat does.
- **The command log was inert.** `log_shanuz_command` was a public export with
  no call sites, so `obj.commands` was always empty where Seurat logs one entry
  per pipeline step. `normalize_data`, `find_variable_features`, `scale_data`,
  `run_pca` and `find_neighbors` now log, keyed as Seurat keys them
  (`NormalizeData.RNA` … `FindNeighbors.RNA.pca`).
- **`orig.ident` was never created.** It is the first column of every Seurat
  object's metadata and the default identity class.
- **`add_meta_data` rejected a plain vector**, which is the form R's
  `AddMetaData` documents and the vignettes pass it.

*Leverage-score sketching: flattened sampling weights and the wrong label transfer (#46)*

- **`leverage_score` whitened against the full rank.** Seurat computes leverage
  from a **rank-50 truncated SVD** — `rowSums(V²)` over the leading 50 right
  singular vectors, so the scores sum to 50. Shanuz used every direction above a
  tolerance, which is the classical hat-matrix definition and equally defensible
  in the abstract, but useless on data of this shape: the scores sum to the rank
  and are capped at 1, so 2000 variable features over a few thousand cells crush
  every score towards `d/n`. On PBMC 3k that meant a max/median of **1.34**
  against Seurat's **6.48**, where uniform sampling scores 1.00 — leverage
  sampling had become an expensive way to sample uniformly, silently. Both
  regimes are now ported: the truncated SVD below `nsketch * 1.5` cells, and
  `CountSketch` → `QR` → `JLEmbed` above it, along with Seurat's `nsketch` bump
  and its "too slow" / "too square" guards. The exact regime now matches R
  per-cell (Spearman **1.000000**, max abs diff 3.4e-6).
- **`leverage_score` read the wrong layer.** The default was `"scale.data"`;
  Seurat scores the log-normalized `"data"`. Changed to match. `sketch_data`
  follows.
- **`project_data` transferred labels through integration anchors.** Seurat's
  `ProjectData` calls `TransferSketchLabels` — a weighted k-nearest-neighbour
  vote *inside the projected reduction*, with the sketch's own rows as the
  reference. The anchor route scored **better** on ifnb (0.936 against Seurat's
  0.905), which is why it survived review; it is still wrong, and it costs
  precisely what sketching exists to remove, so at the scale this API targets it
  is unusable rather than merely different. Now matches Seurat's mechanism and
  its accuracy exactly (**0.9050** each), at 98.1 % per-cell agreement on a
  shared sketch. Seurat's weight kernel is reproduced term-for-term from
  `FindWeightsC`.
- **`sketch_data` gained `method="Uniform"`**, as in Seurat — the control that
  makes "the sketch keeps rare cells" a meaningful claim.

  **Breaking:** `project_data` no longer accepts `seed` (the k-NN vote is
  deterministic), and no longer accepts a raw label array for `refdata` — like
  Seurat it takes a column name on the sketch, or `{new_col: sketch_col}`.
  `leverage_score`'s `eps` changed meaning from an SVD rank tolerance (1e-8) to
  Seurat's Johnson–Lindenstrauss distortion (0.5), and it gained `ndims`.

*JackStraw: a mis-specified permutation null and the wrong significance test (#45)*

- **`jack_straw` built its null against a fixed basis.** R's `JackRandom` permutes
  the selected features and **re-runs the whole PCA**, taking the null loadings
  from that refit basis; shanuz projected the permuted rows onto the *existing*
  embedding. A fixed basis cannot rotate to absorb the scrambled signal, so the
  permuted loadings came out too small and the null was far too tight. On pbmc3k's
  pure-noise PCs 14-20 that put **109-203** of 2000 features below p ≤ 1e-5, where
  R finds **0-5**. Now refits the PCA per replicate, as R does.
- **`score_jackstraw` used a KS test instead of `prop.test`.** R's `ScoreJackStraw`
  tests the count of features below `score.thresh` against the count expected under
  a uniform null; shanuz ran a one-sided Kolmogorov-Smirnov test against
  Uniform(0, 1), which on thousands of features is enormously more sensitive — its
  **largest** score across all 20 pbmc3k PCs was `8.1e-112`, so no PC ever failed
  the threshold. R's `prop.test` is now ported exactly (Yates-corrected two-sample
  chi-square), reproducing R to nine significant figures from 1e-143 to 1.0.
- **Net effect:** shanuz recommended keeping **all 20** PCs where Seurat keeps 13 —
  the function could not do the one thing it exists for. Both now keep **13**. The
  remaining spread is permutation scatter (13-15 across seeds; R fixes each
  replicate's seed to its loop index and is therefore deterministic).
- `JackStrawData.fake_reduction_scores` is now populated, as in R; `jack_straw`
  takes the reduction's stored feature loadings as the observed statistic (matching
  `Loadings(object[[reduction]], projected = FALSE)`) and raises if they are absent,
  rather than silently re-deriving them.
- Both defects were caught by the new R side-by-side, **not** by the test suite,
  which was green throughout: its only JackStraw assertion was that signal features
  score lower than noise features, which stayed true the whole time. Regression
  tests now pin the null's calibration on noise PCs and the aggregation's ability
  to reject, and both were mutation-tested against the old code.

### Added

*Spatial-statistics tutorial — Wave 2's last, and the wave's close*

- `tutorials/xenium_svf_tutorial.py` + `xenium_svf_verify.R` +
  `svf_vignette.md` + `figures_svf/`. Compares the spatial **container**
  (`load_xenium` / `create_fov` / `create_centroids` / `create_segmentation`
  against `LoadXenium` / `CreateFOV` / `CreateCentroids` / `CreateSegmentation`)
  and the one spatial **statistic** never checked against R,
  `find_spatially_variable_features`. The existing Xenium tutorial built its R
  side from `Read10X` plus a coordinate frame, so R never constructed an FOV and
  the whole boundary layer had gone uncompared.
- **38 of 39 anchors match Seurat exactly**, 32 of them with no tolerance at all.
  The one that differs is `GetTissueCoordinates`' shape — R returns `x, y, cell`
  as three columns, shanuz carries the cell as the DataFrame index.
- `weights=` on `find_spatially_variable_features`: `"inverse_square"` (R's, the
  new default) or `"knn"` (the previous behaviour, kept for very large slides).
- `dot_plot` folded into the pbmc3k gallery as `08b_marker_dotplot.png` — the
  last plotting export with no tutorial coverage. Drawing it is what exposed the
  cluster-label transposition fixed above.
- **Wave 2 is complete**: five tutorials, sixteen defects, against Wave 1's four
  and two.

*Object-internals tutorial — the container, side by side with Seurat (#48)*

- `tutorials/pbmc3k_objects_tutorial.py` + `pbmc3k_objects_verify.R` +
  `objects_vignette.md` + `figures_objects/`. The first tutorial to compare the
  **object model** rather than an algorithm: `Cells`/`Features`, the layered
  assay, `Key`, `Embeddings`/`Loadings`/`Stdev`, `Graphs`, `FetchData`,
  `Idents`/`WhichCells`/`RenameIdents`/`subset`, and the command log.
- Nothing in it is stochastic, so **89 of 91 anchors are compared with no
  tolerance at all** — orders, names, dimensions, keys and non-zero counts
  either match Seurat or they do not. The two exceptions are the fields that
  read a PCA, named explicitly rather than covered by a blanket rule.
- 43 tests: 25 on the tutorial's own helpers (including the comparison
  instrument itself — one that always agrees would make every number it prints
  worthless) and 18 regressions on the defects above. All mutation-tested.
- Tutorial coverage of public exports: **36/103 → 81/104**.

*Guards against supported-Python drift (#47)*

- Three tests in `tests/test_packaging.py` cross-checking the four places the
  supported-version decision is written down — `requires-python`, the trove
  classifiers, `[tool.ruff] target-version`, and the CI matrix. Nothing read
  those together before, and each is quiet when wrong in a different way: a
  stale classifier misinforms PyPI without breaking a build; a matrix that has
  moved above the declared floor stops testing the floor, which is the version
  most likely to break; a ruff target below the floor silently disables the lint
  the floor just earned.
- Mutation-tested in all four directions — drop a classifier, raise the floor,
  revert the ruff target, drop the lowest matrix leg — each caught by the
  specific guard(s) it should be, none firing indiscriminately.
- The matrix is parsed from the workflow YAML with a regex rather than PyYAML:
  PyYAML is not a declared dependency, so importing it would pass today and
  begin silently skipping the day that transitive edge disappears — the exact
  drift these tests exist to catch.

*Leverage-score sketching tutorial — side by side with R Seurat (#46)*

- `tutorials/sketch_vignette.md` with `ifnb_sketch_tutorial.py`,
  `ifnb_sketch_verify.R`, and `generate_sketch_plots.py` — `leverage_score`
  (`LeverageScore`), `sketch_data` (`SketchData`) and `project_data`
  (`ProjectData`) on ifnb, on a cell and feature basis shared with the R run.
- **First real-data fidelity result for all three** (synthetic fixtures only
  before), and it found the two defects above. Exercises **both** of Seurat's
  regimes on one dataset by moving `nsketch` rather than the data. Headline:
  exact-regime Spearman **1.000000** against R, leverage tracks cell-type rarity
  at Spearman **−0.929** in both tools (CD4 Naive T 0.76× → Eryth 2.89×), and
  `project_data` matches Seurat's label accuracy exactly.
- ifnb's 13 annotated types — 4,362 cells down to 55 — make the payoff directly
  measurable, against a same-size **uniform** sketch as the control. No synthetic
  fixture reproduces it: several were tried, and R agrees with shanuz on those to
  1e-5 while showing no enrichment either, because real rare types are
  transcriptionally extreme rather than merely scarce.
- The lazy on-disk `LazyMatrix` round-trip is checked too, but reported
  separately and **not** as a side-by-side — R's equivalent is BPCells, which is
  not installed here.

*Dimensional-reduction extras tutorial — side by side with R Seurat (#45)*

- `tutorials/dimreduc_vignette.md` with `pbmc3k_dimreduc_tutorial.py`,
  `pbmc3k_dimreduc_verify.R`, and `generate_dimreduc_plots.py` — `jack_straw` /
  `score_jackstraw` (`JackStraw`/`ScoreJackStraw`), `run_ica` (`RunICA`) and
  `run_tsne` (`RunTSNE`) on PBMC 3k, on a cell and feature basis shared byte-for-byte
  with the R run.
- **First real-data fidelity result for all four** (synthetic fixtures only before).
  After the JackStraw fixes above, both tools keep **13 PCs**. ICA recovers the same
  subspace — components are matched one-to-one by |Pearson r| with the Hungarian
  algorithm, since they are defined only up to sign and order, giving **0.982** mean
  matched |r|. t-SNE is compared on structure rather than coordinates (`Rtsne` is
  Barnes-Hut, shanuz calls scikit-learn): each embedding retains **0.470** / **0.477**
  of its PCA neighbourhoods.
- The comparison reports where the two PCA bases stop matching *in order* (PC 15 on
  pbmc3k) rather than a bare minimum correlation, so a permuted noise tail is not
  mistaken for a disagreeing basis — and so a per-PC finding is only claimed over
  the range where it is like-for-like.

*Cell-cycle & module-score tutorial — side by side with R Seurat (#44)*

- `tutorials/cellcycle_vignette.md` with `thp1_cellcycle_tutorial.py`,
  `thp1_cellcycle_verify.R`, and `generate_cellcycle_plots.py` — `add_module_score`
  (`AddModuleScore`) and `cell_cycle_scoring` (`CellCycleScoring`) on the
  proliferating THP-1 line (GSE153056), compared against Seurat on identical GEO
  counts and the same resolved S / G2M / module gene lists. **Opens Wave 2** of the
  tutorial initiative.
- **First real-data fidelity result for the scoring features** (synthetic fixtures
  only before): per-cell `Phase` concordance with Seurat is **96.62 %** (20,028 of
  20,729 cells), and the `S.Score` / `G2M.Score` / module scores correlate at
  Pearson ≥ 0.998. Both functions sample control genes at random and NumPy's RNG is
  not R's, so the scores are not bit-identical *by construction* — the residual is
  that control-gene RNG (the discrete `Phase` is robust to it), the same documented
  behaviour as `clara` (hashing) and the MULTI-seq KDE. **No defect found.**
- 11 network-free unit tests (`tests/test_cellcycle_tutorial.py`) covering the
  metric helpers and a synthetic run with planted S/G2M populations, plus a gated
  real-data regression in `tests/test_tutorial_smoke.py`.

*Reference mapping tutorial — label transfer, side by side with R Seurat (#43)*

- `tutorials/refmap_vignette.md` with `panc8_reference_mapping_tutorial.py`,
  `panc8_reference_mapping_verify.R`, and `generate_refmap_plots.py` — the
  reference-mapping workflow (`find_transfer_anchors` / `transfer_data` /
  `map_query` / `project_umap`) on the panc8 pancreatic-islet atlas (Baron et al.
  2016), annotating a SMART-seq2 query from a CEL-seq2 reference. Both tools read
  identical exported counts and a shared variable-feature basis; the query's true
  `celltype` is held back as ground truth so the transfer is scored for accuracy,
  not just agreement with R. A single-technology reference is used deliberately, to
  isolate the transfer machinery from the integration path.
- **First real-data fidelity result for the reference-mapping features** (only
  synthetic two-type fixtures before): per-cell label concordance with Seurat is
  **98.71 %** (2,363 of 2,394 query cells get the same `predicted.id`), and each
  tool is ~98.5 % accurate against the held-out cell types (shanuz 0.9845, Seurat
  0.9879). Every abundant cell type is recovered at ≥98 %; the rare types (<10
  reference cells) are noisy in both tools alike — a small single-tech reference's
  honest limit, not a divergence. **No defect found** — the transfer stack ports
  faithfully. Completes Wave 1 of the tutorial initiative.
- 12 network-free unit tests (`tests/test_refmap_tutorial.py`) covering the metric
  helpers and a synthetic two-technology end-to-end run, plus a gated real-data
  accuracy regression in `tests/test_tutorial_smoke.py`.

*Integration tutorial — Harmony / CCA / RPCA, side by side with R Seurat (#41)*

- `tutorials/integration_vignette.md` with `ifnb_integration_tutorial.py`,
  `ifnb_integration_verify.R`, and `generate_integration_plots.py` — the three
  batch-integration paths (`run_harmony` / `integrate_layers(method="cca"|"rpca")`)
  on the ifnb IFN-β benchmark (Kang et al. 2018), compared against Seurat on
  identical exported counts and a shared variable-feature basis. The concordance is
  partition-based (cluster ARI, cell-type-recovery ARI, batch-mixing entropy) since
  integration embeddings are not coordinate-comparable across tools.
- **First real-data fidelity result for the integration features** (v0.2.0; only
  synthetic-fixture tests before): Harmony and CCA reproduce Seurat's batch mixing
  and cell-type recovery to three decimals (batch-mixing entropy py/R 0.991 and
  0.990/0.991). **The first tutorial in the initiative to find real defects** —
  two RPCA bugs (see *Fixed*): a crash on unequal batch sizes, and a ~4×
  under-integration; both are now fixed (the under-integration in follow-up #42).
- Added to the opt-in tutorial smoke suite (`SHANUZ_TUTORIAL_SMOKE=1`) and covered
  by `tests/test_integration_tutorial.py` (network-free: the silhouette/ARI/entropy
  helpers and the load→prep→integrate→score→concordance path on a synthetic
  two-condition dataset with *unequal* batch sizes). Suite 507 → 522. No new `pip`
  dependencies (the R reference uses the already-listed `harmony` package).

*Mixscape tutorial — `CalcPerturbSig` + `RunMixscape` + `MixscapeLDA`, side by side with R Seurat (#40)*

- `tutorials/mixscape_vignette.md` with `thp1_mixscape_tutorial.py`,
  `thp1_mixscape_verify.R`, and `generate_mixscape_plots.py` — the pooled-CRISPR
  Mixscape workflow (`calc_perturb_sig` / `run_mixscape` / `mixscape_lda`) on the
  THP-1 ECCITE-seq screen (GSE153056), compared call-for-call against Seurat on the
  same GEO bytes and a shared variable-feature basis.
- **First real-data fidelity result for the Mixscape features** (all of which
  landed after PR #10 with only synthetic-fixture tests): per-cell class
  concordance is **97.45 %** for both the global class (KO/NP/NT) and the full
  `<gene> KO`/`NP` class. All NT cells agree, the same 14 guides read zero-effect on
  both sides, and the strong interferon-γ hits agree ≥97 %; the residual is
  isolated to the weak boundary guides (MYC/SPI1/BRD4/CUL3) where the EM mixture is
  init-sensitive — a genuine method-level difference (scipy `GaussianMixture` vs R
  `mixtools`, plus per-gene DE tie-breaking), documented in the walkthrough. No
  defect found on a far more stochastic pipeline than the hashing demuxers.
- Added to the opt-in tutorial smoke suite (`SHANUZ_TUTORIAL_SMOKE=1`) and covered
  by `tests/test_mixscape_tutorial.py` (network-free: the perturbation-table /
  concordance helpers and the load→signature→classify→LDA path on a synthetic screen
  with known KO truth). Suite 496 → 507. No new `pip` dependencies (the R reference
  adds the `mixtools` CRAN package for `RunMixscape`).

*Cell-hashing tutorial — `HTODemux` + `MULTIseqDemux`, side by side with R Seurat (#39)*

- `tutorials/hashing_vignette.md` with `pbmc_hashing_tutorial.py`,
  `pbmc_hashing_verify.R`, and `generate_hashing_plots.py` — demultiplexing the
  8-hashtag Cell-Hashing dataset (GSE108313) with `hto_demux` / `multiseq_demux`,
  compared call-for-call against Seurat on byte-identical GEO input.
- **First real-data fidelity result for the hashing features** (all of which
  landed after PR #10 with only synthetic-fixture tests): `HTODemux` is **99.81 %**
  call-concordant with Seurat for both the global class and the sample
  assignment — confirming the CLR-margin fix (#32) and the `clara` default (#34).
  `MULTIseqDemux` agrees on 94.67 %; the residual is a genuine KDE-implementation
  difference (scipy `gaussian_kde` vs R `density()` — bandwidth *and* grid, which
  a single `nrd0` swap makes worse, not better), documented in the walkthrough.
- Added to the opt-in tutorial smoke suite (`SHANUZ_TUTORIAL_SMOKE=1`) and covered
  by `tests/test_hashing_tutorial.py` (network-free: species/concordance helpers
  and the load→demux→figure path on a synthetic barnyard). Suite 489 → 496. No new
  `pip` dependencies.

*Tutorial data infrastructure — the R-side scaffolding for expanding tutorial coverage (#38)*

- `shanuz.datasets.pbmc_hashing` (GSE108313) and `thp1_eccite` (GSE153056) —
  loaders for the Cell-Hashing and ECCITE-seq/Mixscape datasets, parsed straight
  from their original GEO plain-text files so R and Python read identical counts.
  The ECCITE loader also returns the per-cell guide/replicate metadata, so a
  Mixscape tutorial can start from the same annotated state as R's `thp1.eccite`.
- `shanuz.datasets.ifnb` and `panc8`, with `tutorials/export_seuratdata.R` — a
  one-time R bridge that exports the curated SeuratData objects (which have no
  clean cross-language raw source) to a gzipped 10x folder that `read_10x` reads.
  Verified end to end: R exported panc8 and Python read back **51,767,089**
  nonzeros — matching R to the entry — with barcodes and metadata aligned.
- No new `pip` dependencies (the loaders are pure pandas/scipy). R side adds
  `SeuratData` + `harmony`. Very wide count tables are parsed once and memoised to
  a `.npz` sidecar, so a repeat load is ~0.2s rather than minutes.

*Reference mapping and label transfer (milestone v0.3.0)*

- `find_transfer_anchors` and `transfer_data` — atlas-based annotation with
  `pcaproject`/`cca` reduction and both classification and imputation. (#22)
- `map_query` and `project_umap` — place a query dataset in the reference's
  existing UMAP. (#23)

*Integration (milestone v0.2.0, completing it)*

- `integrate_data` and `integrate_layers` — CCA/RPCA anchor-based integration,
  alongside the Harmony path released in 0.2.0. (#21)

*Dimensionality reduction (milestone v0.5.0)*

- `run_spca` — supervised PCA. (#19)
- `glm_pca` — GLM-PCA with Poisson (#19) and negative-binomial (#20) families.
  Pure NumPy/SciPy; the `glmpca-py` dependency proved unnecessary.

*Pseudobulk and differential expression (milestone v0.6.0)*

- `aggregate_expression` and `find_conserved_markers`. (#10)
- `find_markers(test_use="deseq2")` — pseudobulk DESeq2 via optional
  `pydeseq2`. (#11)
- `find_markers(test_use="mast")` — MAST two-part hurdle test. (#12)
- `find_markers(test_use="bimod")` — the McDavid 2013 likelihood-ratio
  test. (#13)

*Spatial transcriptomics (milestone v0.7.0)*

- `load_merscope` — Vizgen MERSCOPE loader. (#14)
- `find_spatially_variable_features` — Moran's I (#15) and markvariogram (#18).
- `VisiumV2` and `load_visium(image=)` — the tissue-image data layer. (#16)
- `spatial_dim_plot` and `spatial_feature_plot` — H&E overlays. (#17)

*Scale and performance (milestone v0.8.0)*

- `sketch_data`, `project_data`, and `leverage_score` — leverage-weighted
  subsampling for million-cell datasets. (#24)
- `LazyMatrix` — BPCells-style on-disk matrices built on NumPy memory-mapping,
  no new dependency. (#25)

*Specialized assays (milestone v0.9.0)*

- `hto_demux` — `HTODemux` cell-hashing demultiplexing. (#26)
- `multiseq_demux` — MULTI-seq demultiplexing. (#27)
- `calc_perturb_sig` and `run_mixscape` — pooled-CRISPR analysis. (#28)
- `mixscape_lda` — supervised map separating guide populations. (#29)
- `plot_perturb_score` and `mixscape_heatmap`. (#30)
- `shanuz._clara` and `hto_demux(kfunc="clara")` — an in-tree port of R's
  `cluster::clara` k-medoids, which is what `HTODemux` actually uses. Needs no
  sklearn. (#33)

  **Known caveat:** R's `clara` is *not* reproducible across CPU architectures —
  it accepts swaps on any improvement below zero, so one ulp in one distance
  flips the whole clustering, and `clara.c` fuses a multiply-add on arm64 that
  it rounds twice on x86-64. "Match R" is therefore not well-defined. This port
  deliberately follows plain IEEE double arithmetic (= `clara.c` on x86-64,
  and what NumPy gives everywhere), and is exact against that reference. It can
  disagree with an arm64 R build, by design.

### Changed

- **Supported Python is now 3.12–3.13; 3.10 and 3.11 are dropped.** The CI matrix
  moves with it, and `requires-python` becomes `>=3.12`.

  The floor tracks [SPEC 0](https://scientific-python.org/specs/spec-0000/) — the
  three-years-past-release window numpy, scipy, pandas and scikit-learn keep for
  themselves — rather than CPython's five-year EOL. By that rule 3.10 lapsed in
  Oct 2024 and 3.11 in Oct 2025, both already past; going by EOL alone would have
  held 3.11 until Oct 2027. Those four packages are what actually constrain this
  library, so theirs is the calendar worth following.

  The full suite passes identically on both legs — 616 passed / 17 skipped — as
  do all 17 tutorial smoke tests, which CI skips on every leg, run here against
  the real datasets.

  **Python 3.14 is not included, though it very nearly works.** Every package in
  the set has a cp314 manylinux wheel except `harmonypy`, which publishes
  cp39–cp313 only. Without a wheel, uv builds it from source, and that needs BLAS
  plus a CMake-fetched armadillo `ubuntu-latest` does not have. Forcing a
  wheels-only resolve is worse: it backtracks to `harmonypy` 0.2.0, which depends
  on torch and pulls in triton and 24 `nvidia-*` packages. Dropping the
  `integration` extra on a 3.14 leg does resolve clean (95 packages, wheels only)
  but would leave 18 harmony tests unrun there. Deferred until `harmonypy` ships
  a cp314 wheel; see `ROADMAP.md`.

  **Nothing breaks retroactively.** `pip` on 3.10 or 3.11 resolves to 0.2.0, the
  last release declaring `>=3.10`.

  *Also removed:* the `tomli` dev dependency (`tomllib` is stdlib from 3.11) and
  its import fallback in `tests/test_packaging.py` — the only version-gated code
  in the repo. `[tool.ruff] target-version` moves to `py312`.

- **`hto_demux` now defaults to `kfunc="clara"`**, matching Seurat; it first
  shipped defaulting to `"kmeans"`. Callers who never passed `kfunc` get
  different output: ~1% of cells change class on synthetic panels, rising with
  tag count to ~3.5% at 12 tags — where `clara` is also the *more* accurate of
  the two. Accuracy is otherwise a wash, so this is a fidelity change, not a
  quality one. `"kmeans"` remains available. Both scale linearly in cells;
  `clara` costs a roughly constant 4× (~1.3 s vs ~0.3 s at 100k cells), so
  choose on fidelity rather than speed. (#34)
- `find_multi_modal_neighbors` is now a full two-stage port of Seurat's WNN
  (`FindModalityWeights` + `MultiModalNN`), replacing an approximation that used
  a linear distance ratio and blended per-modality SNN graphs instead of doing a
  joint neighbour search. The old formula was monotone in the right quantity but
  had no dynamic range, pinning every weight near 0.5 — a weight stuck at 0.5
  cannot say "this cell is decided by protein", which is the one thing WNN
  exists to say. On synthetic data where one group separates only in RNA and
  another only in ADT, the port now gives them ADT weights of 0.073 and 0.993
  (previously 0.482 and 0.575). (#31)
- `hto_demux` and `multiseq_demux` default to `margin=1` instead of `margin=2`.
  **Their behaviour is unchanged** — this compensates for the `_clr_normalize`
  fix below, which would otherwise have silently broken both. Note that
  `margin=1` for hashing is deliberate and correct despite Seurat's *ADT* advice
  being `margin=2`: hashing wants per-hashtag-across-cells, which is what
  Seurat's hashing vignette does at its own default. (#32)

### Fixed

- `integrate_layers(method="rpca")` **under-integrated ~4×** versus Seurat's RPCA
  on real data — on the ifnb benchmark it reached batch-mixing entropy 0.222
  against Seurat's 0.914, with cell-type recovery (0.444) *below* the uncorrected
  baseline, while `reduction="cca"` and `run_harmony` matched Seurat to three
  decimals. Reading the real Seurat source (`ReciprocalProject`, `FindNN`) against
  an anchor-count probe showed the reciprocal-PCA path diverging three ways: it
  scaled the batches **globally** instead of per-object (Seurat's `SplitObject →
  ScaleData` per object), leaving each batch's mean shift in PC1 so reciprocal PCA
  under-found mutual pairs; it searched the **raw** projection instead of Seurat's
  `l2.norm` (standardise each dimension by its SD, then L2-normalise each cell), so
  PC1's variance swamped the neighbour search; and it applied the expression-space
  anchor **filter** Seurat disables for reciprocal PCA. Fixing all three lifts RPCA
  batch-mixing to **0.867** and cell-type recovery to **0.677** (now above
  baseline), with CCA and Harmony unchanged. The residual to 0.914 is the expected
  exact-vs-annoy-neighbour / scikit-learn-vs-irlba-PCA gap. Regression tests: a unit
  test of the embedding normalisation, a check that the RPCA weight embedding is
  L2-normalised (the fix's observable signature — pre-fix rows were 0.79–0.92), and
  a gated ifnb batch-mixing floor, since the *emergent* under-integration reproduces
  on no synthetic fixture (both a 3-type and a hard 6-type unequal-batch fixture
  integrate fine on the pre-fix code). Completes the RPCA pair found in #41. (#42)

- `find_integration_anchors(reduction="rpca")` crashed with `IndexError` on any
  pair of datasets with **unequal cell counts** — i.e. every real dataset (the
  ifnb benchmark is CTRL 6,548 vs STIM 7,451). The reciprocal-PCA branch passed its
  mutual-nearest-neighbour helper the reference/query projections in the wrong
  order, so the query-neighbour list was sized to the reference and indexed past
  its end whenever the query was larger. Balanced synthetic fixtures (equal batch
  sizes) never tripped it. Fixed by restoring the argument order, with two
  regression tests over unequal-size batches (both orderings). Found while building
  the ifnb integration tutorial (#41). *A separate, deeper RPCA under-integration
  (~4× vs Seurat's RPCA) was found at the same time and is **fixed in #42** — see
  the entry above; `reduction="cca"` and `run_harmony` were unaffected throughout.*

- **BREAKING** — `normalize_data`'s CLR `margin` argument was inverted relative
  to Seurat. `margin=1` is now per-feature across cells (Seurat's default) and
  `margin=2` is per-cell across features (what ADT panels want), matching the
  axis R's `CustomNormalize` passes to `apply(data, MARGIN = margin, ...)`.
  Verified against R: shanuz `margin=2` reproduced Seurat `margin=1` and vice
  versa, agreeing to 5e-6. Only the axis was wrong — the per-vector kernel was
  always exact.

  *Who is affected:* callers passing `margin` explicitly to `normalize_data`,
  `hto_demux`, or `multiseq_demux`. Callers using the defaults are unaffected.

  This was the sole cause of the CBMC tutorial's `ADT.weight` gap against
  Seurat; eight of nine cell types now match to 0.02 or better. (#32)

- `sctransform`'s regularized NB model was wrong in four places, and the errors
  compounded into a normalization that erased the fine cell subsets SCTransform
  exists to resolve. A method-of-moments estimator stood in for `theta.ml`; the
  regularization was smoothed against the **arithmetic** gene mean where R uses
  the **geometric** mean, and targeted `log(theta)` rather than the
  overdispersion factor; and residual variance — which ranks the variable
  features — was computed from residuals clipped at `sqrt(N/30)` where Seurat
  clips at `sqrt(N)`, applying the tighter clip only to the stored `scale.data`.
  Verified against a live Seurat 5.5.1 / sctransform 0.4.3 run on PBMC 3k, the
  regularized intercept now matches R at Spearman 1.0000, theta at 0.96, and
  residual variance at 0.9986, with 2,913 of 3,000 variable features shared —
  previously those were 1.0000, **−0.89**, **−0.07** and **414 of 3,000**. (#37)

  *What it looked like:* the SCTransform tutorial resolved **9** clusters against
  R's 12 — and, the real tell, *fewer* than the 11 from plain log-normalization,
  inverting the vignette's whole claim that SCTransform resolves finer subsets.
  It now resolves 13 with 4 T-cell subsets against log-normalization's 11 and 2,
  recovering the pDC, CD8 naive/memory and interferon-response populations. The
  Poisson GLM was never at fault — its intercept and slope always matched R
  exactly; only what was built on top of it was wrong.

  *Also:* `sctransform` now takes `vst_flavor`, defaulting to `"v2"` as Seurat 5
  does (depth slope fixed at `log(10)`, non-overdispersed genes modelled as pure
  Poisson, a variance floor), with `"v1"` for the original 2019 model. Under
  `"v1"` Python and R both resolve 13 clusters; at the `"v2"` default Shanuz
  resolves 13 to R's 12, a real one-cluster difference (R is stable at 12 across
  seeds) left by `vst`'s random step-1 gene sample and the different clustering
  libraries. The assay's `meta_data` now also carries `gmean`.

  *Why it went unseen:* nothing compared the model against R. The tutorial's
  cluster count was documented as an expected implementation difference and
  carried a ⚠️ in `sctransform_vignette.md`, which made a real defect look like a
  known caveat — and `tutorials/README.md` claimed an "exact match" on the 3,000
  variable features when 13.8% of them agreed. `tests/test_sctransform_r_fidelity.py`
  now pins each numerical primitive against R directly, including a port of
  `bw.SJ` (Sheather–Jones; SciPy has no equivalent) that matches R to 3e-7.

- `tutorials/pbmc3k_tutorial.py` — the tutorial the README sends new users to
  first — crashed with `KeyError: 'cluster'` on every fresh install. pandas 3
  stopped passing the grouping column into the callable of
  `groupby(...).apply(...)`, so the top-markers table came back without the
  column the next line filtered on. It now builds the table per cluster and runs
  on pandas 2.0 through 3.x. (#36)

  *Why it went unseen:* the old code works on pandas 2, `pyproject` declares
  `pandas>=2.0`, and no test executed a tutorial — so a developer venv holding
  pandas 2 passed while a fresh `pip install` resolved pandas 3 and broke. The
  full suite passed on the very install where the tutorial died. Tutorials now
  have execution coverage: `tests/test_tutorial_marker_tables.py` runs in CI, and
  `tests/test_tutorial_smoke.py` runs each tutorial end-to-end behind
  `SHANUZ_TUTORIAL_SMOKE=1`.

  Note for anyone touching the plot generators: they use the same
  `groupby(...).apply(...)` idiom and were **not** affected — they never read the
  dropped column. They were rewritten to match anyway, preserving output exactly.
  The obvious rewrite (`sort_values(...).groupby(...).head(n)`) is wrong there: it
  returns the same genes interleaved across clusters, silently scrambling
  `DoHeatmap`'s per-cluster blocks. `test_top_genes_is_cluster_major` pins it.

### Documentation

- CBMC CITE-seq tutorial: Step 8's WNN section written against real figures for
  the first time, which is what exposed both the WNN approximation and the CLR
  margin bug. (#31, #32)
- Side-by-side R Seurat plots added to the PBMC 8k and CBMC tutorials, and
  misleading R/Python figure comparisons corrected. (#8, #9)

## [0.2.0] - 2026-07-05

First release with batch correction, and the last release to date.

### Added

- `run_harmony` and `integrate_layers` — Harmony batch correction. (#6)
- `find_multi_modal_neighbors` and `run_umap(graph=)` — WNN. (#6)
  *Superseded:* see the full WNN port under [Unreleased].
- `run_ica` and `run_tsne` — additional reductions. (#6)

### Fixed

- Spatial tutorial figures are written next to the script rather than into the
  working directory, so the tutorial is safe to run standalone. (#5)
- README links use absolute GitHub URLs, so they resolve on the PyPI page. (#7)

## [0.1.1] - 2026-07-04

First release published to PyPI — `pip install shanuz` works from here on.

### Added

- Spatial Seurat parity: loaders, neighborhood analysis, niches, and
  composition. (#4)
- Xenium spatial tutorial, verified against R Seurat to 8 significant
  figures. (#4)
- [`ROADMAP.md`](ROADMAP.md) — the milestone plan. (#2)
- `shanuz/py.typed` — the package ships PEP 561 type information, so a
  downstream `mypy` reads shanuz's annotations rather than treating it as
  untyped. (#4)

### Fixed

- Pin `numba>=0.59` so the `[analysis]` and `[all]` extras resolve on Python
  3.10+. (#2)
- The README Quick Start example now produces the 500 cells it claims. (#3)

## [0.1.0] - 2026-06-30

Initial release: a port of Seurat's core data structures and analysis pipeline.
Tagged but never published to PyPI; `0.1.1` was the first release on PyPI.

### Added

- Core objects: `Shanuz`, `Assay`, `Assay5`, `StdAssay`, `DimReduc`, `Graph`,
  `Neighbor`, `LogMap`, `JackStrawData`, `ShanuzCommand`.
- Pipeline: `normalize_data`, `find_variable_features`, `scale_data`,
  `percentage_feature_set`, `run_pca`, `find_neighbors`, `find_clusters`,
  `run_umap`, `find_markers`.
- `sctransform`, `module_score`, `jack_straw`.
- Spatial primitives: `FOV`, `Centroids`, `Segmentation`, `Molecules`.
- Plotting, I/O, AnnData interop, and the bundled example datasets.

[Unreleased]: https://github.com/GenomicAI/shanuz/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/GenomicAI/shanuz/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/GenomicAI/shanuz/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/GenomicAI/shanuz/releases/tag/v0.1.0
