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

### Added

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
