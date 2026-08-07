# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

> **Milestones are not releases.** [`ROADMAP.md`](ROADMAP.md) tracks progress in
> `v0.N.0` *milestones*, named after the slice of Seurat they port. Those are
> planning labels and never become version numbers on their own. Before 0.9.0
> the two sequences had drifted a long way apart — the tags stopped at 0.2.0
> while the milestones had run through v0.9.0, with a milestone able to span
> releases (v0.7.0's spatial loaders shipped in 0.1.1 while the rest of it was
> still unreleased). The 0.9.0 release below closes that entire gap in one
> jump, which is why its version number happens to match the milestone it
> completes — a coincidence of this one release, not a policy of matching them
> going forward.

## [Unreleased]

## [1.1.0] - 2026-08-08

Ten weeks of plotting work, a domain-expert review, and the release-readiness
run that found the defect below. The API additions are backwards compatible;
the fixes change output, which is the point of them.

**If you are on 1.0.0, the reasons to move are the corrections**, not the
features: the categorical palette could give two clusters the same colour,
`vln_plot` drew the wrong violin three ways, `aggregate_expression` left raw
sums where Seurat log-normalizes, and `_get_expression_matrix` returned the
wrong layer with mislabelled features.

### Added

- **CI runs the PBMC 3k tutorials against real data.** A new `tutorials` job
  caches the 28 MB dataset and runs nine smoke tests against it — four tutorial
  scripts end to end (guided, SCTransform, dim-reduc, objects), plus the marker
  table, the object-model round trip, the dim-reduc extras and both
  cell-type-map guards. Until now `test_tutorial_smoke.py` ran nowhere but a
  developer's machine, and the entry above is what that cost.

  `lazy_bpcells` and `pbmc3k_de` are deliberately **not** in it. Both pass and
  both have their data cached; they were 65% of the runtime (136s and 102s of
  364s locally, and the first CI run measured **16m33s** for all eleven), and
  neither draws a labelled figure. That is a cost trim rather than a coverage
  judgement, and it is the weaker half of this change: the full opt-in suite is
  the only thing that runs them, which is why the pre-release run stays a
  requirement.

  Two details are the point rather than decoration. **A skip fails the job:**
  every selected test needs only the cached dataset, so a skip can only mean a
  broken cache, a misspelt `TRUECELL_TUTORIAL_SMOKE`, or a new test wanting a
  dataset this job does not fetch — all of which would otherwise produce a green
  tick for a job that verified nothing. And the selection **deselects the eight
  uncached datasets** instead of selecting pbmc3k by name: an allow-list reads
  better but would let a pbmc3k test added later under another name silently
  never run, which is precisely the failure being fixed.

  The other datasets total ~200 MB and stay developer-run, so
  `TRUECELL_TUTORIAL_SMOKE=1 pytest tests/test_tutorial_smoke.py` before a
  release is still the rule, not a formality.

### Fixed

- **The annotated PBMC 3k UMAP captioned the platelet cluster "DC".** Shipped in
  1.0.0. `tutorials/generate_plots.py` carried Seurat's nine `new.cluster.ids`
  while this pipeline resolves eight — DC merges into CD14+ Mono at resolution
  0.5 — and `rename_idents` is positional, so every label from position 7 on
  slid by one. The 14-cell platelet cluster (PPBP 5.85, FCER1A 0.00) was
  labelled DC, and "Platelet" was never applied to anything.

  What makes this one worth reading twice: **the rest of the repository already
  had it right.** `pbmc3k_tutorial.md` prints the corrected eight-entry map
  directly beside the figure drawn from the stale one, and `tutorials/README.md`
  states that all 32 DC cells land in CD14+ Mono. One file lagged, and it was
  the file that draws the pictures. Seven figures are regenerated here — the
  label set also changes the alphabetical palette assignment, so the two
  cluster-coloured panels move too. `10_marker_heatmap.png` correctly does not:
  it restores cluster numbers before drawing.

  Two prose claims counted the same nine and are corrected here: the caption
  directly beneath the figure said *"the biological result — 9 identical cell
  types — is the same"*, and `tutorials/README.md` said the tutorial annotates
  nine. Both contradicted Step 11 and the fidelity table on the same page. The
  caption also attributed the differing cluster numbering to Louvain being
  "non-deterministic by cluster ID"; both tools in fact number clusters by
  descending size (R: 711, 478, 471, 344, 270, 164, 154, 32, 14 — truecell:
  692, 515, 458, 344, 301, 159, 155, 14), so the keys differ because the sizes
  do. R's cluster 7 is the 32-cell DC and its 8 the 14-cell platelets; truecell's
  7 is those same platelets, which is the whole defect in one line.

  The guard test existed and was correct. It never ran: `test_tutorial_smoke.py`
  is opt-in and excluded from CI by design, so nothing failed. It now also
  asserts that the map's keys are exactly the cluster ids produced, which
  catches a length mismatch directly instead of inferring it from a marker
  landing in the wrong place — the old assertion reported this defect as
  "FCER1A is highest in CD14+ Mono", which describes the clustering rather than
  the labelling and reads like a false alarm.

- **The PBMC 3k guided tutorial documented a clustering ARI of 0.938 while
  measuring 0.899**, across six files (`tutorials/README.md` ×2,
  `pbmc3k_tutorial.md` ×2, `docs/fidelity.md`, `docs/quickstart.md`). The
  associated concordance figures were stale too — 2,554/2,638 cells and 0.968,
  against a measured 2,519 and 0.955.

  The likely cause is the graph fixes in #67–#71, which moved cells between
  clusters — the same drift the DE tutorial's `deseq2 top50` band caught at the
  time (25 → 22). This tutorial had **no band on its headline number**, so it
  went stale in six documents instead of failing once. Every swept resolution
  now carries a declared band, checked by `--report`.

  Verified this predates the change: R's vector-form `FindClusters` produces a
  partition identical to the single 0.5 call, and the tutorial's four handoff
  outputs are byte-identical to what `main` produces.

### Changed

- **The guided clustering tutorial scans four resolutions instead of one.**
  Choosing a resolution means running a few and comparing them, so a fidelity
  claim pinned to a single setting says less than it appears to. Both sides now
  use the vector idiom — `find_clusters(pbmc, resolution=[0.4, 0.8, 1.2, 0.5])`
  and `FindClusters(pbmc, resolution = c(...))` — and every resolution is scored
  against R:

  | resolution | truecell | Seurat | ARI |
  |---:|---:|---:|---:|
  | 0.4 | 9 | 9 | 0.8958 |
  | 0.5 | 8 | 9 | 0.8987 |
  | 0.8 | 11 | 11 | 0.8264 |
  | 1.2 | 12 | 12 | 0.7995 |

  Two readings. The **cluster count matches exactly at 0.4, 0.8 and 1.2**, so
  the 8-vs-9 split this tutorial describes is specific to resolution 0.5 rather
  than a standing property of the port. And agreement falls as resolution rises,
  which is expected: finer partitions put more boundaries in play.

  0.5 is given **last** deliberately, since Seurat leaves the object on the last
  resolution in the sequence and every step below the clustering call is written
  against it. The four handoff outputs are byte-identical to before, and a test
  pins the ordering, because reordering the list would silently re-point UMAP,
  the markers, the annotation and the handoff without raising.

- **`add_module_score` is now verified against Seurat as an equality, not a
  correlation.** Every prior comparison was correlation-based and structurally
  had to be: `AddModuleScore` picks control genes with `sample()`, R's RNG is
  not NumPy's, so the two cannot produce the same number. Pearson 0.9995 is a
  real result, but it is the RNG-limited one — it cannot separate a faithful
  port from one whose binning is subtly wrong, because binning error and
  sampling noise land in the same residual.

  `nbin=1` puts every gene in one bin and `ctrl` = pool size draws all of it;
  `sample(n, n)` is a permutation, so the control **set** is forced and only its
  summation order is free. In that regime the two tools agree to
  **6.66e-15 across 20,729 cells** — floating-point associativity over 18,649
  genes, against 1.8e-01 between two R seeds at `ctrl = 8`. The binning, control
  selection and mean subtraction are exact.

  Also adds the first coverage of the multi-program path and of any settings
  other than the defaults: three programs in one call at `nbin=12, ctrl=40`
  (Pearson 0.9972 / 0.9988 / 0.9986), with a guard against column
  transposition — programs are identified by position alone, so a swap would
  leave every correlation high and every label wrong.

  Six mutants, all caught. Two survived the first pass, both because the new
  tests fabricated the R side from Python's own column and so could not detect a
  change in how Python computed it — the same defect shape this repo has hit
  before. Closed by asserting the pipeline's settings produce a column that
  differs from a default-settings recomputation.

### Fixed

- **The R verify scripts now serialize anchors at full float64 precision.**
  Four scripts wrote their JSON references at `digits = NA`, `12` or `15`, each
  of which truncates. Measured against R Seurat 5.5.1 on a 2,000-value spread of
  realistic magnitudes: `NA` and `12` lose about half the values, `15` loses
  ~12%, and **`17` and `22` are exact**. All four now use `22`, matching
  `visium_verify.R`, which had already worked this out.

  Two things about jsonlite's `digits` are worth stating because both are
  counter-intuitive and both were assumed wrong somewhere in this repo:
  it counts **decimal places, not significant digits** — so whether a setting is
  lossy depends on the *magnitude* of the value — and **`NA` is not "max
  precision"**; it round-trips fewer doubles than an explicit `17`.

  What actually moved: the out-of-core reference gained precision on **203 of
  309** values (`1.63587331771851` → `1.6358733177185059`), Xenium spatial on
  **37 of 60** (`24.073884121202` → `24.073884121202223`), and Xenium SVF on 4
  of 81 Moran's I values. The object-model reference moved on **0 of 160**,
  because it already rounds every float to 6 dp on both sides — `NA` was losing
  nothing there, and is changed only so the next unrounded anchor does not
  inherit a trap.

  **No comparison outcome changed**: 91/91 object-model anchors, all Xenium
  deterministic anchors, and every out-of-core and SVF anchor still match, now
  against R's true values rather than truncated ones.

- **The out-of-core report no longer cuts numbers mid-digit.** Once the R
  reference carried full precision, a hard 22-character slice rendered
  `1.0153110547861388e-06` as `1.0153110547861388e-0`, dropping the exponent's
  last digit — a reader would take 1e-06 for 1e-0. Values are now shown whole or
  elided with an ellipsis, and a 1-element list from the R side is unwrapped so
  the number gets the column instead of the brackets.

- **Every tutorial CSV read now round-trips float64.** pandas' default CSV
  *reader* is not correctly rounded — it misparses about a third of random
  doubles by an ULP. `to_csv` was never at fault; it already writes the shortest
  round-trippable form. So a tutorial that wrote a value, handed it to R, read
  both back and reported "these agree exactly" was partly measuring the parser.

  37 call sites across 16 tutorials now pass `float_precision="round_trip"`.
  Three reads are exempt with stated reasons (cell labels, a row count, and the
  DE hex reader, which parses via `float.fromhex`).

  Six reported figures moved, all at the ULP level and most of them downward —
  PBMC 8k's `percent.mt` 5.773e-15 → 5.329e-15, PBMC 3k's VST mean relative
  difference 1.548e-14 → 4.973e-15, SCTransform's `detection_rate` 5.6e-16 →
  5.0e-16. **No declared band moved and all nine `--report` runs still exit
  zero.** The affected vignettes carry the measured values.

  `tests/test_tutorial_csv_precision.py` pins the convention so the next
  `pd.read_csv` cannot quietly reintroduce it, and checks that the three
  exemptions still match something rather than silently widening the lint.

  Not fixed here, because it is not lintable from Python: R's `write.csv`
  renders 15 significant digits and raising it does not help — R's own
  `sprintf("%.17g")` is not correctly rounded either. Where bit-identity is the
  actual question the R script writes a C99 hex-float side table, as
  `pbmc3k_de_verify.R` does.

- **Corrected a marker-agreement bound in the PBMC 3k docs.** They stated
  `avg_log2FC` agreement "to 4.9e-15" for both clusters whose cells match
  exactly. That covers the 151-gene cluster (4.88e-15); the 242-gene one is
  **4.62e-14**, an order of magnitude larger. Both figures are now given. This
  predates the reader fix — the bound was simply the smaller of the two.

- **`aggregate_expression(return_object=True)` left raw sums in the `data`
  layer.** Seurat's `AggregateExpression(return.seurat = TRUE)` runs
  `NormalizeData` over the pseudobulk, so its `data` holds
  `log1p(sums / colSums × 10000)`. Truecell copied the sums across unchanged,
  which every downstream function reading that layer would have taken for
  normalized expression — library-size confounded and off by orders of
  magnitude (14 where Seurat has 6.98).

  `normalization_method` and `scale_factor` now match Seurat's arguments, with
  `normalization_method=None` to opt out.

  This is the one place the two group-summary functions genuinely diverge:
  `average_expression` writes a plain `log1p` of the means with no library-size
  step, and that was already right. Both verified against Seurat 5.5.1.

  Found by pinning `aggregate_expression` against R for the first time. Its nine
  existing tests all compared Python to a Python re-derivation of the same
  formula — the shape of coverage that cannot catch a convention mismatch, and
  the reason the CLR defect once survived its own unit test.

- **`_get_expression_matrix` returned the wrong layer, and mislabelled the right
  one.** Two defects in one function, both reachable from
  `find_markers(layer=...)`, `aggregate_expression` and `average_expression`.

  The Assay5 layer dict is keyed `scale.data`; the Python argument is
  `scale_data`. Only the dotted spelling was matched, so the underscore form
  missed the dict and **fell through to the `data` fallback** — a matrix of the
  right shape, no warning, and normalized values where scaled ones were asked
  for.

  The dotted spelling then hit the second defect: the matrix came back paired
  with the assay's *full* feature list. `scale.data` holds only the scaled
  subset (`scale_data()` defaults to the variable features, as R's `ScaleData`
  does), so a 10-row matrix was labelled with 30 names and every row read as a
  different gene. This is the defect fixed in `reduction.py` under #66; the same
  one survived here. Each layer is now labelled with its own features.

  `reduction.py`, and therefore PCA, was never affected — it has its own
  accessor.

### Added

- **`average_expression`, mirroring Seurat's `AverageExpression`.** Seurat ships
  two group-summary functions and Truecell had only `aggregate_expression`.
  They are not two scalings of one thing: `AggregateExpression` sums raw counts,
  while `AverageExpression` averages the **back-transformed** `data` layer —
  `mean(expm1(x))`, not `mean(x)` and not `expm1(mean(x))`. On a small Poisson
  object the first gene reads 332.84 under one and 3.17 under the other.

  The back-transform applies to the `data` layer alone; `counts` and
  `scale.data` are not log-normalized and are averaged as they stand. All three
  layers, `features=`, multi-column `group_by` and `return_object` were pinned
  against Seurat 5.5.1 — worst absolute difference 6.8e-13 on values of order
  300. `return_object` leaves the averages in `counts` and writes `log1p` of
  them to `data`, as Seurat's `return.seurat = TRUE` does, rather than
  re-normalizing a matrix that is already an average.

  The `expm1` transform is shared with `find_markers`' fold-change path rather
  than copied, since Seurat's fold change is the same back-transform and a
  divergence between the two would be undetectable downstream.

- **`find_clusters` accepts several resolutions, and writes the per-resolution
  column Seurat writes.** `FindClusters(obj, resolution = c(0.4, 0.8, 1.2))` is
  the standard idiom for choosing a resolution — run a few, compare, then pick —
  and it had no Truecell equivalent. `find_clusters(obj, resolution=[0.4, 0.8,
  1.2])` now runs each in turn.

  The larger fix is underneath it and applies at the default settings: Truecell
  wrote **only** `seurat_clusters`, never the `{graph}_res.{resolution}` column,
  so a ported script reading `obj[["RNA_snn_res.0.5"]]` raised `KeyError`. That
  column is now written for a single resolution too.

  Three conventions were pinned against Seurat 5.5.1 rather than assumed, and
  one of them is a trap. The column label is R's number formatting, not
  Python's: `str(1.0)` is `"1.0"` where R's `as.character(1.0)` is `"1"`, so a
  naive implementation names the column `RNA_snn_res.1.0` and the ported script
  still fails. `_res_label` reproduces R's rule — render fixed and scientific,
  keep the shorter, ties to fixed — and is checked against R on 21 values from
  1e-7 to 1234567. The object is left on the **last** resolution given, not the
  largest (`resolution=[1.2, 0.8, 0.4]` ends on 0.4). And each resolution is
  clustered from the same seed rather than a running stream, so a partition does
  not depend on what preceded it or on the order asked for.

  `cluster_name` is supported, matching Seurat's `cluster.name`.

- **`split_by` on `dim_plot`, `feature_plot` and `vln_plot`.** Seurat uses three
  *different* mechanisms for this, which is the part worth getting right; each
  was probed against Seurat 5.5.1 rather than assumed.

  `DimPlot` uses `facet_wrap`: one panel per level holding only that level's
  cells, with `scales = "fixed"`. Both the axis limits and the group colours are
  therefore computed across all cells and shared, because panels are only
  comparable if a position and a colour mean the same thing in each.

  `FeaturePlot` builds a features x levels grid. The colour scale is computed per
  feature over all cells and shared along the row — per-panel scales would make
  two very different levels look alike, each filling its own range.

  `VlnPlot` does not facet at all. `split.plot = FALSE` is its default, and it
  says so itself: "Separate violin plots are now plotted side-by-side". So the
  levels are **dodged** within each group's x position and coloured by level,
  with the group carried by position. The offsets match ggplot's dodge geometry.

  Every existing figure is byte-identical; the feature is purely additive.

### Fixed

- **`vln_plot` drew the wrong violin, three ways.** All three were checked
  against R 4.6.1 / Seurat 5.5.1 running locally, not against recollection of
  what Seurat does.

  **The bandwidth was 2.3x too wide.** The density used scipy's `"scott"` rule,
  which scales the sample standard deviation. R's `stats::density` — and so
  `geom_violin`, and so Seurat — uses `bw.nrd0`, which takes
  `min(sd, IQR/1.34)`. Expression is zero-inflated, so the IQR term is much the
  smaller, and measured against R on zero-inflated draws scipy's bandwidth came
  out 2.0-2.5x wider. The effect is over-smoothing that flattens the spike at
  zero, which is the shape of the distribution. `_bw_nrd0` now reproduces R to
  5e-13 across seven cases including both zero-spread fallback branches.

  (`bw.nrd0` divides the IQR by **1.34**. The neighbouring rule `bw.nrd`, which
  R also ships, uses 1.349 — taking that one is a silent 0.67% error wherever
  the IQR term wins. The first version of this fix had exactly that bug, caught
  by diffing against R.)

  **The violin was not trimmed.** `geom_violin(trim = TRUE)` limits the density
  to the observed range; an untrimmed gaussian KDE tails off past it, so every
  violin extended below zero where expression cannot go.

  **Points were off by default.** `pt_size` defaulted to `0`. Seurat's `VlnPlot`
  passes `pt.size = NULL`, which `ExIPlot` resolves through `AutoPointSize` —
  `min(1583/n, 1)` — so points are shown and shrink as the cell count grows.
  `pt_size=None` is now the default and follows that rule; `pt_size=0` still
  omits them.

  The violin outline now matches `geom_violin(scale = "width", trim = TRUE)`'s
  own computed polygon to 0.2% of full width, with the support equal to the data
  range exactly. Nine tutorial figures are regenerated; they are the nine drawn
  by the six generators that call `vln_plot`, and no others moved.

  Two additions fall out of the rewrite: `violin_width` sets the width of a full
  violin, and `jitter_seed` (default `0`) makes the point jitter reproducible so
  a committed figure redraws identically.

  Still deliberately different from Seurat: the median bar. `geom_violin` draws
  none, and this keeps drawing one.

### Added

- **`ridge_plot` warns when an explicit `figsize` is too small for the group
  labels**, instead of silently returning a figure with one panel's labels
  printed across its neighbour. Ridgeline labels sit outside the axes, so a
  panel is wider than its axes by the length of the longest label; when the
  canvas cannot fit that, matplotlib overlaps the panels rather than erroring.

  The warning names the computed default, which a 162-configuration sweep — up
  to 12 groups, labels as long as `"Haematopoietic stem cell"`, every `ncol` —
  found collision-free, and a test asserts the suggested size actually resolves
  the collision. Recommending it beats scaling the caller's size by some
  function of the overlap, which can still collide.

  No layout engine fixes the underlying case: `constrained_layout` overlaps by
  slightly *more* here and reports "axes sizes collapsed to zero", because the
  space genuinely does not exist. The default path is unaffected and pays no
  cost — the check only runs when `figsize` was passed.

### Fixed

- **The categorical palette repeated colours, so two clusters could render
  identically.** `_PALETTE_36` was named for 36 entries but held only 30
  distinct ones — `#A3A500` sat at both index 14 and index 29, and five more
  colours were doubled in the last five slots. `_palette` sliced that list for
  any `n <= 36` and only fell through to a `tab20` ramp above it, so an object
  with 30 to 36 groups drew two of them in exactly the same colour, silently.
  That is a wrong plot, not an aesthetic one: nothing on the figure said which
  cluster was which.

  The list is gone. `hue_pal(n)` now computes the ramp the way ggplot does —
  `n` evenly spaced HCL hues at c=100, l=65, the polarLUV-to-sRGB conversion
  written out so R's per-channel `fixup = TRUE` clamp is reproduced rather than
  approximated. It cannot run out and it cannot repeat.

  It is also *more* faithful, not just safer. ggplot spreads the hue circle
  across however many groups there are, so the colours for 9 groups are not the
  colours for 8 plus one more. The old list had been built by concatenating
  several such runs, which is where the duplicates came from — and it therefore
  matched Seurat only at n=8. `hue_pal` matches at every n, verified exact
  against R for n = 1-6, 8 and 9.

  **This changes existing figures.** Any plot coloured by group where the group
  count is not 8 gets different (correct) colours. Numbers, ordering and layout
  are untouched.

- **Only `feature_plot` rasterised its points, and it did so unconditionally.**
  Every other scatter drew one vector path per cell, so a PDF or SVG of a
  100k-cell `dim_plot` embedded 100k circles. `dim_plot`, `feature_scatter`,
  `variable_feature_plot`, `image_dim_plot`, `image_feature_plot`,
  `spatial_dim_plot` and `spatial_feature_plot` now take `raster=` alongside
  `feature_plot`, defaulting to `None` — which resolves to Seurat's own rule,
  rasterise above 100,000 points. PNG output is unaffected either way.

### Added

- **A theme layer: `set_theme`, `theme_context`, `get_theme`, `reset_theme`.**
  The module named absolute point sizes at 40 call sites (`fontsize=8` twelve
  times, `fontsize=9` ten times, and so on), so changing the house style meant
  editing the source. Text now scales from one base size through named roles,
  and `set_theme(base_size=13, style="seurat")` moves all of it at once.

  `base_size` also writes `rcParams["font.size"]`, because the roles only cover
  text this module sizes explicitly — axis and tick labels are matplotlib's, and
  without that a larger base grew the titles and left the axis furniture at
  10pt.

  Two style presets ship: `"seurat"` (cowplot's look, which is what Seurat
  draws) and `"minimal"`. Neither is applied unless asked for; the default
  touches no rcParams at all.

  The default `base_size=10` reproduces the previous absolute sizes exactly —
  0.8 x 10 == 8, and so on. Verified by rendering twelve plots on both sides of
  the change with the old palette pinned: all twelve came out byte-identical,
  which is also what isolates the palette above as the only behavioural change.

- `hue_pal(n)` is public — the ggplot/Seurat discrete colour scale, useful for
  matching group colours in a figure drawn outside truecell.

### Changed

- **CI installs from `uv.lock`, and so should you.** The test matrix ran
  `uv pip install -e ".[all]"`, which ignores the lock and resolves fresh
  against the `>=` floors in `pyproject.toml` — so CI, and every contributor,
  could be on a different scientific stack from the one the committed tutorial
  figures were drawn with.

  That is not hypothetical. Rebuilding a development environment moved
  umap-learn, scikit-learn and NumPy forward, and ten CITE-seq figures shifted
  by up to 11% of their pixels while **every** number the vignette asserts
  stayed identical — Pearson 0.9847, 99.29% label concordance, 16 and 21
  clusters, the whole per-cell-type weight table. Clustering runs off the SNN
  graph; the UMAP embedding is only ever displayed. So the figures drift and no
  anchor reports it.

  **`uv.lock` was itself gitignored**, which is the same bug one level down: a
  lock only its author has is not a lock. It is committed now. The usual advice
  to leave a library's lock out of the repository is about the *published*
  package and still holds — what a user installs is set by the `>=` floors in
  `pyproject.toml`, which this does not touch. It constrains the development and
  CI environment, and that environment has an output: the committed figures.

  `uv sync --all-extras --locked` now installs both matrix legs, and `--locked`
  fails rather than re-resolving if `uv.lock` has fallen behind
  `pyproject.toml`. `docs/installation.md` and the `truecell-dev` skill point at
  the same command, so figures are regenerated from the versions CI tested.

- **A weekly job resolves dependencies fresh, and fails when they break us.**
  Pinning CI to the lock gives up what the unpinned install did for free:
  noticing that a new numpy or pandas has broken the package a user installing
  today would get. `freshdeps` is that check, moved rather than dropped — it
  ignores the lock, runs the suite, and **fails loudly** instead of reporting
  advisorily, because a `|| true` here would be a green tick on a check that had
  stopped checking. It runs on a schedule and on demand, not per pull request,
  so upstream breaking is never mistaken for the pull request breaking.

  It resolves under both 3.12 and 3.13, matching the `test` matrix — an upstream
  release can break one interpreter and not the other, and a canary that only
  resolved under 3.12 would stay green through it. The legs do not fail fast,
  because which legs went red is the diagnosis: one is an interpreter problem,
  both is the package.

  The distribution job deliberately still installs the wheel *unlocked*: it
  exists to prove a plain `pip install truecell` works for a real user, and
  pinning it would defeat the point.


## [1.0.0] - 2026-07-29

The rename release. `shanuz` became `truecell`, and the version moves to 1.0.0
rather than continuing the 0.x line: `shanuz` 0.9.0 is the last release under
the old name and `truecell` 0.9.0 was published from the same code before the
rename settled, so neither number was available to carry this. It is the same
codebase, renamed and re-verified, not a maturity claim -- the API is the one
0.9.0 shipped.

### Changed

- **The package is renamed from `shanuz` to `truecell`.** This is a breaking
  change and there is no compatibility shim: `import shanuz` will not work, and
  the top-level class is `Truecell`, not `Shanuz`. Every import, every skill
  directory, the documentation site and the logo assets move with it.

  `shanuz` remains on PyPI at 0.9.0 and is not being withdrawn, but it will not
  receive further releases; new versions are published as `truecell`. The
  GitHub repository was renamed in place rather than forked, so stars, forks,
  issues and pull requests carry over and the old URL redirects — including for
  `git clone` and `pip install git+…`. Clones with an existing `origin` keep
  working through that redirect, though `git remote set-url` is worth running.

  The logo changed with the name. The mark was a lowercase `s`; it is now a
  lowercase `c`, and the wordmark was reset in the new name — which meant
  drawing five letterforms (`t`, `r`, `e`, `c`, `l`) that the six letters of
  `shanuz` had never needed. The two-lobe colour split moved to the `c`, the
  only letter in `truecell` that is a single unbroken arc and so the only one
  whose midpoint is a non-arbitrary place for the colour to change.

### Fixed

- **The sdist shipped the tutorial data, and shipped a different amount of it
  every time.** `docs/tutorials` is a symlink to `tutorials/`, hatchling
  follows it, and the tutorials write their intermediates into the directory
  they live in — the BPCells stores, the `.lazy` matrices, the R handoff CSVs.
  `.gitignore` covers all of those, but an sdist is not built from git, so they
  were being packaged. The published `shanuz` 0.9.0 sdist is 23.5 MB against a
  0.24 MB wheel for that reason, and a build on a machine that had just run the
  tutorials produced 63 MB: the same version number, a different tarball,
  depending on the disk it was built from.

  `[tool.hatch.build.targets.sdist]` now lists its contents explicitly — the
  package, the tests, and the top-level metadata. The sdist is 532 KB. The
  wheel is unaffected; it was always scoped by `packages = ["truecell"]`.

  `tests/test_packaging.py` builds an sdist and looks inside it rather than
  asserting on the configuration, which would pass just as happily against a
  broken artifact.

### Added

- **A logo** (`docs/assets/logo/`), generated by `tools/make_logo.py` rather
  than drawn — the mark is a point cloud whose density traces a lowercase `c`,
  for `cell`. Building a picture out of dots is the nod to Georges Seurat, the
  pointillist the R package is named for. The `c` is two arcs meeting at the
  waist in two colours: the two implementations, one shape. It is also the
  letter the wordmark's own `c` draws, so the mark is a letter of the name
  enlarged rather than an ornament beside it. Ships as a mark, a horizontal lockup, a wordmark,
  a hex sticker, a simplified single-stroke glyph for sizes below about 48px
  where the dots silt up, and favicons; each in a variant for a light ground and
  one for a dark. Colours are the documentation site's, unchanged, so the two
  are one system. The site now carries it as its header logo and favicon, and
  the README as a `<picture>` that follows GitHub's light and dark modes.

- **Agent skills for LLM-assisted work with the package** (`skills/`). Ten
  skills in the [Claude Agent Skills](https://docs.claude.com/en/docs/agents-and-tools/agent-skills/overview)
  format — plain Markdown with YAML frontmatter, so they are equally usable as
  context for any other model. A router (`truecell`) carrying the API contracts
  that break code silently, bundled with a full signature-and-Seurat-equivalent
  API map and an object-model reference; then one skill per area: `workflow`,
  `differential-expression`, `integration`, `multimodal`, `spatial`, `at-scale`,
  `plotting`, `from-seurat`, and `dev` for working on the package itself.

  The content is aimed at the mistakes a model actually makes against this API
  rather than at re-describing the docs: that analysis functions mutate in place
  and return `None` (so `obj = normalize_data(obj)` leaves you holding `None`),
  that `dims` is 0-based, that matrices are features × cells, that the generics
  are not top-level, and — in each domain skill — which differences from Seurat
  are deliberate and must not be "fixed". Every signature, default and measured
  number in them was derived from the installed package or from a recorded R
  comparison; a validation pass resolves each `truecell.*` and `generics.*` name
  in all thirteen files against the live package.

  `.claude/skills` is a symlink to `skills/`, so Claude Code discovers them when
  run in this repo — the same idiom `docs/tutorials` already uses, and for the
  same reason: one copy, not a second one that drifts.

### Fixed

- **The API reference claimed the whole API is exported from the top level.**
  `docs/api/index.md` opened with "Everything on these pages is exported from
  the top level", which is true for eleven of the thirteen pages and false for
  two. 66 of the 73 public generics live only on `truecell.generics`, so
  `truecell.features(obj)` — the call the sentence invites — raises
  `AttributeError`; the loaders on the Loading data page are likewise
  `truecell.io.read_10x`, `truecell.datasets.pbmc3k`,
  `truecell.compat.anndata.as_anndata`. Rewritten to say *most*, name both
  exceptions, and list the seven generics that genuinely are re-exported.
  `tests/test_docs.py` now pins the counts, since prose is the one thing
  `mkdocs build --strict` cannot check.

- **`tutorials/README.md` documented only 10 of the 18 tutorials.** The
  overview table at the top listed all 18, but the detailed walkthrough
  sections stopped after Tutorial 10 and jumped straight to the API quick
  reference — so the eight tutorials added in Waves 2 and 3 (the object
  model, dimensional-reduction extras, sketching, the spatial container and
  Moran's I, the DE test suite, out-of-core `LazyMatrix`, Visium, and anchor
  internals) had no written section, on GitHub or on the
  [docs site](https://genomicai.github.io/truecell/tutorials/). Added all
  eight, in the same format as the rest.

- **Most figures on the docs site 404'd.** MkDocs re-anchors relative paths
  written in Markdown syntax onto the built page, which under
  `use_directory_urls` sits a directory deeper than its source — but it passes
  raw HTML through untouched, and the vignettes write most of their figures as
  `<img>` inside HTML tables so the R and the Python plot sit side by side.
  Those resolved one directory too deep: 110 of the site's 133 figures were
  broken, and the ten vignettes built entirely out of those tables showed no
  images at all. Added a build hook (`tools/mkdocs_html_relpaths.py`) that
  applies the same rewrite to raw `<img src>`, so the one path in the source
  stays correct both on the site and on GitHub, where the vignettes are also
  read. A reference with no file behind it is now a warning, which
  `mkdocs build --strict` turns into a failed CI build.

  The existing figure-existence test only understood Markdown image syntax, so
  it was blind to 126 of the 149 references it was meant to be guarding; it now
  reads both syntaxes.

## [0.9.0] - 2026-07-27

Work from six milestones — reference mapping, extra reductions, pseudobulk DE,
spatial, scale, and the specialized assays — plus one breaking fix, plus the
tutorial fidelity infrastructure that followed (measured bands, staleness
guards on the R-comparison reports). All of it was on `main`; as of this
release, all of it is on PyPI.

### Added

- **A documentation site.** MkDocs + Material + mkdocstrings, built from `main`
  and published to GitHub Pages, at
  [genomicai.github.io/truecell](https://genomicai.github.io/truecell/). The API
  reference is generated from the docstrings — which is the point: eighteen
  tutorials' worth of hard-won fidelity notes lived in them and rendered
  nowhere. Alongside it, four hand-written pages (overview, installation, a
  quickstart whose every line was executed to produce the output it shows, and
  **Fidelity**, which collects how the port is checked against R, what the
  checking has caught, and the differences that are real), plus all eighteen
  vignettes and their figures. `docs/tutorials` is a symlink to `tutorials/`, so
  there is one copy of each vignette rather than a docs copy that drifts.

  **The blocker recorded for this item was not real.** `ROADMAP.md` said
  mkdocstrings resolves annotations and that `typing.get_type_hints()` raising
  `NameError` on the plotting and `graph`/`neighbor` signatures would stop the
  site on day one. It still raises, on 19 callables — but mkdocstrings reads
  annotations statically through griffe, without importing the module, so those
  signatures render and cross-link *because* the `if TYPE_CHECKING:` imports are
  there. Measured rather than assumed, and written up on the Fidelity page.

  Building it surfaced defects the site was the first consumer to notice:

  - **Three `Returns` sections were being parsed as parameters.** `sctransform`,
    `add_module_score` and `composition_test` each ended their parameter list
    with a bare `Returns ...` sentence at the same indentation, so a NumPy
    parser read it as another parameter and the functions documented no return
    value at all. `composition_test` also had its eleven output column names
    parsed as eleven parameters.
  - **Every one of the 370 documented parameters put its description in the
    type slot.** The convention across the package is `name : what it does` —
    prose, never a type, because the types are in the signatures. NumPy format
    reads that slot as the type, so the site printed "min features a cell must
    have to be kept" as `min_features`'s type and left its description empty.
    Fixed in `tools/griffe_sphinx_roles.py`, a build-time translation, rather
    than by rewriting 370 lines of source that are not wrong.
  - **129 Sphinx cross-reference roles rendered as literal `:func:` text.** The
    same extension turns them into real links, but **only** when the target is
    one the site actually publishes; anything else degrades to a code span,
    because `--strict` treats an unresolvable link as a failure.
  - **Fifteen `Slots` blocks collapsed into one run-on line.** Markdown joins
    consecutive lines, so an aligned column block became a paragraph. Rewritten
    as bullet lists in the source, which reads better under `help()` too.
  - **Two heading anchors in `xenium_spatial_tutorial.md` were dead** under
    python-markdown's slugger while being correct on GitHub. Fixed by adopting
    GitHub's slugger site-wide, so one set of links works in both places.

  `tests/test_docs.py` asserts the parts that rot silently: every public export
  appears on an API page, every `:::` directive names something that still
  exists, no symbol is rendered twice, every nav entry and every vignette figure
  resolves, and the whole site builds under `--strict`. All nine mutations of
  those guards were killed.

- **The two R-comparison numbers that were allowed to move are now asserted
  bands, and the reports refuse a stale R reference.** `deseq2`'s overlap with
  Seurat's top 50 and JackStraw's PC cutoff both differ from R for understood
  reasons — pseudobulk against per-cell in one case, a seeded permutation
  against R's loop-indexed one in the other — and both were recorded as prose.
  Prose does not fail, so a genuine regression landing inside the expected
  spread read as ordinary variation. `tutorials/bands.py` gives each number a
  `Band` with a stated reason; `--report` prints the verdicts and **exits
  non-zero** outside them. The bands are measured, not guessed: a 60-seed sweep
  of `jack_straw` (keeps of 12/13/14/15 for 2/28/11/19 seeds against R's
  deterministic 13) and 20 resampled pseudo-replicate splits for `deseq2`
  (20–26, median 22). The seven cell-level DE tests get exact bands, so one
  dropped gene out of the top 50 is now a failure.

  Three defects surfaced while wiring it up, all the same shape — an R
  reference silently older than the handoff it answers:

  - **`pbmc3k_de_tutorial.py --report` compared across runs without noticing.**
    On the working copy the Python tables were from 25 July and the R tables
    from 19 July, taken on a *different* cluster assignment, and the report
    printed a full parity table with `wilcox` at 48/50 and a p-value Spearman of
    0.907 — indistinguishable from a regression in the port. Re-run against a
    matching reference, `wilcox` is 50/50 at Spearman 1.000000. `pct.1`/`pct.2`
    are counts over the shared cells with no statistics in the way, so they now
    gate the comparison: they had differed for 12,491 of 13,712 genes.
  - **`pbmc3k_dimreduc_tutorial.py` NaN-filled instead of failing.** `_read`
    reindexes R's embedding onto the Python cell order, which fills an absent
    barcode with `NaN` rather than raising; every correlation downstream would
    have printed as `nan`. The feature check is now set equality rather than
    one-sided containment — an R run holding *extra* features passed before.
  - **Running `pytest` corrupted the R handoff.** `prep` wrote
    `figures_dimreduc/hvg_features.txt` and `cells.txt` unconditionally, so the
    test suite's synthetic fixtures replaced the 2,000 real HVGs with 100 genes
    named `GENE171` and the 2,700 barcodes with `CELL0..CELL119` — which
    `pbmc3k_dimreduc_verify.R` reads. `prep` now takes `out_dir`.

  Both vignettes are refreshed against re-run R references. The dim-reduction
  basis check improves from `median |r| 0.9759, matched through PC 15` to
  `1.0000, matched through all 20`: the old figure was not PCA disagreement at
  all, but 6 of the 2,000 HVGs having drifted across the selection boundary
  since the reference was built. All eight mutation trials against the new
  guards were killed.

- **The PBMC 3k and PBMC 8k tutorials now have numeric handoffs against R
  Seurat, closing the last two compared entirely by eye.** Every R panel in the
  PBMC 3k walkthrough links a canonical satijalab.org image, so nothing failed
  if the numbers behind those pictures drifted. `pbmc3k_verify.R` now writes
  per-cell QC and PCA keyed by barcode, per-gene VST statistics, the marker
  table and a set of scalar anchors; `pbmc8k_subclustering_verify.R` writes the
  global cell table and the T/NK compartment's membership. Both tutorials gained
  a `--report` that compares them. Neither side is pinned to the other — both
  run their own pipeline from the same 10x bytes, which is what makes the
  comparison mean anything.

  PBMC 3k against Seurat 5.5.1: **the same 2,638 barcodes** survive QC (nCount
  and nFeature exact, percent.mt to 5.3e-15), the VST means agree to 4.8e-14,
  **1,998 of 2,000** variable features are shared, PCA over the 10 dims the
  clustering uses matches at **|r| 0.9988** with no reordering, the kNN graph is
  52,760 edges on both sides, and the clusters agree at **ARI 0.938**
  (8 vs 9 — the one Seurat has and truecell does not is a 32-cell dendritic-cell
  population whose cells land, all 32, in truecell's CD14+ Mono cluster).

  PBMC 8k: **the same 7,475 barcodes**, global clusters at **ARI 0.977**, and
  the T/NK compartment handed to the subclustering stage matches at **Jaccard
  0.9991** — 4,631 of 4,635 cells are the same barcodes. Subclusters then agree
  at ARI 0.916 and the subset labels at 98.2%. The compartment is compared by
  barcode on purpose: everything in stage two is conditioned on which cells
  stage one selected, so a compartment of the right size drawn from the wrong
  clusters would pass a count check and make every later number incomparable.

  Two stale claims fell out. `tutorials/README.md` recorded "Clusters at
  resolution 0.5: 9 ✅" and "All 6 canonical cell types recovered … DC ✅" for
  PBMC 3k; truecell finds 8 and does not separate DC on that dataset. Both, and
  the vignette's "**9 clusters** in both R and Python", are now the measured
  numbers with the provenance traced from the LOESS fit down to the 286 SNN
  edges that move the boundary. PBMC 8k runs the *opposite* way — there truecell
  splits Seurat's merged Platelet/DC cluster in two — so neither run is
  uniformly finer, which the docs now say.

- **The CITE-seq tutorial now compares the CLR transform and the WNN weights,
  per protein and per cell.** `cbmc_citeseq_verify.R` dumps the per-protein CLR
  summary and the per-cell modality weights **keyed by barcode**; `--report`
  compares them. Keying by barcode is the point: a mean-per-cell-type table asks
  two different cluster partitions the same question and cannot separate "the
  weights differ" from "the labels differ".

  Measured against Seurat 5.5.1: ADT CLR agrees to **4.2e-15** on every
  statistic, per-cell ADT weight correlates at Pearson **0.9847** (Spearman
  0.9816, median abs diff 0.0152) over 8,617 shared barcodes, cell-type labels
  are **99.29%** concordant, and cells, genes, proteins, RNA clusters and WNN
  clusters all match exactly.

  **This settles the progenitor question the vignette had left open.** It
  recorded progenitor 0.06 away from Seurat where every other cell type agreed
  to 0.02, and said explicitly that the small-population-noise explanation was
  "a reading, not something this tutorial establishes". Re-grouping truecell's own
  weights by R's labels puts progenitor at **0.283** against R's 0.285 — the
  weights agree; 61 of 8,617 cells sit on the progenitor/erythroid boundary and
  land differently. Not a WNN difference.
- **The SCTransform tutorial now compares the fitted model, not just the
  figures.** `pbmc3k_sctransform_verify.R` dumps Seurat's own
  `SCTModel.list[[1]]@feature.attributes` plus the ranked variable features and
  the clips; the Python side writes the same table and
  `--report` puts them side by side, per gene. This was the highest-risk gap in
  the suite: tutorials 1-4 were compared by eye, and SCTransform is where a
  four-part defect once hid behind a perfectly plausible UMAP.
- `sctransform` now stores the **whole fitted model** on the SCT assay's
  `meta_data`, under Seurat's column names: `residual_mean`, `detection_rate`,
  `(Intercept)` and `log_umi` join the existing `residual_variance`, `theta` and
  `gmean`. The regularized intercept and slope — the model itself — previously
  could not be inspected at all.

  Measured over all 12,572 modelled genes against Seurat 5.5.1: `detection_rate`
  and `gmean` agree to 5.6e-16 and 1.2e-12, `(Intercept)` and `theta` both rank
  identically (Spearman 1.0000), the 3,848 genes v2 declares non-overdispersed
  are exactly the same set (Jaccard 1.0000), and residual variance ranks at
  0.9986. `residual_mean` is the one column that does not track by rank
  (Spearman 0.71, Pearson 0.99) and the one nothing downstream reads.

### Changed

- **BREAKING: `find_variable_features` now writes the column names
  `HVFInfo()` returns.** On an Assay5 the per-feature statistics land in
  `assay.meta_data`, which — unlike SeuratObject, where they sit behind a
  `vf_vst_counts_` prefix that `HVFInfo()` strips — *is* the user-facing table.
  It spelled them `means` / `variances` / `variances.standardized` against
  Seurat's `mean` / `variance` / `variance.standardized`, so the same quantity
  read differently in the two languages and the PBMC 3k handoff had to rename
  all three on the way out for the join to work. Now:

  | method | columns |
  |---|---|
  | `"vst"` | `mean`, `variance`, `variance.expected`, `variance.standardized` |
  | `"mvp"` / `"dispersion"` / `"mean.var.plot"` | `mvp.mean`, `mvp.dispersion`, `mvp.dispersion.scaled` |

  Code reading `meta_data["means"]` and friends must be updated; the old names
  are gone rather than aliased, since keeping both is what let the divergence
  live in the object unnoticed. `highly_variable` is unchanged — it is truecell's
  own flag, and `HVFInfo(status = TRUE)`'s `variable` has no other counterpart.

  Two things fall out of the rename. **`variance.expected` is new**: the vst
  LOESS fit was computed, used to standardize, and then discarded, so nothing
  downstream could tell a genuinely variable gene from one sitting where the fit
  was poor. On PBMC 3k it agrees with Seurat to **2.50e-2 relative** against
  `variance.standardized`'s **2.56e-2** — the two carry the same disagreement,
  as they must, since the standardization is proportional to its reciprocal.
  That is the whole of the VST gap, now attributable: the mean and the observed
  variance match to 1.5e-14 and 5.1e-14, so what the two tools disagree about is
  the LOESS fit and nothing else. And **the mvp path no
  longer writes vst's names**: it had been storing scaled dispersions in a
  column called `variances.standardized` and the raw variance in one called
  `variances`, neither of which is a quantity `HVFInfo(method = "mvp")` reports.
  Those values differ from Seurat's in definition as well as in name — truecell
  computes them on log-normalized data where Seurat round-trips through
  `expm1` — so this fixes the labels, not the numbers.

  `variable_feature_plot` reads the stored statistics and silently falls back to
  recomputing a raw dispersion when it cannot find them, so a consumer left on
  the old spelling would have kept drawing a plot — a different one, on a
  different axis, with no error. `tests/test_hvf_column_names.py` pins the
  column names against a recorded Seurat 5.5.1 probe, pins `variance.expected`'s
  *contents* through the identity `variance.standardized = variance /
  variance.expected` (exact wherever the clip does not bite), and asserts the
  plot takes the stored-statistics branch. The PBMC 3k handoff now carries
  `var.expected` on both sides, which is what separates "the standardization
  drifted" from "the LOESS fit did".
- **The SCTransform vignette's ±1 cluster gap is closed and its docstring
  corrected.** The vignette described truecell resolving 13 clusters against
  Seurat's 12, blamed on the RNG and the differing clustering libraries. Both
  arms now agree exactly — 12 under SCTransform, 11 under LogNormalize — which
  PR #55's graph work closed rather than anything in SCTransform. Separately,
  `truecell/sctransform.py`'s module docstring claimed 99.7% variable-feature
  agreement and theta at 0.96; the reproducible comparison measures 97.1% and
  Spearman 1.0000. The docstring now carries the numbers the tutorial prints.
- **The object-model tutorial's R reference now uses exact neighbours.**
  `pbmc3k_objects_verify.R` called `FindNeighbors` with Seurat's default
  `nn.method = "annoy"`, which is approximate, while truecell's neighbour search
  is exact — so the graph anchors compared two different neighbour tables and
  reported a difference belonging to annoy rather than to either object model.
  It cost 182 SNN edges (199,434 against 199,616). Pinning `nn.method = "rann"`
  makes both sides exact, and the tutorial now matches **91 of 91** anchors,
  up from 90. Confirmed two ways: Seurat under `rann` returns 199,616 directly,
  and feeding Seurat's own `rann` indices into `_build_snn` reproduces both the
  count and the weight sum.
- `objects_vignette.md`, `tutorials/README.md`, `ROADMAP.md` and the tutorial's
  own docstring described the symmetrised kNN graph and the dropped SNN
  self-edge as open. Both were fixed below; the text and the kNN-degree figure
  now record them as closed.

### Fixed

- **One test had been shadowed by a same-named copy and never ran.**
  `test_anchors_seurat_parity.py` defined
  `test_pca_loadings_are_exact_not_randomized` twice; Python keeps the second,
  so the earlier one — which additionally proves sklearn's *randomized* solver
  disagrees on the trailing PCs, i.e. that the exact SVD matters at all — was
  collected by nothing. Renamed and now runs (and passes).
- **`sctransform(vars_to_regress=...)` was tested against no baseline.** The
  guard asserted only that the residuals come back uncorrelated with the
  covariate, which is equally true when there was never any covariate signal to
  remove; it computed a `before` from a *different* object built from
  un-injected counts and then never asserted on it. It now runs SCTransform
  twice over the same injected counts and requires the unregressed residuals to
  carry the signal (0.25) before requiring the regressed ones not to (0.00).
- **`find_markers` named the missing ident in one error and not the other.** An
  empty `ident_2` raised a constant `"No cells found for comparison group."`
  while the `ident_1` branch named what it looked for. Both now name it, and
  the `ident_2=None` case says there is nothing left outside `ident_1` rather
  than reporting a parameter the caller never passed.
- **`ruff check` is clean, at both scopes** — 51 → 0 in `truecell`, 201 → 0 for
  the repo. Mostly unused imports left behind by the Assay5 refactor, plus dead
  locals and semicolon-joined statements. Two of the dead locals were computed
  colour maps in `vln_plot` and `dim_plot` that nothing read, and one was a
  `v3` layout flag in `read_10x` that no branch consumed. The 121 `E402`s in
  `tests/` and `tutorials/` are the deliberate `sys.path` bootstrap those files
  need to run from a clone, and are now ignored by scope in `pyproject.toml`
  rather than by 121 scattered `# noqa`. Note CI resolves **ruff unpinned**,
  which currently means 0.16.0 and a much larger default rule set (PEP 604
  annotations, import sorting); this entry is measured against 0.15.20.

- **Any plot with more than 36 groups raised `AttributeError` on matplotlib
  3.9 or newer.** `_palette` fills the first 36 colours from a fixed list and
  falls back to a colormap past that, and the fallback called
  `plt.cm.get_cmap(name, lut)` — deprecated in 3.7 and **removed in 3.9**. The
  project has supported `matplotlib>=3.7` throughout, so on any current install
  the branch was dead code that crashed the moment it was reached. Nine plotting
  functions route through it; no test had ever asked for that many groups. It
  now uses `.resampled`, the documented replacement, which reaches back to 3.6
  and so works across the whole supported range.
- **A layerless v5 assay came back as `KeyError: None`.** `Assay5.default_layer`
  is `None` when the assay holds no layers, and the two layer getters — in
  `preprocessing` and in `plotting` — indexed the layer dict with it directly.
  Both now raise `ValueError("No layers available.")`, matching what
  `Assay5.layer_data` has always said for the same condition.
- **`multiseq_demux` and `_integrate_anchor_reduction` rebound their own
  parameters**, the same idiom as the dataset loaders below; both use a local
  now. No behaviour change, but neither `multiseq_demux`'s `qrange` nor
  `integrate_layers`' single-element `group_by` list had any test coverage, so
  both are now pinned. This clears the last of the package's type-checker
  errors: **`mypy truecell` reports none**, down from 44.

- **`cbmc_citeseq` reported an unrelated error when the species filter matched
  nothing.** A `species_prefix` that no gene starts with — the wrong one, or a
  file whose row labels are not prefixed — left every chunk empty and surfaced
  as scipy's `blocks must be 2-D, and some must be sparse` from `sp.vstack`,
  which names neither the file nor the prefix that caused it. It now raises a
  `ValueError` naming both.
- **The four oldest dataset loaders rebound their `data_dir` parameter to a
  `Path`.** `pbmc3k`, `pbmc8k`, `cbmc_citeseq` and `xenium_mouse_brain` assigned
  a `Path` back over a parameter declared `Optional[str]`, which was most of the
  package's remaining type-checker noise (33 of 44 errors, now 11). They use a
  local `Path` like the newer loaders alongside them. No behaviour change — but
  nothing had ever checked that an explicit `data_dir` was honoured at all, so
  that is now covered, along with the RNA chunk stitching every real
  `cbmc_citeseq` load depends on and none of the tests reached.

- **The `mvp` statistics were three different quantities under Seurat's
  names.** PR #64 renamed the dispersion path's columns to `mvp.mean`,
  `mvp.dispersion` and `mvp.dispersion.scaled`; the numbers underneath were not
  what `HVFInfo(method = "mvp")` returns. `CalcDispersion` calls `FastExpMean`
  and `FastLogVMR`, and both **undo the log first**:

      mvp.mean       = log1p(mean(expm1(x)))
      mvp.dispersion = log(var(expm1(x)) / mean(expm1(x)))     # sample variance

  truecell took the mean and variance of the log-normalized values directly, with
  an epsilon inside each logarithm. On PBMC 3k its `mvp.mean` column ran 0–2
  where Seurat's runs 1–7. Three further divergences fell out of the same
  function: the bins were **equal-frequency percentiles of `log(mean)`** rather
  than R's 20 **equal-width** bins across the mean's range, which changes every
  z-score because the scaling is per bin; a bin holding a single gene was
  scored **0** where R's `sd` of one value is `NA`; and the within-bin standard
  deviation used `ddof=0` against R's `sd`. Now verified against Seurat 5.5.1
  on all 13,714 PBMC 3k genes — `mvp.mean` 8.9e-16, `mvp.dispersion` 2.2e-15,
  `mvp.dispersion.scaled` 5.3e-15, and the one `NaN` lands on the same gene.
- **`mean.cutoff` and `dispersion.cutoff` were accepted and discarded.** Seurat
  has *two* selectors here, not one: `MVP` (`"mvp"` / `"mean.var.plot"`) keeps
  every gene inside both cutoffs and ignores `nfeatures`, while `DISP`
  (`"dispersion"` / `"disp"`) takes the top `nfeatures` and ignores the cutoffs.
  truecell ran the second for every spelling, so `mean.var.plot` returned exactly
  `nfeatures` genes under a name that promises a cutoff — 2,000 on PBMC 3k
  against Seurat's **1,006**. Both selectors also rank by the **raw**
  dispersion; truecell ranked by the scaled one, which reorders the list. All
  four spellings now return Seurat's features in Seurat's order. `"disp"` is
  newly accepted, and `num_bin` / `binning_method` are exposed.
- **`variable_feature_plot` ignored the mvp columns.** It looked for the vst
  statistics and, not finding them, fell through to recomputing
  `E[x²] − E[x]²` off the data matrix — a third quantity, on a third scale,
  drawn without complaint under a y axis reading "Dispersion". It now plots
  `mvp.mean` against `mvp.dispersion` as `VariableFeaturePlot` does.
- **Re-running `find_variable_features` with another method left the previous
  method's columns behind.** SeuratObject can hold both because it namespaces
  them by method and layer (`vf_vst_counts_mean`) and makes you name one to read
  them back; truecell's `meta_data` has a single flat name per statistic. A vst
  run followed by an mvp run left `variance.standardized` sitting beside an mvp
  `variable_features` list it did not describe — and that is what made
  `variable_feature_plot` draw the vst figure over the mvp genes. Whichever
  method runs now owns the table.
- **`add_module_score` read the expression matrix once per gene, and was 100×
  slower than it needed to be.** Each gene's row was pulled out on its own
  (`mat[i, :]` inside a list comprehension) and the stack handed to `np.mean`.
  Every assay layer here is **CSC**, so slicing a single row walks the whole
  column-major matrix: on the THP-1 ECCITE data (18,381 × 20,729, 69.5M
  nonzeros) that is ~22 ms per gene, and the default `ctrl=100` draws a couple
  of thousand control genes. **50 of 51 profiled seconds** were inside one scipy
  call, `get_csr_submatrix`, invoked once per gene. Now one row selection per
  gene set, transposed to CSR so the rows are summed in the order they were
  asked for:

  | | before | after |
  |---|---|---|
  | `cell_cycle_scoring` (THP-1, 2 programs) | 168.4 s | **0.90 s** |
  | `add_module_score` (30-gene program) | 51.3 s | **0.51 s** |

  R's `CellCycleScoring` takes ~12 s on the same data, so truecell now runs it an
  order of magnitude faster rather than 14× slower. The spelling is load-bearing
  and **bit-identical** to the old arithmetic: an indicator-vector matvec
  (`ind @ mat`) is faster still but lands 7.5e-16 away, and summing the CSC
  selection directly 2.5e-15 away, because each accumulates the columns in a
  different order.
- **The same seed gave a different module score in a different process.** The
  control genes were collected in a `set`, and a mean depends on the order its
  terms are added. Python randomises `str` hashing per process, so iterating
  that set summed the control expression in a different order every run — the
  same object at `seed=1` scored 9.7e-16 apart in two processes. A `dict` now,
  which is also R's semantics (`AddModuleScore` applies `unique()` to the
  sampled names and indexes the matrix with the result). THP-1 `S.Score`,
  `G2M.Score` and the interferon program move by at most 2.1e-15 against the old
  code — which had no fixed value to move from — and all 20,729 `Phase` calls
  are unchanged.
- **A reduction that could not use every feature you asked for labelled its
  loadings with the features you asked for anyway.** `_get_scaled_data` filtered
  the request down to what the layer carried and returned only the matrix, so
  `run_pca` stored *n* names against fewer than *n* rows and every row below the
  first dropped feature answered to the wrong gene. `viz_dim_loadings` and
  `dim_heatmap` read those names straight onto the axis, and `jack_straw` handed
  each per-feature p-value to its neighbour. Nothing raised; every number was
  individually correct. Seurat's `PrepDR5` drops the same features but **warns
  naming them**, and takes its rownames from the subset matrix so the labels
  cannot drift from the rows. Both halves are ported: the drop now warns with the
  gene names (`RuntimeWarning`, first 20 named then a count), requesting nothing
  the layer holds is an error rather than an empty matrix, and `run_pca`,
  `run_spca`, `run_ica` and `jack_straw` all take back the list that labels the
  rows they got. On the ordinary path — PBMC 3k, 2,000 variable features, nothing
  dropped — embeddings and loadings are **bit-identical** to before.
- **`run_pca` on a v3 `Assay` that had been through `scale_data()` read the
  wrong genes, or crashed.** `Assay.scale_data` is a bare ndarray holding only
  the scaled subset — the variable features, by default — with no record of
  *which* subset; R's slot is a matrix and carries its own rownames. Four places
  invented the missing labels, differently: `_get_scaled_data` indexed the layer
  by each gene's position in the *assay*, `features("scale_data")` returned the
  assay's leading *n* features, `subset()` discarded the layer whenever it was
  not full height (i.e. always), and `mapping._scaled_feature_names` returned the
  full feature list. The first either raised `IndexError` or returned another
  gene's row, deterministically, on the standard workflow. `Assay` now carries
  `_scaled_features` alongside `scale_data`, mirroring `Assay5`; the constructor
  refuses to guess labels for a subset rather than inventing them, `subset()`
  subsets the layer by name, and the v3 and v5 paths now return **identical**
  scaled matrices.
- **Integration anchors intersected the feature set against the wrong list.**
  `_integration_features` filtered against the assay's full feature list — while
  its own comment claimed it checked "scale/data" — and `_anchor_feature_matrix`
  filtered against the scaled layer, so a feature one object had never scaled
  survived the first filter and was dropped by the second, **in that object
  alone**. The reference and query matrices are multiplied row-for-row, so that
  is either a shape crash or two matrices describing different genes. The
  intersection is now taken against the layer the anchors are built from, and
  `_anchor_feature_matrix` raises rather than warns if anything is still missing.
- **`run_pca` on an Assay5 that had not been scaled raised "the truth value of
  an array … is ambiguous".** `_get_scaled_data` chose its fallback layer with
  `layers.get("data") or layers.get("counts")`. `or` evaluates `bool()` on its
  left operand, and scipy refuses to answer for any matrix with more than one
  entry — so the fallback raised on **every call that reached it**, which is
  every `run_pca` or `run_ica` on a v5 object that had not been through
  `scale_data()`. The v3 path never had it: it reads `assay_obj.data`
  directly, so the two architectures silently disagreed about whether an
  unscaled object could be reduced at all. Selected with an explicit `is None`
  now, and if neither layer exists the error says which ones it looked for.
  The v5 fallback returns **bit-identical** embeddings to the v3 one it mirrors.
- **Merging a v3-backed and a v5-backed object died with an
  `AttributeError` about a private slot.** `Assay` keeps its cells in
  `_cell_names` and `Assay5` in `_all_cell_names`, so `Truecell.merge` reached
  for whichever the other did not have and failed partway through — after the
  cell names and metadata had already been concatenated, from inside a method
  the caller never named. It now raises a `TypeError` naming the assay and both
  classes, and pointing at `create_truecell_object(..., use_v5=)` as the way to
  make the two sides match. Merging assays of the same class is unaffected.
- **`StdAssay`'s methods returned the abstract base rather than the caller's
  own class.** `join_layers`, `split_layers`, `cast_assay`, `subset`,
  `rename_cells`, `merge` and `_copy` were annotated `-> "StdAssay"` while every
  one of them builds its result from `self.__class__`, so an `Assay5` went in
  and something typed as an ABC came out. Every caller downstream then had to
  widen, which is what put `dict[str, Assay | StdAssay]` where `Truecell` wanted
  `dict[str, Assay | Assay5]`. Annotated `Self`, which is both accurate and what
  the code already did. `tests/test_object_model_typing.py` pins the runtime
  property that makes it true — a subclass gets its own class back from all
  seven — since hardcoding `Assay5(...)` in any of them would leave the
  annotation a promise the code does not keep and no other test would notice.
- **An unreachable second implementation of the object's count columns is
  gone.** `create_truecell_object` dispatched on `hasattr(assay_obj, "calc_n")`,
  which both assay classes satisfy, so the `_calc_n_for_assay5` fallback could
  never run. It derived the `nCount_<ASSAY>` suffix from the assay's *key*
  rather than from the `assay` argument and hardcoded `nCount_RNA` for an assay
  with no default layer — agreeing with the live path for an ordinary RNA
  object, which is why nothing had reason to look at it.

  Together these clear **all 14 mypy errors** in the object-model modules
  (`truecell.py`, `assay5.py`, `reduction.py`, `spatial/fov.py`), taking the
  advisory total from **60 to 46**. The remaining 46 are elsewhere: 40 in
  `datasets.py` (one repeated idiom — a `cache_dir: Optional[str]` parameter
  reassigned to a `Path`), and six spread over plotting, mixscape, multiseq,
  sctransform, preprocessing, module_score and integration.
- **`truecell/command.py` began with a UTF-8 BOM.** CPython decodes source as
  `utf-8-sig`, so the file imported, introspected, linted and type-checked
  without complaint — which is how three bytes survived unnoticed. What they
  broke is reading the source back as text: `Path.read_text()` defaults to plain
  utf-8, keeps the U+FEFF, and `ast.parse` then rejects the file outright with
  `invalid non-printable character`. The AST walk added alongside this had to
  decode around it, and any codemod, doc generator or source-level lint would
  hit the same wall on a file that looks fine. The BOM is gone, and
  `test_no_source_file_starts_with_a_byte_order_mark` scans the whole repo — not
  just the package — so another cannot arrive quietly.

- **Twenty-one annotations named symbols that existed nowhere at module scope.**
  Seventeen plotting functions were declared `-> "plt.Figure"` while `plt` was
  only ever a local inside each body (`plt = _mpl()`); `Graph.as_neighbor`,
  `Neighbor.as_graph` and `from_anndata` named classes imported at call time to
  dodge a circular import. `from __future__ import annotations` makes every
  signature an unevaluated string, so this cost nothing at runtime and no test
  noticed — but it is exactly what a type checker, an IDE, or a documentation
  generator reads. Both mypy (`name-defined`) and ruff (`F821`) had been
  reporting all 21 for as long as the annotations existed.

  Each is now a `if TYPE_CHECKING:` import, which those tools read and the
  interpreter never runs. The plotting signatures say `-> "Figure"` from
  `matplotlib.figure`, which is the class the functions actually return.
  matplotlib stays an optional dependency: `import truecell` eagerly imports
  `truecell.plotting` and still does not pull it in.

  mypy drops **81 → 60** errors and ruff **72 → 51** on `truecell` (222 → 201
  repo-wide) — the same 21 in both counts, so any baseline recorded before this
  is 21 high on both tools. `tests/test_annotations_resolve.py` walks the AST of
  every module in the package and fails on any annotation whose root name is not
  bound at module scope, and separately asserts each deferred import stays
  *inside* its `TYPE_CHECKING` block — hoisting it would satisfy both checkers
  while making matplotlib mandatory.

- **`find_all_markers` was missing Seurat's `return.thresh`, and did not break
  p-value ties.** `FindAllMarkers` defaults to `return.thresh = 1e-2` and
  returns only genes below it; truecell returned everything that survived the
  `min_pct` and `logfc_threshold` pre-filters, including plainly
  non-significant rows. The PBMC 3k handoff found it the decisive way: two
  clusters came out with *identical* cell membership on both sides, and on
  those two truecell returned 190 and 383 genes against Seurat's 151 and 242 —
  every extra row a gene the two agreed about numerically (max `avg_log2FC`
  difference 4.9e-15) and that Seurat simply does not return. Applying the same
  filter reproduces Seurat's gene sets **exactly**, 151/151 and 242/242.

  Rows are now also ordered by `p_val` ascending then `avg_log2FC` descending
  within each cluster, matching `order(gde$p_val, -gde[, 2])`. This is not
  cosmetic: Wilcoxon p-values tie at exactly 0 for the strongest markers — 40
  to 302 genes per cluster on PBMC 3k — so without the tie-break "the top N
  markers", which is what every tutorial prints, was decided by incoming row
  order. Pass `return_thresh=None` for the old unfiltered table.

  Cluster labels are also now iterated numerically rather than
  lexicographically, so a dataset with eleven or more clusters comes back in
  Seurat's factor order instead of 0, 1, 10, 11, 2, …

*The neighbour graphs, against `FindNeighbors` / `FindClusters` (Seurat 5.5.1).
Four defects, found while establishing that the clustering divergence left open
above is **not** a defect — given the same embedding the two tools' KNN indices
are identical and their SNN graphs agree to 2.8e-08, and the partitions differ
only because Seurat's `n.start = 10` restarts find 0.17 % more modularity by
splitting CD14 Mono along the batch axis (73.8 % CTRL / 83.3 % STIM). truecell
scores ARI 0.9195 to `seurat_annotations` against Seurat's 0.7368, so the
Louvain search was deliberately left alone.*

- **`find_neighbors` symmetrised the KNN graph.** Seurat's `nn` is the raw
  ranked neighbour table — `sparseMatrix(i, j, x = 1)`, no symmetrisation — so
  `nnz` is exactly `n * k`, every row sums to `k`, and the column sums vary
  with in-degree (21 to 68 on ifnb). `mat + mat.T` inflated `nnz` and flattened
  that spread, discarding which cells are hubs. Both in-tree consumers
  (`run_umap`, the graph branch of `run_spca`) already symmetrise at the point
  of use, as Seurat's own do.
- **`find_neighbors` deleted the SNN diagonal.** `ComputeSNN` stores
  `SNN[i,i] = k / (2k − k) = 1` for every cell; all 13,999 were present in
  Seurat's graph and none in truecell's. It was invisible to `find_clusters`,
  whose igraph conversion takes the strict upper triangle and so discarded
  precisely the entries that were missing.
- **`run_umap` did not strip that diagonal.** `RunUMAP.Graph` opens with
  `diag(x = object) <- 0`. With the diagonal restored above, omitting this
  would feed the layout `n` zero-length self-edges.
- **`find_neighbors` computed Jaccard weights in float32.** ~3e-08 off Seurat
  on every weight, now 4.4e-16. This never changed which edges were pruned — a
  weak Python `prune_snn` is cast down to float32 for the comparison, so both
  sides round identically — only the stored values.

### Added

- **`find_clusters(group_singletons=True)`**, porting Seurat's
  `GroupSingletons`. Size-1 clusters are absorbed into whichever non-singleton
  cluster they are most connected to, scored by *mean* SNN weight (a sum would
  hand them to the largest cluster instead), with the candidate list fixed
  before the loop so one singleton cannot absorb another. `False` pools them
  into a single `"singleton"` cluster, again matching Seurat.

### Fixed

*The v4 anchor path, against `FindIntegrationAnchors` / `IntegrateData`
(Seurat 5.5.1). Twelve defects; anchor agreement went from 70.0 % to 99.9 %
for CCA and to 100 % for RPCA.*

- **`find_integration_anchors(reduction="cca")` built the wrong shared space.**
  `RunCCA` **standardizes** each cell (z-score) before the cross-covariance SVD;
  truecell L2-normalized it. That is a correlation matrix between cells versus a
  cosine-similarity one — different singular vectors, so every anchor moved.
- **Anchor features constant in either object are now dropped**, as
  `RunCCA`'s `CheckFeatures` does (83 of 2,000 on ifnb). Standardizing works
  down each *cell*, so a constant gene still shifted that cell's mean and SD.
- **The anchor filter now runs on `TopDimFeatures` of the `data` layer.**
  `FilterAnchors` uses at most 200 genes chosen from the CCA loadings and the
  log-normalized values; truecell used all 2,000 anchor features of `scale.data`.
- **The anchor score now uses four neighbour tables.** `ScoreAnchors` gives each
  member `k.score` neighbours within its own dataset **plus** `k.score` in the
  other. A single pooled kNN is nearly all same-batch when a batch effect is
  present, collapsing every score toward the floor.
- **Filtering now happens before scoring.** The score is rescaled against the
  1st/90th percentiles of the set it is handed, so scoring first took those
  percentiles from anchors that were then discarded — shifting every surviving
  value (mean 0.5477 vs Seurat's 0.4971, only 10.5 % identical).
- **`k_filter` now retains every anchor** when either dataset is smaller than
  it, as Seurat does, instead of clamping `k` to the query size.
- **`_pca_loadings` now uses an exact SVD.** sklearn's `PCA` switches to a
  *randomized* solver above `max(shape) > 500`, which drifts in the trailing
  components — only 12–14 of 30 PCs matched irlba above 0.99. Reciprocal PCA
  standardizes each projected dimension by its own SD, which is not
  rotation-invariant, so a drifted axis became a different reciprocal space:
  RPCA anchor recall **44.9 % → 100 %**.
- **`integrate_data` now weights in a PCA of the merged pair.**
  `RunIntegration` re-scales the reference and query together on the anchor
  features and runs a fresh PCA; truecell reused the CCA embedding, a space built
  to make the batches overlap.
- **The weight kernel now matches `FindWeightsC`** —
  `1 − exp(−d̃ · score / (2/sd)²)` over `d̃ = 1 − d/dₖ`, with the score inside
  the exponent — instead of a Gaussian in the raw distance times the score.
  `sd_weight` now widens the kernel as it grows, which is the opposite of what
  the old bandwidth did.
- **`k_weight` counts anchors, not anchor cells.** The neighbour search runs
  over the **unique** query anchor cells and expands each into all of its anchor
  rows, stopping at `k_weight` entries.
- **`integrate_data` raises when there are fewer distinct query anchor cells
  than `k_weight`**, as `FindWeights` does, instead of silently shrinking `k`
  and averaging over far fewer anchors than asked for.
- **`integrate_layers` now corrects onto the larger batch.**
  `PairwiseIntegrateReference` reverses the merge pair when the second object is
  bigger; truecell always took the first. Invisible on an even split, and ifnb is
  CTRL 6,548 vs STIM 7,451.

These also reach `find_transfer_anchors` / `transfer_data`, which share the
scoring, filtering and weighting helpers: `panc8` label-transfer accuracy went
0.9845 → 0.9862 and per-cell concordance with R 0.9871 → 0.9883.

*The v5 `IntegrateLayers` path, against `CCAIntegration` / `RPCAIntegration`
(Seurat 5.5.1). Two more defects, found by treating what had looked like an
"expected implementation gap" on the v4 path (above) as a claim to check.
Embedding agreement on a 2,400-cell probe went from 1/30 PCs above \|r\| = 0.99
to 30/30 for both reductions; on the full 13,999-cell, unequal-batch ifnb,
RPCA batch mixing rose from 0.867 to 0.991 (Seurat: 0.917).*

- **`integrate_layers(method="cca"|"rpca")` was running the wrong algorithm.**
  Seurat's `IntegrateLayers(method = CCAIntegration/RPCAIntegration)` does not
  call `IntegrateData` — it calls **`IntegrateEmbeddings`**, which corrects the
  input PCA embedding directly (transposed into a fake per-dimension assay and
  pushed through the same anchor-weighting machinery) and returns it in the
  same basis, loadings included. truecell was running the v4 workflow behind
  that name — correct expression, re-scale, re-run PCA — landing in a
  different basis with the same shape. Added `integrate_embeddings`
  (`truecell.anchors`), which `integrate_layers` now calls; `k_filter=None` for
  both methods (Seurat forces `k.filter <- NA` on this path) and only
  `rpca` re-scales each batch (`CCAIntegration` slices the object's existing,
  pooled `scale.data`; `RPCAIntegration` runs `ScaleData` per batch).
- **`run_pca` used sklearn's randomized SVD.** It switches solvers once
  `max(shape) > 500`, matching only 15 of 30 PCs above \|r\| = 0.99 against
  Seurat's irlba (one PC at 0.006) — the same drift `_pca_loadings` (above) had
  already needed an exact SVD to avoid. Invisible while only leading PCs were
  read downstream; not invisible once `integrate_embeddings` corrects the
  embedding itself. Replaced with ARPACK (`scipy.sparse.linalg.svds`, seeded
  for determinism), which matches irlba to six decimals and runs 6× faster on
  data this shape; small inputs still take an exact dense SVD, since ARPACK
  misbehaves as `k` approaches the matrix rank. Also fixed to match
  `RunPCA.default`: `sdev = d / sqrt(ncol(object) - 1)` rather than the
  embedding's own SD, and no re-centring of `scale.data`.

### Changed — BREAKING

*The Visium loader, aligned to `Read10X_Image` / `Load10X_Spatial`*

- **`load_visium` now filters to in-tissue spots by default.**
  `filter_by_tissue` defaults to `True`, matching Seurat's
  `filter.matrix = TRUE`. On a `filtered_feature_bc_matrix` bundle nothing
  changes — those barcodes are already tissue-filtered — but on a
  `raw_feature_bc_matrix` the spot count drops from every spot on the capture
  area to the ones under tissue (4,992 → 2,695 on the reference slide).
- **`load_visium` now reads the lowres image by default.** `image_resolution`
  defaults to `"lowres"`, matching `Read10X_Image`'s
  `image.name = "tissue_lowres_image.png"`. The stored array is ~11× smaller;
  `scale_coordinates()` and `spot_radius()` follow the stored resolution, so
  anything reading them in image pixels changes scale accordingly.
- **The image key is now `"slice1"`, not `"spatial"`.** `obj.images["slice1"]`
  is what `Load10X_Spatial` produces and what every ported Seurat script
  indexes; it used to raise `KeyError`. Override with `slice_name=`.
- **`get_tissue_coordinates()` on a `Centroids`, `Segmentation` or FOV now
  returns a `cell` column** alongside `x`/`y`, with the cells still on the
  index — the frame R returns. Callers taking `.values`/`.to_numpy()` over the
  whole frame now get an object array and should select `[["x", "y"]]`.

  The previous behaviour is still available:
  `load_visium(path, image_resolution="hires", filter_by_tissue=False, slice_name="spatial")`.

### Fixed

*The Visium container, audited against Seurat 5.5.1*

- **The tissue image depended on which optional package was installed.**
  `_imread` tried `matplotlib.image` and fell back to `PIL`, which return
  float32 in [0, 1] and uint8 in [0, 255] respectively — **arrays 255× apart,
  with different dtypes, from the same PNG**. Neither library is a declared
  dependency, so `VisiumV2.get_image()` was a function of the environment as
  much as of the file. Plotting never revealed it because `imshow` accepts
  both. Both backends now return the same float array, matching
  `png::readPNG` to float32 epsilon.

  Two Seurat behaviours are **reported and deliberately not matched**: it
  stores `spot_diameter_fullres` in the FOV's `radius` (truecell keeps
  `diameter / 2`; the slide's fixed 100 µm spot pitch shows the field is a
  diameter, since read as a radius the capture areas would overlap by 31 µm),
  and `Radius()` on a `VisiumV2` returns `NULL` because no `Radius.VisiumV2`
  method exists. See [`tutorials/visium_vignette.md`](tutorials/visium_vignette.md).

*The out-of-core path, audited against Seurat 5.5.1 + BPCells 0.3.1*

*The out-of-core path, audited against Seurat 5.5.1 + BPCells 0.3.1*

- **Five functions read an on-disk layer by densifying all of it.**
  `_log_normalize`, `_vst_hvg`, `scale_data`, `find_markers` and
  `add_module_score` each took a dense fallback that ran
  `np.asarray(whole)[idx]` where subsetting first would have been cheap —
  `LazyMatrix.__getitem__` already returns a scipy block. The effect was that
  backing a matrix on disk made things **worse**: `normalize_data` on pbmc3k
  used **4.6× the peak memory** of the sparse path (1169 MB against 253 MB) and
  left a dense `ndarray` **16× larger** than the sparse layer it replaced (296
  MB against 18.3 MB). `col_blocks`, documented as "the primitive for an
  out-of-core reduction", had no callers. Now zero whole-store materialisations
  across the whole pipeline.
- **`percentage_feature_set` raised on an on-disk layer.**
  `sp.issparse(LazyMatrix)` is `False`, so it took the dense branch — but
  indexing a `LazyMatrix` returns scipy, whose `.sum(axis=0)` is a `(1, n)`
  matrix rather than a vector, and the result reached `pd.Series` with the wrong
  shape.
- **Object construction ended laziness before any analysis ran.**
  `create_assay5_object` ran `sp.csc_matrix(np.asarray(matrix))` on anything not
  already scipy, and `calc_n` densified again for `nCount`/`nFeature`. So
  opening a store and building an object on it — the obvious way to use the
  feature — was the one path that could not work. Seurat has a
  `.CalcN.IterableMatrix` for the same reason. `LazyMatrix.nnz_per_row` added to
  support `min_cells` filtering without materialising.
- **`_loess2` was chaotically dependent on sort order.** On pbmc3k **85.5 % of
  genes share their `log10(mean)`** (13,714 genes over 2,837 distinct values,
  largest tied run 627). `np.argsort` defaults to unstable quicksort and the
  LOESS window was chosen by position, so members of a tied run got different
  neighbourhoods and different fitted values — a spread of **1.3e-3** within a
  single tied `x`, where R's `loess` gives 1.8e-15 because a fitted value is a
  function of `x` alone. Perturbing the input by **1e-15** moved fitted values by
  up to **28.8 %**. The fit is now evaluated once per distinct `x`, with
  distance-based neighbourhoods and a lexicographic sort. Checked against R
  before adoption: median \|diff\| vs `loess` 0.000333 → **0.000103**, and HVG
  overlap with Seurat 99.65 % → **99.90 %** (all 6 genes it added are ones
  Seurat picks).
- **Both HVG selectors broke ties the wrong way.** `argsort(v)[::-1]` orders
  ties by descending index; R's `head(order(x, decreasing = TRUE), n)` orders
  them ascending, and unstably at that. Now `argsort(-v, kind="stable")`.
- **The sparse and on-disk paths were two implementations agreeing to 1e-14.**
  That is not enough when a *tie-break* consumes the statistic:
  `variance.standardized` carries exact ties, and genes tied under one summation
  order are not tied under the other. It reordered **147 of 2000** features and
  produced 9 clusters against 8. Both layer types now share one block reduction,
  so the two paths are bit-identical by construction.

*The differential-expression test suite, audited against Seurat 5.5.1*

- **`avg_log2FC` put Seurat's pseudocount on the group mean instead of the group
  sum.** Seurat 5's `log1pdata.mean.fxn` is `log2((sum(expm1(x)) + 1) / n)`, so
  the pseudocount is worth `1/n` on the mean scale; truecell added a whole count to
  the mean, which is Seurat 4's formula. Every fold change was floored toward
  zero — a gene detected in 0 % of one cluster and 24 % of the other read −1.26
  against Seurat's −9.92. **`logfc_threshold` filters on this value**, so the
  error changed which genes were returned, not only what they were labelled:
  2,298 genes against Seurat's 11,931 at a 0.25 threshold (Jaccard 0.193). Now
  matches R to **7.11e-15** across all 13,712 genes.
- **`negbinom` ran a likelihood-ratio test on a moment-estimated dispersion**
  where Seurat's `GLMDETest` fits `MASS::glm.nb` (dispersion by maximum
  likelihood) and takes the Wald p-value. HLA-DRA read 5.5e-128 against R's
  1.1e-321. After the fix, p-values agree exactly for every gene detected above
  5 %; below that the GLM is fitting near-empty rows, and Seurat's
  `min.cells.feature` drops those genes anyway.
- **`test_avg_log2fc_matches_seurat_formula` encoded the same wrong formula** and
  checked truecell agreed with itself, so it was green the whole time under a name
  that claimed Seurat parity. Corrected.
- **The `latent_vars` docstring had MAST backwards** — it advised passing the
  cellular detection rate "to match Seurat's default CDR covariate", but
  `MASTDETest` fits `~ condition` alone and adds no CDR term. Passing CDR is a
  deliberate departure from Seurat, not a way to match it.

*Spatial statistics and the spatial container, audited against Seurat 5.5.1*

- **`find_spatially_variable_features(method="moransi")` used the wrong spatial
  weights.** Seurat builds `1/d²` between every pair of cells and
  `Rfast2::moranI` row-standardises it; truecell used a k-nearest-neighbour graph.
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
  mouse brain). truecell left it `None`, and because `_spot_collection` returns
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
- **`truecell.generics.split_layers` was declared but never registered** for any
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
- **The command log was inert.** `log_truecell_command` was a public export with
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
  singular vectors, so the scores sum to 50. Truecell used every direction above a
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
  from that refit basis; truecell projected the permuted rows onto the *existing*
  embedding. A fixed basis cannot rotate to absorb the scrambled signal, so the
  permuted loadings came out too small and the null was far too tight. On pbmc3k's
  pure-noise PCs 14-20 that put **109-203** of 2000 features below p ≤ 1e-5, where
  R finds **0-5**. Now refits the PCA per replicate, as R does.
- **`score_jackstraw` used a KS test instead of `prop.test`.** R's `ScoreJackStraw`
  tests the count of features below `score.thresh` against the count expected under
  a uniform null; truecell ran a one-sided Kolmogorov-Smirnov test against
  Uniform(0, 1), which on thousands of features is enormously more sensitive — its
  **largest** score across all 20 pbmc3k PCs was `8.1e-112`, so no PC ever failed
  the threshold. R's `prop.test` is now ported exactly (Yates-corrected two-sample
  chi-square), reproducing R to nine significant figures from 1e-143 to 1.0.
- **Net effect:** truecell recommended keeping **all 20** PCs where Seurat keeps 13 —
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

*Differential-expression tutorial — Wave 3's first*

- `tutorials/pbmc3k_de_tutorial.py` + `pbmc3k_de_verify.R` + `de_vignette.md` +
  `figures_de/`. All eight `find_markers` tests against `FindMarkers` on pbmc3k
  clusters 0 vs 1, on a cell assignment exported from the Python side so that no
  clustering difference can be mistaken for a DE difference. **Eight tests, none
  of which had ever been compared to R.**
- After the two fixes, **all seven per-cell tests reproduce Seurat's top 50 genes
  exactly**, with `wilcox` / `t` / `bimod` / `LR` at p-value Spearman ≥ 0.99997
  and `mast` at 0.9993 on genes detected above 5 %.
- Reported rather than changed: `deseq2` is pseudobulk where Seurat tests cells
  as replicates (and requires `sample_col`, so it raises rather than silently
  substituting); Seurat rounds `myAUC` to 3 dp inside `DifferentialAUC`; R's
  `wilcox` returns `NaN` for genes expressed in neither group where truecell
  returns `p = 1`.

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
  as three columns, truecell carries the cell as the DataFrame index.
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
  fixture reproduces it: several were tried, and R agrees with truecell on those to
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
  Barnes-Hut, truecell calls scikit-learn): each embedding retains **0.470** / **0.477**
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
  tool is ~98.5 % accurate against the held-out cell types (truecell 0.9845, Seurat
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
- Added to the opt-in tutorial smoke suite (`TRUECELL_TUTORIAL_SMOKE=1`) and covered
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
- Added to the opt-in tutorial smoke suite (`TRUECELL_TUTORIAL_SMOKE=1`) and covered
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
- Added to the opt-in tutorial smoke suite (`TRUECELL_TUTORIAL_SMOKE=1`) and covered
  by `tests/test_hashing_tutorial.py` (network-free: species/concordance helpers
  and the load→demux→figure path on a synthetic barnyard). Suite 489 → 496. No new
  `pip` dependencies.

*Tutorial data infrastructure — the R-side scaffolding for expanding tutorial coverage (#38)*

- `truecell.datasets.pbmc_hashing` (GSE108313) and `thp1_eccite` (GSE153056) —
  loaders for the Cell-Hashing and ECCITE-seq/Mixscape datasets, parsed straight
  from their original GEO plain-text files so R and Python read identical counts.
  The ECCITE loader also returns the per-cell guide/replicate metadata, so a
  Mixscape tutorial can start from the same annotated state as R's `thp1.eccite`.
- `truecell.datasets.ifnb` and `panc8`, with `tutorials/export_seuratdata.R` — a
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
- `truecell._clara` and `hto_demux(kfunc="clara")` — an in-tree port of R's
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
  Verified against R: truecell `margin=2` reproduced Seurat `margin=1` and vice
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
  `"v1"` Python and R both resolve 13 clusters; at the `"v2"` default Truecell
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
  `TRUECELL_TUTORIAL_SMOKE=1`.

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

First release published to PyPI — `pip install truecell` works from here on.

### Added

- Spatial Seurat parity: loaders, neighborhood analysis, niches, and
  composition. (#4)
- Xenium spatial tutorial, verified against R Seurat to 8 significant
  figures. (#4)
- [`ROADMAP.md`](ROADMAP.md) — the milestone plan. (#2)
- `truecell/py.typed` — the package ships PEP 561 type information, so a
  downstream `mypy` reads truecell's annotations rather than treating it as
  untyped. (#4)

### Fixed

- Pin `numba>=0.59` so the `[analysis]` and `[all]` extras resolve on Python
  3.10+. (#2)
- The README Quick Start example now produces the 500 cells it claims. (#3)

## [0.1.0] - 2026-06-30

Initial release: a port of Seurat's core data structures and analysis pipeline.
Tagged but never published to PyPI; `0.1.1` was the first release on PyPI.

### Added

- Core objects: `Truecell`, `Assay`, `Assay5`, `StdAssay`, `DimReduc`, `Graph`,
  `Neighbor`, `LogMap`, `JackStrawData`, `TruecellCommand`.
- Pipeline: `normalize_data`, `find_variable_features`, `scale_data`,
  `percentage_feature_set`, `run_pca`, `find_neighbors`, `find_clusters`,
  `run_umap`, `find_markers`.
- `sctransform`, `module_score`, `jack_straw`.
- Spatial primitives: `FOV`, `Centroids`, `Segmentation`, `Molecules`.
- Plotting, I/O, AnnData interop, and the bundled example datasets.

[Unreleased]: https://github.com/GenomicAI/truecell/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/GenomicAI/truecell/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/GenomicAI/truecell/compare/v0.9.0...v1.0.0
[0.9.0]: https://github.com/GenomicAI/truecell/compare/v0.2.0...v0.9.0
[0.2.0]: https://github.com/GenomicAI/truecell/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/GenomicAI/truecell/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/GenomicAI/truecell/releases/tag/v0.1.0
