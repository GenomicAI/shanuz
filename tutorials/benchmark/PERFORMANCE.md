# Truecell vs R Seurat — Performance

The companion to the accuracy comparison in [`../README.md`](../README.md).
That work asked whether the two tools produce the same answers. This one asks
what each answer costs.

**Machine.** Apple M4 Pro, 12 cores, 25.8 GB, macOS 26.5.2. Measured 3 Aug 2026.

**Versions.** R 4.6.1 (**linked against Accelerate/vecLib**) · Seurat 5.5.1 ·
Matrix 1.7.5 · presto 1.0.0 · harmony 2.0.5 · uwot 0.2.4 · irlba 2.3.7 ·
RANN 2.6.2 · Rfast2 · data.table 1.18.4 — Python 3.12.13 · truecell 1.0.0 ·
numpy 2.4.6 (Accelerate) · scipy 1.18.0 · umap-learn · scikit-learn · harmonypy.

Reproduce with `bash tutorials/benchmark/sweep.sh` (about 40 minutes).

> **Both stacks now sit on the same BLAS.** An earlier version of this report
> measured an R that linked the unoptimised reference BLAS it ships with, and
> flagged that as its largest caveat. That has been fixed — see
> [section 2.1](#21-both-arms-are-on-accelerate-and-what-that-changed) for what
> it moved, which was more than expected and reversed one of the report's
> conclusions. The reference-BLAS sweep is kept in `results_refblas/` and the
> two are diffable with `compare_blas.py`.

---

## Summary

Across 68 like-for-like step comparisons Truecell is **faster in 46, slower in
21, and within 5% in 1** — but it loses the standard end-to-end workflow by
1.7–2.7x. Three operations account for nearly all of that, and all three have
identifiable causes: marker detection densifies before it filters, umap-learn
single-threads under a seed, and Seurat's Wilcoxon is presto's C++.

| | Truecell | Seurat |
|---|---|---|
| Standard workflow, 2.7k–20.7k cells | | **1.7–2.7x faster** |
| …with `run_umap`'s seed dropped | | **1.2–2.0x faster** |
| Reading counts, building the object | **4–20x faster** | |
| PCA | **1.5–2.2x faster** | |
| Normalisation, QC metrics | **2–13x faster** | |
| VST feature selection, scaling | | **1.3–2.5x faster** |
| Differential expression, 6 of the 7 shared tests | **1.9–9.5x faster** | |
| Wilcoxon — `de_wilcox` and `find_all_markers` | | **3.2–5.9x faster** (presto) |
| Seeded UMAP | | **1.7–4.2x faster** |
| Harmony | **2.1x faster** | |
| CCA / RPCA integration | | **1.1–3.6x faster** |
| SCTransform | **1.3x faster**, half the memory | |
| Moran's I | **87x faster**, and the only one without an *n* limit | |
| Peak memory, standard workflow | | **2.5x lighter** (all of it one step) |
| Peak memory, THP-1 and the spatial/DE benches | **1.04–14x lighter** | |

---

## 1. How this was measured

Both arms run the same pipeline on the same bytes with the same parameters,
step for step. `tutorials/benchmark/bench_truecell.py` and `bench_seurat.R` are
line-for-line counterparts; `run_benchmarks.py` runs them and collects the
numbers.

**Time is measured inside each process. Memory is measured from outside.** R's
`gc()` accounting and Python's `tracemalloc` measure different things, and
neither sees what the other's allocator is holding, so a memory number produced
by either would not be comparable. Instead the parent process samples the
resident set size of the child's whole process tree every 50 ms and afterwards
intersects those samples with the step boundaries the child logged. Both arms
are measured by one instrument, in one unit, on one clock.

Resident set is a high-water mark — neither runtime returns freed pages
promptly, and R's GC returns them more readily than CPython's arenas do — so a
late step's peak partly reports what earlier steps left behind. Per-step figures
are peaks *during* the step; the pipeline row is the process peak, which is the
number that decides whether a machine can run the workload.

**The first repeat of every pair is discarded.** umap-learn's numba kernels
compile on first use and R loads packages lazily; neither is what the benchmark
is asking about. Three timed repeats follow, and the median is reported.

**Every step records an anchor** — cells kept, clusters found, genes tested,
markers returned — printed in the tables beside the timings. A speed comparison
is only worth reading if both sides did the same work, and the anchors are how
you can check rather than take it on trust. They earned their keep when the BLAS
was swapped: all 75 of them were unchanged, which is what made the new timings
safe to believe.

**The instrument is not free.** Sampling costs one `ps` call every 50 ms, a few
percent of one core out of twelve. It is applied identically to both arms, so
the comparisons hold, but treat the absolute seconds as very slightly inflated
on both sides.

**The truecell arm runs first, on purpose.** It writes the cell-to-cluster
assignment and the Xenium cell subset that the R arm reads back. Without that,
two things would have gone wrong: the tools do not always land on the same
number of clusters (9 against Seurat's 12 on PBMC 8k), and one-vs-rest marker
detection costs one test per cluster, so `find_all_markers` would have been
timing a clustering difference. With it, both arms return the same marker table
— 3118 rows against 3118 on PBMC 3k.

---

## 2. Three things that decide most of these numbers

### 2.1 Both arms are on Accelerate, and what that changed

macOS R ships two BLAS builds and symlinks the unoptimised one by default. This
machine was on that default until the sweep below; it is now on Accelerate:

```bash
cd /Library/Frameworks/R.framework/Resources/lib && ln -sf libRblas.vecLib.dylib libRblas.dylib
# revert with: ln -sf libRblas.0.dylib libRblas.dylib
```

`blas_probe` is a control — dense linear algebra on identical inputs, touching
neither Seurat nor truecell:

| Step | Truecell | Truecell (1 thread) | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.46s | 0.47s | — | — | 123 MB | — |  |  |
| `blas_setup` | 0.01s | 0.02s | 0.18s | **12.2x** | 124 MB | 528 MB | 2000 | 2000 |
| `blas_gemm` | 0.01s | 0.02s | 0.13s | **9.3x** | 199 MB | 602 MB | 2002.361 | 2006.143 |
| `blas_svd` | 0.02s | 0.02s | 0.03s | **1.4x** | 199 MB | 605 MB | 66.575 | 67.223 |
| `blas_crossprod_chol` | 0.03s | 0.04s | 0.12s | **4.4x** | 199 MB | 660 MB | 62.883 | 63.745 |
| `library_load` | — | — | 1.38s | — | — | 486 MB |  |  |
| **shared pipeline** | **0.1s** | **0.1s** | **0.5s** | **6.0x** | **219 MB** | **660 MB** |  |  |

**Nothing Seurat computes changed.** Every one of the 75 step anchors across
the whole suite — cluster counts, graph edge counts, PC standard deviations,
marker counts, gene counts — is identical between the two sweeps. The swap is
free of numerical consequence at the precision the suite records.

**What it bought, per step** (`compare_blas.py results_refblas results`):

| Step | Reference BLAS | Accelerate | |
|---|---:|---:|---|
| `integrate_cca` (ifnb) | 88.41s | 6.87s | **12.9x** |
| `pca` (pbmc3k / 8k / ifnb / thp1) | 2.83–13.07s | 0.32–1.79s | **5.5–8.8x** |
| `pca_on_sct` | 1.25s | 0.15s | **8.1x** |
| `hvg_vst` | 0.40–1.74s | 0.15–1.32s | **1.3–2.7x** |
| `prep_to_pca` | 5.73s | 2.52s | **2.3x** |
| `integrate_rpca` | 11.91s | 8.41s | **1.4x** |
| everything else | — | — | unchanged |

Read that last row: neighbours, clustering, UMAP, marker detection, the DE
tests and SCTransform did not move at all, because none of them is BLAS-bound.

**This reversed a conclusion.** The reference-BLAS report called batch
integration "the largest clean win, and the one least contaminated by the
BLAS" — truecell's CCA at 24.3s against Seurat's 88.4s. It was in fact the
*most* contaminated result in the suite: on the same BLAS, Seurat's CCA is
6.87s and truecell is **3.6x slower**. The claim was wrong, and it was wrong in
the direction of flattering this project.

One caveat on the probe itself: `blas_gemm` writes `a %*% t(a)`, which
materialises the transpose, where numpy's `a @ a.T` passes a flag to dgemm and
copies nothing — so that row times an extra 32 MB copy on the R side.
`blas_crossprod_chol` uses `crossprod`, R's flagged form, and is the cleaner of
the two.

### 2.2 presto is installed; glmGamPoi is not

Seurat routes Wilcoxon through presto when it is installed and through a much
slower internal loop when it is not. presto was not present on this machine, so
it was installed (`remotes::install_github("immunogenomics/presto")`, 1.0.0)
before any marker timing was taken. Every `find_all_markers` and `de_wilcox`
figure is Seurat on its fast path — the right comparison, and also where
truecell comes off worst.

The mirror image: **glmGamPoi is absent**, and Seurat says so ("could not find
glmGamPoi installed … falling back to native (slower) implementation"). The
SCTransform figures are Seurat *without* its accelerator; a previous round of
work in this repository found glmGamPoi not safely installable here, so it was
left alone. Treat that row as favourable to truecell by an unmeasured margin.

### 2.3 umap-learn gives up every thread the moment you set a seed

`run_umap(seed=42)` makes umap-learn set `random_state`, and umap-learn then
warns "n_jobs value 1 overridden to 1 by setting random_state" and runs
single-threaded. uwot, which Seurat uses, does not make that trade. The benches
measure both: `umap` is seeded, `umap_unseeded` is the same embedding without
one — 3.4x apart on PBMC 3k and 8.3–8.6x on the three larger sets.

It is no longer enough to turn a total around, but it is still the second
largest thing truecell can fix:

| | Truecell as measured | Truecell, seed dropped | Seurat |
|---|---:|---:|---:|
| pbmc3k | 13.5s | 6.6s | 5.0s |
| pbmc8k | 32.1s | 17.5s | 16.9s |
| ifnb | 43.7s | 30.0s | 18.9s |
| thp1 | 71.6s | 56.5s | 42.0s |

### Two smaller ones

* **Neighbours.** Seurat's default search is approximate (annoy); truecell's is
  exact. The tables compare against `nn.method = "rann"`, Seurat's exact option,
  and report `neighbours_annoy` separately so the default is visible too.
* **UMAP metric.** `metric="cosine"` on both arms, which is `RunUMAP`'s default.
  truecell's own default is euclidean, inherited from umap-learn.


---

## 3. Scaling: the standard workflow, 2.7k to 20.7k cells

| Dataset | Cells | Truecell | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS |
|---|---:|---:|---:|---:|---:|---:|
| `pbmc3k` | 2700 | 13.5s | 5.0s | 2.7x slower | 2938 MB | 1998 MB |
| `pbmc8k` | 8381 | 32.1s | 16.9s | 1.9x slower | 9253 MB | 3295 MB |
| `ifnb` | 13999 | 43.7s | 18.9s | 2.3x slower | 9420 MB | 3796 MB |
| `thp1` | 20729 | 71.6s | 42.0s | 1.7x slower | 10620 MB | 10919 MB |

The ratio is not flat any more. On the reference BLAS it sat at a steady 1.4x;
with both arms on Accelerate it runs 1.7x to 2.7x with no trend in *n* — the
spread now tracks how much of each dataset's total is UMAP and marker
detection, not dataset size. Neither tool is pulling away from the other as
*n* grows on this workload. The one place that is not true is Moran's I
(section 5), where Seurat is quadratic in cells and truecell is not.


## 4. Step by step

#### `pbmc3k_core`

| Step | Truecell | Truecell (1 thread) | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.58s | 0.51s | — | — | 124 MB | — |  |  |
| `read_counts` | 0.06s | 0.05s | 0.69s | **12.5x** | 188 MB | 590 MB | 2700 | 2700 |
| `create_object` | 0.02s | 0.02s | 0.29s | **19.2x** | 222 MB | 644 MB | 13714 | 13714 |
| `qc_metrics` | 0.00s | 0.00s | 0.01s | **4.4x** | 242 MB | 657 MB | 2.216642 | 2.216642 |
| `normalize` | 0.02s | 0.02s | 0.15s | **7.3x** | 242 MB | 696 MB | 4625433.29 | 4625433 |
| `hvg_vst` | 0.25s | 0.26s | 0.15s | 1.7x slower | 407 MB | 723 MB | 2000 | 2000 |
| `scale_hvg` | 0.10s | 0.10s | 0.17s | **1.7x** | 697 MB | 860 MB | 2000 | 2000 |
| `scale_all_genes` | 0.57s | 0.53s | 0.75s | **1.3x** | 2509 MB | 1648 MB | 13714 | 13714 |
| `rescale_hvg` | 0.08s | 0.09s | 0.16s | **1.9x** | 2535 MB | 1648 MB | 2000 | 2000 |
| `pca` | 0.22s | 0.28s | 0.32s | **1.5x** | 2535 MB | 1980 MB | 6.8875 | 6.8737 |
| `neighbours_exact` | 0.24s | 0.21s | 0.07s | 3.4x slower | 2564 MB | 1990 MB | 199616 | 198616 |
| `cluster_louvain` | 0.24s | 0.22s | 0.12s | 1.9x slower | 2627 MB | 1394 MB | 8 | 9 |
| `umap` | 10.00s | 10.49s | 2.41s | 4.2x slower | 2874 MB | 1397 MB | 2700 | 2700 |
| `umap_unseeded` | 2.83s | 2.89s | — | — | 2882 MB | — | 2700 |  |
| `find_all_markers` | 2.36s | 2.51s | 0.65s | 3.6x slower | 2938 MB | 1551 MB | 3118 | 3118 |
| `library_load` | — | — | 1.44s | — | — | 501 MB |  |  |
| `neighbours_annoy` | — | — | 0.38s | — | — | 1998 MB |  | 198484 |
| **shared pipeline** | **13.5s** | **14.2s** | **5.0s** | 2.7x slower | **2938 MB** | **1998 MB** |  |  |

Median of 3 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up and the arm-specific asides (`neighbours_annoy`, `rescale_hvg`, `scale_all_genes`, `umap_unseeded`).

#### `pbmc8k_core`

| Step | Truecell | Truecell (1 thread) | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.59s | 0.65s | — | — | 123 MB | — |  |  |
| `read_counts` | 0.17s | 0.17s | 2.57s | **14.9x** | 506 MB | 1173 MB | 8381 | 8381 |
| `create_object` | 0.07s | 0.07s | 0.59s | **8.3x** | 528 MB | 1185 MB | 18340 | 18340 |
| `qc_metrics` | 0.01s | 0.01s | 0.06s | **4.9x** | 619 MB | 1185 MB | 3.008575 | 3.008575 |
| `normalize` | 0.12s | 0.12s | 0.35s | **2.9x** | 619 MB | 1368 MB | 17833653.695 | 17833654 |
| `hvg_vst` | 0.86s | 0.87s | 0.34s | 2.5x slower | 1091 MB | 1374 MB | 2000 | 2000 |
| `scale_hvg` | 0.39s | 0.41s | 0.27s | 1.4x slower | 1819 MB | 1642 MB | 2000 | 2000 |
| `pca` | 0.51s | 0.71s | 1.11s | **2.2x** | 1820 MB | 2356 MB | 10.6853 | 10.6981 |
| `neighbours_exact` | 0.33s | 0.36s | 0.26s | 1.3x slower | 1930 MB | 2372 MB | 583657 | 585211 |
| `cluster_louvain` | 0.42s | 0.43s | 0.60s | **1.4x** | 1983 MB | 1643 MB | 9 | 12 |
| `umap` | 16.54s | 17.08s | 7.24s | 2.3x slower | 2355 MB | 1742 MB | 8381 | 8381 |
| `umap_unseeded` | 2.02s | 1.98s | — | — | 2399 MB | — | 8381 |  |
| `find_all_markers` | 12.71s | 12.31s | 3.57s | 3.6x slower | 9253 MB | 3295 MB | 7973 | 7974 |
| `library_load` | — | — | 1.49s | — | — | 493 MB |  |  |
| `neighbours_annoy` | — | — | 1.03s | — | — | 2372 MB |  | 585179 |
| **shared pipeline** | **32.1s** | **32.6s** | **16.9s** | 1.9x slower | **9253 MB** | **3295 MB** |  |  |

Median of 3 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up and the arm-specific asides (`neighbours_annoy`, `umap_unseeded`).

#### `ifnb_core`

| Step | Truecell | Truecell (1 thread) | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.59s | 0.66s | — | — | 123 MB | — |  |  |
| `read_counts` | 0.26s | 0.26s | 2.40s | **9.3x** | 474 MB | 1020 MB | 13999 | 13999 |
| `create_object` | 0.05s | 0.05s | 0.44s | **9.4x** | 475 MB | 1105 MB | 13915 | 13915 |
| `qc_metrics` | 0.00s | 0.00s | 0.04s | **12.7x** | 551 MB | 1105 MB | 0.0 | 0 |
| `normalize` | 0.10s | 0.10s | 0.29s | **3.0x** | 552 MB | 1256 MB | 21191453.753 | 21191454 |
| `hvg_vst` | 0.60s | 0.60s | 0.25s | 2.4x slower | 945 MB | 1261 MB | 2000 | 2000 |
| `scale_hvg` | 0.67s | 0.68s | 0.52s | 1.3x slower | 2463 MB | 1941 MB | 2000 | 2000 |
| `pca` | 0.86s | 1.27s | 1.79s | **2.1x** | 2464 MB | 3387 MB | 8.7871 | 8.783 |
| `neighbours_exact` | 0.41s | 0.44s | 0.53s | **1.3x** | 2515 MB | 3430 MB | 955357 | 954109 |
| `cluster_louvain` | 0.73s | 0.77s | 1.40s | **1.9x** | 2643 MB | 1747 MB | 15 | 16 |
| `umap` | 14.70s | 15.04s | 6.92s | 2.1x slower | 2994 MB | 1948 MB | 13999 | 13999 |
| `umap_unseeded` | 1.66s | 1.68s | — | — | 3081 MB | — | 13999 |  |
| `find_all_markers` | 25.37s | 23.23s | 4.30s | 5.9x slower | 9420 MB | 3796 MB | 5895 | 5907 |
| `library_load` | — | — | 1.47s | — | — | 485 MB |  |  |
| `neighbours_annoy` | — | — | 1.81s | — | — | 3611 MB |  | 953731 |
| **shared pipeline** | **43.7s** | **42.4s** | **18.9s** | 2.3x slower | **9420 MB** | **3796 MB** |  |  |

Median of 3 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up and the arm-specific asides (`neighbours_annoy`, `umap_unseeded`).

#### `thp1_core`

| Step | Truecell | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.65s | — | — | 123 MB | — |  |  |
| `read_counts` | 0.61s | 6.04s | **9.9x** | 1219 MB | 9646 MB | 20729 | 20729 |
| `create_object` | 0.40s | 1.61s | **4.0x** | 1691 MB | 10919 MB | 18381 | 18381 |
| `qc_metrics` | 0.05s | 0.31s | **5.7x** | 1755 MB | 10625 MB | 3.675326 | 3.675325 |
| `normalize` | 0.76s | 1.60s | **2.1x** | 2816 MB | 10741 MB | 66869178.271 | 66869178 |
| `hvg_vst` | 2.66s | 1.32s | 2.0x slower | 4084 MB | 10686 MB | 2000 | 2000 |
| `scale_hvg` | 1.07s | 0.57s | 1.9x slower | 5072 MB | 10744 MB | 2000 | 2000 |
| `pca` | 0.70s | 1.51s | **2.2x** | 5077 MB | 10708 MB | 9.2571 | 9.2638 |
| `neighbours_exact` | 0.65s | 1.66s | **2.5x** | 5248 MB | 10678 MB | 1217711 | 1217571 |
| `cluster_louvain` | 2.08s | 2.74s | **1.3x** | 5327 MB | 6543 MB | 7 | 8 |
| `umap` | 17.08s | 10.34s | 1.7x slower | 5700 MB | 6743 MB | 20729 | 20729 |
| `umap_unseeded` | 2.01s | — | — | 5404 MB | — | 20729 |  |
| `find_all_markers` | 45.52s | 14.31s | 3.2x slower | 10620 MB | 10915 MB | 9356 | 9365 |
| `library_load` | — | 1.59s | — | — | 485 MB |  |  |
| `neighbours_annoy` | — | 2.76s | — | — | 10685 MB |  | 1217807 |
| **shared pipeline** | **71.6s** | **42.0s** | 1.7x slower | **10620 MB** | **10919 MB** |  |  |

Median of 3 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up and the arm-specific asides (`neighbours_annoy`, `umap_unseeded`).


## 5. Named operations

#### `blas_probe`

| Step | Truecell | Truecell (1 thread) | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.46s | 0.47s | — | — | 123 MB | — |  |  |
| `blas_setup` | 0.01s | 0.02s | 0.18s | **12.2x** | 124 MB | 528 MB | 2000 | 2000 |
| `blas_gemm` | 0.01s | 0.02s | 0.13s | **9.3x** | 199 MB | 602 MB | 2002.361 | 2006.143 |
| `blas_svd` | 0.02s | 0.02s | 0.03s | **1.4x** | 199 MB | 605 MB | 66.575 | 67.223 |
| `blas_crossprod_chol` | 0.03s | 0.04s | 0.12s | **4.4x** | 199 MB | 660 MB | 62.883 | 63.745 |
| `library_load` | — | — | 1.38s | — | — | 486 MB |  |  |
| **shared pipeline** | **0.1s** | **0.1s** | **0.5s** | **6.0x** | **219 MB** | **660 MB** |  |  |

Median of 3 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up.

#### `pbmc3k_sctransform`

| Step | Truecell | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.51s | — | — | 123 MB | — |  |  |
| `read_counts` | 0.05s | 0.70s | **13.5x** | 138 MB | 604 MB | 2700 | 2700 |
| `create_object` | 0.01s | 0.30s | **19.7x** | 254 MB | 676 MB | 13714 | 13714 |
| `sctransform` | 39.71s | 52.94s | **1.3x** | 1643 MB | 3114 MB | 3000 | 3000 |
| `pca_on_sct` | 0.15s | 0.15s | ~equal | 1001 MB | 3116 MB | 14.3395 | 14.2162 |
| `library_load` | — | 1.48s | — | — | 487 MB |  |  |
| **shared pipeline** | **39.9s** | **54.1s** | **1.4x** | **1643 MB** | **3116 MB** |  |  |

Median of 2 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up.

#### `pbmc3k_de`

| Step | Truecell | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.50s | — | — | 119 MB | — |  |  |
| `prep` | 1.10s | 2.20s | **2.0x** | 818 MB | 985 MB | 8 | 9 |
| `de_wilcox` | 0.44s | 0.08s | 5.9x slower | 952 MB | 1026 MB | 2022 | 2022 |
| `de_t` | 0.52s | 2.12s | **4.1x** | 1032 MB | 1143 MB | 2022 | 2022 |
| `de_bimod` | 0.28s | 2.21s | **7.9x** | 1033 MB | 1062 MB | 2022 | 2022 |
| `de_LR` | 2.79s | 5.24s | **1.9x** | 1084 MB | 1119 MB | 2022 | 2022 |
| `de_negbinom` | 5.00s | 13.38s | **2.7x** | 1084 MB | 1207 MB | 2022 | 2022 |
| `de_roc` | 0.44s | 4.20s | **9.6x** | 1084 MB | 1249 MB | 2022 | 2022 |
| `de_MAST` | 3.22s | 10.34s | **3.2x** | 1086 MB | 2133 MB | 2022 | 2022 |
| `de_DESeq2` | _not available in Truecell_ | | | | | | |
| `library_load` | — | 1.48s | — | — | 495 MB |  |  |
| **shared pipeline** | **14.2s** | **53.6s** | **3.8x** | **1117 MB** | **2673 MB** |  |  |

Median of 2 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up.

#### `ifnb_integration`

| Step | Truecell | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.52s | — | — | 124 MB | — |  |  |
| `read_counts` | 0.26s | 2.43s | **9.3x** | 358 MB | 1063 MB | 13999 | 13999 |
| `prep_to_pca` | 1.87s | 2.52s | **1.3x** | 2574 MB | 3101 MB | 8.7871 | 8.783 |
| `harmony` | 1.06s | 2.22s | **2.1x** | 2580 MB | 3140 MB | 13999 | 13999 |
| `integrate_cca` | 24.76s | 6.87s | 3.6x slower | 5348 MB | 3798 MB | 13999 | 13999 |
| `integrate_rpca` | 9.04s | 8.41s | 1.1x slower | 4414 MB | 5329 MB | 13999 | 13999 |
| `library_load` | — | 1.51s | — | — | 501 MB |  |  |
| **shared pipeline** | **37.0s** | **22.4s** | 1.6x slower | **5348 MB** | **5329 MB** |  |  |

Median of 2 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up.

#### `xenium_spatial`

| Step | Truecell | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.52s | — | — | 124 MB | — |  |  |
| `read_xenium` | 0.12s | 2.31s | **20.0x** | 221 MB | 791 MB | 36602 | 36602 |
| `normalize` | 0.02s | 0.34s | **15.0x** | 255 MB | 812 MB | 36602 | 36602 |
| `morans_i_2k` | 0.03s | 2.55s | **86.9x** | 302 MB | 1042 MB | 248 | 248 |
| `morans_i_full` | 5.62s | — | — | 995 MB | — | 248 |  |
| `library_load` | — | 1.50s | — | — | 496 MB |  |  |
| **shared pipeline** | **0.2s** | **5.2s** | **31.0x** | **1000 MB** | **1042 MB** |  |  |

Median of 2 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up and the arm-specific asides (`morans_i_full`).


**SCTransform** — 39.5s against 52.9s, and half the memory (1.6 GB vs 3.1 GB).
Untouched by the BLAS swap, and read with the caveat in 2.2: Seurat is on its
non-glmGamPoi fallback here, so the true gap is smaller than 1.3x by an
unmeasured margin.

**Differential expression** — the cleanest comparison in the report, and
completely unaffected by the BLAS. Both arms test the same 2,022 genes on the
same two groups of cells (truecell writes the assignment, R adopts it) and both
return 2,022 rows. Truecell wins six of the seven per-cell tests: `roc` 9.6x,
`bimod` 7.9x, `t` 4.1x, `MAST` 3.2x, `negbinom` 2.7x, `LR` 1.9x. It loses
`wilcox` by 5.9x, which is presto's C++ Wilcoxon against
`scipy.stats.mannwhitneyu`.

`DESeq2` has no row because the two tools mean different things by the name:
Seurat's runs per cell, truecell's is a pseudobulk test requiring a replicate
column that PBMC 3k does not have. The DE vignette documents this as **seven**
per-cell tests, not eight, so it is a design divergence rather than a gap.

**Batch integration — the result the BLAS swap reversed.** Harmony still goes
to truecell at 2.1x (harmonypy against R's harmony, neither BLAS-bound). CCA
and RPCA go the other way: Seurat's CCA is **3.6x faster** (6.87s against
24.76s) and its RPCA 1.1x, where on the reference BLAS truecell led CCA by
3.6x. Same 13,999 cells, same PCA, same `k.weight` — the entire 13x that
Seurat's CCA gained came from the BLAS.

**Moran's I — the one asymptotic difference in the report.** On the identical
2,000-cell subset truecell is 87x faster. Then it keeps going: it computes the
full 36,602-cell slide in 5.6s in 1.0 GB, which Seurat cannot do at all.
`RunMoransI` builds `as.matrix(dist(pos))`, a dense n x n distance matrix —
10.7 GB at this n, before any statistic is computed. That is why the spatial
vignette subsets in the first place. Every other gap here is a constant factor;
this one is a wall, and no BLAS changes it.


## 6. The tutorial scripts, end to end

| Tutorial | Python script | R script | Faster by | Python peak RSS | R peak RSS |
|---|---:|---:|---:|---:|---:|
| `pbmc3k` | 16.3s | 8.4s | 1.9x slower | 2883 MB | 1877 MB |
| `sctransform` | 61.4s | 65.2s | **1.1x** | 7386 MB | 3448 MB |
| `de` | 103.8s | 340.7s | **3.3x** | 3551 MB | 5705 MB |
| `dimreduc` | 8.6s | 22.4s | **2.6x** | 1013 MB | 1250 MB |
| `objects` | 2.9s | 4.9s | **1.7x** | 1396 MB | 2020 MB |
| `integration` | 49.8s | 52.6s | **1.1x** | 5443 MB | 6321 MB |
| `cellcycle` | 4.1s | 13.7s | **3.3x** | 3572 MB | 8980 MB |
| `svf` | 1.8s | 6.9s | **3.9x** | 379 MB | 1031 MB |
| `visium` | 3.2s | 7.0s | **2.2x** | 1554 MB | 2058 MB |
| `lazy` | 129.7s | 9.4s | 13.8x slower | 3303 MB | 4099 MB |


These are the tutorial scripts as they ship, not the benches — a different
question, and a noisier one. Each script also prints validation, writes CSVs or
draws figures, and the two sides do not do equal amounts of that. Two rows need
saying out loud:

* **`lazy` is not a 13x loss.** The Python script runs all eight DE tests on an
  out-of-core layer; `FindMarkers` on a BPCells `IterableMatrix` supports
  `wilcox` alone, so the R script attempts eight and completes one. That is a
  capability difference being reported by a stopwatch — the R script finishes
  sooner because there is less it can do.
* **`dimreduc`** was 15.3x on the reference BLAS and is the row the swap hit
  hardest, because JackStraw is 300 permuted PCAs.


---

## 7. What each tool is good at

**Truecell is faster at getting data in.** Reading 10x matrices 9–15x, building
the object 4–20x, normalisation 2–7x, QC metrics 4–13x. None of this is BLAS;
it is sparse I/O and object construction. On the THP-1 dense TSV truecell reads
20,729 cells in 0.60s and 687 MB against Seurat's 6.04s and 9.5 GB, and that
single step sets R's peak for the whole run.

**PCA is still truecell's, but by 1.5–2.2x rather than 12–15x.** Randomized SVD
against irlba, both now on the same Accelerate. That residual gap is the real
one; the rest was the BLAS.

**Truecell wins six of the seven shared DE tests** by 1.9–9.6x, and Moran's I
by 87x with no *n* limit.

**Seurat wins seven operations, and three of them decide the totals.**

| Operation | Gap | Where it comes from |
|---|---|---|
| `find_all_markers` | 3.2–5.9x, and 2.5x the memory | `find_markers` densifies the whole matrix before filtering it |
| `umap` (seeded) | 1.7–4.2x | umap-learn drops to one thread under a `random_state` |
| `de_wilcox` | 5.9x | presto's C++ Wilcoxon vs `scipy.stats.mannwhitneyu` |
| `integrate_cca` | 3.6x | reversed by the BLAS swap; see 2.1 |
| `hvg_vst` | 1.7–2.5x | |
| `scale_hvg` | 1.3–1.9x on 3 of 4 datasets | |
| `neighbours_exact` | 1.3–3.4x on 2 of 4 datasets | |

The bottom three are small in absolute terms — under 1.4s each on every dataset
here — and none changes a total on its own.

**Memory crosses over twice.** Through `hvg_vst` truecell is lighter on every
dataset, by 1.3x on PBMC 3k and 14x on THP-1. From `scale_hvg` onward on the 10x
sets R pulls ahead, partly because its GC returns pages to the OS — visible as
the drop from 3.3 GB to 1.6 GB across `cluster_louvain` on ifnb, which CPython's
arenas never do. Then `find_all_markers` moves truecell from 2.0 GB to 9.4 GB in
one step and settles the process peak: 9.3–9.4 GB against Seurat's 3.3–3.8 GB on
the same cell assignment.

The exception is THP-1, where Seurat is heavier for the whole run — its
`read_counts` peaks at 9.5 GB turning a dense TSV into a sparse matrix, and R
never comes back below it.

## 8. What to do about it

**1. `find_markers` should filter before it densifies.** Now the single largest
thing on this list. `truecell/markers.py` turns the entire gene x cell matrix
into two dense float64 arrays and only then computes the pct and logFC masks
that reduce it to the genes actually tested. Both masks can be computed on the
sparse matrix — the `data` layer is log1p-normalised and `expm1(0) == 0`, so
`expm1` preserves the sparsity pattern exactly — after which only the surviving
genes need to be dense. There is also a redundant `.astype(float)` on an array
`toarray()` already returned as float64. Largest memory step in the pipeline and
the largest timing loss.

**2. Decide what `run_umap`'s seed should cost.** Passing one costs 3.5–8.8x
because umap-learn silently single-threads. Either document it at the call site
or expose the choice, so reproducibility is something a caller opts into
knowingly rather than pays for by default.

**3. Look at `integrate_layers(method="cca")` again.** At 24.8s against Seurat's
6.9s on the same input it is now the third-largest gap, and the reference-BLAS
sweep hid it completely.

**4. Consider a fast path for Wilcoxon.** `scipy.stats.mannwhitneyu` is 5.9x
off presto on the same 2,022 genes, and Wilcoxon is the default test — so this
lands on `find_all_markers` as well, compounding with item 1.

**5. For anyone running the R side: install presto, and glmGamPoi if you can.**
presto is what makes Seurat's marker detection competitive. glmGamPoi does the
same for SCTransform, and Seurat is on its fallback path here.

**6. Keep R on Accelerate.** Already done on this machine. It changed no result
anywhere in the suite and took up to 12.9x off individual steps; if you rebuild
or reinstall R, redo the symlink in 2.1.
