# Truecell vs R Seurat — Performance

The companion to the accuracy comparison in [`../README.md`](../README.md).
That work asked whether the two tools produce the same answers. This one asks
what each answer costs.

**Machine.** Apple M4 Pro, 12 cores, 25.8 GB, macOS 26.5.2. Measured 2 Aug 2026.

**Versions.** R 4.6.1 · Seurat 5.5.1 · Matrix 1.7.5 · presto 1.0.0 ·
harmony 2.0.5 · uwot 0.2.4 · irlba 2.3.7 · RANN 2.6.2 · Rfast2 ·
data.table 1.18.4 — Python 3.12.13 · truecell 1.0.0 · numpy 2.4.6 ·
scipy 1.18.0 · umap-learn · scikit-learn · harmonypy.

Reproduce with `bash tutorials/benchmark/sweep.sh` (about an hour).

---

## Summary

Across 68 like-for-like step comparisons Truecell is **faster in 50, slower in
17, and within 5% in 1**. The 17 losses are only six distinct operations, and
two of them — seeded UMAP and marker detection — are enough to lose the
standard end-to-end workflow by 1.4x on their own. Both have identifiable,
fixable causes.

Every gap larger than 10x in this report is a build or configuration difference
rather than a difference between the two projects. Read section 2 before the
tables.

| | Truecell | Seurat |
|---|---|---|
| Standard workflow, 2.7k–20.7k cells | | **1.4–1.7x faster** |
| …with `run_umap`'s seed dropped | **1.1–1.4x faster** (3 of 4 datasets) | |
| Reading counts, building the object | **4–20x faster** | |
| PCA | **12–15x faster** (mostly BLAS, see 2.1) | |
| Differential expression, 6 of the 7 shared tests | **1.9–9.5x faster** | |
| Wilcoxon — `de_wilcox` and `find_all_markers` | | **3.0–5.8x faster** (presto) |
| Seeded UMAP | | **2.2–4.0x faster** |
| Batch integration (CCA / RPCA / Harmony) | **1.3–3.6x faster** | |
| SCTransform | **1.4x faster**, half the memory | |
| Moran's I | **102x faster**, and the only one without an *n* limit | |
| Peak memory, standard workflow | | **2.5x lighter** (all of it one step) |
| Peak memory, everything else | **1.2–14x lighter** | |

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
late step's peak partly reports what earlier steps left behind. Per-step
figures are peaks *during* the step; the pipeline row is the process peak,
which is the number that decides whether a machine can run the workload.

**The first repeat of every pair is discarded.** umap-learn's numba kernels
compile on first use and R loads packages lazily; neither is what the benchmark
is asking about. Three timed repeats follow, and the median is reported.

**Every step records an anchor** — cells kept, clusters found, genes tested,
markers returned — printed in the tables beside the timings. A speed comparison
is only worth reading if both sides did the same work, and the anchors are how
you can check rather than take it on trust.

**The instrument is not free.** Sampling costs one `ps` call every 50 ms, which
is a few percent of one core out of twelve. It is applied identically to both
arms, so the comparisons hold, but treat the absolute seconds as very slightly
inflated on both sides.

**The truecell arm runs first, on purpose.** It writes the cell-to-cluster
assignment and the Xenium cell subset that the R arm reads back. Without that,
two things would have gone wrong: the tools do not always land on the same
number of clusters (9 against Seurat's 12 on PBMC 8k), and one-vs-rest marker
detection costs one test per cluster, so `find_all_markers` would have been
timing a clustering difference. With it, both arms return the same marker
table — 3118 rows against 3118 on PBMC 3k.

---

## 2. Three things that decide most of these numbers

Read these before the tables. Two of the three largest gaps in this report
belong to build and configuration choices rather than to either project.

### 2.1 This R links the reference BLAS; numpy links Accelerate

`blas_probe` is a control: dense linear algebra on identical inputs, touching
neither Seurat nor truecell.

| Step | Truecell | Truecell (1 thread) | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.46s | 0.47s | — | — | 123 MB | — |  |  |
| `blas_setup` | 0.02s | 0.02s | 0.18s | **11.8x** | 123 MB | 523 MB | 2000 | 2000 |
| `blas_gemm` | 0.02s | 0.02s | 2.10s | **126.8x** | 200 MB | 584 MB | 2002.361 | 2006.143 |
| `blas_svd` | 0.02s | 0.02s | 0.31s | **15.3x** | 200 MB | 615 MB | 66.575 | 67.223 |
| `blas_crossprod_chol` | 0.03s | 0.04s | 3.94s | **133.2x** | 200 MB | 676 MB | 62.883 | 63.745 |
| `library_load` | — | — | 1.38s | — | — | 492 MB |  |  |
| **shared pipeline** | **0.1s** | **0.1s** | **6.5s** | **80.1x** | **305 MB** | **676 MB** |  |  |

The R installed here loads `libRblas.0.dylib`, the unoptimised reference BLAS
that ships with R. numpy loads Apple's Accelerate. **This is not a core-count
effect** — 12 cores cannot produce 127x, and pinning every documented
thread-limit variable to 1 barely moves the truecell column. At minimum a
factor of ten of it is per-core efficiency, which is the quality of the two
BLAS builds. (The one-thread column is a weaker control than it looks:
Accelerate may or may not honour `VECLIB_MAXIMUM_THREADS`. It does not matter
to the conclusion — the arithmetic above holds either way.)

Every `pca` and `scale_data` figure in this report inherits that difference, so
read the 12-15x PCA gaps as mostly BLAS and only partly implementation.

**The fix is one line**, and `libRblas.vecLib.dylib` is already sitting in the
same directory:

```bash
cd /Library/Frameworks/R.framework/Resources/lib && ln -sf libRblas.vecLib.dylib libRblas.dylib
```

That modifies your R installation, so this report measures R as installed. Say
the word and the sweep can be re-run against Accelerate for a second column.

### 2.2 presto was installed for this benchmark

Seurat routes Wilcoxon through the presto package when it is installed and
through a much slower internal loop when it is not. presto was not present on
this machine, so it was installed (`remotes::install_github(
"immunogenomics/presto")`, version 1.0.0) before any marker timing was taken.
Every `find_all_markers` and `de_wilcox` figure below is Seurat on its fast
path — which is the right comparison, and also the reason those are the
numbers where truecell comes off worst.

The mirror image: **glmGamPoi is not installed**, and Seurat says so —
"`vst.flavor` is set to 'v2' but could not find glmGamPoi installed ... falling
back to native (slower) implementation." The SCTransform figures are therefore
Seurat *without* its accelerator. A previous round of work in this repository
found glmGamPoi not safely installable in this environment, so that was left
alone; treat the SCTransform row as favourable to truecell by an unmeasured
margin.

### 2.3 umap-learn gives up every thread the moment you set a seed

`run_umap(seed=42)` makes umap-learn set `random_state`, and umap-learn then
warns "n_jobs value 1 overridden to 1 by setting random_state" and runs
single-threaded. uwot, which Seurat uses, does not make that trade. The benches
measure both: `umap` is seeded, `umap_unseeded` is the same embedding without
one.

The difference is 3.5x on PBMC 3k and 8.4-8.8x on all three larger sets, and it
is most of why truecell's standard workflow comes out behind overall. Swapping
the seeded step for the unseeded one turns three of the four datasets around:

| | Truecell as measured | Truecell, seed dropped | Seurat |
|---|---:|---:|---:|
| pbmc3k | 13.1s | 6.1s | 7.8s |
| pbmc8k | 31.6s | 17.2s | 23.2s |
| ifnb | 41.7s | 28.6s | 30.3s |
| thp1 | 70.2s | 55.5s | 48.8s |

thp1 stays behind because marker detection dominates there, which is section
2.2's story rather than this one.

### Two smaller ones

* **Neighbours.** Seurat's default search is approximate (annoy); truecell's is
  exact. The tables compare against `nn.method = "rann"`, Seurat's exact
  option, and report `neighbours_annoy` separately so the default is visible
  too.
* **UMAP metric.** `metric="cosine"` on both arms, which is `RunUMAP`'s
  default. truecell's own default is euclidean, inherited from umap-learn.


---

## 3. Scaling: the standard workflow, 2.7k to 20.7k cells

| Dataset | Cells | Truecell | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS |
|---|---:|---:|---:|---:|---:|---:|
| `pbmc3k` | 2700 | 13.1s | 7.8s | 1.7x slower | 2931 MB | 2001 MB |
| `pbmc8k` | 8381 | 31.6s | 23.2s | 1.4x slower | 9151 MB | 3243 MB |
| `ifnb` | 13999 | 41.7s | 30.3s | 1.4x slower | 9174 MB | 3669 MB |
| `thp1` | 20729 | 70.2s | 48.8s | 1.4x slower | 9634 MB | 9987 MB |

The ratio is flat: 1.4x from 8,381 cells to 20,729, with the 1.7x at 2,700 an
artifact of how much of that small total is UMAP. Neither tool is pulling away
from the other as *n* grows on this workload — the constant factors differ, the
scaling does not. The one place that is not true is Moran's I (section 5),
where Seurat is quadratic in cells and truecell is not.


## 4. Step by step

#### `pbmc3k_core`

| Step | Truecell | Truecell (1 thread) | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.47s | 0.49s | — | — | 122 MB | — |  |  |
| `read_counts` | 0.05s | 0.05s | 0.69s | **14.0x** | 218 MB | 583 MB | 2700 | 2700 |
| `create_object` | 0.01s | 0.01s | 0.29s | **20.4x** | 218 MB | 640 MB | 13714 | 13714 |
| `qc_metrics` | 0.00s | 0.00s | 0.01s | **4.1x** | 221 MB | 642 MB | 2.216642 | 2.216642 |
| `normalize` | 0.02s | 0.02s | 0.15s | **7.0x** | 292 MB | 692 MB | 4625433.29 | 4625433 |
| `hvg_vst` | 0.25s | 0.25s | 0.40s | **1.6x** | 407 MB | 718 MB | 2000 | 2000 |
| `scale_hvg` | 0.09s | 0.09s | 0.16s | **1.8x** | 592 MB | 884 MB | 2000 | 2000 |
| `scale_all_genes` | 0.54s | 0.54s | 0.75s | **1.4x** | 2437 MB | 1610 MB | 13714 | 13714 |
| `rescale_hvg` | 0.08s | 0.08s | 0.16s | **2.0x** | 2532 MB | 1643 MB | 2000 | 2000 |
| `pca` | 0.21s | 0.27s | 2.83s | **13.7x** | 2533 MB | 1978 MB | 6.8875 | 6.8737 |
| `neighbours_exact` | 0.20s | 0.20s | 0.07s | 2.7x slower | 2561 MB | 1994 MB | 199616 | 198616 |
| `cluster_louvain` | 0.21s | 0.22s | 0.12s | 1.7x slower | 2614 MB | 1391 MB | 8 | 9 |
| `umap` | 9.76s | 10.04s | 2.45s | 4.0x slower | 2867 MB | 1394 MB | 2700 | 2700 |
| `umap_unseeded` | 2.77s | 2.81s | — | — | 2875 MB | — | 2700 |  |
| `find_all_markers` | 2.31s | 2.33s | 0.66s | 3.5x slower | 2930 MB | 1556 MB | 3118 | 3118 |
| `library_load` | — | — | 1.49s | — | — | 493 MB |  |  |
| `neighbours_annoy` | — | — | 0.39s | — | — | 2001 MB |  | 198484 |
| **shared pipeline** | **13.1s** | **13.5s** | **7.8s** | 1.7x slower | **2931 MB** | **2001 MB** |  |  |

Median of 3 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up and the arm-specific asides (`neighbours_annoy`, `rescale_hvg`, `scale_all_genes`, `umap_unseeded`).

#### `pbmc8k_core`

| Step | Truecell | Truecell (1 thread) | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.59s | 0.63s | — | — | 123 MB | — |  |  |
| `read_counts` | 0.17s | 0.17s | 2.64s | **15.5x** | 524 MB | 1174 MB | 8381 | 8381 |
| `create_object` | 0.07s | 0.07s | 0.59s | **8.2x** | 529 MB | 1182 MB | 18340 | 18340 |
| `qc_metrics` | 0.01s | 0.01s | 0.06s | **5.0x** | 620 MB | 1182 MB | 3.008575 | 3.008575 |
| `normalize` | 0.12s | 0.12s | 0.34s | **2.9x** | 621 MB | 1363 MB | 17833653.695 | 17833654 |
| `hvg_vst` | 0.85s | 0.85s | 0.82s | ~equal | 1092 MB | 1390 MB | 2000 | 2000 |
| `scale_hvg` | 0.39s | 0.39s | 0.28s | 1.4x slower | 1820 MB | 1666 MB | 2000 | 2000 |
| `pca` | 0.49s | 0.69s | 6.67s | **13.6x** | 1821 MB | 2344 MB | 10.6853 | 10.6981 |
| `neighbours_exact` | 0.33s | 0.34s | 0.26s | 1.2x slower | 1870 MB | 2368 MB | 583657 | 585211 |
| `cluster_louvain` | 0.41s | 0.42s | 0.60s | **1.5x** | 1997 MB | 1639 MB | 9 | 12 |
| `umap` | 16.36s | 16.26s | 7.28s | 2.2x slower | 2356 MB | 1693 MB | 8381 | 8381 |
| `umap_unseeded` | 1.94s | 1.92s | — | — | 2394 MB | — | 8381 |  |
| `find_all_markers` | 12.44s | 11.88s | 3.63s | 3.4x slower | 9151 MB | 3243 MB | 7973 | 7974 |
| `library_load` | — | — | 1.48s | — | — | 496 MB |  |  |
| `neighbours_annoy` | — | — | 1.05s | — | — | 2423 MB |  | 585179 |
| **shared pipeline** | **31.6s** | **31.2s** | **23.2s** | 1.4x slower | **9151 MB** | **3243 MB** |  |  |

Median of 3 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up and the arm-specific asides (`neighbours_annoy`, `umap_unseeded`).

#### `ifnb_core`

| Step | Truecell | Truecell (1 thread) | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.65s | 0.64s | — | — | 123 MB | — |  |  |
| `read_counts` | 0.27s | 0.26s | 2.38s | **8.9x** | 471 MB | 1093 MB | 13999 | 13999 |
| `create_object` | 0.05s | 0.05s | 0.44s | **9.5x** | 475 MB | 1108 MB | 13915 | 13915 |
| `qc_metrics` | 0.00s | 0.00s | 0.04s | **12.8x** | 549 MB | 1108 MB | 0.0 | 0 |
| `normalize` | 0.10s | 0.10s | 0.29s | **3.0x** | 550 MB | 1258 MB | 21191453.753 | 21191454 |
| `hvg_vst` | 0.60s | 0.59s | 0.51s | 1.2x slower | 945 MB | 1261 MB | 2000 | 2000 |
| `scale_hvg` | 0.67s | 0.67s | 0.52s | 1.3x slower | 2461 MB | 1956 MB | 2000 | 2000 |
| `pca` | 0.88s | 1.26s | 13.07s | **14.8x** | 2463 MB | 3245 MB | 8.7871 | 8.783 |
| `neighbours_exact` | 0.43s | 0.42s | 0.53s | **1.2x** | 2514 MB | 3283 MB | 955357 | 954109 |
| `cluster_louvain` | 0.76s | 0.74s | 1.39s | **1.8x** | 2644 MB | 1599 MB | 15 | 16 |
| `umap` | 14.76s | 14.65s | 6.83s | 2.2x slower | 2989 MB | 1800 MB | 13999 | 13999 |
| `umap_unseeded` | 1.68s | 1.65s | — | — | 3077 MB | — | 13999 |  |
| `find_all_markers` | 23.13s | 23.68s | 4.33s | 5.3x slower | 9174 MB | 3669 MB | 5895 | 5907 |
| `library_load` | — | — | 1.45s | — | — | 493 MB |  |  |
| `neighbours_annoy` | — | — | 1.80s | — | — | 3463 MB |  | 953731 |
| **shared pipeline** | **41.7s** | **42.4s** | **30.3s** | 1.4x slower | **9174 MB** | **3669 MB** |  |  |

Median of 3 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up and the arm-specific asides (`neighbours_annoy`, `umap_unseeded`).

#### `thp1_core`

| Step | Truecell | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.63s | — | — | 122 MB | — |  |  |
| `read_counts` | 0.60s | 5.82s | **9.7x** | 687 MB | 9542 MB | 20729 | 20729 |
| `create_object` | 0.39s | 1.61s | **4.1x** | 1674 MB | 9980 MB | 18381 | 18381 |
| `qc_metrics` | 0.05s | 0.27s | **5.1x** | 1755 MB | 9510 MB | 3.675326 | 3.675325 |
| `normalize` | 0.75s | 1.57s | **2.1x** | 2816 MB | 9782 MB | 66869178.271 | 66869178 |
| `hvg_vst` | 2.63s | 1.74s | 1.5x slower | 4084 MB | 9664 MB | 2000 | 2000 |
| `scale_hvg` | 1.05s | 0.57s | 1.8x slower | 5072 MB | 9612 MB | 2000 | 2000 |
| `pca` | 0.70s | 8.33s | **12.0x** | 5077 MB | 9911 MB | 9.2571 | 9.2638 |
| `neighbours_exact` | 0.65s | 1.61s | **2.5x** | 5245 MB | 9361 MB | 1217711 | 1217571 |
| `cluster_louvain` | 2.03s | 2.63s | **1.3x** | 5363 MB | 9361 MB | 7 | 8 |
| `umap` | 16.70s | 9.96s | 1.7x slower | 5701 MB | 8994 MB | 20729 | 20729 |
| `umap_unseeded` | 1.97s | — | — | 5527 MB | — | 20729 |  |
| `find_all_markers` | 44.62s | 14.74s | 3.0x slower | 9634 MB | 9618 MB | 9356 | 9365 |
| `library_load` | — | 1.53s | — | — | 489 MB |  |  |
| `neighbours_annoy` | — | 2.60s | — | — | 9411 MB |  | 1217807 |
| **shared pipeline** | **70.2s** | **48.8s** | 1.4x slower | **9634 MB** | **9987 MB** |  |  |

Median of 3 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up and the arm-specific asides (`neighbours_annoy`, `umap_unseeded`).


## 5. Named operations

#### `blas_probe`

| Step | Truecell | Truecell (1 thread) | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.46s | 0.47s | — | — | 123 MB | — |  |  |
| `blas_setup` | 0.02s | 0.02s | 0.18s | **11.8x** | 123 MB | 523 MB | 2000 | 2000 |
| `blas_gemm` | 0.02s | 0.02s | 2.10s | **126.8x** | 200 MB | 584 MB | 2002.361 | 2006.143 |
| `blas_svd` | 0.02s | 0.02s | 0.31s | **15.3x** | 200 MB | 615 MB | 66.575 | 67.223 |
| `blas_crossprod_chol` | 0.03s | 0.04s | 3.94s | **133.2x** | 200 MB | 676 MB | 62.883 | 63.745 |
| `library_load` | — | — | 1.38s | — | — | 492 MB |  |  |
| **shared pipeline** | **0.1s** | **0.1s** | **6.5s** | **80.1x** | **305 MB** | **676 MB** |  |  |

Median of 3 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up.

#### `pbmc3k_sctransform`

| Step | Truecell | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.50s | — | — | 119 MB | — |  |  |
| `read_counts` | 0.05s | 0.69s | **13.4x** | 188 MB | 608 MB | 2700 | 2700 |
| `create_object` | 0.01s | 0.30s | **19.9x** | 266 MB | 676 MB | 13714 | 13714 |
| `sctransform` | 39.32s | 54.03s | **1.4x** | 1642 MB | 3104 MB | 3000 | 3000 |
| `pca_on_sct` | 0.15s | 1.25s | **8.3x** | 990 MB | 3105 MB | 14.3395 | 14.2162 |
| `library_load` | — | 1.48s | — | — | 487 MB |  |  |
| **shared pipeline** | **39.5s** | **56.3s** | **1.4x** | **1642 MB** | **3105 MB** |  |  |

Median of 2 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up.

#### `pbmc3k_de`

| Step | Truecell | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.53s | — | — | 124 MB | — |  |  |
| `prep` | 1.13s | 5.16s | **4.6x** | 817 MB | 918 MB | 8 | 9 |
| `de_wilcox` | 0.45s | 0.08s | 5.8x slower | 950 MB | 1000 MB | 2022 | 2022 |
| `de_t` | 0.52s | 2.21s | **4.3x** | 1030 MB | 1154 MB | 2022 | 2022 |
| `de_bimod` | 0.28s | 2.22s | **8.0x** | 1031 MB | 1066 MB | 2022 | 2022 |
| `de_LR` | 2.87s | 5.30s | **1.9x** | 1081 MB | 1142 MB | 2022 | 2022 |
| `de_negbinom` | 5.12s | 13.54s | **2.6x** | 1082 MB | 1185 MB | 2022 | 2022 |
| `de_roc` | 0.44s | 4.16s | **9.5x** | 1082 MB | 1235 MB | 2022 | 2022 |
| `de_MAST` | 3.27s | 10.62s | **3.2x** | 1084 MB | 2094 MB | 2022 | 2022 |
| `de_DESeq2` | _not available in Truecell_ | | | | | | |
| `library_load` | — | 1.49s | — | — | 485 MB |  |  |
| **shared pipeline** | **14.5s** | **57.7s** | **4.0x** | **1115 MB** | **2494 MB** |  |  |

Median of 2 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up.

#### `ifnb_integration`

| Step | Truecell | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.62s | — | — | 122 MB | — |  |  |
| `read_counts` | 0.27s | 2.38s | **8.8x** | 473 MB | 1018 MB | 13999 | 13999 |
| `prep_to_pca` | 1.86s | 5.73s | **3.1x** | 2573 MB | 3078 MB | 8.7871 | 8.783 |
| `harmony` | 1.05s | 2.54s | **2.4x** | 2580 MB | 3129 MB | 13999 | 13999 |
| `integrate_cca` | 24.26s | 88.41s | **3.6x** | 4848 MB | 3610 MB | 13999 | 13999 |
| `integrate_rpca` | 9.00s | 11.91s | **1.3x** | 4364 MB | 5206 MB | 13999 | 13999 |
| `library_load` | — | 1.54s | — | — | 493 MB |  |  |
| **shared pipeline** | **36.4s** | **111.0s** | **3.0x** | **4848 MB** | **5206 MB** |  |  |

Median of 2 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up.

#### `xenium_spatial`

| Step | Truecell | Seurat | Faster by | Truecell peak RSS | Seurat peak RSS | Truecell result | Seurat result |
|---|---:|---:|---:|---:|---:|---:|---:|
| `import` | 0.51s | — | — | 120 MB | — |  |  |
| `read_xenium` | 0.12s | 2.28s | **19.3x** | 229 MB | 788 MB | 36602 | 36602 |
| `normalize` | 0.02s | 0.34s | **14.9x** | 266 MB | 804 MB | 36602 | 36602 |
| `morans_i_2k` | 0.03s | 3.14s | **102.4x** | 379 MB | 1052 MB | 248 | 248 |
| `morans_i_full` | 5.54s | — | — | 1001 MB | — | 248 |  |
| `library_load` | — | 1.50s | — | — | 492 MB |  |  |
| **shared pipeline** | **0.2s** | **5.8s** | **33.6x** | **1001 MB** | **1052 MB** |  |  |

Median of 2 timed repeats per arm, warm-up discarded. *Shared pipeline* excludes interpreter start-up and the arm-specific asides (`morans_i_full`).


**SCTransform** — 39.3s against 54.0s, and half the memory (1.6 GB vs 3.1 GB).
Read with the caveat in 2.2: Seurat is on its non-glmGamPoi fallback path here,
so the true gap is smaller than 1.4x by an unmeasured margin.

**Differential expression** — the cleanest comparison in the report. Both arms
test the same 2,022 genes on the same two groups of cells (truecell writes the
assignment, R adopts it), and both return 2,022 rows. Truecell wins six of the
seven per-cell tests, by 1.9x to 9.5x. It loses `wilcox` by 5.8x, which is
presto's C++ Wilcoxon against truecell's `scipy.stats.mannwhitneyu`.

`DESeq2` has no row because the two tools mean different things by the name:
Seurat's runs per cell, truecell's is a pseudobulk test requiring a replicate
column that PBMC 3k does not have. The DE vignette documents this as **seven**
per-cell tests, not eight, so this is a design divergence rather than a gap.

**Batch integration** — the largest clean win, and the one least contaminated
by the BLAS: Harmony 2.4x, RPCA 1.3x, CCA **3.6x** (24.3s against 88.4s), all
on the same 13,999 cells from the same PCA with the same `k.weight`.

**Moran's I — the one asymptotic difference in the report.** On the identical
2,000-cell subset truecell is 102x faster. Then it keeps going: it computes the
full 36,602-cell slide in 5.5s in 1.0 GB, which Seurat cannot do at all.
`RunMoransI` builds `as.matrix(dist(pos))`, a dense n x n distance matrix —
10.7 GB at this n, before any statistic is computed. That is why the spatial
vignette subsets in the first place. Every other gap here is a constant factor;
this one is a wall.


## 6. The tutorial scripts, end to end

| Tutorial | Python script | R script | Faster by | Python peak RSS | R peak RSS |
|---|---:|---:|---:|---:|---:|
| `pbmc3k` | 16.5s | 11.5s | 1.4x slower | 2882 MB | 1859 MB |
| `sctransform` | 61.9s | 71.8s | **1.2x** | 4993 MB | 3454 MB |
| `de` | 103.4s | 342.6s | **3.3x** | 3553 MB | 5312 MB |
| `dimreduc` | 8.6s | 131.5s | **15.3x** | 1010 MB | 1201 MB |
| `objects` | 2.8s | 5.9s | **2.1x** | 1399 MB | 2021 MB |
| `integration` | 48.8s | 147.8s | **3.0x** | 5212 MB | 5769 MB |
| `cellcycle` | 4.1s | 14.9s | **3.6x** | 3574 MB | 7789 MB |
| `svf` | 1.7s | 7.5s | **4.3x** | 391 MB | 1033 MB |
| `visium` | 3.2s | 8.2s | **2.5x** | 1551 MB | 2069 MB |
| `lazy` | 132.1s | 10.0s | 13.2x slower | 3304 MB | 4101 MB |


These are the tutorial scripts as they ship, not the benches — a different
question, and a noisier one. Each script also prints validation, writes CSVs or
draws figures, and the two sides do not do equal amounts of that. Two rows need
saying out loud:

* **`lazy` is not a 13x loss.** The Python script runs all eight DE tests on an
  out-of-core layer; `FindMarkers` on a BPCells `IterableMatrix` supports
  `wilcox` alone, so the R script attempts eight and completes one. That is a
  capability difference being reported by a stopwatch — the R script finishes
  sooner because there is less it can do.
* **`dimreduc` at 15.3x** is real and is JackStraw: 300 permutation replicates
  of a PCA, which is section 2.1's BLAS gap multiplied 300 times.


---

## 7. What each tool is good at

**Truecell is faster at getting data in and at dense linear algebra.**
Reading 10x matrices is 9-15x faster, building the object 4-20x, normalisation
2-7x, and PCA 12-15x — though most of that last one is the BLAS rather than the
code. Louvain is 1.3-1.8x. On the THP-1 dense TSV, truecell reads 20,729 cells
in 0.60s and 687 MB against Seurat's 5.82s and 9.5 GB, and that single step
sets R's peak for the whole run.

**Truecell is faster at six of the seven differential-expression tests both
tools can run**, on the same 2,022 genes and the same two groups of cells:
`roc` 9.5x, `bimod` 8.0x, `t` 4.3x, `MAST` 3.2x, `negbinom` 2.6x, `LR` 1.9x.
The seventh, `wilcox`, goes the other way — see the table below.

**Seurat wins six operations, and two of them decide the totals.** In order of
how much they matter:

| Operation | Gap | Where it comes from |
|---|---|---|
| `find_all_markers` | 3.0-5.3x, and 2.5x the memory | `find_markers` densifies the whole matrix before filtering it |
| `umap` (seeded) | 1.7-4.0x | umap-learn drops to one thread under a `random_state` |
| `de_wilcox` | 5.8x | presto's C++ Wilcoxon vs `scipy.stats.mannwhitneyu` |
| `scale_hvg` | 1.3-1.8x on 3 of 4 datasets | |
| `hvg_vst` | 1.2-1.5x on 2 of 4 datasets | |
| `neighbours_exact` | 1.2-2.7x on 2 of 4 datasets | |

The last three are small in absolute terms — under a second each on every
dataset here — and none of them changes a total. The first two are not
inherent to either language.

**Memory crosses over twice.** Through `hvg_vst` truecell is lighter on every
dataset, by 1.3x on PBMC 3k and 14x on THP-1. From `scale_hvg` onward on the
10x sets R pulls ahead, partly because its GC returns pages to the OS —
visible as the drop from 3.3 GB to 1.6 GB across `cluster_louvain` on ifnb,
which CPython's arenas never do. Then `find_all_markers` moves truecell from
2.0 GB to 9.2 GB in one step and settles the process peak: 9.2 GB against
Seurat's 3.2-3.7 GB on the same cell assignment.

The exception is THP-1, where Seurat is heavier for the entire run — its
`read_counts` peaks at 9.5 GB turning a dense TSV into a sparse matrix, and R
never comes back below that.

## 8. What to do about it

**1. Point R at Accelerate.** The largest single win available to the Seurat
arm, one line, no code change:

```bash
cd /Library/Frameworks/R.framework/Resources/lib && ln -sf libRblas.vecLib.dylib libRblas.dylib
```

**2. `find_markers` should filter before it densifies.** `truecell/markers.py`
turns the entire gene x cell matrix into two dense float64 arrays and only then
computes the pct and logFC masks that reduce it to the genes actually tested.
Both masks can be computed on the sparse matrix — the `data` layer is
log1p-normalised, `expm1(0) == 0`, so `expm1` preserves the sparsity pattern
exactly — after which only the surviving genes need to be dense. There is also
a redundant `.astype(float)` on an array `toarray()` already returned as
float64. This is the largest memory step in the pipeline and the largest
timing loss.

**3. Decide what `run_umap`'s seed should cost.** Right now passing one costs
3.5-8.5x because umap-learn silently single-threads. Either document it at the
call site or expose the choice, so that reproducibility is something a caller
opts into knowingly rather than something they pay for by default.

**4. For anyone running the R side: install presto, and glmGamPoi if you can.**
presto is what makes Seurat's marker detection competitive; without it that
column would look very different. glmGamPoi does the same for SCTransform, and
Seurat is currently on its fallback path here.

