#!/usr/bin/env Rscript
# The R Seurat arm of the performance benchmark — the mirror of
# bench_truecell.py, step for step, parameter for parameter.
#
# Run through run_benchmarks.py, which measures this process from outside:
#   python tutorials/benchmark/run_benchmarks.py run --bench pbmc3k_core
#
# Two places where the arms genuinely differ get an *extra* step rather than a
# quiet substitution:
#
#   * `neighbours_annoy` — Seurat's default neighbour search is approximate
#     (annoy); truecell's is exact. `neighbours_exact` (nn.method = "rann") is
#     the like-for-like number and `neighbours_annoy` is what a Seurat user
#     actually waits for. Reporting one without the other would be a choice
#     about which tool to flatter.
#
#   * `find_all_markers` — Seurat dispatches Wilcoxon to presto when it is
#     installed and to a much slower internal loop when it is not. Which one
#     ran is printed to the log so the report can say which was measured.
#
# `step(log, "name", { ... })` times the block. The block is a promise, so it
# evaluates in *this* frame and plain `<-` inside it updates the pipeline
# object as usual. Its value, if a single number, becomes the step's anchor.

.T_START <- as.numeric(Sys.time())

.args <- commandArgs(trailingOnly = TRUE)
.all  <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .all[grep("^--file=", .all)])
HERE <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
source(file.path(HERE, "steps.R"))

suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
})
.T_LIBS <- as.numeric(Sys.time())

get_arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), .args)
  if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
BENCH <- get_arg("bench")
STEPS <- get_arg("steps")
if (is.null(BENCH) || is.null(STEPS)) stop("need --bench and --steps")

set.seed(0)
DATA_ROOT <- path.expand("~/.truecell_data")

# `IntegrateLayers` ships work to future workers and refuses to move more than
# 500 MiB by default; ifnb's anchor step is 2 GiB. The tutorials' own verify
# scripts raise it to the same 3 GiB, so the benchmark is not measuring a
# configuration the tutorials do not use.
options(future.globals.maxSize = 3 * 1024^3)

# Pipeline parameters, identical to the Python arm.
N_HVG <- 2000
N_PCS <- 50
DIMS <- 1:10
K_PARAM <- 20
RESOLUTION <- 0.5

HAS_PRESTO <- requireNamespace("presto", quietly = TRUE)

# ---------------------------------------------------------------------------
# Datasets — the same bytes the Python arm reads
# ---------------------------------------------------------------------------

load_counts <- function(bench) {
  key <- sub("_.*$", "", bench)
  if (key == "pbmc3k") {
    list(counts = Read10X(file.path(DATA_ROOT, "pbmc3k",
                                    "filtered_gene_bc_matrices", "hg19")),
         meta = NULL)
  } else if (key == "pbmc8k") {
    list(counts = Read10X(file.path(DATA_ROOT, "pbmc8k",
                                    "filtered_gene_bc_matrices", "GRCh38")),
         meta = NULL)
  } else if (key == "ifnb") {
    d <- file.path(DATA_ROOT, "ifnb")
    counts <- Read10X(d)
    meta <- read.csv(file.path(d, "metadata.csv"), row.names = 1,
                     check.names = FALSE)
    list(counts = counts, meta = meta[colnames(counts), , drop = FALSE])
  } else if (key == "thp1") {
    # A dense 20k-cell TSV, not a 10x folder. data.table::fread is how the
    # tutorial reads it and the only way R gets through it in reasonable time;
    # the Python arm reads the same file. The cost of the format itself lands
    # in `read_counts` on both sides.
    d <- file.path(DATA_ROOT, "thp1_eccite")
    dt <- data.table::fread(cmd = paste("gzcat",
        shQuote(file.path(d, "GSM4633614_ECCITE_cDNA_counts.tsv.gz"))),
        showProgress = FALSE)
    genes <- dt[[1]]
    mat <- as.matrix(dt[, -1]); rownames(mat) <- genes
    rm(dt); invisible(gc())
    counts <- as(Matrix(mat, sparse = TRUE), "CsparseMatrix")
    rm(mat); invisible(gc())
    meta <- read.delim(gzfile(file.path(d, "GSE153056_ECCITE_metadata.tsv.gz")),
                       row.names = 1, check.names = FALSE)
    common <- intersect(colnames(counts), rownames(meta))
    list(counts = counts[, common], meta = meta[common, , drop = FALSE])
  } else {
    stop("unknown dataset for bench ", bench)
  }
}

# ---------------------------------------------------------------------------
# The standard workflow
# ---------------------------------------------------------------------------

# Adopt the cell-to-cluster assignment the Python arm wrote. Untimed: it is
# measurement scaffolding, not pipeline work. Without it the two arms would run
# one-vs-rest markers over different numbers of clusters — 12 against 9 on PBMC
# 8k — and the timing would be reporting a clustering difference.
adopt_python_idents <- function(bench, obj) {
  csv <- file.path(HERE, "results", paste0(bench, "_idents.csv"))
  if (!file.exists(csv)) {
    cat("WARNING: no Python idents file; running on R's own clusters\n")
    return(obj)
  }
  tbl <- read.csv(csv, colClasses = "character")
  shared <- intersect(tbl$cell, colnames(obj))
  obj <- subset(obj, cells = shared)
  Idents(obj) <- factor(tbl$ident[match(colnames(obj), tbl$cell)])
  cat(sprintf("adopted Python's assignment: %d cells, %d idents\n",
              length(shared), length(levels(Idents(obj)))))
  obj
}

run_core <- function(bench, log) {
  step(log, "read_counts", {
    x <- load_counts(bench)
    ncol(x$counts)
  })
  # No meta.data even where the dataset ships some: the Python arm does not
  # attach it either, and a metadata frame neither side computes on would be an
  # asymmetry in the measurement for no gain in what is measured.
  step(log, "create_object", {
    obj <- CreateSeuratObject(x$counts, project = bench, min.cells = 3,
                              min.features = 200)
    nrow(obj)
  })
  rm(x); invisible(gc())

  step(log, "qc_metrics", {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    round(mean(obj[["percent.mt"]][, 1]), 6)
  })
  step(log, "normalize", {
    obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                         scale.factor = 10000, verbose = FALSE)
    round(sum(GetAssayData(obj, layer = "data")), 3)
  })
  step(log, "hvg_vst", {
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                nfeatures = N_HVG, verbose = FALSE)
    length(VariableFeatures(obj))
  })
  hvg <- VariableFeatures(obj)
  step(log, "scale_hvg", {
    obj <- ScaleData(obj, features = hvg, verbose = FALSE)
    length(hvg)
  })

  # Scaling every gene is what the PBMC 3k vignette does and what dominates its
  # memory. Only measured on the small dataset: at 20k cells the dense result is
  # tens of gigabytes in both languages, which measures the machine.
  if (startsWith(bench, "pbmc3k")) {
    step(log, "scale_all_genes", {
      obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
      nrow(obj)
    })
    step(log, "rescale_hvg", {
      obj <- ScaleData(obj, features = hvg, verbose = FALSE)
      length(hvg)
    })
  }

  step(log, "pca", {
    obj <- RunPCA(obj, features = hvg, npcs = N_PCS, verbose = FALSE)
    round(Stdev(obj, reduction = "pca")[1], 4)
  })
  step(log, "neighbours_exact", {
    obj <- FindNeighbors(obj, dims = DIMS, k.param = K_PARAM,
                         nn.method = "rann", verbose = FALSE)
    length(obj@graphs[[paste0(DefaultAssay(obj), "_snn")]]@x)
  })
  # Seurat's default search, measured and discarded: the clusters below come
  # from the exact graph, so both arms cluster the same neighbour table.
  step(log, "neighbours_annoy", {
    g <- FindNeighbors(Embeddings(obj, "pca")[, DIMS], k.param = K_PARAM,
                       nn.method = "annoy", verbose = FALSE)
    length(g$snn@x)
  })
  step(log, "cluster_louvain", {
    obj <- FindClusters(obj, resolution = RESOLUTION, algorithm = 1,
                        random.seed = 0, verbose = FALSE)
    length(levels(Idents(obj)))
  })
  step(log, "umap", {
    obj <- RunUMAP(obj, dims = DIMS, verbose = FALSE, seed.use = 42)
    nrow(Embeddings(obj, "umap"))
  })
  obj <- adopt_python_idents(bench, obj)
  n_markers <- step(log, "find_all_markers", {
    m <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                        logfc.threshold = 0.25, verbose = FALSE)
    nrow(m)
  })
  cat(sprintf("%s: %d cells, %d clusters, %d markers (presto: %s)\n",
              bench, ncol(obj), length(levels(Idents(obj))), n_markers,
              HAS_PRESTO))
}

# ---------------------------------------------------------------------------
# Named heavy operations
# ---------------------------------------------------------------------------

run_sctransform <- function(bench, log) {
  step(log, "read_counts", { x <- load_counts(bench); ncol(x$counts) })
  step(log, "create_object", {
    obj <- CreateSeuratObject(x$counts, project = bench, min.cells = 3,
                              min.features = 200)
    nrow(obj)
  })
  rm(x); invisible(gc())
  step(log, "sctransform", {
    obj <- SCTransform(obj, variable.features.n = 3000, verbose = FALSE)
    length(VariableFeatures(obj))
  })
  step(log, "pca_on_sct", {
    obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
    round(Stdev(obj, reduction = "pca")[1], 4)
  })
}

run_integration <- function(bench, log) {
  n_pcs <- 30
  step(log, "read_counts", { x <- load_counts(bench); ncol(x$counts) })
  step(log, "prep_to_pca", {
    obj <- CreateSeuratObject(x$counts, project = bench, min.cells = 3,
                              meta.data = x$meta)
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                nfeatures = N_HVG, verbose = FALSE)
    obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
    obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = n_pcs,
                  verbose = FALSE)
    round(Stdev(obj, reduction = "pca")[1], 4)
  })
  rm(x); invisible(gc())

  step(log, "harmony", {
    obj <- harmony::RunHarmony(obj, group.by.vars = "stim",
                               reduction.use = "pca",
                               reduction.save = "harmony", verbose = FALSE)
    nrow(Embeddings(obj, "harmony"))
  })
  k_weight <- min(100, min(table(obj$stim)) %/% 2)
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$stim)
  step(log, "integrate_cca", {
    obj <- IntegrateLayers(obj, method = CCAIntegration, orig.reduction = "pca",
                           new.reduction = "cca", k.weight = k_weight,
                           verbose = FALSE)
    nrow(Embeddings(obj, "cca"))
  })
  step(log, "integrate_rpca", {
    obj <- IntegrateLayers(obj, method = RPCAIntegration, orig.reduction = "pca",
                           new.reduction = "rpca", k.weight = k_weight,
                           verbose = FALSE)
    nrow(Embeddings(obj, "rpca"))
  })
}

run_de <- function(bench, log) {
  step(log, "prep", {
    x <- load_counts(bench)
    obj <- CreateSeuratObject(x$counts, project = bench, min.cells = 3,
                              min.features = 200)
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                nfeatures = N_HVG, verbose = FALSE)
    obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
    obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = N_PCS,
                  verbose = FALSE)
    obj <- FindNeighbors(obj, dims = DIMS, k.param = K_PARAM,
                         nn.method = "rann", verbose = FALSE)
    obj <- FindClusters(obj, resolution = RESOLUTION, algorithm = 1,
                        random.seed = 0, verbose = FALSE)
    length(levels(Idents(obj)))
  })

  obj <- adopt_python_idents(bench, obj)

  for (test in c("wilcox", "t", "bimod", "LR", "negbinom", "roc",
                 "MAST", "DESeq2")) {
    step(log, paste0("de_", test), {
      res <- try(suppressWarnings(suppressMessages(
        FindMarkers(obj, ident.1 = "0", ident.2 = "1", test.use = test,
                    logfc.threshold = 0.25, min.pct = 0.1, verbose = FALSE))),
        silent = TRUE)
      if (inherits(res, "try-error")) {
        cat(sprintf("  de_%s unavailable: %s", test, as.character(res)))
        -1
      } else nrow(res)
    })
  }
}

run_spatial <- function(bench, log) {
  d <- file.path(DATA_ROOT, "xenium_mouse_brain")
  step(log, "read_xenium", {
    # molecule.coordinates = FALSE because the cached download is the analysis
    # subset — cell_feature_matrix/ and cells.csv.gz, no transcripts.parquet —
    # and LoadXenium errors rather than skipping it. The svf tutorial's verify
    # script loads it the same way; truecell's load_xenium never reads the
    # molecules, so this is what makes the two arms read the same bytes.
    obj <- LoadXenium(d, fov = "fov", assay = "Xenium",
                      molecule.coordinates = FALSE)
    ncol(obj)
  })
  step(log, "normalize", {
    obj <- NormalizeData(obj, scale.factor = 10000, verbose = FALSE)
    ncol(obj)
  })
  # The matrix form, which is what the svf tutorial's verify script uses. The
  # Seurat-object method goes through FindSpatiallyVariableFeatures.StdAssay
  # and errors on a v5 assay before it computes anything.
  #
  # Only the 2,000-cell subset. RunMoransI does `as.matrix(dist(pos))`, so the
  # full 36,602-cell slide asks for a 10.7 GB dense distance matrix up front;
  # the truecell arm runs both sizes and the report carries the full-slide
  # number on one side only, which is the honest way to show an O(n^2) wall.
  cells_file <- file.path(HERE, "results", "xenium_cells.txt")
  if (!file.exists(cells_file))
    stop("run the truecell arm first; it writes the shared cell subset")
  cells <- readLines(cells_file)
  step(log, "morans_i_2k", {
    dat <- as.matrix(GetAssayData(obj, assay = "Xenium",
                                  layer = "data")[, cells])
    pos <- GetTissueCoordinates(obj[["fov"]])
    rownames(pos) <- pos$cell
    svf <- FindSpatiallyVariableFeatures(dat,
             spatial.location = pos[cells, c("x", "y")],
             selection.method = "moransi", verbose = FALSE)
    nrow(svf)
  })
}

# Dense linear algebra on identical inputs, touching neither library. Not a
# single-cell benchmark: a control. R here links the reference BLAS it ships
# with, numpy links Accelerate. Every PCA, scaling and SVD number in this suite
# inherits that difference, and without measuring it directly there is no way to
# say how much of a gap belongs to the two projects and how much belongs to the
# two BLAS builds underneath them.
run_blas <- function(bench, log) {
  n <- 2000L; k <- 500L
  step(log, "blas_setup", {
    set.seed(0)
    a <- matrix(rnorm(n * n), n, n)
    n
  })
  # `a %*% t(a)` materialises the transpose first; numpy's `a @ a.T` passes a
  # transpose flag to dgemm and copies nothing. So this step times R's extra
  # 32 MB copy alongside the multiply, and overstates the pure-BLAS gap.
  # `blas_crossprod_chol` below uses `crossprod`, R's flagged form, and is the
  # cleaner of the two comparisons.
  step(log, "blas_gemm", {
    g <- a %*% t(a)
    round(g[1, 1], 3)
  })
  step(log, "blas_svd", {
    s <- svd(a[, seq_len(k)], nu = 0, nv = 0)$d
    round(s[1], 3)
  })
  step(log, "blas_crossprod_chol", {
    cc <- chol(crossprod(a) + n * diag(n))
    round(cc[1, 1], 3)
  })
}

# ---------------------------------------------------------------------------

BENCHES <- list(
  blas_probe = run_blas,
  pbmc3k_core = run_core, pbmc8k_core = run_core,
  ifnb_core = run_core, thp1_core = run_core,
  pbmc3k_sctransform = run_sctransform,
  ifnb_integration = run_integration,
  pbmc3k_de = run_de,
  xenium_spatial = run_spatial
)
if (!BENCH %in% names(BENCHES)) stop("unknown bench: ", BENCH)

log <- step_log_open(STEPS)
step_mark(log, "library_load", .T_LIBS - .T_START)
on.exit(step_log_close(log))
BENCHES[[BENCH]](BENCH, log)
