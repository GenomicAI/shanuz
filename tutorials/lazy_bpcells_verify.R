#!/usr/bin/env Rscript
# R/Seurat reference for the out-of-core tutorial: the same pbmc3k pipeline run
# twice, once on a dgCMatrix and once on a BPCells on-disk matrix.
#
# Writes tutorials/figures_lazy/r_anchors.json plus the per-gene tables, which
# `python tutorials/lazy_bpcells_tutorial.py --report` compares against.
#
# Requires: Seurat, BPCells, jsonlite, digest
#   BPCells is NOT on CRAN and needs libhdf5 at build time:
#     brew install hdf5 pkg-config
#     remotes::install_github("bnprks/BPCells/r")
#   It is inert for ordinary input -- every reference to it inside Seurat is
#   gated on `inherits(x, "IterableMatrix")`, never on `requireNamespace`, so
#   unlike glmGamPoi its presence cannot re-flavour a dgCMatrix run. Verified
#   by fingerprinting ten Seurat references before and after installing it:
#   none moved.
#
# Cells come from figures_lazy/cells.txt, written by the Python side, so both
# tools analyse the same matrix.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
for (pkg in c("BPCells", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("%s is required.", pkg), call. = FALSE)
  }
}
suppressPackageStartupMessages(library(BPCells))
set.seed(42)

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
TUT <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG <- file.path(TUT, "figures_lazy")
dir.create(FIG, showWarnings = FALSE)
DATA <- path.expand("~/.shanuz_data/pbmc3k")
STORE <- file.path(FIG, "bpcells_store")

HEAD <- 20L
anchors <- list()
# No rounding, and `digits = NA` on write_json below: at 10 decimal places a
# mean of 0.0034 carries a 3e-8 relative floor, which is larger than several of
# the differences under test and would be mistaken for one.
add <- function(name, value) anchors[[name]] <<- I(as.numeric(value))
addc <- function(name, value) anchors[[name]] <<- I(as.character(value))

cat("Building pbmc3k...\n")
raw <- Read10X(file.path(DATA, "filtered_gene_bc_matrices/hg19"))
obj <- CreateSeuratObject(raw, project = "lazy", min.cells = 3, min.features = 200)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cells_file <- file.path(FIG, "cells.txt")
if (file.exists(cells_file)) {
  keep <- readLines(cells_file)
  obj <- subset(obj, cells = intersect(keep, colnames(obj)))
}
counts <- GetAssayData(obj, layer = "counts")
cat(sprintf("  %d features x %d cells, nnz=%d\n",
            nrow(counts), ncol(counts), length(counts@x)))
writeLines(rownames(counts), file.path(FIG, "r_features.txt"))

# ---------------------------------------------------------------------------
# The store, and what bitpacking is worth
# ---------------------------------------------------------------------------
store_size <- function(d) {
  sum(file.info(list.files(d, recursive = TRUE, full.names = TRUE))$size)
}
unlink(STORE, recursive = TRUE)
bp <- write_matrix_dir(mat = counts, dir = STORE, compress = TRUE)

# Counts arrive from Seurat as doubles. BPCells' bitpacking is built for
# integers and says so; the honest comparison reports both, because the
# integer figure is the one its documentation quotes.
int_dir <- file.path(FIG, "bpcells_store_int")
unlink(int_dir, recursive = TRUE)
bp_int <- write_matrix_dir(convert_matrix_type(counts, "uint32_t"), int_dir,
                           compress = TRUE)

dgc_bytes <- length(counts@x) * 8 + length(counts@i) * 4 + length(counts@p) * 4
add("store.nnz", length(counts@x))
add("store.nrow", nrow(counts))
add("store.ncol", ncol(counts))
add("store.dgc_bytes", dgc_bytes)
add("store.bpcells_double_bytes", store_size(STORE))
add("store.bpcells_uint32_bytes", store_size(int_dir))
addc("store.storage_order", storage_order(bp))
cat(sprintf("  store: %.2f MB double / %.2f MB uint32 (dgCMatrix arrays %.2f MB)\n",
            store_size(STORE) / 1e6, store_size(int_dir) / 1e6, dgc_bytes / 1e6))

# ---------------------------------------------------------------------------
# The same pipeline twice: in memory, and out of core
# ---------------------------------------------------------------------------
run_pipeline <- function(mat, label) {
  o <- CreateSeuratObject(mat, project = label, min.cells = 0, min.features = 0)
  o <- NormalizeData(o, verbose = FALSE)
  o <- FindVariableFeatures(o, selection.method = "vst", nfeatures = 2000,
                            verbose = FALSE)
  o <- ScaleData(o, features = rownames(o), verbose = FALSE)
  o
}

cat("Running the in-memory pipeline...\n")
o_mem <- run_pipeline(counts, "mem")
cat("Running the out-of-core pipeline...\n")
o_bp <- run_pipeline(bp, "bp")

dense <- function(x) if (inherits(x, "IterableMatrix")) as.matrix(x) else as.matrix(x)
sparse_x <- function(o) {
  d <- GetAssayData(o, layer = "data")
  if (inherits(d, "IterableMatrix")) d <- as(d, "dgCMatrix")
  d
}

d_mem <- sparse_x(o_mem); d_bp <- sparse_x(o_bp)
h_mem <- HVFInfo(o_mem);  h_bp <- HVFInfo(o_bp)

# Both runs are emitted. BPCells computes in single precision, so Seurat's
# out-of-core numbers sit ~1e-6 from its own in-memory ones; shanuz stays in
# float64 throughout, so `mem.*` is the series it should be judged against and
# `bp.*` is the series that shows what the precision change costs.
s_bp <- as.matrix(GetAssayData(o_bp, layer = "scale.data"))
s_mem <- as.matrix(GetAssayData(o_mem, layer = "scale.data"))

add("calcn.ncount_head", head(o_bp$nCount_RNA, HEAD))
add("calcn.nfeature_head", head(o_bp$nFeature_RNA, HEAD))

for (tag in c("bp", "mem")) {
  o <- if (tag == "bp") o_bp else o_mem
  dd <- if (tag == "bp") d_bp else d_mem
  hh <- if (tag == "bp") h_bp else h_mem
  ss <- if (tag == "bp") s_bp else s_mem
  add(paste0(tag, ".normalize_head"), head(dd@x, HEAD))
  add(paste0(tag, ".normalize_sum"), sum(dd@x))
  add(paste0(tag, ".vst_mean_head"), head(hh$mean, HEAD))
  add(paste0(tag, ".vst_variance_head"), head(hh$variance, HEAD))
  add(paste0(tag, ".vst_var_std_head"), head(hh$variance.standardized, HEAD))
  addc(paste0(tag, ".vst_selected_head"), head(VariableFeatures(o), HEAD))
  add(paste0(tag, ".scale_head"), head(as.numeric(ss[1, ]), HEAD))
  add(paste0(tag, ".scale_abs_sum"), sum(abs(ss)))
}

# --- the control: Seurat against itself ------------------------------------
add("selfcheck.normalize_max_diff", max(abs(d_mem@x - d_bp@x)))
add("selfcheck.normalize_identical", as.numeric(identical(d_mem@x, d_bp@x)))
add("selfcheck.vst_mean_max_diff", max(abs(h_mem$mean - h_bp$mean)))
add("selfcheck.vst_var_std_max_diff",
    max(abs(h_mem$variance.standardized - h_bp$variance.standardized)))
add("selfcheck.hvg_overlap",
    length(intersect(VariableFeatures(o_mem), VariableFeatures(o_bp))))
common <- intersect(rownames(s_mem), rownames(s_bp))
add("selfcheck.scale_max_diff", max(abs(s_mem[common, ] - s_bp[common, ])))

cat(sprintf("\n  Seurat vs itself: normalize %.3e, var.std %.3e, HVG %d/2000\n",
            max(abs(d_mem@x - d_bp@x)),
            max(abs(h_mem$variance.standardized - h_bp$variance.standardized)),
            length(intersect(VariableFeatures(o_mem), VariableFeatures(o_bp)))))

# ---------------------------------------------------------------------------
# FindMarkers: what the out-of-core path supports
# ---------------------------------------------------------------------------
groups <- rep(c("a", "b"), length.out = ncol(counts))
Idents(o_bp) <- Idents(o_mem) <- groups
supported <- character()
for (tt in c("wilcox", "t", "bimod", "LR", "negbinom", "roc", "MAST", "DESeq2")) {
  r <- try(suppressWarnings(FindMarkers(o_bp, ident.1 = "a", ident.2 = "b",
                                        test.use = tt, logfc.threshold = 0,
                                        min.pct = 0, verbose = FALSE)),
           silent = TRUE)
  if (!inherits(r, "try-error")) supported <- c(supported, tt)
}
addc("markers.supported_tests", supported)
cat(sprintf("  FindMarkers on an IterableMatrix supports: %s\n",
            paste(supported, collapse = ", ")))

mk <- suppressWarnings(FindMarkers(o_bp, ident.1 = "a", ident.2 = "b",
                                   test.use = "wilcox", logfc.threshold = 0,
                                   min.pct = 0, verbose = FALSE))
write.csv(mk, file.path(FIG, "r_markers_wilcox.csv"))
addc("markers.wilcox_top10", head(rownames(mk), 10))
add("markers.n_genes", nrow(mk))

write.csv(h_bp, file.path(FIG, "r_hvf_info.csv"))
writeLines(VariableFeatures(o_bp), file.path(FIG, "r_variable_features.txt"))
jsonlite::write_json(anchors, file.path(FIG, "r_anchors.json"),
                     auto_unbox = TRUE, digits = NA, pretty = TRUE)
cat("\nWrote", file.path(FIG, "r_anchors.json"), "\n")
cat("Compare with: python tutorials/lazy_bpcells_tutorial.py --report\n")
