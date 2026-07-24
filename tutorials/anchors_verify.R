#!/usr/bin/env Rscript
# Seurat reference for the anchor-internals tutorial (anchors_vignette.md).
#
# Runs FindIntegrationAnchors + IntegrateData for BOTH reductions on exactly the
# cells and anchor features the Python side wrote, and records the anchors
# themselves — the pairs and their scores — rather than a downstream clustering.
#
# nn.method = "rann" gives exact nearest neighbours. Seurat's default is annoy,
# which is approximate: re-running the same data with annoy moves ~0.3% of the
# anchors, and that noise would otherwise be indistinguishable from a real
# disagreement with shanuz.
#
# Run order:
#   Rscript tutorials/export_seuratdata.R ifnb    # one-time counts export
#   python  tutorials/anchors_tutorial.py         # writes cells + anchor features
#   Rscript tutorials/anchors_verify.R
#   python  tutorials/anchors_tutorial.py --report
#
# Needs: Seurat, Matrix, jsonlite.
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(jsonlite)
})
set.seed(42)
options(future.globals.maxSize = 3 * 1024^3)

.args   <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
HERE <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG  <- file.path(HERE, "figures_anchors")
DATA <- Sys.getenv("IFNB_DATA", path.expand("~/.shanuz_data/ifnb"))

DIMS <- 1:30

for (f in c("cells_CTRL.txt", "cells_STIM.txt", "anchor_features.txt")) {
  if (!file.exists(file.path(FIG, f)))
    stop(f, " not found in ", FIG,
         "\nRun `python tutorials/anchors_tutorial.py` first (it writes them).")
}
if (!file.exists(file.path(DATA, "matrix.mtx.gz")))
  stop("ifnb counts not found at ", DATA,
       "\nRun `Rscript tutorials/export_seuratdata.R ifnb` first.")

# ---- 1. The same cells and features the Python side used --------------------
cells_ctrl <- readLines(file.path(FIG, "cells_CTRL.txt"))
cells_stim <- readLines(file.path(FIG, "cells_STIM.txt"))
feats      <- readLines(file.path(FIG, "anchor_features.txt"))

counts <- Read10X(DATA)
meta   <- read.csv(file.path(DATA, "metadata.csv"), row.names = 1, check.names = FALSE)
meta   <- meta[colnames(counts), , drop = FALSE]
obj    <- CreateSeuratObject(counts = counts, min.cells = 3, meta.data = meta)
obj    <- subset(obj, cells = c(cells_ctrl, cells_stim))
cat(sprintf("%d genes x %d cells | CTRL %d / STIM %d | %d anchor features\n",
            nrow(obj), ncol(obj), length(cells_ctrl), length(cells_stim), length(feats)))

# Reference first, query second — IntegrateData corrects the query onto the
# reference and copies the reference through untouched.
lst <- list(CTRL = subset(obj, cells = cells_ctrl),
            STIM = subset(obj, cells = cells_stim))
lst <- lapply(lst, function(o) {
  o <- NormalizeData(o, verbose = FALSE)
  VariableFeatures(o) <- feats
  o <- ScaleData(o, features = feats, verbose = FALSE)
  RunPCA(o, features = feats, npcs = max(DIMS), verbose = FALSE)   # rpca needs it
})

# ---- 2. Both reductions ------------------------------------------------------
# The anchor table Seurat hands to IntegrateData is symmetrised and offset into
# a merged index space; grabbing it as FindWeights receives it gives the pairs
# in the reference/query cell order this comparison needs.
.cap <- new.env(); assign(".cap", .cap, envir = globalenv())
trace(Seurat:::FindWeights, exit = quote({
  assign("anchors", Seurat:::GetIntegrationData(returnValue(), "integrated", "anchors"), envir = .cap)
  assign("nbrs",    Seurat:::GetIntegrationData(returnValue(), "integrated", "neighbors"), envir = .cap)
}), print = FALSE)

for (reduction in c("cca", "rpca")) {
  cat("\n", reduction, " ...\n", sep = "")
  a <- FindIntegrationAnchors(object.list = lst, anchor.features = feats,
                              reduction = reduction, dims = DIMS,
                              nn.method = "rann", verbose = FALSE)
  integrated <- IntegrateData(anchorset = a, dims = DIMS, verbose = FALSE)

  AN <- .cap$anchors; NB <- .cap$nbrs
  pairs <- data.frame(cell1 = NB$cells1[AN[, "cell1"]],
                      cell2 = NB$cells2[AN[, "cell2"]],
                      score = AN[, "score"], stringsAsFactors = FALSE)
  write.csv(pairs, file.path(FIG, paste0("r_anchors_", reduction, ".csv")), row.names = FALSE)

  # Only the query half moves; the reference is copied through, so summarising
  # over every cell would average the correction against a block of exact zeros.
  IDATA <- GetAssayData(integrated, assay = "integrated", layer = "data")
  RAW   <- GetAssayData(lst[["STIM"]], layer = "data")[feats, cells_stim]
  DELTA <- as.matrix(IDATA[feats, cells_stim] - RAW)
  stats <- list(n_anchors = nrow(pairs), score_mean = mean(pairs$score),
                delta_mean_abs = mean(abs(DELTA)), delta_max_abs = max(abs(DELTA)),
                delta_frac_nonzero = mean(DELTA != 0),
                seurat_version = as.character(packageVersion("Seurat")))
  write(toJSON(stats, digits = 22, auto_unbox = TRUE, pretty = TRUE),
        file.path(FIG, paste0("r_summary_", reduction, ".json")))
  cat(sprintf("  %d anchors  score mean %.4f  correction mean|d| %.6f  frac nz %.4f\n",
              nrow(pairs), mean(pairs$score), mean(abs(DELTA)), mean(DELTA != 0)))
}
untrace(Seurat:::FindWeights)

cat("\nDONE — wrote r_anchors_*.csv and r_summary_*.json to", FIG, "\n")
cat("Now run: python tutorials/anchors_tutorial.py --report\n")
