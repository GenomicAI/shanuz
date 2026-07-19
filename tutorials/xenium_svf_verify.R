#!/usr/bin/env Rscript
# Reference anchors for the spatial statistics / container tutorial.
#
# Writes tutorials/figures_svf/r_anchors.json, which
# `python tutorials/xenium_svf_tutorial.py --report` compares against.
#
# Requires: Seurat, SeuratObject, Matrix, jsonlite, digest, Rfast2
#   install.packages(c("jsonlite", "digest", "Rfast2"))
# Rfast2 is what Seurat's FindSpatiallyVariableFeatures(selection.method =
# "moransi") calls through RunMoransI; without it Seurat silently falls back to
# ape::Moran.I, which is a different estimator, so the comparison would be
# against something other than Seurat's documented default.
#
# The Moran's I section runs on the 2,000-cell subset the Python side wrote to
# figures_svf/cells.txt. That is not tidiness — RunMoransI does
# `as.matrix(dist(pos))`, so the full 36,602-cell slide needs a 10.7 GB
# allocation before it computes anything.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(jsonlite)
  library(digest)
})

for (pkg in c("Rfast2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("%s is required; install.packages('%s')", pkg, pkg), call. = FALSE)
  }
}

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
TUT <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG <- file.path(TUT, "figures_svf")
dir.create(FIG, showWarnings = FALSE)

DATA <- path.expand("~/.shanuz_data/xenium_mouse_brain")
ASSAY <- "Xenium"
TOP_N <- 10

# Order-sensitive fingerprint, matching the Python `digest()` helper:
# md5 of the names joined by newlines, first 12 hex characters.
name_digest <- function(x) {
  substr(digest(paste(as.character(x), collapse = "\n"),
                algo = "md5", serialize = FALSE), 1, 12)
}

name_anchor <- function(x) {
  x <- as.character(x)
  list(n = length(x), digest = name_digest(x),
       head = utils::head(x, 3), tail = utils::tail(x, 3))
}

cat("Loading Xenium slide...\n")
# molecule.coordinates = FALSE: the public cache ships the cell-feature matrix
# and cells.csv.gz but not transcripts.parquet, and LoadXenium errors rather
# than skipping the molecules on its own.
obj <- LoadXenium(DATA, fov = "fov", assay = ASSAY, molecule.coordinates = FALSE)
fov <- obj[["fov"]]
centroids <- fov[["centroids"]]
coords <- GetTissueCoordinates(fov)

container <- list(
  n_cells = ncol(obj),
  n_features = nrow(obj),
  cells = name_anchor(colnames(obj)),
  features = name_anchor(rownames(obj)),
  n_images = length(Images(obj)),
  # I() keeps a length-1 vector an array: auto_unbox would collapse it to a bare
  # string, and the Python side always emits a list here. That is a difference in
  # the two JSON writers, not in the two libraries, and it does not belong in the
  # parity table.
  boundaries = I(sort(Boundaries(fov))),
  default_boundary = DefaultBoundary(fov),
  fov_assay = DefaultAssay(fov),
  fov_radius_is_none = is.null(Radius(fov)),
  centroids_n = length(Cells(centroids)),
  centroids_nsides = as.integer(slot(centroids, "nsides")),
  centroids_radius = as.numeric(Radius(centroids)),
  coords_shape = as.integer(dim(coords)),
  coords_x_head = as.numeric(utils::head(coords$x, 3)),
  coords_y_head = as.numeric(utils::head(coords$y, 3)),
  coords_cells_head = as.character(utils::head(coords$cell, 3))
)

# ---- the constructors, on input small enough to read ------------------------
cat("Constructors on toy input...\n")
toy_coords <- data.frame(x = c(1, 2, 3, 4), y = c(10, 20, 30, 40),
                         cell = c("a", "b", "c", "d"))
toy_centroids <- CreateCentroids(toy_coords)
toy_fov <- CreateFOV(toy_coords, type = "centroids", assay = "RNA")
square <- data.frame(x = c(0, 1, 1, 0, 5, 6, 6, 5),
                     y = c(0, 0, 1, 1, 5, 5, 6, 6),
                     cell = rep(c("a", "b"), each = 4))
seg <- CreateSegmentation(square)
seg_coords <- GetTissueCoordinates(seg)
ring_a <- seg_coords[seg_coords$cell == "a", c("x", "y")]

toy <- list(
  auto_radius = as.numeric(Radius(toy_centroids)),
  nsides = as.integer(slot(toy_centroids, "nsides")),
  centroid_cells = as.character(Cells(toy_centroids)),
  fov_cells = as.character(Cells(toy_fov)),
  fov_boundaries = I(sort(Boundaries(toy_fov))),
  subset_cells = as.character(Cells(subset(toy_fov, cells = c("b", "d")))),
  segmentation_cells = as.character(Cells(seg)),
  segmentation_rows = nrow(seg_coords),
  segmentation_rows_per_cell = as.integer(
    sapply(Cells(seg), function(c) sum(seg_coords$cell == c))
  ),
  ring_closed = isTRUE(all.equal(unlist(ring_a[1, ]), unlist(ring_a[nrow(ring_a), ]),
                                 check.attributes = FALSE))
)

# ---- Moran's I on the shared subset -----------------------------------------
cells_file <- file.path(FIG, "cells.txt")
if (!file.exists(cells_file)) {
  stop("figures_svf/cells.txt is missing; run the Python tutorial first.",
       call. = FALSE)
}
cells <- as.character(read.csv(cells_file, header = FALSE)[[1]])
cat(sprintf("Moran's I on the shared %d-cell subset...\n", length(cells)))

obj <- NormalizeData(obj, verbose = FALSE)
data <- as.matrix(GetAssayData(obj, assay = ASSAY, layer = "data")[, cells])
pos <- coords[match(cells, coords$cell), c("x", "y")]

svf <- FindSpatiallyVariableFeatures(
  data, spatial.location = pos, selection.method = "moransi", verbose = FALSE
)
svf <- svf[order(-svf$observed), , drop = FALSE]

moransi <- list(
  n_cells = length(cells),
  n_genes = nrow(svf),
  top_genes = rownames(svf)[seq_len(TOP_N)],
  i_head = as.numeric(svf$observed[seq_len(TOP_N)]),
  i_max = max(svf$observed),
  i_min = min(svf$observed),
  ranking = name_digest(rownames(svf))
)

# The full per-gene series, for the figures. The anchors carry only the top 10;
# a scatter drawn from those would be a scatter of the agreement, which is the
# one thing a parity figure must not do.
write.csv(svf, file.path(FIG, "r_moransi.csv"))

anchors <- list(container = container, toy = toy, moransi = moransi)
out <- file.path(FIG, "r_anchors.json")
write(toJSON(anchors, auto_unbox = TRUE, digits = 15, pretty = TRUE), out)
cat("Wrote", out, "\n")
cat(sprintf("  top %d: %s\n", TOP_N, paste(moransi$top_genes, collapse = ", ")))
cat(sprintf("  I range: %.4f .. %.4f\n", moransi$i_min, moransi$i_max))
