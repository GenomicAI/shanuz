#!/usr/bin/env Rscript
# R Seurat reference for the object-internals tutorial (objects_vignette.md).
#
# Builds the same PBMC 3k object as tutorials/pbmc3k_objects_tutorial.py and
# writes the same anchor tree to tutorials/figures_objects/r_anchors.json, so
# that compare_anchors() can put the two side by side field for field.
#
# Run the Python side first — it downloads the data and publishes the shared
# variable-feature and cell lists this script reads back:
#   python tutorials/pbmc3k_objects_tutorial.py
#   Rscript tutorials/pbmc3k_objects_verify.R
#
# Needs: Seurat, SeuratObject, Matrix, jsonlite, digest.
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(jsonlite); library(digest)
})
set.seed(42)

# Run from the repo root or from tutorials/ — both are common.
FIGURES <- if (dir.exists("tutorials")) {
  file.path("tutorials", "figures_objects")
} else {
  "figures_objects"
}
DATA <- Sys.getenv("PBMC3K_DATA",
                   file.path(path.expand("~"), ".shanuz_data", "pbmc3k",
                             "filtered_gene_bc_matrices", "hg19"))

MIN_CELLS <- 3
MIN_FEATURES <- 200
N_PCS <- 20
PROJECT <- "pbmc3k_objects"
IDENT_MARKERS <- list(c("T", "CD3E"), c("B", "MS4A1"), c("Mono", "LYZ"))
IDENT_OTHER <- "Other"

# ---------------------------------------------------------------------------
# Anchor vocabulary — mirrors the Python helpers exactly
# ---------------------------------------------------------------------------

# Order-sensitive fingerprint. `serialize = FALSE` hashes the string's own bytes
# rather than R's serialization envelope, which is what Python's md5 sees.
digest12 <- function(x) {
  substr(digest(paste(as.character(x), collapse = "\n"),
                algo = "md5", serialize = FALSE), 1, 12)
}

# `toJSON(auto_unbox = TRUE)` writes a length-1 vector as a bare scalar, which
# would make a one-assay object's `["RNA"]` compare unequal to R's `"RNA"` for
# no reason. Anything that is semantically a *vector* gets wrapped so it stays
# an array whatever its length.
arr <- function(x) I(as.character(x))

name_anchor <- function(x) {
  x <- as.character(x)
  n <- length(x)
  list(n = n, digest = digest12(x),
       head = arr(if (n >= 3) x[1:3] else x),
       tail = arr(if (n >= 3) x[(n - 2):n] else x))
}

matrix_anchor <- function(m) {
  list(shape = c(nrow(m), ncol(m)),
       nnz = as.integer(sum(m != 0)),
       sum = round(sum(m), 6))
}

layer_anchor <- function(assay) {
  out <- list()
  for (nm in Layers(assay)) {
    m <- LayerData(assay, nm)
    a <- matrix_anchor(m)
    a$cells <- name_anchor(colnames(m))
    out[[nm]] <- a
  }
  out
}

ident_anchor <- function(idents) {
  s <- as.character(idents)
  tab <- table(s)
  lv <- sort(names(tab))
  list(levels = arr(lv),
       counts = setNames(as.list(as.integer(tab[lv])), lv),
       digest = digest12(s))
}

# ---------------------------------------------------------------------------
# Build the object
# ---------------------------------------------------------------------------

cat("Reading", DATA, "\n")
counts <- Read10X(data.dir = DATA)
obj <- CreateSeuratObject(counts = counts, project = PROJECT,
                          min.cells = MIN_CELLS, min.features = MIN_FEATURES)
cat("object:", ncol(obj), "cells x", nrow(obj), "features\n")

# Share the Python side's variable features so both objects rest on one basis.
hvg_path <- file.path(FIGURES, "hvg_features.txt")
stopifnot("run the Python tutorial first" = file.exists(hvg_path))
# Read10X rewrites '_' to '-' in gene symbols (pbmc3k has Y_RNA, RP11-*_* and
# friends); shanuz's loader leaves them alone. Normalise to R's spelling before
# matching — the difference belongs to the two loaders, not to the object model
# this tutorial is auditing. The Python side hashes the same normalised names.
hvg <- gsub("_", "-", readLines(hvg_path))

obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000,
                            verbose = FALSE)
VariableFeatures(obj) <- hvg
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = N_PCS, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)

# Identities from the same deterministic marker gates, applied in the same order.
assign_idents <- function(obj) {
  cells <- colnames(obj)
  labels <- rep(IDENT_OTHER, length(cells))
  for (m in rev(IDENT_MARKERS)) {
    gene <- m[2]
    if (!gene %in% rownames(obj)) next
    x <- as.numeric(LayerData(obj[["RNA"]], "data")[gene, cells])
    labels[x > 0] <- m[1]
  }
  labels
}
Idents(obj) <- assign_idents(obj)

# ---------------------------------------------------------------------------
# Collect the anchors
# ---------------------------------------------------------------------------

assay <- obj[["RNA"]]
pca <- obj[["pca"]]

anchors <- list(
  shape = list(n_cells = ncol(obj), n_features = nrow(obj)),
  cells = name_anchor(colnames(obj)),
  features = name_anchor(rownames(obj)),

  assay = list(
    names = arr(Assays(obj)),
    default = DefaultAssay(obj),
    key = unname(Key(assay)),
    class = class(assay)[1]
  ),
  layers = layer_anchor(assay),

  meta = list(
    columns = arr(sort(colnames(obj[[]]))),
    nCount_RNA_sum = round(sum(obj$nCount_RNA), 6),
    nFeature_RNA_sum = round(sum(obj$nFeature_RNA), 6)
  ),
  idents = ident_anchor(Idents(obj)),

  variable_features = name_anchor(hvg),

  reductions = list(
    names = arr(sort(Reductions(obj))),
    pca_key = unname(Key(pca)),
    pca_dims = dim(Embeddings(pca)),
    pca_n_loadings = dim(Loadings(pca)),
    pca_stdev_head = round(Stdev(pca)[1:5], 6)
  ),

  graphs = {
    g <- list()
    for (nm in sort(Graphs(obj))) {
      m <- obj[[nm]]
      g[[nm]] <- list(shape = c(nrow(m), ncol(m)),
                      nnz = as.integer(sum(m != 0)))
    }
    g
  },

  # Seurat appends one entry per pipeline call. This is the list of what the
  # object believes has been run on it.
  commands = arr(names(obj@commands))
)

# --- WhichCells / RenameIdents / subset ------------------------------------
t_cells <- WhichCells(obj, idents = "T")
renamed <- RenameIdents(obj, "T" = "T_cell")
sub <- subset(obj, idents = "Mono")
anchors$ident_ops <- list(
  which_cells_T = name_anchor(t_cells),
  renamed_levels = arr(sort(unique(as.character(Idents(renamed))))),
  subset_mono = list(
    n_cells = ncol(sub),
    cells = name_anchor(colnames(sub)),
    n_features = nrow(sub)
  )
)

# --- FetchData --------------------------------------------------------------
# All three kinds of variable in one call, each addressed by the name Seurat
# knows it by — `PC_1` comes from the reduction's Key().
fetched <- FetchData(obj, vars = c("nCount_RNA", "nFeature_RNA", "CD3E", "PC_1"))
anchors$fetch <- list(
  columns = arr(colnames(fetched)),
  n_rows = nrow(fetched),
  CD3E_sum = round(sum(fetched$CD3E), 6),
  nCount_RNA_sum = round(sum(fetched$nCount_RNA), 6),
  # Exact, and within one tool: a fetched embedding column must be the object's
  # embedding column, whatever the two PCAs do differently.
  pc1_matches_embeddings = identical(as.numeric(fetched$PC_1),
                                     as.numeric(Embeddings(pca)[, 1])),
  abs_PC_1_sum = round(sum(abs(fetched$PC_1)), 6)
)

# --- split / JoinLayers round trip -----------------------------------------
# On a counts-only assay, matching the Python side: the question is whether the
# round trip is the identity, not how each API spells "split every layer".
batch <- ifelse(seq_len(ncol(obj)) %% 2 == 1, "batch1", "batch2")
fresh <- CreateSeuratObject(counts = counts, project = PROJECT,
                            min.cells = MIN_CELLS, min.features = MIN_FEATURES)
before <- LayerData(fresh[["RNA"]], "counts")
original_cells <- colnames(before)

split_assay <- split(fresh[["RNA"]], f = batch)
joined_assay <- JoinLayers(split_assay)
joined_names <- Layers(joined_assay)
target <- joined_names[1]
after <- LayerData(joined_assay, target)

anchors$split_join <- list(
  layers_before = arr(Layers(fresh[["RNA"]])),
  layers_after_split = arr(Layers(split_assay)),
  layers_after_join = arr(joined_names),
  layer_name_restored = identical(target, "counts"),
  cell_order_restored = identical(colnames(after), original_cells),
  matrix_restored = identical(dim(after), dim(before)) &&
    all(after == before[, colnames(after)]),
  cells_after_join = name_anchor(colnames(after)),
  assay_cells = name_anchor(colnames(fresh[["RNA"]]))
)

# --- JoinLayers on the prepared object -------------------------------------
# The call every real script makes, on an assay that also carries `data` and a
# variable-features-only `scale.data`.
prepared <- obj
prepared[["RNA"]] <- split(prepared[["RNA"]], f = batch)
anchors$join_all_layers <- tryCatch({
  rejoined <- JoinLayers(prepared[["RNA"]])
  list(error = NULL, layers_after_join = arr(Layers(rejoined)))
}, error = function(e) list(error = paste0(class(e)[1], ": ", conditionMessage(e)),
                            layers_after_join = NULL))

# ---------------------------------------------------------------------------
dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)
out <- file.path(FIGURES, "r_anchors.json")
write(toJSON(anchors, auto_unbox = TRUE, digits = NA, null = "null"), out)
cat("wrote", out, "\n")
cat("commands logged by Seurat:", paste(names(obj@commands), collapse = ", "), "\n")
