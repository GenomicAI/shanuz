#!/usr/bin/env Rscript
# R Seurat reference for the integration tutorial (integration_vignette.md).
#
# Mirrors ifnb_integration_tutorial.py: builds the RNA object from the same
# exported ifnb counts, runs the three integration methods (Harmony, CCA, RPCA)
# via Seurat v5's IntegrateLayers, clusters each corrected reduction, and writes
# the per-cell cluster calls (r_calls.csv) that the Python tutorial's
# report_concordance() compares against — the references for shanuz.run_harmony
# / integrate_layers. Also writes the R-side UMAP figures.
#
# Integration embeddings are not coordinate-comparable across tools, so the
# concordance is partition-based (adjusted Rand index of the clusterings, plus
# each tool's cell-type recovery and batch-mixing) — all computable from the
# cluster labels this script writes.
#
# To keep the reduction on ONE shared gene basis, this reads the variable
# features the Python run selected (figures_integration/hvg_features.txt) instead
# of running FindVariableFeatures here, so the only divergences left are the
# genuinely method-level ones (PCA numerics, the integration algorithms, Louvain).
#
# Data: the same counts the Python tutorial reads. Run
#   Rscript tutorials/export_seuratdata.R ifnb          # one-time counts export
#   python  tutorials/ifnb_integration_tutorial.py      # writes the shared HVGs
# then
#   Rscript tutorials/ifnb_integration_verify.R
# Override the data folder with the IFNB_DATA environment variable.
#
# Needs: Seurat, ggplot2, Matrix, harmony.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(Matrix)
})
set.seed(42)
# Seurat v5's IntegrateLayers ships each split layer to a future worker; the
# ifnb layers exceed future's default 500 MiB per-export cap (RPCAIntegration
# alone needs ~1.6 GiB), so raise the ceiling. Not a correctness knob.
options(future.globals.maxSize = 3 * 1024^3)

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
HERE <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG  <- file.path(HERE, "figures_integration")
DATA <- Sys.getenv("IFNB_DATA", path.expand("~/.shanuz_data/ifnb"))
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

HVG_TXT <- file.path(FIG, "hvg_features.txt")
if (!file.exists(file.path(DATA, "matrix.mtx.gz")))
  stop("ifnb counts not found at ", DATA,
       "\nRun `Rscript tutorials/export_seuratdata.R ifnb` first.")
if (!file.exists(HVG_TXT))
  stop("hvg_features.txt not found in ", FIG,
       "\nRun `python tutorials/ifnb_integration_tutorial.py` first (it writes it).")

RESOLUTION <- 0.5
N_PCS <- 30
DIMS <- 1:N_PCS

# ---- 1. Load the same exported counts + metadata ----------------------------
cat("Reading ifnb counts ...\n")
counts <- Read10X(DATA)                                   # features x cells
meta <- read.csv(file.path(DATA, "metadata.csv"), row.names = 1, check.names = FALSE)
meta <- meta[colnames(counts), , drop = FALSE]
obj <- CreateSeuratObject(counts = counts, min.cells = 3, meta.data = meta)
cat(sprintf("RNA %d genes x %d cells | %s\n", nrow(obj), ncol(obj),
            paste(names(table(obj$stim)), table(obj$stim), collapse = " ")))

# ---- 2. Prep to PCA on the SHARED variable-feature basis --------------------
hvg <- readLines(HVG_TXT)
hvg <- hvg[hvg %in% rownames(obj)]
cat(sprintf("Using %d shared variable features from Python\n", length(hvg)))
obj <- NormalizeData(obj, verbose = FALSE)
VariableFeatures(obj) <- hvg
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = N_PCS, verbose = FALSE)

# Split the RNA assay into per-condition layers — the v5 IntegrateLayers API
# integrates across layers.
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$stim)
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = N_PCS, verbose = FALSE)

# ---- 3. The three integrations, each to its own reduction -------------------
cluster_on <- function(object, reduction, key) {
  object <- FindNeighbors(object, reduction = reduction, dims = DIMS, verbose = FALSE)
  object <- FindClusters(object, resolution = RESOLUTION, verbose = FALSE)
  object[[key]] <- as.integer(as.character(object$seurat_clusters))
  object
}

# baseline: cluster the uncorrected PCA
obj <- cluster_on(obj, "pca", "R_pca")

cat("Harmony ...\n")
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca",
                       new.reduction = "harmony", verbose = FALSE)
obj <- cluster_on(obj, "harmony", "R_harmony")

cat("CCA ...\n")
obj <- IntegrateLayers(obj, method = CCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.cca", verbose = FALSE)
obj <- cluster_on(obj, "integrated.cca", "R_cca")

cat("RPCA ...\n")
obj <- IntegrateLayers(obj, method = RPCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.rpca", verbose = FALSE)
obj <- cluster_on(obj, "integrated.rpca", "R_rpca")

obj <- JoinLayers(obj)

# ---- 4. Per-cell calls for the Python concordance report (write first) ------
calls <- data.frame(
  cell               = colnames(obj),
  stim               = as.character(obj$stim),
  seurat_annotations = as.character(obj$seurat_annotations),
  R_pca              = obj$R_pca,
  R_harmony          = obj$R_harmony,
  R_cca              = obj$R_cca,
  R_rpca             = obj$R_rpca,
  stringsAsFactors   = FALSE
)
write.csv(calls, file.path(FIG, "r_calls.csv"), row.names = FALSE)
cat("\nWrote r_calls.csv\n")

# quick R-side batch-mixing sanity print (partition-only, matches Python's metric)
batch_entropy <- function(clusters, batch) {
  levels <- unique(batch); norm <- log(length(levels)); tot <- 0
  for (c in unique(clusters)) {
    m <- clusters == c; p <- sapply(levels, function(l) mean(batch[m] == l))
    p <- p[p > 0]; tot <- tot + sum(m) * (-sum(p * log(p)) / norm)
  }
  tot / length(clusters)
}
for (k in c("R_pca", "R_harmony", "R_cca", "R_rpca"))
  cat(sprintf("  %-10s n_clusters=%2d  batch_mix=%.3f\n", k,
              length(unique(obj[[k]][, 1])), batch_entropy(obj[[k]][, 1], obj$stim)))

# ---- 5. Figures (r_ prefix, side by side with the Python ones) --------------
# Guarded so a plotting hiccup never costs the already-written calls.
try({
  obj <- RunUMAP(obj, reduction = "pca", dims = DIMS,
                 reduction.name = "umap.pca", verbose = FALSE)
  ggsave(file.path(FIG, "r_01_uncorrected_stim.png"),
         DimPlot(obj, reduction = "umap.pca", group.by = "stim") +
           ggtitle("Uncorrected (PCA) — by condition"),
         width = 6, height = 5, dpi = 150)
}, silent = FALSE)

try({
  obj <- RunUMAP(obj, reduction = "harmony", dims = DIMS,
                 reduction.name = "umap.harmony", verbose = FALSE)
  ggsave(file.path(FIG, "r_02_harmony_stim.png"),
         DimPlot(obj, reduction = "umap.harmony", group.by = "stim") +
           ggtitle("Harmony — by condition"),
         width = 6, height = 5, dpi = 150)
  ggsave(file.path(FIG, "r_03_harmony_celltype.png"),
         DimPlot(obj, reduction = "umap.harmony", group.by = "seurat_annotations") +
           ggtitle("Harmony — by cell type"),
         width = 7, height = 5, dpi = 150)
}, silent = FALSE)

cat("\nDONE — wrote r_calls.csv + figures to", FIG, "\n")
