#!/usr/bin/env Rscript
# R Seurat reference for the two SCTransform-tutorial figures the published
# Seurat vignette does not include (sctransform_vignette.md): the annotated
# cell-type UMAP and the SCTransform-vs-LogNormalize comparison. Every other R
# figure there links the canonical satijalab.org vignette image; these two had
# no R counterpart, so this script generates them.
#
# Mirrors pbmc3k_sctransform_tutorial.py: runs BOTH the SCTransform workflow
# (regress percent.mt, dims 1:30) and the standard LogNormalize workflow
# (dims 1:10) on the same cells, annotates by relative marker enrichment, and
# writes into tutorials/figures_sctransform/:
#   * r_02_sct_umap_celltypes.png
#   * r_06_sct_vs_std_umap.png
#
# Data: the same cache the Python tutorial downloads to. Run
#   python tutorials/pbmc3k_sctransform_tutorial.py   # downloads ~24 MB first
# then
#   Rscript tutorials/pbmc3k_sctransform_verify.R
# Needs: Seurat, sctransform, ggplot2, patchwork.
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })
set.seed(0)

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
HERE <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG  <- file.path(HERE, "figures_sctransform")
DATA <- Sys.getenv("PBMC3K_DATA",
                   path.expand("~/.shanuz_data/pbmc3k/filtered_gene_bc_matrices/hg19"))
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(DATA, "matrix.mtx")))
  stop("PBMC 3k data not found at ", DATA,
       "\nRun `python tutorials/pbmc3k_sctransform_tutorial.py` first.")

# Fine-grained lineage panel (identical to pbmc3k_sctransform_tutorial.py).
FINE_MARKERS <- list(
  "Naive CD4 T"   = c("IL7R","CCR7","LEF1","SELL"),
  "Memory CD4 T"  = c("IL7R","S100A4","IL32","ANXA1"),
  "CD8 Naive/Mem" = c("CD8A","CD8B","CCR7"),
  "CD8 Effector"  = c("CD8A","GZMK","CCL5","NKG7"),
  "B"             = c("MS4A1","CD79A","TCL1A","FCER2"),
  "CD14+ Mono"    = c("CD14","LYZ","S100A8","S100A9"),
  "FCGR3A+ Mono"  = c("FCGR3A","MS4A7"),
  "NK"            = c("GNLY","NKG7","KLRD1","XCL1"),
  "DC"            = c("FCER1A","CST3"),
  "pDC"           = c("SERPINF1","ITM2C"),
  "Platelet"      = c("PPBP","PF4")
)

# z-score relative-enrichment annotator (same logic as the other verify scripts).
annotate_clusters <- function(obj, marker_sets, assay) {
  data <- GetAssayData(obj, assay = assay, layer = "data")
  idents <- as.character(Idents(obj))
  clusters <- as.character(sort(unique(as.integer(idents))))
  needed <- intersect(unique(unlist(marker_sets)), rownames(data))
  zmean <- list()
  for (g in needed) {
    expr <- data[g, ]
    per_cluster <- vapply(clusters, function(c) mean(expr[idents == c]), numeric(1))
    sdv <- sd(per_cluster)
    zmean[[g]] <- if (sdv > 1e-9) (per_cluster - mean(per_cluster)) / sdv else rep(0, length(clusters))
  }
  assignment <- character()
  for (ci in seq_along(clusters)) {
    best <- "Unknown"; best_score <- -Inf
    for (lineage in names(marker_sets)) {
      present <- intersect(marker_sets[[lineage]], names(zmean))
      if (length(present) == 0) next
      score <- mean(vapply(present, function(g) zmean[[g]][ci], numeric(1)))
      if (score > best_score) { best_score <- score; best <- lineage }
    }
    assignment[clusters[ci]] <- best
  }
  assignment
}

load_pbmc <- function() {
  o <- CreateSeuratObject(Read10X(DATA), project = "pbmc3k", min.cells = 3, min.features = 200)
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
  o
}

# ---- SCTransform workflow (regress percent.mt, dims 1:30) -------------------
cat("Running SCTransform workflow...\n")
sct <- load_pbmc()
sct <- SCTransform(sct, vars.to.regress = "percent.mt", variable.features.n = 3000,
                   verbose = FALSE)
sct <- RunPCA(sct, npcs = 50, verbose = FALSE)
sct <- FindNeighbors(sct, dims = 1:30, k.param = 20, verbose = FALSE)
sct <- FindClusters(sct, resolution = 0.8, algorithm = 1, verbose = FALSE)
sct <- RunUMAP(sct, dims = 1:30, verbose = FALSE)
n_sct <- length(levels(sct))
sct$sct_clusters <- Idents(sct)
sct_anno <- annotate_clusters(sct, FINE_MARKERS, assay = "SCT")
cat(sprintf("SCTransform: %d clusters\n", n_sct)); print(sct_anno)
Idents(sct) <- sct$sct_clusters
sct <- RenameIdents(sct, sct_anno)
sct$sct_celltype <- Idents(sct)

# ---- Standard LogNormalize workflow (dims 1:10) ----------------------------
cat("Running standard LogNormalize workflow...\n")
std <- load_pbmc()
std <- NormalizeData(std, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
std <- FindVariableFeatures(std, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
std <- ScaleData(std, features = VariableFeatures(std), verbose = FALSE)
std <- RunPCA(std, npcs = 50, verbose = FALSE)
std <- FindNeighbors(std, dims = 1:10, k.param = 20, verbose = FALSE)
std <- FindClusters(std, resolution = 0.8, algorithm = 1, verbose = FALSE)
std <- RunUMAP(std, dims = 1:10, verbose = FALSE)
n_std <- length(levels(std))
std$std_clusters <- Idents(std)
cat(sprintf("LogNormalize: %d clusters\n", n_std))

# ============================ FIGURES =======================================
# 02 SCT UMAP - annotated cell types
p2 <- DimPlot(sct, reduction = "umap", group.by = "sct_celltype", label = TRUE) +
  ggtitle("R Seurat - SCTransform cell types") +
  theme(plot.title = element_text(size = 12, face = "bold"))
ggsave(file.path(FIG, "r_02_sct_umap_celltypes.png"), p2, width = 8.5, height = 6.5, dpi = 150, bg = "white")

# 06 SCTransform vs LogNormalize, side by side
p_sct <- DimPlot(sct, reduction = "umap", group.by = "sct_clusters", label = TRUE) +
  ggtitle(sprintf("SCTransform - %d clusters (dims 1:30)", n_sct)) + NoLegend()
p_std <- DimPlot(std, reduction = "umap", group.by = "std_clusters", label = TRUE) +
  ggtitle(sprintf("LogNormalize - %d clusters (dims 1:10)", n_std)) + NoLegend()
p6 <- (p_sct | p_std) +
  plot_annotation(title = "R Seurat - SCTransform vs standard log-normalization",
                  theme = theme(plot.title = element_text(size = 13, face = "bold")))
ggsave(file.path(FIG, "r_06_sct_vs_std_umap.png"), p6, width = 15, height = 6.5, dpi = 150, bg = "white")

cat("\nWrote r_02_sct_umap_celltypes.png and r_06_sct_vs_std_umap.png to", FIG, "\n")
