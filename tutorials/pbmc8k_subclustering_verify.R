#!/usr/bin/env Rscript
# R Seurat reference for the advanced PBMC 8k subclustering tutorial
# (advanced_pbmc8k_subclustering.md).
#
# Mirrors pbmc8k_subclustering_tutorial.py / generate_advanced_plots.py as
# faithfully as the two implementations allow (same QC, same pipeline
# parameters, and the same two annotation heuristics ported to R), then writes
# the R-side figures for the side-by-side tables into tutorials/figures_advanced/:
#   * r_01_qc_violin.png ... r_11_tnk_markers_heatmap.png
#
# Data: reads the same cache the Python tutorial downloads to. Run
#   python tutorials/pbmc8k_subclustering_tutorial.py   # downloads ~38 MB first
# then
#   Rscript tutorials/pbmc8k_subclustering_verify.R
# Override the data folder with the PBMC8K_DATA environment variable.
#
# Needs: Seurat, ggplot2.
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
set.seed(0)

# --- resolve paths relative to this script ---------------------------------
.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
HERE <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG  <- file.path(HERE, "figures_advanced")
DATA <- Sys.getenv("PBMC8K_DATA",
                   path.expand("~/.shanuz_data/pbmc8k/filtered_gene_bc_matrices/GRCh38"))
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(DATA, "matrix.mtx")) &&
    !file.exists(file.path(DATA, "matrix.mtx.gz")))
  stop("PBMC 8k data not found at ", DATA,
       "\nRun `python tutorials/pbmc8k_subclustering_tutorial.py` first.")

# ---- marker panels (identical to pbmc8k_subclustering_tutorial.py) ----------
BROAD_MARKERS <- list(
  "CD4 T"        = c("IL7R","CD3D","CCR7"),
  "CD8 T"        = c("CD8A","CD8B","GZMK"),
  "NK"           = c("GNLY","NKG7","KLRD1"),
  "B"            = c("MS4A1","CD79A","CD79B"),
  "CD14+ Mono"   = c("CD14","LYZ","S100A8"),
  "FCGR3A+ Mono" = c("FCGR3A","MS4A7"),
  "DC"           = c("FCER1A","CST3"),
  "Platelet"     = c("PPBP","PF4")
)
LYMPHOID_LINEAGES <- c("CD4 T","CD8 T","NK")
TNK_PANEL <- c("CD3D","CD3E","CD8A","CD8B","GNLY","NKG7","KLRD1",
               "CCR7","SELL","LEF1","IL7R","S100A4","GZMK")

# ---- annotation heuristics ported from the Python tutorial ------------------
# z-score each marker's per-cluster mean across clusters; each lineage scores as
# the mean z of its present markers; argmax wins. (np.std ddof=0 vs R sd ddof=1
# is a per-gene constant factor and does not change the argmax.)
annotate_clusters <- function(obj, marker_sets) {
  data <- GetAssayData(obj, assay = "RNA", layer = "data")
  idents <- as.character(Idents(obj))
  clusters <- as.character(sort(unique(as.integer(idents))))
  feats <- rownames(data)
  needed <- intersect(unique(unlist(marker_sets)), feats)
  zmean <- list()
  for (g in needed) {
    expr <- data[g, ]
    per_cluster <- vapply(clusters, function(c) mean(expr[idents == c]), numeric(1))
    sdv <- sd(per_cluster)
    zmean[[g]] <- if (sdv > 1e-9) (per_cluster - mean(per_cluster)) / sdv
                  else rep(0, length(clusters))
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

# Hierarchical lineage-priority gating: NK -> CD8/cytotoxic -> CD4 naive -> CD4 memory.
annotate_tnk_subsets <- function(sub) {
  data <- GetAssayData(sub, assay = "RNA", layer = "data")
  idents <- as.character(Idents(sub))
  clusters <- as.character(sort(unique(as.integer(idents))))
  panel <- intersect(TNK_PANEL, rownames(data))
  expr <- data[panel, , drop = FALSE]
  m <- function(genes, mask) {
    gg <- intersect(genes, panel)
    if (length(gg) == 0) return(0.0)
    mean(vapply(gg, function(g) mean(expr[g, mask]), numeric(1)))
  }
  assignment <- character()
  for (c in clusters) {
    mask <- idents == c
    cd3   <- m(c("CD3D","CD3E"), mask)
    cd8   <- m(c("CD8A","CD8B"), mask)
    nk    <- m(c("GNLY","NKG7"), mask)
    cyto  <- m(c("NKG7","GZMK"), mask)
    naive <- m(c("CCR7","SELL","LEF1"), mask)
    if (cd3 < 0.75 && nk > 1.5)                         assignment[c] <- "NK"
    else if (cd8 > 0.6 || (cd3 >= 0.75 && cyto > 1.5))  assignment[c] <- "CD8 T"
    else if (naive > 0.9)                               assignment[c] <- "CD4 Naive"
    else                                                assignment[c] <- "CD4 Memory"
  }
  assignment
}

top_genes <- function(markers, n) {
  genes <- unlist(lapply(split(markers, markers$cluster), function(df)
    df$gene[order(-df$avg_log2FC)][seq_len(min(n, nrow(df)))]))
  unique(genes)
}

# ---- 1. Load + QC -----------------------------------------------------------
counts <- Read10X(DATA)
pbmc <- CreateSeuratObject(counts = counts, project = "pbmc8k",
                           min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
n_raw <- ncol(pbmc)
raw <- pbmc  # pre-filter snapshot for the QC violin
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cat(sprintf("Loaded %d cells -> %d after QC\n", n_raw, ncol(pbmc)))

# ---- 2. Standard workflow ---------------------------------------------------
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize",
                      scale.factor = 10000, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pbmc <- ScaleData(pbmc, features = rownames(pbmc), verbose = FALSE)
pbmc <- RunPCA(pbmc, npcs = 50, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:10, k.param = 20, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, algorithm = 1, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
cat(sprintf("Global clusters (res=0.5): %d\n", length(levels(pbmc))))

# ---- 3. Markers + broad annotation -----------------------------------------
all_markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25, verbose = FALSE)
pbmc$global_clusters <- Idents(pbmc)
broad <- annotate_clusters(pbmc, BROAD_MARKERS)
cat("Broad lineage assignment:\n"); print(broad)
Idents(pbmc) <- pbmc$global_clusters
pbmc <- RenameIdents(pbmc, broad)
pbmc$broad_celltype <- Idents(pbmc)

# ---- 4. Subcluster the T/NK compartment ------------------------------------
lymphoid_clusters <- names(broad)[broad %in% LYMPHOID_LINEAGES]
sub_cells <- colnames(pbmc)[as.character(pbmc$global_clusters) %in% lymphoid_clusters]
cat(sprintf("Lymphoid clusters %s -> %d T/NK cells\n",
            paste(sort(as.integer(lymphoid_clusters)), collapse = ","), length(sub_cells)))
sub <- subset(pbmc, cells = sub_cells)
# data layer already normalised -> re-run from HVG onward (no NormalizeData)
sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sub <- ScaleData(sub, features = rownames(sub), verbose = FALSE)
sub <- RunPCA(sub, npcs = 30, verbose = FALSE)
sub <- FindNeighbors(sub, dims = 1:10, k.param = 20, verbose = FALSE)
sub <- FindClusters(sub, resolution = 0.6, algorithm = 1, verbose = FALSE)
sub <- RunUMAP(sub, dims = 1:10, verbose = FALSE)
cat(sprintf("T/NK subclusters (res=0.6): %d\n", length(levels(sub))))

sub_markers <- FindAllMarkers(sub, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25, verbose = FALSE)
sub$sub_clusters <- Idents(sub)
sub_anno <- annotate_tnk_subsets(sub)
cat("T/NK subset assignment:\n"); print(sub_anno)
Idents(sub) <- sub$sub_clusters
sub <- RenameIdents(sub, sub_anno)
sub$tnk_subset <- Idents(sub)

# ============================ FIGURES =======================================
sv <- function(p, name, w = 8, h = 6.5)
  ggsave(file.path(FIG, name), p, width = w, height = h, dpi = 150, bg = "white")
titled <- function(p, t) p + ggtitle(t) +
  theme(plot.title = element_text(size = 12, face = "bold"))

# 01 QC violin (pre-filter, like generate_advanced_plots.py)
sv(VlnPlot(raw, c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3, pt.size = 0) &
     theme(plot.title = element_text(size = 11)),
   "r_01_qc_violin.png", 12, 4)

# 02 Elbow
sv(titled(ElbowPlot(pbmc, ndims = 30), "R Seurat - elbow"), "r_02_elbow_plot.png", 7, 4)

# 03 UMAP global clusters
sv(titled(DimPlot(pbmc, reduction = "umap", group.by = "global_clusters", label = TRUE),
          "R Seurat - global clusters"), "r_03_umap_global_clusters.png")

# 04 UMAP broad cell types
sv(titled(DimPlot(pbmc, reduction = "umap", group.by = "broad_celltype", label = TRUE),
          "R Seurat - broad cell types"), "r_04_umap_global_celltypes.png", 8.5, 6.5)

# 05 lineage marker feature plots
lineage_markers <- c("CD3D","CD8A","IL7R","MS4A1","LYZ","FCGR3A","GNLY","FCER1A","PPBP")
sv(FeaturePlot(pbmc, lineage_markers, reduction = "umap", ncol = 3),
   "r_05_lineage_featureplots.png", 12, 10)

# 06 global markers heatmap (top 5 per cluster)
top_global <- top_genes(all_markers, 5)
sv(DoHeatmap(pbmc, features = top_global, group.by = "global_clusters") +
     theme(axis.text.y = element_text(size = 6)),
   "r_06_global_markers_heatmap.png", 14, max(6, length(top_global) * 0.22))

# 07 UMAP T/NK subclusters
sv(titled(DimPlot(sub, reduction = "umap", group.by = "sub_clusters", label = TRUE),
          "R Seurat - T/NK subclusters"), "r_07_umap_tnk_subclusters.png")

# 08 UMAP T/NK annotated subsets
sv(titled(DimPlot(sub, reduction = "umap", group.by = "tnk_subset", label = TRUE),
          "R Seurat - T/NK subsets"), "r_08_umap_tnk_subsets.png", 8.5, 6.5)

# 09 T/NK subset feature plots
tnk_markers <- c("CCR7","SELL","IL7R","S100A4","CD8A","GZMK","GNLY","NKG7")
sv(FeaturePlot(sub, tnk_markers, reduction = "umap", ncol = 4),
   "r_09_tnk_subset_featureplots.png", 15, 7)

# 10 T/NK subset violins
sv(VlnPlot(sub, c("CCR7","S100A4","CD8A","GNLY"), group.by = "tnk_subset", ncol = 2),
   "r_10_tnk_subset_violins.png", 11, 7)

# 11 T/NK markers heatmap (top 6 per subcluster)
top_sub <- top_genes(sub_markers, 6)
sv(DoHeatmap(sub, features = top_sub, group.by = "sub_clusters") +
     theme(axis.text.y = element_text(size = 6)),
   "r_11_tnk_markers_heatmap.png", 12, max(6, length(top_sub) * 0.22))

cat("\nAll R-side advanced figures written to", FIG, "\n")
