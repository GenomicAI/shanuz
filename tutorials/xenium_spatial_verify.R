#!/usr/bin/env Rscript
# R Seurat reference for the Xenium spatial tutorial (xenium_spatial_tutorial.md).
#
# Runs the same pipeline as generate_spatial_plots.py on the public 10x
# mouse-brain Xenium subset and writes, into tutorials/figures_spatial/:
#   * r_reference.json  — deterministic anchors (compared by compare_xenium_anchors.py)
#   * r_*.png           — R-side figures for the side-by-side tables
#
# Data: reads the same cache the Python tutorial downloads to. Run
#   python tutorials/generate_spatial_plots.py   # downloads ~20 MB first
# then
#   Rscript tutorials/xenium_spatial_verify.R
# Override the data folder with the XENIUM_DATA environment variable.
#
# Needs: Seurat, Matrix, FNN, ggplot2, jsonlite.
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(FNN); library(ggplot2); library(jsonlite)
})
set.seed(42)

# --- resolve paths relative to this script ---------------------------------
.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
HERE <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG  <- file.path(HERE, "figures_spatial")
DATA <- Sys.getenv("XENIUM_DATA", path.expand("~/.shanuz_data/xenium_mouse_brain"))
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(DATA, "cell_feature_matrix", "matrix.mtx.gz")))
  stop("Xenium data not found at ", DATA,
       "\nRun `python tutorials/generate_spatial_plots.py` first (it downloads it).")

# Same ordered marker panels as generate_spatial_plots.py (order fixes ties).
CELLTYPE_MARKERS <- list(
  Excitatory      = c("Slc17a7","Slc17a6","Satb2","Fezf2","Neurod6"),
  Inhibitory      = c("Gad1","Gad2","Pvalb","Sst","Vip","Lamp5"),
  Astrocyte       = c("Aqp4","Gfap","Ntsr2"),
  Oligodendrocyte = c("Sox10","Opalin","Gjc3"),
  OPC             = c("Pdgfra","Cspg4"),
  Vascular        = c("Cldn5","Pecam1","Kdr","Emcn","Adgrl4"),
  Immune          = c("Cd68","Trem2","Siglech","Laptm5","Cd53")
)
FOCAL <- "Vascular"

# ---- 1. Load: Read10X (keep Gene Expression) + coordinates ------------------
raw <- Read10X(file.path(DATA, "cell_feature_matrix"))
counts <- if (is.list(raw)) raw[["Gene Expression"]] else raw
cells_csv <- read.csv(gzfile(file.path(DATA, "cells.csv.gz")))
rownames(cells_csv) <- as.character(cells_csv$cell_id)

obj <- CreateSeuratObject(counts = counts, assay = "Xenium", project = "xenium_mb")
n_cells_raw <- ncol(obj); n_genes <- nrow(obj)
cat(sprintf("Loaded %d cells x %d genes\n", n_cells_raw, n_genes))
coords <- cells_csv[colnames(obj), c("x_centroid","y_centroid")]
colnames(coords) <- c("x","y")
obj$x <- coords$x; obj$y <- coords$y

# ---- 2. QC filter ----------------------------------------------------------
obj <- subset(obj, subset = nCount_Xenium >= 10)
n_cells_qc <- ncol(obj)
cat(sprintf("QC nCount_Xenium>=10: %d cells retained\n", n_cells_qc))

# ---- 3. Deterministic marker cell types (argmax of summed raw counts) ------
cmat <- GetAssayData(obj, assay = "Xenium", layer = "counts")
score <- sapply(names(CELLTYPE_MARKERS), function(t) {
  g <- intersect(CELLTYPE_MARKERS[[t]], rownames(cmat))
  if (length(g) == 0) return(rep(0, ncol(cmat)))
  Matrix::colSums(cmat[g, , drop = FALSE])
})
best <- max.col(score, ties.method = "first")                  # == np.argmax
cell_type <- names(CELLTYPE_MARKERS)[best]
cell_type[apply(score, 1, max) == 0] <- "Other"
obj$cell_type <- cell_type
ct_counts <- table(obj$cell_type)
cat("cell_type counts:\n"); print(ct_counts)

# ---- 4. Deterministic spatial region split (median y) ----------------------
ymed <- median(obj$y)
obj$region <- ifelse(obj$y >= ymed, "ventral", "dorsal")

# ---- 5. Standard unsupervised pipeline (structure view) --------------------
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = n_genes, verbose = FALSE)
obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
obj <- RunPCA(obj, features = rownames(obj), npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.3, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
n_clusters <- length(unique(obj$seurat_clusters))
cat(sprintf("clusters (res=0.3): %d\n", n_clusters))

# ---- 6. Spatial statistics on focal type (FNN == shanuz idiom) -------------
xy <- as.matrix(obj@meta.data[, c("x","y")])
is_focal <- obj$cell_type == FOCAL
foc_xy <- xy[is_focal, , drop = FALSE]
nn <- FNN::get.knn(foc_xy, k = 1)$nn.dist[, 1]                 # nearest same-type
nn_median <- median(nn); nn_mean <- mean(nn)
knn <- FNN::get.knnx(data = xy, query = foc_xy, k = 11)$nn.index[, -1]  # drop self
prop_focal <- rowMeans(matrix(is_focal[knn], nrow = nrow(knn)))
dens_mean <- mean(prop_focal)
cat(sprintf("%s: n=%d  NN median=%.4f  local prop mean=%.5f\n",
            FOCAL, sum(is_focal), nn_median, dens_mean))

# ---- 7. Composition test across region split (mirrors composition_test) ----
tab <- table(obj$cell_type, obj$region)
ref <- "dorsal"; test <- "ventral"
n_ref <- sum(tab[, ref]); n_test <- sum(tab[, test])
comp <- data.frame(group = rownames(tab), stringsAsFactors = FALSE)
comp$log2_ratio <- NA_real_; comp$odds_ratio <- NA_real_; comp$p <- NA_real_
for (i in seq_len(nrow(tab))) {
  a <- tab[i, test]; b <- tab[i, ref]; c <- n_test - a; d <- n_ref - b
  ft <- fisher.test(matrix(c(a, c, b, d), nrow = 2))
  comp$log2_ratio[i] <- log2((a / n_test) / (b / n_ref))
  comp$odds_ratio[i] <- unname(ft$estimate)
  comp$p[i] <- ft$p.value
}
comp$padj <- p.adjust(comp$p, method = "BH")
comp$enriched_in <- ifelse(comp$log2_ratio > 0, test, ref)
comp <- comp[order(-comp$log2_ratio), ]
chisq_p <- suppressWarnings(chisq.test(tab)$p.value)
cat(sprintf("composition_test region (chi2 p=%.3g):\n", chisq_p))
print(comp[, c("group","log2_ratio","odds_ratio","padj","enriched_in")], row.names = FALSE)

# ---- 8. Niches (knn composition + kmeans; mirrors build_niche_assay) --------
knn20 <- FNN::get.knn(xy, k = 20)$nn.index
ctypes_sorted <- sort(unique(obj$cell_type))
comp_mat <- sapply(ctypes_sorted, function(t)
  rowMeans(matrix(obj$cell_type[knn20] == t, nrow = nrow(xy))))
set.seed(0)
obj$niches <- paste0("niche_", kmeans(comp_mat, centers = 6, nstart = 10)$cluster)
cat(sprintf("niches: %d\n", length(unique(obj$niches))))

# ============================ FIGURES (manual ggplot) =======================
# ImageDimPlot renders blank under ggplot2 4.x, so draw centroids directly.
theme_sp <- theme_void() + theme(legend.position = "right",
                                 plot.title = element_text(size = 10, face = "bold"))
df <- data.frame(x = obj$x, y = obj$y, cell_type = obj$cell_type,
                 cluster = factor(obj$seurat_clusters), niches = obj$niches,
                 region = obj$region, Slc17a7 = as.numeric(GetAssayData(
                   obj, assay = "Xenium", layer = "data")["Slc17a7", ]))
sp_plot <- function(colvar, title, discrete = TRUE) {
  p <- ggplot(df, aes(x = x, y = -y, color = .data[[colvar]])) +
    geom_point(size = 0.35) + coord_fixed() + ggtitle(title) + theme_sp +
    guides(color = guide_legend(override.aes = list(size = 3)))
  if (!discrete) p <- p + scale_color_viridis_c() + guides(color = "colorbar")
  p
}
sv <- function(p, name, w = 6.5, h = 5.5)
  ggsave(file.path(FIG, name), p, width = w, height = h, dpi = 130, bg = "white")
sv(sp_plot("cell_type", "R Seurat - marker cell types"), "r_03_image_celltype.png")
sv(sp_plot("cluster",   "R Seurat - clusters (res=0.3)"), "r_04_image_clusters.png")
sv(sp_plot("Slc17a7",   "R Seurat - Slc17a7", discrete = FALSE), "r_05_image_feature_Slc17a7.png", 6)
sv(sp_plot("niches",    "R Seurat - niches"), "r_06_image_niches.png")
df$is_focal <- ifelse(df$cell_type == FOCAL, FOCAL, "other")
sv(ggplot(df, aes(x = x, y = -y, color = is_focal)) + geom_point(size = 0.35) +
     coord_fixed() + scale_color_manual(values = c(Vascular = "#d62728", other = "#dddddd")) +
     ggtitle(paste0("R Seurat - ", FOCAL, " cells")) + theme_sp +
     guides(color = guide_legend(override.aes = list(size = 3))), "r_07_image_focal.png")
umap <- as.data.frame(Embeddings(obj, "umap")); colnames(umap) <- c("UMAP1","UMAP2")
umap$cell_type <- obj$cell_type
sv(ggplot(umap, aes(UMAP1, UMAP2, color = cell_type)) + geom_point(size = 0.3) +
     ggtitle("R Seurat - UMAP marker cell types") + theme_bw() +
     theme(plot.title = element_text(size = 10, face = "bold")) +
     guides(color = guide_legend(override.aes = list(size = 3))), "r_02_umap_celltype.png", 6.5, 5)
qcl <- rbind(data.frame(cell_type = obj$cell_type, metric = "nCount_Xenium", value = obj$nCount_Xenium),
             data.frame(cell_type = obj$cell_type, metric = "nFeature_Xenium", value = obj$nFeature_Xenium))
sv(ggplot(qcl, aes(cell_type, value, fill = cell_type)) + geom_violin(scale = "width") +
     facet_wrap(~metric, scales = "free_y") + theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none"),
   "r_01_qc_violin.png", 12, 4)

# ============================ ANCHORS =======================================
comp_json <- setNames(lapply(seq_len(nrow(comp)), function(i) list(
  log2_ratio = comp$log2_ratio[i], odds_ratio = comp$odds_ratio[i],
  p = comp$p[i], padj = comp$padj[i], enriched_in = comp$enriched_in[i])),
  comp$group)
anchors <- list(
  n_cells_raw = n_cells_raw, n_genes = n_genes, n_cells_qc = n_cells_qc,
  celltype_counts = as.list(setNames(as.integer(ct_counts), names(ct_counts))),
  focal_type = FOCAL, n_focal = sum(is_focal),
  focal_nn_median = nn_median, focal_nn_mean = nn_mean,
  focal_local_density_mean = dens_mean, region_ymed = ymed,
  composition_chisq_p = chisq_p, composition = comp_json,
  n_clusters = n_clusters, n_niches = length(unique(obj$niches)))
write_json(anchors, file.path(FIG, "r_reference.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 12)
cat("\nwrote", file.path(FIG, "r_reference.json"), "\n")
