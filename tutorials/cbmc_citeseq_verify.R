#!/usr/bin/env Rscript
# R Seurat reference for the multimodal CITE-seq tutorial (multimodal_citeseq.md).
#
# Mirrors cbmc_citeseq_tutorial.py / generate_multimodal_plots.py: builds the
# RNA object, attaches the CLR-normalised ADT assay, runs the RNA workflow,
# annotates clusters by surface protein (annotate_cells ported to R), and runs
# WNN joint clustering. Writes the R-side figures for the side-by-side tables
# into tutorials/figures_multimodal/:
#   * r_01_rna_umap_clusters.png ... r_10_adt_weight_by_celltype.png
#
# Data: reads the same cache the Python tutorial downloads to. Run
#   python tutorials/cbmc_citeseq_tutorial.py   # downloads ~15 MB first
# then
#   Rscript tutorials/cbmc_citeseq_verify.R
# Override the data folder with the CBMC_DATA environment variable.
#
# Needs: Seurat, ggplot2, patchwork.
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })
set.seed(0)

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
HERE <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG  <- file.path(HERE, "figures_multimodal")
DATA <- Sys.getenv("CBMC_DATA", path.expand("~/.shanuz_data/cbmc"))
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
RNA_CSV <- file.path(DATA, "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz")
ADT_CSV <- file.path(DATA, "GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz")
if (!file.exists(RNA_CSV))
  stop("CBMC data not found at ", DATA,
       "\nRun `python tutorials/cbmc_citeseq_tutorial.py` first.")

# ---- 1. Load RNA (keep human genes) + attach CLR-normalised ADT -------------
cat("Reading RNA csv (large, ~1 min)...\n")
rna_raw <- as.sparse(read.csv(gzfile(RNA_CSV), header = TRUE, row.names = 1))
hmask <- startsWith(rownames(rna_raw), "HUMAN_")          # mirror the Python loader
rna <- rna_raw[hmask, ]
rownames(rna) <- make.unique(sub("^HUMAN_", "", rownames(rna)))
adt <- as.sparse(read.csv(gzfile(ADT_CSV), header = TRUE, row.names = 1))
common <- colnames(rna)[colnames(rna) %in% colnames(adt)]  # shared barcodes, RNA order
rna <- rna[, common]; adt <- adt[, common]
cat(sprintf("RNA %d human genes x %d cells | ADT %d proteins\n",
            nrow(rna), length(common), nrow(adt)))

obj <- CreateSeuratObject(counts = rna, project = "cbmc", min.cells = 3, min.features = 0)
adt_assay <- CreateAssayObject(counts = adt[, colnames(obj)])
obj[["ADT"]] <- adt_assay
Key(obj[["ADT"]]) <- "adt_"
obj <- NormalizeData(obj, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = FALSE)
cat(sprintf("Proteins: %s\n", paste(rownames(obj[["ADT"]]), collapse = ", ")))

# ---- 2. RNA workflow --------------------------------------------------------
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:15, k.param = 20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.6, algorithm = 1, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:15, verbose = FALSE)
cat(sprintf("RNA clusters (res=0.6): %d\n", length(levels(obj))))

# per-cluster mean CLR for the marker panel (mirrors Python tutorial Step 3)
{
  ad <- GetAssayData(obj, assay = "ADT", layer = "data")
  ids <- as.character(Idents(obj)); cls <- as.character(sort(unique(as.integer(ids))))
  panel <- c("CD3","CD4","CD8","CD19","CD14","CD16","CD56","CD11c","CD34")
  tab <- sapply(cls, function(c) round(rowMeans(ad[panel, ids == c, drop = FALSE]), 2))
  rownames(tab) <- panel; colnames(tab) <- paste0("c", cls)
  cat("Per-cluster mean CLR:\n"); print(tab)
}

# ---- 3. Annotate by surface protein (annotate_cells ported from Python) -----
# Same protein-gated priority logic and the same CLR cut-offs as
# cbmc_citeseq_tutorial.py. Seurat CLR (margin=2) puts T clusters at CD3~1.5-2.3
# and non-T at ~0.3-0.75; Shanuz reproduces that transform exactly, so the
# thresholds are shared verbatim between the two scripts.
RNA_FALLBACK <- list(Platelet = c("PPBP","PF4"), Erythroid = c("HBB","HBA1"),
                     pDC = c("IGJ","PLD4","SERPINF1"), Cycling = c("STMN1","MKI67","TUBB"))
annotate_cells <- function(obj) {
  idents <- as.character(Idents(obj))
  clusters <- as.character(sort(unique(as.integer(idents))))
  adt_data <- GetAssayData(obj, assay = "ADT", layer = "data")
  prot <- rownames(adt_data)
  rna_data <- GetAssayData(obj, assay = "RNA", layer = "data")
  rna_feats <- rownames(rna_data)
  pm <- function(p, mask) if (p %in% prot) mean(adt_data[p, mask]) else -Inf
  rmean <- function(genes, mask) {
    present <- intersect(genes, rna_feats)
    if (length(present) == 0) return(0.0)
    mean(vapply(present, function(g) mean(rna_data[g, mask]), numeric(1)))
  }
  rna_fallback <- function(mask) {
    best <- "Other"; best_score <- 0.30
    for (label in names(RNA_FALLBACK)) {
      present <- intersect(RNA_FALLBACK[[label]], rna_feats)
      if (length(present) == 0) next
      score <- mean(vapply(present, function(g) mean(rna_data[g, mask]), numeric(1)))
      if (score > best_score) { best_score <- score; best <- label }
    }
    best
  }
  assignment <- character()
  for (c in clusters) {
    mask <- idents == c
    cd3 <- pm("CD3", mask); cd4 <- pm("CD4", mask); cd8 <- pm("CD8", mask)
    cd19 <- pm("CD19", mask); cd14 <- pm("CD14", mask)
    cd16 <- pm("CD16", mask); cd56 <- pm("CD56", mask)
    cd11c <- pm("CD11c", mask); cd34 <- pm("CD34", mask)
    if (rmean(c("PPBP","PF4"), mask) > 2.0)                assignment[c] <- "Platelet"
    else if (rmean(c("HBB","HBA1"), mask) > 2.5)          assignment[c] <- "Erythroid"
    else if (cd3 > 1.0)                                   assignment[c] <- if (cd8 > 1.0) "CD8 T" else "CD4 T"
    else if (cd16 > 0.8 && cd56 > 0.8)                    assignment[c] <- "NK"
    else if (cd19 > 1.5)                                  assignment[c] <- "B"
    else if (cd34 > 1.0)                                  assignment[c] <- "Progenitor"
    else if (cd14 > 0.9)                                  assignment[c] <- "CD14+ Mono"
    else if (cd11c > 1.0)                                 assignment[c] <- "DC / Mono"
    else                                                  assignment[c] <- rna_fallback(mask)
  }
  assignment
}
obj$rna_clusters <- Idents(obj)
anno <- annotate_cells(obj)
cat("Protein cell-type assignment:\n"); print(anno)
Idents(obj) <- obj$rna_clusters
obj <- RenameIdents(obj, anno)
obj$protein_celltype <- Idents(obj)

# ---- 4. Weighted Nearest Neighbor (mirrors run_wnn in the Python tutorial) ---
# The ADT panel is 13 proteins, so apca must stay under that: npcs = 12, and
# approx = FALSE because irlba needs npcs well below the feature count.
DefaultAssay(obj) <- "ADT"
VariableFeatures(obj) <- rownames(obj[["ADT"]])
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, reduction.name = "apca", reduction.key = "apca_",
              npcs = 12, approx = FALSE, verbose = FALSE)
obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "apca"),
                               dims.list = list(1:15, 1:12), k.nn = 20, verbose = FALSE)
obj <- FindClusters(obj, graph.name = "wsnn", resolution = 0.6, verbose = FALSE)
obj$wnn_clusters <- Idents(obj)
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap",
               reduction.key = "wnnUMAP_", verbose = FALSE)
DefaultAssay(obj) <- "RNA"
Idents(obj) <- obj$protein_celltype   # FindClusters left the WNN ids active
cat(sprintf("WNN clusters (res=0.6): %d | mean weight RNA %.2f | ADT %.2f\n",
            length(unique(obj$wnn_clusters)), mean(obj$RNA.weight), mean(obj$ADT.weight)))
cat("Mean ADT weight by protein cell type:\n")
print(round(sort(tapply(obj$ADT.weight, obj$protein_celltype, mean), decreasing = TRUE), 3))

# ============================ FIGURES =======================================
sv <- function(p, name, w = 8, h = 6.5)
  ggsave(file.path(FIG, name), p, width = w, height = h, dpi = 150, bg = "white")
brand <- function(p, t) p + plot_annotation(title = t,
  theme = theme(plot.title = element_text(size = 13, face = "bold")))
titled <- function(p, t) p + ggtitle(t) +
  theme(plot.title = element_text(size = 12, face = "bold"))

# 01 RNA UMAP - clusters
sv(titled(DimPlot(obj, reduction = "umap", group.by = "rna_clusters", label = TRUE),
          "R Seurat - RNA clusters"), "r_01_rna_umap_clusters.png")
# 02 RNA UMAP - protein cell types
sv(titled(DimPlot(obj, reduction = "umap", group.by = "protein_celltype", label = TRUE),
          "R Seurat - cell types (protein + RNA)"), "r_02_rna_umap_celltypes.png", 8.5, 6.5)

# 03 surface-protein feature plots on the RNA UMAP
DefaultAssay(obj) <- "ADT"
sv(brand(FeaturePlot(obj, c("CD3","CD4","CD8","CD19","CD14","CD16","CD56","CD11c"),
                     reduction = "umap", ncol = 4, min.cutoff = "q05", max.cutoff = "q95"),
         "R Seurat - surface proteins on the RNA UMAP"),
   "r_03_adt_featureplots.png", 15, 7)
DefaultAssay(obj) <- "RNA"

# 04 protein (ADT) vs RNA for the same marker
sv(brand(FeaturePlot(obj, c("adt_CD19","rna_CD19","adt_CD3","rna_CD3E",
                            "adt_CD8","rna_CD8A","adt_CD14","rna_CD14"),
                     reduction = "umap", ncol = 2, min.cutoff = "q05", max.cutoff = "q95"),
         "R Seurat - protein (left) vs RNA (right)"),
   "r_04_protein_vs_rna.png", 9, 16)

# 05 ADT ridge plots by cell type
DefaultAssay(obj) <- "ADT"
sv(brand(RidgePlot(obj, features = c("CD3","CD19","CD14","CD56"),
                   group.by = "protein_celltype", ncol = 2),
         "R Seurat - ADT ridge plots"), "r_05_adt_ridgeplots.png", 12, 9)
DefaultAssay(obj) <- "RNA"

# 06 / 07 ADT feature scatter — protein bivariates
sv(titled(FeatureScatter(obj, "adt_CD4", "adt_CD8", group.by = "protein_celltype"),
          "R Seurat - ADT CD4 vs CD8"), "r_06_adt_scatter_CD4_CD8.png", 7, 5.5)
sv(titled(FeatureScatter(obj, "adt_CD19", "adt_CD3", group.by = "protein_celltype"),
          "R Seurat - ADT CD19 vs CD3"), "r_07_adt_scatter_CD19_CD3.png", 7, 5.5)

# 08 WNN joint clusters on the joint embedding
sv(titled(DimPlot(obj, reduction = "wnn.umap", group.by = "wnn_clusters", label = TRUE),
          "R Seurat - WNN joint clusters"), "r_08_wnn_umap_clusters.png")

# 09 RNA-only vs joint embedding, same cells and labels
sv(brand((DimPlot(obj, reduction = "umap", group.by = "protein_celltype", label = TRUE) +
            ggtitle("RNA-only UMAP") + NoLegend()) |
           (DimPlot(obj, reduction = "wnn.umap", group.by = "protein_celltype", label = TRUE) +
              ggtitle("WNN joint UMAP")),
         "R Seurat - RNA-only vs WNN joint embedding"),
   "r_09_wnn_vs_rna_umap.png", 14, 6.5)

# 10 learned per-cell modality weights
sv(titled(VlnPlot(obj, features = "ADT.weight", group.by = "protein_celltype",
                  pt.size = 0, sort = FALSE) + NoLegend(),
          "R Seurat - ADT weight by cell type"), "r_10_adt_weight_by_celltype.png", 9, 5)

cat("\nAll R-side multimodal figures written to", FIG, "\n")
