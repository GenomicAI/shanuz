#!/usr/bin/env Rscript
# R Seurat reference for the one PBMC 3k figure the published Seurat vignette
# does not include: the RidgePlot (pbmc3k_tutorial.md, final step). Every other
# R figure in that tutorial links the canonical satijalab.org vignette image;
# only the RidgePlot had no R counterpart, so this script generates it.
#
# Runs the standard PBMC 3k workflow (mirrors pbmc3k_tutorial.py) and writes
#   tutorials/figures/r_12_ridge_plot.png
#
# Data: the same cache the Python tutorial downloads to. Run
#   python tutorials/pbmc3k_tutorial.py     # downloads ~24 MB first
# then
#   Rscript tutorials/pbmc3k_verify.R
# Needs: Seurat, ggplot2, ggridges.
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
set.seed(0)

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
HERE <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG  <- file.path(HERE, "figures")
DATA <- Sys.getenv("PBMC3K_DATA",
                   path.expand("~/.shanuz_data/pbmc3k/filtered_gene_bc_matrices/hg19"))
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(DATA, "matrix.mtx")))
  stop("PBMC 3k data not found at ", DATA,
       "\nRun `python tutorials/pbmc3k_tutorial.py` first.")

# ---- standard workflow (mirrors pbmc3k_tutorial.py) -------------------------
pbmc <- CreateSeuratObject(Read10X(DATA), project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pbmc <- ScaleData(pbmc, features = rownames(pbmc), verbose = FALSE)
pbmc <- RunPCA(pbmc, npcs = 50, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:10, k.param = 20, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, algorithm = 1, verbose = FALSE)
cat(sprintf("PBMC 3k: %d cells, %d clusters\n", ncol(pbmc), length(levels(pbmc))))

# ---- RidgePlot (LYZ / NKG7 / MS4A1 / CD8A across clusters) -------------------
p <- RidgePlot(pbmc, features = c("LYZ", "NKG7", "MS4A1", "CD8A"), ncol = 2) &
  theme(plot.title = element_text(size = 11))
ggsave(file.path(FIG, "r_12_ridge_plot.png"), p, width = 11, height = 8, dpi = 150, bg = "white")
cat("wrote", file.path(FIG, "r_12_ridge_plot.png"), "\n")
