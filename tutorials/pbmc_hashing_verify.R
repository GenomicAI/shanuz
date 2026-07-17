#!/usr/bin/env Rscript
# R Seurat reference for the cell-hashing tutorial (hashing_vignette.md).
#
# Mirrors pbmc_hashing_tutorial.py: attaches the hashtag counts as an "HTO"
# assay, CLR-normalises them (margin 1), and runs HTODemux + MULTIseqDemux —
# the references for shanuz.hto_demux and shanuz.multiseq_demux. Writes the
# per-cell calls (r_calls.csv) that the Python tutorial's report_concordance()
# reads, plus the R-side figures for the side-by-side tables into
# tutorials/figures_hashing/:
#   * r_01_ridge.png  r_02_scatter.png  r_03_ncount_violin.png
#
# Data: reads the same cache the Python tutorial downloads to. Run
#   python tutorials/pbmc_hashing_tutorial.py   # downloads ~34 MB first
# then
#   Rscript tutorials/pbmc_hashing_verify.R
# Override the data folder with the HASHING_DATA environment variable.
#
# Needs: Seurat, ggplot2.
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
set.seed(42)

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
HERE <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG  <- file.path(HERE, "figures_hashing")
DATA <- Sys.getenv("HASHING_DATA", path.expand("~/.shanuz_data/pbmc_hashing"))
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
RNA_TSV <- file.path(DATA, "GSM2895282_Hashtag-RNA.umi.txt.gz")
HTO_CSV <- file.path(DATA, "GSM2895283_Hashtag-HTO-count.csv.gz")
if (!file.exists(HTO_CSV))
  stop("hashing data not found at ", DATA,
       "\nRun `python tutorials/pbmc_hashing_tutorial.py` first.")

# ---- 1. Load the HTO matrix (RNA only supplies the barcode set) -------------
# HTODemux needs only the HTO assay. The RNA text is ~41k genes x 50k cells,
# dense (~13 GB if read in full), so read just its header for the barcode list
# the Python loader intersects the HTO matrix down to — same cell set, cheaply.
hdr <- readLines(gzfile(RNA_TSV), n = 1)
rna_cells <- strsplit(hdr, "\t")[[1]][-1]              # drop the "GENE" index corner

hto <- read.csv(gzfile(HTO_CSV), row.names = 1, check.names = FALSE)
hto <- hto[!rownames(hto) %in% c("bad_struct", "no_match", "total_reads"), ]
common <- intersect(rna_cells, colnames(hto))         # RNA order == Python cell set
hto <- as.matrix(hto[, common])
rownames(hto) <- sub("-.*$", "", rownames(hto))       # BatchA-AGG... -> BatchA (match Python)
cat(sprintf("HTO %d hashtags x %d cells\n", nrow(hto), ncol(hto)))

obj <- CreateSeuratObject(counts = as.sparse(hto), assay = "HTO")
obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR",
                     margin = 1, verbose = FALSE)

# ---- 2. HTODemux + MULTIseqDemux --------------------------------------------
obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99, verbose = FALSE)
obj <- MULTIseqDemux(obj, assay = "HTO", quantile = 0.7)

ms_glob <- ifelse(obj$MULTI_ID %in% c("Doublet", "Negative"),
                  as.character(obj$MULTI_ID), "Singlet")
cat("\nHTODemux global:\n"); print(table(obj$HTO_classification.global))
cat("\nSinglets per sample (hash.ID):\n"); print(table(obj$hash.ID))
cat("\nMULTIseqDemux global:\n"); print(table(ms_glob))
cat("\nHTODemux (rows) x MULTIseqDemux (cols):\n")
print(table(obj$HTO_classification.global, ms_glob))

# ---- 3. Per-cell calls for the Python concordance report --------------------
calls <- data.frame(
  cell         = colnames(obj),
  R_HTO_global = as.character(obj$HTO_classification.global),
  R_hash_ID    = as.character(obj$hash.ID),
  R_MULTI_ID   = as.character(obj$MULTI_ID),
  stringsAsFactors = FALSE
)
write.csv(calls, file.path(FIG, "r_calls.csv"), row.names = FALSE)

# ---- 4. Figures (r_ prefix, side by side with the Python ones) --------------
Idents(obj) <- "hash.ID"
ggsave(file.path(FIG, "r_01_ridge.png"),
       RidgePlot(obj, assay = "HTO", features = rownames(obj), ncol = 3),
       width = 11, height = 7, dpi = 150)
ggsave(file.path(FIG, "r_02_scatter.png"),
       FeatureScatter(obj, feature1 = rownames(obj)[1], feature2 = rownames(obj)[2],
                      group.by = "HTO_classification.global"),
       width = 6, height = 5, dpi = 150)
ggsave(file.path(FIG, "r_03_ncount_violin.png"),
       VlnPlot(obj, features = "nCount_HTO", group.by = "HTO_classification.global",
               log = TRUE, pt.size = 0),
       width = 6, height = 5, dpi = 150)

cat("\nDONE — wrote r_calls.csv + 3 figures to", FIG, "\n")
