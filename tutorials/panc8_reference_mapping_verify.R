#!/usr/bin/env Rscript
# R Seurat reference for the reference-mapping tutorial (refmap_vignette.md).
#
# Mirrors panc8_reference_mapping_tutorial.py: builds the reference (celseq2) and
# query (smartseq2) objects from the same exported panc8 counts, runs Seurat's
# FindTransferAnchors(reduction = "pcaproject") + TransferData to annotate the
# query, and writes the per-cell predicted labels (r_calls.csv) that the Python
# tutorial's report_concordance() compares against — the reference for shanuz's
# find_transfer_anchors / transfer_data. Also writes the R-side projection figure.
#
# Transferred labels ARE comparable per cell (both tools pick a predicted.id from
# the same reference label set), so the concordance is a direct per-cell label
# agreement plus each tool's accuracy against the query's known cell types — no
# coordinate comparison, which the projected embeddings would not support.
#
# To keep the projection on ONE shared gene basis, this reads the variable
# features the Python run selected (figures_refmap/hvg_features.txt) instead of
# running FindVariableFeatures here, so the only divergences left are the
# genuinely method-level ones (PCA numerics, the anchor/weight kernels, kNN ties).
#
# Data: the same counts the Python tutorial reads. Run
#   Rscript tutorials/export_seuratdata.R panc8               # one-time counts export
#   python  tutorials/panc8_reference_mapping_tutorial.py     # writes the shared HVGs
# then
#   Rscript tutorials/panc8_reference_mapping_verify.R
# Override the data folder with the PANC8_DATA environment variable.
#
# Needs: Seurat, ggplot2, Matrix.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(Matrix)
})
set.seed(42)

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
HERE <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG  <- file.path(HERE, "figures_refmap")
DATA <- Sys.getenv("PANC8_DATA", path.expand("~/.shanuz_data/panc8"))
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

REFERENCE_TECH <- "celseq2"
QUERY_TECH     <- "smartseq2"
CELLTYPE       <- "celltype"
N_PCS <- 30
DIMS  <- 1:N_PCS

HVG_TXT <- file.path(FIG, "hvg_features.txt")
if (!file.exists(file.path(DATA, "matrix.mtx.gz")))
  stop("panc8 counts not found at ", DATA,
       "\nRun `Rscript tutorials/export_seuratdata.R panc8` first.")
if (!file.exists(HVG_TXT))
  stop("hvg_features.txt not found in ", FIG,
       "\nRun `python tutorials/panc8_reference_mapping_tutorial.py` first (it writes it).")

# ---- 1. Load the same exported counts + metadata, split by technology --------
cat("Reading panc8 counts ...\n")
counts <- Read10X(DATA)                                   # features x cells
meta <- read.csv(file.path(DATA, "metadata.csv"), row.names = 1, check.names = FALSE)
meta <- meta[colnames(counts), , drop = FALSE]
# One object over the full gene universe, then subset — so reference and query
# share an identical feature set (min.cells filtered once, like the Python run).
full <- CreateSeuratObject(counts = counts, min.cells = 3, meta.data = meta)
reference <- subset(full, subset = tech == REFERENCE_TECH)
query     <- subset(full, subset = tech == QUERY_TECH)
cat(sprintf("reference %s %d cells | query %s %d cells | %d cell types\n",
            REFERENCE_TECH, ncol(reference), QUERY_TECH, ncol(query),
            length(unique(query[[CELLTYPE]][, 1]))))

# ---- 2. Prep on the SHARED variable-feature basis from Python -----------------
hvg <- readLines(HVG_TXT)
hvg <- hvg[hvg %in% rownames(full)]
cat(sprintf("Using %d shared variable features from Python\n", length(hvg)))

reference <- NormalizeData(reference, verbose = FALSE)
VariableFeatures(reference) <- hvg
reference <- ScaleData(reference, features = hvg, verbose = FALSE)
reference <- RunPCA(reference, features = hvg, npcs = N_PCS, verbose = FALSE)

query <- NormalizeData(query, verbose = FALSE)
query <- ScaleData(query, features = hvg, verbose = FALSE)

# ---- 3. Transfer anchors (project the query into the reference PCA) + labels --
cat("FindTransferAnchors (pcaproject) ...\n")
anchors <- FindTransferAnchors(
  reference = reference, query = query, reference.reduction = "pca",
  reduction = "pcaproject", features = hvg, dims = DIMS, k.filter = 200,
  verbose = FALSE
)
cat(sprintf("  %d transfer anchors\n", nrow(slot(anchors, "anchors"))))

predictions <- TransferData(
  anchorset = anchors, refdata = reference[[CELLTYPE]][, 1], dims = DIMS,
  k.weight = 50, verbose = FALSE
)

# ---- 4. Per-cell calls for the Python concordance report (write first) --------
truth <- as.character(query[[CELLTYPE]][, 1])
calls <- data.frame(
  cell        = colnames(query),
  tech        = as.character(query$tech),
  celltype    = truth,
  R_predicted = as.character(predictions$predicted.id),
  R_score_max = as.numeric(predictions$prediction.score.max),
  stringsAsFactors = FALSE
)
write.csv(calls, file.path(FIG, "r_calls.csv"), row.names = FALSE)
cat("\nWrote r_calls.csv\n")

acc <- mean(calls$R_predicted == calls$celltype)
cat(sprintf("  R transfer accuracy: %.4f (%d/%d query cells)\n",
            acc, sum(calls$R_predicted == calls$celltype), nrow(calls)))
cat("  R per-class recall:\n")
for (ct in names(sort(table(truth), decreasing = TRUE))) {
  m <- truth == ct
  cat(sprintf("    %-20s support=%4d  recall=%.4f\n", ct, sum(m),
              mean(calls$R_predicted[m] == ct)))
}

# ---- 5. Figures (r_ prefix, side by side with the Python ones) ----------------
# Guarded so a plotting/version hiccup never costs the already-written calls.
# MapQuery projects the query into the reference's UMAP (Seurat's ProjectUMAP).
try({
  reference <- RunUMAP(reference, dims = DIMS, reduction = "pca",
                       return.model = TRUE, verbose = FALSE)
  ggsave(file.path(FIG, "r_01_reference_umap_celltype.png"),
         DimPlot(reference, reduction = "umap", group.by = CELLTYPE, label = TRUE) +
           ggtitle(sprintf("Reference (%s) — by cell type", REFERENCE_TECH)),
         width = 7, height = 5, dpi = 150)

  # Carry the TransferData predicted.id onto the query, then MapQuery purely for
  # the projection (refdata = NULL), so the figure colours by the same labels
  # written to r_calls.csv rather than relying on MapQuery's transfer naming.
  query$predicted_id <- as.character(predictions[colnames(query), "predicted.id"])
  query <- MapQuery(anchorset = anchors, reference = reference, query = query,
                    refdata = NULL, reference.reduction = "pca",
                    reduction.model = "umap")
  ggsave(file.path(FIG, "r_02_query_projected_predicted.png"),
         DimPlot(query, reduction = "ref.umap", group.by = "predicted_id", label = TRUE) +
           ggtitle(sprintf("Query (%s) projected — predicted labels", QUERY_TECH)),
         width = 7, height = 5, dpi = 150)
  ggsave(file.path(FIG, "r_03_query_projected_truth.png"),
         DimPlot(query, reduction = "ref.umap", group.by = CELLTYPE, label = TRUE) +
           ggtitle(sprintf("Query (%s) projected — true labels", QUERY_TECH)),
         width = 7, height = 5, dpi = 150)
}, silent = FALSE)

cat("\nDONE — wrote r_calls.csv + figures to", FIG, "\n")
