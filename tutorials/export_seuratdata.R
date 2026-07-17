#!/usr/bin/env Rscript
# Wave 0 data bridge: export a SeuratData dataset to a language-neutral 10x-style
# matrix folder that Python's shanuz.datasets can read, guaranteeing R and Python
# analyse *byte-identical* counts.
#
# Needed only for the two Wave-1 datasets with no clean raw source: ifnb and
# panc8 are curated SeuratData `.rda` objects (R binaries). Every other tutorial
# dataset is fetched from its original GEO/10x files by shanuz.datasets directly.
#
# This is the mirror image of the existing verify scripts, where R depends on the
# Python download (e.g. pbmc3k_verify.R needs `python pbmc3k_tutorial.py` first).
# Here the Python tutorial depends on this one-time R export.
#
# Usage:
#   Rscript tutorials/export_seuratdata.R ifnb
#   Rscript tutorials/export_seuratdata.R panc8
# Writes into ~/.shanuz_data/<name>/ :
#   matrix.mtx  features.tsv  barcodes.tsv   (read by shanuz.io.read_10x, v3 plain)
#   metadata.csv                             (per-cell labels: stim/tech/celltype…)
# Needs: Seurat, SeuratData, Matrix.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(Matrix)
})
# The SeuratData packages are large (ifnb ~394 MB, panc8 ~117 MB); R's default
# 60s download timeout is far too short for InstallData's first-time fetch.
options(timeout = 1800)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop("usage: Rscript export_seuratdata.R <ifnb|panc8> [out_dir]")
name <- args[[1]]
out  <- if (length(args) >= 2) args[[2]] else file.path(path.expand("~/.shanuz_data"), name)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

log <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "-", ..., "\n")

# --- load (installing the dataset package on first use) ----------------------
installed <- tryCatch(rownames(InstalledData()), error = function(e) character(0))
if (!name %in% installed) {
  log("InstallData('", name, "') — first-time download ...", sep = "")
  InstallData(name)
}
log("LoadData('", name, "') ...", sep = "")
obj <- LoadData(name)
obj <- UpdateSeuratObject(obj)
DefaultAssay(obj) <- "RNA"

# Seurat v5 may split counts across per-batch layers; join them back so we
# export a single counts matrix in one cell order.
obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
counts <- tryCatch(
  SeuratObject::LayerData(obj, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
)
counts <- as(counts, "CsparseMatrix")           # dgCMatrix (features x cells)
stopifnot(nrow(counts) == nrow(obj), ncol(counts) == ncol(obj))
log(sprintf("counts: %d genes x %d cells (nnz=%d)",
            nrow(counts), ncol(counts), length(counts@x)))

# --- write the 10x-style trio (gzipped v3 layout) + metadata ----------------
writeMM(counts, file.path(out, "matrix.mtx"))
# features.tsv: id <tab> symbol <tab> type. These datasets carry only symbols,
# so id == symbol; shanuz.io.read_10x(var_names="gene_symbols") reads column 1.
feats <- data.frame(id = rownames(counts), symbol = rownames(counts),
                    type = "Gene Expression")
write.table(feats, file.path(out, "features.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeLines(colnames(counts), file.path(out, "barcodes.tsv"))
# Gzip the trio to the 10x v3 layout read_10x detects (matrix.mtx.gz +
# features.tsv.gz + barcodes.tsv.gz). An ASCII MatrixMarket file of tens of
# millions of nonzeros is ~900 MB; gzip cuts it to ~150 MB. All three must be
# gzipped together — read_10x keys the .gz layout off features.tsv.gz.
for (f in c("matrix.mtx", "features.tsv", "barcodes.tsv"))
  system2("gzip", c("-f", shQuote(file.path(out, f))))
# metadata: per-cell labels the tutorials group by (stim / tech / celltype / …).
write.csv(obj@meta.data, file.path(out, "metadata.csv"), row.names = TRUE)

log("wrote:", normalizePath(out))
log("  meta columns:", paste(colnames(obj@meta.data), collapse = ", "))
log("DONE")
