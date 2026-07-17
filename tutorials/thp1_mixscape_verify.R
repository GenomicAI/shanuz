#!/usr/bin/env Rscript
# R Seurat reference for the Mixscape tutorial (mixscape_vignette.md).
#
# Mirrors thp1_mixscape_tutorial.py: builds the RNA object from the same GEO
# bytes, computes the perturbation signature (CalcPerturbSig), classifies each
# guide's cells into KO / NP (RunMixscape), and fits the guide-separating LDA
# (MixscapeLDA) — the references for shanuz.calc_perturb_sig / run_mixscape /
# mixscape_lda. Writes the per-cell calls (r_calls.csv) that the Python
# tutorial's report_concordance() reads, plus the R-side figures for the
# side-by-side tables into tutorials/figures_mixscape/:
#   * r_01_perturb_score.png  r_02_lda.png  r_03_heatmap.png
#
# To keep the perturbation signature on ONE shared gene basis across the two
# tools, this reads the variable-feature list the Python run selected
# (figures_mixscape/hvg_features.txt) instead of running FindVariableFeatures
# here — so the only divergences left are the genuinely method-level ones
# (PCA numerics, kNN ties, the DE test, the EM mixture).
#
# Data: reads the same cache the Python tutorial downloads to. Run
#   python tutorials/thp1_mixscape_tutorial.py   # downloads ~66 MB, writes HVGs
# then
#   Rscript tutorials/thp1_mixscape_verify.R
# Override the data folder with the MIXSCAPE_DATA environment variable.
#
# Needs: Seurat, ggplot2, data.table, Matrix.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(data.table); library(Matrix)
})
set.seed(42)

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
HERE <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG  <- file.path(HERE, "figures_mixscape")
DATA <- Sys.getenv("MIXSCAPE_DATA", path.expand("~/.shanuz_data/thp1_eccite"))
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

RNA_TSV  <- file.path(DATA, "GSM4633614_ECCITE_cDNA_counts.tsv.gz")
META_TSV <- file.path(DATA, "GSE153056_ECCITE_metadata.tsv.gz")
HVG_TXT  <- file.path(FIG, "hvg_features.txt")
if (!file.exists(RNA_TSV))
  stop("data not found at ", DATA,
       "\nRun `python tutorials/thp1_mixscape_tutorial.py` first.")
if (!file.exists(HVG_TXT))
  stop("hvg_features.txt not found in ", FIG,
       "\nRun `python tutorials/thp1_mixscape_tutorial.py` first (it writes it).")

# ---- 1. Load counts + metadata (the same GEO bytes Python reads) -------------
# The cDNA table is a dense ~18.6k gene x 20.7k cell matrix; fread reads it fast,
# then it is sparsified immediately so only one dense copy is ever held.
cat("Reading cDNA counts ...\n")
dt <- data.table::fread(cmd = paste("gzcat", shQuote(RNA_TSV)), showProgress = FALSE)
genes <- dt[[1]]
mat <- as.matrix(dt[, -1])
rownames(mat) <- genes
rm(dt); invisible(gc())
counts <- as(Matrix(mat, sparse = TRUE), "CsparseMatrix")
rm(mat); invisible(gc())

meta <- read.delim(gzfile(META_TSV), row.names = 1, check.names = FALSE)
common <- intersect(colnames(counts), rownames(meta))       # same cell set as Python
counts <- counts[, common]
meta   <- meta[common, ]
cat(sprintf("RNA %d genes x %d cells\n", nrow(counts), ncol(counts)))

obj <- CreateSeuratObject(counts = counts, min.cells = 3, meta.data = meta)

# ---- 2. Prep to PCA on the SHARED variable-feature basis ---------------------
hvg <- readLines(HVG_TXT)
hvg <- hvg[hvg %in% rownames(obj)]
cat(sprintf("Using %d shared variable features from Python\n", length(hvg)))
obj <- NormalizeData(obj, verbose = FALSE)
VariableFeatures(obj) <- hvg
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = 50, verbose = FALSE)

# ---- 3. Perturbation signature + Mixscape classification --------------------
obj <- CalcPerturbSig(
  object = obj, assay = "RNA", features = hvg, slot = "data",
  gd.class = "gene", nt.cell.class = "NT", reduction = "pca",
  ndims = 40, num.neighbors = 20, split.by = "replicate",
  new.assay.name = "PRTB"
)
DefaultAssay(obj) <- "PRTB"
VariableFeatures(obj) <- hvg
obj <- ScaleData(obj, features = hvg, do.scale = FALSE, do.center = TRUE, verbose = FALSE)
obj <- RunMixscape(
  object = obj, assay = "PRTB", slot = "scale.data", labels = "gene",
  nt.class.name = "NT", de.assay = "RNA", min.de.genes = 5, min.cells = 5,
  logfc.threshold = 0.25, iter.num = 10, prtb.type = "KO", verbose = FALSE
)

gclass <- as.character(obj$mixscape_class.global)
mclass <- as.character(obj$mixscape_class)
cat("\nMixscape global:\n"); print(table(gclass))
cat("\nPer-gene KO rate:\n")
ko_tbl <- as.data.frame.matrix(table(obj$gene, gclass))
ko_tbl$ko_rate <- round(ko_tbl$KO / (ko_tbl$KO + ko_tbl$NP), 3)
print(ko_tbl[order(-ko_tbl$ko_rate), ])

# ---- 4. Per-cell calls for the Python concordance report (write first) ------
calls <- data.frame(
  cell              = colnames(obj),
  R_mixscape_global = gclass,
  R_mixscape_class  = mclass,
  stringsAsFactors  = FALSE
)
write.csv(calls, file.path(FIG, "r_calls.csv"), row.names = FALSE)
cat("\nWrote r_calls.csv\n")

# ---- 5. Figures (r_ prefix, side by side with the Python ones) --------------
# Guarded so a plotting hiccup never costs the already-written calls.
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "mixscape_class"

try({
  ggsave(file.path(FIG, "r_01_perturb_score.png"),
         PlotPerturbScore(obj, target.gene.ident = "IFNGR2",
                          mixscape.class = "mixscape_class", col = "orange2") +
           ggtitle("Perturbation score — IFNGR2"),
         width = 6, height = 5, dpi = 150)
}, silent = FALSE)

try({
  ggsave(file.path(FIG, "r_03_heatmap.png"),
         MixscapeHeatmap(obj, ident.1 = "NT", ident.2 = "IFNGR2 KO",
                         balanced = TRUE, assay = "RNA", max.genes = 20,
                         angle = 0, size = 6) + NoLegend(),
         width = 7, height = 6, dpi = 150)
}, silent = FALSE)

try({
  obj <- MixscapeLDA(obj, assay = "RNA", pc.assay = "PRTB", labels = "gene",
                     nt.label = "NT", npcs = 10, logfc.threshold = 0.25,
                     verbose = FALSE)
  lda_key <- if ("lda" %in% Reductions(obj)) "lda" else Reductions(obj)[length(Reductions(obj))]
  ggsave(file.path(FIG, "r_02_lda.png"),
         DimPlot(obj, reduction = lda_key, group.by = "mixscape_class.global",
                 pt.size = 0.3) + ggtitle("MixscapeLDA (global class)"),
         width = 7, height = 5, dpi = 150)
}, silent = FALSE)

cat("\nDONE — wrote r_calls.csv + figures to", FIG, "\n")
