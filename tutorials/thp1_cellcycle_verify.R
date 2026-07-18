#!/usr/bin/env Rscript
# R Seurat reference for the cell-cycle & module-score tutorial (cellcycle_vignette.md).
#
# Mirrors thp1_cellcycle_tutorial.py: builds the RNA object from the same GEO
# bytes, scores the cell cycle (CellCycleScoring) and an interferon program
# (AddModuleScore) on the SAME resolved gene lists the Python run wrote, and
# writes the per-cell calls (r_calls.csv) the Python tutorial's
# report_concordance() reads — the references for shanuz.cell_cycle_scoring /
# add_module_score. Also writes the R-side figures.
#
# Both AddModuleScore and CellCycleScoring sample control genes at random, and
# R's RNG is not NumPy's, so the per-cell SCORES are not expected to be identical
# across tools — only to correlate tightly, with the discrete Phase concordant.
# To rule out gene-list drift as a second source of divergence, this reads the
# exact S / G2M / IFN gene symbols the Python run resolved against the assay
# (figures_cellcycle/*.txt) rather than Seurat's cc.genes.updated.2019 directly.
#
# Data: reads the same cache the Python tutorial (and the Mixscape tutorial)
# download to. Run
#   python tutorials/thp1_cellcycle_tutorial.py   # downloads ~66 MB, writes gene lists
# then
#   Rscript tutorials/thp1_cellcycle_verify.R
# Override the data folder with the CELLCYCLE_DATA environment variable.
#
# Needs: Seurat, ggplot2, data.table, Matrix.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(data.table); library(Matrix)
})
set.seed(42)

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
HERE <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG  <- file.path(HERE, "figures_cellcycle")
DATA <- Sys.getenv("CELLCYCLE_DATA", path.expand("~/.shanuz_data/thp1_eccite"))
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

RNA_TSV  <- file.path(DATA, "GSM4633614_ECCITE_cDNA_counts.tsv.gz")
META_TSV <- file.path(DATA, "GSE153056_ECCITE_metadata.tsv.gz")
S_TXT    <- file.path(FIG, "s_genes.txt")
G2M_TXT  <- file.path(FIG, "g2m_genes.txt")
IFN_TXT  <- file.path(FIG, "ifn_genes.txt")
if (!file.exists(RNA_TSV))
  stop("data not found at ", DATA,
       "\nRun `python tutorials/thp1_cellcycle_tutorial.py` first.")
if (!all(file.exists(S_TXT, G2M_TXT, IFN_TXT)))
  stop("gene lists not found in ", FIG,
       "\nRun `python tutorials/thp1_cellcycle_tutorial.py` first (it writes them).")

# ---- 1. Load counts + metadata (the same GEO bytes Python reads) -------------
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
obj <- NormalizeData(obj, verbose = FALSE)

# ---- 2. Score on the SHARED resolved gene lists from Python -------------------
s_genes   <- readLines(S_TXT)
g2m_genes <- readLines(G2M_TXT)
ifn_genes <- readLines(IFN_TXT)
cat(sprintf("Using %d S + %d G2M + %d IFN genes from Python\n",
            length(s_genes), length(g2m_genes), length(ifn_genes)))

obj <- CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes,
                        set.ident = FALSE)
obj <- AddModuleScore(obj, features = list(ifn_genes), name = "IFN.Response",
                      seed = 1)
# AddModuleScore appends the program index to `name`; a single program -> "...1".
ifn_col <- if ("IFN.Response1" %in% colnames(obj[[]])) "IFN.Response1" else "IFN.Response"

# ---- 3. Per-cell calls for the Python concordance report (write first) --------
calls <- data.frame(
  cell        = colnames(obj),
  R_Phase     = as.character(obj$Phase),
  R_S_Score   = as.numeric(obj$S.Score),
  R_G2M_Score = as.numeric(obj$G2M.Score),
  R_IFN       = as.numeric(obj[[ifn_col]][, 1]),
  stringsAsFactors = FALSE
)
write.csv(calls, file.path(FIG, "r_calls.csv"), row.names = FALSE)
cat("\nWrote r_calls.csv\n")

cat("  R phase distribution:\n")
for (p in c("G1", "S", "G2M"))
  cat(sprintf("    %-4s %5d  (%.3f)\n", p, sum(calls$R_Phase == p),
              mean(calls$R_Phase == p)))

# ---- 4. Figures (r_ prefix, side by side with the Python ones) ----------------
# Guarded so a plotting hiccup never costs the already-written calls. The classic
# cell-cycle read: S vs G2M score, coloured by the assigned phase.
try({
  ggsave(file.path(FIG, "r_01_score_scatter.png"),
         ggplot(calls, aes(R_S_Score, R_G2M_Score, colour = R_Phase)) +
           geom_point(size = 0.4, alpha = 0.6) +
           labs(x = "S.Score", y = "G2M.Score", colour = "Phase",
                title = "R Seurat — cell-cycle scores by phase") +
           theme_bw(),
         width = 6, height = 5, dpi = 150)

  prop <- as.data.frame(table(Phase = factor(calls$R_Phase, c("G1", "S", "G2M"))))
  ggsave(file.path(FIG, "r_02_phase_bar.png"),
         ggplot(prop, aes(Phase, Freq, fill = Phase)) +
           geom_col() +
           labs(y = "cells", title = "R Seurat — phase distribution") +
           theme_bw(),
         width = 5, height = 4, dpi = 150)
}, silent = FALSE)

cat("\nDONE — wrote r_calls.csv + figures to", FIG, "\n")
