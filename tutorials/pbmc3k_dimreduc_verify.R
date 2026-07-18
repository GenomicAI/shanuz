#!/usr/bin/env Rscript
# R Seurat reference for the dimensional-reduction extras tutorial
# (dimreduc_vignette.md).
#
# Mirrors pbmc3k_dimreduc_tutorial.py: builds the PBMC 3k object from the same
# 10x bytes, restricted to the SAME cells and the SAME variable features the
# Python run used, then runs JackStraw / ScoreJackStraw / RunICA / RunTSNE and
# writes the per-PC, per-feature and per-cell references the Python tutorial's
# report_concordance() reads.
#
# Why the shared cell + feature lists: JackStraw's null is built from the scaled
# matrix and the PCA basis, so if the two tools disagree about which cells or
# genes are in play, nothing downstream is interpretable. Reading Python's lists
# back pins both, leaving PCA numerics as the only pre-JackStraw difference —
# and r_pca.csv exists so that can be checked directly (step 0 of the report).
#
# Writes, into tutorials/figures_dimreduc/:
#   * r_pca.csv            cell x PC1..20 embeddings      (the shared-basis check)
#   * r_jackstraw_p.csv    feature x PC1..20 empirical p  (the per-feature null)
#   * r_jackstraw_pcs.csv  PC, R_Score, R_n_sig_features  (ScoreJackStraw)
#   * r_ica.csv            cell x ICA1..20                (order/sign arbitrary)
#   * r_tsne.csv           cell x tSNE_1, tSNE_2          (coordinates arbitrary)
#   * r_*.png              R-side figures for the side-by-side tables
#
# Data: reads the same cache the Python tutorial downloads to. Run
#   python tutorials/pbmc3k_dimreduc_tutorial.py   # downloads ~8 MB, writes the lists
# then
#   Rscript tutorials/pbmc3k_dimreduc_verify.R
# Override the data folder with the PBMC3K_DATA environment variable.
#
# SLOW: JackStraw re-runs a full PCA for every one of the 100 replicates, so
# this step alone takes a few minutes. Set JS_REPLICATES to trade accuracy for
# time when iterating.
#
# Needs: Seurat, ggplot2, ica, Rtsne.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2)
})
set.seed(42)

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
HERE <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG  <- file.path(HERE, "figures_dimreduc")
DATA <- Sys.getenv("PBMC3K_DATA",
                   path.expand("~/.shanuz_data/pbmc3k/filtered_gene_bc_matrices/hg19"))
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# Mirror the Python constants (pbmc3k_dimreduc_tutorial.py).
N_PCS        <- 50
JS_DIMS      <- 20
JS_REPLICATES <- as.integer(Sys.getenv("JS_REPLICATES", "100"))
JS_PROP_FREQ <- 0.01
SCORE_THRESH <- 1e-5
N_ICS        <- 20
TSNE_DIMS    <- 10

CELLS_TXT <- file.path(FIG, "cells.txt")
HVG_TXT   <- file.path(FIG, "hvg_features.txt")
if (!file.exists(file.path(DATA, "matrix.mtx")))
  stop("PBMC 3k data not found at ", DATA,
       "\nRun `python tutorials/pbmc3k_dimreduc_tutorial.py` first.")
if (!all(file.exists(CELLS_TXT, HVG_TXT)))
  stop("cell / feature lists not found in ", FIG,
       "\nRun `python tutorials/pbmc3k_dimreduc_tutorial.py` first (it writes them).")

# ---- 1. Object on the SHARED cell + feature basis ----------------------------
cells <- readLines(CELLS_TXT)
hvg   <- readLines(HVG_TXT)
# Read10X() rewrites underscores in feature names to dashes ("Y_RNA" -> "Y-RNA"),
# so Python's symbols have to be mapped through the same rule before they can
# index the R object. The Python side normalises the same way when it reads the
# per-feature matrix back, so the two still line up by name.
hvg <- gsub("_", "-", hvg)
cat(sprintf("Shared basis from Python: %d cells, %d variable features\n",
            length(cells), length(hvg)))

pbmc <- CreateSeuratObject(Read10X(DATA), project = "pbmc3k_dimreduc",
                           min.cells = 3, min.features = 200)
missing <- setdiff(cells, colnames(pbmc))
if (length(missing))
  stop(length(missing), " of Python's cells are absent from the R object — the ",
       "two loaders disagree about QC; fix that before reading anything below.")
missing_f <- setdiff(hvg, rownames(pbmc))
if (length(missing_f))
  stop(length(missing_f), " of Python's variable features are absent from the R ",
       "object (e.g. ", paste(head(missing_f, 3), collapse = ", "), ")")
pbmc <- subset(pbmc, cells = cells)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize",
                      scale.factor = 10000, verbose = FALSE)
VariableFeatures(pbmc) <- hvg                 # Python's list, not a fresh vst run
pbmc <- ScaleData(pbmc, features = hvg, verbose = FALSE)
pbmc <- RunPCA(pbmc, features = hvg, npcs = N_PCS, verbose = FALSE)
cat(sprintf("PCA: %d cells x %d PCs\n", ncol(pbmc), ncol(Embeddings(pbmc, "pca"))))

emb <- Embeddings(pbmc, "pca")[, 1:JS_DIMS, drop = FALSE]
write.csv(data.frame(cell = rownames(emb), emb, check.names = FALSE),
          file.path(FIG, "r_pca.csv"), row.names = FALSE)

# ---- 2. JackStraw + ScoreJackStraw ------------------------------------------
# The step this tutorial exists for. R's JackRandom permutes prop.freq of the
# rows and re-runs a full PCA per replicate, taking the null loadings from that
# refit basis; shanuz projects the permuted rows onto the fixed original basis.
# Both matrices are written out so the difference (if any) can be located.
cat(sprintf("JackStraw: %d replicates x %d dims (re-runs a PCA per replicate) ...\n",
            JS_REPLICATES, JS_DIMS))
t0 <- Sys.time()
pbmc <- JackStraw(pbmc, reduction = "pca", dims = JS_DIMS,
                  num.replicate = JS_REPLICATES, prop.freq = JS_PROP_FREQ,
                  verbose = FALSE)
pbmc <- ScoreJackStraw(pbmc, dims = 1:JS_DIMS, score.thresh = SCORE_THRESH)
cat(sprintf("  done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

emp <- JS(object = pbmc[["pca"]], slot = "empirical.p.values")
write.csv(data.frame(feature = rownames(emp), emp, check.names = FALSE),
          file.path(FIG, "r_jackstraw_p.csv"), row.names = FALSE)

overall <- JS(object = pbmc[["pca"]], slot = "overall")
score <- as.numeric(overall[, "Score"])
n_sig <- apply(emp[, 1:JS_DIMS, drop = FALSE], 2,
               function(p) sum(p <= SCORE_THRESH))
write.csv(data.frame(PC = 1:JS_DIMS, R_Score = score,
                     R_n_sig_features = as.integer(n_sig)),
          file.path(FIG, "r_jackstraw_pcs.csv"), row.names = FALSE)
cat("  significant PCs (Score <= 0.05):",
    paste(which(score <= 0.05), collapse = ", "), "\n")

# ---- 3. ICA + t-SNE ----------------------------------------------------------
# Neither is compared coordinate-wise (ICA components are sign- and
# order-arbitrary; t-SNE coordinates differ by implementation), so these are
# written raw and matched/scored on the Python side.
pbmc <- RunICA(pbmc, features = hvg, nics = N_ICS, verbose = FALSE, seed.use = 42)
ica <- Embeddings(pbmc, "ica")
write.csv(data.frame(cell = rownames(ica), ica, check.names = FALSE),
          file.path(FIG, "r_ica.csv"), row.names = FALSE)

pbmc <- RunTSNE(pbmc, dims = 1:TSNE_DIMS, seed.use = 42)
tsne <- Embeddings(pbmc, "tsne")
write.csv(data.frame(cell = rownames(tsne), tsne, check.names = FALSE),
          file.path(FIG, "r_tsne.csv"), row.names = FALSE)
cat(sprintf("ICA %d comps, t-SNE on PC 1-%d written\n", ncol(ica), TSNE_DIMS))

# ---- 4. Figures (r_ prefix, side by side with the Python ones) ---------------
# Guarded so a plotting hiccup never costs the already-written references.
try({
  js_df <- data.frame(PC = 1:JS_DIMS, Score = score)
  ggsave(file.path(FIG, "r_01_jackstraw_scores.png"),
         ggplot(js_df, aes(PC, -log10(pmax(Score, 1e-300)))) +
           geom_point(size = 2) + geom_line(alpha = 0.4) +
           geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                      colour = "grey40") +
           labs(y = expression(-log[10](Score)),
                title = "R Seurat — ScoreJackStraw per PC") +
           theme_bw(),
         width = 6, height = 4, dpi = 150)

  ggsave(file.path(FIG, "r_02_elbow.png"),
         ElbowPlot(pbmc, ndims = 30) +
           labs(title = "R Seurat — elbow plot") + theme_bw(),
         width = 6, height = 4, dpi = 150)

  # Coloured by a marker rather than by cluster: cluster *labels* are arbitrary
  # and would not correspond across the two tools, whereas LYZ expression is the
  # same number on the same cell either side — so the two panels are directly
  # comparable even though the coordinates are not.
  ggsave(file.path(FIG, "r_03_tsne.png"),
         FeaturePlot(pbmc, features = "LYZ", reduction = "tsne") +
           labs(title = "R Seurat — t-SNE (PC 1-10), LYZ") + theme_bw(),
         width = 6, height = 5, dpi = 150)

  ica_df <- as.data.frame(ica[, 1:2])
  names(ica_df) <- c("IC1", "IC2")
  ggsave(file.path(FIG, "r_04_ica.png"),
         ggplot(ica_df, aes(IC1, IC2)) +
           geom_point(size = 0.4, alpha = 0.5) +
           labs(title = "R Seurat — ICA components 1 & 2") + theme_bw(),
         width = 6, height = 5, dpi = 150)
}, silent = FALSE)

cat("\nDONE — wrote references + figures to", FIG, "\n")
