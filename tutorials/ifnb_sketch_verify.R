#!/usr/bin/env Rscript
# R Seurat reference for the leverage-score sketching tutorial (sketch_vignette.md).
#
# Mirrors ifnb_sketch_tutorial.py: builds the ifnb object from the same exported
# counts, restricted to the SAME cells and scored on the SAME variable features
# the Python run used, then runs LeverageScore in both of Seurat's regimes,
# SketchData (leverage and uniform) and ProjectData, and writes the per-cell
# references the Python tutorial's report_concordance() reads.
#
# Why the shared cell + feature lists: leverage is a property of a matrix. If the
# two tools disagree about which cells or genes are in play, every number below
# is uninterpretable.
#
# ONE GOTCHA WORTH THE WARNING. `features` is passed to LeverageScore
# EXPLICITLY. Setting VariableFeatures(obj) <- hvg does NOT register for
# layer = "data": LeverageScore looks up variable features for that layer, finds
# none, falls back, and silently scores ALL ~13,700 genes instead of the 2,000
# intended. No error, no warning, a completely different matrix, and a divergence
# that looks like a port bug. Do not remove the `features =` argument.
#
# The two regimes are selected by moving nsketch, not the data: LeverageScore
# takes a truncated-SVD path when ncells < nsketch * 1.5 and a
# CountSketch/QR/Johnson-Lindenstrauss path otherwise. ifnb has 13,999 cells, so
#   nsketch = 10000 -> exact      (13999 < 15000)
#   nsketch =  5000 -> sketched   (13999 > 7500)
# The exact path is reproducible to ~1e-6; the sketched path draws random
# matrices from R's RNG and is only ever comparable statistically.
#
# Writes, into tutorials/figures_sketch/:
#   * r_leverage.csv             cell, r_leverage_exact, r_leverage_sketched
#   * r_sketch_composition.csv   per cell type: share of full data vs of each sketch
#   * r_sketch_cells.csv         the cells each method drew (leverage / uniform)
#   * r_projection.csv           cell, r_predicted, truth   (ProjectData labels)
#   * r_*.png                    R-side figures for the side-by-side
#
# Data: reads the counts exported by tutorials/export_seuratdata.R. Run
#   Rscript tutorials/export_seuratdata.R ifnb     # one-time
#   python  tutorials/ifnb_sketch_tutorial.py      # writes the shared lists
# then
#   Rscript tutorials/ifnb_sketch_verify.R
# Override the data folder with the IFNB_DATA environment variable.
#
# Needs: Seurat (>= 5), ggplot2, irlba.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2)
})
set.seed(42)

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
HERE <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG  <- file.path(HERE, "figures_sketch")
DATA <- Sys.getenv("IFNB_DATA", path.expand("~/.shanuz_data/ifnb"))
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# Mirror the Python constants (ifnb_sketch_tutorial.py).
NSKETCH_EXACT    <- 10000L
NSKETCH_SKETCHED <- 5000L
SKETCH_CELLS     <- 2000L
N_PCS            <- 30L
SEED             <- 123L
CELLTYPE         <- "seurat_annotations"

CELLS_TXT <- file.path(FIG, "cells.txt")
HVG_TXT   <- file.path(FIG, "hvg_features.txt")
if (!file.exists(file.path(DATA, "matrix.mtx")) &&
    !file.exists(file.path(DATA, "matrix.mtx.gz")))
  stop("ifnb counts not found at ", DATA,
       "\nRun `Rscript tutorials/export_seuratdata.R ifnb` first.")
if (!all(file.exists(CELLS_TXT, HVG_TXT)))
  stop("cell / feature lists not found in ", FIG,
       "\nRun `python tutorials/ifnb_sketch_tutorial.py` first (it writes them).")

# ---- 1. Object on the SHARED cell + feature basis ----------------------------
cells <- readLines(CELLS_TXT)
hvg   <- readLines(HVG_TXT)
# Read10X() rewrites underscores in feature names to dashes; normalise Python's
# symbols the same way before they are used to index the R object.
hvg <- gsub("_", "-", hvg)
cat(sprintf("Shared basis from Python: %d cells, %d variable features\n",
            length(cells), length(hvg)))

counts <- Read10X(DATA)
ifnb <- CreateSeuratObject(counts, project = "ifnb_sketch", min.cells = 3)
missing <- setdiff(cells, colnames(ifnb))
if (length(missing))
  stop(length(missing), " of Python's cells are absent from the R object — the ",
       "two loaders disagree about QC; fix that before reading anything below.")
missing_f <- setdiff(hvg, rownames(ifnb))
if (length(missing_f))
  stop(length(missing_f), " of Python's variable features are absent from the R ",
       "object (e.g. ", paste(head(missing_f, 3), collapse = ", "), ")")

ifnb <- subset(ifnb, cells = cells)
ifnb <- NormalizeData(ifnb, normalization.method = "LogNormalize",
                      scale.factor = 10000, verbose = FALSE)
VariableFeatures(ifnb) <- hvg

# The cell-type annotations the metrics are measured against. The exporter
# writes them alongside the counts; without them nothing below is meaningful.
meta_path <- file.path(DATA, "metadata.csv")
if (!file.exists(meta_path))
  stop("metadata.csv not found in ", DATA,
       "\nRe-run `Rscript tutorials/export_seuratdata.R ifnb`.")
meta <- read.csv(meta_path, row.names = 1, check.names = FALSE)
if (!CELLTYPE %in% colnames(meta))
  stop("metadata.csv has no '", CELLTYPE, "' column")
ifnb[[CELLTYPE]] <- meta[colnames(ifnb), CELLTYPE]
cat(sprintf("Cell types: %d, sizes %d..%d\n", length(unique(ifnb[[CELLTYPE]][, 1])),
            min(table(ifnb[[CELLTYPE]][, 1])), max(table(ifnb[[CELLTYPE]][, 1]))))

# ---- 2. LeverageScore, both regimes -----------------------------------------
# `features = hvg` is load-bearing — see the header.
scores <- list()
for (regime in c("exact", "sketched")) {
  nsketch <- if (regime == "exact") NSKETCH_EXACT else NSKETCH_SKETCHED
  t0 <- Sys.time()
  tmp <- Seurat::LeverageScore(ifnb, features = hvg, nsketch = nsketch,
                               seed = SEED, verbose = FALSE)
  s <- tmp[["leverage.score"]][, 1]
  scores[[regime]] <- s
  cat(sprintf("LeverageScore %-9s nsketch=%5d  sum=%10.1f  CV=%.4f  max/median=%5.2f  (%.1fs)\n",
              regime, nsketch, sum(s), sd(s) / mean(s), max(s) / median(s),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
write.csv(data.frame(cell = colnames(ifnb),
                     r_leverage_exact = as.numeric(scores[["exact"]]),
                     r_leverage_sketched = as.numeric(scores[["sketched"]])),
          file.path(FIG, "r_leverage.csv"), row.names = FALSE)

# Does leverage track rarity? The metric the method is actually sold on.
ct <- as.character(ifnb[[CELLTYPE]][, 1])
sizes <- table(ct)
per_type <- tapply(scores[["exact"]], ct, mean)
cat(sprintf("  Spearman(mean leverage, cell-type size) = %+.3f  (want strongly negative)\n",
            cor(as.numeric(per_type), as.numeric(sizes[names(per_type)]),
                method = "spearman")))

# ---- 3. SketchData — leverage vs the uniform control -------------------------
# Uniform is not decoration: a 2000-cell sketch contains rare cells by accident,
# so "the sketch has rare cells" means nothing without the same-size control.
sketch_cells <- list()
for (method in c("LeverageScore", "Uniform")) {
  sk <- SketchData(object = ifnb, ncells = SKETCH_CELLS, method = method,
                   sketched.assay = "sketch", features = hvg,
                   seed = SEED, verbose = FALSE)
  drawn <- colnames(sk[["sketch"]])
  sketch_cells[[method]] <- drawn
  cat(sprintf("SketchData %-14s drew %d cells\n", method, length(drawn)))
  if (method == "LeverageScore") sketched_obj <- sk
}
write.csv(data.frame(
  method = rep(names(sketch_cells), lengths(sketch_cells)),
  cell = unlist(sketch_cells, use.names = FALSE)),
  file.path(FIG, "r_sketch_cells.csv"), row.names = FALSE)

full_frac <- as.numeric(sizes) / length(ct)
names(full_frac) <- names(sizes)
comp <- data.frame(celltype = names(sizes), n_full = as.integer(sizes),
                   frac_full = full_frac)
for (method in names(sketch_cells)) {
  tab <- table(factor(ct[match(sketch_cells[[method]], colnames(ifnb))],
                      levels = names(sizes)))
  frac <- as.numeric(tab) / sum(tab)
  comp[[paste0("frac_", method)]] <- frac
  comp[[paste0("fold_", method)]] <- frac / full_frac
}
comp <- comp[order(-comp$n_full), ]
write.csv(comp, file.path(FIG, "r_sketch_composition.csv"), row.names = FALSE)
cat("  fold-enrichment of the three rarest types (leverage / uniform):\n")
for (i in seq_len(min(3, nrow(comp)))) {
  row <- comp[nrow(comp) - i + 1, ]
  cat(sprintf("    %-22s n=%5d  %.2fx / %.2fx\n", row$celltype, row$n_full,
              row$fold_LeverageScore, row$fold_Uniform))
}

# ---- 4. ProjectData — the sketch's analysis back onto every cell -------------
# Guarded: ProjectData is the heaviest step and the leverage numbers above are
# the tutorial's point, so a failure here must not cost them.
ok <- try({
  DefaultAssay(sketched_obj) <- "sketch"
  # Python's shared feature list, NOT a fresh FindVariableFeatures on the sketch:
  # recomputing here would give the two tools different PCA bases and make the
  # accuracy below a comparison of feature selection rather than of ProjectData.
  sketch_hvg <- intersect(hvg, rownames(sketched_obj[["sketch"]]))
  VariableFeatures(sketched_obj) <- sketch_hvg
  sketched_obj <- ScaleData(sketched_obj, features = sketch_hvg, verbose = FALSE)
  sketched_obj <- RunPCA(sketched_obj, features = sketch_hvg, npcs = N_PCS,
                         verbose = FALSE)
  sketched_obj <- FindNeighbors(sketched_obj, dims = 1:N_PCS, verbose = FALSE)
  sketched_obj <- FindClusters(sketched_obj, resolution = 0.6, verbose = FALSE)
  sketched_obj <- RunUMAP(sketched_obj, dims = 1:N_PCS, return.model = TRUE,
                          verbose = FALSE)
  sketched_obj <- ProjectData(
    object = sketched_obj, assay = "RNA", full.reduction = "pca.full",
    sketched.assay = "sketch", sketched.reduction = "pca",
    umap.model = "umap", dims = 1:N_PCS,
    refdata = list(projected.celltype = CELLTYPE), verbose = FALSE)
  pred <- sketched_obj[["projected.celltype"]][, 1]
  truth <- as.character(ifnb[[CELLTYPE]][colnames(sketched_obj), 1])
  write.csv(data.frame(cell = colnames(sketched_obj),
                       r_predicted = as.character(pred), truth = truth),
            file.path(FIG, "r_projection.csv"), row.names = FALSE)
  cat(sprintf("ProjectData: label accuracy vs held-out annotations = %.4f\n",
              mean(as.character(pred) == truth, na.rm = TRUE)))
}, silent = FALSE)
if (inherits(ok, "try-error"))
  cat("  ProjectData failed; r_projection.csv not written (leverage results stand).\n")

# ---- 5. Figures (r_ prefix, side by side with the Python ones) ---------------
try({
  lev_df <- data.frame(leverage = as.numeric(scores[["exact"]]), celltype = ct)
  lev_df$celltype <- factor(lev_df$celltype,
                            levels = names(sort(sizes, decreasing = TRUE)))
  ggsave(file.path(FIG, "r_01_leverage_by_type.png"),
         ggplot(lev_df, aes(celltype, leverage)) +
           geom_boxplot(outlier.size = 0.3) +
           scale_y_log10() +
           labs(title = "R Seurat — leverage by cell type (commonest first)",
                x = NULL, y = "leverage (log scale)") +
           theme_bw() +
           theme(axis.text.x = element_text(angle = 45, hjust = 1)),
         width = 7.5, height = 4.5, dpi = 150)

  fold_df <- data.frame(
    celltype = factor(comp$celltype, levels = comp$celltype),
    fold = c(comp$fold_LeverageScore, comp$fold_Uniform),
    method = rep(c("LeverageScore", "Uniform"), each = nrow(comp)))
  ggsave(file.path(FIG, "r_02_sketch_enrichment.png"),
         ggplot(fold_df, aes(celltype, fold, fill = method)) +
           geom_col(position = "dodge") +
           geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
           labs(title = "R Seurat — cell-type share in the sketch vs the full data",
                x = NULL, y = "fold enrichment") +
           theme_bw() +
           theme(axis.text.x = element_text(angle = 45, hjust = 1)),
         width = 7.5, height = 4.5, dpi = 150)
}, silent = FALSE)

cat("\nDONE — wrote references + figures to", FIG, "\n")
