#!/usr/bin/env Rscript
# R Seurat reference for the differential-expression test suite.
#
# Writes tutorials/figures_de/r_<test>.csv, which
# `python tutorials/pbmc3k_de_tutorial.py --report` compares against.
#
# Requires: Seurat, MAST, DESeq2
#   BiocManager::install(c("MAST", "DESeq2"))
#
# NOT glmGamPoi. It is a Suggests of DESeq2 rather than an Imports, FindMarkers
# never calls it, and installing it flips sctransform's vst onto a different
# backend — which would move the SCTransform R reference that
# pbmc3k_sctransform_verify.R is pinned against. Installing with default
# dependencies leaves Suggests alone, which is what you want here.
#
# Groups come from figures_de/groups.csv, written by the Python side, so both
# tools test the same cells. Louvain cluster numbering is not guaranteed to agree
# across implementations and a clustering difference would look exactly like a DE
# difference.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

for (pkg in c("MAST", "DESeq2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("%s is required: BiocManager::install('%s')", pkg, pkg),
         call. = FALSE)
  }
}

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
TUT <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
FIG <- file.path(TUT, "figures_de")
dir.create(FIG, showWarnings = FALSE)

DATA <- path.expand("~/.truecell_data/pbmc3k")
TESTS <- c("wilcox", "t", "bimod", "LR", "negbinom", "roc", "MAST", "DESeq2")

cat("Building pbmc3k...\n")
raw <- Read10X(file.path(DATA, "filtered_gene_bc_matrices/hg19"))
obj <- CreateSeuratObject(raw, project = "pbmc3k_de", min.cells = 3, min.features = 200)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
obj <- NormalizeData(obj, verbose = FALSE)

groups_file <- file.path(FIG, "groups.csv")
if (!file.exists(groups_file)) {
  stop("figures_de/groups.csv is missing; run the Python tutorial first.",
       call. = FALSE)
}
g <- read.csv(groups_file, row.names = 1)
obj <- subset(obj, cells = rownames(g))
Idents(obj) <- factor(as.character(g[colnames(obj), "group"]))
cat(sprintf("Shared cells: %d  (%s)\n", ncol(obj),
            paste(names(table(Idents(obj))), table(Idents(obj)),
                  sep = "=", collapse = " ")))

# A second, exact copy of each table, for the comparisons that need bit-identity.
#
# `write.csv` renders doubles at 15 significant digits, which does not round-trip
# float64 — so a Python side reading it is comparing against a value R has
# already perturbed by a few ULP, and "identical" can never be measured. Raising
# the digit count does not fix it: R's own `sprintf("%.17g", x)` is not correctly
# rounded, and emits digits denoting a *different* double (0.1234567890123456789
# comes back as ...571 where the shortest exact form is ...568). Verified in both
# languages.
#
# `%a` is C99 hex float — a transcription of the IEEE-754 bits rather than a
# decimal approximation of them — so it is exact by construction, and Python
# reads it with `float.fromhex`. The human-readable CSV above is kept as it is,
# because the vignette quotes from it.
write_exact <- function(df, path) {
  out <- as.data.frame(lapply(df, function(col) {
    if (is.numeric(col)) sprintf("%a", col) else as.character(col)
  }), stringsAsFactors = FALSE)
  rownames(out) <- rownames(df)
  write.csv(out, path)
}

# logfc.threshold = 0 and min.pct = 0 so the comparison sees every gene rather
# than only the ones that survive a filter whose *input* is the number under
# test. Seurat's defaults would hide exactly the disagreement worth seeing.
for (tt in TESTS) {
  t0 <- Sys.time()
  res <- try(FindMarkers(obj, ident.1 = "0", ident.2 = "1", test.use = tt,
                         logfc.threshold = 0, min.pct = 0, verbose = FALSE),
             silent = TRUE)
  if (inherits(res, "try-error")) {
    cat(sprintf("  %-9s FAILED: %s\n", tt,
                substr(gsub("\n", " ", as.character(res)), 1, 100)))
    next
  }
  write.csv(res, file.path(FIG, sprintf("r_%s.csv", tolower(tt))))
  write_exact(res, file.path(FIG, sprintf("r_%s_exact.csv", tolower(tt))))
  cat(sprintf("  %-9s %6d genes  %6.1fs\n", tt, nrow(res),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

cat("\nWrote reference tables to", FIG, "\n")
cat("Compare with: python tutorials/pbmc3k_de_tutorial.py --report\n")
