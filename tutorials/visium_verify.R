#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# R side of the Visium tutorial: Seurat 5.5.1 on the 10x V1_Mouse_Brain_Sagittal
# _Anterior Space Ranger 1.1.0 bundle, the same files shanuz reads.
#
# Emits tutorials/figures_visium/r_anchors.json.
#
# On precision: `digits = 22`, not `digits = NA`. jsonlite's NA writes 15
# significant digits, which does not round-trip a double — 89.47199235723474
# comes back as 89.4719923572347. That 4e-14 gap is larger than several of the
# residuals under test, and in an earlier tutorial it invented three findings
# that were not real. Never round on the way out.
#
# Usage:  Rscript tutorials/visium_verify.R
# ---------------------------------------------------------------------------
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))
suppressMessages(library(jsonlite))

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
TUT <- if (length(.script)) dirname(normalizePath(.script)) else getwd()
OUT <- file.path(TUT, "figures_visium")
dir.create(OUT, showWarnings = FALSE)
DATA <- path.expand("~/.shanuz_data/visium_mouse_brain")

anchors <- list()
add <- function(name, value) anchors[[name]] <<- I(value)
addn <- function(name, value) anchors[[name]] <<- I(as.numeric(value))

if (!dir.exists(file.path(DATA, "spatial"))) {
  stop("Visium bundle missing. Run in Python first:\n",
       "  from shanuz.datasets import visium_mouse_brain; visium_mouse_brain()")
}

# ---------------------------------------------------------------------------
# The container: Read10X_Image on its own defaults
# ---------------------------------------------------------------------------
img <- Read10X_Image(image.dir = file.path(DATA, "spatial"))

add("load.default_image_name", as.character(formals(Read10X_Image)$image.name))
add("load.default_filter_matrix", as.logical(formals(Read10X_Image)$filter.matrix))
addn("load.n_spots", length(Cells(img)))
addn("load.n_spots_unfiltered",
     length(Cells(Read10X_Image(file.path(DATA, "spatial"), filter.matrix = FALSE))))
add("load.image_class", class(img)[1])

# Radius: reported, not matched. Seurat passes scale.factors[["spot"]] --
# spot_diameter_fullres -- straight into CreateFOV's `radius`, and there is no
# Radius.VisiumV2 method, so the accessor on the class itself returns NULL.
addn("radius.centroids", Radius(img@boundaries[["centroids"]]))
add("radius.visium_is_null", is.null(Radius(img)))
add("radius.has_visiumv2_method",
    "Radius.VisiumV2" %in% as.character(utils::methods("Radius")))

sf <- ScaleFactors(img)
addn("sf.spot", sf$spot); addn("sf.fiducial", sf$fiducial)
addn("sf.hires", sf$hires); addn("sf.lowres", sf$lowres)

co <- GetTissueCoordinates(img)
co <- co[order(co$cell), ]
add("coords.colnames", colnames(GetTissueCoordinates(img)))
addn("coords.n", nrow(co))
addn("coords.x_head", head(co$x, 10)); addn("coords.y_head", head(co$y, 10))
addn("coords.x_sum", sum(co$x));       addn("coords.y_sum", sum(co$y))
add("coords.cells_head", head(co$cell, 5))

im <- GetImage(img, mode = "raw")
addn("image.dim", dim(im))
addn("image.corner", im[1, 1, ])
addn("image.range", range(im))
addn("image.mean", mean(im))

# ---------------------------------------------------------------------------
# The object, then the standard pipeline
# ---------------------------------------------------------------------------
# Load10X_Spatial reads an .h5; this bundle ships the MTX triplet, so this is
# the documented manual equivalent it wraps.
counts <- Read10X(file.path(DATA, "filtered_feature_bc_matrix"))
obj <- CreateSeuratObject(counts = counts, assay = "Spatial")
obj[["slice1"]] <- img[Cells(obj)]

addn("obj.n_cells", ncol(obj)); addn("obj.n_features", nrow(obj))
add("obj.image_names", Images(obj))
add("obj.assay", DefaultAssay(obj))

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst",
                            nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

norm <- GetAssayData(obj, assay = "Spatial", layer = "data")
addn("norm.head", norm[1:10, 1]); addn("norm.sum", sum(norm))
addn("qc.ncount_head", head(obj$nCount_Spatial, 10))
addn("qc.nfeature_head", head(obj$nFeature_Spatial, 10))

hvi <- HVFInfo(obj)
addn("vst.mean_head", head(hvi$mean, 10))
addn("vst.variance_head", head(hvi$variance, 10))
addn("vst.var_std_head", head(hvi$variance.standardized, 10))
writeLines(VariableFeatures(obj), file.path(OUT, "r_variable_features.txt"))

# PCA sign is arbitrary; Stdev is not.
addn("pca.stdev_head", head(Stdev(obj, reduction = "pca"), 10))

writeLines(toJSON(anchors, auto_unbox = TRUE, digits = 22),
           file.path(OUT, "r_anchors.json"))
cat("Wrote", file.path(OUT, "r_anchors.json"), "\n")
cat("  spots", length(Cells(img)), " radius", Radius(img@boundaries[["centroids"]]),
    " Radius(VisiumV2) NULL:", is.null(Radius(img)), "\n")
