"""Visium — the spatial container Seurat's own accessor cannot read.

shanuz's ``load_visium`` / ``VisiumV2`` against Seurat 5.5.1's ``Read10X_Image`` /
``Load10X_Spatial``, on the 10x V1_Mouse_Brain_Sagittal_Anterior Space Ranger
1.1.0 bundle — 2,695 in-tissue spots x 32,285 genes, the anterior half of
``stxBrain``.

Why this one is different from the other sixteen
------------------------------------------------
Every earlier tutorial treated Seurat as the reference and a difference as
shanuz's defect. That was the right default and it found 29 of them. Here it
would have introduced one.

Seurat builds the Visium FOV with ``radius = scale.factors[["spot"]]``, and
``Read10X_ScaleFactors`` sets that field from ``spot_diameter_fullres``. A
diameter goes into a slot named radius. Matching it would have been wrong, so
the tutorial reports the number instead of asserting it, and proves which
reading is correct from the *slide*, not from either tool: Visium spots sit on
a fixed 100 um centre-to-centre grid, so the pixel spacing gives px/um, and
``spot_diameter_fullres`` read as a radius describes overlapping spots.

The second finding needs no arithmetic: ``Radius()`` on a ``VisiumV2`` returns
``NULL``. There are ``Radius`` methods for Centroids, STARmap, SlideSeq,
SpatialImage and VisiumV1 — but none for VisiumV2, Seurat 5's own current
Visium class. The value is stored, and reachable on the centroids underneath.

Usage
-----
    python tutorials/visium_tutorial.py            # writes py_anchors.json
    Rscript tutorials/visium_verify.R              # writes r_anchors.json
    python tutorials/visium_tutorial.py --report   # compares them
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from shanuz import (  # noqa: E402
    find_variable_features,
    load_visium,
    normalize_data,
    run_pca,
    scale_data,
)
from shanuz.datasets import visium_mouse_brain  # noqa: E402

FIGURES = Path(__file__).parent / "figures_visium"
ASSAY = "Spatial"
HEAD = 10
NFEATURES = 2000

# Visium slide geometry, fixed by the hardware and independent of both tools.
SPOT_PITCH_UM = 100.0        # centre-to-centre distance between spots
REFERENCE_SPOT_UM = 65.0     # what 10x defines spot_diameter_fullres against

# Tolerances are each a measured residual rounded up one decade. The floor for
# every one of them is jsonlite's serialisation, which is why visium_verify.R
# writes with digits = 22 rather than digits = NA -- NA emits 15 significant
# digits and does not round-trip a double.
FLOAT_TOLERANCES = {
    # Coordinates are integers read from the same CSV: exact or wrong.
    "coords.x_head": 0.0,
    "coords.y_head": 0.0,
    "coords.x_sum": 0.0,
    "coords.y_sum": 0.0,
    "coords.n": 0.0,
    "load.n_spots": 0.0,
    "load.n_spots_unfiltered": 0.0,
    "obj.n_cells": 0.0,
    "obj.n_features": 0.0,
    # Scale factors are parsed straight from the json.
    "sf.spot": 0.0,
    "sf.fiducial": 0.0,
    "sf.hires": 0.0,
    "sf.lowres": 0.0,
    # Counts.
    "qc.ncount_head": 0.0,
    "qc.nfeature_head": 0.0,
    # The image: R keeps png::readPNG's double, shanuz normalises to float32,
    # so the gap is float32 epsilon on values in [0, 1].
    "image.dim": 0.0,
    "image.corner": 1e-6,
    "image.range": 1e-6,
    "image.mean": 1e-6,
    # Summation order over 2,695 spots / 32,285 genes.
    "norm.head": 1e-13,
    "norm.sum": 1e-11,
    "vst.mean_head": 1e-13,
    "vst.variance_head": 1e-13,
}

# Reported rather than matched — each is a finding or a deliberate difference.
REPORTED_ONLY = (
    "radius.centroids",             # Seurat stores the diameter here
    "radius.visium_is_null",        # no Radius.VisiumV2 method exists
    "radius.has_visiumv2_method",
    "load.default_image_name",      # a filename in R, a resolution word in shanuz
    "load.default_filter_matrix",
    "load.image_class",
    "obj.image_names",
    "obj.assay",
    "coords.cells_head",
    # The known LOESS residual, unchanged from the out-of-core tutorial:
    # variance.standardized is proportional to 1/10^loess_fitted, and shanuz's
    # NumPy local-quadratic fit differs from R's cloess Fortran.
    "vst.var_std_head",
    # Downstream of that residual, not independent of it: the two tools select
    # 1995 of the same 2000 features, so the PCA runs on slightly different
    # matrices and cannot agree exactly. Asserting a loose tolerance here would
    # hide whether the PCA itself is right, so `--report` re-runs it on Seurat's
    # own feature list instead and prints that residual.
    "pca.stdev_head",
)


# ----------------------------------------------------------------------
# pipeline
# ----------------------------------------------------------------------


def build(data_dir=None):
    """Load the bundle on load_visium's own defaults — those are under test."""
    path = data_dir or visium_mouse_brain()
    return load_visium(path)


def run_pipeline(obj):
    normalize_data(obj, assay=ASSAY)
    find_variable_features(obj, assay=ASSAY, selection_method="vst", nfeatures=NFEATURES)
    scale_data(obj, assay=ASSAY)
    run_pca(obj, assay=ASSAY, n_pcs=30)
    return obj


def _positions_row_count(path) -> int:
    """Rows in spatial/tissue_positions*.csv — every spot on the capture area."""
    import pandas as pd

    spatial = Path(path) / "spatial"
    for name in ("tissue_positions.csv", "tissue_positions_list.csv"):
        f = spatial / name
        if f.exists():
            header = 0 if name == "tissue_positions.csv" else None
            return len(pd.read_csv(f, header=header))
    raise FileNotFoundError(f"no tissue_positions file in {spatial}")


def pca_on_shared_features(obj, r_features):
    """Re-run scale+PCA on Seurat's variable features, to separate the PCA from
    the feature selection feeding it."""
    present = [f for f in r_features if f in set(obj.assays[ASSAY].features())]
    scale_data(obj, assay=ASSAY, features=present)
    run_pca(obj, assay=ASSAY, n_pcs=30, features=present,
            reduction_name="pca_shared")
    return np.asarray(obj.reductions["pca_shared"].stdev[:HEAD], dtype=float), len(present)


def spot_geometry(fov):
    """Decide radius-vs-diameter from the slide, not from either tool.

    Visium spots are on a fixed 100 um centre-to-centre grid, so the median
    nearest-neighbour distance in pixels fixes the px/um scale. Reading
    ``spot_diameter_fullres`` as a radius then implies spots wider than the
    pitch, which would mean overlapping capture areas.
    """
    from scipy.spatial import cKDTree

    coords = fov.get_tissue_coordinates()
    xy = coords[["x", "y"]].to_numpy(float)
    nn = float(np.median(cKDTree(xy).query(xy, k=2)[0][:, 1]))
    px_per_um = nn / SPOT_PITCH_UM
    spot_px = float(fov.scale_factors.spot)
    return {
        "geometry.nn_spacing_px": nn,
        "geometry.px_per_um": px_per_um,
        "geometry.spot_um_if_diameter": spot_px / px_per_um,
        "geometry.spot_um_if_radius": 2.0 * spot_px / px_per_um,
        "geometry.pitch_um": SPOT_PITCH_UM,
    }


def collect_anchors(obj, fov):
    a = {}
    import inspect

    sig = inspect.signature(load_visium).parameters
    a["load.default_image_name"] = sig["image_resolution"].default
    a["load.default_filter_matrix"] = bool(sig["filter_by_tissue"].default)
    a["load.n_spots"] = len(fov.cells())
    # Spots on the slide, before either the in-tissue filter or the intersection
    # with the matrix. Seurat's Read10X_Image(filter.matrix = FALSE) reads the
    # positions file alone, so this has to be the positions file too --
    # load_visium(filter_by_tissue=False) still meets filtered_feature_bc_matrix,
    # which is already tissue-filtered, and would answer 2695 for another reason.
    a["load.n_spots_unfiltered"] = _positions_row_count(visium_mouse_brain())
    a["load.image_class"] = type(fov).__name__

    a["radius.centroids"] = float(fov.radius())
    a["radius.visium_is_null"] = fov.radius() is None
    a["radius.has_visiumv2_method"] = hasattr(type(fov), "radius")

    sf = fov.scale_factors
    a["sf.spot"] = float(sf.spot)
    a["sf.fiducial"] = float(sf.fiducial)
    a["sf.hires"] = float(sf.hires)
    a["sf.lowres"] = float(sf.lowres)

    co = fov.get_tissue_coordinates().sort_index()
    a["coords.colnames"] = list(co.columns)
    a["coords.n"] = len(co)
    a["coords.x_head"] = co["x"].to_numpy(float)[:HEAD].tolist()
    a["coords.y_head"] = co["y"].to_numpy(float)[:HEAD].tolist()
    a["coords.x_sum"] = float(co["x"].sum())
    a["coords.y_sum"] = float(co["y"].sum())
    a["coords.cells_head"] = [str(c) for c in co.index[:5]]

    img = fov.get_image()
    a["image.dim"] = list(img.shape)
    a["image.corner"] = img[0, 0].astype(float).tolist()
    a["image.range"] = [float(img.min()), float(img.max())]
    a["image.mean"] = float(img.mean())

    a["obj.n_cells"] = len(obj.cell_names())
    a["obj.n_features"] = len(obj.feature_names())
    a["obj.image_names"] = list(obj.images)
    a["obj.assay"] = obj.active_assay

    assay = obj.assays[ASSAY]
    data = assay.layers["data"]
    dense_col = np.asarray(data[:, 0].todense()).ravel() if hasattr(data, "todense") \
        else np.asarray(data)[:, 0]
    a["norm.head"] = dense_col[:HEAD].tolist()
    a["norm.sum"] = float(data.sum())
    a["qc.ncount_head"] = obj.meta_data[f"nCount_{ASSAY}"].to_numpy(float)[:HEAD].tolist()
    a["qc.nfeature_head"] = obj.meta_data[f"nFeature_{ASSAY}"].to_numpy(float)[:HEAD].tolist()

    md = assay.meta_data
    a["vst.mean_head"] = md["means"].to_numpy(float)[:HEAD].tolist()
    a["vst.variance_head"] = md["variances"].to_numpy(float)[:HEAD].tolist()
    a["vst.var_std_head"] = md["variances.standardized"].to_numpy(float)[:HEAD].tolist()

    a["pca.stdev_head"] = list(map(float, obj.reductions["pca"].stdev[:HEAD]))
    a.update(spot_geometry(fov))
    return a


# ----------------------------------------------------------------------
# comparison
# ----------------------------------------------------------------------


def _match(py, r, tol):
    if isinstance(py, (str, bool)) or isinstance(r, (str, bool)):
        return str(py) == str(r)
    py_arr = np.atleast_1d(np.asarray(py, dtype=object))
    r_arr = np.atleast_1d(np.asarray(r, dtype=object))
    if py_arr.shape != r_arr.shape:
        return False
    if any(isinstance(v, str) for v in py_arr.ravel()):
        return [str(v) for v in py_arr.ravel()] == [str(v) for v in r_arr.ravel()]
    a = np.asarray(py, dtype=float).ravel()
    b = np.asarray(r, dtype=float).ravel()
    if tol == 0.0:
        return bool(np.array_equal(a, b))
    scale = np.maximum(np.abs(a), np.abs(b))
    scale[scale == 0] = 1.0
    return bool((np.abs(a - b) / scale <= tol).all())


def compare(py_anchors, r_anchors):
    matched, differed, reported = [], [], []
    for name in sorted(set(py_anchors) & set(r_anchors)):
        if name in REPORTED_ONLY:
            reported.append((name, py_anchors[name], r_anchors[name]))
            continue
        tol = FLOAT_TOLERANCES.get(name, 0.0)
        if _match(py_anchors[name], r_anchors[name], tol):
            matched.append((name, tol))
        else:
            differed.append((name, py_anchors[name], r_anchors[name], tol))
    return matched, differed, reported


def report_concordance(py_anchors=None, r_anchors=None):
    py_anchors = py_anchors or json.loads((FIGURES / "py_anchors.json").read_text())
    r_path = FIGURES / "r_anchors.json"
    if not r_path.exists():
        raise SystemExit("figures_visium/r_anchors.json missing — run "
                         "`Rscript tutorials/visium_verify.R` first.")
    r_anchors = r_anchors or json.loads(r_path.read_text())
    # jsonlite's I() boxes every scalar in a 1-element array; unwrap so a string
    # or bool prints as itself rather than ['Spatial'].
    r_anchors = {k: (v[0] if isinstance(v, list) and len(v) == 1 else v)
                 for k, v in r_anchors.items()}

    matched, differed, reported = compare(py_anchors, r_anchors)
    total = len(matched) + len(differed)
    print(f"\n{'=' * 74}\nshanuz vs Seurat 5.5.1 — Visium\n{'=' * 74}")
    print(f"  {len(matched)} of {total} compared anchors match "
          f"({len([m for m, t in matched if t == 0.0])} exactly)\n")
    for name, py, r, tol in differed:
        print(f"  DIFFERS  {name}  (tol {tol:g})")
        print(f"           shanuz {np.asarray(py).ravel()[:4]}")
        print(f"           Seurat {np.asarray(r).ravel()[:4]}")

    overlap = None
    r_hvg_path, py_hvg_path = FIGURES / "r_variable_features.txt", \
        FIGURES / "py_variable_features.txt"
    if r_hvg_path.exists() and py_hvg_path.exists():
        r_hvg = {ln.strip() for ln in r_hvg_path.read_text().splitlines() if ln.strip()}
        py_hvg = {ln.strip().replace("_", "-")
                  for ln in py_hvg_path.read_text().splitlines() if ln.strip()}
        overlap = len(py_hvg & r_hvg)
        print(f"\n  variable features: {overlap}/{len(r_hvg)} shared with Seurat")

    print(f"\n{'-' * 74}\nReported, not matched\n{'-' * 74}")
    for name, py, r in reported:
        print(f"  {name:<32} shanuz {str(py)[:24]:<26} Seurat {str(r)[:24]}")

    if overlap is not None and overlap < len(r_hvg):
        print(f"\n{'-' * 74}\nIs the PCA gap the PCA, or the {len(r_hvg) - overlap} "
              f"features upstream of it?\n{'-' * 74}")
        obj = run_pipeline(build())
        r_feats = [ln.strip() for ln in r_hvg_path.read_text().splitlines() if ln.strip()]
        shared, n = pca_on_shared_features(obj, r_feats)
        rstd = np.asarray(r_anchors["pca.stdev_head"], dtype=float)
        own = np.asarray(py_anchors["pca.stdev_head"], dtype=float)
        rel = lambda v: float(np.abs(v - rstd).max() / np.abs(rstd).max())  # noqa: E731
        print(f"  PCA on shanuz's own {len(py_hvg)} features : "
              f"max rel diff vs Seurat {rel(own):.3g}")
        print(f"  PCA on Seurat's {n} features        : "
              f"max rel diff vs Seurat {rel(shared):.3g}")
        print("  -> the difference is the feature selection, not the decomposition.")

    g = {k: v for k, v in py_anchors.items() if k.startswith("geometry.")}
    if g:
        print(f"\n{'-' * 74}\nWhich tool is right about the radius — from the slide"
              f"\n{'-' * 74}")
        print(f"  spot spacing {g['geometry.nn_spacing_px']:.3f} px over a "
              f"{g['geometry.pitch_um']:.0f} um pitch "
              f"=> {g['geometry.px_per_um']:.4f} px/um")
        print(f"  spot_diameter_fullres as a DIAMETER -> "
              f"{g['geometry.spot_um_if_diameter']:.2f} um spot  "
              f"({g['geometry.spot_um_if_diameter'] / REFERENCE_SPOT_UM:.4f}x "
              f"10x's {REFERENCE_SPOT_UM:.0f} um reference)")
        print(f"  spot_diameter_fullres as a RADIUS   -> "
              f"{g['geometry.spot_um_if_radius']:.2f} um spot  "
              f"-> overlaps its neighbour by "
              f"{g['geometry.spot_um_if_radius'] - g['geometry.pitch_um']:.0f} um")
    return matched, differed, reported, overlap


# ----------------------------------------------------------------------
# entry points
# ----------------------------------------------------------------------


def run_full(data_dir=None, verbose=True):
    FIGURES.mkdir(exist_ok=True)
    if verbose:
        print("Loading the Visium bundle...")
    obj = build(data_dir)
    fov = obj.images["slice1"]
    if verbose:
        print(f"  {len(obj.feature_names())} genes x {len(obj.cell_names())} spots, "
              f"{fov.image_resolution} image {fov.get_image().shape}")
        print("Running normalize -> HVG -> scale -> PCA...")
    run_pipeline(obj)

    anchors = collect_anchors(obj, fov)
    (FIGURES / "py_variable_features.txt").write_text(
        "\n".join(obj.assays[ASSAY].variable_features) + "\n")
    (FIGURES / "py_anchors.json").write_text(json.dumps(anchors, indent=2))
    if verbose:
        print(f"\nWrote {FIGURES / 'py_anchors.json'}")
        print("Now run: Rscript tutorials/visium_verify.R")
    return obj, anchors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visium side-by-side")
    parser.add_argument("--report", action="store_true",
                        help="compare against the R anchors")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    if args.report:
        report_concordance()
    else:
        run_full(data_dir=args.data_dir)
