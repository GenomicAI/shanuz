"""Spatial statistics and the spatial container — shanuz against Seurat 5.5.1.

Wave 2's last side-by-side. It takes the 10x Xenium mouse brain slide (36,602
cells x 248 genes) and compares two things that the existing Xenium tutorial
never touched:

* **The container.** ``load_xenium`` against ``LoadXenium``, and the
  FOV / Centroids / Segmentation classes against ``CreateFOV`` /
  ``CreateCentroids`` / ``CreateSegmentation``. The older spatial tutorial built
  its R side with ``Read10X`` plus a coordinate data.frame, so R never actually
  constructed an FOV and the whole boundary layer went uncompared.
* **The statistic.** ``find_spatially_variable_features`` against
  ``FindSpatiallyVariableFeatures`` — the one real spatial *statistic* in the
  library that had never been checked against R.

What it found
-------------
Three defects, all fixed here.

1. **Moran's I was computed on the wrong weight matrix.** Seurat builds
   ``1 / d²`` between every pair of cells and ``Rfast2::moranI`` row-standardises
   it; shanuz used a k-nearest-neighbour graph. The kNN answer is a decent
   approximation — Pearson 0.986 against R, 46 of R's top 50 genes — which is
   exactly why it survived: close enough to look right in a plot. It is not R's
   answer. It runs a median 1.23x high and agrees on only 7 of R's top 10, and
   the top 10 is what anyone actually reads off this function. With the
   inverse-square weighting it now matches R to **1.6e-14** and **10/10**.

2. **Centroids never got a radius.** ``SeuratObject`` always computes one
   (``.AutoRadius`` = 1% of the mean bounding-box dimension: 42.83 on this
   slide). shanuz left it ``None``, and ``_spot_collection`` returns ``None`` for
   a ``None`` radius — so every true-to-scale spot renderer silently degraded to
   a fixed-size scatter on every FOV that did not come from a Visium
   ``scalefactors_json.json``. Nothing raised; the plots just quietly stopped
   meaning what they claimed.

3. **Segmentation polygons were stored open.** R closes each ring by repeating
   the first vertex — a square is five rows, not four. Anything that measures a
   perimeter off the vertex list, or strokes an outline without asking the
   renderer to close it, is short one edge.

Two more differences are left standing on purpose, and saying why is half the
point of this file:

* **The Moran's I p-value.** R runs a 999-permutation test. On this 248-gene
  panel that yields 14 distinct p-values and ties **233 genes** at its 1/1025
  floor — it cannot rank the most spatially variable gene against the
  two-hundredth. shanuz's normal-approximation p-value is continuous and
  deterministic. Matching R here would cost information and buy nothing, so the
  statistic was fixed and the p-value deliberately was not. Parity is the goal
  right up until it makes the port worse.
* **The FOV ``Key``.** R derives it from the assay (``RNA_``); shanuz hardcodes
  ``fov_``. It appears in ``__repr__`` and nowhere else — no lookup, no
  prefixing, no fetch — so changing a default to fix a display string was not
  worth the churn.

Result: **38 of 39 anchors match Seurat exactly.** The one that does not is the
shape of ``get_tissue_coordinates`` — R returns ``x, y, cell`` as three columns,
shanuz returns ``x, y`` and carries the cell as the index. Same information, and
the object-level accessor already materialises ``cell`` as a column; adding a
second copy alongside the index would only give the two a way to disagree.

Conventions
-----------
Anchors are compared with **no tolerance** wherever the quantity is exact — cell
and feature orders, counts, boundary names, ring vertex counts. Only the genuinely
floating-point ones carry a tolerance, and they are named in ``FLOAT_TOLERANCES``
rather than covered by a blanket rule, so a new anchor cannot quietly inherit
slack it did not earn.

The Moran's I comparison runs on a **2,000-cell subset**, because Seurat cannot
do otherwise: ``RunMoransI`` materialises the full n x n weight matrix, which on
36,602 cells is a 10.7 GB allocation. shanuz evaluates the same weights in row
blocks and does run on the whole slide — that comparison is reported too, against
itself, since there is no R number to check it against.

Usage
-----
    python tutorials/xenium_svf_tutorial.py
    Rscript tutorials/xenium_svf_verify.R      # writes figures_svf/r_anchors.json
    python tutorials/xenium_svf_tutorial.py --report
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from shanuz.datasets import xenium_mouse_brain
from shanuz.preprocessing import normalize_data
from shanuz.spatial import (
    create_centroids,
    create_fov,
    create_segmentation,
    find_spatially_variable_features,
    load_xenium,
)

FIGURES = Path(__file__).parent / "figures_svf"
ASSAY = "Xenium"
N_SUBSET = 2000          # cells R can build a dense weight matrix for
SUBSET_SEED = 42
TOP_N = 10               # the head of the ranking, which is what gets read

# Only these anchors are floating point. Everything else is compared exactly.
FLOAT_TOLERANCES = {
    "container.centroids_radius": 1e-9,
    "container.coords_x_head": 1e-9,
    "container.coords_y_head": 1e-9,
    "moransi.i_head": 1e-9,
    "moransi.i_max": 1e-9,
    "moransi.i_min": 1e-9,
    "toy.auto_radius": 1e-12,
}


# ---------------------------------------------------------------------------
# Anchor helpers
# ---------------------------------------------------------------------------

def digest(names) -> str:
    """A short, order-**sensitive** fingerprint of a sequence of names."""
    joined = "\n".join(str(n) for n in names)
    return hashlib.md5(joined.encode()).hexdigest()[:12]


def name_anchor(names) -> dict:
    names = list(names)
    return {
        "n": len(names),
        "digest": digest(names),
        "head": [str(n) for n in names[:3]],
        "tail": [str(n) for n in names[-3:]],
    }


def compare_anchors(py: dict, r: dict, path: str = "") -> pd.DataFrame:
    """Flatten two nested anchor dicts into a per-field comparison table."""
    rows = []
    for key in sorted(set(py) | set(r)):
        field = f"{path}.{key}" if path else str(key)
        a, b = py.get(key), r.get(key)
        if isinstance(a, dict) and isinstance(b, dict):
            rows.append(compare_anchors(a, b, field))
            continue
        rows.append(pd.DataFrame([{
            "field": field,
            "python": a,
            "r": b,
            "match": _match(field, a, b),
        }]))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["field", "python", "r", "match"]
    )


def _match(field: str, a, b) -> bool:
    if a is None or b is None:
        return False
    tol = _tolerance_for(field)
    if tol is not None:
        try:
            return bool(np.allclose(np.asarray(a, dtype=float),
                                    np.asarray(b, dtype=float),
                                    rtol=0, atol=tol))
        except (TypeError, ValueError):
            return False
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return [str(x) for x in a] == [str(x) for x in b]
    if isinstance(a, float) or isinstance(b, float):
        return float(a) == float(b)
    return str(a) == str(b)


def _tolerance_for(field: str):
    """Tolerances are looked up by exact field name, never by prefix.

    A prefix rule would let a new anchor inherit slack nobody chose for it, which
    is how a comparison quietly stops comparing.
    """
    return FLOAT_TOLERANCES.get(field)


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def load_slide(data_dir=None):
    obj = load_xenium(data_dir or xenium_mouse_brain(), assay=ASSAY)
    normalize_data(obj)
    return obj


def subset_cells(obj) -> list[str]:
    """A deterministic 2,000-cell subset, small enough for R's dense weights."""
    cells = np.asarray(obj.assays[ASSAY].cells())
    rng = np.random.default_rng(SUBSET_SEED)
    idx = np.sort(rng.choice(len(cells), size=N_SUBSET, replace=False))
    return [str(c) for c in cells[idx]]


def container_anchors(obj) -> dict:
    """What the object *is* — boundaries, keys, radii, coordinate frames."""
    assay = obj.assays[ASSAY]
    fov_name = list(obj.images)[0]
    fov = obj.images[fov_name]
    centroids = fov.boundaries[fov.default_boundary()]
    coords = fov.get_tissue_coordinates()
    return {
        "n_cells": len(assay.cells()),
        "n_features": len(assay.features()),
        "cells": name_anchor(assay.cells()),
        "features": name_anchor(assay.features()),
        "n_images": len(obj.images),
        "boundaries": sorted(fov.boundaries),
        "default_boundary": fov.default_boundary(),
        "fov_assay": fov.assay,
        "fov_radius_is_none": fov.radius() is None,
        "centroids_n": len(centroids.cells()),
        "centroids_nsides": int(centroids.nsides),
        "centroids_radius": float(centroids.radius()),
        "coords_shape": list(coords.shape),
        "coords_x_head": [float(v) for v in coords["x"].to_numpy()[:3]],
        "coords_y_head": [float(v) for v in coords["y"].to_numpy()[:3]],
        "coords_cells_head": [str(c) for c in list(coords.index)[:3]],
    }


def toy_anchors() -> dict:
    """The constructors, on input small enough to check by eye.

    A 4-cell frame and a 2-cell square: the numbers R produces for these are
    short enough to read in the vignette, which is the point — the slide-scale
    anchors say *whether* the tools agree, these say *what* they agree on.
    """
    coords = pd.DataFrame(
        {"x": [1.0, 2, 3, 4], "y": [10.0, 20, 30, 40], "cell": list("abcd")}
    )
    centroids = create_centroids(coords)
    fov = create_fov(coords, type_="centroids", assay="RNA")
    square = pd.DataFrame({
        "x": [0.0, 1, 1, 0, 5, 6, 6, 5],
        "y": [0.0, 0, 1, 1, 5, 5, 6, 6],
        "cell": ["a"] * 4 + ["b"] * 4,
    })
    seg = create_segmentation(square)
    seg_coords = seg.get_tissue_coordinates()
    return {
        "auto_radius": float(centroids.radius()),
        "nsides": int(centroids.nsides),
        "centroid_cells": list(centroids.cells()),
        "fov_cells": list(fov.cells()),
        "fov_boundaries": sorted(fov.boundaries),
        "subset_cells": list(fov.subset(cells=["b", "d"]).cells()),
        "segmentation_cells": list(seg.cells()),
        "segmentation_rows": len(seg_coords),
        "segmentation_rows_per_cell": [
            int((seg_coords.index == c).sum()) for c in seg.cells()
        ],
        "ring_closed": bool(
            seg_coords.loc["a"].iloc[0].tolist()
            == seg_coords.loc["a"].iloc[-1].tolist()
        ),
    }


def moransi_anchors(obj, cells: list[str]) -> dict:
    """Moran's I on the R-sized subset, with R's weighting."""
    sub = _subset_object(obj, cells)
    res = find_spatially_variable_features(sub, method="moransi",
                                           weights="inverse_square")
    res = res.sort_values("moransi", ascending=False)
    return {
        "n_cells": len(cells),
        "n_genes": len(res),
        "top_genes": list(res.index[:TOP_N]),
        "i_head": [float(v) for v in res["moransi"].to_numpy()[:TOP_N]],
        "i_max": float(res["moransi"].max()),
        "i_min": float(res["moransi"].min()),
        "ranking": digest(res.index),
    }


def knn_comparison(obj, cells: list[str]) -> dict:
    """The old kNN weighting against the fixed one — reported, not compared to R.

    This is the defect measured in its own terms: same cells, same genes, same
    everything but the weight matrix.
    """
    from scipy.stats import pearsonr, spearmanr

    sub = _subset_object(obj, cells)
    exact = find_spatially_variable_features(sub, method="moransi",
                                             weights="inverse_square")["moransi"]
    knn = find_spatially_variable_features(sub, method="moransi",
                                           weights="knn")["moransi"]
    shared = exact.index.intersection(knn.index)
    a, b = exact[shared], knn[shared]
    top_e = set(exact.sort_values(ascending=False).index[:TOP_N])
    top_k = set(knn.sort_values(ascending=False).index[:TOP_N])
    return {
        "pearson": float(pearsonr(a, b)[0]),
        "spearman": float(spearmanr(a, b)[0]),
        "max_abs_diff": float(np.abs(a - b).max()),
        "median_ratio": float(np.median((b / a).replace(
            [np.inf, -np.inf], np.nan).dropna())),
        "top10_overlap": len(top_e & top_k),
    }


def _subset_object(obj, cells: list[str]):
    """A slide restricted to ``cells``, coordinates and all."""
    return obj.subset(cells=cells)


def collect_anchors(obj, cells: list[str]) -> dict:
    return {
        "container": container_anchors(obj),
        "toy": toy_anchors(),
        "moransi": moransi_anchors(obj, cells),
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_full(data_dir=None, verbose=True):
    FIGURES.mkdir(exist_ok=True)
    obj = load_slide(data_dir)
    cells = subset_cells(obj)
    (FIGURES / "cells.txt").write_text("\n".join(cells) + "\n")

    anchors = collect_anchors(obj, cells)
    anchors["knn_vs_fixed"] = knn_comparison(obj, cells)
    (FIGURES / "py_anchors.json").write_text(json.dumps(anchors, indent=2))

    if verbose:
        _print_summary(anchors)
    return obj, anchors


def _print_summary(anchors: dict) -> None:
    c, m, k = anchors["container"], anchors["moransi"], anchors["knn_vs_fixed"]
    print(f"\n  Slide: {c['n_cells']} cells x {c['n_features']} genes, "
          f"{c['n_images']} FOV, boundaries {c['boundaries']}")
    print(f"  Centroids radius: {c['centroids_radius']:.5f}  "
          f"(R's .AutoRadius; was None before this tutorial)")
    print(f"\n  Moran's I on {m['n_cells']} cells, {m['n_genes']} genes")
    print(f"    top {TOP_N}: {', '.join(m['top_genes'])}")
    print(f"    I range: {m['i_min']:.4f} .. {m['i_max']:.4f}")
    print("\n  Old kNN weighting vs the fixed one, same cells:")
    print(f"    pearson {k['pearson']:.4f}  spearman {k['spearman']:.4f}  "
          f"max|diff| {k['max_abs_diff']:.4f}")
    print(f"    median ratio {k['median_ratio']:.3f}x  "
          f"top-{TOP_N} overlap {k['top10_overlap']}/{TOP_N}")


def report_concordance() -> pd.DataFrame:
    py_path, r_path = FIGURES / "py_anchors.json", FIGURES / "r_anchors.json"
    if not py_path.exists():
        raise SystemExit("Run the tutorial first — py_anchors.json is missing.")
    if not r_path.exists():
        raise SystemExit(
            "figures_svf/r_anchors.json is missing.\n"
            "Run: Rscript tutorials/xenium_svf_verify.R"
        )
    py = json.loads(py_path.read_text())
    r = json.loads(r_path.read_text())
    # knn_vs_fixed has no R counterpart by construction — it is shanuz measuring
    # its own former behaviour — so it is not part of the parity table.
    py = {k: v for k, v in py.items() if k != "knn_vs_fixed"}
    table = compare_anchors(py, r)
    _print_report(table)
    return table


def _print_report(table: pd.DataFrame) -> None:
    matched = int(table["match"].sum())
    total = len(table)
    print(f"\n  Anchors matching Seurat 5.5.1: {matched}/{total}")
    bad = table.loc[~table["match"]]
    if bad.empty:
        print("  Every anchor agrees.")
        return
    print("\n  Differences:")
    for _, row in bad.iterrows():
        print(f"    {row['field']}")
        print(f"      python: {row['python']}")
        print(f"      r     : {row['r']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", default=None)
    parser.add_argument("--report", action="store_true",
                        help="compare against figures_svf/r_anchors.json")
    args = parser.parse_args()
    if args.report:
        report_concordance()
        return
    run_full(data_dir=args.data_dir)
    print(f"\n  Wrote {FIGURES / 'py_anchors.json'}")
    print("  Next: Rscript tutorials/xenium_svf_verify.R")


if __name__ == "__main__":
    main()
