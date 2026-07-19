"""Generate the figures for the spatial statistics / container tutorial.

Renders to tutorials/figures_svf/:
  * py_01_moransi_vs_r.png    both weightings against Seurat's Moran's I
  * py_02_top_n_overlap.png   where the kNN approximation actually costs you
  * py_03_spot_radius.png     the radius defect, drawn

01 is the headline: the fixed weighting lands on the identity line, the old one
does not. 02 exists because 01 is not damning enough on its own — a Pearson of
0.986 looks like agreement until you ask which genes come out on top, which is
the only thing this function is used for. 03 is the defect you cannot see in any
number at all, only in a picture.

Usage
-----
    python tutorials/generate_svf_plots.py [--data-dir PATH]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from tutorials.xenium_svf_tutorial import (  # noqa: E402
    FIGURES,
    _subset_object,
    load_slide,
    subset_cells,
)
from shanuz.spatial import find_spatially_variable_features  # noqa: E402

FIGURES.mkdir(exist_ok=True)

_EXACT = "#48a9a6"    # the fixed weighting
_KNN = "#c1666b"      # the old one
_R = "0.35"


def _save(fig, name):
    import matplotlib.pyplot as plt

    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path.name}")


def _r_moransi():
    """R's Moran's I per gene, from the verify script's anchors, if present."""
    path = FIGURES / "r_moransi.csv"
    if path.exists():
        return pd.read_csv(path, index_col=0)["observed"]
    return None


def moransi_vs_r(exact, knn, r):
    """Both weightings against R, on the identity line."""
    import matplotlib.pyplot as plt

    shared = exact.index.intersection(r.index)
    fig, ax = plt.subplots(figsize=(6.2, 6))
    lim = [min(r[shared].min(), knn[shared].min()) - 0.02,
           max(r[shared].max(), knn[shared].max()) + 0.02]
    ax.plot(lim, lim, color=_R, lw=1.0, ls="--", zorder=1,
            label="identity (perfect agreement)")
    ax.scatter(r[shared], knn[shared], s=14, color=_KNN, alpha=0.75, zorder=2,
               label="before — kNN graph (k=10)")
    ax.scatter(r[shared], exact[shared], s=14, color=_EXACT, alpha=0.9, zorder=3,
               label="after — R's 1/d², row-standardised")
    ax.set_xlabel("Seurat 5.5.1  ·  FindSpatiallyVariableFeatures(moransi)")
    ax.set_ylabel("shanuz  ·  find_spatially_variable_features")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_title("Moran's I per gene, 248 genes on 2,000 shared cells",
                 fontsize=11)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    return fig


def top_n_overlap(exact, knn, r):
    """How many of R's top N each weighting recovers, as N grows.

    The point the scatter understates: the approximation is good on average and
    worst exactly where the function is read.
    """
    import matplotlib.pyplot as plt

    ns = np.arange(5, 105, 5)
    r_rank = r.sort_values(ascending=False)
    e_rank = exact.sort_values(ascending=False)
    k_rank = knn.sort_values(ascending=False)

    def frac(series_rank, n):
        return len(set(series_rank.index[:n]) & set(r_rank.index[:n])) / n

    fig, ax = plt.subplots(figsize=(7, 4.2))
    # Plotted as a *fraction*: the raw count rises with N either way, which makes
    # a 3-of-10 miss at the head of the list look like a rounding error next to a
    # 7-of-100 miss in the tail. It is the opposite.
    ax.axhline(1.0, ls="--", lw=1.0, color=_R, zorder=1,
               label="perfect recovery")
    ax.plot(ns, [frac(e_rank, n) for n in ns], "-o", ms=4, color=_EXACT,
            zorder=3, label="after — R's weighting (exact at every N)")
    ax.plot(ns, [frac(k_rank, n) for n in ns], "-o", ms=4, color=_KNN,
            zorder=2, label="before — kNN graph")
    ax.set_xlabel("N — size of the top-N gene list")
    ax.set_ylabel("fraction of Seurat's top N recovered")
    ax.set_ylim(0.6, 1.03)
    ax.set_title(
        "Recovery of Seurat's ranking — the approximation is worst\n"
        "at the head of the list, which is the part anyone reads",
        fontsize=11)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    return fig


def spot_radius(obj):
    """The radius defect: fixed-size dots against true-to-scale spots.

    Drawn as outlines on a tight crop. Filled circles at this radius merge into a
    solid block — which is not a drawing mistake to hide but the honest result:
    ``.AutoRadius`` is 1% of the *bounding box*, so on a densely packed section it
    lands well above the cell spacing. Seurat draws this slide the same way. The
    defect being fixed is that shanuz had no radius at all, not that R's heuristic
    is the right size for every slide.
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import EllipseCollection

    fov = obj.images[list(obj.images)[0]]
    coords = fov.get_tissue_coordinates()
    radius = fov.boundaries[fov.default_boundary()].radius()
    x, y = coords["x"].to_numpy(), coords["y"].to_numpy()
    # ~200 µm across: a few dozen cells, so individual outlines stay legible.
    cx, cy, half = np.median(x), np.median(y), 100.0
    m = (np.abs(x - cx) < half) & (np.abs(y - cy) < half)
    xs, ys = x[m], y[m]

    # Half the median nearest-neighbour spacing: what a reader would actually
    # pass for a section this dense, and only expressible now that the slot is
    # populated at all.
    from shanuz.spatial import spatial_knn

    tuned = 0.5 * float(np.median(spatial_knn(np.column_stack([x, y]), k=1)[0][:, 0]))

    fig, axes = plt.subplots(1, 3, figsize=(14, 5.2), sharex=True, sharey=True)
    axes[0].scatter(xs, ys, s=18, color=_KNN, linewidths=0)
    axes[0].set_title(
        f"before — radius None\n{len(xs)} cells, fixed-size dots at any zoom",
        fontsize=10)
    for ax, r_, title in (
        (axes[1], radius,
         f"R's .AutoRadius — {radius:.2f} µm\ntrue-to-scale, and overlapping as R's does"),
        (axes[2], tuned,
         f"an explicit radius — {tuned:.2f} µm\nhalf the median cell spacing"),
    ):
        ax.add_collection(EllipseCollection(
            widths=2 * r_, heights=2 * r_, angles=0.0, units="xy",
            offsets=np.column_stack([xs, ys]), offset_transform=ax.transData,
            facecolors="none", edgecolors=_EXACT, linewidths=1.1,
        ))
        ax.scatter(xs, ys, s=3, color=_EXACT, linewidths=0)
        ax.set_title(title, fontsize=10)
    for ax in axes:
        ax.set_xlim(cx - half, cx + half)
        ax.set_ylim(cy - half, cy + half)
        ax.set_aspect("equal")
        ax.set_xlabel("x (µm)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("y (µm)")
    fig.suptitle(
        "Centroids carried no radius, so every true-to-scale renderer silently "
        "fell back to dots (left).\n"
        "The fix restores R's default (middle) and, just as importantly, makes "
        "the slot mean something you can set (right).",
        fontsize=10, y=1.02)
    return fig


def main(data_dir=None):
    obj = load_slide(data_dir)
    cells = subset_cells(obj)
    sub = _subset_object(obj, cells)

    exact = find_spatially_variable_features(
        sub, method="moransi", weights="inverse_square")["moransi"]
    knn = find_spatially_variable_features(
        sub, method="moransi", weights="knn")["moransi"]

    r = _r_moransi()
    if r is not None:
        _save(moransi_vs_r(exact, knn, r), "py_01_moransi_vs_r.png")
        _save(top_n_overlap(exact, knn, r), "py_02_top_n_overlap.png")
    else:
        print("\n  NOTE: no R Moran's I available — skipping panels 01 and 02.")
        print("  Run `Rscript tutorials/xenium_svf_verify.R` first.\n")

    _save(spot_radius(obj), "py_03_spot_radius.png")
    print(f"\n  Wrote figures to {FIGURES}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="spatial-statistics figures")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    main(data_dir=args.data_dir)
