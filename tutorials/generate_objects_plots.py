"""Generate all figures for the object-internals tutorial.

Renders to tutorials/figures_objects/:
  * py_01_split_join.png    the split/join round trip, broken against fixed
  * py_02_concordance.png   anchors matching R, by section
  * py_03_nn_degree.png     kNN degree distribution against Seurat's flat k

01 is the point of the tutorial. The defect it draws is invisible in every
summary statistic — same shape, same values, same column sums — and shows up
only when you look at *which* column each cell's counts landed in. A picture is
the honest way to report it.

03 illustrates the one finding this tutorial deliberately does not fix, so that
"left for a separate comparison" does not quietly become "forgotten".

Usage
-----
    python tutorials/generate_objects_plots.py [--data-dir PATH]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from tutorials.pbmc3k_objects_tutorial import (  # noqa: E402
    FIGURES,
    batch_labels,
    compare_anchors,
    run_full,
)
from shanuz.plotting import _palette  # noqa: E402

FIGURES.mkdir(exist_ok=True)

_PY_COLOR, _R_COLOR = _palette(2)
_BROKEN_COLOR = "#c1666b"
_FIXED_COLOR = "#48a9a6"


def _save(fig, name):
    import matplotlib.pyplot as plt

    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path.name}")


# ---------------------------------------------------------------------------


def split_join_diagram(n_cells=12):
    """The round trip on a toy matrix, before the fix and after.

    Deliberately drawn on twelve cells rather than the tutorial's 2,700: the
    whole point is that you can read the column labels. Each cell's column is
    coloured by its *original* position, so a permutation is visible as a
    scrambled gradient rather than as a number you have to trust.
    """
    import matplotlib.pyplot as plt

    # The same alternating assignment the tutorial splits on.
    batches = batch_labels(n_cells)
    order = np.arange(n_cells)
    even = order[::2]
    odd = order[1::2]
    # What the old code produced: the split layers concatenated as they came.
    broken = np.concatenate([even, odd])
    fixed = order

    fig, axes = plt.subplots(4, 1, figsize=(9, 5.2),
                             gridspec_kw={"hspace": 0.8})
    # (title, what is in each column, what the *labels* on that row say, colour)
    #
    # The split row labels itself honestly, because a split layer records its
    # own cells. The two join rows are labelled with the assay's cell vector —
    # c0…c11, untouched by any of this — because that is what a consumer
    # indexes by. On the broken row the colours therefore run out of step with
    # the labels above them, which is exactly the bug.
    rows = [
        ("original `counts`", order, order, "0.55"),
        (f"after `split(f = batch)` → {batches[0]} | {batches[1]}",
         broken, broken, "0.75"),
        ("after `join_layers()` — before the fix  (labels = `assay.cells()`)",
         broken, order, _BROKEN_COLOR),
        ("after `join_layers()` — after the fix",
         fixed, order, _FIXED_COLOR),
    ]
    for ax, (title, contents, labels, edge) in zip(axes, rows):
        ax.imshow(contents[np.newaxis, :], aspect="auto",
                  cmap="viridis", vmin=0, vmax=n_cells - 1)
        ax.set_yticks([])
        ax.set_xticks(range(n_cells))
        ax.set_xticklabels([f"c{i}" for i in labels], fontsize=7)
        ax.set_title(title, fontsize=9, loc="left", pad=4)
        for spine in ax.spines.values():
            spine.set_edgecolor(edge)
            spine.set_linewidth(2.0)
    # The split row gets a divider so the two layers read as two layers.
    axes[1].axvline(len(even) - 0.5, color="white", lw=2.5)

    fig.suptitle(
        "Colour = the cell whose counts are in that column.\n"
        "Before the fix the join returned every value, in the wrong column, "
        "while the assay's own cell vector never moved.",
        fontsize=9, y=1.06,
    )
    return fig


def concordance_by_section(anchors):
    """How many anchors match R, grouped by the section they belong to."""
    import matplotlib.pyplot as plt

    path = FIGURES / "r_anchors.json"
    if not path.exists():
        return None
    table = compare_anchors(anchors, json.loads(path.read_text()))
    table["section"] = table["field"].str.split(".").str[0]
    grouped = table.groupby("section")["match"].agg(["sum", "count"])
    grouped = grouped.sort_values("count", ascending=True)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    y = np.arange(len(grouped))
    # Teal for agreement, grey for the shortfall: the default palette's first
    # colour is a red that reads as a warning, which is the wrong signal for
    # the bar that means "matches".
    ax.barh(y, grouped["count"], color="0.88", label="anchors compared")
    ax.barh(y, grouped["sum"], color=_FIXED_COLOR, label="matching R exactly")
    ax.set_yticks(y)
    ax.set_yticklabels(grouped.index, fontsize=9)
    ax.set_xlabel("anchors")
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    total, matched = int(grouped["count"].sum()), int(grouped["sum"].sum())
    ax.set_title(f"Agreement with Seurat 5.5.1 — {matched}/{total} anchors",
                 fontsize=11)
    for i, (m, c) in enumerate(zip(grouped["sum"], grouped["count"])):
        if m < c:
            ax.text(c + 0.15, i, f"{c - m} differ", va="center", fontsize=8,
                    color=_BROKEN_COLOR)
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    return fig


def nn_degree(obj):
    """Neighbours per cell, against the flat k Seurat's directed graph gives.

    The finding this tutorial reports but does not fix.
    """
    import matplotlib.pyplot as plt

    graph = obj.graphs.get("RNA_nn")
    if graph is None:
        return None
    degrees = np.diff(graph.tocsr().indptr)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(degrees, bins=range(int(degrees.min()), int(degrees.max()) + 2),
            color=_PY_COLOR, label="shanuz (symmetrized)")
    ax.axvline(20, color=_R_COLOR, lw=2.5,
               label="Seurat (directed): exactly 20 for every cell")
    ax.set_xlabel("neighbours per cell in the `RNA_nn` graph")
    ax.set_ylabel("cells")
    ax.set_title(
        f"kNN degree — shanuz min {degrees.min()}, max {degrees.max()}, "
        f"mean {degrees.mean():.1f}", fontsize=11)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    return fig


def main(data_dir=None):
    # The diagram is schematic and needs no data, so draw it first — it is the
    # figure worth having even if the download is unavailable.
    _save(split_join_diagram(), "py_01_split_join.png")

    obj, anchors = run_full(data_dir=data_dir, verbose=False)

    fig = concordance_by_section(anchors)
    if fig is not None:
        _save(fig, "py_02_concordance.png")
    else:
        print("\n  NOTE: r_anchors.json absent — skipping the concordance "
              "panel.\n  Run `Rscript tutorials/pbmc3k_objects_verify.R` first.")

    fig = nn_degree(obj)
    if fig is not None:
        _save(fig, "py_03_nn_degree.png")

    print(f"\n  Wrote figures to {FIGURES}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="object-internals tutorial figures")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    main(data_dir=args.data_dir)
