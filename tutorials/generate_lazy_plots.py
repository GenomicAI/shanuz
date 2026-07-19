"""Generate the figures for the out-of-core tutorial.

Renders to tutorials/figures_lazy/:
  * py_01_self_consistency.png  each tool's on-disk path against its own in-memory one
  * py_02_memory.png            what the densification cost, and what removing it saved
  * py_03_storage.png           what the two on-disk formats cost

01 is the tutorial's finding: whether a tool returns the same answer once the
matrix moves to disk. 02 is why the fixes mattered — before them, backing a
matrix on disk *raised* peak memory. 03 is where shanuz loses.

The R-side numbers are constants read off a live Seurat 5.5.1 / BPCells 0.3.1
session rather than recomputed here, so these figures render without R
installed. `Rscript tutorials/lazy_bpcells_verify.R` reproduces them.

Usage
-----
    python tutorials/generate_lazy_plots.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from tutorials.lazy_bpcells_tutorial import FIGURES  # noqa: E402

_SHANUZ = "#48a9a6"
_SEURAT = "#c1666b"
_GREY = "0.35"
_BEFORE = "#c1666b"

# R: Seurat 5.5.1 + BPCells 0.3.1, pbmc3k 13,714 x 2,638, in-memory run against
# out-of-core run. See tutorials/lazy_bpcells_verify.R.
SEURAT_SELF = {
    "normalise": 1.0153110547861e-06,
    "scale.data": 1.7167773513904e-06,
    "variance.standardized": 2.05553676092025e-02,
}
# shanuz's are all exactly zero; plotted at a floor so a log axis can show them.
FLOOR = 1e-17

# Measured on this machine, pbmc3k, separate processes. See the vignette.
MEMORY = {
    "in-memory (sparse)": 253,
    "on-disk, before the fixes": 1169,
    "on-disk, after": 335,
}
# bytes
STORAGE = {
    "dgCMatrix arrays": 26_875_340,
    "shanuz LazyMatrix": 26_875_830,
    "BPCells (float64)": 20_710_000,
    "BPCells (uint32)": 4_501_708,
}


def _save(fig, name):
    import matplotlib.pyplot as plt

    FIGURES.mkdir(exist_ok=True)
    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path.name}")


def self_consistency():
    """Does moving the matrix to disk change the answer?"""
    import matplotlib.pyplot as plt

    labels = list(SEURAT_SELF)
    seurat = [SEURAT_SELF[k] for k in labels]
    y = np.arange(len(labels))
    left = 1e-9

    fig, ax = plt.subplots(figsize=(8.6, 3.4))
    ax.barh(y, seurat, height=0.42, color=_SEURAT, left=left,
            label="Seurat 5.5.1 + BPCells")
    # shanuz's differences are exactly zero, which has no length on a log axis.
    # Drawing a stub at some floor would read as a measurement; a marker on the
    # zero rule says what is true.
    ax.plot([left] * len(labels), y, "o", color=_SHANUZ, ms=8, zorder=3,
            label="shanuz — exactly zero")
    ax.axvline(left, color=_SHANUZ, lw=1.2, ls=":", zorder=1)

    ax.set_xscale("log")
    ax.set_xlim(left / 3, 3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("max |on-disk − in-memory|, same input matrix")
    ax.set_title("A tool against itself: what moving the matrix to disk costs\n"
                 "pbmc3k, 13,714 × 2,638 — identical input to both paths",
                 fontsize=11)
    for i, v in enumerate(seurat):
        ax.text(v * 2.2, i, f"{v:.2e}", va="center", fontsize=8.5, color=_SEURAT)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(0.5, -0.14,
             "Bars start at the zero rule; shanuz has no bar because its two "
             "paths return the same bits. BPCells computes in single precision, "
             "so Seurat's\nout-of-core run also selects a different variable "
             "feature — 1999/2000 agree with its own in-memory run.",
             ha="center", fontsize=8.5, color=_GREY)
    return fig


def memory():
    """Peak RSS for `normalize_data`, and the layer it leaves behind."""
    import matplotlib.pyplot as plt

    labels = list(MEMORY)
    values = [MEMORY[k] for k in labels]
    colours = ["0.78", _BEFORE, _SHANUZ]

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.0),
                                  gridspec_kw={"width_ratios": [1.35, 1]})
    bars = ax.bar(range(len(labels)), values, color=colours, width=0.62)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(["in-memory\n(sparse)", "on-disk\nbefore", "on-disk\nafter"],
                       fontsize=9)
    ax.set_ylabel("peak RSS, MB")
    ax.set_title("`normalize_data` on pbmc3k", fontsize=11)
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, v + 25, f"{v} MB",
                ha="center", fontsize=9)
    ax.set_ylim(0, max(values) * 1.18)
    ax.spines[["top", "right"]].set_visible(False)

    layers = [18.3, 296.2, 18.3]
    bars2 = ax2.bar(range(3), layers, color=colours, width=0.62)
    ax2.set_xticks(range(3))
    ax2.set_xticklabels(["sparse", "before", "after"], fontsize=9)
    ax2.set_ylabel("resulting `data` layer, MB")
    ax2.set_title("...and what it left behind", fontsize=11)
    for bar, v, t in zip(bars2, layers, ("csc_matrix", "ndarray", "csc_matrix")):
        ax2.text(bar.get_x() + bar.get_width() / 2, v + 8,
                 f"{v:.1f}\n{t}", ha="center", fontsize=8)
    ax2.set_ylim(0, max(layers) * 1.25)
    ax2.spines[["top", "right"]].set_visible(False)

    fig.suptitle("Before the fixes, moving the matrix on disk cost 4.6× the peak "
                 "memory and produced a dense layer 16× larger.",
                 fontsize=10, y=1.03)
    return fig


def storage():
    """What each on-disk format costs — the axis shanuz loses on."""
    import matplotlib.pyplot as plt

    labels = list(STORAGE)
    values = [STORAGE[k] / 1e6 for k in labels]
    colours = ["0.78", _SHANUZ, _SEURAT, _SEURAT]

    fig, ax = plt.subplots(figsize=(7.6, 3.8))
    bars = ax.barh(range(len(labels)), values, color=colours, height=0.6)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("on disk, MB  ·  pbmc3k counts, 2,238,732 non-zeros")
    ax.set_title("shanuz stores the CSC arrays verbatim; BPCells bitpacks",
                 fontsize=11)
    for bar, v in zip(bars, values):
        ax.text(v + 0.4, bar.get_y() + bar.get_height() / 2, f"{v:.2f} MB",
                va="center", fontsize=9)
    ax.set_xlim(0, max(values) * 1.22)
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(0.5, -0.06,
             "BPCells' headline compression needs the integer conversion — on "
             "the doubles Seurat hands it, it gets 1.30×, not 5.97×.",
             ha="center", fontsize=8.5, color=_GREY)
    return fig


def main():
    _save(self_consistency(), "py_01_self_consistency.png")
    _save(memory(), "py_02_memory.png")
    _save(storage(), "py_03_storage.png")
    print(f"\n  Wrote figures to {FIGURES}")


if __name__ == "__main__":
    argparse.ArgumentParser(description="out-of-core tutorial figures").parse_args()
    main()
