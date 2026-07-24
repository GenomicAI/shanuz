#!/usr/bin/env python
"""Figures for the anchor-internals vignette.

Reads only the artifacts the tutorial and the R verify script already wrote —
no recomputation, so a figure can never disagree with the numbers in the text.

    python tutorials/anchors_tutorial.py
    Rscript tutorials/anchors_verify.R
    python tutorials/generate_anchors_plots.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures_anchors"

BLUE, ORANGE, GREY = "#3b6ea5", "#d1873b", "#8a8a8a"


def _load(reduction):
    r = pd.read_csv(FIG / f"r_anchors_{reduction}.csv")
    p = pd.read_csv(FIG / f"py_anchors_{reduction}.csv")
    return r, p


def fig_agreement():
    """Anchor overlap and score agreement, one column per reduction."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for col, reduction in enumerate(("cca", "rpca")):
        r, p = _load(reduction)
        rp = set(zip(r["cell1"], r["cell2"]))
        pp = set(zip(p["cell1"], p["cell2"]))
        shared = rp & pp

        ax = axes[0, col]
        only_r, only_p = len(rp - pp), len(pp - rp)
        ax.bar(["Seurat only", "both", "shanuz only"], [only_r, len(shared), only_p],
               color=[ORANGE, BLUE, GREY])
        ax.set_title(f"{reduction.upper()} — anchor pairs\n"
                     f"{len(shared):,} of Seurat's {len(rp):,} recovered "
                     f"({100*len(shared)/len(rp):.1f}%)")
        ax.set_ylabel("anchors")
        for i, v in enumerate([only_r, len(shared), only_p]):
            ax.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=9)

        ax = axes[1, col]
        rs = r.set_index(["cell1", "cell2"])["score"]
        ps = p.set_index(["cell1", "cell2"])["score"]
        order = list(shared)
        x = np.array([ps.loc[k] for k in order])
        y = np.array([rs.loc[k] for k in order])
        ax.scatter(y, x, s=4, alpha=0.25, color=BLUE, edgecolors="none")
        ax.plot([0, 1], [0, 1], color=GREY, lw=1, ls="--")
        ax.set_xlabel("Seurat anchor score")
        ax.set_ylabel("shanuz anchor score")
        ax.set_title(f"scores on the {len(order):,} shared anchors\n"
                     f"r = {np.corrcoef(x, y)[0, 1]:.5f}, "
                     f"{100*np.mean(np.isclose(x, y, atol=1e-9)):.1f}% identical")
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.03, 1.03)
    fig.tight_layout()
    fig.savefig(FIG / "py_01_anchor_agreement.png", dpi=150)
    plt.close(fig)


def fig_correction():
    """The corrected expression of the query half, both tools side by side."""
    fig, ax = plt.subplots(figsize=(8, 4.6))
    labels, py_vals, r_vals = [], [], []
    for reduction in ("cca", "rpca"):
        pj = json.loads((FIG / f"py_summary_{reduction}.json").read_text())
        rj = json.loads((FIG / f"r_summary_{reduction}.json").read_text())
        labels += [f"{reduction.upper()}\nmean |Δ|", f"{reduction.upper()}\nfrac nonzero"]
        py_vals += [pj["delta_mean_abs"], pj["delta_frac_nonzero"]]
        r_vals += [rj["delta_mean_abs"], rj["delta_frac_nonzero"]]

    x = np.arange(len(labels))
    ax.bar(x - 0.19, r_vals, 0.38, label="Seurat 5.5.1", color=ORANGE)
    ax.bar(x + 0.19, py_vals, 0.38, label="shanuz", color=BLUE)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_yscale("log")
    ax.set_ylabel("value (log scale)")
    ax.set_title("Batch correction applied to the query half\n"
                 "(the reference is copied through untouched by both tools)")
    ax.legend()
    for xi, (a, b) in enumerate(zip(r_vals, py_vals)):
        ax.text(xi - 0.19, a, f"{a:.4f}", ha="center", va="bottom", fontsize=7.5)
        ax.text(xi + 0.19, b, f"{b:.4f}", ha="center", va="bottom", fontsize=7.5)
    fig.tight_layout()
    fig.savefig(FIG / "py_02_correction.png", dpi=150)
    plt.close(fig)


def main():
    need = [FIG / f"r_anchors_{r}.csv" for r in ("cca", "rpca")]
    missing = [p.name for p in need if not p.exists()]
    if missing:
        print(f"Missing {missing}. Run the tutorial and anchors_verify.R first.")
        return 1
    fig_agreement()
    fig_correction()
    print(f"Wrote py_01_anchor_agreement.png and py_02_correction.png to {FIG.name}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
