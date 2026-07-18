"""Generate all figures for the leverage-score sketching tutorial.

Runs ifnb_sketch_tutorial.run_full() and renders the shanuz-side figures to
tutorials/figures_sketch/, alongside the r_*.png the R verify script writes:
  * py_01_leverage_by_type.png   leverage per cell type, commonest first
  * py_02_sketch_enrichment.png  cell-type share in the sketch vs uniform
  * py_03_leverage_vs_r.png      shanuz against R, both regimes
  * py_04_rarity.png             mean leverage against population size
  * py_05_projected_umap.png     the sketch's UMAP, extended to every cell

01 and 02 are the point of the tutorial: 01 shows leverage rising as a
population gets rarer, 02 shows the sketch acting on it against a same-size
uniform control. 03 is the fidelity check that caught the defect — before the
fix the shanuz series was a flat band with no relationship to R's.

Usage
-----
    python tutorials/generate_sketch_plots.py [--data-dir PATH]
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

from tutorials.ifnb_sketch_tutorial import (
    run_full, FIGURES, CELLTYPE, SKETCH_CELLS, _read_r_scores,
)
from shanuz.plotting import _palette

FIGURES.mkdir(exist_ok=True)

_PY_COLOR, _R_COLOR = _palette(2)
_UNIFORM_COLOR = "0.65"


def _save(fig, name):
    import matplotlib.pyplot as plt

    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path.name}")


def _r_leverage(cells, regime):
    """R's per-cell leverage, if the verify script has run."""
    path = FIGURES / "r_leverage.csv"
    if not path.exists():
        return None
    return _read_r_scores(path, f"r_leverage_{regime}", cells)


def leverage_by_type(summary):
    """Leverage per cell type, ordered commonest to rarest.

    The headline: the boxes climb steadily from left to right, which is the whole
    claim of the method — rare populations carry more influence per cell, so
    sampling by leverage keeps them.
    """
    import matplotlib.pyplot as plt

    scores = summary["leverage"]["exact"]
    groups = summary["celltype"]
    order = list(summary["enrichment"].index)          # already commonest-first

    fig, ax = plt.subplots(figsize=(8, 4.5))
    # Tick labels set separately: boxplot's own keyword was renamed from
    # `labels` to `tick_labels` in matplotlib 3.9, and this way works on both.
    ax.boxplot([scores[groups == g] for g in order],
               showfliers=False, patch_artist=True,
               boxprops=dict(facecolor=_PY_COLOR, alpha=0.6),
               medianprops=dict(color="black"))
    ax.set_xticks(np.arange(1, len(order) + 1))
    ax.set_xticklabels(order)
    ax.set_yscale("log")
    ax.axhline(scores.mean(), color="grey", lw=0.8, ls=":", label="overall mean")
    ax.set_ylabel("leverage (log scale)")
    ax.set_title("Shanuz — leverage by cell type (commonest first)")
    ax.tick_params(axis="x", rotation=45)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    return fig


def sketch_enrichment(summary):
    """Cell-type share in the sketch vs the full data — leverage against uniform.

    The uniform bars are the control that makes this readable: they scatter
    around 1.0, so any systematic lift in the leverage bars is the weighting
    doing work rather than a large sketch containing rare cells by accident.
    """
    import matplotlib.pyplot as plt

    lev = summary["sketch_composition"]["LeverageScore"]
    uni = summary["sketch_composition"]["Uniform"]
    names = list(lev.index)
    x = np.arange(len(names))
    width = 0.4

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(x - width / 2, lev["fold"], width, color=_PY_COLOR, label="LeverageScore")
    ax.bar(x + width / 2, uni.loc[names, "fold"], width, color=_UNIFORM_COLOR,
           label="Uniform (control)")
    ax.axhline(1.0, color="grey", lw=0.8, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha="right")
    ax.set_ylabel("fold enrichment in the sketch")
    ax.set_title(f"Shanuz — cell-type share, {SKETCH_CELLS}-cell sketch vs full data")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    return fig


def leverage_vs_r(summary):
    """shanuz against R, per cell, in both regimes.

    The exact regime is the fidelity check the defect failed: it should be a
    clean diagonal. Before the fix shanuz's scores were a flat band spanning a
    third of a decade with essentially no relationship to R's (Spearman 0.24).
    """
    import matplotlib.pyplot as plt

    cells = summary["cells"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, regime in zip(axes, ("exact", "sketched")):
        py = summary["leverage"][regime]
        r = _r_leverage(cells, regime)
        if r is None:
            ax.text(0.5, 0.5, "R reference not found\nrun ifnb_sketch_verify.R",
                    ha="center", va="center", transform=ax.transAxes, fontsize=9)
            ax.set_title(f"{regime} regime")
            continue
        ax.scatter(r, py, s=2, alpha=0.25, color=_PY_COLOR, edgecolors="none")
        lo = min(r.min(), py.min())
        hi = max(r.max(), py.max())
        ax.plot([lo, hi], [lo, hi], color=_R_COLOR, lw=1, ls="--", label="y = x")
        from scipy.stats import spearmanr
        ax.set_title(f"{regime} regime — Spearman {spearmanr(r, py).statistic:.3f}")
        ax.set_xlabel("R Seurat leverage")
        ax.set_ylabel("shanuz leverage")
        ax.legend(fontsize=8, frameon=False)
    fig.suptitle("Per-cell leverage, shanuz vs R Seurat", y=1.02)
    fig.tight_layout()
    return fig


def rarity(summary):
    """Mean leverage against population size — the method's claim in one panel."""
    import matplotlib.pyplot as plt

    table = summary["enrichment"]
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.scatter(table["n"], table["mean_leverage"], s=40, color=_PY_COLOR,
               edgecolors="none")
    for name, row in table.iterrows():
        ax.annotate(str(name), (row["n"], row["mean_leverage"]), fontsize=7,
                    xytext=(4, 3), textcoords="offset points")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("cells of this type (log scale)")
    ax.set_ylabel("mean leverage (log scale)")
    rho = summary["rarity_spearman"]["exact"]
    ax.set_title(f"Rarer types carry more leverage — Spearman {rho:+.3f}")
    fig.tight_layout()
    return fig


def projected_umap(obj, summary):
    """The sketch's UMAP, extended to every cell by project_data."""
    import matplotlib.pyplot as plt

    if "ref.umap" not in obj.reductions:
        return None
    coords = obj.reductions["ref.umap"].cell_embeddings
    groups = summary["celltype"]
    names = list(summary["enrichment"].index)
    colors = _palette(len(names))

    fig, ax = plt.subplots(figsize=(7, 5.5))
    for name, color in zip(names, colors):
        mask = groups == name
        ax.scatter(coords[mask, 0], coords[mask, 1], s=2, alpha=0.5, color=color,
                   edgecolors="none", label=f"{name} ({int(mask.sum())})")
    ax.set_xlabel("refUMAP_1")
    ax.set_ylabel("refUMAP_2")
    ax.set_title(f"Shanuz — {SKETCH_CELLS}-cell sketch's UMAP, projected to all "
                 f"{len(obj)} cells")
    ax.legend(fontsize=6, frameon=False, markerscale=3, loc="center left",
              bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout()
    return fig


def main(data_dir=None):
    obj, summary = run_full(data_dir=data_dir, verbose=False)

    _save(leverage_by_type(summary), "py_01_leverage_by_type.png")
    _save(sketch_enrichment(summary), "py_02_sketch_enrichment.png")
    _save(leverage_vs_r(summary), "py_03_leverage_vs_r.png")
    _save(rarity(summary), "py_04_rarity.png")
    umap_fig = projected_umap(obj, summary)
    if umap_fig is not None:
        _save(umap_fig, "py_05_projected_umap.png")

    if not (FIGURES / "r_leverage.csv").exists():
        print("\n  NOTE: r_leverage.csv absent — the comparison panel shows")
        print("  only the shanuz series. Run ifnb_sketch_verify.R for both.")
    print(f"\n  Wrote figures to {FIGURES}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ifnb sketching tutorial figures")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    main(data_dir=args.data_dir)
