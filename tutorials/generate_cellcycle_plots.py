"""Generate all figures for the cell-cycle & module-score tutorial.

Runs thp1_cellcycle_tutorial.run_full() and renders the shanuz-side figures to
tutorials/figures_cellcycle/, alongside the r_*.png the R verify script writes:
  * py_01_score_scatter.png   S.Score vs G2M.Score, coloured by assigned phase
  * py_02_phase_bar.png       cell count per phase (G1/S/G2M)
  * py_03_ifn_hist.png        the IFN-response module score distribution (AddModuleScore)

Usage
-----
    python tutorials/generate_cellcycle_plots.py [--data-dir PATH]
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

from tutorials.thp1_cellcycle_tutorial import (
    run_full, PHASE_COL, S_COL, G2M_COL, PHASES, IFN_NAME,
)
from shanuz.plotting import _palette

FIGURES = Path(__file__).parent / "figures_cellcycle"
FIGURES.mkdir(exist_ok=True)

_PHASE_COLOR = dict(zip(PHASES, _palette(3)))


def _save(fig, name):
    import matplotlib.pyplot as plt

    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path.name}")


def score_scatter(meta):
    """S.Score vs G2M.Score, coloured by the assigned phase — the classic read."""
    import matplotlib.pyplot as plt

    s = np.asarray(meta[S_COL], dtype=float)
    g2m = np.asarray(meta[G2M_COL], dtype=float)
    phase = np.asarray(meta[PHASE_COL]).astype(str)

    fig, ax = plt.subplots(figsize=(6, 5))
    for p in PHASES:
        m = phase == p
        ax.scatter(s[m], g2m[m], s=4, alpha=0.6, label=f"{p} (n={int(m.sum())})",
                   color=_PHASE_COLOR[p], edgecolors="none")
    ax.axhline(0, color="grey", lw=0.6)
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_xlabel("S.Score")
    ax.set_ylabel("G2M.Score")
    ax.set_title("Shanuz — cell-cycle scores by phase")
    ax.legend(fontsize=8, frameon=False, markerscale=2)
    fig.tight_layout()
    return fig


def phase_bar(summary):
    """Cell count per phase (G1/S/G2M)."""
    import matplotlib.pyplot as plt

    dist = summary["phase_distribution"]
    fig, ax = plt.subplots(figsize=(5, 4))
    bars = ax.bar(dist["phase"], dist["n"],
                  color=[_PHASE_COLOR[p] for p in dist["phase"]])
    ax.set_ylabel("cells")
    ax.set_title("Shanuz — phase distribution")
    ax.bar_label(bars, fmt="%d", fontsize=8, padding=2)
    fig.tight_layout()
    return fig


def ifn_hist(meta):
    """Distribution of the interferon-response module score (add_module_score)."""
    import matplotlib.pyplot as plt

    ifn = np.asarray(meta[IFN_NAME], dtype=float)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(ifn, bins=60, color=_palette(1)[0], alpha=0.85)
    ax.axvline(0, color="grey", lw=0.8, ls="--")
    ax.set_xlabel(f"{IFN_NAME} module score")
    ax.set_ylabel("cells")
    ax.set_title("Shanuz — add_module_score: interferon-response program")
    fig.tight_layout()
    return fig


def main(data_dir=None):
    obj, summary = run_full(data_dir=data_dir, verbose=False)
    meta = obj.meta_data

    _save(score_scatter(meta), "py_01_score_scatter.png")
    _save(phase_bar(summary), "py_02_phase_bar.png")
    _save(ifn_hist(meta), "py_03_ifn_hist.png")

    print(f"\n  Wrote figures to {FIGURES}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="THP-1 cell-cycle tutorial figures")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    main(data_dir=args.data_dir)
