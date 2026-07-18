"""Generate all figures for the dimensional-reduction extras tutorial.

Runs pbmc3k_dimreduc_tutorial.run_full() and renders the shanuz-side figures to
tutorials/figures_dimreduc/, alongside the r_*.png the R verify script writes:
  * py_01_jackstraw_scores.png  per-PC ScoreJackStraw, both tools
  * py_02_nsig_compare.png      features below the significance threshold, per PC
  * py_03_elbow.png             the elbow plot — the other way to pick a cutoff
  * py_04_tsne.png              t-SNE coloured by LYZ (comparable across tools)
  * py_05_ica.png               independent components 1 & 2

The first two are the point of the tutorial: 01 shows both tools' per-PC scores
falling off the same cliff after PC 13, and 02 shows the per-feature counts
underneath them. Both were badly divergent before this tutorial found and fixed
two defects in ``jack_straw`` / ``score_jackstraw`` — see dimreduc_vignette.md.

Usage
-----
    python tutorials/generate_dimreduc_plots.py [--data-dir PATH]
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

from tutorials.pbmc3k_dimreduc_tutorial import (
    run_full, FIGURES, JS_DIMS, ALPHA, SCORE_THRESH, TSNE_DIMS,
)
from shanuz.plotting import _palette

FIGURES.mkdir(exist_ok=True)

_PY_COLOR, _R_COLOR = _palette(2)
# A score of exactly 0 has no logarithm; clamp to the smallest positive double so
# a saturated PC plots at the top of the axis instead of vanishing from it.
_FLOOR = np.nextafter(0, 1)


def _save(fig, name):
    import matplotlib.pyplot as plt

    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path.name}")


def _r_pcs():
    """R's per-PC JackStraw reference, if the verify script has run."""
    path = FIGURES / "r_jackstraw_pcs.csv"
    return pd.read_csv(path) if path.exists() else None


def jackstraw_scores(summary):
    """Per-PC significance score, both tools, on a -log10 axis.

    The headline figure: the two curves fall off the same cliff after PC 13, so
    both tools tell an analyst to keep the same PCs. Before the JackStraw fix
    shanuz's curve never came down at all — it sat above 110 on every PC, noise
    included, because the null was too tight and the KS-test aggregation was far
    too sensitive.
    """
    import matplotlib.pyplot as plt

    py = np.maximum(summary["pc_scores"], _FLOOR)
    pcs = np.arange(1, py.size + 1)
    r = _r_pcs()

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(pcs, -np.log10(py), "o-", color=_PY_COLOR, label="shanuz")
    if r is not None:
        ax.plot(r["PC"], -np.log10(np.maximum(r["R_Score"], _FLOOR)), "s--",
                color=_R_COLOR, label="R Seurat")
    ax.axhline(-np.log10(ALPHA), color="grey", lw=0.8, ls=":",
               label=f"alpha = {ALPHA}")
    ax.set_xlabel("PC")
    ax.set_ylabel(r"$-\log_{10}$(score)")
    ax.set_title("ScoreJackStraw per PC — higher is more significant")
    ax.set_xticks(pcs[::2])
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    return fig


def nsig_compare(summary):
    """Features below the significance threshold, per PC — where the gap starts.

    Computed the same way from either tool's per-feature p-value matrix, so
    unlike the aggregate score this is genuinely like-for-like. The two agree on
    the leading PCs (real signal) and part company on the trailing ones, which
    locates the difference in the null rather than in the aggregation.
    """
    import matplotlib.pyplot as plt

    table = summary["pc_table"]
    r = _r_pcs()

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    pcs = table["PC"].to_numpy()
    width = 0.4
    ax.bar(pcs - width / 2, table["py_n_sig_features"], width,
           color=_PY_COLOR, label="shanuz")
    if r is not None:
        ax.bar(pcs + width / 2, r["R_n_sig_features"], width,
               color=_R_COLOR, label="R Seurat")
    ax.set_xlabel("PC")
    ax.set_ylabel(f"features with p <= {SCORE_THRESH:g}")
    ax.set_title("JackStraw: per-feature significance, of 2000 variable features")
    ax.set_xticks(pcs[::2])
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    return fig


def elbow(obj, ndims=30):
    """Standard deviation per PC — the cheap heuristic JackStraw is meant to replace."""
    import matplotlib.pyplot as plt

    stdev = np.asarray(obj.reductions["pca"].stdev, dtype=float)[:ndims]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(np.arange(1, stdev.size + 1), stdev, "o-", color=_PY_COLOR, ms=4)
    ax.set_xlabel("PC")
    ax.set_ylabel("standard deviation")
    ax.set_title("Shanuz — elbow plot")
    fig.tight_layout()
    return fig


def tsne_marker(obj, gene="LYZ"):
    """t-SNE coloured by a marker gene.

    Coloured by expression rather than by cluster on purpose: cluster labels are
    arbitrary integers that would not correspond between the two tools, whereas
    this gene's value on this cell is the same number on both sides — so the
    panels compare even though the coordinates do not.
    """
    import matplotlib.pyplot as plt

    coords = obj.reductions["tsne"].cell_embeddings
    assay = obj.assays["RNA"]
    row = assay.features().index(gene)
    expr = np.asarray(assay.layer_data("data")[row, :].todense()).ravel()

    fig, ax = plt.subplots(figsize=(6, 5))
    order = np.argsort(expr)              # high-expressing cells drawn on top
    sc = ax.scatter(coords[order, 0], coords[order, 1], c=expr[order], s=4,
                    cmap="viridis", edgecolors="none")
    fig.colorbar(sc, ax=ax, label=f"{gene} (log-normalized)")
    ax.set_xlabel("tSNE_1")
    ax.set_ylabel("tSNE_2")
    ax.set_title(f"Shanuz — t-SNE (PC 1-{TSNE_DIMS}), {gene}")
    fig.tight_layout()
    return fig


def ica_scatter(obj):
    """Independent components 1 & 2."""
    import matplotlib.pyplot as plt

    ica = obj.reductions["ica"].cell_embeddings
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(ica[:, 0], ica[:, 1], s=4, alpha=0.5, color=_PY_COLOR,
               edgecolors="none")
    ax.set_xlabel("ICA_1")
    ax.set_ylabel("ICA_2")
    ax.set_title("Shanuz — ICA components 1 & 2")
    fig.tight_layout()
    return fig


def main(data_dir=None):
    obj, summary = run_full(data_dir=data_dir, verbose=False)

    _save(jackstraw_scores(summary), "py_01_jackstraw_scores.png")
    _save(nsig_compare(summary), "py_02_nsig_compare.png")
    _save(elbow(obj), "py_03_elbow.png")
    _save(tsne_marker(obj), "py_04_tsne.png")
    _save(ica_scatter(obj), "py_05_ica.png")

    if _r_pcs() is None:
        print("\n  NOTE: r_jackstraw_pcs.csv absent — the comparison figures show")
        print("  only the shanuz series. Run pbmc3k_dimreduc_verify.R for both.")
    print(f"\n  Wrote figures to {FIGURES}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PBMC 3k dim-reduction tutorial figures")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    main(data_dir=args.data_dir)
