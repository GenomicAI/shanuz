"""Generate all figures for the cell-hashing tutorial.

Runs pbmc_hashing_tutorial.run_full() and renders the shanuz-side figures to
tutorials/figures_hashing/, alongside the r_*.png the R verify script writes:
  * py_01_ridge.png          hashtag CLR enrichment per assigned sample
  * py_02_scatter.png        two hashtags, coloured by HTODemux global class
  * py_03_ncount_violin.png  total hashtag counts per global class
  * py_04_species_scatter.png  human vs mouse UMIs — the barnyard ground truth

Usage
-----
    python tutorials/generate_hashing_plots.py [--data-dir PATH]
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

from tutorials.pbmc_hashing_tutorial import run_full
from shanuz.plotting import ridge_plot, feature_scatter, vln_plot, _palette

FIGURES = Path(__file__).parent / "figures_hashing"
FIGURES.mkdir(exist_ok=True)

# Fixed colour for each global class, shared across the scatter and violin.
_GLOBAL_ORDER = ["Singlet", "Doublet", "Negative"]
_GLOBAL_COLORS = dict(zip(_GLOBAL_ORDER, _palette(3)))


def _save(fig, name):
    import matplotlib.pyplot as plt

    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path.name}")


def species_scatter(obj):
    """Human vs mouse UMIs per droplet, coloured by HTODemux global class.

    The barnyard the hashtags never saw: cells hugging an axis are single-species,
    cells off both axes are species-mixed multiplets. If demultiplexing works, the
    Doublet cloud should lean into the interior more than the Singlet cloud does.
    """
    import matplotlib.pyplot as plt

    meta = obj.meta_data
    human = np.log10(meta["nCount_human"].to_numpy(float) + 1)
    mouse = np.log10(meta["nCount_mouse"].to_numpy(float) + 1)
    glob = meta["HTO_classification.global"].to_numpy()

    fig, ax = plt.subplots(figsize=(6, 5))
    for cls in _GLOBAL_ORDER:
        m = glob == cls
        ax.scatter(human[m], mouse[m], s=3, alpha=0.4,
                   color=_GLOBAL_COLORS[cls], label=f"{cls} ({int(m.sum())})",
                   rasterized=True)
    ax.set_xlabel("log10 human UMIs + 1")
    ax.set_ylabel("log10 mouse UMIs + 1")
    ax.set_title("Cross-species ground truth (HTODemux colours)")
    leg = ax.legend(markerscale=3, frameon=False, fontsize=8)
    for lh in leg.legend_handles:
        lh.set_alpha(1.0)
    return fig


def main(data_dir=None):
    obj, tags, _ = run_full(data_dir=data_dir, verbose=False)

    # 1. Ridgeline of each hashtag's CLR enrichment, grouped by assigned sample.
    _save(ridge_plot(obj, features=tags, group_by="hash.ID", assay="HTO", ncol=3),
          "py_01_ridge.png")

    # 2. Two hashtags against each other, coloured by global class — singlets on
    #    the axes, doublets in the interior.
    _save(feature_scatter(obj, tags[0], tags[1], assay="HTO",
                          group_by="HTO_classification.global"),
          "py_02_scatter.png")

    # 3. Total hashtag counts per global class (doublets carry more signal).
    _save(vln_plot(obj, "nCount_HTO", group_by="HTO_classification.global"),
          "py_03_ncount_violin.png")

    # 4. The barnyard ground truth.
    _save(species_scatter(obj), "py_04_species_scatter.png")

    print(f"\n  Wrote figures to {FIGURES}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cell-hashing tutorial figures")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    main(data_dir=args.data_dir)
