"""Generate all figures for the Mixscape tutorial.

Runs thp1_mixscape_tutorial.run_full() and renders the shanuz-side figures to
tutorials/figures_mixscape/, alongside the r_*.png the R verify script writes:
  * py_01_perturb_score.png  IFNGR2 perturbation-score density (KO vs NP vs NT)
  * py_02_lda.png            MixscapeLDA map, coloured by global class
  * py_03_heatmap.png        DE heatmap, NT vs IFNGR2 KO, cells by KO probability
  * py_04_ko_rate.png        per-guide knockout rate — which edits actually took

Usage
-----
    python tutorials/generate_mixscape_plots.py [--data-dir PATH]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from tutorials.thp1_mixscape_tutorial import run_full, perturbation_table
from shanuz.plotting import plot_perturb_score, dim_plot, mixscape_heatmap, _palette

FIGURES = Path(__file__).parent / "figures_mixscape"
FIGURES.mkdir(exist_ok=True)

FOCUS_GENE = "IFNGR2"          # the vignette's canonical strong knockout


def _save(fig, name):
    import matplotlib.pyplot as plt

    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path.name}")


def ko_rate_bar(obj):
    """Per-guide knockout rate as a horizontal bar — the headline Mixscape result.

    Each guide's cells split into KO (the edit took) and NP (escaped); the bar is
    the KO fraction. Strong-phenotype guides sit near the top; guides whose effect
    is invisible at the RNA level fall to zero — not a failure of the guide, but of
    RNA to report it (several are checkpoint proteins the screen also measured).
    """
    import matplotlib.pyplot as plt

    tbl = perturbation_table(obj.meta_data)
    tbl = tbl.iloc[::-1]                      # highest rate at the top of the bar
    colors = [_palette(3)[0] if r >= 0.5 else _palette(3)[2] for r in tbl["ko_rate"]]

    fig, ax = plt.subplots(figsize=(6, 8))
    ax.barh(tbl["gene"], tbl["ko_rate"], color=colors)
    ax.axvline(0.5, color="black", linewidth=0.8, linestyle="--")
    ax.set_xlabel("knockout rate (KO / guide cells)")
    ax.set_xlim(0, 1)
    ax.set_title("Per-guide knockout rate (Mixscape)")
    for y, (rate, n) in enumerate(zip(tbl["ko_rate"], tbl["n_cells"])):
        ax.text(min(rate + 0.02, 0.98), y, f"{rate:.0%} (n={n})",
                va="center", fontsize=7)
    fig.tight_layout()
    return fig


def main(data_dir=None):
    obj, _ = run_full(data_dir=data_dir, verbose=False, do_lda=True)

    # 1. Why Mixscape split IFNGR2 the way it did: the score it thresholded on.
    _save(plot_perturb_score(obj, FOCUS_GENE, mixscape_class="mixscape_class"),
          "py_01_perturb_score.png")

    # 2. The supervised map — genuine knockouts pull away from the NT/NP cloud.
    _save(dim_plot(obj, reduction="lda", group_by="mixscape_class.global",
                   label=False, pt_size=3),
          "py_02_lda.png")

    # 3. The genes underneath the score, cells ordered by KO probability.
    _save(mixscape_heatmap(obj, ident_1="NT", ident_2=f"{FOCUS_GENE} KO",
                           max_genes=20),
          "py_03_heatmap.png")

    # 4. The headline: which guides actually knocked their gene out.
    _save(ko_rate_bar(obj), "py_04_ko_rate.png")

    print(f"\n  Wrote figures to {FIGURES}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mixscape tutorial figures")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    main(data_dir=args.data_dir)
