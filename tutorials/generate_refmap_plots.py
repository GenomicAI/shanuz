"""Generate all figures for the reference-mapping tutorial.

Runs panc8_reference_mapping_tutorial.run_full(do_umap=True) and renders the
shanuz-side figures to tutorials/figures_refmap/, alongside the r_*.png the R
verify script writes:
  * py_01_reference_umap_celltype.png    the reference atlas UMAP, by cell type
  * py_02_query_projected_predicted.png  the query projected into the reference
                                         UMAP, coloured by transferred label
  * py_03_query_projected_truth.png      the same projection, by the true label
  * py_04_perclass_recall.png            per-cell-type transfer recall vs support

Usage
-----
    python tutorials/generate_refmap_plots.py [--data-dir PATH]
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

from tutorials.panc8_reference_mapping_tutorial import (
    run_full, CELLTYPE, REFERENCE_TECH, QUERY_TECH,
)
from shanuz.plotting import dim_plot, _palette

FIGURES = Path(__file__).parent / "figures_refmap"
FIGURES.mkdir(exist_ok=True)


def _save(fig, name):
    import matplotlib.pyplot as plt

    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path.name}")


def perclass_recall_bars(summary):
    """Per-cell-type transfer recall, ordered by support.

    The honest story of a single-technology reference in one panel: the abundant
    cell types (alpha, beta, ductal, …) are recovered near-perfectly, while the
    rare types with only a handful of reference cells (epsilon, schwann) are
    annotated noisily. Bars are recall; each is labelled with the cell type's
    support so the low-recall tail is visibly the low-support tail.
    """
    import matplotlib.pyplot as plt

    per_class = summary["per_class"]
    types = list(per_class["celltype"])
    recall = per_class["recall"].to_numpy()
    support = per_class["support"].to_numpy()

    x = np.arange(len(types))
    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.bar(x, recall, 0.7, color=_palette(1)[0])
    ax.set_xticks(x)
    ax.set_xticklabels(types, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("recall (fraction of true type recovered)")
    ax.set_ylim(0, 1.08)
    ax.set_title(f"Label transfer {REFERENCE_TECH} → {QUERY_TECH} — recall per cell type")
    ax.bar_label(bars, labels=[f"n={s}" for s in support], fontsize=7, padding=2)
    fig.tight_layout()
    return fig


def main(data_dir=None):
    reference, query, _anchors, predictions, summary = run_full(
        data_dir=data_dir, verbose=False, do_umap=True)
    # Carry the transferred label onto the query so the projection can be
    # coloured by it (as Seurat's MapQuery writes predicted.id onto the query).
    query.meta_data["predicted.id"] = (
        predictions["predicted.id"].reindex(query.cell_names()).to_numpy())

    # 1. The reference atlas: the annotated space the query is mapped into.
    _save(dim_plot(reference, reduction="umap", group_by=CELLTYPE, label=True,
                   pt_size=2, label_size=7,
                   title=f"Reference ({REFERENCE_TECH}) — by cell type"),
          "py_01_reference_umap_celltype.png")

    # 2. The query projected into the reference UMAP, by transferred label.
    _save(dim_plot(query, reduction="ref.umap", group_by="predicted.id", label=True,
                   pt_size=2, label_size=7,
                   title=f"Query ({QUERY_TECH}) projected — predicted labels"),
          "py_02_query_projected_predicted.png")

    # 3. The same projection by the query's true label — the two should match
    #    wherever the transfer succeeded.
    _save(dim_plot(query, reduction="ref.umap", group_by=CELLTYPE, label=True,
                   pt_size=2, label_size=7,
                   title=f"Query ({QUERY_TECH}) projected — true labels"),
          "py_03_query_projected_truth.png")

    # 4. The headline: where transfer works (abundant types) vs where it is noisy
    #    (rare, low-support types).
    _save(perclass_recall_bars(summary), "py_04_perclass_recall.png")

    print(f"\n  Wrote figures to {FIGURES}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="panc8 reference-mapping figures")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    main(data_dir=args.data_dir)
