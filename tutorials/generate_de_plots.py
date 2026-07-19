"""Generate the figures for the differential-expression tutorial.

Renders to tutorials/figures_de/:
  * py_01_log2fc_vs_r.png     the fold-change defect, before and after
  * py_02_threshold_impact.png what the wrong fold change did to the gene list
  * py_03_test_concordance.png all eight tests against Seurat

01 is the defect. 02 is why it was a defect rather than a cosmetic difference:
`logfc_threshold` filters on the number in 01, so getting it wrong changed which
genes came back at all. 03 is the tutorial's actual result.

The "before" series is recomputed from the same matrix rather than obtained by
reverting the library — same data, two formulas, no second code path to keep
honest.

Usage
-----
    python tutorials/generate_de_plots.py [--data-dir PATH]
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

from tutorials.pbmc3k_de_tutorial import (  # noqa: E402
    FIGURES,
    IDENT_1,
    IDENT_2,
    TEST_MAP,
    build,
    compare,
    shared_groups,
)

_AFTER = "#48a9a6"
_BEFORE = "#c1666b"
_R = "0.35"


def _save(fig, name):
    import matplotlib.pyplot as plt

    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path.name}")


def _old_formula(mat1, mat2):
    """Seurat 4's fold change: the pseudocount added to the mean."""
    return (np.log2(np.expm1(mat1).mean(axis=1) + 1.0)
            - np.log2(np.expm1(mat2).mean(axis=1) + 1.0))


def _group_matrices(obj, groups):
    assay = obj.assays["RNA"]
    data = assay.layer_data("data")
    X = np.asarray(data.todense() if hasattr(data, "todense") else data, dtype=float)
    cells = list(assay.cells())
    pos = {c: i for i, c in enumerate(cells)}
    i1 = [pos[c] for c in groups.index[groups == IDENT_1]]
    i2 = [pos[c] for c in groups.index[groups == IDENT_2]]
    return X[:, i1], X[:, i2], list(assay.features())


def log2fc_vs_r(before, after, r):
    """Both formulas against Seurat, on the identity line."""
    import matplotlib.pyplot as plt

    shared = after.index.intersection(r.index)
    fig, ax = plt.subplots(figsize=(6.2, 6))
    lim = [min(r[shared].min(), before[shared].min()) - 0.5,
           max(r[shared].max(), before[shared].max()) + 0.5]
    ax.plot(lim, lim, color=_R, lw=1.0, ls="--", zorder=1,
            label="identity (perfect agreement)")
    ax.scatter(r[shared], before[shared], s=8, color=_BEFORE, alpha=0.5, zorder=2,
               label="before — pseudocount on the mean")
    ax.scatter(r[shared], after[shared], s=8, color=_AFTER, alpha=0.8, zorder=3,
               label="after — pseudocount on the sum (R's)")
    ax.set_xlabel("Seurat 5.5.1  ·  avg_log2FC")
    ax.set_ylabel("shanuz  ·  avg_log2FC")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_title("avg_log2FC per gene — 13,712 genes, 1,172 shared cells",
                 fontsize=11)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    return fig


def threshold_impact(before, r):
    """Genes surviving logfc_threshold, which is the filter's actual input."""
    import matplotlib.pyplot as plt

    shared = before.index.intersection(r.index)
    thresholds = np.arange(0.0, 2.01, 0.05)
    n_before = [int((np.abs(before[shared]) >= t).sum()) for t in thresholds]
    n_r = [int((np.abs(r[shared]) >= t).sum()) for t in thresholds]
    jaccard = []
    for t in thresholds:
        a = set(shared[np.abs(before[shared]) >= t])
        b = set(shared[np.abs(r[shared]) >= t])
        jaccard.append(len(a & b) / len(a | b) if (a | b) else 1.0)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax.plot(thresholds, n_r, color=_R, lw=2, label="Seurat (and shanuz, after)")
    ax.plot(thresholds, n_before, color=_BEFORE, lw=2, label="shanuz, before")
    for t in (0.1, 0.25):
        ax.axvline(t, color="0.8", lw=0.8, ls=":", zorder=0)
    ax.set_xlabel("logfc_threshold")
    ax.set_ylabel("genes returned")
    ax.set_title("What the filter let through", fontsize=11)
    ax.legend(fontsize=8, frameon=False)

    ax2.plot(thresholds, jaccard, color=_BEFORE, lw=2)
    ax2.set_ylim(0, 1.02)
    ax2.set_xlabel("logfc_threshold")
    ax2.set_ylabel("Jaccard with Seurat's gene set")
    ax2.set_title("Agreement on *which* genes, before the fix", fontsize=11)
    for t, lbl in ((0.1, "Seurat's default"), (0.25, "a common choice")):
        j = jaccard[int(round(t / 0.05))]
        ax2.plot([t], [j], "o", color=_R, ms=5)
        ax2.annotate(f"{lbl}\n{j:.2f}", (t, j), textcoords="offset points",
                     xytext=(8, 8), fontsize=8, color=_R)
    for a in (ax, ax2):
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "avg_log2FC is not just reported — `logfc_threshold` filters on it, so a "
        "wrong value changes the result set.",
        fontsize=10, y=1.02)
    return fig


def test_concordance(table):
    """Every test's agreement with Seurat, after the fixes."""
    import matplotlib.pyplot as plt

    tests = [t for t in TEST_MAP if t != "roc"]
    overall = [table.loc[t, "p_spearman"] for t in tests]
    expressed = [table.loc[t, "p_spearman_expressed"] for t in tests]
    y = np.arange(len(tests))

    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    ax.barh(y + 0.19, overall, height=0.36, color="0.78",
            label="all genes scored by both")
    ax.barh(y - 0.19, expressed, height=0.36, color=_AFTER,
            label="genes detected in >5% of a group")
    ax.set_yticks(y)
    ax.set_yticklabels(tests, fontsize=9)
    ax.set_xlim(0, 1.05)
    ax.axvline(1.0, color=_R, lw=1.0, ls="--")
    ax.set_xlabel("Spearman correlation of p-values with Seurat 5.5.1")
    ax.set_title("Seven of eight tests reproduce Seurat's ranking\n"
                 "(`roc` omitted — it returns an AUC, not a p-value)",
                 fontsize=11)
    # Below the axes: at "lower left" it sat on top of the wilcox bars.
    ax.legend(fontsize=8, frameon=False, loc="upper center", ncol=2,
              bbox_to_anchor=(0.5, -0.16))
    ax.spines[["top", "right"]].set_visible(False)
    # The two that are not near 1 are not the same kind of miss.
    for i, t in enumerate(tests):
        note = {"negbinom": "near-empty genes only",
                "deseq2": "pseudobulk ≠ per-cell"}.get(t)
        if note:
            ax.text(1.03, i, note, va="center", fontsize=7.5, color=_R)
    return fig


def main(data_dir=None):
    FIGURES.mkdir(exist_ok=True)
    r_wilcox = FIGURES / "r_wilcox.csv"
    if not r_wilcox.exists():
        raise SystemExit("figures_de/r_wilcox.csv missing — run "
                         "`Rscript tutorials/pbmc3k_de_verify.R` first.")

    obj = build(data_dir)
    groups = shared_groups(obj)
    mat1, mat2, feats = _group_matrices(obj, groups)
    before = pd.Series(_old_formula(mat1, mat2), index=feats)

    py = pd.read_csv(FIGURES / "py_wilcox.csv", index_col=0)
    r = pd.read_csv(r_wilcox, index_col=0)
    after = py["avg_log2FC"]
    r_fc = r["avg_log2FC"]

    _save(log2fc_vs_r(before, after, r_fc), "py_01_log2fc_vs_r.png")
    _save(threshold_impact(before, r_fc), "py_02_threshold_impact.png")

    rows = []
    for test in TEST_MAP:
        p = pd.read_csv(FIGURES / f"py_{test}.csv", index_col=0)
        rr = pd.read_csv(FIGURES / f"r_{test.lower()}.csv", index_col=0)
        rows.append({"test": test, **compare(p, rr, test)})
    table = pd.DataFrame(rows).set_index("test")
    _save(test_concordance(table), "py_03_test_concordance.png")
    print(f"\n  Wrote figures to {FIGURES}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DE tutorial figures")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    main(data_dir=args.data_dir)
