"""Cell-cycle & module scoring tutorial — Shanuz vs R Seurat on THP-1.

A Python port of Seurat's [cell-cycle vignette](https://satijalab.org/seurat/articles/cell_cycle_vignette)
and `AddModuleScore`, run on the THP-1 ECCITE-seq dataset (**GSE153056**,
Papalexi et al. 2021). THP-1 is a proliferating monocytic-leukemia line, so —
unlike the resting PBMCs of the earlier tutorials — it carries substantial S and
G2/M populations, which is what makes cell-cycle scoring a meaningful test rather
than "everything is G1". (We reuse the Mixscape tutorial's counts; the
perturbation labels are irrelevant here — only the raw RNA matters.)

The task these two functions solve: **score a gene program per cell against a
matched control set**, so the score reflects the program and not a cell's overall
depth or the genes' baseline abundance.

  * :func:`shanuz.add_module_score` ↔ ``AddModuleScore`` — the general primitive:
    a program's score is its genes' mean expression minus the mean of control
    genes drawn from the *same average-expression bins*, so highly- and
    lowly-expressed programs are put on the same footing.
  * :func:`shanuz.cell_cycle_scoring` ↔ ``CellCycleScoring`` — module scoring
    applied to the Tirosh 2016 S and G2/M gene sets (Seurat's
    ``cc.genes.updated.2019``), then a discrete ``Phase`` per cell: ``G1`` when
    both scores are ≤ 0, else whichever of S / G2M is larger.

Why this tutorial exists
------------------------
``add_module_score`` / ``cell_cycle_scoring`` landed with only synthetic
fixtures. This is the first time they meet a real dataset with a Seurat
reference. The wrinkle worth stating up front: both functions **sample control
genes at random** (binned by expression), and NumPy's RNG is not R's — so on
identical counts and identical gene lists the per-cell *scores* will not be
byte-identical. They will, however, correlate extremely tightly (the algorithm is
the same; only the control draw differs), and the discrete ``Phase`` — a
thresholding of those scores — is robust to the small score wobble. So the
comparison targets are **per-cell Phase concordance** and **score correlation**,
not exact score equality — the same "faithful port, RNG-driven residual" story as
``clara`` (hashing) and the KDE step (MULTI-seq).

Note on the data and the R comparison
-------------------------------------
Both tools read the *same* GEO counts the Mixscape tutorial caches, and — to rule
out gene-list drift as a source of divergence — the Python run writes the exact
S / G2M / module gene symbols it resolved against the assay to
``figures_cellcycle/*.txt``, which the R script reads back. So the only thing
left to differ is the control-gene RNG. The dataset also ships Papalexi's own
published ``Phase`` (from their Seurat run); we keep it as ``published_Phase`` for
a bonus external sanity check, but the controlled comparison is shanuz vs a fresh
R ``CellCycleScoring`` on identical input.

Usage
-----
    python tutorials/thp1_cellcycle_tutorial.py    # downloads ~66 MB (shared with Mixscape), writes gene lists

Then, for the side-by-side numbers and figures:

    Rscript tutorials/thp1_cellcycle_verify.R
    python  tutorials/generate_cellcycle_plots.py

References
----------
Tirosh I, Izar B, Prakadan SM, et al. (2016) **Dissecting the multicellular
ecosystem of metastatic melanoma by single-cell RNA-seq.** Science 352, 189-196.
https://doi.org/10.1126/science.aad0501

Papalexi E, Mimitou EP, Butler AW, et al. (2021) **Characterizing the molecular
regulation of inhibitory immune checkpoints with multimodal single-cell
screens.** Nature Genetics 53, 322-331. https://doi.org/10.1038/s41588-021-00778-2
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from shanuz.datasets import thp1_eccite
from shanuz.shanuz import create_shanuz_object
from shanuz.preprocessing import normalize_data
from shanuz.module_score import add_module_score, cell_cycle_scoring, CC_GENES

FIGURES = Path(__file__).parent / "figures_cellcycle"

# The columns cell_cycle_scoring writes, and the discrete phases.
PHASE_COL = "Phase"
S_COL = "S.Score"
G2M_COL = "G2M.Score"
PHASES = ("G1", "S", "G2M")

# A compact interferon-stimulated-gene program to demonstrate the general
# add_module_score primitive on something other than the cell cycle. Apt for
# THP-1, whose ECCITE-seq screen targets interferon-γ regulators.
IFN_PROGRAM = [
    "STAT1", "IRF1", "GBP1", "GBP2", "GBP4", "GBP5", "ISG15", "IFI6", "IFIT1",
    "IFIT3", "MX1", "OAS1", "OASL", "IFI44", "IFI44L", "IFITM1", "IFITM3",
    "PSMB8", "PSMB9", "TAP1", "B2M", "WARS", "UBE2L6", "IRF7", "STAT2",
]
IFN_NAME = "IFN.Response"


def section(title: str) -> None:
    print(f"\n{'=' * 70}\n{title}\n{'=' * 70}")


# ---------------------------------------------------------------------------
# Pure helpers (unit-tested without network in tests/test_cellcycle_tutorial.py)
# ---------------------------------------------------------------------------

def phase_concordance(a, b) -> float:
    """Fraction of cells assigned the same discrete phase by two tools.

    The headline cell-cycle metric — like Mixscape's class concordance, it is a
    per-cell agreement over a discretised call, robust to the small score wobble
    the control-gene RNG introduces. Order-sensitive.
    """
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape or a.size == 0:
        raise ValueError("both phase vectors must be the same non-zero length")
    return float(np.mean(a == b))


def score_correlation(a, b) -> dict:
    """Pearson and Spearman correlation between two per-cell score vectors.

    Because both tools compute the *same* module-score algorithm and differ only
    in the random control set, the continuous scores should correlate almost
    perfectly even though they are not identical — this is the number that shows
    the port is faithful. Returns ``{"pearson": .., "spearman": ..}``.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape or a.size == 0:
        raise ValueError("both score vectors must be the same non-zero length")
    from scipy.stats import pearsonr, spearmanr

    return {
        "pearson": float(pearsonr(a, b)[0]),
        "spearman": float(spearmanr(a, b)[0]),
    }


def phase_distribution(phases) -> pd.DataFrame:
    """Per-phase cell count and fraction, ordered G1 → S → G2M.

    A partition-only summary (needs just the phase labels), so R and Python are
    compared on identical terms. Phases absent from the data still appear, at 0.
    """
    phases = np.asarray(phases).astype(str)
    n = len(phases)
    rows = []
    for p in PHASES:
        c = int((phases == p).sum())
        rows.append({"phase": p, "n": c, "fraction": (c / n) if n else float("nan")})
    return pd.DataFrame(rows, columns=["phase", "n", "fraction"])


def build_scoreboard(rows: list[dict]) -> pd.DataFrame:
    """Assemble per-comparison score metrics into a tidy, ordered table."""
    cols = ["metric", "pearson", "spearman", "concordance"]
    df = pd.DataFrame(rows)
    return df[[c for c in cols if c in df.columns]]


def resolved_genes(obj, genes, assay: str = "RNA") -> list[str]:
    """The subset of ``genes`` actually present in the assay (order preserved).

    Both cell_cycle_scoring and add_module_score silently drop program genes the
    assay lacks; writing out exactly what survived lets the R run score the
    identical gene set, so the comparison isolates the control-gene RNG.
    """
    present = set(obj.assays[assay].features())
    return [g for g in genes if g in present]


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def load_thp1_object(data_dir=None, min_cells=3):
    """Load the THP-1 RNA counts and build the object, preserving the published phase.

    Reads the same GEO cache the Mixscape tutorial downloads. The dataset ships
    Papalexi's own ``Phase`` / ``S.Score`` / ``G2M.Score``; those column names
    would collide with what :func:`cell_cycle_scoring` writes, so they are kept
    under a ``published_`` prefix for the bonus external comparison.
    """
    rna_mat, rna_genes, _adt, _adt_names, meta, cells = thp1_eccite(data_dir=data_dir)
    keep = pd.DataFrame(index=list(cells))
    for col, out in ((PHASE_COL, "published_Phase"),
                     (S_COL, "published_S.Score"),
                     (G2M_COL, "published_G2M.Score")):
        if col in meta.columns:
            keep[out] = meta.loc[cells, col].to_numpy()
    obj = create_shanuz_object(
        counts=rna_mat, assay="RNA", min_cells=min_cells, min_features=0,
        project="thp1_cellcycle", feature_names=list(rna_genes),
        cell_names=list(cells), meta_data=keep,
    )
    return obj


def run_scoring(obj, seed=1):
    """Normalize, then score the cell cycle and the interferon program.

    ``NormalizeData`` → :func:`cell_cycle_scoring` (S.Score / G2M.Score / Phase) →
    :func:`add_module_score` for the IFN program. Writes the exact resolved gene
    lists to ``figures_cellcycle/`` for the R reference. Returns the resolved
    ``(s_genes, g2m_genes, ifn_genes)``.
    """
    normalize_data(obj, assay="RNA")

    s_genes = resolved_genes(obj, CC_GENES["s_genes"])
    g2m_genes = resolved_genes(obj, CC_GENES["g2m_genes"])
    ifn_genes = resolved_genes(obj, IFN_PROGRAM)

    cell_cycle_scoring(obj, s_features=s_genes, g2m_features=g2m_genes, seed=seed)
    add_module_score(obj, features={IFN_NAME: ifn_genes}, seed=seed)

    FIGURES.mkdir(exist_ok=True)
    (FIGURES / "s_genes.txt").write_text("\n".join(s_genes) + "\n")
    (FIGURES / "g2m_genes.txt").write_text("\n".join(g2m_genes) + "\n")
    (FIGURES / "ifn_genes.txt").write_text("\n".join(ifn_genes) + "\n")
    return s_genes, g2m_genes, ifn_genes


def summarize(obj, verbose=True) -> dict:
    """Summarise the Python-side scoring: phase distribution + published agreement."""
    phase = np.asarray(obj.meta_data[PHASE_COL]).astype(str)
    dist = phase_distribution(phase)

    out = {
        "n_cells": int(len(obj.cell_names())),
        "phase_distribution": dist,
        "mean_S": float(np.mean(obj.meta_data[S_COL])),
        "mean_G2M": float(np.mean(obj.meta_data[G2M_COL])),
    }
    # Bonus: agreement with Papalexi's published phase, if it came through.
    if "published_Phase" in obj.meta_data.columns:
        pub = np.asarray(obj.meta_data["published_Phase"]).astype(str)
        out["published_concordance"] = phase_concordance(phase, pub)

    if not verbose:
        return out

    section("Cell-cycle phase distribution (THP-1, a proliferating line)")
    print("    " + dist.to_string(index=False).replace("\n", "\n    "))
    if "published_concordance" in out:
        print(f"\n  Phase vs Papalexi's published call: "
              f"{out['published_concordance']:.4f} "
              f"(context — their pipeline, not a controlled comparison).")
    return out


def report_concordance(obj, r_calls_path=None, verbose=True) -> dict | None:
    """Compare Python vs R scoring, if ``thp1_cellcycle_verify.R`` has run.

    Reads ``figures_cellcycle/r_calls.csv`` (R's Phase + S.Score + G2M.Score +
    IFN score per cell) and reports per-cell Phase concordance and the score
    correlations. Returns ``None`` (with a hint) when the R calls are absent.
    """
    path = Path(r_calls_path) if r_calls_path else FIGURES / "r_calls.csv"
    if not path.exists():
        if verbose:
            print(f"\n  [concordance] {path.name} not found — run "
                  "`Rscript tutorials/thp1_cellcycle_verify.R` first.")
        return None

    r = pd.read_csv(path).set_index("cell")
    cells = obj.cell_names()
    r = r.reindex(cells)
    meta = obj.meta_data

    out = {
        "phase_concordance": phase_concordance(
            np.asarray(meta[PHASE_COL]).astype(str),
            r["R_Phase"].astype(str).to_numpy()),
        "s_score": score_correlation(meta[S_COL].to_numpy(), r["R_S_Score"].to_numpy()),
        "g2m_score": score_correlation(meta[G2M_COL].to_numpy(), r["R_G2M_Score"].to_numpy()),
        "ifn_score": score_correlation(meta[IFN_NAME].to_numpy(), r["R_IFN"].to_numpy()),
    }

    if verbose:
        section("R-vs-Python concordance (shared counts & gene lists)")
        board = build_scoreboard([
            {"metric": "S.Score", "pearson": round(out["s_score"]["pearson"], 4),
             "spearman": round(out["s_score"]["spearman"], 4)},
            {"metric": "G2M.Score", "pearson": round(out["g2m_score"]["pearson"], 4),
             "spearman": round(out["g2m_score"]["spearman"], 4)},
            {"metric": IFN_NAME, "pearson": round(out["ifn_score"]["pearson"], 4),
             "spearman": round(out["ifn_score"]["spearman"], 4)},
        ])
        print("    " + board.to_string(index=False).replace("\n", "\n    "))
        print(f"\n  Per-cell Phase concordance: {out['phase_concordance']:.4f} "
              f"({int(out['phase_concordance'] * len(cells))}/{len(cells)} cells).")
        print("  pearson/spearman: score agreement — high despite the control-gene")
        print("                    RNG differing between NumPy and R (scores are not")
        print("                    expected to be identical, only to track).")
    return out


def run_full(data_dir=None, verbose=True, seed=1):
    """Load, normalize, score cell cycle + IFN module, summarize, compare to R."""
    t0 = time.time()
    if verbose:
        section("Loading THP-1 ECCITE-seq RNA (cell-cycle & module scoring)")
    obj = load_thp1_object(data_dir=data_dir)
    if verbose:
        print(f"  {len(obj.assays['RNA'].features())} genes x "
              f"{len(obj.cell_names())} cells")

    t1 = time.time()
    s_genes, g2m_genes, ifn_genes = run_scoring(obj, seed=seed)
    if verbose:
        print(f"  scored: {len(s_genes)} S + {len(g2m_genes)} G2M cell-cycle genes, "
              f"{len(ifn_genes)}-gene IFN program ({time.time() - t1:.1f}s)")
        print("  gene lists written to figures_cellcycle/ for the R reference")

    summary = summarize(obj, verbose=verbose)
    report_concordance(obj, verbose=verbose)

    if verbose:
        section("Summary")
        print(f"  Total runtime: {time.time() - t0:.1f}s")
    return obj, summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="THP-1 cell-cycle & module-score tutorial")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    run_full(data_dir=args.data_dir)
