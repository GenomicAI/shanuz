"""Pooled CRISPR screen tutorial — Mixscape with Shanuz (CalcPerturbSig + RunMixscape).

A Python port of Seurat's Mixscape vignette
(https://satijalab.org/seurat/articles/mixscape_vignette) on the THP-1 ECCITE-seq
screen (Papalexi, Mimitou et al. 2021, GSE153056): 20,729 stimulated THP-1 cells,
each carrying one guide RNA against one of 25 immune-regulatory genes (plus
non-targeting controls), profiled with paired RNA + protein.

The problem Mixscape solves: **carrying a guide is not the same as being
perturbed.** A cell can pick up a ``STAT1`` guide and still escape the knockout —
the edit fails, one allele survives — so the cells labelled ``STAT1`` are a
*mixture* of true knockouts (KO) and non-perturbed escapers (NP) that look just
like controls. Averaging over that mixture dilutes, and can entirely mask, the
phenotype. Mixscape separates the two so downstream analysis runs on genuinely
perturbed cells.

It demonstrates Shanuz's three Mixscape functions side by side with Seurat's:
  * :func:`shanuz.calc_perturb_sig` — Seurat's ``CalcPerturbSig``: subtract from
    each cell the mean of its nearest non-targeting (NT) neighbours, cancelling
    the technical variation guide assignment is confounded by, and leaving the
    **local perturbation signature**,
  * :func:`shanuz.run_mixscape` — Seurat's ``RunMixscape``: per target gene, an
    iterative two-component Gaussian mixture over the perturbation score splits
    KO from NP, so each guide gets a KO *rate* rather than an assumption,
  * :func:`shanuz.mixscape_lda` — Seurat's ``MixscapeLDA``: one supervised map on
    which each guide population forms its own cloud, the classic Mixscape LDA plot.

Why this tutorial exists
------------------------
Every Mixscape feature landed after PR #10 and had only been checked against
synthetic fixtures with a *known* KO/NP truth. This is the first time
``calc_perturb_sig`` / ``run_mixscape`` / ``mixscape_lda`` meet real screen data
with a Seurat reference. Mixscape is a stochastic, multi-stage pipeline (kNN,
per-gene differential expression, an EM-refined mixture), so the comparison
target is **not** byte-identical calls but *do the two tools recover the same
biology*: the same guides called strong-effect, KO/NP/NT proportions in step, and
per-cell class concordance, all reported by :func:`report_concordance` once
``thp1_mixscape_verify.R`` has run.

Note on the data and the R comparison
-------------------------------------
These are the raw GEO count matrices (cDNA + ADT + the published per-cell
metadata), the same bytes Shanuz and Seurat both read — *not* SeuratData's
pre-built ``thp1.eccite`` object, whose internal processing has no clean
cross-language form. The metadata already carries each cell's guide assignment
(``gene`` = target or ``NT``, ``guide_ID``, ``replicate``, cell-cycle ``Phase``),
so both tools start from the identical annotated state. To keep the perturbation
signature on one shared gene basis, the Python run writes the variable-feature set
it selected to ``figures_mixscape/hvg_features.txt`` and the R script reads it
back — so the only divergences left are the genuinely method-level ones.

Usage
-----
    python tutorials/thp1_mixscape_tutorial.py [--data-dir PATH]

The dataset (~66 MB) downloads automatically to ~/.shanuz_data/thp1_eccite.
Then, for the side-by-side numbers and figures:

    Rscript tutorials/thp1_mixscape_verify.R
    python  tutorials/generate_mixscape_plots.py

References
----------
Papalexi E, Mimitou EP, Butler AW, et al. (2021)
**Characterizing the molecular regulation of inhibitory immune checkpoints with
multimodal single-cell screens.** Nature Genetics 53, 322-331.
https://doi.org/10.1038/s41588-021-00778-2
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
from shanuz.assay5 import create_assay5_object
from shanuz.preprocessing import normalize_data, find_variable_features, scale_data
from shanuz.reduction import run_pca
from shanuz.mixscape import calc_perturb_sig, run_mixscape, mixscape_lda

FIGURES = Path(__file__).parent / "figures_mixscape"

# Metadata columns carried onto the object — the guide annotation the screen
# ships with, plus cell-cycle for optional confounder checks.
_META_COLS = [
    "gene", "guide_ID", "replicate", "crispr", "Phase",
    "S.Score", "G2M.Score", "percent.mito",
]

# Mixscape prep mirrors the vignette: signature over the first 40 PCs, each cell
# referenced to its 20 nearest NT controls, found within its own replicate.
NDIMS = 40
NUM_NEIGHBORS = 20
SPLIT_BY = "replicate"


def section(title: str) -> None:
    print(f"\n{'=' * 70}\n{title}\n{'=' * 70}")


# ---------------------------------------------------------------------------
# Pure helpers (unit-tested without network in tests/test_mixscape_tutorial.py)
# ---------------------------------------------------------------------------

def perturbation_table(
    meta,
    labels: str = "gene",
    global_col: str = "mixscape_class.global",
    nt_class: str = "NT",
    prtb_type: str = "KO",
) -> pd.DataFrame:
    """Per-target-gene knockout summary from mixscape metadata.

    Returns a DataFrame with one row per target gene (``nt_class`` excluded) and
    columns ``gene`` / ``n_cells`` / ``n_ko`` / ``n_np`` / ``ko_rate``, sorted by
    ``ko_rate`` descending — the "which guides actually knocked their gene out"
    table that separates a strong perturbation from one whose cells nearly all
    escaped.
    """
    g = pd.Series(np.asarray(meta[labels]), dtype=object)
    gc = pd.Series(np.asarray(meta[global_col]), dtype=object)
    rows = []
    for gene in sorted(x for x in set(g) if isinstance(x, str) and x != nt_class):
        m = (g == gene).to_numpy()
        n = int(m.sum())
        n_ko = int((gc[m] == prtb_type).sum())
        n_np = int((gc[m] == "NP").sum())
        rows.append((gene, n, n_ko, n_np, n_ko / n if n else np.nan))
    tbl = pd.DataFrame(rows, columns=["gene", "n_cells", "n_ko", "n_np", "ko_rate"])
    return tbl.sort_values("ko_rate", ascending=False).reset_index(drop=True)


def responsive_genes(table: pd.DataFrame, min_ko_rate: float = 0.5) -> list[str]:
    """Genes whose knockout 'took' in at least ``min_ko_rate`` of their cells."""
    return list(table.loc[table["ko_rate"] >= min_ko_rate, "gene"])


def call_concordance(a, b) -> float:
    """Fraction of positions where two per-cell call vectors agree."""
    a = np.asarray(a, dtype=object)
    b = np.asarray(b, dtype=object)
    if a.shape != b.shape or a.size == 0:
        raise ValueError("call vectors must be the same non-zero length")
    return float(np.mean(a == b))


def variable_feature_list(obj, assay: str = "RNA") -> list[str]:
    """The assay's selected variable features, as a plain list of names."""
    a = obj.assays[assay]
    vf = getattr(a, "variable_features", None)
    if vf is None:
        vf = getattr(a, "var_features", None)
    return list(vf) if vf is not None else []


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def load_screen_object(data_dir=None, min_cells=3):
    """Load GSE153056, build the RNA object, attach ADT and the guide metadata.

    Returns the :class:`~shanuz.Shanuz` object carrying the ``gene`` / ``guide_ID``
    / ``replicate`` / ``Phase`` metadata mixscape reads, plus an ``ADT`` assay of
    the four surface proteins (the screen is CITE-seq) and ``nCount_ADT``.
    """
    rna, genes, adt, adt_names, meta, cells = thp1_eccite(data_dir=data_dir)

    keep_meta = meta[[c for c in _META_COLS if c in meta.columns]].copy()
    obj = create_shanuz_object(
        counts=rna, assay="RNA", min_cells=min_cells, min_features=0,
        project="thp1_eccite", feature_names=genes, cell_names=cells,
        meta_data=keep_meta,
    )
    kept = obj.cell_names()
    cpos = {c: i for i, c in enumerate(cells)}
    idx = [cpos[c] for c in kept]

    adt_kept = adt[:, idx].tocsc()
    obj.assays["ADT"] = create_assay5_object(
        counts=adt_kept, feature_names=list(adt_names), cell_names=kept, key="adt_",
    )
    obj.meta_data["nCount_ADT"] = np.asarray(adt_kept.sum(axis=0)).ravel()
    return obj


def prep_reduction(obj, n_hvg=2000, n_pcs=50, verbose=False):
    """Standard RNA prep up to PCA — the reduction ``calc_perturb_sig`` needs.

    Normalize (LogNormalize) → variable features → scale → PCA, exactly the
    vignette's ``NormalizeData %>% FindVariableFeatures %>% ScaleData %>% RunPCA``.
    Returns the selected variable-feature list (the shared signature basis).
    """
    normalize_data(obj, assay="RNA")
    find_variable_features(obj, assay="RNA", nfeatures=n_hvg)
    scale_data(obj, assay="RNA")
    run_pca(obj, assay="RNA", n_pcs=n_pcs)
    return variable_feature_list(obj, assay="RNA")


def run_signature(obj, features, verbose=False):
    """Local perturbation signature over the first ``NDIMS`` PCs (``CalcPerturbSig``)."""
    calc_perturb_sig(
        obj, assay="RNA", features=features, labels="gene", nt_class="NT",
        split_by=SPLIT_BY, num_neighbors=NUM_NEIGHBORS, reduction="pca",
        ndims=NDIMS, new_assay="PRTB",
    )
    return obj


def run_classification(obj, verbose=False):
    """Split each guide's cells into KO / NP (``RunMixscape``)."""
    run_mixscape(
        obj, assay="PRTB", labels="gene", nt_class="NT", de_assay="RNA",
        min_de_genes=5, iter_num=10, prtb_type="KO", verbose=verbose,
    )
    return obj


def run_lda(obj, npcs=10, verbose=False):
    """Supervised map separating the guide populations (``MixscapeLDA``)."""
    mixscape_lda(
        obj, labels="gene", nt_class="NT", assay="PRTB", de_assay="RNA",
        npcs=npcs, verbose=verbose,
    )
    return obj


def summarize(obj, verbose=True) -> dict:
    """Print the Mixscape report and return the headline numbers as a dict."""
    meta = obj.meta_data
    gclass = pd.Series(np.asarray(meta["mixscape_class.global"]), dtype=object)
    table = perturbation_table(meta)
    globals_ = gclass.value_counts().to_dict()

    out = {
        "n_cells": int(len(gclass)),
        "global": globals_,
        "n_genes_tested": int(len(table)),
        "responsive": responsive_genes(table, 0.5),
        "table": table,
    }
    if not verbose:
        return out

    section("1. Perturbation signature (CalcPerturbSig → PRTB assay)")
    prtb = obj.assays["PRTB"]
    print(f"    PRTB assay: {len(prtb.features())} genes x {len(prtb.cells())} cells")
    print(f"    (each cell minus its {NUM_NEIGHBORS} nearest NT controls, "
          f"within {SPLIT_BY}, over {NDIMS} PCs)")

    section("2. Mixscape classification (RunMixscape) — KO vs escaper vs control")
    for k in ("KO", "NP", "NT"):
        print(f"    {k:<4} {globals_.get(k, 0):>7}")
    print("\n  Per-guide knockout rate (how often the edit 'took'):")
    print("    " + table.to_string(index=False).replace("\n", "\n    "))
    print(f"\n  {len(out['responsive'])} guides with KO rate >= 50%: "
          f"{', '.join(out['responsive'])}")

    if "lda_assignments" in meta.columns and "lda" in obj.reductions:
        section("3. MixscapeLDA — guide-population separation (visualization)")
        lda = obj.reductions["lda"]
        used = list(lda.misc.get("genes_used", []))
        out["lda_axes"] = int(lda.cell_embeddings.shape[1])
        out["lda_genes_used"] = used
        print(f"    LDA reduction: {lda.cell_embeddings.shape[1]} discriminant axes "
              f"from {len(used)} guides with a detectable signature")
        print(f"    ({', '.join(used)})")
        print("    Escaper (NP) cells fall onto the NT cloud by design, so this map is")
        print("    read visually (see figures), not as a per-cell accuracy.")
    return out


def report_concordance(obj, r_calls_path=None, verbose=True) -> dict | None:
    """Compare Python vs R per-cell calls, if ``thp1_mixscape_verify.R`` has run.

    Reads ``figures_mixscape/r_calls.csv``, aligns it to the object's cells, and
    reports the agreement fraction for the global class (KO/NP/NT), the full
    ``mixscape_class`` (``<gene> KO`` / ``<gene> NP`` / ``NT``), and the LDA guide
    assignment. Returns ``None`` (and prints a hint) when the R calls are absent.
    """
    path = Path(r_calls_path) if r_calls_path else FIGURES / "r_calls.csv"
    if not path.exists():
        if verbose:
            print(f"\n  [concordance] {path.name} not found — run "
                  "`Rscript tutorials/thp1_mixscape_verify.R` first.")
        return None

    r = pd.read_csv(path).set_index("cell")
    cells = obj.cell_names()
    r = r.reindex(cells)
    meta = obj.meta_data
    pairs = {
        "R_mixscape_global": "mixscape_class.global",
        "R_mixscape_class": "mixscape_class",
    }
    if "R_lda" in r.columns and "lda_assignments" in meta.columns:
        pairs["R_lda"] = "lda_assignments"

    agree = {}
    for rcol, pycol in pairs.items():
        py = pd.Series(np.asarray(meta[pycol]), index=cells)
        agree[rcol] = call_concordance(py.values, r[rcol].values)

    if verbose:
        section("4. R-vs-Python concordance (shared inputs & variable-feature basis)")
        print(f"    Mixscape global class (KO/NP/NT)     : {agree['R_mixscape_global']:.4f}")
        print(f"    Mixscape full class (<gene> KO/NP)   : {agree['R_mixscape_class']:.4f}")
        if "R_lda" in agree:
            print(f"    LDA guide assignment                 : {agree['R_lda']:.4f}")
        print("\n  Mixscape global — Python (rows) x R (cols):")
        py = pd.Series(np.asarray(meta["mixscape_class.global"]), index=cells)
        print(pd.crosstab(py.values, r["R_mixscape_global"].values).to_string())
    return agree


def run_full(data_dir=None, verbose=True, do_lda=True):
    """Load, build the signature, classify, (optionally) LDA, report, compare to R."""
    t0 = time.time()
    if verbose:
        section("Loading GSE153056 THP-1 ECCITE-seq screen")
    obj = load_screen_object(data_dir=data_dir)
    genes = [g for g in pd.unique(obj.meta_data["gene"]) if g != "NT"]
    if verbose:
        print(f"  {len(obj.assays['RNA'].features())} genes x "
              f"{len(obj.cell_names())} cells | {len(genes)} target guides + NT")

    hvg = prep_reduction(obj, verbose=verbose)
    FIGURES.mkdir(exist_ok=True)
    (FIGURES / "hvg_features.txt").write_text("\n".join(hvg) + "\n")
    if verbose:
        print(f"  variable features: {len(hvg)} (written to "
              f"figures_mixscape/hvg_features.txt for the R reference)")

    run_signature(obj, hvg, verbose=verbose)
    run_classification(obj, verbose=verbose)
    if do_lda:
        run_lda(obj, verbose=verbose)

    summary = summarize(obj, verbose=verbose)
    report_concordance(obj, verbose=verbose)

    if verbose:
        section("Summary")
        print(f"  Total runtime: {time.time() - t0:.1f}s")
    return obj, summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="THP-1 ECCITE-seq Mixscape tutorial")
    parser.add_argument("--data-dir", default=None)
    parser.add_argument("--no-lda", action="store_true", help="skip the (slow) LDA step")
    args = parser.parse_args()
    run_full(data_dir=args.data_dir, do_lda=not args.no_lda)
