"""Cell-hashing demultiplexing tutorial — HTODemux + MULTIseqDemux with Shanuz.

A Python port of Seurat's hashing vignette
(https://satijalab.org/seurat/articles/hashing_vignette) using the original
Cell-Hashing dataset (Stoeckius et al. 2018, GSE108313): eight samples, each
tagged with a distinct antibody hashtag (BatchA–H), pooled and sequenced on one
lane. Demultiplexing recovers, per droplet, which sample the cell came from —
and which droplets caught two cells (doublets) or none (negatives).

It demonstrates Shanuz's two demultiplexers side by side with Seurat's:
  * build the object, attach the hashtag counts as an ``"HTO"`` assay, and
    CLR-normalise them (margin 1 — per hashtag across cells, Seurat's default),
  * :func:`shanuz.hto_demux` — Seurat's ``HTODemux``: a negative-binomial cutoff
    per hashtag, learned from its background cluster,
  * :func:`shanuz.multiseq_demux` — Seurat's ``MULTIseqDemux``: a kernel-density
    cutoff per barcode instead, and how the two methods differ at the margins,
  * an independent ground-truth check that only this dataset offers: the RNA is
    aligned to a *combined human+mouse* genome, so a droplet with both human and
    mouse transcripts is a species-mixed multiplet regardless of what the
    hashtags say.

Why this tutorial exists
------------------------
Every hashing feature landed after PR #10 and had only been checked against
synthetic fixtures. This is the first time ``hto_demux`` / ``multiseq_demux``
meet real data with a Seurat reference — a direct test of the CLR margin fix
(#32) and the ``clara`` k-medoids default (#34). The comparison target is
**R-vs-Python call concordance on byte-identical inputs**, reported by
:func:`report_concordance` once ``pbmc_hashing_verify.R`` has written its calls.

Note on the data
----------------
These are the raw GEO matrices (unfiltered barcodes), not the pre-filtered
``pbmc_umi_mtx.rds`` the vignette downloads from Dropbox (an R binary with no
clean cross-language form). So the singlet/doublet/negative totals here will
*not* match the vignette's headline numbers — the raw barcode list is dominated
by empty-ish droplets, hence a high Negative rate. That is expected and does not
affect the R-vs-Python comparison, which is the point.

Usage
-----
    python tutorials/pbmc_hashing_tutorial.py [--data-dir PATH]

The dataset (~34 MB) downloads automatically to ~/.shanuz_data/pbmc_hashing.
Then, for the side-by-side numbers and figures:

    Rscript tutorials/pbmc_hashing_verify.R
    python  tutorials/generate_hashing_plots.py

References
----------
Stoeckius M, Zheng S, Houck-Loomis B, et al. (2018)
**Cell Hashing with barcoded antibodies enables multiplexing and doublet
detection for single cell genomics.** Genome Biology 19, 224.
https://doi.org/10.1186/s13059-018-1603-1
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

from shanuz.datasets import pbmc_hashing
from shanuz.shanuz import create_shanuz_object
from shanuz.assay5 import create_assay5_object
from shanuz.preprocessing import normalize_data
from shanuz.hto import hto_demux
from shanuz.multiseq import multiseq_demux

FIGURES = Path(__file__).parent / "figures_hashing"


def section(title: str) -> None:
    print(f"\n{'=' * 70}\n{title}\n{'=' * 70}")


# ---------------------------------------------------------------------------
# Pure helpers (unit-tested without network in tests/test_hashing_tutorial.py)
# ---------------------------------------------------------------------------

def short_tags(names) -> list[str]:
    """``BatchA-AGGACCATCCAA`` → ``BatchA`` — the oligo sequence is display noise.

    Applied identically on the R side so hashtag names, ``hash.ID`` values and
    plot labels line up between the two tools.
    """
    return [str(n).split("-", 1)[0] for n in names]


def species_mask(genes) -> np.ndarray:
    """Boolean mask, True where a gene symbol is human.

    The combined reference names human genes in all upper case (``MALAT1``,
    ``CD3E``) and mouse genes in title case (``Xkr4``, ``mt-Nd1``, ``…Rik``), so
    "no lower-case letters" cleanly separates the two species.
    """
    return np.array([str(g) == str(g).upper() for g in genes])


def species_labels(counts, genes, human_hi: float = 0.9, human_lo: float = 0.1):
    """Per-cell ``(human_fraction, label)`` from a combined human+mouse count matrix.

    ``label`` is ``"human"`` / ``"mouse"`` above/below the fraction cutoffs and
    ``"mixed"`` in between — a species-mixed droplet, i.e. an independent doublet
    call that owes nothing to the hashtags. ``human_fraction`` is ``nan`` for
    cells with no human+mouse counts at all.
    """
    import scipy.sparse as sp

    mask = species_mask(genes)
    csr = counts.tocsr() if sp.issparse(counts) else np.asarray(counts)
    human = np.asarray(csr[mask].sum(axis=0)).ravel().astype(float)
    mouse = np.asarray(csr[~mask].sum(axis=0)).ravel().astype(float)
    total = human + mouse
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(total > 0, human / total, np.nan)
    label = np.where(
        np.isnan(frac), "empty",
        np.where(frac > human_hi, "human",
                 np.where(frac < human_lo, "mouse", "mixed")),
    )
    return frac, label


def call_concordance(a, b) -> float:
    """Fraction of positions where two per-cell call vectors agree."""
    a = np.asarray(a, dtype=object)
    b = np.asarray(b, dtype=object)
    if a.shape != b.shape or a.size == 0:
        raise ValueError("call vectors must be the same non-zero length")
    return float(np.mean(a == b))


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def load_hashing_object(data_dir=None, min_cells=3, short_names=True):
    """Load GSE108313, build the RNA object, attach the CLR-normalised HTO assay.

    Returns ``(obj, tags)`` where ``tags`` is the (optionally shortened) hashtag
    name list. The object carries per-cell ``human_frac`` / ``species`` metadata
    (the cross-species ground truth) and ``nCount_HTO``.
    """
    rna, genes, hto, hto_names, cells = pbmc_hashing(data_dir=data_dir)

    obj = create_shanuz_object(
        counts=rna, assay="RNA", min_cells=min_cells, min_features=0,
        project="hashing", feature_names=genes, cell_names=cells,
    )
    kept = obj.cell_names()
    cpos = {c: i for i, c in enumerate(cells)}
    idx = [cpos[c] for c in kept]

    tags = short_tags(hto_names) if short_names else list(hto_names)
    hto_kept = hto[:, idx].tocsc()
    obj.assays["HTO"] = create_assay5_object(
        counts=hto_kept, feature_names=tags, cell_names=kept, key="hto_",
    )
    # Margin 1 (per hashtag across cells) is Seurat's default and what the
    # hashing vignette normalises with; hto_demux/multiseq_demux then read this
    # `data` layer directly (normalize=False below).
    normalize_data(obj, assay="HTO", normalization_method="CLR", margin=1)

    obj.meta_data["nCount_HTO"] = np.asarray(hto_kept.sum(axis=0)).ravel()

    rna_kept = rna[:, idx].tocsr()
    mask = species_mask(genes)
    obj.meta_data["nCount_human"] = np.asarray(rna_kept[mask].sum(axis=0)).ravel()
    obj.meta_data["nCount_mouse"] = np.asarray(rna_kept[~mask].sum(axis=0)).ravel()
    frac, label = species_labels(rna_kept, genes)
    obj.meta_data["human_frac"] = frac
    obj.meta_data["species"] = label
    return obj, tags


def run_demux(obj, positive_quantile=0.99, verbose=False):
    """Run both demultiplexers on the object's ``HTO`` assay (in place)."""
    hto_demux(obj, assay="HTO", positive_quantile=positive_quantile,
              normalize=False, verbose=verbose)
    multiseq_demux(obj, assay="HTO", quantile=0.7, normalize=False, verbose=verbose)
    return obj


def summarize(obj, tags, verbose=True) -> dict:
    """Print the demux report and return the headline counts as a dict."""
    meta = obj.meta_data
    ht_glob = pd.Series(meta["HTO_classification.global"].values)
    hash_id = pd.Series(meta["hash.ID"].values)
    multi = pd.Series(meta["MULTI_ID"].values)
    ms_glob = multi.map(lambda s: s if s in ("Doublet", "Negative") else "Singlet")

    out = {
        "n_cells": int(len(ht_glob)),
        "hto_global": ht_glob.value_counts().to_dict(),
        "multi_global": ms_glob.value_counts().to_dict(),
    }
    if not verbose:
        return out

    section("1. HTODemux — negative-binomial threshold per hashtag")
    for k in ("Singlet", "Doublet", "Negative"):
        print(f"    {k:<9} {out['hto_global'].get(k, 0):>7}")
    print("\n  Singlets assigned per sample (hash.ID):")
    for tag in tags:
        print(f"    {tag:<8} {int((hash_id == tag).sum()):>7}")

    section("2. MULTIseqDemux — kernel-density threshold per barcode")
    for k in ("Singlet", "Doublet", "Negative"):
        print(f"    {k:<9} {out['multi_global'].get(k, 0):>7}")
    print("\n  HTODemux (rows) x MULTIseqDemux (cols):")
    print(pd.crosstab(ht_glob, ms_glob).to_string())

    section("3. Cross-species ground truth (combined human+mouse reference)")
    sp = pd.Series(meta["species"].values)
    print("  A species-'mixed' droplet is a multiplet no hashtag was consulted for.")
    print(pd.crosstab(ht_glob, sp).to_string())
    n_mixed = int((sp == "mixed").sum())
    dbl_mixed = int(((sp == "mixed") & (ht_glob == "Doublet")).sum())
    print(f"\n  {n_mixed} species-mixed cells; {dbl_mixed} of them HTODemux calls Doublet.")
    print("  (At raw-barcode depth the mixed pool is mostly low-count ambient "
          "droplets,\n   so this is a weak sanity check, not a clean doublet "
          "benchmark.)")
    return out


def report_concordance(obj, r_calls_path=None, verbose=True) -> dict | None:
    """Compare Python vs R per-cell calls, if ``pbmc_hashing_verify.R`` has run.

    Reads ``figures_hashing/r_calls.csv`` (or ``r_calls_path``), aligns it to the
    object's cells, and reports the agreement fraction for the global class, the
    ``hash.ID`` singlet assignment, and the ``MULTI_ID`` call. Returns ``None``
    (and prints a hint) when the R calls are not present yet.
    """
    path = Path(r_calls_path) if r_calls_path else FIGURES / "r_calls.csv"
    if not path.exists():
        if verbose:
            print(f"\n  [concordance] {path.name} not found — run "
                  "`Rscript tutorials/pbmc_hashing_verify.R` first.")
        return None

    r = pd.read_csv(path).set_index("cell")
    cells = obj.cell_names()
    r = r.reindex(cells)
    meta = obj.meta_data
    py = {
        "R_HTO_global": pd.Series(meta["HTO_classification.global"].values, index=cells),
        "R_hash_ID": pd.Series(meta["hash.ID"].values, index=cells),
        "R_MULTI_ID": pd.Series(meta["MULTI_ID"].values, index=cells),
    }
    agree = {col: call_concordance(py[col].values, r[col].values) for col in py}

    if verbose:
        section("4. R-vs-Python call concordance (byte-identical inputs)")
        print(f"    HTODemux global (Singlet/Doublet/Negative) : {agree['R_HTO_global']:.4f}")
        print(f"    HTODemux sample assignment (hash.ID)       : {agree['R_hash_ID']:.4f}")
        print(f"    MULTIseqDemux call (MULTI_ID)              : {agree['R_MULTI_ID']:.4f}")
        print("\n  HTODemux global — Python (rows) x R (cols):")
        print(pd.crosstab(py["R_HTO_global"].values, r["R_HTO_global"].values).to_string())
    return agree


def run_full(data_dir=None, verbose=True, positive_quantile=0.99):
    """Load, demultiplex, report, and (if available) compare against R."""
    t0 = time.time()
    if verbose:
        section("Loading GSE108313 cell-hashing data")
    obj, tags = load_hashing_object(data_dir=data_dir)
    if verbose:
        print(f"  {len(obj.assays['RNA']._all_feature_names)} genes x "
              f"{len(obj.cell_names())} cells | {len(tags)} hashtags: "
              f"{', '.join(tags)}")

    run_demux(obj, positive_quantile=positive_quantile, verbose=False)
    summary = summarize(obj, tags, verbose=verbose)
    report_concordance(obj, verbose=verbose)

    if verbose:
        section("Summary")
        print(f"  Total runtime: {time.time() - t0:.1f}s")
    return obj, tags, summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cell-hashing demultiplexing tutorial")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    run_full(data_dir=args.data_dir)
