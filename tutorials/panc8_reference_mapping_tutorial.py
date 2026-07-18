"""Reference mapping tutorial — annotate a query by label transfer with Shanuz on panc8.

A Python port of Seurat's *Mapping and annotating query datasets* vignette
(https://satijalab.org/seurat/articles/integration_mapping) on the pancreatic
islet dataset ``panc8``: ~14,900 human islet cells profiled across five
technologies (CEL-seq, CEL-seq2, SMART-seq2, Fluidigm C1, inDrop). Each
technology is its own batch with its own capture chemistry, so the same cell
type looks measurably different from one to the next — which is exactly what
makes cross-technology *annotation transfer* a real test.

The task reference mapping solves: **you have an annotated atlas you trust (the
reference) and a new, unlabelled dataset (the query); borrow the atlas's labels
to annotate the query — without ever moving or re-clustering the reference.**
Unlike integration (``ifnb_integration_tutorial.py``), which corrects several
datasets *onto each other*, mapping is deliberately asymmetric: the reference is
fixed, the query is projected into the reference's own PCA, and labels are
carried across mutual-nearest-neighbour anchors.

Reference vs query
------------------
Seurat's vignette builds a multi-technology *integrated* reference. This
tutorial deliberately uses a **single-technology reference** — CEL-seq2 (2,285
cells, carrying all 13 annotated cell types) — and maps the **SMART-seq2** query
(2,394 cells, also all 13 types) onto it. Two reasons: it isolates the label-
transfer machinery (:func:`shanuz.find_transfer_anchors` /
:func:`shanuz.transfer_data`) from the integration machinery the ifnb tutorial
already checks, and it is still a genuine cross-chemistry transfer (tag-based
CEL-seq2 → full-length SMART-seq2). Because the query's own ``celltype`` labels
ship with the data, the transfer can be scored directly for **accuracy** against
ground truth, not merely for agreement with R.

It demonstrates Shanuz's reference-mapping API side by side with Seurat's:
  * :func:`shanuz.find_transfer_anchors` — Seurat's ``FindTransferAnchors``:
    project the query into the reference's PCA (``reduction="pcaproject"``) and
    find scored mutual-nearest-neighbour anchors,
  * :func:`shanuz.transfer_data` — Seurat's ``TransferData``: turn the anchors
    into a per-query-cell weighted vote over the reference labels
    (``predicted.id`` + ``prediction.score.*``),
  * :func:`shanuz.map_query` / :func:`shanuz.project_umap` — Seurat's
    ``MapQuery`` / ``ProjectUMAP``: place the query in the reference's own UMAP
    so the new cells land on the atlas you already know how to read (figures).

Why this tutorial exists
------------------------
The whole reference-mapping stack (``find_transfer_anchors`` / ``transfer_data``
/ ``map_query`` / ``project_umap``) landed with only synthetic ``default_rng``
fixtures for tests — two cell types, a planted batch block, balanced sizes. This
is the first time it meets a real atlas with thirteen cell types (several of
them rare), genuine cross-technology batch structure, and a Seurat reference to
match. The comparison target is **not** byte-identical coordinates — the query's
projected embedding depends on irlba-vs-sklearn PCA sign and umap-learn-vs-uwot
transform — but the *labels*, which are a robust weighted argmax, are directly
comparable per cell: does shanuz assign each query cell the same ``predicted.id``
as Seurat, and does each tool recover the query's true cell type.

Note on the data and the R comparison
-------------------------------------
``panc8`` is a curated SeuratData object with no clean raw source, so both tools
read the *same* counts exported once from R by ``tutorials/export_seuratdata.R``
(a 10x-style matrix folder), guaranteeing byte-identical input and cell order.
To keep the projection on one shared gene basis, the Python run writes the
reference variable features it selected to ``figures_refmap/hvg_features.txt``
and the R script reads them back — so the only divergences left are the
genuinely method-level ones (PCA numerics, the anchor/weight kernels, kNN ties).

Usage
-----
    Rscript tutorials/export_seuratdata.R panc8            # one-time, writes the counts
    python  tutorials/panc8_reference_mapping_tutorial.py  # writes the shared HVGs

Then, for the side-by-side numbers and figures:

    Rscript tutorials/panc8_reference_mapping_verify.R
    python  tutorials/generate_refmap_plots.py

References
----------
Baron M, Veres A, Wolock SL, et al. (2016) **A single-cell transcriptomic map of
the human and mouse pancreas reveals inter- and intra-cell population
structure.** Cell Systems 3, 346-360. https://doi.org/10.1016/j.cels.2016.08.011

Stuart T, Butler A, Hoffman P, et al. (2019) **Comprehensive integration of
single-cell data.** Cell 177, 1888-1902. https://doi.org/10.1016/j.cell.2019.05.031
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

from shanuz.datasets import panc8
from shanuz.shanuz import create_shanuz_object
from shanuz.preprocessing import normalize_data, find_variable_features, scale_data
from shanuz.reduction import run_pca
from shanuz.umap import run_umap
from shanuz.transfer import find_transfer_anchors, transfer_data
from shanuz.mapping import project_umap

FIGURES = Path(__file__).parent / "figures_refmap"

# The two labels the tutorial turns on: the technology that defines the
# reference/query split (the batch) and the annotation being transferred.
TECH = "tech"
CELLTYPE = "celltype"

# A single-technology reference and a single-technology query (see the module
# docstring): CEL-seq2 → SMART-seq2, both carrying all 13 annotated cell types.
REFERENCE_TECH = "celseq2"
QUERY_TECH = "smartseq2"

# Params mirror Seurat's FindTransferAnchors / TransferData defaults: 2000
# variable features, project on 30 reference PCs, k.anchor 5, k.filter 200,
# k.weight 50.
N_HVG = 2000
N_PCS = 30
K_ANCHOR = 5
K_FILTER = 200
K_WEIGHT = 50


def section(title: str) -> None:
    print(f"\n{'=' * 70}\n{title}\n{'=' * 70}")


# ---------------------------------------------------------------------------
# Pure helpers (unit-tested without network in tests/test_refmap_tutorial.py)
# ---------------------------------------------------------------------------

def transfer_accuracy(predicted, truth) -> float:
    """Fraction of query cells whose transferred label matches the truth.

    The headline reference-mapping metric: ``mean(predicted.id == celltype)``.
    Order-sensitive — ``predicted`` and ``truth`` must be aligned to the same
    query cells.
    """
    predicted = np.asarray(predicted)
    truth = np.asarray(truth)
    if predicted.shape != truth.shape or predicted.size == 0:
        raise ValueError("predicted and truth must be the same non-zero length")
    return float(np.mean(predicted == truth))


def label_concordance(pred_a, pred_b) -> float:
    """Fraction of query cells the two tools assign the *same* label.

    Distinct from :func:`transfer_accuracy`: this compares two *predictors* to
    each other (shanuz vs Seurat's ``predicted.id``), not a predictor to ground
    truth. 1.0 means the two tools annotate every cell identically.
    """
    pred_a = np.asarray(pred_a)
    pred_b = np.asarray(pred_b)
    if pred_a.shape != pred_b.shape or pred_a.size == 0:
        raise ValueError("both label vectors must be the same non-zero length")
    return float(np.mean(pred_a == pred_b))


def per_class_recall(predicted, truth) -> pd.DataFrame:
    """Per-cell-type recall of the transfer: how often each true type is recovered.

    For every class present in ``truth``, the fraction of its cells that got the
    correct ``predicted.id`` (recall), with the class's support (cell count).
    Sorted by support descending, so the abundant, statistically meaningful cell
    types lead and the rare ones (which a small single-tech reference annotates
    noisily) sit at the bottom where they belong.
    """
    predicted = np.asarray(predicted)
    truth = np.asarray(truth)
    rows = []
    for c in sorted({str(x) for x in truth.tolist()}):
        mask = truth.astype(str) == c
        n = int(mask.sum())
        recall = float((predicted.astype(str)[mask] == c).mean()) if n else float("nan")
        rows.append({"celltype": c, "support": n, "recall": recall})
    df = pd.DataFrame(rows, columns=["celltype", "support", "recall"])
    return df.sort_values("support", ascending=False).reset_index(drop=True)


def macro_recall(predicted, truth) -> float:
    """Unweighted mean of the per-class recalls (every cell type counts equally).

    Complements :func:`transfer_accuracy` (which is cell-weighted, so the
    abundant alpha/beta cells dominate): the macro average exposes whether the
    rare types are being recovered at all.
    """
    rec = per_class_recall(predicted, truth)["recall"].to_numpy()
    rec = rec[np.isfinite(rec)]
    return float(rec.mean()) if rec.size else float("nan")


def build_scoreboard(rows: list[dict]) -> pd.DataFrame:
    """Assemble per-tool transfer metrics into a tidy, ordered table."""
    cols = ["tool", "n_query", "n_anchors", "accuracy", "macro_recall", "mean_score"]
    df = pd.DataFrame(rows)
    return df[[c for c in cols if c in df.columns]]


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

def load_panc8_split(reference_tech=REFERENCE_TECH, query_tech=QUERY_TECH,
                     data_dir=None, min_cells=3):
    """Load ``panc8`` and split it into a reference and a query by technology.

    Reads the counts exported by ``tutorials/export_seuratdata.R panc8`` (raises
    a helpful error if that one-time export has not run), builds one object over
    the full gene universe (so reference and query share an identical feature
    set), then subsets it to the two technologies. Both carry the :data:`TECH`
    and :data:`CELLTYPE` metadata; the query's ``celltype`` is held back as
    ground truth and never fed to the transfer.
    """
    counts, genes, cells, meta = panc8(data_dir=data_dir)
    keep = meta[[c for c in (TECH, CELLTYPE) if c in meta.columns]].copy()
    full = create_shanuz_object(
        counts=counts, assay="RNA", min_cells=min_cells, min_features=0,
        project="panc8", feature_names=list(genes), cell_names=list(cells),
        meta_data=keep,
    )

    cell_arr = np.asarray(full.cell_names())
    tech_vals = full.meta_data.loc[cell_arr, TECH].to_numpy().astype(str)
    ref_cells = cell_arr[tech_vals == reference_tech].tolist()
    query_cells = cell_arr[tech_vals == query_tech].tolist()
    if not ref_cells:
        raise ValueError(f"no cells for reference tech {reference_tech!r}")
    if not query_cells:
        raise ValueError(f"no cells for query tech {query_tech!r}")

    return full.subset(cells=ref_cells), full.subset(cells=query_cells)


def prep_reference(reference, n_hvg=N_HVG, n_pcs=N_PCS, do_umap=False, seed=42):
    """Prep the reference atlas: normalize → variable features → scale → PCA (→ UMAP).

    The reference's standard ``NormalizeData %>% FindVariableFeatures %>%
    ScaleData %>% RunPCA``. When ``do_umap`` is set, also fit a returnable UMAP
    (``run_umap`` stashes the model) so :func:`shanuz.project_umap` can place the
    query in the reference's embedding for the figures. Returns the selected
    variable-feature list — the shared basis written for the R run.
    """
    normalize_data(reference, assay="RNA")
    find_variable_features(reference, assay="RNA", nfeatures=n_hvg)
    scale_data(reference, assay="RNA")
    run_pca(reference, assay="RNA", n_pcs=n_pcs)
    if do_umap:
        dims = range(min(n_pcs, reference.reductions["pca"].cell_embeddings.shape[1]))
        run_umap(reference, reduction="pca", dims=dims, seed=seed)
    return variable_feature_list(reference, assay="RNA")


def prep_query(query, features):
    """Prep the query on the reference's variable-feature basis.

    ``NormalizeData`` then ``ScaleData(features = <reference HVGs>)`` — the query
    is scaled on exactly the reference's anchor features so the projection into
    the reference PCA and the anchor search line up gene-for-gene. The query is
    never given its own variable features or PCA; it only ever lives in the
    reference's space.
    """
    normalize_data(query, assay="RNA")
    scale_data(query, assay="RNA", features=list(features))


def run_mapping(reference, query, features, do_umap=False, k_anchor=K_ANCHOR,
                k_filter=K_FILTER, k_weight=K_WEIGHT, n_pcs=N_PCS, seed=42):
    """Find transfer anchors, transfer the cell-type labels, (optionally) project.

    Mirrors ``FindTransferAnchors(reduction="pcaproject") %>% TransferData``:
    project the query into the reference PCA on the shared ``features``, score
    the mutual-nearest-neighbour anchors, and carry the reference's
    :data:`CELLTYPE` across them to a per-query-cell ``predicted.id`` +
    ``prediction.score.*``. When ``do_umap`` is set, also
    :func:`shanuz.project_umap` the query into the reference UMAP (Seurat's
    ``MapQuery``/``ProjectUMAP``) for the atlas-placement figure.

    Returns ``(anchors, predictions)``.
    """
    anchors = find_transfer_anchors(
        reference, query, anchor_features=list(features), reduction="pcaproject",
        dims=n_pcs, k_anchor=k_anchor, k_filter=k_filter, seed=seed,
    )
    predictions = transfer_data(anchors, refdata=CELLTYPE, k_weight=k_weight)
    if do_umap:
        project_umap(query, reference)
    return anchors, predictions


def _score_columns(predictions):
    """The ``prediction.score.<class>`` columns (excluding ``.max``)."""
    return [c for c in predictions.columns
            if c.startswith("prediction.score.") and c != "prediction.score.max"]


def summarize(reference, query, anchors, predictions, verbose=True) -> dict:
    """Score the Python transfer against the query's ground-truth cell types."""
    truth = np.asarray(query.meta_data[CELLTYPE]).astype(str)
    pred = predictions["predicted.id"].to_numpy().astype(str)
    acc = transfer_accuracy(pred, truth)
    per_class = per_class_recall(pred, truth)

    board = build_scoreboard([{
        "tool": "shanuz (pcaproject)",
        "n_query": int(len(query.cell_names())),
        "n_anchors": int(len(anchors.anchors)),
        "accuracy": round(acc, 4),
        "macro_recall": round(macro_recall(pred, truth), 4),
        "mean_score": round(float(predictions["prediction.score.max"].mean()), 4),
    }])

    out = {
        "n_ref": int(len(reference.cell_names())),
        "n_query": int(len(query.cell_names())),
        "n_celltypes": int(query.meta_data[CELLTYPE].nunique()),
        "n_anchors": int(len(anchors.anchors)),
        "accuracy": acc,
        "scoreboard": board,
        "per_class": per_class,
    }
    if not verbose:
        return out

    section(f"Label transfer {REFERENCE_TECH} → {QUERY_TECH} "
            f"(predicted.id vs the query's true cell type)")
    print("    " + board.to_string(index=False).replace("\n", "\n    "))
    print(f"\n  {int(acc * out['n_query'])}/{out['n_query']} query cells "
          f"annotated correctly ({acc:.1%}).")
    print("\n  Per cell type (recall = fraction of that true type recovered):")
    shown = per_class.round({"recall": 4})
    print("    " + shown.to_string(index=False).replace("\n", "\n    "))
    return out


def report_concordance(query, predictions, r_calls_path=None, verbose=True) -> dict | None:
    """Compare Python vs R label transfer, if ``panc8_reference_mapping_verify.R`` has run.

    Reads ``figures_refmap/r_calls.csv`` (Seurat's ``predicted.id`` +
    ``prediction.score.max`` per query cell) and reports each tool's accuracy
    against the ground-truth ``celltype`` and — the direct comparison — the
    fraction of query cells the two tools annotate *identically*. Returns
    ``None`` (with a hint) when the R calls are absent.
    """
    path = Path(r_calls_path) if r_calls_path else FIGURES / "r_calls.csv"
    if not path.exists():
        if verbose:
            print(f"\n  [concordance] {path.name} not found — run "
                  "`Rscript tutorials/panc8_reference_mapping_verify.R` first.")
        return None

    r = pd.read_csv(path).set_index("cell")
    cells = query.cell_names()
    r = r.reindex(cells)
    truth = np.asarray(query.meta_data[CELLTYPE]).astype(str)
    py_pred = predictions["predicted.id"].reindex(cells).to_numpy().astype(str)
    r_pred = r["R_predicted"].to_numpy().astype(str)

    out = {
        "py_accuracy": transfer_accuracy(py_pred, truth),
        "r_accuracy": transfer_accuracy(r_pred, truth),
        "concordance": label_concordance(py_pred, r_pred),
        "py_mean_score": float(predictions["prediction.score.max"].reindex(cells).mean()),
        "r_mean_score": float(np.nanmean(r["R_score_max"].to_numpy())),
    }

    if verbose:
        section("R-vs-Python concordance (shared counts & variable-feature basis)")
        tbl = pd.DataFrame([
            {"tool": "shanuz", "accuracy": round(out["py_accuracy"], 4),
             "mean_score": round(out["py_mean_score"], 4)},
            {"tool": "Seurat R", "accuracy": round(out["r_accuracy"], 4),
             "mean_score": round(out["r_mean_score"], 4)},
        ])
        print("    " + tbl.to_string(index=False).replace("\n", "\n    "))
        print(f"\n  Label concordance (same predicted.id per cell): "
              f"{out['concordance']:.4f} "
              f"({int(out['concordance'] * len(cells))}/{len(cells)} query cells).")
        print("  accuracy  : each tool's predicted.id vs the known cell types.")
        print("  mean_score: mean prediction.score.max (transfer confidence).")
    return out


def run_full(data_dir=None, verbose=True, do_umap=False,
             reference_tech=REFERENCE_TECH, query_tech=QUERY_TECH):
    """Load, split, prep reference + query, transfer labels, score, compare to R."""
    t0 = time.time()
    if verbose:
        section(f"Loading panc8 — reference {reference_tech} → query {query_tech}")
    reference, query = load_panc8_split(
        reference_tech=reference_tech, query_tech=query_tech, data_dir=data_dir)
    if verbose:
        print(f"  reference {reference_tech}: {len(reference.cell_names())} cells | "
              f"query {query_tech}: {len(query.cell_names())} cells | "
              f"{query.meta_data[CELLTYPE].nunique()} cell types")

    hvg = prep_reference(reference, do_umap=do_umap)
    prep_query(query, hvg)
    FIGURES.mkdir(exist_ok=True)
    (FIGURES / "hvg_features.txt").write_text("\n".join(hvg) + "\n")
    if verbose:
        print(f"  variable features: {len(hvg)} (written to "
              f"figures_refmap/hvg_features.txt for the R reference)")

    t1 = time.time()
    anchors, predictions = run_mapping(reference, query, hvg, do_umap=do_umap)
    if verbose:
        print(f"  {len(anchors.anchors)} transfer anchors + labels "
              f"transferred ({time.time() - t1:.1f}s)")

    summary = summarize(reference, query, anchors, predictions, verbose=verbose)
    report_concordance(query, predictions, verbose=verbose)

    if verbose:
        section("Summary")
        print(f"  Total runtime: {time.time() - t0:.1f}s")
    return reference, query, anchors, predictions, summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="panc8 reference-mapping (label transfer) tutorial")
    parser.add_argument("--data-dir", default=None)
    parser.add_argument("--reference-tech", default=REFERENCE_TECH)
    parser.add_argument("--query-tech", default=QUERY_TECH)
    parser.add_argument("--umap", action="store_true",
                        help="also project the query into the reference UMAP (figures)")
    args = parser.parse_args()
    run_full(data_dir=args.data_dir, do_umap=args.umap,
             reference_tech=args.reference_tech, query_tech=args.query_tech)
