"""The object model itself — Shanuz vs R Seurat on pbmc3k.

Every other tutorial in this series compares an *algorithm*: does ``run_pca``
land where ``RunPCA`` lands, does ``hto_demux`` call the same cells. This one
compares the **container** — the accessors, the layer machinery and the
bookkeeping that Seurat's own "command cheat sheet" vignette is made of:

  * ``Cells`` / ``Features`` / ``dim`` — who is in the object, in what order;
  * ``Layers`` / ``LayerData`` / ``split`` / ``JoinLayers`` — the v5 layered
    assay, the single biggest change from Seurat v4;
  * ``DefaultAssay`` / ``GetAssayData`` / ``Key`` / ``VariableFeatures``;
  * ``Embeddings`` / ``Loadings`` / ``Stdev`` on a ``DimReduc``;
  * ``Graphs`` and ``as_graph`` / ``as_neighbor``;
  * ``FetchData``, ``Idents`` / ``WhichCells`` / ``RenameIdents`` / ``subset``;
  * ``Command`` — the log of what has been run on the object;
  * ``Misc`` / ``Tool`` — the user-extensible slots.

Why this is the sharpest net in the series
------------------------------------------
Almost nothing here is stochastic. Cell orders, layer names, matrix dimensions,
key strings, non-zero counts and command names either match R exactly or they
are wrong — 89 of the 91 anchors are compared with no tolerance at all. The
other tutorials have to argue about whether a 2 % gap is RNG or a bug; this one
mostly does not have that conversation. The exceptions are the two anchors that
read a PCA (``pca_stdev_head`` and ``abs_PC_1_sum``), where shanuz's randomized
SVD and Seurat's ``irlba`` differ in the fourth decimal; they are named
explicitly in :data:`FLOAT_TOLERANCES` rather than covered by a blanket rule.

That matters because the object layer is where a defect is quietest. A broken
``run_pca`` shows up in a plot. A broken ``JoinLayers`` returns a matrix of
exactly the right shape, full of exactly the right numbers, in the wrong order.

Coverage gap this closes
------------------------
``Assay5``, ``StdAssay``, ``Graph``, ``ShanuzCommand`` and ``log_shanuz_command``
are public exports that appeared in **no** tutorial before this one, and
``join_layers`` / ``split_layers`` had **zero call sites and zero tests**
anywhere in the package — the defining feature of the v5 object model, and
nothing had ever run it. (The legacy v3 ``Assay`` and ``create_assay_object``
are still uncovered: this tutorial builds ``Assay5`` objects throughout, which
is what ``create_shanuz_object`` returns.)

What it found
-------------
**Eleven** fidelity defects across four areas, all fixed in the pull request
that brought this tutorial: five in the layered assay, three in ``FetchData``,
and three in the object's bookkeeping. Most were in code no test touched; the
rest were guarded by tests that could not fail.

**The layered assay — ``split`` / ``JoinLayers`` was not a round trip.** After a
split and a join, shanuz returned a layer named ``joined`` (Seurat restores the
original name), holding the right numbers **in the split's order** rather than
the assay's. The assay's own cell vector never moves during a split, so the
matrix came back silently transposed against the metadata that indexes it: ask
for cell ``c1``'s column and get ``c2``'s. Nothing raised, the shape was right,
and the checksums were right, because every value was still present. The
no-argument ``join_layers()`` — R's idiom, and the only call a real script makes
— additionally raised ``ValueError`` on any prepared assay, because it hstacked
``counts``, ``data`` and a variable-features-only ``scale.data`` together
regardless of their differing feature counts. And ``shanuz.generics.split_layers``
was declared but never registered for any type, so the documented generic raised
``NotImplementedError`` while the method it should have dispatched to worked.

**``FetchData`` returned objects instead of numbers.** ``np.asarray`` on a sparse
matrix yields a 0-d *object* array wrapping it, not its contents, so
``.flatten()`` broadcast one ``csc_matrix`` down all 2,700 rows: every cell's
"expression" was a copy of the whole matrix. This is the most-called accessor in
Seurat, on the default assay class. It had a test — ``test_fetch_data`` — which
asserted the column name and the row count, both of which that satisfies
perfectly. Two smaller gaps came with it: embedding columns could not be
addressed by their ``Key`` (``PC_1`` raised ``KeyError``; R names them that way
precisely so they can be), and an unqualified fetch read ``counts`` where R reads
``data``, returning raw integers where every vignette shows normalized values.

**The command log was inert.** ``log_shanuz_command`` was a public export with no
call sites, so ``obj.commands`` was always empty. Seurat logs five entries for
this pipeline, keyed ``NormalizeData.RNA`` through ``FindNeighbors.RNA.pca``.
Also: ``orig.ident`` — the first column of every Seurat object's metadata — was
never created, and ``add_meta_data`` rejected the plain vector that R's
``AddMetaData`` documents and every vignette passes it.

Two differences are left standing on purpose. ``FindNeighbors`` builds a
*symmetrized* kNN graph where Seurat's is directed (nnz 75,740 against exactly
2,700 × 20 = 54,000), and the SNN drops the self-edge Seurat keeps. Both are
``find_neighbors``, not the object model, and both change what clustering
consumes, so they belong to their own comparison rather than to this one. They
are recorded rather than quietly absorbed.

Conventions
-----------
The batch column that drives the layer split is assigned **deterministically**
(alternating cells), not with ``sample()`` as Seurat's vignette does, so that R
and Python split on identical membership and any difference downstream belongs
to the layer code rather than to two different random draws.

Identities come from marker-gene thresholds rather than from clustering for the
same reason: Louvain drifts by a cluster between the two tools, which would put
noise into ``WhichCells`` / ``subset`` / ``RenameIdents`` comparisons that are
otherwise exact.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz import generics as G  # noqa: E402
from shanuz.datasets import pbmc3k  # noqa: E402
from shanuz.graph import as_graph  # noqa: E402
from shanuz.neighbors import find_neighbors  # noqa: E402
from shanuz.preprocessing import (  # noqa: E402
    find_variable_features,
    normalize_data,
    scale_data,
)
from shanuz.reduction import run_pca  # noqa: E402
from shanuz.shanuz import create_shanuz_object  # noqa: E402

FIGURES = Path(__file__).parent / "figures_objects"

MIN_CELLS = 3
MIN_FEATURES = 200
N_HVG = 2000
N_PCS = 20
PROJECT = "pbmc3k_objects"

#: Marker gates for the deterministic identity assignment. Order matters — the
#: first gate a cell passes wins, so the mapping is a function of the data and
#: nothing else. These are the canonical PBMC lineage markers; the point is a
#: reproducible partition, not a defensible cell-type call.
IDENT_MARKERS = (("T", "CD3E"), ("B", "MS4A1"), ("Mono", "LYZ"))
IDENT_OTHER = "Other"


# ---------------------------------------------------------------------------
# Pure helpers — the comparison vocabulary
# ---------------------------------------------------------------------------

def digest(names) -> str:
    """A short, order-**sensitive** fingerprint of a sequence of names.

    Cell and feature vectors here run to thousands of entries, and the anchor
    file has to stay readable. A hash compares them in one field; order matters
    because half of what this tutorial checks is that R and Python hold their
    cells in the *same* order, not merely the same set.
    """
    joined = "\n".join(str(n) for n in names)
    return hashlib.md5(joined.encode()).hexdigest()[:12]


def r_safe_names(names) -> list[str]:
    """Feature names as R's ``Read10X`` would spell them.

    ``Read10X`` rewrites ``_`` to ``-`` in gene symbols — pbmc3k carries
    ``Y_RNA`` and a long tail of ``RP11-*_*`` — while shanuz's loader leaves them
    alone. That difference belongs to the two file readers, not to the object
    model this tutorial audits, so feature-name anchors are compared in R's
    spelling on both sides rather than being allowed to fail every field that
    mentions a gene.
    """
    return [str(n).replace("_", "-") for n in names]


def name_anchor(names) -> dict:
    """Fingerprint plus enough plain text to debug a mismatch by eye.

    A bare hash tells you *that* two lists differ and nothing else, which is
    useless at 2,700 cells. The head and tail localise the usual suspects — a
    different filter cut, a sort, a reversed order — without storing the lot.
    """
    names = [str(n) for n in names]
    return {
        "n": len(names),
        "digest": digest(names),
        "head": names[:3],
        "tail": names[-3:],
    }


def matrix_anchor(mat) -> dict:
    """Shape, sparsity and total of a layer, as exactly comparable scalars.

    ``nnz`` is the one that catches a silently densified layer: a stored
    explicit zero changes it while leaving both the shape and the sum intact.
    """
    dense = mat.toarray() if sp.issparse(mat) else np.asarray(mat)
    return {
        "shape": list(dense.shape),
        "nnz": int(np.count_nonzero(dense)),
        "sum": round(float(dense.sum()), 6),
    }


def layer_anchor(assay) -> dict:
    """Every layer of an assay: its name, dimensions, and the cells it holds."""
    out = {}
    for name in assay.layers_list():
        anchor = matrix_anchor(assay.layer_data(name))
        anchor["cells"] = name_anchor(assay.cells(name))
        out[name] = anchor
    return out


def split_join_roundtrip(assay, groups, layer="counts") -> dict:
    """Split an assay into per-group layers, join it back, and audit the result.

    Seurat's contract for this pair is that it is a **round trip**: after
    ``obj[["RNA"]] <- split(obj[["RNA"]], f = ...)`` and then ``JoinLayers``,
    you are back to the layer you started from — same name, same cells, in the
    same order. Anything else silently decouples the matrix from the metadata
    that indexes it, because the object's own cell vector never moved.

    Run on a **counts-only** assay in both tools, deliberately: on a prepared
    object R's ``split`` divides every layer at once while shanuz's takes one
    layer at a time, and comparing those two would be comparing the API shapes
    rather than the round trip. Stripped to one layer, both tools are being
    asked exactly the same question.

    Returns the observed layer names at each stage and the three facts that
    decide whether the join honoured the contract.
    """
    before_names = assay.layers_list()
    before = assay.layer_data(layer)
    before = before.toarray() if sp.issparse(before) else np.asarray(before)
    original_cells = list(assay.cells(layer))

    split = assay.split_layers(list(groups), layer=layer)
    split_names = split.layers_list()
    joined = split.join_layers()

    joined_names = joined.layers_list()
    target = joined_names[0] if joined_names else None
    after = joined.layer_data(target) if target else None
    after = (after.toarray() if sp.issparse(after) else np.asarray(after)) \
        if after is not None else np.empty((0, 0))
    joined_cells = list(joined.cells(target)) if target else []

    return {
        "layers_before": before_names,
        "layers_after_split": split_names,
        "layers_after_join": joined_names,
        # Does the rejoined layer carry the name it had before the split?
        "layer_name_restored": target == layer,
        # Are the cells back in their original order?
        "cell_order_restored": joined_cells == original_cells,
        # And does the matrix itself come back unchanged?
        "matrix_restored": bool(
            after.shape == before.shape and np.array_equal(after, before)
        ),
        "cells_after_join": name_anchor(joined_cells),
        # The object's own cell vector is untouched by either operation, so if
        # the layer's order moved, this is what it is now out of step with.
        "assay_cells": name_anchor(assay.cells()),
    }


def graph_neighbor_roundtrip(obj, graph_name="RNA_nn") -> dict:
    """``Graph`` ↔ ``Neighbor`` conversion — a within-tool invariant.

    Not compared against R, and kept out of the anchor tree for that reason:
    ``as.Neighbor`` reads a *directed* kNN graph, and shanuz's is symmetrized
    (one of the two divergences this tutorial reports but does not fix), so the
    two tools are not being asked the same question here.

    What is checked is the round trip's own consistency: a graph converted to a
    neighbour list and back must describe the same edges it started with.
    """
    graph = obj.graphs.get(graph_name)
    if graph is None:
        return {"available": False, "reason": f"no graph named {graph_name!r}"}
    neighbor = G.as_neighbor(graph)
    back = as_graph(neighbor, cell_names=graph.cells())
    original = graph.tocsr()
    restored = back.tocsr()
    degrees = np.diff(original.indptr)
    return {
        "available": True,
        "n_cells": original.shape[0],
        # Seurat's kNN graph gives every cell exactly `k.param` neighbours, so
        # these three numbers would all read 20. They do not, because the graph
        # is symmetrized — and `as_neighbor` pads its index matrix out to the
        # widest row, so the resulting Neighbor's second dimension is the max
        # degree rather than any meaningful k.
        "degree_min": int(degrees.min()),
        "degree_max": int(degrees.max()),
        "degree_mean": round(float(degrees.mean()), 2),
        "neighbor_width": int(neighbor.indices().shape[1]),
        # Every edge the neighbour list knows about must survive the trip back.
        "edges_preserved": bool(
            (restored.multiply(original) != 0).nnz == restored.count_nonzero()
        ),
    }


def join_all_layers(assay, groups, layer="counts") -> dict:
    """``JoinLayers`` on a *prepared* assay — the call every real script makes.

    After the standard pipeline an assay holds ``counts``, ``data`` and
    ``scale.data``, and ``scale.data`` has only the variable features. Seurat's
    ``JoinLayers(obj)`` handles that by rejoining each split stem separately and
    leaving everything else alone. Whether shanuz's no-argument call survives
    the same object is the question.
    """
    split = assay.split_layers(list(groups), layer=layer)
    try:
        joined = split.join_layers()
    except Exception as exc:  # noqa: BLE001 — the failure mode *is* the finding
        return {"error": f"{type(exc).__name__}: {exc}", "layers_after_join": None}
    return {"error": None, "layers_after_join": joined.layers_list()}


def assign_idents(obj, markers=IDENT_MARKERS, other=IDENT_OTHER) -> np.ndarray:
    """Label each cell by the first marker gate it passes.

    Deterministic by construction: a fixed gene order, a fixed ``> 0`` gate on
    normalized data, and a fixed fallback. The R side applies the identical rule
    to its own normalized matrix, so the two partitions can be compared cell by
    cell.

    Reads the layer directly rather than through :func:`shanuz.generics.fetch_data`
    — the obvious way to write this — because ``fetch_data`` is one of the things
    under test here, and a tutorial that assigns its identities with the function
    it is auditing cannot report on it. See :func:`fetch_anchor`.
    """
    assay = obj.get_assay("RNA")
    cells = obj.cell_names()
    labels = np.full(len(cells), other, dtype=object)
    # Reverse order so that earlier gates overwrite later ones — "first gate
    # wins" without a per-cell loop.
    for label, gene in reversed(markers):
        row = assay.layer_data(layer="data", features=[gene], cells=cells)
        row = row.toarray() if sp.issparse(row) else np.asarray(row)
        labels[row.ravel() > 0] = label
    return labels


def ident_anchor(idents) -> dict:
    """Levels and per-level counts, in a stable order."""
    series = pd.Series([str(i) for i in idents])
    counts = series.value_counts()
    return {
        "levels": sorted(counts.index.tolist()),
        "counts": {k: int(counts[k]) for k in sorted(counts.index)},
        "digest": digest(series),
    }


#: The only fields not required to match exactly, each with the reason it is
#: here. Kept deliberately short and explicit: the value of this tutorial is
#: that almost everything in it is an exact comparison, and a general tolerance
#: would quietly buy back all the sensitivity that makes it worth running.
#:
#: Both entries are the same cause — shanuz's PCA uses a randomized SVD where
#: Seurat's `RunPCA` uses `irlba`, and the two agree to about four decimals on
#: this data. Neither is a property of the object model.
#: Measured on pbmc3k: stdevs agree to ~4e-4 relative, and the summed absolute
#: PC_1 to 1.1e-3 (15804.19 against 15822.31) — accumulated over 2,700 cells.
#: PCA fidelity is T-dr's subject, not this tutorial's; what is checked exactly
#: here is that `fetch_data` returns *the object's own* embedding column, which
#: `fetch.pc1_matches_embeddings` pins with no tolerance at all.
FLOAT_TOLERANCES = {
    "reductions.pca_stdev_head": 1e-3,
    "fetch.abs_PC_1_sum": 5e-3,
}

#: Fields whose *contents* are compared but whose order is not. Seurat's own
#: order here is an artifact of how `JoinLayers` rebuilds the assay — it returns
#: `data, counts, scale.data`, which is not the order it started in either — so
#: requiring a match would be requiring shanuz to reproduce an accident.
UNORDERED_FIELDS = {"join_all_layers.layers_after_join"}


def compare_anchors(py: dict, r: dict, path: str = "") -> pd.DataFrame:
    """Walk two anchor trees in parallel and report every leaf as match/mismatch.

    The output is deliberately a flat table of booleans rather than a distance:
    for this tutorial there is no such thing as *nearly* the same layer name or
    *approximately* the same cell order, so a tolerance would only hide things.
    The handful of genuine exceptions are named in :data:`FLOAT_TOLERANCES` and
    :data:`UNORDERED_FIELDS` rather than applied across the board.
    """
    rows = []
    for key in sorted(set(py) | set(r)):
        here = f"{path}.{key}" if path else key
        if key not in py or key not in r:
            rows.append({"field": here, "python": py.get(key, "—"),
                         "r": r.get(key, "—"), "match": False})
            continue
        pv, rv = py[key], r[key]
        if isinstance(pv, dict) and isinstance(rv, dict):
            rows.extend(compare_anchors(pv, rv, here).to_dict("records"))
            continue
        rows.append({"field": here, "python": pv, "r": rv,
                     "match": _leaf_matches(here, pv, rv)})
    return pd.DataFrame(rows, columns=["field", "python", "r", "match"])


def _leaf_matches(field: str, pv, rv) -> bool:
    """Compare one leaf, honouring the two documented exception lists."""
    if field in UNORDERED_FIELDS:
        return sorted(map(str, pv or [])) == sorted(map(str, rv or []))
    rtol = FLOAT_TOLERANCES.get(field)
    if rtol is not None:
        return bool(np.allclose(np.atleast_1d(np.asarray(pv, dtype=float)),
                                np.atleast_1d(np.asarray(rv, dtype=float)),
                                rtol=rtol, atol=0))
    if isinstance(pv, float) or isinstance(rv, float):
        return bool(np.isclose(float(pv), float(rv), rtol=1e-6, atol=0))
    return bool(pv == rv)


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def load_pbmc3k_object(data_dir=None):
    """The PBMC 3k object at the plain gates, matching the other pbmc3k tutorials.

    No mitochondrial or upper-nFeature filtering: this tutorial is about the
    container, and every extra QC step is one more way the two tools could hold
    different cells before the accessors under test are ever called.
    """
    counts, genes, cells = pbmc3k(data_dir=data_dir)
    return create_shanuz_object(
        counts=counts, assay="RNA", min_cells=MIN_CELLS, min_features=MIN_FEATURES,
        project=PROJECT, feature_names=list(genes), cell_names=list(cells),
    )


def prep(obj, n_hvg=N_HVG, n_pcs=N_PCS):
    """Normalize → HVG → scale → PCA → neighbours, then publish the shared basis.

    The R script reads the HVG list and the cell list back so that both tools
    describe an object built on identical membership — the same trick every
    other tutorial in the series uses.
    """
    normalize_data(obj, normalization_method="LogNormalize", scale_factor=10000)
    find_variable_features(obj, selection_method="vst", nfeatures=n_hvg)
    hvg = list(obj.assays["RNA"].variable_features)
    scale_data(obj, features=hvg)
    run_pca(obj, n_pcs=n_pcs, features=hvg, reduction_name="pca")
    find_neighbors(obj, dims=list(range(10)))

    FIGURES.mkdir(exist_ok=True)
    (FIGURES / "hvg_features.txt").write_text("\n".join(hvg) + "\n")
    (FIGURES / "cells.txt").write_text("\n".join(obj.cell_names()) + "\n")
    return hvg


def batch_labels(n_cells: int) -> list[str]:
    """Alternating batch assignment — deterministic on both sides by design."""
    return ["batch1" if i % 2 == 0 else "batch2" for i in range(n_cells)]


def collect_anchors(obj, hvg) -> dict:
    """Every structural fact this tutorial compares, as one JSON-able tree."""
    assay = obj.get_assay("RNA")
    cells = obj.cell_names()

    anchors = {
        # --- shape and membership -------------------------------------
        "shape": {"n_cells": len(cells), "n_features": len(obj.feature_names())},
        "cells": name_anchor(cells),
        "features": name_anchor(r_safe_names(obj.feature_names())),

        # --- assays ---------------------------------------------------
        "assay": {
            "names": list(obj.assays.keys()),
            "default": obj.default_assay,
            "key": G.key(assay),
            "class": type(assay).__name__,
        },
        "layers": layer_anchor(assay),

        # --- metadata and identities ----------------------------------
        "meta": {
            "columns": sorted(obj.meta_data.columns.tolist()),
            "nCount_RNA_sum": round(float(obj.meta_data["nCount_RNA"].sum()), 6),
            "nFeature_RNA_sum": round(float(obj.meta_data["nFeature_RNA"].sum()), 6),
        },
        "idents": ident_anchor(obj.idents),

        # --- variable features ----------------------------------------
        "variable_features": name_anchor(r_safe_names(hvg)),

        # --- reductions -----------------------------------------------
        "reductions": {
            "names": sorted(obj.reductions.keys()),
            "pca_key": G.key(obj.reductions["pca"]),
            "pca_dims": list(G.embeddings(obj.reductions["pca"]).shape),
            "pca_n_loadings": list(G.loadings(obj.reductions["pca"]).shape),
            # Standard deviations are sign-independent, unlike the embeddings
            # themselves, so they are the part of a PCA worth comparing here.
            "pca_stdev_head": [round(float(v), 6)
                               for v in G.stdev(obj.reductions["pca"])[:5]],
        },

        # --- graphs ---------------------------------------------------
        # ``count_nonzero``, not ``nnz``: the two tools need not agree on which
        # zeros they bother to store, only on which edges exist.
        "graphs": {
            name: {"shape": list(g.shape),
                   "nnz": int(g.tocsr().count_nonzero())}
            for name, g in sorted(obj.graphs.items())
        },

        # --- the command log ------------------------------------------
        # Seurat appends an entry per pipeline call; this is the list of what
        # the object believes has been run on it.
        "commands": [c.key for c in obj.commands],
    }
    return anchors


def ident_operations(obj) -> dict:
    """``WhichCells`` / ``RenameIdents`` / ``subset``, exercised and measured."""
    t_cells = G.which_cells(obj, ident="T")
    renamed = obj.rename_idents({"T": "T_cell"})
    subset = obj.subset(idents="Mono")
    return {
        "which_cells_T": name_anchor(t_cells),
        "renamed_levels": sorted({str(i) for i in renamed.idents}),
        "subset_mono": {
            "n_cells": len(subset.cell_names()),
            "cells": name_anchor(subset.cell_names()),
            "n_features": len(subset.feature_names()),
        },
    }


def _column_sum(frame, column) -> float | str:
    """Sum a fetched column, or say why it cannot be summed.

    ``FetchData`` is under test, so this cannot assume the column holds numbers.
    A column that comes back as ``object`` — holding, say, sparse matrices
    rather than values — is reported as such and compared against R like any
    other anchor, instead of raising and taking the whole run down with it.
    """
    if column not in frame.columns:
        return "missing"
    series = frame[column]
    if series.dtype == object:
        return f"non-numeric: {type(series.iloc[0]).__name__}"
    return round(float(series.sum()), 6)


def fetch_anchor(obj) -> dict:
    """``FetchData`` across the three kinds of variable it has to reach.

    Metadata, a gene, and a reduction column — the three branches of R's
    ``FetchData``, which addresses every one of them by name in a single call.
    ``PC_1`` is named the way R names it: by the reduction's ``Key()``, which is
    what ``Key()`` is for.
    """
    wanted = ["nCount_RNA", "nFeature_RNA", "CD3E", "PC_1"]
    try:
        frame = G.fetch_data(obj, wanted)
    except KeyError as exc:
        return {"error": f"KeyError: {exc.args[0]}", "columns": [], "n_rows": 0}
    return {
        "columns": list(frame.columns),
        "n_rows": int(frame.shape[0]),
        "CD3E_sum": _column_sum(frame, "CD3E"),
        "nCount_RNA_sum": _column_sum(frame, "nCount_RNA"),
        # Exact, and within one tool: whatever the two PCAs do differently, a
        # fetched embedding column must *be* the object's embedding column. This
        # is the anchor that actually guards the fetch, and it needs no
        # tolerance because both sides are asking the same object.
        "pc1_matches_embeddings": bool(
            "PC_1" in frame.columns
            and frame["PC_1"].dtype != object
            and np.array_equal(
                frame["PC_1"].to_numpy(),
                G.embeddings(obj.reductions["pca"])[:, 0],
            )
        ),
        # PC_1's sign is arbitrary between the two tools, so its *absolute* sum
        # is the comparable quantity; the other two are directly comparable.
        "abs_PC_1_sum": (
            round(float(np.abs(frame["PC_1"]).sum()), 6)
            if "PC_1" in frame.columns and frame["PC_1"].dtype != object
            else _column_sum(frame, "PC_1")
        ),
    }


def run_full(data_dir=None, verbose=True):
    """The whole tutorial. Returns ``(obj, anchors)``."""
    obj = load_pbmc3k_object(data_dir=data_dir)
    hvg = prep(obj)
    if verbose:
        print(f"pbmc3k: {len(obj.cell_names())} cells x "
              f"{len(obj.feature_names())} features")

    obj.set_ident(obj.cell_names(), "unassigned")
    labels = assign_idents(obj)
    for label in sorted(set(labels)):
        obj.set_ident([c for c, lab in zip(obj.cell_names(), labels) if lab == label],
                      label)

    anchors = collect_anchors(obj, hvg)
    anchors["ident_ops"] = ident_operations(obj)
    anchors["fetch"] = fetch_anchor(obj)

    # The round trip runs on a pristine counts-only assay (see the function's
    # docstring for why), so reload rather than reuse the prepared one.
    batches = batch_labels(len(obj.cell_names()))
    fresh = load_pbmc3k_object(data_dir=data_dir).get_assay("RNA")
    anchors["split_join"] = split_join_roundtrip(fresh, batches, layer="counts")
    anchors["join_all_layers"] = join_all_layers(
        obj.get_assay("RNA"), batches, layer="counts",
    )

    FIGURES.mkdir(exist_ok=True)
    (FIGURES / "py_anchors.json").write_text(json.dumps(anchors, indent=2) + "\n")

    # Reported alongside, but not part of the compared tree — see the docstring.
    roundtrip = graph_neighbor_roundtrip(obj)

    if verbose:
        _print_report(anchors, report_concordance(anchors))
        if roundtrip["available"]:
            print(f"\n  Graph → Neighbor → Graph (no R counterpart): "
                  f"edges preserved={roundtrip['edges_preserved']}")
            print(f"    degree min={roundtrip['degree_min']} "
                  f"max={roundtrip['degree_max']} "
                  f"mean={roundtrip['degree_mean']}  "
                  f"(Seurat: exactly k=20 for every cell)")
    return obj, anchors


def report_concordance(anchors, figures=FIGURES):
    """Compare against the R reference, if ``pbmc3k_objects_verify.R`` has run."""
    path = figures / "r_anchors.json"
    if not path.exists():
        print("\n  R reference not found — run "
              "`Rscript tutorials/pbmc3k_objects_verify.R`")
        print("  for the side-by-side table. Python-side results above stand alone.")
        return None
    r = json.loads(path.read_text())
    return compare_anchors(anchors, r)


def _print_report(anchors, concordance):
    def section(title):
        print(f"\n{title}\n{'-' * len(title)}")

    section("Object")
    print(f"  assays      {anchors['assay']['names']}  "
          f"default={anchors['assay']['default']}  key={anchors['assay']['key']!r}")
    print(f"  layers      {list(anchors['layers'])}")
    print(f"  reductions  {anchors['reductions']['names']}  "
          f"graphs {list(anchors['graphs'])}")
    print(f"  idents      {anchors['idents']['counts']}")
    print(f"  commands    {anchors['commands'] or '(empty)'}")

    section("split / join round trip")
    sj = anchors["split_join"]
    print(f"  layers  {sj['layers_before']} → {sj['layers_after_split']} "
          f"→ {sj['layers_after_join']}")
    for label, key in (("layer name restored", "layer_name_restored"),
                       ("cell order restored", "cell_order_restored"),
                       ("matrix restored", "matrix_restored")):
        print(f"  {label:<22} {sj[key]}")
    jal = anchors["join_all_layers"]
    print(f"  join on prepared assay  "
          f"{jal['error'] or jal['layers_after_join']}")

    section("FetchData")
    for key, value in anchors["fetch"].items():
        print(f"  {key:<16} {value}")

    if concordance is None:
        return
    section("Concordance with R Seurat")
    failed = concordance[~concordance["match"]]
    print(f"  {len(concordance) - len(failed)}/{len(concordance)} anchors match")
    if len(failed):
        print("\n  mismatches:")
        for row in failed.itertuples():
            print(f"    {row.field}\n      python={row.python!r}\n      r     ={row.r!r}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="pbmc3k object-internals tutorial")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    run_full(data_dir=args.data_dir)
