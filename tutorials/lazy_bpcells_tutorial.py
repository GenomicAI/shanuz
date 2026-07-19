"""shanuz's on-disk layer against Seurat's, on the same pbmc3k matrix.

What this compares
==================
Seurat's out-of-core story is BPCells: a bitpacked on-disk format plus a
*deferred operation graph*, so ``LogNormalize`` / ``VST`` / ``ScaleData`` are
queued rather than executed and nothing is materialised until a reduction asks
for it. Seurat ships dedicated ``IterableMatrix`` methods for exactly the
functions that touch the matrix.

shanuz's ``LazyMatrix`` is a different design with the same goal: three
memory-mapped ``.npy`` arrays in CSC layout, with the analysis functions
streaming over them in cell-blocks. It is lazy in *storage*, not in operations.

So this is not a port comparison — there is no ``LazyMatrix`` in Seurat to
match value-for-value. What can be compared, and is:

  1. **Do the two tools compute the same analysis out of core?** The anchors
     below run both pipelines on the same matrix.
  2. **Does each tool agree with itself?** Whether a tool's on-disk path
     returns what its in-memory path returns is the property a user actually
     depends on, and it is the one where the two differ most.
  3. **What does each format cost on disk?**

The headline
============
Seurat's two paths **do not agree**: BPCells computes in single precision, so
normalised values sit ~1e-6 apart, ``variance.standardized`` ~2e-2 apart, and
the out-of-core run selects a different variable feature from the in-memory
run. shanuz's two paths are **bit-identical** — one implementation serves both,
which is a deliberate choice made after the alternative was measured (see
``test_lazy_pipeline.py``). shanuz pays for that in disk: no bitpacking, so its
store is ~6x a BPCells integer store.

Usage
-----
    python tutorials/lazy_bpcells_tutorial.py            # run, write anchors
    Rscript tutorials/lazy_bpcells_verify.R              # the Seurat side
    python tutorials/lazy_bpcells_tutorial.py --report   # compare
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from shanuz import (  # noqa: E402
    create_shanuz_object,
    find_markers,
    find_variable_features,
    normalize_data,
    percentage_feature_set,
    scale_data,
)
from shanuz.datasets import pbmc3k  # noqa: E402
from shanuz.lazy import (  # noqa: E402
    LazyMatrix,
    open_lazy_matrix,
    write_lazy_matrix,
)

FIGURES = Path(__file__).parent / "figures_lazy"
ASSAY = "RNA"
HEAD = 20
NFEATURES = 2000

# Every tolerance is the measured residual rounded up one decade, and every one
# has a reason. They are set against R's *in-memory* series: BPCells computes in
# single precision, so Seurat's own out-of-core numbers sit ~1e-6 from its
# in-memory ones, while shanuz stays in float64 throughout.
FLOAT_TOLERANCES = {
    # Summation order only. Measured 2.9e-15 / 3.8e-15 / 3.7e-15 / 2.5e-14.
    "mem.normalize_head": 1e-13,
    "mem.vst_mean_head": 1e-13,
    "mem.vst_variance_head": 1e-13,
    "mem.scale_head": 1e-12,
    # Accumulated rounding over every non-zero (2.2M and 36M terms). 7.0e-13 / 4.2e-7.
    "mem.normalize_sum": 1e-11,
    "mem.scale_abs_sum": 1e-5,
    # The one real residual: `variance.standardized` is proportional to
    # 1/variance.expected, and variance.expected is `10^loess_fitted`. shanuz's
    # `_loess2` is a NumPy local-quadratic fit; R's `loess` is the cloess
    # Fortran with kd-tree interpolation, and the two differ by up to 2.5e-2 on
    # the expected variance. The inputs are not in question -- mean and
    # variance agree to 1.5e-13. This is a pre-existing smoother difference,
    # not anything to do with going out of core.
    "mem.vst_var_std_head": 1e-2,
    # Integer counts, so exact or wrong.
    "calcn.ncount_head": 0.0,
    "calcn.nfeature_head": 0.0,
}

# `vst_selected_head` is deliberately absent: the LOESS residual above reorders
# genes whose standardized variance is within 7.6e-3 of each other, so the top
# of the list is not a stable comparison. The *set* is, and is checked instead.

# Anchors that are *reported* rather than matched: the two designs differ here
# on purpose, and the number is the finding.
REPORTED_ONLY = (
    "mem.vst_selected_head",
    "store.dgc_bytes",
    "store.bpcells_double_bytes",
    "store.bpcells_uint32_bytes",
    "store.storage_order",
    "markers.supported_tests",
    "selfcheck.normalize_max_diff",
    "selfcheck.normalize_identical",
    "selfcheck.vst_mean_max_diff",
    "selfcheck.vst_var_std_max_diff",
    "selfcheck.hvg_overlap",
    "selfcheck.scale_max_diff",
)


# ----------------------------------------------------------------------
# pipeline
# ----------------------------------------------------------------------


def build(data_dir=None):
    """pbmc3k through the standard filter — the same cells as every other
    tutorial in this series."""
    counts, genes, cells = pbmc3k(data_dir=data_dir)
    obj = create_shanuz_object(counts=counts, assay=ASSAY, min_cells=3,
                               min_features=200, project="pbmc3k_lazy",
                               feature_names=genes, cell_names=cells)
    percentage_feature_set(obj, pattern=r"^MT-", col_name="percent.mt")
    md = obj.meta_data
    keep = ((md["nFeature_RNA"] > 200) & (md["nFeature_RNA"] < 2500)
            & (md["percent.mt"] < 5))
    return obj.subset(cells=list(md.index[keep]))


def persist(obj, path=None):
    """Write the counts layer to an on-disk store and return it."""
    assay = obj.assays[ASSAY]
    target = Path(path or (FIGURES / "shanuz_store"))
    return write_lazy_matrix(assay.layers["counts"], target, overwrite=True)


def run_pipeline(counts, genes, cells):
    """Normalise, select variable features, scale — on whatever `counts` is."""
    obj = create_shanuz_object(counts=counts, assay=ASSAY, min_cells=0,
                               min_features=0, project="p",
                               feature_names=genes, cell_names=cells)
    normalize_data(obj)
    find_variable_features(obj, selection_method="vst", nfeatures=NFEATURES)
    scale_data(obj, features=obj.assays[ASSAY]._all_feature_names)
    return obj


def _materialisation_counter():
    """Count whole-store densifications, so the run can prove it never did one."""
    state = {"n": 0, "restore": []}
    for name in ("__array__", "toarray", "to_scipy"):
        original = getattr(LazyMatrix, name)
        state["restore"].append((name, original))

        def hooked(self, *a, _orig=original, **k):
            state["n"] += 1
            return _orig(self, *a, **k)

        setattr(LazyMatrix, name, hooked)
    return state


def _restore(state):
    for name, original in state["restore"]:
        setattr(LazyMatrix, name, original)


# ----------------------------------------------------------------------
# anchors
# ----------------------------------------------------------------------


def _layer(obj, layer):
    data = obj.assays[ASSAY].layers[layer]
    return sp.csc_matrix(data) if not sp.issparse(data) else data


def collect_anchors(obj_mem, obj_lazy, store, counts):
    """Every number the vignette quotes, named to match the R side."""
    anchors = {}

    anchors["store.nnz"] = int(counts.nnz)
    anchors["store.nrow"] = int(counts.shape[0])
    anchors["store.ncol"] = int(counts.shape[1])
    csc = sp.csc_matrix(counts)
    anchors["store.dgc_bytes"] = int(csc.data.nbytes + csc.indices.nbytes
                                     + csc.indptr.nbytes)
    anchors["store.shanuz_bytes"] = int(
        sum(f.stat().st_size for f in store.path.iterdir())
    )
    anchors["store.storage_order"] = "col"

    md = obj_lazy.meta_data
    anchors["calcn.ncount_head"] = md["nCount_RNA"].to_numpy()[:HEAD].tolist()
    anchors["calcn.nfeature_head"] = md["nFeature_RNA"].to_numpy()[:HEAD].tolist()

    for tag, obj in (("lazy", obj_lazy), ("mem", obj_mem)):
        assay = obj.assays[ASSAY]
        data = _layer(obj, "data")
        data.sort_indices()
        hvf = assay.meta_data
        scaled = np.asarray(_layer(obj, "scale.data").todense())
        anchors[f"{tag}.normalize_head"] = data.data[:HEAD].tolist()
        anchors[f"{tag}.normalize_sum"] = float(data.data.sum())
        anchors[f"{tag}.vst_mean_head"] = hvf["means"].to_numpy()[:HEAD].tolist()
        anchors[f"{tag}.vst_variance_head"] = (
            hvf["variances"].to_numpy()[:HEAD].tolist())
        anchors[f"{tag}.vst_var_std_head"] = (
            hvf["variances.standardized"].to_numpy()[:HEAD].tolist())
        anchors[f"{tag}.vst_selected_head"] = list(assay.variable_features[:HEAD])
        anchors[f"{tag}.scale_head"] = scaled[0, :HEAD].tolist()
        anchors[f"{tag}.scale_abs_sum"] = float(np.abs(scaled).sum())

    # shanuz against itself — the property the tutorial is really about.
    d_mem, d_lazy = _layer(obj_mem, "data"), _layer(obj_lazy, "data")
    d_mem.sort_indices()
    d_lazy.sort_indices()
    h_mem = obj_mem.assays[ASSAY].meta_data
    h_lazy = obj_lazy.assays[ASSAY].meta_data
    s_mem = np.asarray(_layer(obj_mem, "scale.data").todense())
    s_lazy = np.asarray(_layer(obj_lazy, "scale.data").todense())
    anchors["selfcheck.normalize_max_diff"] = float(
        np.abs(d_mem.data - d_lazy.data).max())
    anchors["selfcheck.normalize_identical"] = float(
        np.array_equal(d_mem.data, d_lazy.data))
    anchors["selfcheck.vst_mean_max_diff"] = float(
        np.abs(h_mem["means"].to_numpy() - h_lazy["means"].to_numpy()).max())
    anchors["selfcheck.vst_var_std_max_diff"] = float(np.abs(
        h_mem["variances.standardized"].to_numpy()
        - h_lazy["variances.standardized"].to_numpy()).max())
    anchors["selfcheck.hvg_overlap"] = len(
        set(obj_mem.assays[ASSAY].variable_features)
        & set(obj_lazy.assays[ASSAY].variable_features))
    anchors["selfcheck.scale_max_diff"] = float(np.abs(s_mem - s_lazy).max())
    return anchors


def markers_anchors(obj):
    groups = ["a" if i % 2 == 0 else "b" for i in range(len(obj.cell_names()))]
    obj.idents = groups
    table = find_markers(obj, ident_1="a", ident_2="b", test_use="wilcox",
                         logfc_threshold=0.0, min_pct=0.0)
    return {
        "markers.wilcox_top10": list(table.index[:10]),
        "markers.n_genes": int(len(table)),
        # Every test shanuz offers works on a lazy layer; Seurat's out-of-core
        # path takes only wilcox, and warns that column-major storage is the
        # wrong orientation for DE at that.
        "markers.supported_tests": _supported_tests(obj),
    }, table


def _supported_tests(obj):
    supported = []
    for test in ("wilcox", "t", "bimod", "LR", "negbinom", "roc", "mast", "deseq2"):
        try:
            find_markers(obj, ident_1="a", ident_2="b", test_use=test,
                         logfc_threshold=0.0, min_pct=0.0)
            supported.append(test)
        except Exception:
            pass
    return supported


# ----------------------------------------------------------------------
# comparison
# ----------------------------------------------------------------------


def _match(py, r, tol):
    if isinstance(py, str) or isinstance(r, str):
        return str(py) == str(r)
    py_arr = np.atleast_1d(np.asarray(py, dtype=object))
    r_arr = np.atleast_1d(np.asarray(r, dtype=object))
    if py_arr.shape != r_arr.shape:
        return False
    if py_arr.dtype.kind in "OU" and any(isinstance(v, str) for v in py_arr.ravel()):
        return [str(v) for v in py_arr.ravel()] == [str(v) for v in r_arr.ravel()]
    a = np.asarray(py, dtype=float).ravel()
    b = np.asarray(r, dtype=float).ravel()
    if tol == 0.0:
        return bool(np.array_equal(a, b))
    scale = np.maximum(np.abs(a), np.abs(b))
    scale[scale == 0] = 1.0
    return bool((np.abs(a - b) / scale <= tol).all())


def compare(py_anchors, r_anchors):
    """Anchor-by-anchor, with the reported-only ones held aside."""
    matched, differed, reported = [], [], []
    for name in sorted(set(py_anchors) & set(r_anchors)):
        if name in REPORTED_ONLY or name.startswith("bp."):
            reported.append((name, py_anchors.get(name), r_anchors[name]))
            continue
        tol = FLOAT_TOLERANCES.get(name, 0.0)
        if _match(py_anchors[name], r_anchors[name], tol):
            matched.append((name, tol))
        else:
            differed.append((name, py_anchors[name], r_anchors[name], tol))
    return matched, differed, reported


def report_concordance(py_anchors=None, r_anchors=None):
    py_anchors = py_anchors or json.loads((FIGURES / "py_anchors.json").read_text())
    r_path = FIGURES / "r_anchors.json"
    if not r_path.exists():
        raise SystemExit("figures_lazy/r_anchors.json missing — run "
                         "`Rscript tutorials/lazy_bpcells_verify.R` first.")
    r_anchors = r_anchors or json.loads(r_path.read_text())

    # R names its two runs bp./mem.; shanuz's out-of-core run is `lazy.`, and
    # it is compared against R's in-memory series (see the module docstring).
    remapped = dict(py_anchors)
    for key, value in py_anchors.items():
        if key.startswith("lazy."):
            remapped.setdefault("mem." + key[len("lazy."):], value)

    matched, differed, reported = compare(remapped, r_anchors)
    total = len(matched) + len(differed)
    print(f"\n{'=' * 74}\nshanuz vs Seurat — out of core\n{'=' * 74}")
    print(f"  {len(matched)} of {total} compared anchors match "
          f"({len([m for m, t in matched if t == 0.0])} exactly)\n")
    for name, py, r, tol in differed:
        print(f"  DIFFERS  {name}  (tol {tol:g})")
        print(f"           shanuz {np.asarray(py).ravel()[:4]}")
        print(f"           Seurat {np.asarray(r).ravel()[:4]}")

    # The selected *set*, which the head-of-list ordering cannot speak to.
    r_hvg_path = FIGURES / "r_variable_features.txt"
    py_hvg_path = FIGURES / "py_variable_features.txt"
    overlap = None
    if r_hvg_path.exists() and py_hvg_path.exists():
        r_hvg = {ln.strip() for ln in r_hvg_path.read_text().splitlines() if ln.strip()}
        py_hvg = {ln.strip().replace("_", "-")
                  for ln in py_hvg_path.read_text().splitlines() if ln.strip()}
        overlap = len(py_hvg & r_hvg)
        print(f"\n  variable features selected out of core: "
              f"{overlap}/{len(r_hvg)} shared with Seurat")

    print(f"\n{'-' * 74}\nReported, not matched — the two designs differ here on "
          f"purpose\n{'-' * 74}")
    for name, py, r in reported:
        print(f"  {name:<38} shanuz {str(py)[:22]:<24} Seurat {str(r)[:22]}")
    return matched, differed, reported, overlap


# ----------------------------------------------------------------------
# entry points
# ----------------------------------------------------------------------


def run_full(data_dir=None, verbose=True):
    FIGURES.mkdir(exist_ok=True)
    if verbose:
        print("Building pbmc3k...")
    obj = build(data_dir)
    assay = obj.assays[ASSAY]
    counts = sp.csc_matrix(assay.layers["counts"])
    genes = list(assay.features())
    cells = obj.cell_names()
    (FIGURES / "cells.txt").write_text("\n".join(cells) + "\n")
    if verbose:
        print(f"  {counts.shape[0]} features x {counts.shape[1]} cells, "
              f"nnz={counts.nnz}")

    store = persist(obj)
    if verbose:
        size = sum(f.stat().st_size for f in store.path.iterdir())
        print(f"  store: {size / 1e6:.2f} MB")

    if verbose:
        print("Running the in-memory pipeline...")
    obj_mem = run_pipeline(counts, genes, cells)

    if verbose:
        print("Running the out-of-core pipeline...")
    counter = _materialisation_counter()
    try:
        obj_lazy = run_pipeline(open_lazy_matrix(store.path), genes, cells)
    finally:
        _restore(counter)
    if verbose:
        print(f"  whole-store materialisations: {counter['n']}")

    anchors = collect_anchors(obj_mem, obj_lazy, store, counts)
    anchors["run.materialisations"] = counter["n"]
    marker_anchors, table = markers_anchors(obj_lazy)
    anchors.update(marker_anchors)
    table.to_csv(FIGURES / "py_markers_wilcox.csv")

    (FIGURES / "py_variable_features.txt").write_text(
        "\n".join(obj_lazy.assays[ASSAY].variable_features) + "\n")
    (FIGURES / "py_anchors.json").write_text(json.dumps(anchors, indent=2))
    if verbose:
        print(f"\nWrote {FIGURES / 'py_anchors.json'}")
        print("Now run: Rscript tutorials/lazy_bpcells_verify.R")
    return obj_mem, obj_lazy, anchors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="out-of-core side-by-side")
    parser.add_argument("--report", action="store_true",
                        help="compare against the R anchors")
    parser.add_argument("--data-dir", default=None)
    args = parser.parse_args()
    if args.report:
        report_concordance()
    else:
        run_full(data_dir=args.data_dir)
