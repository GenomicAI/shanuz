"""What the BPCells side-by-side established, pinned so it cannot quietly rot.

Constants marked "R:" were read off a live Seurat 5.5.1 / BPCells 0.3.1 session
on pbmc3k (13,714 x 2,638), not derived from shanuz's own output. The tests
themselves need neither R nor BPCells — CI has no R — so the R side is
transcribed rather than invoked, the same way `test_de_parity.py` does it.

The finding worth protecting: **shanuz's on-disk and in-memory paths are
bit-identical, and Seurat's are not.** BPCells computes in single precision, so
Seurat's out-of-core run returns different numbers from its own in-memory run
and selects a different variable feature. That is the property a user actually
leans on when they move a matrix to disk, and it is the one shanuz should not
lose.
"""
import json
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

from shanuz import (
    create_shanuz_object,
    find_variable_features,
    normalize_data,
    scale_data,
)
from shanuz.lazy import write_lazy_matrix

TUTORIALS = Path(__file__).resolve().parent.parent / "tutorials"

# R: Seurat 5.5.1 with BPCells 0.3.1, in-memory run vs out-of-core run,
# identical input matrix. `Rscript tutorials/lazy_bpcells_verify.R` reproduces.
SEURAT_SELF_NORMALIZE_MAX_DIFF = 1.0153110547861e-06
SEURAT_SELF_SCALE_MAX_DIFF = 1.7167773513904e-06
SEURAT_SELF_VAR_STD_MAX_DIFF = 0.0205553676092025
SEURAT_SELF_HVG_OVERLAP = 1999           # out of 2000
# R: FindMarkers on an IterableMatrix, all eight tests attempted.
SEURAT_OUT_OF_CORE_TESTS = ("wilcox",)
# R: store sizes on the same matrix, bytes.
BPCELLS_UINT32_BYTES = 4_501_708
DGC_ARRAY_BYTES = 26_875_340


def _counts(n_genes=90, n_cells=60, seed=0):
    rng = np.random.default_rng(seed)
    dense = rng.poisson(0.5, size=(n_genes, n_cells)).astype(float)
    dense[0, :] = rng.poisson(9.0, size=n_cells)
    m = sp.csc_matrix(dense)
    m.eliminate_zeros()
    return m


def _pipeline(counts, genes, cells):
    obj = create_shanuz_object(counts=counts, assay="RNA", feature_names=genes,
                               cell_names=cells, min_cells=0, min_features=0,
                               project="p")
    normalize_data(obj)
    find_variable_features(obj, selection_method="vst", nfeatures=20)
    scale_data(obj, features=obj.assays["RNA"]._all_feature_names)
    return obj


# ---------------------------------------------------------------------------
# the thesis
# ---------------------------------------------------------------------------


def test_shanuz_on_disk_and_in_memory_paths_are_bit_identical(tmp_path):
    """Seurat's differ by 1.0e-06 on normalise and 2.1e-02 on
    variance.standardized; shanuz's must differ by nothing at all."""
    counts = _counts()
    genes = [f"g{i}" for i in range(counts.shape[0])]
    cells = [f"c{i}" for i in range(counts.shape[1])]
    store = write_lazy_matrix(counts, tmp_path / "s", overwrite=True)

    mem = _pipeline(counts, genes, cells)
    lazy = _pipeline(store, genes, cells)

    d_mem = sp.csc_matrix(mem.assays["RNA"].layers["data"])
    d_lazy = sp.csc_matrix(lazy.assays["RNA"].layers["data"])
    d_mem.sort_indices()
    d_lazy.sort_indices()
    np.testing.assert_array_equal(d_mem.data, d_lazy.data)

    h_mem = mem.assays["RNA"].meta_data
    h_lazy = lazy.assays["RNA"].meta_data
    for column in ("means", "variances", "variances.standardized"):
        np.testing.assert_array_equal(h_mem[column].to_numpy(),
                                      h_lazy[column].to_numpy())

    assert (mem.assays["RNA"].variable_features
            == lazy.assays["RNA"].variable_features)
    np.testing.assert_array_equal(
        np.asarray(sp.csc_matrix(mem.assays["RNA"].layers["scale.data"]).todense()),
        np.asarray(sp.csc_matrix(lazy.assays["RNA"].layers["scale.data"]).todense()),
    )


def test_seurats_two_paths_are_known_to_disagree():
    """Documents the comparison's point. If Seurat ever closes these, the
    tutorial's headline needs rewriting rather than silently going stale."""
    assert SEURAT_SELF_NORMALIZE_MAX_DIFF > 1e-7, (
        "BPCells' single precision is what makes Seurat's two paths differ"
    )
    assert SEURAT_SELF_VAR_STD_MAX_DIFF > 1e-3
    assert SEURAT_SELF_HVG_OVERLAP < 2000, (
        "Seurat's out-of-core run selected a different variable feature"
    )


def test_shanuz_supports_more_de_tests_out_of_core_than_seurat():
    """Seurat's IterableMatrix path takes wilcox alone — `FindMarkers.default`
    raises for anything else, and warns that column-major storage is the wrong
    orientation for DE at that."""
    from shanuz.markers import find_markers  # noqa: F401

    assert SEURAT_OUT_OF_CORE_TESTS == ("wilcox",)


# ---------------------------------------------------------------------------
# the storage gap, stated as a number rather than a vibe
# ---------------------------------------------------------------------------


def test_shanuz_store_is_uncompressed_and_that_costs_what_it_costs(tmp_path):
    """shanuz writes the three CSC arrays verbatim, so the store is the size of
    the matrix. BPCells bitpacks: 4.50 MB against 26.88 MB of dgCMatrix arrays
    on pbmc3k. Closing that means implementing BP128 delta encoding, which is
    not a tolerance to loosen — it is a format."""
    counts = _counts()
    store = write_lazy_matrix(counts, tmp_path / "s", overwrite=True)
    on_disk = sum(f.stat().st_size for f in store.path.iterdir())
    in_memory = counts.data.nbytes + counts.indices.nbytes + counts.indptr.nbytes

    # Fixed overhead only: three .npy headers plus meta.json, ~490 bytes,
    # independent of the matrix. A ratio would pass for the wrong reason on a
    # large fixture and fail for the wrong reason on a small one.
    overhead = on_disk - in_memory
    assert 0 <= overhead < 2048, (
        f"store is {overhead} bytes over the in-memory arrays; it is expected "
        "to be those arrays verbatim plus a small header. If this changed, a "
        "compression scheme was added and the tutorial's storage comparison "
        "needs remeasuring"
    )
    assert BPCELLS_UINT32_BYTES * 4 < DGC_ARRAY_BYTES, (
        "BPCells' integer store is several times smaller — the gap is real"
    )


# ---------------------------------------------------------------------------
# the tutorial's own comparison machinery
# ---------------------------------------------------------------------------


def _tutorial():
    import sys

    root = str(TUTORIALS.parent)
    if root not in sys.path:
        sys.path.insert(0, root)
    from tutorials import lazy_bpcells_tutorial

    return lazy_bpcells_tutorial


def test_every_tolerance_is_justified_and_bounded():
    """A tolerance loose enough to hide a defect is not a comparison. The
    loosest here is 1e-2, on the one anchor with a named cause."""
    t = _tutorial()
    assert t.FLOAT_TOLERANCES, "tolerances must be declared, not implicit"
    for name, tol in t.FLOAT_TOLERANCES.items():
        assert 0.0 <= tol <= 1e-2, f"{name} tolerance {tol:g} is too loose to mean anything"
    assert t.FLOAT_TOLERANCES["calcn.ncount_head"] == 0.0, (
        "integer counts must match exactly"
    )
    assert t.FLOAT_TOLERANCES["mem.vst_var_std_head"] == 1e-2, (
        "the LOESS residual is 7.6e-3; the tolerance should sit just above it"
    )


def test_tolerance_lookup_is_exact_not_prefix():
    """`mem.vst_var_std_head` is the loosest tolerance in the table. A prefix
    match would hand it to every other `mem.vst_*` anchor."""
    t = _tutorial()
    matched, differed, reported = t.compare(
        {"mem.vst_mean_head": [1.0]}, {"mem.vst_mean_head": [1.0 + 5e-3]}
    )
    assert differed, "a 5e-3 error on vst_mean must not inherit vst_var_std's 1e-2"


def test_reported_anchors_are_not_counted_as_matches():
    """The store sizes and self-check numbers differ on purpose; counting them
    as matches or as failures would both be wrong."""
    t = _tutorial()
    py = {"store.dgc_bytes": 1, "selfcheck.normalize_max_diff": 0.0}
    r = {"store.dgc_bytes": 999, "selfcheck.normalize_max_diff": 1e-6}
    matched, differed, reported = t.compare(py, r)
    assert not matched and not differed
    assert len(reported) == 2


def test_exact_anchors_use_exact_comparison():
    t = _tutorial()
    _, differed, _ = t.compare({"calcn.ncount_head": [10.0]},
                               {"calcn.ncount_head": [10.000000001]})
    assert differed, "a zero tolerance must reject any difference at all"


@pytest.mark.skipif(not (TUTORIALS / "figures_lazy" / "r_anchors.json").exists(),
                    reason="R anchors absent; run tutorials/lazy_bpcells_verify.R")
def test_recorded_run_matches_every_compared_anchor():
    """Guards the committed result, when a local R run has produced anchors."""
    t = _tutorial()
    fig = TUTORIALS / "figures_lazy"
    py = json.loads((fig / "py_anchors.json").read_text())
    r = json.loads((fig / "r_anchors.json").read_text())
    remapped = dict(py)
    for key, value in py.items():
        if key.startswith("lazy."):
            remapped.setdefault("mem." + key[len("lazy."):], value)
    matched, differed, reported = t.compare(remapped, r)
    assert not differed, f"anchors drifted: {[d[0] for d in differed]}"
    assert len(matched) >= 14
