"""Tests for BPCells-style lazy, on-disk matrices (v0.8.0 Scale).

  * write_lazy_matrix / open_lazy_matrix  (lazy.py)
  * LazyMatrix indexing, streaming reductions, materialisation
  * LazyMatrix as a drop-in Assay5 layer

Every LazyMatrix operation is checked against the equivalent in-memory
``scipy.sparse`` result, so the on-disk store is verified to be a faithful,
byte-for-byte drop-in for the sparse layer it replaces. Network-free; all
matrices are written to pytest's ``tmp_path``.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz._sparse import as_dense  # noqa: E402
from shanuz.lazy import (  # noqa: E402
    LazyMatrix,
    is_lazy,
    open_lazy_matrix,
    write_lazy_matrix,
)
from shanuz.shanuz import create_shanuz_object  # noqa: E402


def _rand_sparse(nrow=40, ncol=30, density=0.2, seed=0):
    rng = np.random.default_rng(seed)
    m = sp.random(nrow, ncol, density=density, format="csc", random_state=rng)
    m.data = np.round(m.data * 10.0, 3)  # nicer values, still float64
    return m


# ----------------------------------------------------------------------
# round-trip & metadata
# ----------------------------------------------------------------------


def test_roundtrip_preserves_values(tmp_path):
    m = _rand_sparse()
    lazy = write_lazy_matrix(m, tmp_path / "m.mat")
    assert isinstance(lazy, LazyMatrix)
    assert np.array_equal(lazy.toarray(), m.toarray())
    # reopening the store gives the same values without re-writing
    again = open_lazy_matrix(tmp_path / "m.mat")
    assert np.array_equal(again.toarray(), m.toarray())


def test_reports_shape_dtype_nnz(tmp_path):
    m = _rand_sparse(nrow=17, ncol=23)
    lazy = write_lazy_matrix(m, tmp_path / "m.mat")
    assert lazy.shape == (17, 23)
    assert lazy.nrow == 17 and lazy.ncol == 23
    assert lazy.ndim == 2
    assert lazy.nnz == m.nnz
    assert lazy.dtype == m.data.dtype
    assert len(lazy) == 17


def test_is_lazy(tmp_path):
    lazy = write_lazy_matrix(_rand_sparse(), tmp_path / "m.mat")
    assert is_lazy(lazy)
    assert not is_lazy(_rand_sparse())
    assert not is_lazy(np.zeros((3, 3)))


# ----------------------------------------------------------------------
# indexing — every idiom matches scipy
# ----------------------------------------------------------------------


def test_column_subset_matches_scipy(tmp_path):
    m = _rand_sparse()
    lazy = write_lazy_matrix(m, tmp_path / "m.mat")
    cols = [5, 0, 17, 3, 3]  # out of order, with a repeat
    assert np.array_equal(lazy[:, cols].toarray(), m[:, cols].toarray())


def test_row_subset_matches_scipy(tmp_path):
    m = _rand_sparse()
    lazy = write_lazy_matrix(m, tmp_path / "m.mat")
    rows = [10, 2, 39, 2]
    assert np.array_equal(lazy[rows, :].toarray(), m[rows, :].toarray())


def test_ix_cross_product_matches_scipy(tmp_path):
    m = _rand_sparse()
    lazy = write_lazy_matrix(m, tmp_path / "m.mat")
    r = [1, 8, 20, 33]
    c = [0, 4, 9, 15, 29]
    assert np.array_equal(
        lazy[np.ix_(r, c)].toarray(), m[np.ix_(r, c)].toarray()
    )


def test_contiguous_column_slice_matches_scipy(tmp_path):
    m = _rand_sparse()
    lazy = write_lazy_matrix(m, tmp_path / "m.mat")
    assert np.array_equal(lazy[:, 5:19].toarray(), m[:, 5:19].toarray())
    assert np.array_equal(lazy[:, :].toarray(), m.toarray())


def test_boolean_mask_matches_scipy(tmp_path):
    m = _rand_sparse()
    lazy = write_lazy_matrix(m, tmp_path / "m.mat")
    col_mask = np.zeros(m.shape[1], dtype=bool)
    col_mask[[2, 7, 11]] = True
    assert np.array_equal(lazy[:, col_mask].toarray(), m[:, col_mask].toarray())


# ----------------------------------------------------------------------
# streaming reductions
# ----------------------------------------------------------------------


def test_sum_axes_match_scipy(tmp_path):
    m = _rand_sparse()
    lazy = write_lazy_matrix(m, tmp_path / "m.mat")
    assert lazy.sum() == pytest.approx(float(m.sum()))
    assert np.allclose(lazy.sum(axis=0), np.asarray(m.sum(axis=0)).ravel())
    assert np.allclose(lazy.sum(axis=1), np.asarray(m.sum(axis=1)).ravel())


def test_mean_axes_match_scipy(tmp_path):
    m = _rand_sparse()
    dense = m.toarray()
    lazy = write_lazy_matrix(m, tmp_path / "m.mat")
    assert lazy.mean() == pytest.approx(dense.mean())
    assert np.allclose(lazy.mean(axis=0), dense.mean(axis=0))
    assert np.allclose(lazy.mean(axis=1), dense.mean(axis=1))


def test_nnz_per_col_matches_scipy(tmp_path):
    m = _rand_sparse()
    lazy = write_lazy_matrix(m, tmp_path / "m.mat")
    assert np.array_equal(lazy.nnz_per_col(), np.diff(m.indptr))


def test_col_blocks_reconstruct_matrix(tmp_path):
    m = _rand_sparse(ncol=30)
    lazy = write_lazy_matrix(m, tmp_path / "m.mat")
    seen = 0
    blocks = []
    for start, stop, block in lazy.col_blocks(block_size=7):
        assert block.shape == (m.shape[0], stop - start)
        seen += stop - start
        blocks.append(block)
    assert seen == m.shape[1]
    rebuilt = sp.hstack(blocks).toarray()
    assert np.array_equal(rebuilt, m.toarray())


# ----------------------------------------------------------------------
# materialisation escape hatches
# ----------------------------------------------------------------------


def test_np_asarray_and_as_dense_materialise(tmp_path):
    m = _rand_sparse()
    lazy = write_lazy_matrix(m, tmp_path / "m.mat")
    assert np.array_equal(np.asarray(lazy), m.toarray())
    assert np.array_equal(as_dense(lazy), m.toarray())
    assert np.array_equal(lazy.to_scipy().toarray(), m.toarray())


# ----------------------------------------------------------------------
# write guards & lifecycle
# ----------------------------------------------------------------------


def test_write_refuses_existing_without_overwrite(tmp_path):
    m = _rand_sparse()
    write_lazy_matrix(m, tmp_path / "m.mat")
    with pytest.raises(FileExistsError):
        write_lazy_matrix(m, tmp_path / "m.mat")
    # overwrite replaces cleanly with a different matrix
    m2 = _rand_sparse(nrow=40, ncol=30, seed=99)
    lazy = write_lazy_matrix(m2, tmp_path / "m.mat", overwrite=True)
    assert np.array_equal(lazy.toarray(), m2.toarray())


def test_open_missing_store_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        open_lazy_matrix(tmp_path / "nope")


def test_write_accepts_dense_and_lazy(tmp_path):
    dense = _rand_sparse().toarray()
    lazy = write_lazy_matrix(dense, tmp_path / "d.mat")
    assert np.array_equal(lazy.toarray(), dense)
    # a LazyMatrix can itself be re-persisted
    copy = write_lazy_matrix(lazy, tmp_path / "d2.mat")
    assert np.array_equal(copy.toarray(), dense)


def test_context_manager_closes(tmp_path):
    m = _rand_sparse()
    with write_lazy_matrix(m, tmp_path / "m.mat") as lazy:
        assert np.array_equal(lazy[:, [1, 2]].toarray(), m[:, [1, 2]].toarray())
    assert lazy._data is None  # memmaps released


# ----------------------------------------------------------------------
# drop-in Assay5 layer
# ----------------------------------------------------------------------


def test_lazy_layer_roundtrips_through_assay5(tmp_path):
    rng = np.random.default_rng(0)
    G, N = 25, 40
    counts = sp.csc_matrix(rng.poisson(0.5, size=(G, N)).astype(float))
    features = [f"g{i}" for i in range(G)]
    cells = [f"c{i}" for i in range(N)]
    obj = create_shanuz_object(
        counts=counts, assay="RNA",
        feature_names=features, cell_names=cells,
        meta_data=pd.DataFrame(index=cells),
    )
    assay = obj.get_assay()
    original = assay.layers["counts"].toarray()

    lazy = write_lazy_matrix(assay.layers["counts"], tmp_path / "counts.mat")
    assay.set_layer_data("counts", lazy)
    assert is_lazy(assay.layers["counts"])

    # Full layer still reads back identically...
    assert np.array_equal(as_dense(assay.layer_data("counts")), original)

    # ...and a feature × cell block subset matches the in-memory answer.
    fsub = ["g3", "g10", "g20"]
    csub = ["c1", "c5", "c30", "c39"]
    block = assay.layer_data("counts", cells=csub, features=fsub)
    fidx = [features.index(f) for f in fsub]
    cidx = [cells.index(c) for c in csub]
    assert np.array_equal(as_dense(block), original[np.ix_(fidx, cidx)])
