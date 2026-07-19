"""The out-of-core path through the *pipeline*, not the matrix in isolation.

``test_lazy.py`` checks ``LazyMatrix`` as an object — indexing, reductions,
round-trip — and every one of those tests passed while the feature was
inoperative: five library functions read a lazy layer by densifying all of it
first, so backing a matrix on disk cost more peak memory than leaving it in
RAM, and ``percentage_feature_set`` raised outright. Nothing here is about the
matrix; it is about what the analysis functions do when handed one.

Two properties, held for every step:

  * **Nothing is materialised whole.** The three escape hatches
    (``__array__`` / ``toarray`` / ``to_scipy``) are counted, and the count is
    zero. This is the test the feature never had.
  * **The answer does not change.** A lazy layer and the sparse layer it
    replaces give the same result, bit-for-bit where the arithmetic is
    identical and to a stated tolerance where the summation order differs.
"""
import sys
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import shanuz  # noqa: E402
from shanuz.lazy import LazyMatrix, write_lazy_matrix  # noqa: E402
from shanuz.preprocessing import _log_normalize, _loess2, _vst_hvg  # noqa: E402
from shanuz.shanuz import create_shanuz_object  # noqa: E402


# ----------------------------------------------------------------------
# fixtures & instrumentation
# ----------------------------------------------------------------------


def _counts(n_genes=120, n_cells=80, seed=0):
    """A small counts matrix with the tie structure real data has.

    Integer counts over few cells produce many genes sharing an exact mean,
    which is what the HVG tie handling has to cope with — on pbmc3k 85 % of
    genes share a ``log10(mean)`` with another gene.
    """
    rng = np.random.default_rng(seed)
    dense = rng.poisson(0.4, size=(n_genes, n_cells)).astype(float)
    dense[0, :] = rng.poisson(8.0, size=n_cells)  # one clearly variable gene
    m = sp.csc_matrix(dense)
    m.eliminate_zeros()
    return m


def _object(counts, lazy_path=None):
    genes = [f"g{i}" for i in range(counts.shape[0])]
    cells = [f"c{i}" for i in range(counts.shape[1])]
    obj = create_shanuz_object(
        counts=counts, assay="RNA", feature_names=genes, cell_names=cells,
        min_cells=0, min_features=0, project="lazypipe",
    )
    if lazy_path is not None:
        assay = obj.assays["RNA"]
        assay.layers["counts"] = write_lazy_matrix(
            assay.layers["counts"], lazy_path, overwrite=True
        )
    return obj


class _MaterialisationCounter:
    """Counts whole-store materialisations for the duration of a `with` block.

    ``__array__`` delegates to ``toarray`` which delegates to ``to_scipy``, so
    one logical densification trips three hooks; the count is only ever
    compared against zero, which is unambiguous either way.
    """

    def __init__(self):
        self.count = 0
        self._saved = {}

    def __enter__(self):
        for name in ("__array__", "toarray", "to_scipy"):
            original = getattr(LazyMatrix, name)
            self._saved[name] = original

            def hooked(inner_self, *args, _orig=original, **kwargs):
                self.count += 1
                return _orig(inner_self, *args, **kwargs)

            setattr(LazyMatrix, name, hooked)
        return self

    def __exit__(self, *exc):
        for name, original in self._saved.items():
            setattr(LazyMatrix, name, original)


def test_the_counter_actually_fires():
    """Guard the guard: a counter that never fires would pass every test below."""
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        lazy = write_lazy_matrix(_counts(20, 10), Path(d) / "m", overwrite=True)
        with _MaterialisationCounter() as counter:
            np.asarray(lazy)
        assert counter.count > 0


# ----------------------------------------------------------------------
# the property the feature exists for
# ----------------------------------------------------------------------


def test_building_an_object_from_a_store_never_materialises_it(tmp_path):
    """The obvious way to use the feature: open a store, build an object on it.

    It was also the one path that could not work. `create_assay5_object` ran
    `sp.csc_matrix(np.asarray(matrix))` on anything not already scipy, and
    `calc_n` densified again for `nCount`/`nFeature` — so laziness ended in the
    constructor, before a single analysis function was called. Every test that
    swapped a layer in *after* construction missed it.
    """
    lazy = write_lazy_matrix(_counts(), tmp_path / "counts.lazy", overwrite=True)
    genes = [f"g{i}" for i in range(lazy.shape[0])]
    cells = [f"c{i}" for i in range(lazy.shape[1])]

    with _MaterialisationCounter() as counter:
        obj = create_shanuz_object(
            counts=lazy, assay="RNA", feature_names=genes, cell_names=cells,
            min_cells=0, min_features=0, project="fromstore",
        )

    assert counter.count == 0
    assert isinstance(obj.assays["RNA"].layers["counts"], LazyMatrix)


def test_qc_metadata_matches_between_a_store_and_a_sparse_matrix(tmp_path):
    """`calc_n` has a streaming branch now; it must agree with the dense one."""
    counts = _counts()
    lazy = write_lazy_matrix(counts, tmp_path / "counts.lazy", overwrite=True)
    genes = [f"g{i}" for i in range(counts.shape[0])]
    cells = [f"c{i}" for i in range(counts.shape[1])]

    built = [
        create_shanuz_object(counts=m, assay="RNA", feature_names=genes,
                             cell_names=cells, min_cells=0, min_features=0,
                             project="p")
        for m in (counts, lazy)
    ]
    for column in ("nCount_RNA", "nFeature_RNA"):
        np.testing.assert_array_equal(built[0].meta_data[column].to_numpy(),
                                      built[1].meta_data[column].to_numpy())


def test_min_cells_filtering_agrees_between_a_store_and_a_sparse_matrix(tmp_path):
    """Filtering cannot leave the result on disk, but it must still select the
    same features and cells — the counts it filters on come from a different
    code path for a lazy layer (`nnz_per_row` rather than a CSR indptr diff).

    The matrix is built so the filter actually bites, and asymmetrically: with
    a dense-ish block plus deliberately near-empty genes and cells, using the
    column counts where the row counts belong selects the wrong features
    instead of silently keeping everything.
    """
    rng = np.random.default_rng(3)
    dense = np.zeros((40, 25))
    dense[:12, :] = rng.poisson(1.5, size=(12, 25))     # genes that survive
    dense[12:20, :2] = rng.poisson(2.0, size=(8, 2))    # genes in too few cells
    dense[20:, :] = 0                                   # genes in no cell
    dense[:, 20:] = 0                                   # cells with no genes
    dense[:6, 20:23] = rng.poisson(2.0, size=(6, 3))    # ...but not quite empty
    counts = sp.csc_matrix(dense)
    counts.eliminate_zeros()

    lazy = write_lazy_matrix(counts, tmp_path / "counts.lazy", overwrite=True)
    genes = [f"g{i}" for i in range(counts.shape[0])]
    cells = [f"c{i}" for i in range(counts.shape[1])]

    built = [
        create_shanuz_object(counts=m, assay="RNA", feature_names=genes,
                             cell_names=cells, min_cells=3, min_features=5,
                             project="p")
        for m in (counts, lazy)
    ]
    kept_features = list(built[0].assays["RNA"].features())
    kept_cells = built[0].cell_names()
    assert 0 < len(kept_features) < counts.shape[0], (
        "the fixture must make min_cells actually drop features, or swapping "
        "row counts for column counts would keep everything either way"
    )
    assert 0 < len(kept_cells) < counts.shape[1]

    assert kept_features == list(built[1].assays["RNA"].features())
    assert kept_cells == built[1].cell_names()
    np.testing.assert_array_equal(
        np.asarray(built[0].assays["RNA"].layer_data("counts").todense()),
        np.asarray(built[1].assays["RNA"].layer_data("counts").todense()),
    )


def test_nnz_per_row_matches_scipy(tmp_path):
    counts = _counts()
    lazy = write_lazy_matrix(counts, tmp_path / "counts.lazy", overwrite=True)
    expected = np.diff(counts.tocsr().indptr).astype(np.int64)
    np.testing.assert_array_equal(expected, lazy.nnz_per_row())
    np.testing.assert_array_equal(expected, lazy.nnz_per_row(block_size=3))


def test_whole_pipeline_never_materialises_the_store(tmp_path):
    """The regression test for the defect: every step, zero densifications."""
    obj = _object(_counts(), tmp_path / "counts.lazy")
    assay = obj.assays["RNA"]
    obj.idents = ["a"] * 40 + ["b"] * 40

    with _MaterialisationCounter() as counter:
        shanuz.percentage_feature_set(obj, pattern=r"^g1", col_name="pct")
        shanuz.normalize_data(obj)
        shanuz.find_variable_features(obj, nfeatures=20)
        shanuz.scale_data(obj, features=assay.variable_features)
        shanuz.run_pca(obj, n_pcs=5, features=assay.variable_features,
                       reduction_name="pca")
        shanuz.add_module_score(obj, features=[list(assay.features())[:10]],
                                name="mod")
        shanuz.find_markers(obj, ident_1="a", ident_2="b")

    assert counter.count == 0, (
        f"{counter.count} whole-store materialisations during the pipeline; "
        "an on-disk layer that gets densified has defeated its own purpose"
    )


def _lazy_data_layer(tmp_path, counts=None):
    """An object whose **data** layer is on disk.

    `scale_data`, `find_markers` and `add_module_score` all read `data`, not
    `counts`, so a test that only persists `counts` cannot exercise them —
    `normalize_data` has replaced `data` with an in-memory sparse matrix long
    before they run, and reverting their fix changes nothing observable.
    """
    obj = _object(counts if counts is not None else _counts())
    shanuz.normalize_data(obj)
    assay = obj.assays["RNA"]
    assay.layers["data"] = write_lazy_matrix(
        assay.layers["data"], tmp_path / "data.lazy", overwrite=True
    )
    return obj


def test_downstream_steps_never_materialise_a_lazy_data_layer(tmp_path):
    """The other half of the pipeline, where the lazy layer is `data`."""
    obj = _lazy_data_layer(tmp_path)
    obj.idents = ["a"] * 40 + ["b"] * 40

    with _MaterialisationCounter() as counter:
        shanuz.find_variable_features(obj, nfeatures=20)
        shanuz.scale_data(obj, features=[f"g{i}" for i in range(30)])
        shanuz.add_module_score(obj, features=[[f"g{i}" for i in range(5)]],
                                name="mod", seed=0)
        shanuz.find_markers(obj, ident_1="a", ident_2="b")

    assert counter.count == 0, (
        f"{counter.count} whole-store materialisations reading a lazy data layer"
    )


def test_normalizing_a_lazy_layer_leaves_a_sparse_data_layer(tmp_path):
    """The dense fallback produced an ndarray — 16x the sparse layer's bytes."""
    obj = _object(_counts(), tmp_path / "counts.lazy")
    shanuz.normalize_data(obj)
    data = obj.assays["RNA"].layers["data"]
    assert sp.issparse(data), f"data layer came back as {type(data).__name__}"


def test_counts_layer_survives_the_pipeline_still_lazy(tmp_path):
    obj = _object(_counts(), tmp_path / "counts.lazy")
    shanuz.normalize_data(obj)
    shanuz.find_variable_features(obj, nfeatures=20)
    assert isinstance(obj.assays["RNA"].layers["counts"], LazyMatrix)


# ----------------------------------------------------------------------
# lazy and sparse must agree
# ----------------------------------------------------------------------


def test_log_normalize_is_bit_identical_on_a_lazy_layer(tmp_path):
    """Not `allclose` — the streaming path forms the same products in the
    same order, so anything short of equality is a real difference."""
    counts = _counts()
    lazy = write_lazy_matrix(counts, tmp_path / "m", overwrite=True)

    a = sp.csc_matrix(_log_normalize(counts, 10000.0))
    b = sp.csc_matrix(_log_normalize(lazy, 10000.0))
    a.sort_indices()
    b.sort_indices()

    assert np.array_equal(a.indices, b.indices)
    assert np.array_equal(a.indptr, b.indptr)
    assert np.array_equal(a.data, b.data)


@pytest.mark.parametrize("block_size", [1, 3, 17, 10_000])
def test_streaming_normalize_is_invariant_to_block_size(tmp_path, block_size):
    from shanuz.preprocessing import _log_normalize_lazy

    counts = _counts()
    lazy = write_lazy_matrix(counts, tmp_path / "m", overwrite=True)
    reference = sp.csc_matrix(_log_normalize(counts, 10000.0))
    reference.sort_indices()

    got = sp.csc_matrix(_log_normalize_lazy(lazy, 10000.0, block_size=block_size))
    got.sort_indices()
    assert np.array_equal(reference.data, got.data)


def test_percentage_feature_set_works_on_a_lazy_layer(tmp_path):
    """It raised ValueError: the dense branch ran, but indexing a lazy layer
    returns scipy, whose `.sum(axis=0)` is a (1, n) matrix rather than a vector."""
    counts = _counts()
    sparse_obj = _object(counts)
    lazy_obj = _object(counts, tmp_path / "counts.lazy")

    shanuz.percentage_feature_set(sparse_obj, pattern=r"^g1", col_name="pct")
    shanuz.percentage_feature_set(lazy_obj, pattern=r"^g1", col_name="pct")

    expected = sparse_obj.meta_data["pct"].to_numpy()
    got = lazy_obj.meta_data["pct"].to_numpy()
    assert got.shape == (counts.shape[1],)
    assert np.array_equal(expected, got)


def test_variable_features_are_bit_identical_between_lazy_and_sparse(tmp_path):
    """Bit-identical, not merely close.

    An earlier version ran a streaming reduction for lazy layers and scipy's
    row reduction for sparse ones. They agreed to 1e-14 — which sounds like
    enough, and is not: `variance.standardized` carries exact ties, a tie-break
    decides which of two tied genes is selected, and genes tied under one
    summation order are not tied under the other. On pbmc3k that reordered 147
    of 2000 features, which changed the PCA row order and, downstream, gave
    9 clusters against 8. One implementation now serves both.
    """
    counts = _counts()
    lazy = write_lazy_matrix(counts, tmp_path / "m", overwrite=True)

    idx_sparse, mean_s, var_s, vst_s = _vst_hvg(counts, nfeatures=20)
    idx_lazy, mean_l, var_l, vst_l = _vst_hvg(lazy, nfeatures=20)

    np.testing.assert_array_equal(idx_sparse, idx_lazy)
    np.testing.assert_array_equal(mean_s, mean_l)
    np.testing.assert_array_equal(var_s, var_l)
    np.testing.assert_array_equal(vst_s, vst_l)


def test_scale_data_agrees_between_lazy_and_sparse(tmp_path):
    counts = _counts()
    sparse_obj = _object(counts)
    shanuz.normalize_data(sparse_obj)
    lazy_obj = _lazy_data_layer(tmp_path, counts)

    features = [f"g{i}" for i in range(30)]
    for obj in (sparse_obj, lazy_obj):
        shanuz.scale_data(obj, features=features)

    a = sparse_obj.assays["RNA"].layer_data("scale.data")
    b = lazy_obj.assays["RNA"].layer_data("scale.data")
    np.testing.assert_array_equal(np.asarray(a.todense()), np.asarray(b.todense()))


def test_find_markers_agrees_between_lazy_and_sparse(tmp_path):
    counts = _counts()
    sparse_obj = _object(counts)
    shanuz.normalize_data(sparse_obj)
    lazy_obj = _lazy_data_layer(tmp_path, counts)

    out = []
    for obj in (sparse_obj, lazy_obj):
        obj.idents = ["a"] * 40 + ["b"] * 40
        out.append(shanuz.find_markers(obj, ident_1="a", ident_2="b"))

    assert list(out[0].index) == list(out[1].index)
    np.testing.assert_array_equal(out[0]["avg_log2FC"].to_numpy(),
                                  out[1]["avg_log2FC"].to_numpy())


def test_add_module_score_agrees_between_lazy_and_sparse(tmp_path):
    counts = _counts()
    sparse_obj = _object(counts)
    shanuz.normalize_data(sparse_obj)
    lazy_obj = _lazy_data_layer(tmp_path, counts)

    program = [f"g{i}" for i in range(5)]
    for obj in (sparse_obj, lazy_obj):
        shanuz.add_module_score(obj, features=[program], name="mod", seed=0)

    np.testing.assert_array_equal(sparse_obj.meta_data["mod1"].to_numpy(),
                                  lazy_obj.meta_data["mod1"].to_numpy())


# ----------------------------------------------------------------------
# LOESS determinism — the amplifier that made the above hard to test
# ----------------------------------------------------------------------


def _loess_input(seed=0, n=400):
    """x with heavy exact ties, as counts-derived log10(mean) always has."""
    rng = np.random.default_rng(seed)
    x = np.round(rng.normal(0, 1, size=n), 2)  # rounding manufactures the ties
    y = 0.7 * x + rng.normal(0, 0.15, size=n)
    return x, y


def test_loess_gives_one_fitted_value_per_distinct_x():
    """R's loess cannot return two fitted values for one x — its within-x
    spread on pbmc3k is 1.8e-15, pure interpolation noise. shanuz's was
    1.3e-3, because the window was chosen by position in an unstably sorted
    array, so members of a tied run got different neighbourhoods."""
    x, y = _loess_input()
    assert len(np.unique(x)) < len(x), "fixture must contain ties to be a test"

    fitted = _loess2(x, y, frac=0.3)
    for value in np.unique(x):
        spread = np.ptp(fitted[x == value])
        assert spread == 0.0, f"x={value} got fitted values spanning {spread:.3e}"


def test_loess_is_invariant_to_input_row_order():
    x, y = _loess_input()
    rng = np.random.default_rng(1)
    perm = rng.permutation(len(x))
    inverse = np.empty_like(perm)
    inverse[perm] = np.arange(len(x))

    np.testing.assert_array_equal(
        _loess2(x, y, frac=0.3), _loess2(x[perm], y[perm], frac=0.3)[inverse]
    )


def test_loess_does_not_amplify_a_rounding_level_perturbation():
    """A 1e-15 nudge to x moved fitted values by up to 29 %, because it
    reshuffled tied runs under an unstable sort. Any 1e-14-scale difference —
    a different BLAS, a different summation order, a different platform —
    was enough to move the HVG statistics."""
    x, y = _loess_input()
    rng = np.random.default_rng(2)
    nudged = x * (1 + rng.choice([-1.0, 1.0], size=x.shape) * 1e-15)

    base = _loess2(x, y, frac=0.3)
    moved = _loess2(nudged, y, frac=0.3)
    scale = np.maximum(np.abs(base), np.abs(moved))
    nonzero = scale > 0
    worst = float((np.abs(base - moved)[nonzero] / scale[nonzero]).max())
    assert worst < 1e-9, f"1e-15 input change moved the fit by {worst:.3e}"


def test_hvg_ties_break_by_ascending_index_like_r_order():
    """`head(order(x, decreasing = TRUE), n)` breaks ties by ascending original
    index; `argsort(x)[::-1]` breaks them descending, and unstably."""
    counts = _counts()
    idx, _, _, var_standardized = _vst_hvg(counts, nfeatures=len(counts.toarray()))

    ranked = var_standardized[idx]
    assert np.all(np.diff(ranked) <= 0), "must be ordered by descending score"
    for value in np.unique(ranked):
        tied = idx[ranked == value]
        if len(tied) > 1:
            assert list(tied) == sorted(tied), (
                f"ties at {value} came out as {list(tied)}, not ascending index"
            )
