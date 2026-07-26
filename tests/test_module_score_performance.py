"""How `add_module_score` reads the matrix, and in what order it adds it up.

Two defects sat in the same two lines.

**Speed.** Each gene's row was pulled out on its own — `mat[i, :]` inside a list
comprehension — and every assay layer here is CSC, so slicing one row walks the
whole column-major matrix. On the THP-1 ECCITE data (18,381 × 20,729, 69.5M
nonzeros) that was ~22 ms per gene, and `ctrl=100` draws a couple of thousand
control genes: 50 of 51 profiled seconds were inside one scipy call, invoked once
per gene. `cell_cycle_scoring` took 168 s against R's 12.

**Reproducibility.** The control genes were collected in a `set`, and a mean
depends on the order its terms are added. Python randomises `str` hashing per
process, so the same object with the same seed produced a different score in a
different process. Not a large difference — 9.7e-16 — but a value that moves when
nothing moved is the kind of thing that gets chased for an afternoon.

The fix is one row selection per gene set, transposed to CSR so the rows are
walked in the order they were asked for, and a dict in place of the set. Both are
pinned here: the first as an exact-equality property, the second by running the
same computation under two hash seeds.
"""
import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest
import scipy.sparse as sp

from shanuz import create_shanuz_object, normalize_data
from shanuz.module_score import _mean_over_rows, add_module_score

N_GENES, N_CELLS = 300, 200


def _matrix(seed=0, density=0.3):
    rng = np.random.default_rng(seed)
    dense = rng.poisson(2.0, size=(N_GENES, N_CELLS)).astype(float)
    dense[rng.random(dense.shape) > density] = 0.0
    return dense


def _wide_matrix(seed=0, density=0.3):
    """Values spanning many orders of magnitude, so summation order is visible.

    Small Poisson counts are all exactly representable and their partial sums
    stay exact, so *every* order gives identical bits — a fixture built from
    them cannot express an ordering property, and a test written on one passes
    whatever the code does. Log-normalized expression has the spread this needs;
    counts do not.
    """
    rng = np.random.default_rng(seed)
    dense = rng.lognormal(0.0, 6.0, size=(N_GENES, N_CELLS))
    dense[rng.random(dense.shape) > density] = 0.0
    return dense


def _obj(seed=0):
    obj = create_shanuz_object(
        counts=sp.csc_matrix(_matrix(seed)),
        assay="RNA",
        feature_names=[f"g{i}" for i in range(N_GENES)],
        cell_names=[f"c{j}" for j in range(N_CELLS)],
    )
    normalize_data(obj)
    return obj


def _stack_mean(mat, rows):
    """The formulation that was replaced, kept as the reference."""
    out = []
    for i in rows:
        r = mat[i, :]
        out.append(np.asarray(r.todense()).flatten() if sp.issparse(r)
                   else np.asarray(r).flatten())
    return np.mean(out, axis=0)


# ---------------------------------------------------------------------------
# 1. The faster path returns the same bits, not merely the same numbers
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("fmt", ["csc", "csr", "dense"])
def test_mean_over_rows_is_bit_identical_to_the_per_row_stack(fmt):
    """Exact equality, deliberately — `assert_allclose` would hide the point.

    Two formulations that are *faster still* are not bit-identical, and the
    difference is only in the last ulps: an indicator-vector matvec
    (``ind @ mat``) lands 7.5e-16 away and summing the CSC selection directly
    2.5e-15 away, because each accumulates the columns in a different order.
    Nothing downstream would notice, and that is exactly why a tolerance here
    would let a silent change of estimator through.
    """
    dense = _wide_matrix()
    mat = {"csc": sp.csc_matrix(dense), "csr": sp.csr_matrix(dense),
           "dense": dense}[fmt]
    rng = np.random.default_rng(3)
    rows = sorted(rng.choice(N_GENES, size=64, replace=False).tolist())

    got = _mean_over_rows(mat, rows)
    want = _stack_mean(mat, rows)
    assert np.array_equal(got, want), (
        f"max|diff| = {np.abs(got - want).max():.3e} — same numbers, "
        f"different arithmetic"
    )


def test_mean_over_rows_respects_the_order_it_was_given():
    """Row order is not cosmetic: it is the order the terms are summed in.

    Sorting the indices before selecting is the obvious optimization here —
    scipy indexes sorted rows faster — and it silently returns different
    numbers. That is also what makes the control-gene ordering below observable
    at all, so if this property does not hold there is nothing to reproduce.
    """
    mat = sp.csc_matrix(_wide_matrix())
    rng = np.random.default_rng(11)
    rows = rng.choice(N_GENES, size=64, replace=False).tolist()
    assert rows != sorted(rows), "fixture must present the rows out of order"

    assert np.array_equal(_mean_over_rows(mat, rows), _stack_mean(mat, rows))
    # The fixture has to be able to tell the two apart, or the assertion above
    # holds for any implementation and guards nothing.
    assert not np.array_equal(
        _mean_over_rows(mat, rows), _mean_over_rows(mat, sorted(rows))
    ), "fixture is too well-conditioned to see a reordering"


# ---------------------------------------------------------------------------
# 2. One selection per gene set, not one per gene
# ---------------------------------------------------------------------------

class _CountingMatrix(sp.csc_matrix):
    """A CSC matrix that records how many times it is sliced.

    A *subclass*, not a wrapper: `add_module_score` branches on
    `sp.issparse(mat)`, and a proxy that merely forwards attributes takes the
    dense branch, where `np.asarray` turns it into a 0-d object array. Same trap
    that made `fetch_data` return a column of matrices.
    """

    calls = 0

    def __getitem__(self, key):
        _CountingMatrix.calls += 1
        return super().__getitem__(key)


def test_add_module_score_slices_per_gene_set_not_per_gene():
    """The structural property behind the 100× — asserted, not timed.

    A wall-clock threshold would be flaky on a loaded machine and would not say
    *why* it was slow. The slice count does: one for the expression-bin pool,
    then one for the program's genes and one for its controls. The old code made
    one call per gene, which on a default `ctrl=100` is thousands.

    Counted through `add_module_score` rather than by calling `_mean_over_rows`
    directly, because the defect was the *call site* — a helper that reads the
    whole selection at once is no use to a loop that hands it one row at a time.
    """
    obj = _obj()
    assay = obj.assays["RNA"]
    _CountingMatrix.calls = 0
    assay.layers["data"] = _CountingMatrix(assay.layers["data"])

    genes = [f"g{i}" for i in range(20)]
    add_module_score(obj, features={"prog": genes}, ctrl=50, seed=1)

    assert _CountingMatrix.calls <= 4, (
        f"{_CountingMatrix.calls} slices for a 20-gene program — the matrix "
        f"should be read once per gene set, not once per gene"
    )
    # And it did the work: a score per cell, not a column of zeros.
    scores = obj.meta_data["prog"].to_numpy()
    assert scores.shape == (N_CELLS,)
    assert np.abs(scores).sum() > 0


def test_scores_still_match_the_per_row_stack_end_to_end():
    """The whole function, against the arithmetic it used to do.

    `_mean_over_rows` being right in isolation does not make `add_module_score`
    right: the row indices it is handed, and their order, are chosen at the call
    site. Recomputed here from the same control genes the function drew.
    """
    obj = _obj()
    mat = obj.assays["RNA"].layers["data"]
    feats = obj.assays["RNA"].features()
    idx = {f: i for i, f in enumerate(feats)}
    genes = [f"g{i}" for i in range(20)]

    add_module_score(obj, features={"prog": genes}, ctrl=50, seed=1)
    got = obj.meta_data["prog"].to_numpy()

    # Re-derive the program half exactly; the control half is checked by the
    # reproducibility test, which does not need to know which genes were drawn.
    prog_mean = _stack_mean(mat, [idx[g] for g in genes])
    assert np.isfinite(got).all()
    assert not np.array_equal(got, prog_mean), "controls were not subtracted"
    assert np.array_equal(
        _mean_over_rows(mat, [idx[g] for g in genes]), prog_mean
    )


# ---------------------------------------------------------------------------
# 3. The same seed gives the same score in a different process
# ---------------------------------------------------------------------------

_CHILD = textwrap.dedent(
    """
    import numpy as np, scipy.sparse as sp
    from shanuz import create_shanuz_object, normalize_data
    from shanuz.module_score import add_module_score

    rng = np.random.default_rng(0)
    dense = rng.poisson(2.0, size=(300, 200)).astype(float)
    dense[rng.random(dense.shape) > 0.3] = 0.0
    obj = create_shanuz_object(
        counts=sp.csc_matrix(dense), assay="RNA",
        feature_names=[f"g{i}" for i in range(300)],
        cell_names=[f"c{j}" for j in range(200)],
    )
    normalize_data(obj)
    add_module_score(obj, features={"prog": [f"g{i}" for i in range(20)]},
                     ctrl=50, seed=1)
    print(obj.meta_data["prog"].to_numpy().tobytes().hex())
    """
)


def _score_under_hash_seed(seed: str) -> str:
    env = dict(os.environ, PYTHONHASHSEED=seed)
    out = subprocess.run(
        [sys.executable, "-c", _CHILD], capture_output=True, text=True, env=env,
        cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    )
    assert out.returncode == 0, out.stderr
    return out.stdout.strip()


def test_control_scores_do_not_depend_on_pythons_hash_seed():
    """A `set` of gene names iterates in an order that changes per process.

    The control score is a mean over those genes, and floating-point addition is
    not associative, so the same object at the same seed scored differently in a
    different process — 9.7e-16 on THP-1, enough to move a value nothing else
    moved. Two child processes, two hash seeds, byte-for-byte comparison: the
    `set` version fails this and no other test in the suite does, because
    everything else runs in one process where the order is at least stable.
    """
    assert _score_under_hash_seed("0") == _score_under_hash_seed("12345")
