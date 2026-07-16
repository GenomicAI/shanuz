"""Tests for cell-hashing demultiplexing (v0.9.0 Specialized Assays).

  * hto_demux  (hto.py)

A synthetic hashtag matrix is built with *known* ground truth — clean singlets
for each of four hashtags, doublets spanning tag pairs, and empty-ish negatives —
and ``hto_demux`` is asked to recover it. The Seurat metadata columns, the
``hash.ID`` identity, and the singlet/doublet/negative accuracy are all checked
against that ground truth. Network-free and deterministic (fixed RNG + seed).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz._clara import clara, clara_sampsize  # noqa: E402
from shanuz.hto import hto_demux  # noqa: E402
from shanuz.preprocessing import normalize_data  # noqa: E402
from shanuz.shanuz import create_shanuz_object  # noqa: E402

TAGS = ["HTO-A", "HTO-B", "HTO-C", "HTO-D"]
DOUBLET_PAIRS = [(0, 1), (1, 2), (2, 3), (0, 2), (1, 3)]


def _hashing_counts(seed=0):
    """A hashtag count matrix (4 × N) plus per-cell ground-truth labels."""
    rng = np.random.default_rng(seed)
    H = len(TAGS)
    cols, truth = [], []

    # Background staining is ~uniform per antibody across all droplets; a cell is
    # "positive" for a tag only when its on-target signal towers over that floor.
    bg = 1.0
    for h in range(H):                       # 30 singlets per hashtag
        for _ in range(30):
            col = rng.poisson(bg, size=H).astype(float)
            col[h] = rng.poisson(100.0)
            cols.append(col)
            truth.append(TAGS[h])

    for _ in range(4):                       # 20 doublets across tag pairs
        for a, b in DOUBLET_PAIRS:
            col = rng.poisson(bg, size=H).astype(float)
            col[a] = rng.poisson(70.0)
            col[b] = rng.poisson(70.0)
            cols.append(col)
            truth.append("Doublet")

    for _ in range(15):                      # 15 negatives (only background)
        cols.append(rng.poisson(bg, size=H).astype(float))
        truth.append("Negative")

    counts = np.asarray(cols).T              # H × N
    return counts, np.asarray(truth)


def _hashing_object(seed=0):
    counts, truth = _hashing_counts(seed)
    cells = [f"c{i}" for i in range(counts.shape[1])]
    obj = create_shanuz_object(
        counts=sp.csc_matrix(counts), assay="HTO",
        feature_names=TAGS, cell_names=cells,
        meta_data=pd.DataFrame(index=cells),
    )
    return obj, truth


# ----------------------------------------------------------------------
# metadata columns & identity
# ----------------------------------------------------------------------


def test_writes_seurat_metadata_columns():
    obj, _ = _hashing_object()
    hto_demux(obj)
    for col in (
        "HTO_maxID", "HTO_secondID", "HTO_margin",
        "HTO_classification", "HTO_classification.global", "hash.ID",
    ):
        assert col in obj.meta_data.columns


def test_global_classification_domain():
    obj, _ = _hashing_object()
    hto_demux(obj)
    vals = set(obj.meta_data["HTO_classification.global"])
    assert vals <= {"Negative", "Singlet", "Doublet"}


def test_idents_set_to_hash_id():
    obj, _ = _hashing_object()
    hto_demux(obj)
    assert len(obj.idents) == len(obj.cell_names())
    assert list(obj.idents) == list(obj.meta_data["hash.ID"])


def test_margin_nonnegative_and_ids_distinct():
    obj, _ = _hashing_object()
    hto_demux(obj)
    margin = obj.meta_data["HTO_margin"].to_numpy()
    assert (margin >= -1e-9).all()
    max_id = obj.meta_data["HTO_maxID"].to_numpy()
    second_id = obj.meta_data["HTO_secondID"].to_numpy()
    assert (max_id != second_id).all()


# ----------------------------------------------------------------------
# recovery against ground truth
# ----------------------------------------------------------------------


def test_singlets_recovered():
    obj, truth = _hashing_object()
    hto_demux(obj)
    gclass = obj.meta_data["HTO_classification.global"].to_numpy()
    hash_id = obj.meta_data["hash.ID"].to_numpy()

    singlet = np.isin(truth, TAGS)
    assert (gclass[singlet] == "Singlet").mean() >= 0.9
    # correct singlets are assigned to their true hashtag
    assert (hash_id[singlet] == truth[singlet]).mean() >= 0.9


def test_doublets_recovered():
    obj, truth = _hashing_object()
    hto_demux(obj)
    gclass = obj.meta_data["HTO_classification.global"].to_numpy()
    doublet = truth == "Doublet"
    assert (gclass[doublet] == "Doublet").mean() >= 0.75
    # a called doublet's classification names its two tags, sorted + joined
    hto_class = obj.meta_data["HTO_classification"].to_numpy()
    for lab in hto_class[doublet & (gclass == "Doublet")]:
        assert "_" in lab and lab == "_".join(sorted(lab.split("_")))


def test_negatives_recovered():
    obj, truth = _hashing_object()
    hto_demux(obj)
    gclass = obj.meta_data["HTO_classification.global"].to_numpy()
    negative = truth == "Negative"
    assert (gclass[negative] == "Negative").mean() >= 0.8


# ----------------------------------------------------------------------
# bookkeeping, options & guards
# ----------------------------------------------------------------------


def test_cutoffs_stored_in_misc():
    obj, _ = _hashing_object()
    hto_demux(obj)
    info = obj.misc["hto_demux"]["HTO"]
    assert set(info["cutoffs"]) == set(TAGS)
    assert info["ncenters"] == len(TAGS) + 1
    assert all(np.isfinite(c) for c in info["cutoffs"].values())


def test_init_overrides_center_count():
    obj, _ = _hashing_object()
    hto_demux(obj, init=6)
    assert obj.misc["hto_demux"]["HTO"]["ncenters"] == 6


def test_normalize_false_uses_data_layer():
    obj, truth = _hashing_object()
    normalize_data(obj, normalization_method="CLR", margin=1, assay="HTO")
    hto_demux(obj, normalize=False)
    gclass = obj.meta_data["HTO_classification.global"].to_numpy()
    singlet = np.isin(truth, TAGS)
    assert (gclass[singlet] == "Singlet").mean() >= 0.9


def test_deterministic():
    obj_a, _ = _hashing_object()
    obj_b, _ = _hashing_object()
    hto_demux(obj_a, seed=7)
    hto_demux(obj_b, seed=7)
    assert list(obj_a.meta_data["hash.ID"]) == list(obj_b.meta_data["hash.ID"])


def test_requires_two_hashtags():
    rng = np.random.default_rng(0)
    counts = sp.csc_matrix(rng.poisson(5.0, size=(1, 20)).astype(float))
    cells = [f"c{i}" for i in range(20)]
    obj = create_shanuz_object(
        counts=counts, assay="HTO",
        feature_names=["HTO-A"], cell_names=cells,
        meta_data=pd.DataFrame(index=cells),
    )
    with pytest.raises(ValueError):
        hto_demux(obj)


def test_unsupported_kfunc_raises():
    obj, _ = _hashing_object()
    with pytest.raises(NotImplementedError):
        hto_demux(obj, kfunc="pam")


# ----------------------------------------------------------------------
# clara (k-medoids) -- Seurat's default kfunc
# ----------------------------------------------------------------------

# Ground truth captured from R, not re-derived in Python:
#   library(cluster); clara(x, k = <k>, samples = <s>, sampsize = <ss>)$clustering
# cluster 2.1.8.2, at the defaults Seurat's HTODemux relies on (euclidean,
# rngR = FALSE, pamLike = FALSE, stand = FALSE). The four cases cover the three
# branches of clara's sampling loop plus the degenerate one:
#   FULL   n == sampsize  -> one sample is the whole dataset
#   LRG    n <  2*sampsize -> clara draws the *complement* and inverts it
#   NORMAL n >  2*sampsize -> the usual medoid-seeded resampling
#   K1     k == 1          -> BUILD only, SWAP never runs
#
# Every value is a multiple of 1/4, so the squared-distance sums are exact and
# the fixtures cannot drift with the floating-point rounding discussed in
# shanuz/_clara.py. Verified identical from clara.c built for arm64 and x86_64.
_R_FULL_X = np.array([[2,1,1.5], [0.75,0.25,1.75], [1.5,2.25,1.25], [0,2,2.25], [0.25,2.75,1.5], [3,0,0.5], [1.5,0.75,2.25], [2.5,0.5,1.25], [0.25,1.25,1.75], [2.5,2.25,0.25], [0.5,2.25,3], [0,1.25,0.25], [1,0.75,1.25], [1,2.75,1.25], [2.25,0.75,0], [1.25,2.25,0.5], [2.25,2.75,0.5], [1.5,2,1.75], [2,1.5,1.25], [1,1.25,1.5], [1,2,2.75], [2,1.75,1.25], [2,2.75,1.75], [1,2,1.5]])
_R_FULL_CL = np.array([1, 1, 2, 2, 2, 1, 1, 1, 2, 3, 2, 2, 1, 2, 1, 2, 3, 2, 1, 2, 2, 1, 2, 2])
_R_LRG_X = np.array([[1,0.5,1.25], [1.25,2,2.75], [1.25,1.5,0.75], [1.75,1.75,2.75], [0,3,0.75], [3,3,2], [0,1.25,2.25], [2.75,0.25,0.25], [3,1.5,1.25], [2,0.25,1.25], [0.25,0.5,0.5], [2.5,0.75,3], [0,0.5,0.5], [0.5,0,1.75], [1.25,1.5,1.25], [0.25,2,3], [0.5,3,3], [1.5,0,0], [1.75,0.25,1], [1.5,1.75,1.25], [0,0.75,0], [1.25,3,2.25], [2,1,2], [0.75,2.75,1], [2.5,0.25,1.75], [1.25,1,0.75], [2,1.25,0], [1.75,3,0], [3,1.5,1], [1.25,0.25,2.5]])
_R_LRG_CL = np.array([1, 2, 3, 2, 3, 4, 2, 1, 4, 1, 3, 2, 3, 1, 3, 2, 2, 1, 1, 3, 3, 2, 1, 3, 1, 3, 3, 3, 4, 1])
_R_NORMAL_X = np.array([[1,1,3,2.25], [2.25,2.5,0.75,0], [2.75,2.25,1.5,2], [1.5,0.75,0.25,1.25], [0.75,2,0.75,0.5], [2.25,1,1.75,2], [1.75,1.25,1.75,1.5], [2.5,1.5,0.75,0.25], [1.75,1.5,0.75,0.5], [0.75,0.5,2.5,2.75], [2.25,1.75,2.25,3], [1.5,2.75,2.5,2.25], [1.75,0.75,3,1.5], [1.75,2.5,2,2.5], [1.75,3,3,1], [1,1,0.5,1.5], [0.25,2.25,2,2.25], [3,1.75,1.25,0.5], [2.75,1.5,0.5,3], [3,2.75,0.25,1.25], [1,2,1,3], [1.75,2.5,2.25,1.25], [1,0.25,0.75,1], [2,2,3,1], [2,2.5,0.25,0.25], [1.75,2.5,0.25,0], [1.25,1,2.5,2.25], [0.25,2,0.25,0.25], [2,2.25,0,1.75], [2.75,2.25,3,1.5], [1.75,1.5,2.25,0.75], [1.25,0,1.75,2.75], [2.75,3,2.25,0.5], [0.5,2.75,1.75,0], [1.25,2.25,0.25,0], [1.75,1.75,1.5,0], [2.25,1.75,2,0], [2.5,3,1.25,2.5], [3,2,2.75,3], [1.5,0.25,1.25,1.75], [0.5,0.25,2.25,2.25], [1,1.5,0,0.25], [1.5,3,1.75,2.75], [1.25,3,2.25,0.75], [1.75,0.5,1.5,1.5], [2.25,2.75,0.5,1.25], [0.75,3,2,2.25], [0,2.75,0,2.75], [0.5,2.75,0.75,2.75], [1.25,0,2.75,1], [0.75,1.25,1.25,0.5], [2,1.25,2.75,3], [2.25,3,0.75,1.5], [3,2.75,1,3], [1.5,1.75,1.25,2], [2.25,0,0.5,3], [0,2,0.5,3], [2,2.25,1.25,0], [0,1,1.5,1.25], [2,1.5,0,0], [1,3,1,1.75], [0.25,3,1.5,0.5], [1.25,0.25,0.5,2.75], [3,1,2,0.25], [1.75,0.5,0.25,1], [2.5,2.75,0.25,3], [2.25,0.75,2.5,1.5], [0.75,0.5,0.5,1.5], [2,2.25,1,1], [0.5,0.5,2.25,1.5], [0.75,3,2.25,2.75], [2.75,2,1.5,0.5], [0.25,0,1,0.5], [0,0.75,2.5,0.5], [0.75,2.25,2.5,1.25], [1.25,1.75,2.25,2.75], [0,0.25,2,1], [3,0.5,3,1.75], [0.25,2,2,0.5], [1,1.25,1.75,2], [0.25,0.75,3,1.5], [1.75,2.25,1,2.25], [1.5,2.75,0.75,1.5], [1.5,0,1.5,0.25], [2.25,1.5,1.25,0], [1.25,0.75,0,3], [0.5,3,0,2.75], [0,0.5,0.25,3], [2.25,1,0.5,1.75], [0.25,1.25,0.5,1.25]])
_R_NORMAL_CL = np.array([1, 2, 3, 4, 2, 1, 1, 2, 2, 1, 1, 3, 1, 3, 2, 4, 1, 2, 3, 2, 3, 2, 4, 1, 2, 2, 1, 4, 3, 1, 2, 1, 2, 2, 2, 2, 2, 3, 1, 4, 1, 4, 3, 2, 4, 2, 3, 3, 3, 1, 4, 1, 2, 3, 3, 4, 3, 2, 4, 2, 3, 2, 4, 2, 4, 3, 1, 4, 2, 1, 3, 2, 4, 1, 1, 1, 1, 1, 2, 1, 1, 3, 2, 4, 2, 4, 3, 4, 4, 4])
_R_K1_X = np.array([[1.75,3], [2.5,1], [0.5,3], [0.5,0.25], [1.5,2.75], [2.75,1], [0.5,3], [1.25,1.25], [1,1], [2.25,0], [0.5,1.25], [1.75,2.5], [1.25,1], [0.25,0.75], [1.75,0], [0,1.75], [0.75,2.75], [2.5,1.25], [0.5,2.5], [2.5,3]])
_R_K1_CL = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

_R_CLARA_CASES = [
    ("FULL", _R_FULL_X, _R_FULL_CL, 3, 5, 24),
    ("LRG", _R_LRG_X, _R_LRG_CL, 4, 10, 20),
    ("NORMAL", _R_NORMAL_X, _R_NORMAL_CL, 4, 20, 30),
    ("K1", _R_K1_X, _R_K1_CL, 1, 5, 12),
]


@pytest.mark.parametrize(
    "name,x,expected,k,samples,sampsize",
    _R_CLARA_CASES,
    ids=[c[0] for c in _R_CLARA_CASES],
)
def test_clara_matches_r_ground_truth(name, x, expected, k, samples, sampsize):
    """clara must reproduce R's cluster assignments, not merely cluster sensibly.

    The whole point of porting clara rather than reaching for a k-medoids library
    is that HTODemux's cutoffs depend on which cells land in each cluster, and
    clara's answer is decided by details no library shares: its own 16-bit RNG,
    a swap rule that is not pam()'s, and tie-breaks that go opposite ways in
    BUILD and SWAP. Only fixed output from a real cluster::clara run pins those.
    """
    got = clara(x, k=k, samples=samples, sampsize=sampsize) + 1  # R is 1-based
    assert np.array_equal(got, expected)


def test_clara_takes_no_seed_and_is_deterministic():
    """clara's sampling is not seedable -- in R either.

    ``HTODemux`` calls ``set.seed(seed)`` before clustering, but ``clara`` at
    ``rngR = FALSE`` (its default, and Seurat's) draws from a built-in generator
    that ``set.seed`` cannot reach, so that call is a no-op on this path. Anyone
    who "fixes" _clara by threading a seed through it will silently stop matching R.
    """
    x = _R_NORMAL_X
    a = clara(x, k=4, samples=20, sampsize=30)
    b = clara(x, k=4, samples=20, sampsize=30)
    assert np.array_equal(a, b)
    import inspect

    assert "seed" not in inspect.signature(clara).parameters


def test_clara_full_sample_ignores_sample_count():
    """When sampsize == n every draw is the whole dataset, so samples cannot matter.

    C short-circuits this with `if (full_sample) break;`. Worth pinning because a
    port that kept resampling would still *pass* the ground-truth tests -- it would
    just redo identical work.
    """
    x = _R_FULL_X
    one = clara(x, k=3, samples=1, sampsize=len(x))
    many = clara(x, k=3, samples=500, sampsize=len(x))
    assert np.array_equal(one, many)


def test_clara_default_sampsize_matches_r():
    """R's default is sampsize = min(n, 40 + 2 * k)."""
    assert clara_sampsize(1000, 9) == 58
    assert clara_sampsize(30, 9) == 30  # n smaller than the formula -> full sample
    x = _R_NORMAL_X
    assert np.array_equal(
        clara(x, k=4, samples=20),
        clara(x, k=4, samples=20, sampsize=clara_sampsize(len(x), 4)),
    )


@pytest.mark.parametrize(
    "kwargs,match",
    [
        (dict(k=0), "at least 1"),
        (dict(k=24), "at most n-1"),          # n == 24 here, so k must be <= 23
        (dict(k=3, sampsize=3), "at least 4"),  # needs max(2, k+1)
        (dict(k=3, sampsize=99), "larger than"),
        (dict(k=3, samples=0), "at least 1"),
    ],
)
def test_clara_rejects_bad_arguments(kwargs, match):
    with pytest.raises(ValueError, match=match):
        clara(_R_FULL_X, **kwargs)


def test_clara_rejects_non_finite():
    """R's clara has a whole missing-data code path that this port does not have;
    refuse rather than quietly disagree with it."""
    x = _R_FULL_X.copy()
    x[0, 0] = np.nan
    with pytest.raises(ValueError, match="missing or non-finite"):
        clara(x, k=3, samples=5, sampsize=24)


def test_hto_demux_clara_recovers_the_truth():
    """The clara path must demultiplex as well as the k-means one."""
    obj, truth = _hashing_object()
    hto_demux(obj, kfunc="clara")
    gclass = obj.meta_data["HTO_classification.global"].to_numpy()
    hash_id = obj.meta_data["hash.ID"].to_numpy()
    singlet = np.isin(truth, TAGS)
    assert (gclass[singlet] == "Singlet").mean() >= 0.9
    assert (hash_id[singlet] == truth[singlet]).mean() >= 0.9
    assert (gclass[truth == "Doublet"] == "Doublet").mean() >= 0.8


def test_hto_demux_clara_and_kmeans_agree():
    """The two kfuncs pick the same background clusters, so the calls agree.

    Only each tag's *least-expressing* cluster feeds the negative-binomial fit, and
    that is robust to where exactly the two algorithms draw their boundaries.
    """
    a, _ = _hashing_object()
    b, _ = _hashing_object()
    hto_demux(a, kfunc="kmeans")
    hto_demux(b, kfunc="clara")
    assert list(a.meta_data["hash.ID"]) == list(b.meta_data["hash.ID"])


def test_hto_demux_clara_stores_cutoffs():
    obj, _ = _hashing_object()
    hto_demux(obj, kfunc="clara", nsamples=20)
    info = obj.misc["hto_demux"]["HTO"]
    assert set(info["cutoffs"]) == set(TAGS)
    assert info["ncenters"] == len(TAGS) + 1


# ----------------------------------------------------------------------
# clara's rounding sensitivity
# ----------------------------------------------------------------------

# Unlike the fixtures above, this one is NOT pinned to the R on any particular
# machine, because for an input like this there is no machine-independent R
# answer to pin to. clara accepts a swap on any improvement below zero -- R really
# does take swaps worth -2.2e-16 -- so a one-ulp change in a single distance can
# flip a swap, hence the winning sub-sample, hence the whole clustering. That
# makes the arithmetic's *order* load-bearing, and order is a compiler's choice:
# clara.c built for arm64 contracts `clk += d*d` into a fused multiply-add (one
# rounding) while the same source built for baseline x86_64, which has no FMA in
# its ISA, rounds twice. The two builds return materially different clusterings
# for this input, so R disagrees with R across architectures.
#
# shanuz follows plain IEEE double arithmetic, which is what numpy gives on every
# platform and what clara.c gives on x86_64. The value below was checked against
# clara.c compiled for x86_64. Its job is to fail if someone "simplifies" the
# summation order -- swapping the cumsum in _selec for np.sum is enough to break
# it -- not to assert that this is the One True Answer.
_IEEE_X = np.array([
    0.3, -1.2, 2.3, 1.2, -0.3, 0.1, 0.1, -0.4, -1.6, -0.6, -1.5, 0.2, -1.4,
    -1.4, 0.6, -0.9, -0.4, 0.3, -1.9, -0.1, -1, -0, 1.7, -0.8, 1.4, -0.9, 0.2,
    -0.7, -1, 1.8, -0.5, -1.6, -1, -0.4, 0.4, -0.5, -1.6, -0.2, -0.1, 0.4,
    -1.2, -1.5, -0.9, -0.1, 1.2, 0.1, -1.1, -0.2, -2, 0.9, 0.2, -0.7, -2.9,
    -0.8, 0.3, 1.4, -0.4, 0.3, 0, -1.7, 0, 0.4, 1.1, -0.3, 0, -0, -1.2, 1.2,
    1.1, 1.5, 0.4, 0.2, -0.5, -0.1, 1.1, -0.5, -0.4, -1.2, 0.2, -0.2, -2.6,
    1.3, -1.7, 0.6, -0.7, -2.2, 0.9, 0.1, 0.9, 3, -0.9, -0.3, 0, -0.5, -1.8,
    -0.6, 1.3, -0.1, 0.6, 0.2, -1, -1.4
])[:, None]
_IEEE_CL = np.array([
    1, 2, 3, 4, 5, 6, 6, 5, 7, 8, 7, 1, 2, 2, 1, 9, 5, 1, 7, 6, 9, 6, 4, 9, 4,
    9, 1, 8, 9, 4, 8, 7, 9, 5, 1, 8, 7, 5, 6, 1, 2, 7, 9, 6, 4, 6, 2, 5, 7, 4,
    1, 8, 10, 9, 1, 4, 5, 1, 6, 7, 6, 1, 4, 5, 6, 6, 2, 4, 4, 4, 1, 1, 8, 6,
    4, 8, 5, 2, 1, 5, 10, 4, 7, 1, 8, 10, 4, 6, 4, 3, 9, 5, 6, 8, 7, 8, 4, 6,
    1, 1, 9, 2
])


def test_clara_summation_order_is_load_bearing():
    """Pins the IEEE reference on a rounding-sensitive input.

    See the comment above: on inputs like this clara is chaotic, and shanuz
    deliberately tracks plain IEEE arithmetic (== clara.c on x86_64) rather than
    drifting with whatever summation order numpy finds convenient.
    """
    got = clara(_IEEE_X, k=10, samples=50, sampsize=60) + 1
    assert np.array_equal(got, _IEEE_CL)
