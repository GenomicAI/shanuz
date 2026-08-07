"""``find_clusters`` at one resolution and at several, against Seurat's conventions.

Seurat's ``FindClusters`` writes each resolution to its own metadata column named
``{graph}_res.{resolution}`` and leaves the object sitting on the **last**
resolution given. Truecell wrote only ``seurat_clusters`` and accepted only a
scalar, so a ported script reading ``obj[["RNA_snn_res.0.5"]]`` raised
``KeyError`` — at the default settings, with nothing in the suite noticing.

Every expected value here was pinned against R Seurat 5.5.1 rather than derived
from the Python side, because what is being tested is a naming and ordering
*convention*, which is exactly the class of thing a Python-side re-derivation
cannot check. See ``REVIEW_PLAN.md``.
"""
import numpy as np
import pandas as pd
import pytest

from truecell import (
    create_truecell_object,
    find_clusters,
    find_neighbors,
    find_variable_features,
    normalize_data,
    run_pca,
    scale_data,
)
from truecell.clustering import _res_label


@pytest.fixture
def clustered(small_counts, feature_names, cell_names):
    """A small object carried as far as the SNN graph."""
    import scipy.sparse as sp
    rng = np.random.default_rng(1)
    counts = sp.csr_matrix(rng.poisson(2.0, size=(200, 300)).astype(float))
    obj = create_truecell_object(
        counts=counts,
        feature_names=[f"gene_{i}" for i in range(200)],
        cell_names=[f"cell_{i}" for i in range(300)],
    )
    normalize_data(obj)
    find_variable_features(obj)
    scale_data(obj)
    run_pca(obj, n_pcs=10)
    find_neighbors(obj, dims=list(range(1, 6)))
    return obj


# ---------------------------------------------------------------------------
# The column label is R's number formatting, not Python's
# ---------------------------------------------------------------------------

# Pinned against R: `as.character(x)` for each value, R 4.x at the default
# `scipen = 0`. The pair that matters most in practice is 1.0 -> "1": Python's
# own `str(1.0)` gives "1.0", which would name the column `RNA_snn_res.1.0` and
# leave a ported script reading `RNA_snn_res.1` with a KeyError.
@pytest.mark.parametrize("value,expected", [
    (0.4, "0.4"),
    (1.0, "1"),
    (0.50, "0.5"),
    (1.25, "1.25"),
    (2, "2"),
    (0.125, "0.125"),
    (1e-2, "0.01"),
    (10, "10"),
    (0.8, "0.8"),
    (1234.5, "1234.5"),
    (123456, "123456"),
    (1234567, "1234567"),
    (0.1, "0.1"),
    (0.02, "0.02"),
    # R switches to scientific when it is strictly shorter, ties going to fixed:
    # "0.001" and "1e-03" are both 5 characters, so 0.001 stays fixed, while
    # "0.0001" (6) loses to "1e-04" (5).
    (0.001, "0.001"),
    (0.0001, "1e-04"),
    (0.00012, "0.00012"),
    (1e-5, "1e-05"),
    (1e-7, "1e-07"),
    (1e5, "1e+05"),
    (1e6, "1e+06"),
])
def test_resolution_label_matches_r_as_character(value, expected):
    assert _res_label(value) == expected


# ---------------------------------------------------------------------------
# The column exists at all — the gap that predates the multi-resolution work
# ---------------------------------------------------------------------------

def test_a_single_resolution_still_writes_its_own_column(clustered):
    find_clusters(clustered, resolution=0.5)
    assert "RNA_snn_res.0.5" in clustered.meta_data.columns
    assert (clustered.meta_data["RNA_snn_res.0.5"].astype(str)
            == clustered.meta_data["seurat_clusters"].astype(str)).all()


def test_a_whole_number_resolution_drops_its_trailing_zero(clustered):
    """The regression the label helper exists for, exercised end to end."""
    find_clusters(clustered, resolution=1.0)
    assert "RNA_snn_res.1" in clustered.meta_data.columns
    assert "RNA_snn_res.1.0" not in clustered.meta_data.columns


def test_the_column_name_follows_the_graph_it_clustered(clustered):
    find_clusters(clustered, resolution=0.5, graph_name="RNA_nn")
    assert "RNA_nn_res.0.5" in clustered.meta_data.columns


# ---------------------------------------------------------------------------
# Several resolutions at once
# ---------------------------------------------------------------------------

def test_each_resolution_gets_its_own_column(clustered):
    find_clusters(clustered, resolution=[0.4, 0.8, 1.2])
    for name in ("RNA_snn_res.0.4", "RNA_snn_res.0.8", "RNA_snn_res.1.2"):
        assert name in clustered.meta_data.columns


def test_the_object_is_left_on_the_last_resolution_given(clustered):
    """Last as *given*, not largest — Seurat takes the last column of its frame."""
    find_clusters(clustered, resolution=[1.2, 0.8, 0.4])
    last = clustered.meta_data["RNA_snn_res.0.4"].astype(str)
    assert (clustered.meta_data["seurat_clusters"].astype(str) == last).all()
    assert (pd.Series(clustered.idents).astype(str).values == last.values).all()
    # ...and specifically not the largest, which is the plausible wrong answer.
    largest = clustered.meta_data["RNA_snn_res.1.2"].astype(str)
    if not (largest == last).all():
        assert not (clustered.meta_data["seurat_clusters"].astype(str)
                    == largest).all()


def test_a_partition_does_not_depend_on_what_preceded_it(clustered, small_counts,
                                                         feature_names, cell_names):
    """Verified against Seurat 5.5.1, which reseeds rather than running the stream on.

    If the RNG ran on between resolutions, 0.8's partition would differ between
    these three calls, and a user comparing resolutions would be comparing
    partitions that also differ by seed.
    """
    find_clusters(clustered, resolution=[0.4, 0.8, 1.2])
    forward = clustered.meta_data["RNA_snn_res.0.8"].astype(str).values.copy()

    find_clusters(clustered, resolution=[1.2, 0.8, 0.4])
    reverse = clustered.meta_data["RNA_snn_res.0.8"].astype(str).values.copy()

    find_clusters(clustered, resolution=0.8)
    alone = clustered.meta_data["RNA_snn_res.0.8"].astype(str).values

    assert (forward == reverse).all()
    assert (forward == alone).all()


def test_resolutions_are_not_all_the_same_partition(clustered):
    """Guards the loop actually using its own `res` rather than a fixed one.

    Without this, passing the first resolution to every call would satisfy every
    other test in this file: the columns would all exist and all be consistent.
    """
    find_clusters(clustered, resolution=[0.01, 2.0])
    coarse = clustered.meta_data["RNA_snn_res.0.01"].astype(str)
    fine = clustered.meta_data["RNA_snn_res.2"].astype(str)
    assert coarse.nunique() < fine.nunique()


# ---------------------------------------------------------------------------
# cluster_name
# ---------------------------------------------------------------------------

def test_cluster_name_replaces_the_generated_names(clustered):
    find_clusters(clustered, resolution=[0.4, 0.8], cluster_name=["lo", "hi"])
    assert {"lo", "hi"} <= set(clustered.meta_data.columns)
    assert "RNA_snn_res.0.4" not in clustered.meta_data.columns
    # Seurat writes seurat_clusters either way.
    assert "seurat_clusters" in clustered.meta_data.columns


def test_a_single_cluster_name_is_accepted_as_a_string(clustered):
    find_clusters(clustered, resolution=0.5, cluster_name="my_clusters")
    assert "my_clusters" in clustered.meta_data.columns


def test_a_mismatched_cluster_name_count_is_rejected(clustered):
    with pytest.raises(ValueError, match="one per resolution"):
        find_clusters(clustered, resolution=[0.4, 0.8], cluster_name=["only_one"])


# ---------------------------------------------------------------------------
# Validation happens before anything is written
# ---------------------------------------------------------------------------

def test_an_empty_resolution_list_is_rejected(clustered):
    with pytest.raises(ValueError, match="at least one"):
        find_clusters(clustered, resolution=[])


@pytest.mark.parametrize("algorithm,exc", [(3, NotImplementedError),
                                           (9, ValueError)])
def test_a_bad_algorithm_writes_no_columns_at_all(clustered, algorithm, exc):
    """A rejected `algorithm` leaves the object exactly as it found it.

    Note what this does *not* establish. Moving the check from above the loop to
    inside it passes this test unchanged — mutation testing confirmed it — because
    the dispatch is the first statement of the loop body, so the first resolution
    raises before anything is stored. The eager check is tidiness, not a guard
    against a partial write, and the comment in ``clustering.py`` says so.
    """
    before = set(clustered.meta_data.columns)
    with pytest.raises(exc):
        find_clusters(clustered, resolution=[0.4, 0.8], algorithm=algorithm)
    assert set(clustered.meta_data.columns) == before
