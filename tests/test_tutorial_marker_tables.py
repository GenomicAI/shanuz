"""Tests for the tutorials' marker-table presentation helpers.

These exist because the tutorials had no execution coverage at all: the full
suite passed on an install where ``tutorials/pbmc3k_tutorial.py`` crashed
outright. The other tutorial tests import annotation *helpers* and run them on
synthetic data, so nothing here reached the code that formats marker tables.

The bug: ``all_markers.groupby("cluster").apply(...)`` returned a frame with no
``cluster`` column under pandas 3 (the grouping column stopped being passed into
the callable), and the next line filtered on it -> ``KeyError: 'cluster'``.

Note the version asymmetry, because it explains why this shipped: on pandas 2
the old code *works*. No unit test can fail on pandas 2 for a bug that only
exists on pandas 3. What these guard is every pandas-3 environment — which is
CI and, since ``pyproject`` declares ``pandas>=2.0``, every fresh install.
"""
import importlib
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tutorials.generate_advanced_plots import _top_genes  # noqa: E402
from tutorials.pbmc3k_tutorial import top_markers_per_cluster  # noqa: E402
from tutorials.pbmc8k_subclustering_tutorial import top_markers_table  # noqa: E402

TUTORIAL_MODULES = [
    "tutorials.pbmc3k_tutorial",
    "tutorials.pbmc8k_subclustering_tutorial",
    "tutorials.cbmc_citeseq_tutorial",
    "tutorials.pbmc3k_sctransform_tutorial",
    "tutorials.generate_plots",
    "tutorials.generate_advanced_plots",
    "tutorials.generate_multimodal_plots",
    "tutorials.generate_sctransform_plots",
    "tutorials.generate_spatial_plots",
]


def _markers(n_clusters=12, per_cluster=4):
    """A find_all_markers-shaped frame.

    12 clusters on purpose: cluster labels are strings, so anything sorting them
    lexicographically orders 0, 1, 10, 11, 2 ... and a 9-cluster fixture (pbmc3k's
    real shape) would hide it.

    Within each cluster, p_val ranks genes in the reverse order of avg_log2FC, so
    a helper that sorts by the wrong column picks visibly wrong genes rather than
    coincidentally right ones.
    """
    rows = []
    for c in range(n_clusters):
        for i in range(per_cluster):
            rows.append({
                "cluster": str(c),
                "gene": f"g{c}_{i}",
                "p_val": 10.0 ** (-(per_cluster - i) - c),
                "avg_log2FC": float(i + 1),
                "pct.1": 0.9,
                "pct.2": 0.1,
                "p_val_adj": 10.0 ** (-(per_cluster - i) - c),
            })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# pbmc3k — the helper that replaced the crashing groupby.apply
# ----------------------------------------------------------------------


def test_top_markers_per_cluster_keys_every_cluster():
    """The regression itself: the result must carry cluster identity.

    The old code dropped the ``cluster`` column on pandas 3 and raised KeyError
    one line later. Anything that reintroduces ``groupby(...).apply(...)`` here
    fails this on pandas 3.
    """
    out = top_markers_per_cluster(_markers(), n=3)
    assert set(out) == {str(c) for c in range(12)}
    assert all(len(genes) == 3 for genes in out.values())


def test_top_markers_per_cluster_picks_the_most_significant_genes():
    """Smallest p_val wins — not largest, and not avg_log2FC.

    The fixture ranks p_val opposite to avg_log2FC, so a helper that sorted by
    the wrong column would return the exact reverse of this (``g0_3, g0_2``) —
    which is what ``_top_genes`` correctly returns for its own log2FC ranking.
    """
    out = top_markers_per_cluster(_markers(), n=2)
    assert out["0"] == ["g0_0", "g0_1"]
    assert out["11"] == ["g11_0", "g11_1"]


def test_top_markers_per_cluster_orders_clusters_numerically():
    """0, 1, 2 ... 11 — not the lexicographic 0, 1, 10, 11, 2 ...

    The helper sorts with ``key=int``; a plain ``sorted()`` on these string
    labels would pass every other assertion here while printing the table in a
    confusing order.
    """
    assert list(top_markers_per_cluster(_markers(), n=1)) == [str(c) for c in range(12)]


def test_top_markers_per_cluster_handles_a_cluster_with_fewer_than_n_markers():
    """Real marker frames are ragged — a sparse cluster must not raise."""
    df = _markers(n_clusters=3)
    df = df[~((df["cluster"] == "1") & (df["gene"] != "g1_0"))]
    out = top_markers_per_cluster(df, n=3)
    assert out["1"] == ["g1_0"]
    assert len(out["0"]) == 3


# ----------------------------------------------------------------------
# pbmc8k — the helper pbmc3k's fix was modelled on (was already correct)
# ----------------------------------------------------------------------


def test_top_markers_table_renders_every_cluster():
    text = top_markers_table(_markers(), n=2)
    lines = text.splitlines()
    assert len(lines) == 12
    assert lines[0].startswith("    Cluster 0: ")
    # numeric ordering, same as pbmc3k's helper
    assert lines[-1].startswith("    Cluster 11: ")


# ----------------------------------------------------------------------
# plot generators — never broken, but one "simplification" away from wrong
# ----------------------------------------------------------------------


def test_top_genes_is_cluster_major():
    """Gene order is load-bearing: it sets the heatmap's row order.

    This is the guard that matters here. ``_top_genes`` was *not* broken by
    pandas 3 — the old ``groupby(...).apply(...)`` still returned the right
    genes. The hazard is the obvious-looking rewrite,
    ``sort_values("avg_log2FC", ascending=False).groupby(...).head(n)``, which
    returns the same *set* of genes interleaved across clusters, silently
    turning the DoHeatmap's per-cluster blocks into noise. Verified: that
    rewrite fails this test, on both pandas 2 and 3.
    """
    genes = _top_genes(_markers(n_clusters=3), n=2)
    # groupby sorts string labels lexicographically; with 3 clusters that is also
    # numeric order. Each cluster's genes must be contiguous and log2FC-ranked.
    assert genes == ["g0_3", "g0_2", "g1_3", "g1_2", "g2_3", "g2_2"]


def test_top_genes_deduplicates_but_keeps_first_position():
    """A gene topping two clusters appears once, at its first occurrence."""
    df = _markers(n_clusters=2, per_cluster=2)
    df.loc[df["gene"] == "g1_1", "gene"] = "g0_1"  # now shared by both clusters
    genes = _top_genes(df, n=2)
    assert genes.count("g0_1") == 1
    assert genes.index("g0_1") < genes.index("g1_0")


# ----------------------------------------------------------------------
# every tutorial module at least imports
# ----------------------------------------------------------------------


@pytest.mark.parametrize("module", TUTORIAL_MODULES)
def test_tutorial_module_imports(module):
    """Cheap, but it is more than the suite had before.

    Would not have caught the pandas 3 crash (that was mid-function), so it is a
    floor, not a substitute for running the tutorials. The end-to-end runs live
    in test_tutorial_smoke.py and need the cached datasets.
    """
    assert importlib.import_module(module) is not None
