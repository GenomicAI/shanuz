"""Panel layout in multi-panel plots.

Ridgeline group labels are drawn outside the axes, so a panel is wider than its
axes by the length of the longest label. When the figure is too small for that,
matplotlib draws the panels over each other rather than complaining, and one
panel's labels end up printed across its neighbour.

The computed default leaves room; an explicit `figsize` can not. These tests fix
both halves: the default must never collide, and a too-small explicit size must
say so.
"""
import itertools
import warnings

import numpy as np
import pytest
import scipy.sparse as sp

import truecell as tc

# Long enough to actually crowd a panel. The real tutorial labels
# ("FCGR3A+ Mono") are in this range; the last two are deliberately worse.
LONG_LABELS = [
    "Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B cells", "CD8 T",
    "FCGR3A+ Mono", "NK", "Dendritic", "Platelet", "Plasmacytoid DC",
    "Erythrocyte precursor", "Haematopoietic stem cell",
]


@pytest.fixture(scope="module")
def obj():
    pytest.importorskip("matplotlib")
    rng = np.random.default_rng(0)
    n_genes, n_cells = 24, 600
    genes = [f"G{i:02d}" for i in range(n_genes)]
    o = tc.create_truecell_object(
        counts=sp.csc_matrix(rng.poisson(3.0, size=(n_genes, n_cells)).astype(float)),
        assay="RNA", feature_names=genes,
        cell_names=[f"C{i:03d}" for i in range(n_cells)], project="layout",
    )
    tc.normalize_data(o)
    o.meta_data["ct"] = [LONG_LABELS[i % 12] for i in range(n_cells)]
    return o


def worst_overlap(fig) -> float:
    """Largest overlap between any two panels' *tight* boxes, as a figure
    fraction. The tight box includes tick labels and the title, so this measures
    exactly what "labels printed over the neighbouring panel" means."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    boxes = [ax.get_tightbbox(renderer).transformed(inv)
             for ax in fig.axes if ax.get_visible()]
    worst = 0.0
    for a, b in itertools.combinations(boxes, 2):
        ox = min(a.x1, b.x1) - max(a.x0, b.x0)
        oy = min(a.y1, b.y1) - max(a.y0, b.y0)
        if ox > 0 and oy > 0:
            worst = max(worst, min(ox, oy))
    return worst


@pytest.mark.parametrize("n_features", [2, 4, 6])
@pytest.mark.parametrize("ncol", [None, 2, 3])
def test_default_figsize_never_collides(obj, n_features, ncol):
    """The property the computed default exists to have. Twelve groups with the
    longest labels in the set is the crowded end of realistic usage."""
    plt = pytest.importorskip("matplotlib.pyplot")
    genes = obj.assays["RNA"]._all_feature_names[:n_features]
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)  # must not warn either
        fig = tc.ridge_plot(obj, genes, group_by="ct", ncol=ncol)
    assert worst_overlap(fig) <= 0.005
    plt.close(fig)


def test_too_small_explicit_figsize_warns(obj):
    plt = pytest.importorskip("matplotlib.pyplot")
    genes = obj.assays["RNA"]._all_feature_names[:4]
    with pytest.warns(UserWarning, match="panels overlap"):
        fig = tc.ridge_plot(obj, genes, group_by="ct", figsize=(8, 5))
    plt.close(fig)


def test_the_suggested_figsize_actually_resolves_it(obj):
    """A suggestion that still collides is worse than none — it costs the reader
    a second round trip. Parse the size out of the warning and check it."""
    import re
    plt = pytest.importorskip("matplotlib.pyplot")
    genes = obj.assays["RNA"]._all_feature_names[:4]

    with pytest.warns(UserWarning, match="panels overlap") as rec:
        plt.close(tc.ridge_plot(obj, genes, group_by="ct", figsize=(8, 5)))
    msg = str(rec[0].message)

    m = re.search(r"Use figsize=\(([\d.]+), ([\d.]+)\)", msg)
    assert m, f"warning carries no parseable figsize: {msg}"
    suggested = (float(m.group(1)), float(m.group(2)))

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        fig = tc.ridge_plot(obj, genes, group_by="ct", figsize=suggested)
    assert worst_overlap(fig) <= 0.005
    plt.close(fig)


def test_a_roomy_explicit_figsize_stays_silent(obj):
    plt = pytest.importorskip("matplotlib.pyplot")
    genes = obj.assays["RNA"]._all_feature_names[:4]
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        fig = tc.ridge_plot(obj, genes, group_by="ct", figsize=(20, 8))
    plt.close(fig)


def test_short_labels_fit_in_a_small_figure(obj):
    """The warning must key on whether the panels actually collide, not on the
    figure being small — a small figure with short labels is fine."""
    plt = pytest.importorskip("matplotlib.pyplot")
    obj.meta_data["short"] = [["A", "B", "C", "D"][i % 4]
                              for i in range(obj.meta_data.shape[0])]
    genes = obj.assays["RNA"]._all_feature_names[:2]
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        fig = tc.ridge_plot(obj, genes, group_by="short", figsize=(8, 5))
    assert worst_overlap(fig) <= 0.005
    plt.close(fig)
