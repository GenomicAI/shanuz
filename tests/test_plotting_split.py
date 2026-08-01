"""`split_by` faceting, matching what Seurat actually does for each plot.

The three are not the same mechanism, which is the thing worth pinning down.
Probed against Seurat 5.5.1:

  DimPlot     facet_wrap  — one panel per level, only that level's cells,
                            `scales = "fixed"` so the ranges are shared
  FeaturePlot patchwork   — a features x levels grid
  VlnPlot     dodge       — `split.plot = FALSE` is the default, and its own
                            message says "plotted side-by-side"; the levels sit
                            within each group's x position, not in new panels
"""
import numpy as np
import pytest
import scipy.sparse as sp

import truecell as tc


@pytest.fixture(scope="module")
def obj():
    pytest.importorskip("matplotlib")
    rng = np.random.default_rng(0)
    n_genes, n_cells = 20, 300
    counts = rng.poisson(2.0, size=(n_genes, n_cells)).astype(float)
    counts[0, :100] = rng.poisson(25, size=100)
    o = tc.create_truecell_object(
        counts=sp.csc_matrix(counts), assay="RNA",
        feature_names=[f"G{i:02d}" for i in range(n_genes)],
        cell_names=[f"C{i:03d}" for i in range(n_cells)], project="split",
    )
    tc.normalize_data(o)
    tc.find_variable_features(o, nfeatures=15)
    tc.scale_data(o, features=o.assays["RNA"]._all_feature_names)
    tc.run_pca(o, n_pcs=6)
    tc.run_umap(o, dims=range(5), seed=42)
    # The split columns must be independent of `cl`, or the cross is degenerate.
    # `cl = i % 4` with `stim = i % 2` confounds them perfectly — every c0 cell
    # is CTRL, every c1 is STIM — so half the group x level cells are empty and
    # a dodge test silently sees four violins where it expects eight. Stepping
    # the split every four cells decorrelates them.
    o.meta_data["cl"] = [["c0", "c1", "c2", "c3"][i % 4] for i in range(n_cells)]
    o.meta_data["stim"] = [["CTRL", "STIM"][(i // 4) % 2] for i in range(n_cells)]
    o.meta_data["batch"] = [["b1", "b2", "b3"][(i // 4) % 3] for i in range(n_cells)]

    cross = o.meta_data.groupby(["cl", "stim"], observed=True).size()
    assert len(cross) == 8 and cross.min() > 2, "fixture is degenerate"
    return o


def _scatters(ax):
    from matplotlib.collections import PathCollection
    return [c for c in ax.collections if isinstance(c, PathCollection)]


def _visible(fig):
    return [a for a in fig.axes if a.get_visible() and a.get_position().width > 0.02]


# ---------------------------------------------------------------------------
# dim_plot — facet_wrap
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("col,n_levels", [("stim", 2), ("batch", 3), ("cl", 4)])
def test_dim_plot_makes_one_panel_per_level(obj, col, n_levels):
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = tc.dim_plot(obj, group_by="cl", split_by=col, label=False)
    titles = [a.get_title() for a in _visible(fig) if a.get_title()]
    levels = sorted(set(obj.meta_data[col].astype(str)))
    assert sorted(titles) == levels
    assert len(titles) == n_levels
    plt.close(fig)


def test_dim_plot_panels_hold_only_their_own_cells(obj):
    """facet_wrap subsets; it does not grey the others out."""
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = tc.dim_plot(obj, group_by="cl", split_by="stim", label=False)
    panels = [a for a in _visible(fig) if a.get_title()]
    counts = [sum(len(c.get_offsets()) for c in _scatters(a)) for a in panels]
    expected = obj.meta_data["stim"].value_counts().sort_index().tolist()
    assert counts == expected
    assert sum(counts) == obj.meta_data.shape[0]
    plt.close(fig)


def test_dim_plot_shares_axis_limits_across_panels(obj):
    """`scales = "fixed"`. Per-panel limits would silently rescale each one, and
    the whole point of splitting is that a position means the same thing."""
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = tc.dim_plot(obj, group_by="cl", split_by="batch", label=False)
    panels = [a for a in _visible(fig) if a.get_title()]
    xlims = {tuple(np.round(a.get_xlim(), 9)) for a in panels}
    ylims = {tuple(np.round(a.get_ylim(), 9)) for a in panels}
    assert len(xlims) == 1, f"panels disagree on x range: {xlims}"
    assert len(ylims) == 1, f"panels disagree on y range: {ylims}"
    plt.close(fig)


def test_dim_plot_colours_are_consistent_across_panels(obj):
    """A group must be the same colour in every panel, or the panels cannot be
    read against each other. Colours come from all groups, not the panel's."""
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = tc.dim_plot(obj, group_by="cl", split_by="stim", label=False)
    panels = [a for a in _visible(fig) if a.get_title()]
    per_panel = [[tuple(np.round(c.get_facecolor()[0], 6)) for c in _scatters(a)]
                 for a in panels]
    assert per_panel[0] == per_panel[1]
    assert len(set(per_panel[0])) == 4, "four clusters should be four colours"
    plt.close(fig)


def test_dim_plot_split_legend_lists_every_group_once(obj):
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = tc.dim_plot(obj, group_by="cl", split_by="stim", label=False)
    assert len(fig.legends) == 1
    texts = [t.get_text() for t in fig.legends[0].get_texts()]
    assert texts == ["c0", "c1", "c2", "c3"]
    plt.close(fig)


# ---------------------------------------------------------------------------
# feature_plot — a features x levels grid
# ---------------------------------------------------------------------------

def test_feature_plot_lays_out_features_by_levels(obj):
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = tc.feature_plot(obj, ["G00", "G01"], split_by="batch")
    titles = [a.get_title() for a in fig.axes if a.get_title()]
    assert titles == ["G00 — b1", "G00 — b2", "G00 — b3",
                      "G01 — b1", "G01 — b2", "G01 — b3"]
    plt.close(fig)


def test_feature_plot_shares_the_colour_scale_along_a_row(obj):
    """Computed over all cells, not the panel's subset. Per-panel scales would
    make two very different levels look alike — each would fill its own range."""
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = tc.feature_plot(obj, ["G00"], split_by="stim")
    panels = [a for a in fig.axes if a.get_title()]
    clims = set()
    for a in panels:
        for c in _scatters(a):
            if c.get_array() is not None:
                clims.add(tuple(np.round(c.get_clim(), 9)))
    assert len(clims) == 1, f"colour scale differs across the row: {clims}"
    plt.close(fig)


def test_feature_plot_split_panels_hold_only_their_own_cells(obj):
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = tc.feature_plot(obj, ["G00"], split_by="stim")
    panels = [a for a in fig.axes if a.get_title()]
    counts = [sum(len(c.get_offsets()) for c in _scatters(a)) for a in panels]
    assert counts == obj.meta_data["stim"].value_counts().sort_index().tolist()
    plt.close(fig)


# ---------------------------------------------------------------------------
# vln_plot — a dodge, not a facet
# ---------------------------------------------------------------------------

def _violin_centres(ax):
    from matplotlib.collections import PolyCollection
    out = []
    for c in ax.collections:
        if isinstance(c, PolyCollection) and c.get_paths():
            xs = c.get_paths()[0].to_polygons()[0][:, 0]
            # Midrange, not mean: the violin is symmetric about its centre by
            # construction, so (min + max) / 2 is exact, while the mean of the
            # outline's vertices is pulled around by how they are distributed.
            out.append(round((float(xs.min()) + float(xs.max())) / 2, 6))
    return sorted(out)


def test_vln_plot_dodges_rather_than_facets(obj):
    """One panel, and the levels share each group's x position."""
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = tc.vln_plot(obj, ["G00"], group_by="cl", split_by="stim", pt_size=0)
    assert len([a for a in fig.axes if a.has_data()]) == 1, "split_by faceted"
    centres = _violin_centres(fig.axes[0])
    assert len(centres) == 8, "4 groups x 2 levels"
    plt.close(fig)


def test_vln_plot_dodge_offsets_match_ggplots_geometry(obj):
    """ggplot dodges within one group's width: with two levels the violins sit
    at +/- width/4 of the group's centre. Seurat's own build put them at 0.775
    and 1.225 around x = 1, i.e. +/- 0.225 for a width of 0.9."""
    plt = pytest.importorskip("matplotlib.pyplot")
    width = 0.8
    fig = tc.vln_plot(obj, ["G00"], group_by="cl", split_by="stim",
                      pt_size=0, violin_width=width)
    centres = np.array(_violin_centres(fig.axes[0]))
    expected = np.sort(np.concatenate([np.arange(4) - width / 4,
                                       np.arange(4) + width / 4]))
    assert centres == pytest.approx(expected, abs=1e-6)
    plt.close(fig)


def test_vln_plot_split_colours_by_level_not_by_group(obj):
    """The dodge carries the split in colour; the group is carried by position."""
    plt = pytest.importorskip("matplotlib.pyplot")
    from matplotlib.collections import PolyCollection
    fig = tc.vln_plot(obj, ["G00"], group_by="cl", split_by="stim", pt_size=0)
    cols = [tuple(np.round(c.get_facecolor()[0], 6))
            for c in fig.axes[0].collections if isinstance(c, PolyCollection)]
    assert len(set(cols)) == 2, "should be one colour per split level"
    plt.close(fig)


def test_vln_plot_split_violins_actually_hold_different_data(obj):
    """Position and colour say a split happened; only the shape proves it did.

    Dropping the `split_labels == lv` mask — so both levels draw every cell in
    the group — left the dodge, the colours and the tick labels exactly as they
    should be, and every other test in this file passed. The violins were
    identical twins and nothing noticed. So build a column where the two levels
    have disjoint value ranges within each group, and require the drawn shapes
    to differ.
    """
    plt = pytest.importorskip("matplotlib.pyplot")
    from matplotlib.collections import PolyCollection

    # Its own object, not the shared fixture: this one needs the *values* rigged
    # per level, and mutating a module-scoped fixture would leak into the rest.
    rng = np.random.default_rng(11)
    n = 320
    # `(i // 4) % 2` against `cl = i % 4`, for the same reason as the fixture —
    # `i % 2` would confound them and half the cells would be empty.
    cl = np.array([["c0", "c1", "c2", "c3"][i % 4] for i in range(n)])
    lohi = np.array([["lo", "hi"][(i // 4) % 2] for i in range(n)])
    counts = rng.poisson(1.0, size=(2, n)).astype(float)
    counts[1] = np.where(lohi == "hi",
                         rng.normal(40, 2, n), rng.normal(4, 2, n)).clip(0)

    o2 = tc.create_truecell_object(
        counts=sp.csc_matrix(counts), assay="RNA", feature_names=["G00", "G01"],
        cell_names=[f"X{i:03d}" for i in range(n)], project="lohi",
    )
    o2.meta_data["cl"] = cl
    o2.meta_data["lohi"] = lohi

    # Raw counts on purpose. LogNormalize divides by each cell's library size,
    # and with two genes G01 *is* most of it — normalising pulls the two levels
    # back on top of each other and the fixture stops testing anything.
    fig = tc.vln_plot(o2, ["G01"], group_by="cl", split_by="lohi",
                      layer="counts", pt_size=0)
    by_centre = {}
    for c in fig.axes[0].collections:
        if isinstance(c, PolyCollection) and c.get_paths():
            poly = c.get_paths()[0].to_polygons()[0]
            centre = (float(poly[:, 0].min()) + float(poly[:, 0].max())) / 2
            by_centre[round(centre, 4)] = (float(poly[:, 1].min()),
                                           float(poly[:, 1].max()))
    plt.close(fig)

    centres = sorted(by_centre)
    assert len(centres) == 8, "4 groups x 2 levels"
    # Within each group the two violins must cover disjoint ranges. Which side
    # holds which level follows the level order — "hi" sorts before "lo", so it
    # takes the left slot — and the claim here is only that they differ.
    for a, b in zip(centres[0::2], centres[1::2]):
        ra, rb = by_centre[a], by_centre[b]
        disjoint = rb[0] > ra[1] or ra[0] > rb[1]
        assert disjoint, (
            f"violins at x={a} and x={b} overlap in y ({ra} vs {rb}) "
            "— both levels drew the same cells"
        )


def test_vln_plot_split_x_ticks_still_label_the_groups(obj):
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = tc.vln_plot(obj, ["G00"], group_by="cl", split_by="stim", pt_size=0)
    labels = [t.get_text() for t in fig.axes[0].get_xticklabels()]
    assert labels == ["c0", "c1", "c2", "c3"]
    plt.close(fig)


def test_vln_plot_split_three_levels_stay_inside_the_group_slot(obj):
    """Adjacent groups must not collide: every violin has to stay within half a
    unit of its group centre however many levels there are."""
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = tc.vln_plot(obj, ["G00"], group_by="cl", split_by="batch",
                      pt_size=0, violin_width=0.8)
    from matplotlib.collections import PolyCollection
    for c in fig.axes[0].collections:
        if isinstance(c, PolyCollection) and c.get_paths():
            xs = c.get_paths()[0].to_polygons()[0][:, 0]
            nearest = round(float(xs.mean()))
            assert np.abs(xs - nearest).max() < 0.5
    plt.close(fig)


# ---------------------------------------------------------------------------
# shared
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("fn,kw", [
    (tc.dim_plot, {}),
    (tc.feature_plot, {"features": ["G00"]}),
    (tc.vln_plot, {"features": ["G00"]}),
])
def test_unknown_split_column_is_reported_clearly(obj, fn, kw):
    with pytest.raises(KeyError, match="nope"):
        fn(obj, split_by="nope", **kw)


@pytest.mark.parametrize("fn,kw", [
    (tc.dim_plot, {"label": False}),
    (tc.feature_plot, {"features": ["G00"]}),
    (tc.vln_plot, {"features": ["G00"], "pt_size": 0}),
])
def test_omitting_split_by_draws_a_single_panel(obj, fn, kw):
    """The unsplit path must not have grown a grid."""
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = fn(obj, **kw)
    # feature_plot attaches a colourbar, which is itself an Axes with data.
    data_axes = [a for a in fig.axes
                 if a.has_data() and a.get_label() != "<colorbar>"]
    assert len(data_axes) == 1
    plt.close(fig)
    plt.close(fig)
