"""Violin geometry: bandwidth, trim, width scaling, and the jittered points.

A Seurat violin is `geom_violin(scale = "width", adjust = adjust, trim = TRUE)`
smoothed with R's `nrd0` bandwidth, plus jittered points whose size comes from
`AutoPointSize`. Each of those four is asserted here against values taken from
R 4.6.1 / Seurat 5.5.1 rather than from a reimplementation of the same idea.
"""
import numpy as np
import pytest
import scipy.sparse as sp

import truecell as tc
from truecell.plotting import _auto_point_size, _bw_nrd0

# `stats::bw.nrd0(x)` in R 4.6.1, printed to 12 places.
R_BW_NRD0 = {
    "arange10": (np.arange(10.0), 1.719286404692),
    "zeros_then_tail": (np.array([0.0] * 10 + [1, 2, 3, 8.0]), 0.297148380291),
    "mostly_zero": (np.array([0, 0, 0, 0, 1.0]), 0.291718187405),
    "zero_inflated": (
        np.array([0.1, 0.4, 0.9, 1.6, 2.5, 3.6, 4.9, 6.4, 8.1, 10.0, 0, 0, 0, 0, 0]),
        1.660765786387,
    ),
    "two_points": (np.array([1.0, 2.0]), 0.292349069764),
    # Zero spread takes bw.nrd0's fallback chain: sd and IQR are both 0, so it
    # falls to abs(x[1]) and then to 1. Both branches are exercised.
    "all_same_nonzero": (np.array([5.0, 5, 5, 5]), 3.410362274648),
    "all_zero": (np.array([0.0, 0, 0, 0]), 0.682072454930),
}


@pytest.mark.parametrize("name", sorted(R_BW_NRD0))
def test_bw_nrd0_matches_r(name):
    x, expected = R_BW_NRD0[name]
    assert _bw_nrd0(x) == pytest.approx(expected, abs=1e-11)


def test_bw_nrd0_uses_the_right_divisor():
    """`bw.nrd0` divides the IQR by 1.34; the neighbouring rule `bw.nrd` uses
    1.349. Taking the wrong one is a silent 0.67% error in every bandwidth where
    the IQR term wins — which, on zero-inflated expression, is most of them."""
    # The IQR term has to win *and* be non-zero. Too much zero inflation puts
    # both quartiles at zero, which sends bw.nrd0 down its fallback chain
    # instead — so a tight bulk with a heavy tail, not 80% zeros.
    x = np.array([0.0] * 20 + [1.0] * 20 + [50.0] * 5)
    hi = float(np.std(x, ddof=1))
    q75, q25 = np.percentile(x, [75, 25])
    iqr_term = (q75 - q25) / 1.34
    assert 0 < iqr_term < hi, "fixture does not exercise the IQR branch"
    assert _bw_nrd0(x) == pytest.approx(0.9 * iqr_term * x.size ** -0.2)
    # And the wrong divisor would be visibly different, not a rounding matter.
    assert _bw_nrd0(x) != pytest.approx(
        0.9 * (q75 - q25) / 1.349 * x.size ** -0.2, rel=1e-4
    )


def test_bw_nrd0_is_narrower_than_scipy_scott_on_expression_like_data():
    """The reason this exists. Scott scales the sd; nrd0 takes min(sd, IQR/1.34),
    and zero inflation makes the IQR term much the smaller. scipy's default
    over-smooths, flattening the spike at zero that is the shape of the data."""
    rng = np.random.default_rng(0)
    x = np.concatenate([np.zeros(1400), rng.gamma(2, 1, 600)])
    scott = float(np.std(x, ddof=1)) * x.size ** -0.2
    assert scott > 2 * _bw_nrd0(x)


def test_bw_nrd0_rejects_a_single_point():
    with pytest.raises(ValueError, match="at least 2"):
        _bw_nrd0(np.array([1.0]))


def test_auto_point_size_follows_seurats_rule():
    """Seurat: min(1583 / n, 1), capped so small objects get full-size points."""
    assert _auto_point_size(100) == pytest.approx(_auto_point_size(1583))
    assert _auto_point_size(1583) > _auto_point_size(2638) > _auto_point_size(10_000)
    # The cap, and the 1/n^2 falloff in matplotlib's area units.
    ratio = _auto_point_size(20_000) / _auto_point_size(10_000)
    assert ratio == pytest.approx(0.25, rel=1e-6)


@pytest.fixture(scope="module")
def obj():
    pytest.importorskip("matplotlib")
    rng = np.random.default_rng(0)
    n_genes, n_cells = 12, 400
    counts = rng.poisson(0.4, size=(n_genes, n_cells)).astype(float)
    counts[0, :120] = rng.poisson(30, size=120)  # one clearly expressed gene
    genes = [f"G{i:02d}" for i in range(n_genes)]
    o = tc.create_truecell_object(
        counts=sp.csc_matrix(counts), assay="RNA", feature_names=genes,
        cell_names=[f"C{i:03d}" for i in range(n_cells)], project="vln",
    )
    tc.normalize_data(o)
    o.meta_data["grp"] = [["a", "b", "c", "d"][i % 4] for i in range(n_cells)]
    return o


def _violin_paths(ax):
    from matplotlib.collections import PolyCollection
    return [c for c in ax.collections if isinstance(c, PolyCollection)]


def test_violins_are_trimmed_to_the_observed_range(obj):
    """trim = TRUE. Expression cannot be negative, so a violin must not tail off
    below the data — an untrimmed gaussian KDE does exactly that."""
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = tc.vln_plot(obj, ["G00"], group_by="grp", pt_size=0)
    ax = fig.axes[0]
    data_min = tc.plotting._get_expression(obj, "G00").min()
    for coll in _violin_paths(ax):
        ys = np.concatenate([p[:, 1] for p in coll.get_paths()[0].to_polygons()])
        assert ys.min() >= data_min - 1e-9, "violin extends below the data"
    assert data_min >= 0
    plt.close(fig)


def test_every_violin_reaches_the_same_maximum_width(obj):
    """scale = "width" — groups are compared on shape, not on cell count."""
    plt = pytest.importorskip("matplotlib.pyplot")
    width = 0.8
    fig = tc.vln_plot(obj, ["G00"], group_by="grp", pt_size=0, violin_width=width)
    widths = []
    for coll in _violin_paths(fig.axes[0]):
        poly = coll.get_paths()[0].to_polygons()[0]
        centre = np.round(poly[:, 0].mean())
        widths.append(2 * np.abs(poly[:, 0] - centre).max())
    assert len(widths) == 4
    assert all(w == pytest.approx(width, rel=1e-6) for w in widths)
    plt.close(fig)


def test_the_rendered_violin_uses_nrd0_not_scipys_default(obj):
    """The one that matters, and the one component-level bandwidth tests miss.

    Asserting `_bw_nrd0` matches R says nothing about whether `vln_plot` calls
    it: swapping the KDE to `bw_method="scott"` left every other test in this
    file green. So compare the *drawn* outline against both candidate densities
    and require it to be the nrd0 one.
    """
    from scipy.stats import gaussian_kde
    plt = pytest.importorskip("matplotlib.pyplot")

    # Zero-inflated, where the two rules diverge most.
    rng = np.random.default_rng(3)
    n = 600
    vals = np.concatenate([np.zeros(420), rng.gamma(2.0, 1.0, n - 420)])
    o = tc.create_truecell_object(
        counts=sp.csc_matrix(np.vstack([vals, vals])),
        assay="RNA", feature_names=["F0", "F1"],
        cell_names=[f"C{i:03d}" for i in range(n)], project="bw",
    )
    o.assays["RNA"].layers["data"] = o.assays["RNA"].layers["counts"].copy()
    o.meta_data["one"] = ["g"] * n

    width = 0.8
    fig = tc.vln_plot(o, ["F0"], group_by="one", pt_size=0, violin_width=width)
    poly = _violin_paths(fig.axes[0])[0].get_paths()[0].to_polygons()[0]
    plt.close(fig)

    # Right edge of the drawn outline, as half-width against y.
    right = poly[poly[:, 0] >= poly[:, 0].mean()]
    order = np.argsort(right[:, 1])
    ys, half = right[order, 1], right[order, 0] - poly[:, 0].mean()

    grid = np.linspace(vals.min(), vals.max(), 512)
    sd = float(np.std(vals, ddof=1))
    nrd0 = gaussian_kde(vals, bw_method=_bw_nrd0(vals) / sd)(grid)
    scott = gaussian_kde(vals, bw_method="scott")(grid)
    pred_nrd0 = np.interp(ys, grid, nrd0 / nrd0.max()) * width / 2
    pred_scott = np.interp(ys, grid, scott / scott.max()) * width / 2

    err_nrd0 = np.abs(half - pred_nrd0).max()
    err_scott = np.abs(half - pred_scott).max()

    # Measured separation is ~590x (nrd0 4e-4 against scott 2.3e-1). The residual
    # on the nrd0 side is polygon-extraction noise — the closing edges of the
    # filled path and interpolation onto its vertices — not a bandwidth error,
    # so the bound is loose in absolute terms and tight relative to the gap.
    assert err_nrd0 < 5e-3, f"outline does not match nrd0 (err {err_nrd0:.2e})"
    assert err_scott > 50 * err_nrd0, (
        f"nrd0 (err {err_nrd0:.2e}) and scott (err {err_scott:.2e}) are not "
        "distinguishable here — the fixture is too weak to prove anything"
    )


def test_points_are_drawn_by_default(obj):
    """Seurat's VlnPlot passes pt.size = NULL, which ExIPlot resolves through
    AutoPointSize — so points are on unless asked otherwise."""
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = tc.vln_plot(obj, ["G00"], group_by="grp")
    from matplotlib.collections import PathCollection
    scatters = [c for c in fig.axes[0].collections if isinstance(c, PathCollection)]
    assert scatters, "no jittered points drawn"
    assert sum(len(c.get_offsets()) for c in scatters) == obj.meta_data.shape[0]
    plt.close(fig)


def test_pt_size_zero_still_suppresses_points(obj):
    plt = pytest.importorskip("matplotlib.pyplot")
    fig = tc.vln_plot(obj, ["G00"], group_by="grp", pt_size=0)
    from matplotlib.collections import PathCollection
    assert not [c for c in fig.axes[0].collections
                if isinstance(c, PathCollection)]
    plt.close(fig)


def test_jitter_is_reproducible_by_default(obj):
    """A figure that redraws differently every call cannot be diffed, and the
    tutorial figures are committed."""
    plt = pytest.importorskip("matplotlib.pyplot")
    from matplotlib.collections import PathCollection

    def offsets(**kw):
        fig = tc.vln_plot(obj, ["G00"], group_by="grp", **kw)
        got = np.concatenate([c.get_offsets() for c in fig.axes[0].collections
                              if isinstance(c, PathCollection)])
        plt.close(fig)
        return got

    assert np.array_equal(offsets(), offsets())
    assert not np.array_equal(offsets(jitter_seed=1), offsets(jitter_seed=2))


def test_a_constant_group_still_appears(obj):
    """A group whose values are all identical has no density. It must not vanish
    from the panel — that would read as "no cells" rather than "no spread"."""
    plt = pytest.importorskip("matplotlib.pyplot")
    obj.meta_data["const"] = ["only"] * obj.meta_data.shape[0]
    fig = tc.vln_plot(obj, ["G11"], group_by="const", pt_size=0)
    ax = fig.axes[0]
    assert _violin_paths(ax) or ax.lines, "constant group drew nothing"
    plt.close(fig)
