"""The plotting theme, the ggplot palette, and rasterisation.

Three things are asserted here that the module could not previously state:

  * `hue_pal` reproduces `scales::hue_pal()` exactly, so group colours match
    Seurat at every n rather than only at the one length a fixed list was
    copied from;
  * the theme scales every text element from one base size, and leaves output
    unchanged at its default;
  * `raster=None` resolves to Seurat's own rule.
"""
import io
import re

import numpy as np
import pytest

from truecell.plotting import (
    _FONT_ROLES,
    _RASTER_THRESHOLD,
    _fs,
    _palette,
    _should_raster,
    get_theme,
    hue_pal,
    reset_theme,
    set_theme,
    theme_context,
)

_HEX = re.compile(r"^#[0-9A-F]{6}$")


@pytest.fixture(autouse=True)
def _clean_theme():
    """No test may leak a theme into the next one — it is module-global state."""
    reset_theme()
    yield
    reset_theme()


# ---------------------------------------------------------------------------
# hue_pal against R
# ---------------------------------------------------------------------------

# Ground truth: `scales::hue_pal()(n)` in R, which is what ggplot — and so
# Seurat — assigns to n discrete groups. These are the values to beat; the
# implementation is a from-scratch polarLUV -> sRGB conversion, so an error
# anywhere in it shows up as a wrong hex here.
_R_HUE_PAL = {
    1: ["#F8766D"],
    2: ["#F8766D", "#00BFC4"],
    3: ["#F8766D", "#00BA38", "#619CFF"],
    4: ["#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"],
    5: ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"],
    6: ["#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"],
    8: ["#F8766D", "#CD9600", "#7CAE00", "#00BE67",
        "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC"],
    9: ["#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F",
        "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3"],
}


@pytest.mark.parametrize("n,expected", sorted(_R_HUE_PAL.items()))
def test_hue_pal_matches_r_exactly(n, expected):
    assert hue_pal(n) == expected


def test_hue_pal_depends_on_n_not_just_position():
    """The property a fixed list cannot have.

    ggplot spreads the hue circle across however many groups there are, so the
    colours for 9 groups are not the colours for 8 plus one more. The old
    36-entry list was built by concatenating several such runs, which is how it
    ended up with duplicates.
    """
    assert hue_pal(9)[1] != hue_pal(8)[1]
    assert hue_pal(5)[1] != hue_pal(8)[1]
    # Only the first is stable across n — it is always hue 15.
    assert {hue_pal(n)[0] for n in range(1, 20)} == {"#F8766D"}


@pytest.mark.parametrize("n", [1, 2, 5, 29, 30, 31, 36, 37, 64, 200])
def test_palette_never_repeats_a_colour(n):
    """The defect this replaced: `_palette` sliced a 36-name list holding 30
    distinct colours, so 30-36 groups silently drew two clusters identically."""
    colors = _palette(n)
    assert len(colors) == n
    assert len(set(colors)) == n, "two groups would render in the same colour"
    assert all(_HEX.match(c) for c in colors)


def test_hue_pal_rejects_nothing_and_returns_empty_below_one():
    assert hue_pal(0) == []
    assert hue_pal(-3) == []


# ---------------------------------------------------------------------------
# Theme
# ---------------------------------------------------------------------------

def test_default_base_size_reproduces_the_previous_absolute_sizes():
    """The refactor replaced 40 hard-coded sizes with multipliers. At the
    default base of 10 they must come out as the exact numbers they replaced,
    or every committed figure shifts for no reason."""
    assert _fs("tiny") == 7
    assert _fs("small") == 8
    assert _fs("label") == 9
    assert _fs("title") == 11
    assert _fs("large") == 12
    assert _fs("suptitle") == 13


def test_base_size_scales_every_role_together():
    set_theme(base_size=20)
    assert _fs("small") == 16
    assert _fs("suptitle") == 26
    # Ordering is preserved, which is what makes it a hierarchy rather than a
    # set of unrelated numbers.
    sizes = [_fs(r) for r in ("tiny", "small", "label", "title", "large", "suptitle")]
    assert sizes == sorted(sizes)


def test_roles_are_ordered_by_construction():
    mults = [_FONT_ROLES[r] for r in
             ("tiny", "small", "label", "title", "large", "suptitle")]
    assert mults == sorted(mults)


def test_set_theme_only_changes_what_it_is_given():
    set_theme(base_size=14)
    set_theme(font="DejaVu Sans")
    assert get_theme()["base_size"] == 14
    assert get_theme()["font"] == "DejaVu Sans"


def test_reset_theme_restores_defaults():
    set_theme(base_size=30, palette=["#000000"], style="seurat")
    reset_theme()
    t = get_theme()
    assert t["base_size"] == 10.0
    assert t["palette"] is None
    assert t["style"] is None
    assert _fs("label") == 9


def test_theme_context_restores_on_exit_and_on_error():
    with theme_context(base_size=18):
        assert _fs("label") == pytest.approx(16.2)
    assert _fs("label") == 9

    with pytest.raises(RuntimeError), theme_context(base_size=18):
        raise RuntimeError("boom")
    assert _fs("label") == 9, "an exception inside the block must still restore"


def test_set_theme_rejects_bad_input():
    with pytest.raises(ValueError, match="base_size must be positive"):
        set_theme(base_size=0)
    with pytest.raises(ValueError, match="Unknown style"):
        set_theme(style="nonexistent")


def test_theme_palette_is_used_and_falls_back_rather_than_repeating():
    set_theme(palette=["#111111", "#222222", "#333333"])
    assert _palette(2) == ["#111111", "#222222"]
    assert _palette(3) == ["#111111", "#222222", "#333333"]
    # More groups than the palette has entries: cycling would repeat a colour,
    # so hue_pal takes over instead.
    got = _palette(6)
    assert len(set(got)) == 6
    assert got == hue_pal(6)


def test_base_size_reaches_text_this_module_does_not_size_itself():
    """The roles only cover text with an explicit ``fontsize=``. Axis labels and
    tick labels are left to matplotlib, so `base_size` has to land in
    ``rcParams["font.size"]`` too — otherwise raising it grew the titles and left
    the axis furniture at 10pt, which is what it did before this was added.
    """
    plt = pytest.importorskip("matplotlib.pyplot")
    set_theme(base_size=20)
    fig, ax = plt.subplots()
    ax.set_xlabel("UMAP_1")
    ax.plot([1, 2], [1, 2])
    assert ax.xaxis.label.get_fontsize() == 20
    assert ax.get_xticklabels()[0].get_fontsize() == 20
    plt.close(fig)


def test_default_base_size_leaves_matplotlibs_font_size_where_it_was():
    """10.0 is matplotlib's own default, so the theme writing it is a no-op —
    which is what keeps every committed figure unchanged."""
    plt = pytest.importorskip("matplotlib.pyplot")
    set_theme(base_size=10)
    assert plt.rcParams["font.size"] == 10.0


def test_style_preset_writes_rcparams():
    plt = pytest.importorskip("matplotlib.pyplot")
    set_theme(style="seurat")
    assert plt.rcParams["axes.grid"] is False
    assert plt.rcParams["axes.spines.top"] is False
    reset_theme()


# ---------------------------------------------------------------------------
# Rasterisation
# ---------------------------------------------------------------------------

def test_should_raster_follows_seurats_rule():
    assert _RASTER_THRESHOLD == 100_000
    assert _should_raster(None, 1_000) is False
    assert _should_raster(None, _RASTER_THRESHOLD) is False
    assert _should_raster(None, _RASTER_THRESHOLD + 1) is True
    # An explicit choice always wins over the automatic one.
    assert _should_raster(True, 1) is True
    assert _should_raster(False, 10_000_000) is False


def test_raster_flag_reaches_the_pdf(tmp_path):
    """The property that actually matters, checked on real output rather than
    on the flag: an unrasterised scatter emits one path per point into the PDF,
    a rasterised one emits a single image. Saving both and comparing size is
    the cheapest way to tell them apart without parsing PDF internals.
    """
    plt = pytest.importorskip("matplotlib.pyplot")
    rng = np.random.default_rng(0)
    xy = rng.normal(size=(5_000, 2))

    sizes = {}
    for flag in (False, True):
        fig, ax = plt.subplots()
        ax.scatter(xy[:, 0], xy[:, 1], s=3, linewidths=0, rasterized=flag)
        buf = io.BytesIO()
        fig.savefig(buf, format="pdf")
        sizes[flag] = len(buf.getvalue())
        plt.close(fig)

    assert sizes[True] < sizes[False], (
        f"rasterised PDF ({sizes[True]} B) should be smaller than vector "
        f"({sizes[False]} B) at 5,000 points"
    )
