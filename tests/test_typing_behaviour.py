"""Behaviour behind the last of the mypy fixes.

Companion to `test_object_model_typing.py`. Most of the annotation work is
advisory — CI runs mypy with `|| true` — but two of the sites it pointed at were
doing something wrong at runtime, and those are what is pinned here:

  * `_palette` reached for `plt.cm.get_cmap`, which matplotlib **removed** in
    3.9. Every plot with more than 36 groups raised `AttributeError` on any
    modern matplotlib, and no test ever asked for that many.
  * `_get_layer` / `_get_data_matrix` indexed an Assay5's layers with
    `default_layer`, which is `None` when the assay has none — so a layerless
    assay came back as `KeyError: None` instead of a sentence.

The annotations themselves are mypy's business and are not asserted.
"""
import re

import numpy as np
import pytest
import scipy.sparse as sp

from truecell.assay5 import Assay5
from truecell.plotting import _get_data_matrix, _palette
from truecell.preprocessing import _get_layer

_HEX = re.compile(r"^#[0-9a-f]{6}$", re.IGNORECASE)


@pytest.mark.parametrize("n", [1, 2, 8, 30, 36, 37, 40, 64])
def test_palette_is_well_formed_and_distinct_at_every_n(n):
    # Distinctness is the property, at every n rather than only the small ones.
    # It did not hold before: `_palette` sliced a 36-entry list holding 30
    # distinct colours, so 30-36 groups drew duplicates, and past 36 the
    # `tab20` fallback repeated hues by construction. Two clusters rendered in
    # the same colour is a wrong plot, not a cosmetic one.
    colors = _palette(n)
    assert len(colors) == n
    assert all(_HEX.match(c) for c in colors)
    assert len(set(colors)) == n


def test_vln_plot_survives_more_groups_than_the_fixed_palette():
    # `_palette` is private; this is the public path that used to crash. 40
    # groups, one over the point where the colormap branch takes over.
    plt = pytest.importorskip("matplotlib.pyplot")
    import pandas as pd

    from truecell.plotting import vln_plot
    from truecell.truecell import create_truecell_object

    n_groups, per = 40, 3
    n_cells = n_groups * per
    rng = np.random.default_rng(0)
    counts = sp.csc_matrix(rng.poisson(5, size=(6, n_cells)).astype(np.float64))
    cells = [f"c{i}" for i in range(n_cells)]
    meta = pd.DataFrame({"grp": [f"g{i // per}" for i in range(n_cells)]}, index=cells)
    obj = create_truecell_object(
        counts=counts, feature_names=[f"gene{i}" for i in range(6)],
        cell_names=cells, meta_data=meta, min_cells=0, min_features=0,
    )

    fig = vln_plot(obj, features="gene0", group_by="grp")
    plt.close("all")
    assert fig is not None


@pytest.mark.parametrize("getter", [_get_layer, _get_data_matrix])
def test_a_layerless_assay5_says_so(getter):
    # `default_layer` is None for an assay with no layers, so both getters used
    # to index the dict with None. The message matches Assay5.layer_data's.
    empty = Assay5(layers={}, feature_names=["g1"], cell_names=["c1"])
    assert empty.default_layer is None
    with pytest.raises(ValueError, match="No layers available"):
        getter(empty, None)


@pytest.mark.parametrize("getter", [_get_layer, _get_data_matrix])
def test_a_populated_assay5_still_resolves_without_a_named_layer(getter):
    # The guard sits after the data/counts search, so it must not intercept the
    # ordinary path — nor the third case, a layer under neither of those names.
    mat = sp.csc_matrix(np.array([[1.0, 2.0]]))
    both = Assay5(
        layers={"counts": mat, "data": mat * 2},
        feature_names=["g1"], cell_names=["c1", "c2"],
    )
    assert (getter(both, None).toarray() == (mat * 2).toarray()).all()

    odd = Assay5(
        layers={"scale.data": mat}, feature_names=["g1"], cell_names=["c1", "c2"],
    )
    assert (getter(odd, None).toarray() == mat.toarray()).all()
