"""Tests for the Visium tissue-image plots (spatial_dim_plot / spatial_feature_plot)."""
import json

import matplotlib
import numpy as np
import pandas as pd
import pytest
import scipy.io as sio
import scipy.sparse as sp

matplotlib.use("Agg")

from matplotlib.collections import EllipseCollection  # noqa: E402

from shanuz import load_visium, spatial_dim_plot, spatial_feature_plot  # noqa: E402
from shanuz.preprocessing import normalize_data  # noqa: E402

SPOT_DIAMETER = 20.0
HIRES_SCALEF = 0.1
LOWRES_SCALEF = 0.03
# Spots reach x = 1000 fullres px, so the images must be big enough to hold
# them once scaled: 1000 * 0.1 = 100 < 120 hires, 1000 * 0.03 = 30 < 40 lowres.
HIRES_PX, LOWRES_PX = 120, 40

GENES = ["Gad1", "Slc17a7", "Sox9"]
N_SPOTS = 6


def _write_visium(tmp_path, with_image=True, with_scalefactors=True):
    """A minimal Visium bundle: MTX triplet + spatial/ (mirrors test_visium_image)."""
    mtx = tmp_path / "filtered_feature_bc_matrix"
    mtx.mkdir()
    counts = sp.csr_matrix(np.arange(len(GENES) * N_SPOTS, dtype=float).reshape(
        len(GENES), N_SPOTS))
    sio.mmwrite(mtx / "matrix.mtx", counts)
    barcodes = [f"AAAC{i}-1" for i in range(N_SPOTS)]
    pd.DataFrame({"bc": barcodes}).to_csv(
        mtx / "barcodes.tsv", sep="\t", header=False, index=False)
    pd.DataFrame({"id": [f"ENSG{i}" for i in range(len(GENES))], "sym": GENES}).to_csv(
        mtx / "genes.tsv", sep="\t", header=False, index=False)

    spatial = tmp_path / "spatial"
    spatial.mkdir()
    pd.DataFrame({
        "barcode": barcodes,
        "in_tissue": [1] * N_SPOTS,
        "array_row": np.arange(N_SPOTS),
        "array_col": np.arange(N_SPOTS),
        "pxl_row_in_fullres": np.arange(N_SPOTS, dtype=float) * 100.0,
        "pxl_col_in_fullres": np.arange(N_SPOTS, dtype=float) * 200.0,
    }).to_csv(spatial / "tissue_positions.csv", index=False)

    if with_scalefactors:
        with open(spatial / "scalefactors_json.json", "w") as fh:
            json.dump({
                "spot_diameter_fullres": SPOT_DIAMETER,
                "fiducial_diameter_fullres": 30.0,
                "tissue_hires_scalef": HIRES_SCALEF,
                "tissue_lowres_scalef": LOWRES_SCALEF,
            }, fh)

    if with_image:
        import matplotlib.image as mpimg
        rng = np.random.default_rng(0)
        for res, size in (("hires", HIRES_PX), ("lowres", LOWRES_PX)):
            mpimg.imsave(spatial / f"tissue_{res}_image.png", rng.random((size, size, 3)))
    return barcodes


@pytest.fixture
def visium(tmp_path):
    barcodes = _write_visium(tmp_path)
    obj = load_visium(tmp_path)
    normalize_data(obj, assay="Spatial")
    obj.add_meta_data(
        pd.Series(["a", "b"] * (N_SPOTS // 2), index=barcodes), col_name="region")
    return obj


def _spots(fig):
    """The spot collection of the first panel, if the spots were drawn to scale."""
    for coll in fig.axes[0].collections:
        if isinstance(coll, EllipseCollection):
            return coll
    return None


def _diameters(coll):
    """Spot diameters in data units. EllipseCollection keeps half-widths, privately."""
    return np.asarray(coll._widths) * 2.0


# ---------------------------------------------------------------------------
# The image actually gets drawn
# ---------------------------------------------------------------------------

def test_spatial_dim_plot_draws_the_tissue_image(visium):
    fig = spatial_dim_plot(visium, group_by="region")
    ax = fig.axes[0]
    assert len(ax.images) == 1
    assert ax.images[0].get_array().shape[:2] == (LOWRES_PX, LOWRES_PX)


def test_spatial_feature_plot_draws_the_tissue_image(visium):
    fig = spatial_feature_plot(visium, "Gad1")
    ax = fig.axes[0]
    assert len(ax.images) == 1
    # A colourbar axis is added alongside the panel.
    assert len(fig.axes) >= 2


def test_image_alpha_is_applied(visium):
    fig = spatial_dim_plot(visium, group_by="region", image_alpha=0.3)
    assert fig.axes[0].images[0].get_alpha() == pytest.approx(0.3)


# ---------------------------------------------------------------------------
# Spots land on the image, at the right size
# ---------------------------------------------------------------------------

def test_spots_are_scaled_into_image_pixel_space(visium):
    fig = spatial_dim_plot(visium, group_by="region", crop=False)
    offsets = _spots(fig).get_offsets()
    # Fullres x = col * 200, y = row * 100; the lowres image is 0.03× that.
    expected = np.column_stack([
        np.arange(N_SPOTS) * 200.0 * LOWRES_SCALEF,
        np.arange(N_SPOTS) * 100.0 * LOWRES_SCALEF,
    ])
    assert np.allclose(np.asarray(offsets), expected)
    # ...and they land inside the image, which the raw fullres coords would not.
    assert np.asarray(offsets).max() <= LOWRES_PX


def test_spot_diameter_matches_the_scale_factor(visium):
    fig = spatial_dim_plot(visium, group_by="region", pt_size_factor=1.0)
    coll = _spots(fig)
    assert coll.get_offset_transform() is fig.axes[0].transData    # sized in data units
    assert np.allclose(_diameters(coll), SPOT_DIAMETER * LOWRES_SCALEF)


def test_pt_size_factor_scales_the_spots(visium):
    small = _diameters(_spots(spatial_dim_plot(visium, group_by="region",
                                               pt_size_factor=1.0)))
    big = _diameters(_spots(spatial_dim_plot(visium, group_by="region",
                                             pt_size_factor=2.0)))
    assert np.allclose(big, small * 2.0)


def test_hires_resolution_rescales_spots_and_image(tmp_path):
    """The non-default resolution. lowres is what load_visium now picks on its own,
    so asking for hires is the case that proves the choice is wired through."""
    _write_visium(tmp_path)
    obj = load_visium(tmp_path, image_resolution="hires")
    fig = spatial_dim_plot(obj, crop=False)
    assert fig.axes[0].images[0].get_array().shape[:2] == (HIRES_PX, HIRES_PX)
    offsets = np.asarray(_spots(fig).get_offsets())
    assert np.allclose(offsets[:, 0], np.arange(N_SPOTS) * 200.0 * HIRES_SCALEF)


# ---------------------------------------------------------------------------
# Framing
# ---------------------------------------------------------------------------

def test_y_axis_points_down_like_the_image(visium):
    ax = spatial_dim_plot(visium, group_by="region").axes[0]
    lo, hi = ax.get_ylim()
    assert lo > hi          # inverted: row 0 at the top, as in the photo


def test_crop_false_shows_the_whole_slide(visium):
    ax = spatial_dim_plot(visium, group_by="region", crop=False).axes[0]
    assert ax.get_xlim() == pytest.approx((-0.5, LOWRES_PX - 0.5))
    assert ax.get_ylim() == pytest.approx((LOWRES_PX - 0.5, -0.5))


def test_crop_zooms_to_the_spots(visium):
    ax = spatial_dim_plot(visium, group_by="region", crop=True).axes[0]
    x0, x1 = ax.get_xlim()
    # Spots span x = 0..30 in lowres px; cropping must be tighter than the slide.
    assert x1 < LOWRES_PX - 0.5
    assert x0 > -LOWRES_PX


# ---------------------------------------------------------------------------
# Graceful degradation — no image stored
# ---------------------------------------------------------------------------

def test_falls_back_to_a_scatter_without_an_image(tmp_path):
    _write_visium(tmp_path)
    obj = load_visium(tmp_path, image=False)          # plain FOV, no photo
    normalize_data(obj, assay="Spatial")
    fig = spatial_dim_plot(obj)
    ax = fig.axes[0]
    assert len(ax.images) == 0                        # nothing to draw underneath
    assert len(ax.collections) == 1                   # but the spots are still there
    lo, hi = ax.get_ylim()
    assert lo > hi                                    # same orientation as with an image


def test_fallback_keeps_fullres_coordinates(tmp_path):
    _write_visium(tmp_path)
    obj = load_visium(tmp_path, image=False)
    ax = spatial_dim_plot(obj).axes[0]
    offsets = np.asarray(ax.collections[0].get_offsets())
    assert np.allclose(offsets[:, 0], np.arange(N_SPOTS) * 200.0)   # unscaled


def test_feature_plot_falls_back_without_an_image(tmp_path):
    _write_visium(tmp_path)
    obj = load_visium(tmp_path, image=False)
    normalize_data(obj, assay="Spatial")
    fig = spatial_feature_plot(obj, "Gad1")
    assert len(fig.axes[0].images) == 0
    assert len(fig.axes) >= 2                          # colourbar still drawn


def test_no_scalefactors_means_plain_scatter_points(tmp_path):
    _write_visium(tmp_path, with_scalefactors=False)
    obj = load_visium(tmp_path)
    fig = spatial_dim_plot(obj)
    # The photo is drawn, but with no spot_diameter_fullres we cannot size spots
    # to scale, so they degrade to fixed-size scatter points.
    assert len(fig.axes[0].images) == 1
    assert _spots(fig) is None
    assert len(fig.axes[0].collections) == 1


# ---------------------------------------------------------------------------
# Colouring and errors
# ---------------------------------------------------------------------------

def test_group_by_colours_and_legend(visium):
    fig = spatial_dim_plot(visium, group_by="region",
                           cols={"a": "#ff0000", "b": "#0000ff"})
    faces = _spots(fig).get_facecolor()
    assert np.allclose(faces[0][:3], (1.0, 0.0, 0.0))     # spot 0 is region "a"
    assert np.allclose(faces[1][:3], (0.0, 0.0, 1.0))     # spot 1 is region "b"
    assert [t.get_text() for t in fig.legends[0].get_texts()] == ["a", "b"]


def test_feature_plot_colours_by_expression(visium):
    fig = spatial_feature_plot(visium, "Gad1")
    arr = np.asarray(_spots(fig).get_array(), dtype=float)
    row = visium.assays["Spatial"].layer_data("data")[GENES.index("Gad1"), :]
    expected = np.asarray(row.todense()).ravel() if sp.issparse(row) else np.ravel(row)
    assert np.allclose(arr, expected)


def test_unknown_image_name_raises(visium):
    with pytest.raises(KeyError, match="No such image"):
        spatial_dim_plot(visium, image="nope")


def test_object_without_images_raises(visium):
    visium.images = {}
    with pytest.raises(ValueError, match="no spatial images"):
        spatial_dim_plot(visium)


def test_subset_object_still_plots(visium):
    keep = visium.cell_names()[:3]
    sub = visium[keep]
    fig = spatial_dim_plot(sub, group_by="region")
    assert len(fig.axes[0].images) == 1                  # image survives subsetting
    assert len(_spots(fig).get_offsets()) == 3
