"""Tests for Visium tissue-image support (load_visium image=..., VisiumV2)."""
import json
import sys

import numpy as np
import pandas as pd
import pytest
import scipy.io as sio
import scipy.sparse as sp

from shanuz import VisiumV2, load_visium
from shanuz.spatial.visium import ScaleFactors, read_scale_factors, read_tissue_image

SPOT_DIAMETER = 20.0
HIRES_SCALEF = 0.1
LOWRES_SCALEF = 0.03

GENES = ["Gad1", "Slc17a7", "Sox9"]
N_SPOTS = 6


def _write_visium(tmp_path, with_image=True, with_scalefactors=True,
                  resolutions=("hires", "lowres"), all_in_tissue=True):
    """Write a minimal Visium bundle: MTX triplet + spatial/."""
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
    in_tissue = [1] * N_SPOTS if all_in_tissue else [1, 1, 1, 1, 0, 0]
    pd.DataFrame({
        "barcode": barcodes,
        "in_tissue": in_tissue,
        "array_row": np.arange(N_SPOTS),
        "array_col": np.arange(N_SPOTS),
        # y = row, x = col — in FULL-RESOLUTION pixels.
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
        for res, size in (("hires", 40), ("lowres", 12)):
            if res in resolutions:
                mpimg.imsave(spatial / f"tissue_{res}_image.png",
                             rng.random((size, size, 3)))
    return barcodes


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def test_read_scale_factors(tmp_path):
    _write_visium(tmp_path)
    sf = read_scale_factors(tmp_path / "spatial")            # dir or file both work
    assert isinstance(sf, ScaleFactors)
    assert sf.spot == SPOT_DIAMETER
    assert sf.hires == HIRES_SCALEF
    assert sf.scale_factor("lowres") == LOWRES_SCALEF
    with pytest.raises(ValueError, match="hires"):
        sf.scale_factor("fullres")


def test_read_tissue_image_falls_back_to_other_resolution(tmp_path):
    _write_visium(tmp_path, resolutions=("lowres",))
    img, res = read_tissue_image(tmp_path / "spatial", resolution="hires")
    assert res == "lowres"                                   # hires absent -> fall back
    assert img.shape[:2] == (12, 12)


def test_read_tissue_image_returns_none_without_png(tmp_path):
    _write_visium(tmp_path, with_image=False)
    assert read_tissue_image(tmp_path / "spatial") is None


# ---------------------------------------------------------------------------
# The image must be a function of the file, not of the environment.
#
# matplotlib returns float in [0, 1] for an 8-bit PNG; Pillow returns uint8 in
# [0, 255]. Neither is a declared dependency, so before normalising, the same
# bundle gave arrays 255x apart depending on which happened to be installed.
# Asserting the dtype on one backend proves nothing — the two have to be read
# and compared against each other.
# ---------------------------------------------------------------------------

def _no_matplotlib(monkeypatch):
    """Force _imread's Pillow branch. sys.modules[name] = None makes import raise."""
    monkeypatch.setitem(sys.modules, "matplotlib.image", None)


def test_both_image_backends_return_the_same_array(tmp_path, monkeypatch):
    _write_visium(tmp_path)
    via_mpl, _ = read_tissue_image(tmp_path / "spatial", resolution="lowres")
    with monkeypatch.context() as m:
        _no_matplotlib(m)
        via_pil, _ = read_tissue_image(tmp_path / "spatial", resolution="lowres")

    assert via_pil.dtype == via_mpl.dtype
    np.testing.assert_array_equal(via_pil, via_mpl)


@pytest.mark.parametrize("pillow_only", [False, True])
def test_image_is_unit_float_on_either_backend(tmp_path, monkeypatch, pillow_only):
    _write_visium(tmp_path)
    with monkeypatch.context() as m:
        if pillow_only:
            _no_matplotlib(m)
        img, _ = read_tissue_image(tmp_path / "spatial", resolution="lowres")

    assert np.issubdtype(img.dtype, np.floating), f"got {img.dtype}"
    assert 0.0 <= img.min() and img.max() <= 1.0, f"range [{img.min()}, {img.max()}]"


def test_the_pillow_branch_is_actually_reached(tmp_path, monkeypatch):
    """Guard the guard: if the monkeypatch stopped working the two tests above
    would silently both exercise matplotlib and could never fail."""
    _write_visium(tmp_path)
    with monkeypatch.context() as m:
        _no_matplotlib(m)
        with pytest.raises(ImportError):
            import matplotlib.image  # noqa: F401


@pytest.mark.parametrize("pillow_only", [False, True])
def test_load_visium_image_is_unit_float(tmp_path, monkeypatch, pillow_only):
    """Through the real entry point, not just the helper — and on both backends,
    since with matplotlib installed the Pillow branch is never reached and the
    assertion below cannot fail."""
    _write_visium(tmp_path)                      # writes the PNG; needs matplotlib
    with monkeypatch.context() as m:
        if pillow_only:
            _no_matplotlib(m)
        img = load_visium(tmp_path).images["slice1"].get_image()
    assert np.issubdtype(img.dtype, np.floating)
    assert img.max() <= 1.0


# ---------------------------------------------------------------------------
# load_visium
# ---------------------------------------------------------------------------

def test_load_visium_defaults_match_seurat(tmp_path):
    """The defaults themselves, not just the helpers behind them.

    `Read10X_Image` reads tissue_lowres_image.png, filters to in-tissue spots,
    and `Load10X_Spatial` keys the image 'slice1'. A test that passes every
    argument explicitly would keep passing if any of these drifted.
    """
    _write_visium(tmp_path, all_in_tissue=False)             # 2 of 6 off tissue
    obj = load_visium(tmp_path)

    assert list(obj.images) == ["slice1"]                    # not 'spatial'
    (fov,) = obj.images.values()
    assert fov.image_resolution == "lowres"                  # not 'hires'
    assert len(obj.cell_names()) == 4                        # off-tissue spots dropped


def test_load_visium_attaches_image_and_scalefactors(tmp_path):
    barcodes = _write_visium(tmp_path)

    obj = load_visium(tmp_path)

    assert obj.cell_names() == barcodes
    (fov,) = obj.images.values()
    assert isinstance(fov, VisiumV2)
    assert fov.image_resolution == "lowres"
    assert fov.get_image().shape[:2] == (12, 12)
    assert fov.scale_factors.spot == SPOT_DIAMETER
    # Spot radius is half the diameter, in fullres pixels.
    assert fov.radius() == SPOT_DIAMETER / 2
    assert fov.spot_radius() == SPOT_DIAMETER / 2 * LOWRES_SCALEF


def test_coordinates_stay_fullres_and_scale_on_demand(tmp_path):
    _write_visium(tmp_path)
    obj = load_visium(tmp_path)
    (fov,) = obj.images.values()

    coords = fov.get_tissue_coordinates()
    np.testing.assert_allclose(coords["x"].to_numpy(),
                               np.arange(N_SPOTS) * 200.0)     # pxl_col_in_fullres
    np.testing.assert_allclose(coords["y"].to_numpy(),
                               np.arange(N_SPOTS) * 100.0)     # pxl_row_in_fullres

    scaled = fov.scale_coordinates()
    np.testing.assert_allclose(scaled["x"].to_numpy(),
                               coords["x"].to_numpy() * LOWRES_SCALEF)
    np.testing.assert_allclose(scaled["y"].to_numpy(),
                               coords["y"].to_numpy() * LOWRES_SCALEF)
    # An explicit resolution overrides the stored one.
    hi = fov.scale_coordinates(resolution="hires")
    np.testing.assert_allclose(hi["x"].to_numpy(),
                               coords["x"].to_numpy() * HIRES_SCALEF)


def test_load_visium_hires_request(tmp_path):
    _write_visium(tmp_path)
    obj = load_visium(tmp_path, image_resolution="hires")
    (fov,) = obj.images.values()
    assert fov.image_resolution == "hires"
    assert fov.get_image().shape[:2] == (40, 40)
    assert fov.spot_radius() == SPOT_DIAMETER / 2 * HIRES_SCALEF


def test_load_visium_without_image_is_a_plain_fov(tmp_path):
    _write_visium(tmp_path)
    obj = load_visium(tmp_path, image=False)
    (fov,) = obj.images.values()
    assert not isinstance(fov, VisiumV2)
    assert fov.get_image() is None                            # base class behaviour


def test_bundle_with_no_spatial_extras_still_loads(tmp_path):
    """Backwards compatible: a positions-only bundle loads as it always did."""
    barcodes = _write_visium(tmp_path, with_image=False, with_scalefactors=False)
    obj = load_visium(tmp_path)
    (fov,) = obj.images.values()
    assert not isinstance(fov, VisiumV2)
    assert obj.cell_names() == barcodes
    assert len(obj.get_tissue_coordinates()) == N_SPOTS


def test_scalefactors_without_png_still_gives_visium_fov(tmp_path):
    """Radius/scale factors are useful even when the PNG is missing."""
    _write_visium(tmp_path, with_image=False)
    obj = load_visium(tmp_path)
    (fov,) = obj.images.values()
    assert isinstance(fov, VisiumV2)
    assert fov.get_image() is None
    assert fov.radius() == SPOT_DIAMETER / 2


def test_filter_by_tissue_drops_offtissue_spots(tmp_path):
    _write_visium(tmp_path, all_in_tissue=False)

    obj = load_visium(tmp_path)                          # filtering is now the default
    assert len(obj.cell_names()) == 4                    # 2 spots had in_tissue == 0
    assert len(obj.assays["Spatial"].cells()) == 4       # dropped from the matrix too
    assert len(obj.get_tissue_coordinates()) == 4

    # Opting out keeps the off-tissue spots.
    kept = load_visium(tmp_path, filter_by_tissue=False)
    assert len(kept.cell_names()) == N_SPOTS


# ---------------------------------------------------------------------------
# VisiumV2 behaviour
# ---------------------------------------------------------------------------

def test_subset_preserves_image_and_scalefactors(tmp_path):
    barcodes = _write_visium(tmp_path)
    obj = load_visium(tmp_path)

    sub = obj.subset(cells=barcodes[:3])
    (fov,) = sub.images.values()
    assert isinstance(fov, VisiumV2)
    assert fov.cells() == barcodes[:3]
    assert fov.get_image().shape[:2] == (12, 12)
    assert fov.radius() == SPOT_DIAMETER / 2


def test_rename_cells_preserves_image(tmp_path):
    barcodes = _write_visium(tmp_path)
    obj = load_visium(tmp_path)
    (fov,) = obj.images.values()

    renamed = fov.rename_cells([f"new-{i}" for i in range(len(barcodes))])
    assert isinstance(renamed, VisiumV2)
    assert renamed.cells() == [f"new-{i}" for i in range(len(barcodes))]
    assert renamed.get_image().shape[:2] == (12, 12)
    assert renamed.scale_factors.spot == SPOT_DIAMETER


def test_get_image_generic_dispatches(tmp_path):
    """shanuz.generics.get_image(image) works, as R's GetImage(obj[['slice1']]) does."""
    from shanuz.generics import get_image, radius

    _write_visium(tmp_path)
    (visium,) = load_visium(tmp_path).images.values()
    (plain,) = load_visium(tmp_path, image=False).images.values()

    assert get_image(visium).shape[:2] == (12, 12)
    assert radius(visium) == SPOT_DIAMETER / 2
    assert get_image(plain) is None          # non-Visium images have no photo


def test_visium_fov_without_scalefactors_has_no_radius():
    fov = VisiumV2()
    assert fov.radius() is None
    assert fov.spot_radius() is None
    assert fov.scale_factor() == 1.0        # unscaled: coords pass through unchanged
    assert fov.get_image() is None
