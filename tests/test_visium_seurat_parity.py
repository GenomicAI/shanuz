"""What the Visium side-by-side established, pinned so it cannot quietly rot.

Constants marked "R:" were read off a live Seurat 5.5.1 session on the 10x
V1_Mouse_Brain_Sagittal_Anterior bundle, not derived from shanuz's own output.
The tests need neither R nor the dataset — CI has no R — so the R side is
transcribed rather than invoked, as `test_lazy_bpcells_parity.py` does.

The finding worth protecting is unusual for this suite: **here shanuz is right
and Seurat is not.** Seurat builds the Visium FOV with
``radius = scale.factors[["spot"]]``, and that field holds
``spot_diameter_fullres`` — a diameter in a slot named radius. So the tests
below pin a *deliberate* 2x divergence from Seurat, and pin the arithmetic that
licenses it, so that a later "parity fix" cannot quietly introduce the bug.
"""
import json
from pathlib import Path

import numpy as np
import pytest

TUTORIALS = Path(__file__).resolve().parent.parent / "tutorials"
FIGURES = TUTORIALS / "figures_visium"

# R: Seurat 5.5.1, Read10X_Image on its own defaults.
# `Rscript tutorials/visium_verify.R` reproduces every one of these.
SEURAT_N_SPOTS = 2695                     # filter.matrix = TRUE
SEURAT_N_SPOTS_UNFILTERED = 4992          # every spot on the capture area
SEURAT_N_FEATURES = 32285
SEURAT_RADIUS = 89.47199235723474         # == spot_diameter_fullres, the diameter
SEURAT_RADIUS_ON_VISIUMV2_IS_NULL = True  # no Radius.VisiumV2 method exists
SEURAT_HVG_OVERLAP = 1995                 # out of 2000
# R: PCA stdev, max relative difference against Seurat's, measured two ways.
SEURAT_PCA_REL_OWN_FEATURES = 0.00194     # shanuz picks 1995 of the same 2000
SEURAT_PCA_REL_SHARED_FEATURES = 2.49e-05  # same features -> the PCA itself

# 10x Visium slide geometry. Fixed by the hardware, not by either tool — this is
# the evidence that decides which reading of spot_diameter_fullres is correct.
SPOT_PITCH_UM = 100.0
REFERENCE_SPOT_UM = 65.0
MEASURED_NN_SPACING_PX = 137.0            # median, this slide


def _tutorial():
    import sys

    if str(TUTORIALS.parent) not in sys.path:
        sys.path.insert(0, str(TUTORIALS.parent))
    from tutorials import visium_tutorial

    return visium_tutorial


# ---------------------------------------------------------------------------
# The R-side finding, and the arithmetic that establishes it
# ---------------------------------------------------------------------------

def test_reading_the_spot_diameter_as_a_radius_implies_overlapping_spots():
    """Pure arithmetic on published slide geometry — no tool involved.

    Visium capture spots are physically distinct wells on a 100 um grid. Any
    reading of `spot_diameter_fullres` that makes them wider than the pitch is
    wrong on its face, and reading it as a radius does exactly that.
    """
    px_per_um = MEASURED_NN_SPACING_PX / SPOT_PITCH_UM
    as_diameter = SEURAT_RADIUS / px_per_um
    as_radius = 2.0 * SEURAT_RADIUS / px_per_um

    assert as_diameter < SPOT_PITCH_UM, "spots must fit inside their own pitch"
    assert as_radius > SPOT_PITCH_UM, (
        "if this stops being true the overlap argument no longer holds and the "
        "radius finding needs re-deriving")
    # And read as a diameter it lands on 10x's reference spot size.
    assert as_diameter == pytest.approx(REFERENCE_SPOT_UM, rel=0.01)


def test_shanuz_radius_is_half_seurats_on_purpose():
    """A 2x divergence that must not be 'fixed' into agreement."""
    from shanuz.spatial.visium import ScaleFactors, VisiumV2

    sf = ScaleFactors(spot=SEURAT_RADIUS, fiducial=144.5316799616869,
                      hires=0.17211704, lowres=0.051635113)
    fov = VisiumV2(scale_factors=sf)

    assert fov.radius() == SEURAT_RADIUS / 2.0
    assert SEURAT_RADIUS / fov.radius() == pytest.approx(2.0)


def test_shanuz_exposes_a_radius_where_seurats_visiumv2_returns_null():
    from shanuz.spatial.visium import ScaleFactors, VisiumV2

    assert SEURAT_RADIUS_ON_VISIUMV2_IS_NULL, "Seurat 5.5.1 has no Radius.VisiumV2"
    fov = VisiumV2(scale_factors=ScaleFactors(spot=10.0, fiducial=20.0,
                                              hires=0.1, lowres=0.03))
    assert fov.radius() == 5.0          # reachable on the class, not just the centroids


def test_radius_is_none_rather_than_a_guess_without_scale_factors():
    from shanuz.spatial.visium import VisiumV2

    assert VisiumV2().radius() is None


# ---------------------------------------------------------------------------
# The defaults, which changed to match Seurat
# ---------------------------------------------------------------------------

def test_load_visium_defaults_are_seurats():
    import inspect

    from shanuz import load_visium

    p = inspect.signature(load_visium).parameters
    assert p["filter_by_tissue"].default is True       # Read10X_Image filter.matrix
    assert p["image_resolution"].default == "lowres"   # tissue_lowres_image.png
    assert p["slice_name"].default == "slice1"         # Load10X_Spatial's key


def test_tissue_coordinates_have_seurats_columns():
    """R: GetTissueCoordinates(image) returns x, y, cell — and the rownames are
    the cells as well. Keeping cell only on the index is the pandas idiom, but
    it is not the same frame, and `coords$cell` ported from R would fail on it.
    """
    import pandas as pd

    from shanuz import create_centroids

    cells = ["a", "b", "c"]
    cen = create_centroids(pd.DataFrame({"x": [1.0, 2.0, 3.0],
                                         "y": [4.0, 5.0, 6.0],
                                         "cell": cells}))
    df = cen.get_tissue_coordinates()
    assert list(df.columns) == ["x", "y", "cell"]
    assert list(df["cell"]) == cells
    assert list(df.index) == cells          # rownames too, as in R
    # x/y stay numeric — a caller taking the whole frame must not silently get
    # an object array.
    assert df[["x", "y"]].to_numpy().dtype.kind == "f"


# ---------------------------------------------------------------------------
# The comparison instrument itself
# ---------------------------------------------------------------------------

def test_every_tolerance_is_justified_and_bounded():
    t = _tutorial()
    assert t.FLOAT_TOLERANCES, "the tolerance table must not be empty"
    for name, tol in t.FLOAT_TOLERANCES.items():
        assert 0.0 <= tol <= 1e-5, f"{name} tolerance {tol} is too loose to mean anything"


def test_reported_anchors_are_never_counted_as_matches():
    t = _tutorial()
    fake = {n: 1.0 for n in t.REPORTED_ONLY}
    other = {n: 999.0 for n in t.REPORTED_ONLY}
    matched, differed, reported = t.compare(fake, other)
    assert not matched and not differed
    assert len(reported) == len(t.REPORTED_ONLY)


def test_a_reported_anchor_is_not_also_given_a_tolerance():
    t = _tutorial()
    overlap = set(t.REPORTED_ONLY) & set(t.FLOAT_TOLERANCES)
    assert not overlap, f"{overlap} is both excused and toleranced"


def test_exact_anchors_use_exact_comparison():
    """tol == 0.0 must mean array_equal, not 'approximately'."""
    t = _tutorial()
    assert t._match([1.0, 2.0], [1.0, 2.0], 0.0)
    assert not t._match([1.0, 2.0], [1.0, 2.0 + 1e-15], 0.0)


def test_the_comparator_can_actually_fail():
    """Guard the guard: a comparator that always matched would make every
    number this tutorial prints worthless."""
    t = _tutorial()
    assert not t._match([1.0], [2.0], 1e-6)
    assert not t._match([1.0, 2.0], [1.0], 0.0)      # shape mismatch
    assert not t._match("a", "b", 0.0)


# ---------------------------------------------------------------------------
# The recorded run, when it is present
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not (FIGURES / "r_anchors.json").exists(),
                    reason="figures_visium/r_anchors.json not generated (needs R)")
def test_recorded_run_matches_every_compared_anchor():
    t = _tutorial()
    py = json.loads((FIGURES / "py_anchors.json").read_text())
    r = json.loads((FIGURES / "r_anchors.json").read_text())
    r = {k: (v[0] if isinstance(v, list) and len(v) == 1 else v) for k, v in r.items()}
    matched, differed, _ = t.compare(py, r)
    assert not differed, [d[0] for d in differed]
    assert len(matched) >= 20


@pytest.mark.skipif(not (FIGURES / "r_anchors.json").exists(),
                    reason="figures_visium/r_anchors.json not generated (needs R)")
def test_recorded_run_reproduces_the_transcribed_r_constants():
    """If R and the constants above ever disagree, one of them is stale."""
    r = json.loads((FIGURES / "r_anchors.json").read_text())
    r = {k: (v[0] if isinstance(v, list) and len(v) == 1 else v) for k, v in r.items()}
    assert r["load.n_spots"] == SEURAT_N_SPOTS
    assert r["load.n_spots_unfiltered"] == SEURAT_N_SPOTS_UNFILTERED
    assert r["obj.n_features"] == SEURAT_N_FEATURES
    assert r["radius.centroids"] == SEURAT_RADIUS
    assert bool(r["radius.visium_is_null"]) is SEURAT_RADIUS_ON_VISIUMV2_IS_NULL
    assert bool(r["radius.has_visiumv2_method"]) is False


@pytest.mark.skipif(not (FIGURES / "py_anchors.json").exists(),
                    reason="figures_visium/py_anchors.json not generated")
def test_recorded_geometry_supports_the_overlap_argument():
    py = json.loads((FIGURES / "py_anchors.json").read_text())
    assert py["geometry.nn_spacing_px"] == pytest.approx(MEASURED_NN_SPACING_PX)
    assert py["geometry.spot_um_if_diameter"] < SPOT_PITCH_UM
    assert py["geometry.spot_um_if_radius"] > SPOT_PITCH_UM
    assert py["geometry.spot_um_if_diameter"] == pytest.approx(REFERENCE_SPOT_UM, rel=0.01)
    assert np.isclose(py["radius.centroids"], SEURAT_RADIUS / 2.0)
