"""Generate the figures for the Visium tutorial.

Renders to tutorials/figures_visium/:
  * py_01_radius_geometry.png  the slide geometry that decides radius vs diameter
  * py_02_spots_on_tissue.png  the loaded object drawn on its own H&E image
  * py_03_pca_isolation.png    the PCA gap, traced to feature selection

01 is the tutorial's finding, drawn: Seurat passes `spot_diameter_fullres` into
`CreateFOV`'s `radius`, and read that way the capture spots overlap. 02 is the
container working end to end — coordinates, scale factors and image agreeing
well enough to land 2,695 spots on the tissue. 03 separates shanuz's PCA from
the feature selection feeding it.

R-side numbers are constants read off a live Seurat 5.5.1 session, so these
render without R installed. `Rscript tutorials/visium_verify.R` reproduces them.

Usage
-----
    python tutorials/generate_visium_plots.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np  # noqa: E402

_ROOT = Path(__file__).parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

FIGURES = Path(__file__).parent / "figures_visium"

_SHANUZ = "#48a9a6"
_SEURAT = "#c1666b"
_GREY = "0.35"

# R: Seurat 5.5.1, Read10X_Image on the V1_Mouse_Brain_Sagittal_Anterior bundle.
SEURAT_RADIUS = 89.47199235723474      # == spot_diameter_fullres
# Measured on this slide by tutorials/visium_tutorial.py.
NN_SPACING_PX = 137.0
PITCH_UM = 100.0                        # fixed by the Visium slide
REFERENCE_SPOT_UM = 65.0                # what 10x defines the field against
# PCA stdev, max relative difference vs Seurat, from `--report`.
PCA_REL_OWN = 0.00194
PCA_REL_SHARED = 2.49e-05
HVG_OVERLAP, HVG_TOTAL = 1995, 2000


def _save(fig, name):
    import matplotlib.pyplot as plt

    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / name, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {name}")


def radius_geometry():
    """Draw both readings of spot_diameter_fullres on the real spot lattice."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    px_per_um = NN_SPACING_PX / PITCH_UM
    as_diameter = SEURAT_RADIUS / px_per_um          # um
    as_radius = 2.0 * SEURAT_RADIUS / px_per_um      # um

    # A patch of the Visium lattice: 100 um centre-to-centre, rows offset by half.
    centres = [(col * PITCH_UM + (row % 2) * PITCH_UM / 2, row * PITCH_UM * np.sqrt(3) / 2)
               for row in range(4) for col in range(4)]

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.0))
    for ax, dia, label, colour in (
        (axes[0], as_diameter, "read as a DIAMETER — shanuz", _SHANUZ),
        (axes[1], as_radius, "read as a RADIUS — Seurat", _SEURAT),
    ):
        for cx, cy in centres:
            ax.add_patch(Circle((cx, cy), dia / 2, facecolor=colour, alpha=0.45,
                                edgecolor=colour, lw=1.0))
            ax.plot(cx, cy, ".", color="0.2", ms=2.5, zorder=3)
        ax.set_aspect("equal")
        ax.set_xlim(-80, 430)
        ax.set_ylim(-80, 350)
        ax.set_title(f"{label}\nspot {dia:.1f} µm on a {PITCH_UM:.0f} µm pitch",
                     fontsize=10.5)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)

    axes[0].text(0.5, -0.06, f"fits — {as_diameter / PITCH_UM * 100:.0f}% of the pitch",
                 transform=axes[0].transAxes, ha="center", fontsize=9.5, color=_SHANUZ)
    axes[1].text(0.5, -0.06,
                 f"overlaps by {as_radius - PITCH_UM:.0f} µm — physically impossible",
                 transform=axes[1].transAxes, ha="center", fontsize=9.5, color=_SEURAT)

    fig.suptitle("Seurat stores spot_diameter_fullres in the radius slot", fontsize=12)
    fig.text(0.5, -0.03,
             f"Spot spacing measured at {NN_SPACING_PX:.0f} px over the Visium slide's "
             f"fixed {PITCH_UM:.0f} µm pitch, so {NN_SPACING_PX / PITCH_UM:.2f} px/µm. "
             f"Read as a diameter the field lands on 10x's {REFERENCE_SPOT_UM:.0f} µm\n"
             f"reference spot ({as_diameter / REFERENCE_SPOT_UM:.4f}×); read as a radius "
             f"the capture areas would run into each other. Neither tool is consulted "
             f"for this — it is the slide.",
             ha="center", fontsize=8.5, color=_GREY)
    return fig


def spots_on_tissue(data_dir=None):
    """2,695 spots placed on the H&E image, straight out of load_visium."""
    import matplotlib.pyplot as plt

    from shanuz import load_visium
    from shanuz.datasets import visium_mouse_brain

    obj = load_visium(data_dir or visium_mouse_brain())
    fov = obj.images["slice1"]
    img = fov.get_image()
    coords = fov.scale_coordinates()

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 5.2))
    axes[0].imshow(img)
    axes[0].set_title(f"tissue_lowres_image.png\n{img.shape[1]}×{img.shape[0]}, "
                      f"float {img.min():.2f}–{img.max():.2f}", fontsize=10)
    axes[1].imshow(img, alpha=0.55)
    axes[1].scatter(coords["x"], coords["y"], s=3.0, c=_SHANUZ, linewidths=0)
    axes[1].set_title(f"{len(coords):,} in-tissue spots, scaled by "
                      f"tissue_lowres_scalef", fontsize=10)
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle("The container works: coordinates, scale factors and image agree",
                 fontsize=12)
    fig.text(0.5, -0.02,
             "Spot positions are stored in full-resolution pixels, as "
             "tissue_positions_list.csv gives them, and scaled on demand. Every "
             "coordinate matches Seurat exactly — max|dx| = max|dy| = 0 over all "
             "2,695 spots.",
             ha="center", fontsize=8.5, color=_GREY)
    return fig


def pca_isolation():
    """Is the PCA gap the PCA, or the features feeding it?"""
    import matplotlib.pyplot as plt

    labels = [f"shanuz's own\n{HVG_TOTAL} features",
              f"Seurat's\n{HVG_TOTAL} features"]
    values = [PCA_REL_OWN, PCA_REL_SHARED]

    fig, ax = plt.subplots(figsize=(7.4, 3.9))
    bars = ax.bar(range(2), values, color=[_SEURAT, _SHANUZ], width=0.5)
    ax.set_yscale("log")
    ax.set_xticks(range(2))
    ax.set_xticklabels(labels, fontsize=9.5)
    ax.set_ylabel("max relative difference in PCA stdev vs Seurat")
    ax.set_ylim(values[1] / 6, values[0] * 6)
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, v * 1.35, f"{v:.3g}",
                ha="center", fontsize=10)
    ax.set_title("Running shanuz's PCA on Seurat's feature list closes the gap "
                 f"{values[0] / values[1]:.0f}×", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(0.5, -0.10,
             f"The two tools select {HVG_OVERLAP} of the same {HVG_TOTAL} variable "
             f"features; the {HVG_TOTAL - HVG_OVERLAP} that differ trace to the LOESS "
             f"smoother, not to the decomposition.\nGiven the same input the "
             f"decompositions agree to {PCA_REL_SHARED:.3g}, so the gap is feature "
             f"selection — which is why the tutorial reports this anchor instead of "
             f"asserting a loose tolerance on it.",
             ha="center", fontsize=8.5, color=_GREY)
    return fig


def main(data_dir=None):
    _save(radius_geometry(), "py_01_radius_geometry.png")
    _save(spots_on_tissue(data_dir), "py_02_spots_on_tissue.png")
    _save(pca_isolation(), "py_03_pca_isolation.png")
    anchors = FIGURES / "py_anchors.json"
    if anchors.exists():
        a = json.loads(anchors.read_text())
        assert np.isclose(a["geometry.nn_spacing_px"], NN_SPACING_PX), (
            "the measured spot spacing moved; the constants in this file are stale")
    print(f"\n  Wrote figures to {FIGURES}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Visium tutorial figures")
    p.add_argument("--data-dir", default=None)
    main(p.parse_args().data_dir)
