"""Visium tissue image support — scale factors and the ``VisiumV2`` FOV.

Mirrors Seurat v5's ``VisiumV2`` class: an :class:`~shanuz.spatial.fov.FOV` that
additionally carries the H&E tissue image and the ``scalefactors_json.json``
values needed to place spots on it.

Coordinate convention
---------------------
Everything stored in the FOV — spot centroids and :meth:`VisiumV2.radius` — is in
**full-resolution pixel space**, the space ``tissue_positions.csv`` uses. The
scale factors convert that to the pixel space of the *stored* (hires or lowres)
image: multiply by :meth:`VisiumV2.scale_factor`, or just call
:meth:`VisiumV2.scale_coordinates`. Keeping the FOV in fullres means the spatial
analysis functions (``spatial_knn``, ``nearest_neighbor_distance``, Moran's I)
see real, image-independent distances.
"""
from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd

from .centroids import Centroids
from .fov import FOV
from .molecules import Molecules


@dataclass
class ScaleFactors:
    """The four values in a Visium ``spatial/scalefactors_json.json``.

    spot     : spot diameter in full-resolution pixels
    fiducial : fiducial-marker diameter in full-resolution pixels
    hires    : multiply a fullres coordinate by this to land on tissue_hires_image.png
    lowres   : ditto for tissue_lowres_image.png
    """

    spot: float
    fiducial: float
    hires: float
    lowres: float

    def scale_factor(self, resolution: str) -> float:
        if resolution == "hires":
            return self.hires
        if resolution == "lowres":
            return self.lowres
        raise ValueError(f"resolution must be 'hires' or 'lowres', got {resolution!r}.")


def read_scale_factors(path: Union[str, Path]) -> ScaleFactors:
    """Read a ``scalefactors_json.json`` (the file itself, or its parent dir)."""
    path = Path(path)
    if path.is_dir():
        path = path / "scalefactors_json.json"
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    return ScaleFactors(
        spot=float(d.get("spot_diameter_fullres", np.nan)),
        fiducial=float(d.get("fiducial_diameter_fullres", np.nan)),
        hires=float(d.get("tissue_hires_scalef", np.nan)),
        lowres=float(d.get("tissue_lowres_scalef", np.nan)),
    )


def read_tissue_image(
    path: Union[str, Path],
    resolution: str = "lowres",
) -> Optional[tuple[np.ndarray, str]]:
    """Read ``tissue_{hires,lowres}_image.png`` from a Visium ``spatial/`` dir.

    Falls back to the other resolution when the requested one is missing, so a
    lowres-only bundle still yields an image. Returns ``(image, resolution)`` —
    the resolution actually read, which may differ from the one asked for — or
    ``None`` when no image file exists.

    Reading a PNG needs matplotlib or Pillow; neither is a core dependency, so if
    both are missing this warns and returns ``None`` rather than failing the load.
    """
    if resolution not in ("hires", "lowres"):
        raise ValueError(f"resolution must be 'hires' or 'lowres', got {resolution!r}.")
    path = Path(path)
    order = [resolution] + [r for r in ("hires", "lowres") if r != resolution]
    for res in order:
        png = path / f"tissue_{res}_image.png"
        if png.exists():
            img = _imread(png)
            return None if img is None else (img, res)
    return None


def _to_unit_float(img: np.ndarray) -> np.ndarray:
    """Integer pixel values → float in [0, 1]; an already-float image passes through."""
    if np.issubdtype(img.dtype, np.integer):
        return img.astype(np.float32) / float(np.iinfo(img.dtype).max)
    return img


def _imread(png: Path) -> Optional[np.ndarray]:
    """Read a PNG as float in [0, 1], whichever backend happens to be installed.

    The two backends do not agree on their own: matplotlib returns float in
    [0, 1] for an 8-bit PNG, Pillow returns uint8 in [0, 255]. Both are "the
    image", but they are not the same array, so without normalising here
    ``get_image()`` would be a function of the environment as much as of the
    file — and neither matplotlib nor Pillow is a declared dependency. R's
    ``png::readPNG`` is always double in [0, 1]; match that.
    """
    img = None
    try:
        import matplotlib.image as mpimg

        img = np.asarray(mpimg.imread(png))
    except ImportError:
        try:
            from PIL import Image

            with Image.open(png) as im:
                # matplotlib expands a palette PNG to RGB(A); Pillow hands back
                # palette indices unless asked, which is the same disagreement.
                if im.mode == "P":
                    im = im.convert("RGBA" if "transparency" in im.info else "RGB")
                img = np.asarray(im)
        except ImportError:
            warnings.warn(
                "Reading the Visium tissue image needs matplotlib or Pillow; "
                "loading without it. Install shanuz[analysis] to enable tissue plots.",
                stacklevel=3,
            )
            return None
    return _to_unit_float(img)


class VisiumV2(FOV):
    """An FOV that also carries the Visium tissue image and its scale factors.

    Mirrors Seurat v5's ``VisiumV2``. Behaves as a normal FOV everywhere (its
    coordinates stay in fullres pixels); the extra slots are what
    ``spatial_dim_plot`` / ``spatial_feature_plot`` need to draw spots on tissue.

    Slots
    -----
    image            : (H, W[, C]) array, or None
    scale_factors    : ScaleFactors, or None
    image_resolution : 'hires' | 'lowres' — which image is stored
    """

    __slots__ = ("image", "scale_factors", "image_resolution")

    def __init__(
        self,
        boundaries: Optional[dict] = None,
        molecules: Optional[dict[str, Molecules]] = None,
        image: Optional[np.ndarray] = None,
        scale_factors: Optional[ScaleFactors] = None,
        image_resolution: str = "lowres",
        coords_x_orientation: str = "horizontal",
        assay: str = "",
        key: str = "slice1_",
        misc: Optional[dict] = None,
    ) -> None:
        super().__init__(
            boundaries=boundaries,
            molecules=molecules,
            coords_x_orientation=coords_x_orientation,
            assay=assay,
            key=key,
            misc=misc,
        )
        self.image = image
        self.scale_factors = scale_factors
        self.image_resolution = image_resolution

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_fov(
        cls,
        fov: FOV,
        image: Optional[np.ndarray] = None,
        scale_factors: Optional[ScaleFactors] = None,
        image_resolution: str = "lowres",
    ) -> "VisiumV2":
        """Upgrade a plain FOV in place-ish: same boundaries, plus the image."""
        obj = cls(
            boundaries=dict(fov.boundaries),
            molecules=dict(fov.molecules),
            image=image,
            scale_factors=scale_factors,
            image_resolution=image_resolution,
            coords_x_orientation=fov.coords_x_orientation,
            assay=fov.assay,
            key=fov._key,
            misc=dict(fov.misc),
        )
        # Spot boundaries know their own radius, in the same fullres pixel space.
        r = obj.radius()
        if r is not None:
            for b in obj.boundaries.values():
                if isinstance(b, Centroids):
                    b.radius_ = r
        return obj

    # ------------------------------------------------------------------
    # SpatialImage interface
    # ------------------------------------------------------------------

    def get_image(self) -> Optional[np.ndarray]:
        """The stored tissue image, or None when the bundle had none."""
        return self.image

    def radius(self) -> Optional[float]:
        """Spot radius in full-resolution pixels (half the spot diameter)."""
        if self.scale_factors is None or not np.isfinite(self.scale_factors.spot):
            return None
        return float(self.scale_factors.spot) / 2.0

    def rename_cells(self, new_names: list[str]) -> "VisiumV2":
        return self._with(super().rename_cells(new_names))

    def subset(self, cells: list[str]) -> "VisiumV2":
        return self._with(super().subset(cells))

    def _with(self, fov: FOV) -> "VisiumV2":
        """Re-attach image + scale factors to an FOV produced by the base class."""
        return VisiumV2(
            boundaries=fov.boundaries,
            molecules=fov.molecules,
            image=self.image,
            scale_factors=self.scale_factors,
            image_resolution=self.image_resolution,
            coords_x_orientation=fov.coords_x_orientation,
            assay=fov.assay,
            key=fov._key,
            misc=dict(fov.misc),
        )

    # ------------------------------------------------------------------
    # Image-space helpers
    # ------------------------------------------------------------------

    def scale_factor(self, resolution: Optional[str] = None) -> float:
        """Fullres → image-pixel multiplier for the stored (or given) resolution."""
        if self.scale_factors is None:
            return 1.0
        return self.scale_factors.scale_factor(resolution or self.image_resolution)

    def scale_coordinates(
        self,
        cells: Optional[list[str]] = None,
        resolution: Optional[str] = None,
    ) -> pd.DataFrame:
        """Tissue coordinates rescaled into the pixel space of the stored image.

        Same frame as :meth:`get_tissue_coordinates`, with ``x``/``y`` multiplied
        by the scale factor — i.e. ready to overlay on :meth:`get_image`.
        """
        coords = self.get_tissue_coordinates(cells=cells).copy()
        f = self.scale_factor(resolution)
        for col in ("x", "y"):
            if col in coords.columns:
                coords[col] = coords[col].astype(float) * f
        return coords

    def spot_radius(self, resolution: Optional[str] = None) -> Optional[float]:
        """Spot radius in the pixel space of the stored image (None if unknown)."""
        r = self.radius()
        return None if r is None else r * self.scale_factor(resolution)

    def __repr__(self) -> str:
        if self.image is None:
            img = "no image"
        else:
            h, w = self.image.shape[:2]
            img = f"{self.image_resolution} image {w}×{h}"
        return (
            f"VisiumV2: {len(self.cells())} spots  Key: {self._key!r}  "
            f"Assay: {self.assay!r}\n"
            f"  Boundaries: {list(self.boundaries)}\n"
            f"  Image:      {img}"
        )
