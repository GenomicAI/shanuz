from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from .base import SpatialImage


def _close_rings(coords: pd.DataFrame) -> pd.DataFrame:
    """Repeat each cell's first vertex at the end of its ring, as R does.

    ``CreateSegmentation`` stores closed rings — a square arrives as four
    vertices and comes back as five. That is the ``sf``/GEOS convention R
    inherits, and it is not decoration: code that reads the vertex list to
    measure a perimeter, or to draw an outline without asking matplotlib to close
    it, is off by one edge on an open ring.

    Already-closed rings are left alone, so this is idempotent.
    """
    if coords.empty:
        return coords
    out = []
    for _, ring in coords.groupby("cell", sort=False):
        first, last = ring.iloc[0], ring.iloc[-1]
        if len(ring) > 2 and (first["x"], first["y"]) != (last["x"], last["y"]):
            ring = pd.concat([ring, ring.iloc[[0]]], ignore_index=False)
        out.append(ring)
    return pd.concat(out) if out else coords


class Segmentation(SpatialImage):
    """Cell boundary polygon coordinates.

    Mirrors R's Segmentation class from segmentation.R.
    Each cell may have multiple (x, y) polygon vertices.

    Slots
    -----
    _coords : pd.DataFrame  columns: x, y, cell   (multiple rows per cell)
    """

    __slots__ = ("_coords", "assay", "misc", "_key")

    def __init__(
        self,
        coords: pd.DataFrame,
        assay: str = "",
        key: str = "segmentation_",
        misc: Optional[dict] = None,
    ) -> None:
        super().__init__(assay=assay, key=key, misc=misc)
        for col in ("x", "y", "cell"):
            if col not in coords.columns:
                raise ValueError(f"coords must have a '{col}' column.")
        self._coords = _close_rings(coords[["x", "y", "cell"]].copy())

    # ------------------------------------------------------------------
    # SpatialImage interface
    # ------------------------------------------------------------------

    def cells(self) -> list[str]:
        return list(self._coords["cell"].unique())

    def dim(self) -> tuple[int, int]:
        return (self._coords["cell"].nunique(), 2)

    def get_tissue_coordinates(
        self,
        cells: Optional[list[str]] = None,
    ) -> pd.DataFrame:
        df = self._coords.set_index("cell")[["x", "y"]]
        if cells is not None:
            df = df.loc[df.index.isin(cells)]
        return df.copy()

    def rename_cells(self, new_names: list[str]) -> "Segmentation":
        old_names = self.cells()
        if len(new_names) != len(old_names):
            raise ValueError("new_names length must match number of unique cells.")
        mapping = dict(zip(old_names, new_names))
        new_coords = self._coords.copy()
        new_coords["cell"] = new_coords["cell"].map(mapping)
        return Segmentation(
            coords=new_coords,
            assay=self.assay,
            key=self._key,
            misc=dict(self.misc),
        )

    def subset(self, cells: list[str]) -> "Segmentation":
        mask = self._coords["cell"].isin(cells)
        return Segmentation(
            coords=self._coords[mask].copy(),
            assay=self.assay,
            key=self._key,
            misc=dict(self.misc),
        )

    # ------------------------------------------------------------------
    # Simplify polygons
    # ------------------------------------------------------------------

    def simplify(self, tol: float = 0.5) -> "Segmentation":
        """Reduce polygon vertex count by removing vertices closer than tol.

        Simple Douglas–Peucker–style approximation per cell.
        """
        groups = []
        for cell_id, group in self._coords.groupby("cell", sort=False):
            pts = group[["x", "y"]].values
            if len(pts) <= 3:
                groups.append(group)
                continue
            # keep[i] = True means retain pts[i]; always keep first and last
            keep = np.zeros(len(pts), dtype=bool)
            keep[0] = True
            keep[-1] = True
            diffs = np.linalg.norm(np.diff(pts, axis=0), axis=1)  # len(pts)-1
            keep[1:-1] = diffs[:-1] > tol
            groups.append(group.iloc[np.where(keep)[0]])

        new_coords = pd.concat(groups, ignore_index=True)
        return Segmentation(
            coords=new_coords,
            assay=self.assay,
            key=self._key,
            misc=dict(self.misc),
        )


def create_segmentation(
    coords: pd.DataFrame,
    assay: str = "",
    key: str = "segmentation_",
) -> Segmentation:
    return Segmentation(coords=coords, assay=assay, key=key)
