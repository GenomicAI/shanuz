from __future__ import annotations

from typing import Optional

import pandas as pd

from .base import SpatialImage


class Molecules(SpatialImage):
    """Spatially-resolved molecule (FISH) data.

    Mirrors R's Molecules class from molecules.R.
    Each row represents one molecule detection event.

    Slots
    -----
    _coords : pd.DataFrame  columns: x, y, gene, [cell]
                             'cell' is optional (not all FISH protocols assign cells)
    """

    __slots__ = ("_coords", "assay", "misc", "_key")

    def __init__(
        self,
        coords: pd.DataFrame,
        assay: str = "",
        key: str = "molecules_",
        misc: Optional[dict] = None,
    ) -> None:
        super().__init__(assay=assay, key=key, misc=misc)
        for col in ("x", "y", "gene"):
            if col not in coords.columns:
                raise ValueError(f"coords must have a '{col}' column.")
        keep = ["x", "y", "gene"]
        if "cell" in coords.columns:
            keep.append("cell")
        self._coords = coords[keep].copy()

    # ------------------------------------------------------------------
    # SpatialImage interface
    # ------------------------------------------------------------------

    def cells(self) -> list[str]:
        if "cell" not in self._coords.columns:
            return []
        return list(self._coords["cell"].dropna().unique())

    def features(self) -> list[str]:
        return list(self._coords["gene"].unique())

    def dim(self) -> tuple[int, int]:
        return (len(self._coords), 2)

    def get_tissue_coordinates(
        self,
        cells: Optional[list[str]] = None,
    ) -> pd.DataFrame:
        df = self._coords[["x", "y"]].copy()
        if cells is not None and "cell" in self._coords.columns:
            mask = self._coords["cell"].isin(cells)
            df = df[mask]
        return df

    def rename_cells(self, new_names: list[str]) -> "Molecules":
        if "cell" not in self._coords.columns:
            raise ValueError("This Molecules object has no 'cell' column to rename.")
        old_names = self.cells()
        if len(new_names) != len(old_names):
            raise ValueError("new_names length must match number of unique cells.")
        mapping = dict(zip(old_names, new_names))
        new_coords = self._coords.copy()
        new_coords["cell"] = new_coords["cell"].map(mapping)
        return Molecules(coords=new_coords, assay=self.assay, key=self._key, misc=dict(self.misc))

    def subset(self, cells: list[str]) -> "Molecules":
        if "cell" not in self._coords.columns:
            return Molecules(
                coords=self._coords.copy(), assay=self.assay, key=self._key, misc=dict(self.misc)
            )
        mask = self._coords["cell"].isin(cells)
        return Molecules(
            coords=self._coords[mask].copy(),
            assay=self.assay,
            key=self._key,
            misc=dict(self.misc),
        )


def create_molecules(
    coords: pd.DataFrame,
    assay: str = "",
    key: str = "molecules_",
) -> Molecules:
    return Molecules(coords=coords, assay=assay, key=key)
