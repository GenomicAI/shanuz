"""Spatial technology loaders — Xenium / Visium / CosMx.

Mirror Seurat's ``LoadXenium`` / ``Load10X_Spatial`` / ``LoadNanostring``: read a
platform's on-disk output into a Shanuz object with the expression assay AND a
populated ``seurat.images`` (per-FOV centroids), so the spatial accessors and
``shanuz.spatial.analysis`` functions work immediately.
"""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..io import read_10x
from .fov import create_fovs


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _read_table(path: Path) -> pd.DataFrame:
    """Read a cell/metadata table (.parquet, .csv, or .csv.gz)."""
    s = str(path)
    if s.endswith(".parquet"):
        return pd.read_parquet(path)
    return pd.read_csv(path)


def _feature_types(mtx_dir: Path) -> Optional[list[str]]:
    """Read the 3rd column (feature type) of a 10x features.tsv[.gz], if present.

    Xenium/Visium feature tables tag every row with a type — ``Gene Expression``
    for real genes and ``Negative Control Probe`` / ``Negative Control
    Codeword`` / ``Blank Codeword`` / ``Deprecated Codeword`` for QC controls.
    Returns one type per matrix row (file order), or ``None`` if unavailable.
    """
    feat = _first_existing(mtx_dir, ["features.tsv.gz", "features.tsv"])
    if feat is None:
        return None
    opn = gzip.open if str(feat).endswith(".gz") else open
    types: list[str] = []
    with opn(feat, "rt", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            types.append(parts[2] if len(parts) > 2 else "Gene Expression")
    return types


def _first_existing(base: Path, names: list[str]) -> Optional[Path]:
    for n in names:
        p = base / n
        if p.exists():
            return p
    return None


def _cell_id_column(df: pd.DataFrame) -> str:
    """Name of the cell-identifier column.

    MERSCOPE tables key on an unnamed leading index column in some exports and on
    an explicit ``cell``/``EntityID`` column in others; fall back to the first
    column, which is the cell id in every Vizgen layout.
    """
    for c in ("cell", "cell_id", "cell_ID", "EntityID"):
        if c in df.columns:
            return c
    return str(df.columns[0])


def _build_spatial_object(
    counts: sp.spmatrix,
    feature_names: list[str],
    cell_names: list[str],
    coords: pd.DataFrame,
    assay: str,
    project: str,
    fov: Optional[Union[str, np.ndarray]] = None,
    meta_data: Optional[pd.DataFrame] = None,
):
    """Assemble a Shanuz object + images from expression + coordinate parts."""
    from ..shanuz import create_shanuz_object

    obj = create_shanuz_object(
        counts, assay=assay, project=project,
        feature_names=feature_names, cell_names=cell_names,
    )
    kept = obj.cell_names()
    coords = coords.set_index("cell").reindex(kept)
    coords["cell"] = kept
    coords = coords.dropna(subset=["x", "y"])
    fov_labels = coords[fov].to_numpy() if isinstance(fov, str) and fov in coords else None
    obj.images = create_fovs(coords[["x", "y", "cell"]], fov=fov_labels, assay=assay,
                             default_name=assay.lower())
    if meta_data is not None:
        md = meta_data.reindex(kept)
        for c in md.columns:
            if c not in obj.meta_data.columns:
                obj.meta_data[c] = md[c].values
    return obj


# ---------------------------------------------------------------------------
# Xenium
# ---------------------------------------------------------------------------

def load_xenium(
    path: Union[str, Path],
    assay: str = "Xenium",
    fov_column: Optional[str] = None,
    project: str = "Xenium",
    keep_controls: bool = False,
):
    """Load a 10x Xenium output bundle into a Shanuz object with images.

    Expects (from the Xenium output folder):
      * ``cell_feature_matrix/`` — 10x MTX triplet (barcodes/features/matrix)
      * ``cells.parquet`` or ``cells.csv[.gz]`` — with ``cell_id``,
        ``x_centroid``, ``y_centroid`` (and optionally ``fov`` / transcript QC)

    By default only ``Gene Expression`` features are kept in the assay (matching
    Seurat's ``LoadXenium``, which routes Negative Control / Blank codewords to
    separate assays); set ``keep_controls=True`` to retain every feature row.

    ``fov_column``, if given and present in the cells table, splits the object
    into one image per FOV; otherwise a single image is created.
    """
    path = Path(path)
    mtx_dir = path / "cell_feature_matrix"
    if not mtx_dir.exists():
        raise FileNotFoundError(
            f"{mtx_dir} not found. load_xenium expects the unpacked "
            "cell_feature_matrix/ MTX directory."
        )
    counts, feats, cells = read_10x(mtx_dir)

    if not keep_controls:
        types = _feature_types(mtx_dir)
        if types is not None and len(types) == len(feats) and "Gene Expression" in types:
            mask = np.array([t == "Gene Expression" for t in types])
            counts = counts[mask, :]
            feats = [f for f, m in zip(feats, mask) if m]

    cell_candidates = [path / n for n in ("cells.parquet", "cells.csv.gz", "cells.csv")
                       if (path / n).exists()]
    if not cell_candidates:
        raise FileNotFoundError(f"No cells.parquet/csv found in {path}.")
    cdf = None
    for cf in cell_candidates:                      # prefer parquet, fall back to csv
        try:
            cdf = _read_table(cf)
            break
        except ImportError:                         # no parquet engine → try next
            continue
    if cdf is None:
        raise ImportError(
            "cells.parquet found but no parquet engine is installed. Install "
            "pyarrow, or provide cells.csv[.gz] alongside it."
        )
    rename = {"cell_id": "cell", "x_centroid": "x", "y_centroid": "y"}
    cdf = cdf.rename(columns={k: v for k, v in rename.items() if k in cdf.columns})
    if not {"cell", "x", "y"} <= set(cdf.columns):
        raise ValueError("cells table must contain cell_id, x_centroid, y_centroid.")
    cdf["cell"] = cdf["cell"].astype(str)

    fov = fov_column if (fov_column and fov_column in cdf.columns) else None
    coords = cdf[["cell", "x", "y"] + ([fov] if fov else [])]
    meta = cdf.set_index("cell").drop(columns=["x", "y"], errors="ignore")
    return _build_spatial_object(counts, feats, [str(c) for c in cells], coords,
                                 assay, project, fov=fov, meta_data=meta)


# ---------------------------------------------------------------------------
# Visium
# ---------------------------------------------------------------------------

def load_visium(
    path: Union[str, Path],
    assay: str = "Spatial",
    project: str = "Visium",
):
    """Load a 10x Visium output into a Shanuz object with spot coordinates.

    Expects:
      * ``filtered_feature_bc_matrix/`` — 10x MTX triplet
      * ``spatial/tissue_positions.csv`` (or ``tissue_positions_list.csv``) with
        barcode, in_tissue, array row/col and pixel row/col columns
    """
    path = Path(path)
    mtx_dir = _first_existing(path, ["filtered_feature_bc_matrix", "raw_feature_bc_matrix"])
    if mtx_dir is None:
        raise FileNotFoundError(f"No filtered_feature_bc_matrix/ in {path}.")
    counts, feats, cells = read_10x(mtx_dir)

    pos_file = _first_existing(path / "spatial",
                               ["tissue_positions.csv", "tissue_positions_list.csv"])
    if pos_file is None:
        raise FileNotFoundError(f"No spatial/tissue_positions.csv in {path}.")
    header = 0 if pos_file.name == "tissue_positions.csv" else None
    pos = pd.read_csv(pos_file, header=header)
    if header is None:
        pos.columns = ["barcode", "in_tissue", "array_row", "array_col",
                       "pxl_row_in_fullres", "pxl_col_in_fullres"]
    pos = pos.rename(columns={"barcode": "cell", "pxl_col_in_fullres": "x",
                              "pxl_row_in_fullres": "y"})
    pos["cell"] = pos["cell"].astype(str)
    coords = pos[["cell", "x", "y"]]
    return _build_spatial_object(counts, feats, [str(c) for c in cells], coords,
                                 assay, project)


# ---------------------------------------------------------------------------
# CosMx / Nanostring
# ---------------------------------------------------------------------------

def load_cosmx(
    path: Union[str, Path],
    expr_file: Optional[str] = None,
    meta_file: Optional[str] = None,
    assay: str = "Nanostring",
    fov_column: str = "fov",
    project: str = "CosMx",
):
    """Load NanoString CosMx output (exprMat + metadata CSVs) into a Shanuz object.

    ``expr_file`` is a cell×gene CSV (rows = cells, first columns cell/fov ids);
    ``meta_file`` carries ``CenterX_global_px`` / ``CenterY_global_px`` and a FOV
    column. If names are omitted, ``*exprMat_file.csv`` / ``*metadata_file.csv``
    are auto-detected in ``path``.
    """
    path = Path(path)
    expr = Path(expr_file) if expr_file else next(iter(path.glob("*exprMat_file.csv")), None)
    meta = Path(meta_file) if meta_file else next(iter(path.glob("*metadata_file.csv")), None)
    if expr is None or meta is None:
        raise FileNotFoundError("Could not locate exprMat_file.csv / metadata_file.csv.")

    edf = pd.read_csv(expr)
    mdf = pd.read_csv(meta)
    id_cols = [c for c in ("fov", "cell_ID", "cell_id", "cell") if c in edf.columns]
    gene_cols = [c for c in edf.columns if c not in id_cols]

    def _cid(df):
        fovc = fov_column if fov_column in df.columns else id_cols[0]
        cidc = next((c for c in ("cell_ID", "cell_id", "cell") if c in df.columns), None)
        return (df[fovc].astype(str) + "_" + df[cidc].astype(str)).to_numpy()

    cell_ids = _cid(edf)
    counts = sp.csc_matrix(edf[gene_cols].to_numpy(dtype=float).T)   # genes × cells

    mdf = mdf.copy()
    mdf["cell"] = _cid(mdf)
    mcoord = mdf.rename(columns={"CenterX_global_px": "x", "CenterY_global_px": "y"})
    coords = mcoord[["cell", "x", "y"] + ([fov_column] if fov_column in mcoord else [])]
    return _build_spatial_object(counts, gene_cols, list(cell_ids), coords,
                                 assay, project,
                                 fov=fov_column if fov_column in coords else None,
                                 meta_data=mdf.set_index("cell"))


# ---------------------------------------------------------------------------
# MERSCOPE / Vizgen
# ---------------------------------------------------------------------------

def load_merscope(
    path: Union[str, Path],
    expr_file: Optional[str] = None,
    meta_file: Optional[str] = None,
    assay: str = "Vizgen",
    fov_column: str = "fov",
    project: str = "MERSCOPE",
    keep_controls: bool = False,
):
    """Load a Vizgen MERSCOPE output into a Shanuz object with images.

    Mirrors Seurat's ``LoadVizgen``. Expects, in ``path``:
      * ``cell_by_gene.csv`` — cell × gene counts (leading column = cell id)
      * ``cell_metadata.csv`` — with ``center_x`` / ``center_y`` (and usually
        ``fov``, ``volume``)

    Blank/control barcodes (``Blank-*`` columns) are dropped by default, matching
    ``LoadVizgen``; set ``keep_controls=True`` to retain them.

    ``fov_column``, if present in the metadata, splits the object into one image
    per FOV; otherwise a single image is created.
    """
    path = Path(path)
    expr = Path(expr_file) if expr_file else _first_existing(
        path, ["cell_by_gene.csv", "cell_by_gene.csv.gz"])
    meta = Path(meta_file) if meta_file else _first_existing(
        path, ["cell_metadata.csv", "cell_metadata.csv.gz"])
    if expr is None or meta is None:
        raise FileNotFoundError(
            f"Could not locate cell_by_gene.csv / cell_metadata.csv in {path}."
        )

    edf = pd.read_csv(expr)
    mdf = pd.read_csv(meta)

    ecid = _cell_id_column(edf)
    gene_cols = [c for c in edf.columns if c != ecid]
    if not keep_controls:
        gene_cols = [g for g in gene_cols if not str(g).lower().startswith("blank")]
    if not gene_cols:
        raise ValueError(f"No gene columns found in {expr}.")
    cell_ids = edf[ecid].astype(str).to_numpy()
    counts = sp.csc_matrix(edf[gene_cols].to_numpy(dtype=float).T)   # genes × cells

    mdf = mdf.copy()
    mdf["cell"] = mdf[_cell_id_column(mdf)].astype(str)
    mcoord = mdf.rename(columns={"center_x": "x", "center_y": "y"})
    if not {"x", "y"} <= set(mcoord.columns):
        raise ValueError("cell_metadata must contain center_x / center_y columns.")
    fov = fov_column if fov_column in mcoord.columns else None
    coords = mcoord[["cell", "x", "y"] + ([fov] if fov else [])]
    return _build_spatial_object(counts, gene_cols, list(cell_ids), coords,
                                 assay, project, fov=fov,
                                 meta_data=mdf.set_index("cell"))
