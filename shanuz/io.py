"""10X Genomics and other single-cell data I/O utilities.

Mirrors R's Read10X() from Seurat.
"""
from __future__ import annotations

import gzip
import os
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def read_10x(
    data_dir: Union[str, Path],
    var_names: str = "gene_symbols",
    make_unique: bool = True,
) -> tuple[sp.csc_matrix, list[str], list[str]]:
    """Read 10X Genomics output directory.

    Mirrors R's Read10X(). Supports both v2 (genes.tsv) and v3
    (features.tsv.gz) directory layouts.

    Parameters
    ----------
    data_dir    : path to the 10X output directory
    var_names   : 'gene_symbols' (default) or 'gene_ids'
    make_unique : append suffix to duplicate gene names

    Returns
    -------
    (matrix, feature_names, cell_names)
    matrix is (features × cells) csc_matrix.
    """
    data_dir = Path(data_dir)

    # Detect v2 vs v3 layout
    if (data_dir / "features.tsv.gz").exists():
        matrix_file = data_dir / "matrix.mtx.gz"
        barcodes_file = data_dir / "barcodes.tsv.gz"
        features_file = data_dir / "features.tsv.gz"
        v3 = True
    elif (data_dir / "genes.tsv").exists():
        matrix_file = data_dir / "matrix.mtx"
        barcodes_file = data_dir / "barcodes.tsv"
        features_file = data_dir / "genes.tsv"
        v3 = False
    elif (data_dir / "features.tsv").exists():
        matrix_file = data_dir / "matrix.mtx"
        barcodes_file = data_dir / "barcodes.tsv"
        features_file = data_dir / "features.tsv"
        v3 = True
    else:
        raise FileNotFoundError(
            f"No 10X matrix files found in {data_dir}. "
            "Expected genes.tsv/features.tsv[.gz], barcodes.tsv[.gz], matrix.mtx[.gz]."
        )

    # Read matrix
    mat = _read_mtx(matrix_file)
    mat = sp.csc_matrix(mat)

    # Read barcodes
    cell_names = _read_tsv_column(barcodes_file, col=0)

    # Read features
    gene_ids = _read_tsv_column(features_file, col=0)
    gene_symbols = _read_tsv_column(features_file, col=1)

    feature_names = gene_symbols if var_names == "gene_symbols" else gene_ids

    if make_unique:
        feature_names = _make_unique(feature_names)

    # Validate shapes
    if mat.shape[0] != len(feature_names):
        raise ValueError(
            f"Matrix rows ({mat.shape[0]}) != feature count ({len(feature_names)})."
        )
    if mat.shape[1] != len(cell_names):
        raise ValueError(
            f"Matrix cols ({mat.shape[1]}) != barcode count ({len(cell_names)})."
        )

    return mat, feature_names, cell_names


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def _open_file(path: Path):
    """Open plain or gzipped file."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _read_mtx(path: Path) -> sp.coo_matrix:
    """Read a (possibly gzipped) Matrix Market file."""
    if str(path).endswith(".gz"):
        import tempfile, shutil
        with gzip.open(path, "rb") as src:
            with tempfile.NamedTemporaryFile(suffix=".mtx", delete=False) as dst:
                shutil.copyfileobj(src, dst)
                tmp_path = dst.name
        mat = scipy.io.mmread(tmp_path)
        os.unlink(tmp_path)
        return mat
    return scipy.io.mmread(str(path))


def _read_tsv_column(path: Path, col: int = 0) -> list[str]:
    """Read one column from a TSV file (plain or gzipped)."""
    values = []
    with _open_file(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            values.append(parts[col] if col < len(parts) else "")
    return values


def _make_unique(names: list[str]) -> list[str]:
    """Append .1, .2 … suffixes to duplicate names (mirrors make.unique in R)."""
    seen: dict[str, int] = {}
    result = []
    for n in names:
        if n in seen:
            seen[n] += 1
            result.append(f"{n}.{seen[n]}")
        else:
            seen[n] = 0
            result.append(n)
    return result
