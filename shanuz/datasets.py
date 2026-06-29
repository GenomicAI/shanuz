"""Built-in dataset loaders.

Provides automatic download of commonly used benchmark datasets.
"""
from __future__ import annotations

import os
import tarfile
import urllib.request
from pathlib import Path
from typing import Optional

import scipy.sparse as sp


# PBMC 3k dataset download sources (tried in order)
_PBMC3K_URLS = [
    "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
]
_PBMC3K_DIR = "filtered_gene_bc_matrices/hg19"

# PBMC 8k dataset (10x Genomics, GRCh38, v2 chemistry) — ~8,400 cells.
_PBMC8K_URLS = [
    "https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/"
    "pbmc8k_filtered_gene_bc_matrices.tar.gz",
]
_PBMC8K_DIR = "filtered_gene_bc_matrices/GRCh38"


def pbmc3k(
    data_dir: Optional[str] = None,
    force_download: bool = False,
) -> tuple[sp.csc_matrix, list[str], list[str]]:
    """Download (if needed) and load the PBMC 3k dataset.

    Returns (counts_matrix, gene_names, cell_barcodes).
    matrix is (genes × cells) csc_matrix with raw counts.

    Parameters
    ----------
    data_dir        : directory to cache the raw files
                      (defaults to ~/.shanuz_data/pbmc3k)
    force_download  : re-download even if files exist
    """
    from .io import read_10x

    if data_dir is None:
        data_dir = Path.home() / ".shanuz_data" / "pbmc3k"
    else:
        data_dir = Path(data_dir)

    matrix_dir = data_dir / _PBMC3K_DIR

    if force_download or not (matrix_dir / "matrix.mtx").exists():
        _download_10x(_PBMC3K_URLS, data_dir, label="PBMC3k", size_mb=24)

    mat, genes, cells = read_10x(matrix_dir, var_names="gene_symbols")
    return mat, genes, cells


def pbmc8k(
    data_dir: Optional[str] = None,
    force_download: bool = False,
) -> tuple[sp.csc_matrix, list[str], list[str]]:
    """Download (if needed) and load the 10x Genomics PBMC 8k dataset.

    ~8,400 peripheral blood mononuclear cells (GRCh38, v2 chemistry). Larger
    than :func:`pbmc3k` and used by the advanced subclustering tutorial.

    Returns (counts_matrix, gene_names, cell_barcodes); matrix is
    (genes x cells) csc_matrix with raw counts.

    Parameters
    ----------
    data_dir        : directory to cache the raw files
                      (defaults to ~/.shanuz_data/pbmc8k)
    force_download  : re-download even if files exist
    """
    from .io import read_10x

    if data_dir is None:
        data_dir = Path.home() / ".shanuz_data" / "pbmc8k"
    else:
        data_dir = Path(data_dir)

    matrix_dir = data_dir / _PBMC8K_DIR

    if force_download or not (matrix_dir / "matrix.mtx").exists():
        _download_10x(_PBMC8K_URLS, data_dir, label="PBMC8k", size_mb=38)

    mat, genes, cells = read_10x(matrix_dir, var_names="gene_symbols")
    return mat, genes, cells


def _download_10x(urls: list[str], dest_dir: Path, label: str, size_mb: int) -> None:
    """Download and extract a 10x Genomics ``*.tar.gz`` matrix bundle."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    tar_path = dest_dir / "download.tar.gz"

    print(f"Downloading {label} dataset (~{size_mb} MB)...")
    print(f"  Dest: {tar_path}")

    for url in urls:
        print(f"  Trying: {url}")
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": (
                        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                        "AppleWebKit/537.36 (KHTML, like Gecko) "
                        "Chrome/120.0.0.0 Safari/537.36"
                    ),
                    "Accept": "*/*",
                    "Referer": "https://www.10xgenomics.com/",
                },
            )
            with urllib.request.urlopen(req, timeout=120) as response:
                total = int(response.headers.get("Content-Length", 0))
                downloaded = 0
                chunk = 65536
                with open(tar_path, "wb") as fh:
                    while True:
                        block = response.read(chunk)
                        if not block:
                            break
                        fh.write(block)
                        downloaded += len(block)
                        if total > 0:
                            pct = min(downloaded / total * 100, 100)
                            print(f"\r  {pct:.1f}%", end="", flush=True)
            print()  # newline after progress
            break  # success
        except Exception as e:
            print(f"\n  Failed ({e}), trying next source...")
    else:
        raise RuntimeError(
            f"Could not download {label} data from any known source.\n"
            "Please download manually from:\n  " + "\n  ".join(urls) +
            "\nand extract to: " + str(dest_dir)
        )

    print("Extracting...")
    with tarfile.open(tar_path, "r:gz") as tf:
        tf.extractall(dest_dir)

    os.unlink(tar_path)
    print("Done.")
