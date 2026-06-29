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
        _download_pbmc3k(data_dir)

    mat, genes, cells = read_10x(matrix_dir, var_names="gene_symbols")
    return mat, genes, cells


def _download_pbmc3k(dest_dir: Path) -> None:
    """Download and extract the PBMC3k tarball into dest_dir."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    tar_path = dest_dir / "pbmc3k.tar.gz"

    print(f"Downloading PBMC3k dataset (~24 MB)...")
    print(f"  Dest: {tar_path}")

    for url in _PBMC3K_URLS:
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
            "Could not download PBMC3k data from any known source.\n"
            "Please download manually from:\n"
            "  https://cf.10xgenomics.com/samples/cell/pbmc3k/"
            "pbmc3k_filtered_gene_bc_matrices.tar.gz\n"
            "and extract to: " + str(dest_dir)
        )

    print("Extracting...")
    with tarfile.open(tar_path, "r:gz") as tf:
        tf.extractall(dest_dir)

    os.unlink(tar_path)
    print("Done.")
