"""Built-in dataset loaders.

Provides automatic download of commonly used benchmark datasets.
"""
from __future__ import annotations

import os
import tarfile
import urllib.request
from pathlib import Path
from typing import Optional

import numpy as np
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


# CBMC CITE-seq dataset (GSE100866) — cord-blood mononuclear cells with paired
# RNA + 13-antibody surface-protein (ADT) measurements. Used by the multimodal
# tutorial. Files are gzipped CSVs (features x cells) hosted on NCBI GEO.
_CBMC_BASE = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/"
)
_CBMC_RNA = "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz"
_CBMC_ADT = "GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz"
_CBMC_SPECIES_PREFIX = "HUMAN_"


def cbmc_citeseq(
    data_dir: Optional[str] = None,
    force_download: bool = False,
    species_prefix: str = _CBMC_SPECIES_PREFIX,
):
    """Download (if needed) and load the CBMC CITE-seq dataset (GSE100866).

    ~8,600 cord-blood mononuclear cells profiled for both RNA and 13 surface
    proteins (ADT). The RNA matrix mixes human and mouse spike-in genes; this
    loader keeps the human genes and strips the ``HUMAN_`` prefix (mirroring
    Seurat's CollapseSpeciesExpressionMatrix). RNA and ADT are aligned to their
    shared cell barcodes.

    Returns
    -------
    (rna_counts, rna_genes, adt_counts, adt_proteins, cell_names)
      rna_counts : (genes x cells) csc_matrix
      adt_counts : (proteins x cells) csc_matrix, same cell order as rna_counts
    """
    import pandas as pd

    from .io import _make_unique

    if data_dir is None:
        data_dir = Path.home() / ".shanuz_data" / "cbmc"
    else:
        data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)

    rna_path = data_dir / _CBMC_RNA
    adt_path = data_dir / _CBMC_ADT
    for fname, path in ((_CBMC_RNA, rna_path), (_CBMC_ADT, adt_path)):
        if force_download or not path.exists():
            _download_file(_CBMC_BASE + fname, path, label=fname)

    # ADT is tiny — read directly (proteins x cells).
    adt_df = pd.read_csv(adt_path, index_col=0)
    adt_cells = list(adt_df.columns)

    # RNA is larger — read in row-chunks, keeping only human genes.
    rna_blocks, rna_genes, rna_cells = [], [], None
    for chunk in pd.read_csv(rna_path, index_col=0, chunksize=4000):
        if rna_cells is None:
            rna_cells = list(chunk.columns)
        mask = np.asarray(chunk.index.str.startswith(species_prefix))
        sub = chunk[mask]
        if len(sub):
            rna_genes.extend(g[len(species_prefix):] for g in sub.index)
            rna_blocks.append(sp.csr_matrix(sub.values.astype(np.float32)))
    rna_mat = sp.vstack(rna_blocks, format="csc")
    rna_genes = _make_unique(rna_genes)

    # Align to shared barcodes, ordered by the RNA matrix.
    adt_set = set(adt_cells)
    rna_pos = {c: i for i, c in enumerate(rna_cells)}
    common = [c for c in rna_cells if c in adt_set]
    rna_cols = [rna_pos[c] for c in common]
    rna_mat = rna_mat[:, rna_cols].tocsc()
    adt_mat = sp.csc_matrix(adt_df[common].values.astype(np.float32))

    return rna_mat, rna_genes, adt_mat, list(adt_df.index), common


# Xenium mouse brain (coronal, CTX+HP subset) — 10x Genomics public dataset,
# the same section used in Seurat's spatial vignette. 36,602 cells x 248 genes.
# Only the lightweight analysis components are fetched (~20 MB), not the multi-GB
# morphology-image bundle — enough for LoadXenium/load_xenium.
_XENIUM_MB_BASE = (
    "https://cf.10xgenomics.com/samples/xenium/1.0.2/"
    "Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/"
    "Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP"
)
_XENIUM_MB_FILES = {
    "cell_feature_matrix.tar.gz": "_cell_feature_matrix.tar.gz",
    "cells.csv.gz": "_cells.csv.gz",
}


def xenium_mouse_brain(
    data_dir: Optional[str] = None,
    force_download: bool = False,
) -> Path:
    """Download (if needed) the 10x Xenium mouse-brain coronal subset.

    Fetches only the analysis components (``cell_feature_matrix/`` + ``cells``),
    ~20 MB, into ``data_dir`` (default ``~/.shanuz_data/xenium_mouse_brain``) and
    returns the folder path — ready to pass to :func:`shanuz.load_xenium`.

    This is the public section featured in Seurat's Xenium spatial vignette, so
    the same analysis runs in R (``LoadXenium``) and Python (``load_xenium``).
    """
    if data_dir is None:
        data_dir = Path.home() / ".shanuz_data" / "xenium_mouse_brain"
    else:
        data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)

    mtx_dir = data_dir / "cell_feature_matrix"
    need = force_download or not (mtx_dir / "matrix.mtx.gz").exists()
    if need:
        tar_dest = data_dir / "cell_feature_matrix.tar.gz"
        _download_file(_XENIUM_MB_BASE + _XENIUM_MB_FILES["cell_feature_matrix.tar.gz"],
                       tar_dest, label="Xenium mouse brain matrix (~11 MB)")
        with tarfile.open(tar_dest, "r:gz") as tf:
            tf.extractall(data_dir)
        os.unlink(tar_dest)

    cells_dest = data_dir / "cells.csv.gz"
    if force_download or not cells_dest.exists():
        _download_file(_XENIUM_MB_BASE + _XENIUM_MB_FILES["cells.csv.gz"],
                       cells_dest, label="Xenium mouse brain cells (~2 MB)")
    return data_dir


def _download_file(url: str, dest: Path, label: str) -> None:
    """Download a single file (with a simple progress print)."""
    print(f"Downloading {label} ...")
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Mozilla/5.0", "Accept": "*/*"},
    )
    with urllib.request.urlopen(req, timeout=300) as response:
        total = int(response.headers.get("Content-Length", 0))
        downloaded = 0
        with open(dest, "wb") as fh:
            while True:
                block = response.read(65536)
                if not block:
                    break
                fh.write(block)
                downloaded += len(block)
                if total:
                    print(f"\r  {min(downloaded / total * 100, 100):.1f}%",
                          end="", flush=True)
    print()


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
