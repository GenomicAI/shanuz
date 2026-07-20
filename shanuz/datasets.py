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


# Visium mouse brain (sagittal anterior, section 1) — 10x Genomics public dataset,
# the anterior half of Seurat's `stxBrain`. 2,695 in-tissue spots x 32,285 genes.
#
# Fetched as the raw Space Ranger output rather than through SeuratData, because
# a loader can only be tested against the files it is supposed to read: the
# curated .rda carries a built object, not `spatial/`. This is Space Ranger 1.1.0,
# so positions arrive as the headerless `tissue_positions_list.csv` — the older
# of the two layouts both `Read10X_Image` and `load_visium` have to handle.
_VISIUM_MB_BASE = (
    "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/"
    "V1_Mouse_Brain_Sagittal_Anterior/"
    "V1_Mouse_Brain_Sagittal_Anterior"
)
_VISIUM_MB_FILES = {
    "spatial": ("_spatial.tar.gz", "spatial images + scale factors (~9 MB)"),
    "filtered_feature_bc_matrix": ("_filtered_feature_bc_matrix.tar.gz",
                                   "filtered spot matrix (~55 MB)"),
}


def visium_mouse_brain(
    data_dir: Optional[str] = None,
    force_download: bool = False,
) -> Path:
    """Download (if needed) the 10x Visium mouse-brain sagittal-anterior section.

    Fetches the Space Ranger bundle (~64 MB) into ``data_dir`` (default
    ``~/.shanuz_data/visium_mouse_brain``) and returns the folder path — ready to
    pass to :func:`shanuz.load_visium`, and to R's ``Read10X_Image`` /
    ``Load10X_Spatial``, so the same slide runs in both languages.
    """
    # A local Path rather than rebinding the str|None parameter — the older
    # loaders in this module do the latter and it is most of their mypy noise.
    root = (Path(data_dir) if data_dir is not None
            else Path.home() / ".shanuz_data" / "visium_mouse_brain")
    root.mkdir(parents=True, exist_ok=True)

    # Each component is one tarball that unpacks to a directory of the same name.
    for name, (suffix, label) in _VISIUM_MB_FILES.items():
        if not force_download and (root / name).is_dir():
            continue
        tar_dest = root / f"{name}.tar.gz"
        _download_file(_VISIUM_MB_BASE + suffix, tar_dest,
                       label=f"Visium mouse brain {label}")
        with tarfile.open(tar_dest, "r:gz") as tf:
            tf.extractall(root)
        os.unlink(tar_dest)
    return root


# PBMC "Cell Hashing" dataset (Stoeckius et al. 2018, GSE108313) — the 8-hashtag
# experiment used by Seurat's hashing vignette. Cells are labelled with 8 HTOs
# (BatchA–H) and the RNA is aligned to a *combined human+mouse* reference, so
# cross-species doublets are an independent ground truth for the HTO doublet
# calls. Files are gzipped plain text on NCBI GEO (per-sample suppl paths), so R
# (HTODemux) and Python (hto_demux) read byte-identical inputs.
#
# Note: this is the raw GEO matrix (unfiltered barcodes), *not* the pre-filtered
# `pbmc_umi_mtx.rds` the vignette downloads from Dropbox — that is an R binary
# with no clean cross-language form. So counts here will not reproduce the
# vignette's headline singlet/doublet totals; the comparison target is R-vs-Python
# parity on these same files (plus the cross-species doublet ground truth).
_HASHING_RNA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2895nnn/GSM2895282/suppl/"
    "GSM2895282_Hashtag-RNA.umi.txt.gz"
)
_HASHING_HTO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2895nnn/GSM2895283/suppl/"
    "GSM2895283_Hashtag-HTO-count.csv.gz"
)
# Rows in the HTO matrix that are QC tallies, not hashtags — dropped on load.
_HASHING_HTO_SKIP = frozenset({"bad_struct", "no_match", "total_reads"})


def pbmc_hashing(
    data_dir: Optional[str] = None,
    force_download: bool = False,
):
    """Download (if needed) and load the PBMC Cell-Hashing dataset (GSE108313).

    The 8-HTO experiment from Stoeckius et al. (2018), as used by Seurat's
    hashing vignette. Returns raw RNA counts aligned to the HTO counts on their
    shared cell barcodes; the three non-hashtag QC rows (``bad_struct``,
    ``no_match``, ``total_reads``) are dropped from the HTO matrix, leaving the
    8 hashtags (``BatchA``–``BatchH``).

    The RNA reference is a combined human+mouse genome (both ``MT-`` and ``mt-``
    genes are present), which is deliberate: cross-species doublets validate the
    HTO doublet calls.

    Returns
    -------
    (rna_counts, rna_genes, hto_counts, hto_names, cell_names)
      rna_counts : (genes x cells) csc_matrix, raw counts
      hto_counts : (8 x cells) csc_matrix, same cell order as rna_counts
    """
    from .io import _make_unique

    base = (Path(data_dir) if data_dir is not None
            else Path.home() / ".shanuz_data" / "pbmc_hashing")
    base.mkdir(parents=True, exist_ok=True)

    rna_path = base / "GSM2895282_Hashtag-RNA.umi.txt.gz"
    hto_path = base / "GSM2895283_Hashtag-HTO-count.csv.gz"
    if force_download or not rna_path.exists():
        _download_file(_HASHING_RNA_URL, rna_path, label="Cell Hashing RNA (~32 MB)")
    if force_download or not hto_path.exists():
        _download_file(_HASHING_HTO_URL, hto_path, label="Cell Hashing HTO (~1 MB)")

    rna_mat, rna_genes, rna_cells = _read_table_cached(rna_path, sep="\t")
    hto_mat, hto_names, hto_cells = _read_table_cached(hto_path, sep=",")

    # Drop the QC-tally rows, keeping only the hashtags.
    keep = [i for i, n in enumerate(hto_names) if n not in _HASHING_HTO_SKIP]
    hto_mat = hto_mat[keep, :].tocsr()
    hto_names = [hto_names[i] for i in keep]

    rna_genes = _make_unique(rna_genes)
    rna_mat, hto_mat, cells = _align_on_cells(rna_mat, rna_cells, hto_mat, hto_cells)
    return rna_mat, rna_genes, hto_mat, hto_names, cells


# THP-1 ECCITE-seq dataset (Papalexi et al. 2021, GSE153056) — the pooled CRISPR
# screen behind Seurat's Mixscape vignette (SeuratData's `thp1.eccite`). 111 gRNAs
# targeting cell-surface / regulatory genes in stimulated THP-1 cells, with paired
# RNA + ADT (protein) + HTO + GDO (guide) assays. The series-level metadata table
# carries the per-cell guide assignment already, so the Python tutorial starts
# from the same annotated state as R's LoadData("thp1.eccite").
_ECCITE_RNA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4633nnn/GSM4633614/suppl/"
    "GSM4633614_ECCITE_cDNA_counts.tsv.gz"
)
_ECCITE_ADT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4633nnn/GSM4633615/suppl/"
    "GSM4633615_ECCITE_ADT_counts.tsv.gz"
)
_ECCITE_META_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE153nnn/GSE153056/suppl/"
    "GSE153056_ECCITE_metadata.tsv.gz"
)


def thp1_eccite(
    data_dir: Optional[str] = None,
    force_download: bool = False,
):
    """Download (if needed) and load the THP-1 ECCITE-seq dataset (GSE153056).

    The pooled-CRISPR screen from Papalexi et al. (2021) used by Seurat's
    Mixscape vignette. Returns RNA + ADT counts and the per-cell metadata
    (guide assignment, targeted gene, replicate, cell-cycle phase), all aligned
    to their shared barcodes. The ``gene`` / ``guide_ID`` / ``NT`` columns of the
    metadata are the perturbation labels ``run_mixscape`` needs; ``NT`` marks the
    non-targeting controls.

    Returns
    -------
    (rna_counts, rna_genes, adt_counts, adt_names, meta, cell_names)
      rna_counts : (genes x cells) csc_matrix, raw counts
      adt_counts : (proteins x cells) csc_matrix, same cell order
      meta       : pandas.DataFrame indexed by cell barcode (guide_ID, gene,
                   NT, crispr, replicate, Phase, S.Score, G2M.Score, ...)
    """
    import pandas as pd

    from .io import _make_unique

    base = (Path(data_dir) if data_dir is not None
            else Path.home() / ".shanuz_data" / "thp1_eccite")
    base.mkdir(parents=True, exist_ok=True)

    rna_path = base / "GSM4633614_ECCITE_cDNA_counts.tsv.gz"
    adt_path = base / "GSM4633615_ECCITE_ADT_counts.tsv.gz"
    meta_path = base / "GSE153056_ECCITE_metadata.tsv.gz"
    if force_download or not rna_path.exists():
        _download_file(_ECCITE_RNA_URL, rna_path, label="ECCITE RNA (~64 MB)")
    if force_download or not adt_path.exists():
        _download_file(_ECCITE_ADT_URL, adt_path, label="ECCITE ADT (~1 MB)")
    if force_download or not meta_path.exists():
        _download_file(_ECCITE_META_URL, meta_path, label="ECCITE metadata (~1 MB)")

    rna_mat, rna_genes, rna_cells = _read_table_cached(rna_path, sep="\t")
    adt_mat, adt_names, adt_cells = _read_table_cached(adt_path, sep="\t")
    meta = pd.read_csv(meta_path, sep="\t", index_col=0)

    rna_genes = _make_unique(rna_genes)
    rna_mat, adt_mat, cells = _align_on_cells(rna_mat, rna_cells, adt_mat, adt_cells)
    # Restrict to cells that also have metadata, preserving the RNA order.
    meta_set = set(meta.index)
    keep = [i for i, c in enumerate(cells) if c in meta_set]
    cells = [cells[i] for i in keep]
    rna_mat = rna_mat[:, keep].tocsc()
    adt_mat = adt_mat[:, keep].tocsc()
    meta = meta.loc[cells]
    return rna_mat, rna_genes, adt_mat, adt_names, meta, cells


def _load_seuratdata_export(name: str, data_dir: Optional[str]):
    """Read a dataset exported by ``tutorials/export_seuratdata.R``.

    Returns ``(counts, genes, cells, meta)`` where ``meta`` is a
    :class:`pandas.DataFrame` indexed by (and aligned to) the matrix's cell
    order. Raises a helpful error if the one-time R export has not been run.
    """
    import pandas as pd

    from .io import read_10x

    base = Path(data_dir) if data_dir is not None else Path.home() / ".shanuz_data" / name

    has_mtx = (base / "matrix.mtx").exists() or (base / "matrix.mtx.gz").exists()
    if not has_mtx:
        raise FileNotFoundError(
            f"{name!r} is a curated SeuratData object with no clean cross-language "
            f"raw source, so it must be exported from R first:\n"
            f"    Rscript tutorials/export_seuratdata.R {name}\n"
            f"Expected 10x-style files in {base}."
        )
    counts, genes, cells = read_10x(base, var_names="gene_symbols")
    meta = pd.read_csv(base / "metadata.csv", index_col=0)
    meta = meta.loc[cells]  # align metadata to the matrix cell order
    return counts, genes, cells, meta


def ifnb(data_dir: Optional[str] = None):
    """Load the IFNB-stimulated PBMC dataset (Kang et al. 2018), via SeuratData.

    ~14,000 human PBMCs, half stimulated with interferon-beta and half control —
    the standard benchmark for batch integration (correcting the stim/ctrl shift
    while preserving cell type). Curated as SeuratData's ``ifnb``, so it is loaded
    through the R export bridge (see :func:`_load_seuratdata_export`).

    Returns ``(counts, genes, cells, meta)``; ``meta`` carries ``stim``
    (CTRL/STIM — the batch) and ``seurat_annotations`` (the cell types).
    """
    return _load_seuratdata_export("ifnb", data_dir)


def panc8(data_dir: Optional[str] = None):
    """Load the human pancreatic-islet dataset ``panc8`` (8 techs), via SeuratData.

    ~14,900 cells profiled across five/eight technologies (CEL-seq, CEL-seq2,
    Fluidigm C1, SMART-seq2, inDrop) — a cross-technology integration and
    reference-mapping benchmark. Loaded through the R export bridge.

    Returns ``(counts, genes, cells, meta)``; ``meta`` carries ``tech`` (the
    batch / technology) and ``celltype`` (the reference annotation).
    """
    return _load_seuratdata_export("panc8", data_dir)


def _read_dense_table_sparse(
    path: Path,
    sep: str,
    chunksize: int = 4000,
) -> tuple[sp.csc_matrix, list[str], list[str]]:
    """Read a gzipped dense ``features x cells`` text table into a sparse matrix.

    Row 0 is the header (its first field is the — possibly empty — index name,
    the rest are cell barcodes); every later row is a feature label followed by
    its counts. The table is read in row-chunks and sparsified as it goes, so a
    file whose dense form would not fit in memory still loads at bounded cost.
    Quoted fields (R's ``write.table`` default) are handled.

    Returns ``(matrix, feature_names, cell_names)`` with ``matrix`` a
    ``(features x cells)`` :class:`scipy.sparse.csc_matrix`.
    """
    import pandas as pd

    blocks: list[sp.csr_matrix] = []
    features: list[str] = []
    cells: Optional[list[str]] = None
    # na_filter=False skips pandas' NA scan — a ~2.6x speedup on these very wide
    # (tens of thousands of columns) count tables, where NA detection dominates.
    for chunk in pd.read_csv(path, sep=sep, index_col=0, chunksize=chunksize,
                             na_filter=False):
        if cells is None:
            cells = [str(c) for c in chunk.columns]
        features.extend(str(f) for f in chunk.index)
        blocks.append(sp.csr_matrix(chunk.to_numpy(dtype=np.float32)))
    if cells is None:  # empty file → no header row
        raise ValueError(f"{path} has no data rows")
    mat = sp.vstack(blocks, format="csc") if blocks else sp.csc_matrix((0, len(cells)))
    return mat, features, cells


def _read_table_cached(
    path: Path,
    sep: str,
    chunksize: int = 4000,
) -> tuple[sp.csc_matrix, list[str], list[str]]:
    """:func:`_read_dense_table_sparse`, memoised to a ``.parsed.npz`` sidecar.

    These GEO count tables are dense text tens of thousands of columns wide;
    parsing one takes a couple of minutes. The sparse result is cached next to
    the source file so only the first load pays that cost — every later load is
    a fast :func:`scipy.sparse.load_npz`. The cache is ignored if the source is
    newer, so re-downloading transparently reparses.
    """
    mat_cache = Path(str(path) + ".parsed.npz")
    names_cache = Path(str(path) + ".parsed.names.npz")
    if (mat_cache.exists() and names_cache.exists()
            and mat_cache.stat().st_mtime >= path.stat().st_mtime):
        mat = sp.load_npz(mat_cache).tocsc()
        nz = np.load(names_cache)
        return mat, list(map(str, nz["features"])), list(map(str, nz["cells"]))

    mat, features, cells = _read_dense_table_sparse(path, sep, chunksize)
    sp.save_npz(mat_cache, mat)
    np.savez(names_cache, features=np.array(features), cells=np.array(cells))
    return mat, features, cells


def _align_on_cells(
    mat_a: sp.spmatrix,
    cells_a: list[str],
    mat_b: sp.spmatrix,
    cells_b: list[str],
) -> tuple[sp.csc_matrix, sp.csc_matrix, list[str]]:
    """Subset two feature×cell matrices to shared barcodes, ordered by ``mat_a``."""
    b_pos = {c: j for j, c in enumerate(cells_b)}
    common = [c for c in cells_a if c in b_pos]
    a_pos = {c: j for j, c in enumerate(cells_a)}
    cols_a = [a_pos[c] for c in common]
    cols_b = [b_pos[c] for c in common]
    return mat_a.tocsc()[:, cols_a], mat_b.tocsc()[:, cols_b], common


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
