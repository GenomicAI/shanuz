"""Unit tests for the tutorial dataset loaders' parsing layer.

These exercise the pure parse/cache/align helpers and the loader plumbing on
tiny *synthetic* gzipped tables — no network, so they run in CI. The real GEO
downloads are covered by the opt-in tutorial smoke tests, not here.

The formats mirror the real files (verified 2026-07 against GEO):
  * GSE108313 HTO — a comma-separated ``features x cells`` table whose non-hashtag
    QC rows (``bad_struct`` / ``no_match`` / ``total_reads``) must be dropped;
  * GSE153056 ECCITE — tab-separated tables with R ``write.table`` quoting;
  * ifnb / panc8 — a 10x-style folder written by ``export_seuratdata.R``.
"""
from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

from shanuz import datasets
from shanuz.datasets import (
    _align_on_cells,
    _read_dense_table_sparse,
    _read_table_cached,
    cbmc_citeseq,
    ifnb,
    pbmc3k,
    pbmc8k,
    pbmc_hashing,
    xenium_mouse_brain,
)


def _write_gzip_table(path: Path, header: list[str], rows: list[list], sep: str,
                      quote: bool = False) -> None:
    """Write a ``features x cells`` table (row 0 = header) as gzipped text."""
    q = csv.QUOTE_NONNUMERIC if quote else csv.QUOTE_MINIMAL
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter=sep, quoting=q)
        w.writerow(header)
        for r in rows:
            w.writerow(r)


# --------------------------------------------------------------------------- #
# _read_dense_table_sparse
# --------------------------------------------------------------------------- #

def test_read_dense_table_sparse_round_trips(tmp_path):
    p = tmp_path / "t.tsv.gz"
    _write_gzip_table(
        p,
        header=["", "cellA", "cellB", "cellC"],
        rows=[["g1", 0, 3, 0], ["g2", 1, 0, 0], ["g3", 0, 0, 5]],
        sep="\t",
    )
    mat, feats, cells = _read_dense_table_sparse(p, sep="\t")
    assert isinstance(mat, sp.csc_matrix)
    assert feats == ["g1", "g2", "g3"]
    assert cells == ["cellA", "cellB", "cellC"]
    np.testing.assert_array_equal(
        mat.toarray(), np.array([[0, 3, 0], [1, 0, 0], [0, 0, 5]], dtype=np.float32)
    )


def test_read_dense_table_sparse_handles_quoting_and_chunk_boundaries(tmp_path):
    # R's write.table quotes strings; the chunked reader must stitch chunks that
    # split the rows, so force chunksize below the row count.
    p = tmp_path / "q.tsv.gz"
    rows = [[f"gene{i}", i % 3, 0, (i + 1) % 2] for i in range(10)]
    _write_gzip_table(p, header=["", "c1", "c2", "c3"], rows=rows, sep="\t", quote=True)
    mat, feats, cells = _read_dense_table_sparse(p, sep="\t", chunksize=3)
    assert feats == [f"gene{i}" for i in range(10)]
    assert cells == ["c1", "c2", "c3"]
    assert mat.shape == (10, 3)
    assert mat[0, 1] == 0 and mat[1, 0] == 1 and mat[2, 0] == 2


def test_read_dense_table_sparse_comma_separator(tmp_path):
    p = tmp_path / "c.csv.gz"
    _write_gzip_table(p, header=["", "x", "y"], rows=[["r1", 2, 0], ["r2", 0, 4]], sep=",")
    mat, feats, cells = _read_dense_table_sparse(p, sep=",")
    assert feats == ["r1", "r2"] and cells == ["x", "y"]
    assert mat.toarray().tolist() == [[2, 0], [0, 4]]


# --------------------------------------------------------------------------- #
# _read_table_cached
# --------------------------------------------------------------------------- #

def test_read_table_cached_writes_and_reuses_sidecar(tmp_path):
    p = tmp_path / "t.tsv.gz"
    _write_gzip_table(p, header=["", "a", "b"], rows=[["g1", 1, 2], ["g2", 0, 3]], sep="\t")
    m1, f1, c1 = _read_table_cached(p, sep="\t")
    assert (tmp_path / "t.tsv.gz.parsed.npz").exists()
    assert (tmp_path / "t.tsv.gz.parsed.names.npz").exists()
    # second call must come from the cache and be identical
    m2, f2, c2 = _read_table_cached(p, sep="\t")
    assert f2 == f1 == ["g1", "g2"] and c2 == c1 == ["a", "b"]
    assert (m1 != m2).nnz == 0


def test_read_table_cached_reparses_when_source_is_newer(tmp_path):
    import os
    import time

    p = tmp_path / "t.tsv.gz"
    _write_gzip_table(p, header=["", "a"], rows=[["g1", 1]], sep="\t")
    m1, _, _ = _read_table_cached(p, sep="\t")
    assert m1.toarray().tolist() == [[1]]
    # rewrite the source with different data and a strictly newer mtime
    time.sleep(0.01)
    _write_gzip_table(p, header=["", "a"], rows=[["g1", 9]], sep="\t")
    os.utime(p, (time.time() + 5, time.time() + 5))
    m2, _, _ = _read_table_cached(p, sep="\t")
    assert m2.toarray().tolist() == [[9]], "stale cache was served after source changed"


# --------------------------------------------------------------------------- #
# _align_on_cells
# --------------------------------------------------------------------------- #

def test_align_on_cells_intersects_and_orders_by_first():
    a = sp.csc_matrix(np.array([[1, 2, 3]], dtype=np.float32))  # cells p,q,r
    b = sp.csc_matrix(np.array([[7, 8, 9]], dtype=np.float32))  # cells r,q,z
    aa, bb, common = _align_on_cells(a, ["p", "q", "r"], b, ["r", "q", "z"])
    assert common == ["q", "r"]                     # ordered by a, intersection only
    assert aa.toarray().tolist() == [[2, 3]]
    assert bb.toarray().tolist() == [[8, 7]]         # b reordered to match


# --------------------------------------------------------------------------- #
# pbmc_hashing plumbing (QC-row drop + alignment), no network
# --------------------------------------------------------------------------- #

def test_pbmc_hashing_drops_qc_rows_and_aligns(tmp_path):
    # RNA over cells c1..c3; HTO over c2..c4 (+ QC rows). Shared = c2,c3.
    _write_gzip_table(
        tmp_path / "GSM2895282_Hashtag-RNA.umi.txt.gz",
        header=["GENE", "c1", "c2", "c3"],
        rows=[["Xkr4", 0, 5, 1], ["MALAT1", 2, 0, 3]],
        sep="\t",
    )
    _write_gzip_table(
        tmp_path / "GSM2895283_Hashtag-HTO-count.csv.gz",
        header=["", "c2", "c3", "c4"],
        rows=[
            ["BatchA-AAA", 10, 0, 1],
            ["BatchB-CCC", 0, 20, 2],
            ["bad_struct", 3, 3, 3],
            ["no_match", 4, 4, 4],
            ["total_reads", 99, 99, 99],
        ],
        sep=",",
    )
    rna, genes, hto, hto_names, cells = pbmc_hashing(data_dir=str(tmp_path))
    assert cells == ["c2", "c3"]                       # RNA order, shared with HTO
    assert hto_names == ["BatchA-AAA", "BatchB-CCC"]   # 3 QC rows dropped
    assert hto.shape == (2, 2)
    assert rna.shape == (2, 2)
    assert hto.toarray().tolist() == [[10, 0], [0, 20]]
    assert rna.toarray().tolist() == [[5, 1], [0, 3]]  # cols c2,c3 of RNA


# --------------------------------------------------------------------------- #
# SeuratData export reader (ifnb / panc8)
# --------------------------------------------------------------------------- #

def test_ifnb_errors_when_bridge_not_run(tmp_path):
    with pytest.raises(FileNotFoundError, match="export_seuratdata.R ifnb"):
        ifnb(data_dir=str(tmp_path))


# --------------------------------------------------------------------------- #
# data_dir is honoured (no network, and the default cache is made unreachable)
# --------------------------------------------------------------------------- #

@pytest.fixture
def offline(monkeypatch, tmp_path):
    """Point ``Path.home()`` at an empty dir and make any download fail loudly.

    Both halves are load-bearing. Redirecting home means a loader that ignored
    ``data_dir`` cannot quietly succeed off a developer's populated
    ``~/.shanuz_data``; failing the download means it announces itself here
    instead of reaching the network in CI.
    """
    fake_home = tmp_path / "home"
    fake_home.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: fake_home))

    def _no_download(*args, **kwargs):
        raise AssertionError("loader tried to download — it did not use data_dir")

    monkeypatch.setattr(datasets, "_download_10x", _no_download)
    monkeypatch.setattr(datasets, "_download_file", _no_download)
    return fake_home


def _write_10x_v2(directory: Path, mat: sp.csc_matrix, genes: list[str],
                  cells: list[str]) -> None:
    """Write a v2 (``genes.tsv``) 10x trio, the layout both PBMC bundles use."""
    import scipy.io

    directory.mkdir(parents=True, exist_ok=True)
    with open(directory / "matrix.mtx", "wb") as fh:
        scipy.io.mmwrite(fh, mat)
    (directory / "genes.tsv").write_text(
        "".join(f"ENSG{i:05d}\t{g}\n" for i, g in enumerate(genes))
    )
    (directory / "barcodes.tsv").write_text("".join(f"{c}\n" for c in cells))


@pytest.mark.parametrize("loader, subdir", [
    (pbmc3k, "filtered_gene_bc_matrices/hg19"),
    (pbmc8k, "filtered_gene_bc_matrices/GRCh38"),
])
def test_pbmc_loaders_read_the_given_data_dir(loader, subdir, tmp_path, offline):
    # The reference lives at data_dir/<bundle subdir>, so this pins the join as
    # well as the directory: a loader reading data_dir itself finds nothing.
    mat = sp.csc_matrix(np.array([[1, 0], [0, 2], [3, 4]], dtype=np.float32))
    _write_10x_v2(tmp_path / subdir, mat, ["GeneA", "GeneB", "GeneC"], ["c1", "c2"])

    counts, genes, cells = loader(data_dir=str(tmp_path))
    assert genes == ["GeneA", "GeneB", "GeneC"]
    assert cells == ["c1", "c2"]
    assert counts.toarray().tolist() == [[1, 0], [0, 2], [3, 4]]
    assert not list(offline.iterdir()), "loader wrote to the default cache"


def test_xenium_returns_the_given_data_dir(tmp_path, offline):
    # Both files present, so nothing is fetched and the folder is returned as is.
    (tmp_path / "cell_feature_matrix").mkdir()
    (tmp_path / "cell_feature_matrix" / "matrix.mtx.gz").write_bytes(b"")
    (tmp_path / "cells.csv.gz").write_bytes(b"")

    got = xenium_mouse_brain(data_dir=str(tmp_path))
    assert got == tmp_path
    assert isinstance(got, Path), "load_xenium is handed this path directly"
    assert not list(offline.iterdir()), "loader wrote to the default cache"


def test_cbmc_keeps_human_genes_and_aligns_to_shared_barcodes(tmp_path, offline):
    # RNA over c1..c3 with one mouse spike-in; ADT over c3,c2,c9. Shared = c2,c3,
    # which must come back in *RNA* order and drop MOUSE_Actb with its prefix.
    _write_gzip_table(
        tmp_path / "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",
        header=["", "c1", "c2", "c3"],
        rows=[["HUMAN_CD3E", 1, 5, 0], ["MOUSE_Actb", 9, 9, 9], ["HUMAN_MS4A1", 2, 0, 7]],
        sep=",",
    )
    _write_gzip_table(
        tmp_path / "GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz",
        header=["", "c3", "c2", "c9"],
        rows=[["CD3", 30, 20, 90], ["CD19", 33, 22, 99]],
        sep=",",
    )
    rna, genes, adt, prots, cells = cbmc_citeseq(data_dir=str(tmp_path))
    assert genes == ["CD3E", "MS4A1"], "mouse row kept, or prefix left on"
    assert cells == ["c2", "c3"]
    assert prots == ["CD3", "CD19"]
    assert rna.toarray().tolist() == [[5, 0], [0, 7]]
    assert adt.toarray().tolist() == [[20, 30], [22, 33]]  # reordered to the RNA
    assert not list(offline.iterdir()), "loader wrote to the default cache"


def test_cbmc_stitches_chunks_and_skips_an_all_mouse_one(tmp_path, offline):
    # The loader reads the RNA table in 4,000-row chunks and the real file is
    # ~36,000 rows, so every genuine load crosses several boundaries — the tests
    # above all fit in one chunk and never exercise the accumulation at all.
    # Three chunks, the middle one entirely mouse, so *two* human blocks have to
    # be stitched across a skipped one. First column carries the source row
    # number, which is what ties each matrix row back to its gene.
    rows = [[f"HUMAN_G{i}", i + 1, 1] for i in range(4000)]          # chunk 1
    rows += [[f"MOUSE_m{i}", 7, 7] for i in range(4000)]             # chunk 2, dropped
    rows += [[f"HUMAN_G{4000 + i}", 8001 + i, 1] for i in range(500)]  # chunk 3
    _write_gzip_table(tmp_path / "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",
                      header=["", "c1", "c2"], rows=rows, sep=",")
    _write_gzip_table(tmp_path / "GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz",
                      header=["", "c1", "c2"], rows=[["CD3", 4, 5]], sep=",")

    rna, genes, adt, prots, cells = cbmc_citeseq(data_dir=str(tmp_path))
    assert rna.shape == (4500, 2), "chunks not concatenated, or mouse chunk kept"
    assert genes[:2] == ["G0", "G1"] and genes[-1] == "G4499"  # source order held
    assert cells == ["c1", "c2"]
    # Rows must still line up with `genes` after the stitch: G0 came from source
    # row 1 and G4499 from row 8500, so a reordered block shows up here.
    assert rna[0, 0] == 1 and rna[4499, 0] == 8500
    # Column 2 is all ones, so a dropped or duplicated block shows up here.
    assert rna[:, 1].toarray().sum() == 4500


@pytest.mark.parametrize("rna_rows", [
    [],                              # header only
    [["MOUSE_Actb", 1]],             # rows, none matching the prefix
])
def test_cbmc_names_the_species_prefix_when_nothing_survives(rna_rows, tmp_path, offline):
    # Both cases used to reach sp.vstack([]) and report "blocks must be 2-D",
    # which says nothing about the prefix that actually filtered everything out.
    _write_gzip_table(tmp_path / "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",
                      header=["", "c1"], rows=rna_rows, sep=",")
    _write_gzip_table(tmp_path / "GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz",
                      header=["", "c1"], rows=[["CD3", 1]], sep=",")
    with pytest.raises(ValueError, match="HUMAN_"):
        cbmc_citeseq(data_dir=str(tmp_path))


def test_ifnb_reads_bridge_output(tmp_path):
    # Minimal 10x plain-v3 trio + metadata, as export_seuratdata.R writes.
    import scipy.io

    mat = sp.csc_matrix(np.array([[1, 0], [0, 2], [3, 4]], dtype=np.float32))
    with open(tmp_path / "matrix.mtx", "wb") as fh:
        scipy.io.mmwrite(fh, mat)
    (tmp_path / "features.tsv").write_text(
        "GeneA\tGeneA\tGene Expression\nGeneB\tGeneB\tGene Expression\n"
        "GeneC\tGeneC\tGene Expression\n"
    )
    (tmp_path / "barcodes.tsv").write_text("cell1\ncell2\n")
    (tmp_path / "metadata.csv").write_text(
        ",stim,seurat_annotations\ncell1,CTRL,B\ncell2,STIM,T\n"
    )
    counts, genes, cells, meta = ifnb(data_dir=str(tmp_path))
    assert genes == ["GeneA", "GeneB", "GeneC"]
    assert cells == ["cell1", "cell2"]
    assert counts.shape == (3, 2)
    # metadata is aligned to the matrix cell order
    assert list(meta.index) == ["cell1", "cell2"]
    assert list(meta["stim"]) == ["CTRL", "STIM"]
