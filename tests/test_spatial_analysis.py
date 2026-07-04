"""Tests for the spatial-analysis / composition / loader features."""
import gzip

import numpy as np
import pandas as pd
import pytest
import scipy.io
import scipy.sparse as sp

from shanuz import (
    build_niche_assay,
    composition_test,
    create_shanuz_object,
    get_tissue_coordinates,
    local_neighborhood,
    nearest_neighbor_distance,
    spatial_knn,
)
from shanuz.spatial import create_fovs, load_xenium


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def spatial_seurat():
    """40 cells over 2 FOVs; cell_type (Mast/Other) and condition (A/B)."""
    rng = np.random.default_rng(0)
    n, g = 40, 15
    counts = sp.csc_matrix(rng.poisson(1.0, size=(g, n)).astype(float))
    feats = [f"gene_{i}" for i in range(g)]
    cells = [f"cell_{i}" for i in range(n)]
    obj = create_shanuz_object(counts, assay="Xenium", feature_names=feats,
                               cell_names=cells)
    # Two 20-cell FOVs laid on a jittered grid.
    xs, ys, fov = [], [], []
    for f in range(2):
        for i in range(20):
            xs.append((i % 5) * 10 + rng.uniform(-0.5, 0.5))
            ys.append((i // 5) * 10 + rng.uniform(-0.5, 0.5))
            fov.append(f"fov{f}")
    coords = pd.DataFrame({"x": xs, "y": ys, "cell": cells})
    obj.images = create_fovs(coords, fov=np.array(fov))
    # Mast cells over-represented in condition B.
    cell_type = (["Mast"] * 4 + ["Other"] * 16) + (["Mast"] * 12 + ["Other"] * 8)
    condition = ["A"] * 20 + ["B"] * 20
    obj.meta_data["cell_type"] = cell_type
    obj.meta_data["condition"] = condition
    return obj


# ---------------------------------------------------------------------------
# Low-level KNN (FNN replacement)
# ---------------------------------------------------------------------------

def test_spatial_knn_self_excluded_matches_bruteforce():
    rng = np.random.default_rng(1)
    coords = rng.uniform(0, 100, size=(30, 2))
    dist, idx = spatial_knn(coords, k=3)
    assert dist.shape == (30, 3) and idx.shape == (30, 3)
    # brute-force check for point 0
    d0 = np.linalg.norm(coords - coords[0], axis=1)
    d0[0] = np.inf                       # exclude self
    expected = np.sort(d0)[:3]
    np.testing.assert_allclose(np.sort(dist[0]), expected, rtol=1e-6)
    assert 0 not in idx[0]               # self never returned


def test_spatial_knn_query_mode():
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [5.0, 0.0]])
    query = np.array([[0.1, 0.0]])
    dist, idx = spatial_knn(coords, k=1, query=query)
    assert idx[0, 0] == 0                # nearest to (0.1,0) is (0,0)


# ---------------------------------------------------------------------------
# Coordinate access
# ---------------------------------------------------------------------------

def test_get_tissue_coordinates(spatial_seurat):
    coords = get_tissue_coordinates(spatial_seurat)
    assert list(coords.columns) == ["x", "y", "cell", "image"]
    assert len(coords) == 40
    assert set(coords["image"]) == {"fov0", "fov1"}
    # object-level accessor agrees
    assert spatial_seurat.image_names() == ["fov0", "fov1"]
    assert len(spatial_seurat.get_tissue_coordinates("fov0")) == 20


# ---------------------------------------------------------------------------
# Nearest-neighbour distance / local neighbourhood
# ---------------------------------------------------------------------------

def test_nearest_neighbor_distance_same_type(spatial_seurat):
    df = nearest_neighbor_distance(spatial_seurat, "cell_type", "Mast")
    assert set(df.columns) == {"cell", "image", "reference", "target", "distance"}
    assert (df["distance"] >= 0).all()
    assert (df["reference"] == "Mast").all() and (df["target"] == "Mast").all()


def test_local_neighborhood_density(spatial_seurat):
    nb = local_neighborhood(spatial_seurat, "cell_type", reference="Mast", k=5)
    assert "prop_Mast" in nb.columns and "n_Other" in nb.columns
    props = nb["prop_Mast"]
    assert ((props >= 0) & (props <= 1)).all()
    # each row sums k neighbours
    ncols = [c for c in nb.columns if c.startswith("n_")]
    assert (nb[ncols].sum(axis=1) == 5).all()


# ---------------------------------------------------------------------------
# Niche assay
# ---------------------------------------------------------------------------

def test_build_niche_assay(spatial_seurat):
    build_niche_assay(spatial_seurat, "cell_type", k=5, niches=2)
    assert "niche" in spatial_seurat.assays
    assert "niches" in spatial_seurat.meta_data.columns
    niche_assay = spatial_seurat.assays["niche"]
    assert sorted(niche_assay.features()) == ["Mast", "Other"]


# ---------------------------------------------------------------------------
# Composition test
# ---------------------------------------------------------------------------

def test_composition_test_direction_and_stats(spatial_seurat):
    res = composition_test(spatial_seurat, group_by="cell_type",
                           split_by="condition", reference="A")
    assert {"group", "log2_ratio", "p", "padj", "sig", "enriched_in"} <= set(res.columns)
    assert "chisq_p" in res.attrs
    # Mast were seeded as enriched in condition B (16 in B vs 4 in A)
    mast = res.set_index("group").loc["Mast"]
    assert mast["log2_ratio"] > 0
    assert mast["enriched_in"] == "B"
    assert (res["padj"] >= res["p"] - 1e-12).all()


def test_composition_test_requires_two_levels(small_seurat):
    small_seurat.meta_data["grp"] = ["a"] * 10 + ["b"] * 10
    small_seurat.meta_data["cond3"] = (["x"] * 7 + ["y"] * 7 + ["z"] * 6)
    with pytest.raises(ValueError):
        composition_test(small_seurat, "grp", "cond3")


# ---------------------------------------------------------------------------
# add_module_score(search=)
# ---------------------------------------------------------------------------

def test_resolve_symbols_matching():
    from shanuz.module_score import _resolve_symbols
    feats = ["HLA-A", "cd8a", "MS4A1", "CD3D"]     # mixed case / punctuation
    # search resolves case + punctuation differences
    assert _resolve_symbols(["HLA.A", "CD8A"], feats, search=True) == ["HLA-A", "cd8a"]
    # without search, non-verbatim aliases resolve to nothing
    assert _resolve_symbols(["HLA.A", "CD8A"], feats, search=False) == []
    # verbatim matches always pass through
    assert _resolve_symbols(["MS4A1"], feats, search=False) == ["MS4A1"]


def test_add_module_score_search_writes_column():
    from shanuz import add_module_score
    from shanuz.preprocessing import normalize_data
    rng = np.random.default_rng(3)
    counts = sp.csc_matrix(rng.poisson(1.0, size=(6, 40)).astype(float))
    feats = ["HLA-A", "cd8a", "MS4A1", "CD3D", "GENE5", "GENE6"]
    obj = create_shanuz_object(counts, feature_names=feats,
                               cell_names=[f"c{i}" for i in range(40)])
    normalize_data(obj)
    # aliases resolve only when search=True; missing → all-zero column
    add_module_score(obj, {"hit": ["HLA.A", "CD8A"]}, nbin=3, ctrl=2, search=True, seed=1)
    add_module_score(obj, {"miss": ["HLA.A", "CD8A"]}, nbin=3, ctrl=2, search=False, seed=1)
    assert "hit" in obj.meta_data.columns
    assert np.allclose(obj.meta_data["miss"].to_numpy(), 0)


# ---------------------------------------------------------------------------
# Spatial plots
# ---------------------------------------------------------------------------

def test_image_plots_return_figures(spatial_seurat):
    import matplotlib
    matplotlib.use("Agg")
    from shanuz import image_dim_plot, image_feature_plot
    from shanuz.preprocessing import normalize_data
    normalize_data(spatial_seurat, assay="Xenium")
    fig1 = image_dim_plot(spatial_seurat, group_by="cell_type")
    fig2 = image_feature_plot(spatial_seurat, "gene_0")
    assert fig1 is not None and fig2 is not None
    assert len(fig1.axes) >= 2          # one panel per FOV (+ maybe colourbar)


# ---------------------------------------------------------------------------
# from_anndata spatial round-trip (the P0 regression)
# ---------------------------------------------------------------------------

def test_from_anndata_builds_images_not_reduction():
    anndata = pytest.importorskip("anndata")
    rng = np.random.default_rng(5)
    n, g = 30, 10
    X = sp.csr_matrix(rng.poisson(1.0, size=(n, g)).astype(float))
    obs = pd.DataFrame({"fov": ["s1"] * 15 + ["s2"] * 15},
                       index=[f"cell{i}" for i in range(n)])
    var = pd.DataFrame(index=[f"g{i}" for i in range(g)])
    xy = rng.uniform(0, 50, size=(n, 2))
    adata = anndata.AnnData(X=X, obs=obs, var=var, obsm={"spatial": xy})

    from shanuz.compat.anndata import from_anndata
    obj = from_anndata(adata, assay="Xenium")
    # spatial must become images, NOT a reduction
    assert "spatial" not in obj.reductions
    assert obj.image_names() == ["s1", "s2"]
    coords = obj.get_tissue_coordinates()
    assert len(coords) == 30
    # coordinates round-trip
    got = coords.set_index("cell").loc["cell0", ["x", "y"]].to_numpy(dtype=float)
    np.testing.assert_allclose(got, xy[0], rtol=1e-6)


# ---------------------------------------------------------------------------
# Xenium loader (synthetic on-disk bundle)
# ---------------------------------------------------------------------------

def test_load_xenium(tmp_path):
    genes = ["KIT", "TPSAB1", "CD3D"]
    cells = [f"cell-{i}" for i in range(6)]
    mat = sp.csc_matrix(np.array([[1, 0, 2, 0, 1, 3],
                                  [0, 1, 1, 0, 2, 0],
                                  [2, 2, 0, 1, 0, 1]], dtype=float))
    mtx_dir = tmp_path / "cell_feature_matrix"
    mtx_dir.mkdir()
    with gzip.open(mtx_dir / "matrix.mtx.gz", "wb") as fh:
        scipy.io.mmwrite(fh, mat)
    with gzip.open(mtx_dir / "features.tsv.gz", "wt") as fh:
        for gname in genes:
            fh.write(f"{gname}\t{gname}\tGene Expression\n")
    with gzip.open(mtx_dir / "barcodes.tsv.gz", "wt") as fh:
        for c in cells:
            fh.write(c + "\n")
    pd.DataFrame({
        "cell_id": cells,
        "x_centroid": np.arange(6, dtype=float),
        "y_centroid": np.arange(6, dtype=float)[::-1],
        "fov": ["A", "A", "A", "B", "B", "B"],
    }).to_csv(tmp_path / "cells.csv", index=False)

    obj = load_xenium(tmp_path, fov_column="fov")
    assert obj.active_assay == "Xenium"
    assert set(obj.feature_names()) == set(genes)
    assert obj.image_names() == ["A", "B"]
    assert len(obj.get_tissue_coordinates()) == 6


def test_load_xenium_drops_control_features(tmp_path):
    # Two real genes + two control-probe rows; controls dropped by default,
    # kept with keep_controls=True (mirrors Seurat's LoadXenium assay split).
    rows = [("KIT", "Gene Expression"), ("TPSAB1", "Gene Expression"),
            ("NegControlProbe_1", "Negative Control Probe"),
            ("BLANK_1", "Blank Codeword")]
    cells = [f"cell-{i}" for i in range(4)]
    mat = sp.csc_matrix(np.array([[1, 0, 2, 1],
                                  [0, 1, 1, 0],
                                  [5, 5, 5, 5],       # control noise
                                  [7, 7, 7, 7]], dtype=float))
    mtx_dir = tmp_path / "cell_feature_matrix"
    mtx_dir.mkdir()
    with gzip.open(mtx_dir / "matrix.mtx.gz", "wb") as fh:
        scipy.io.mmwrite(fh, mat)
    with gzip.open(mtx_dir / "features.tsv.gz", "wt") as fh:
        for name, ftype in rows:
            fh.write(f"{name}\t{name}\t{ftype}\n")
    with gzip.open(mtx_dir / "barcodes.tsv.gz", "wt") as fh:
        fh.write("\n".join(cells) + "\n")
    pd.DataFrame({"cell_id": cells,
                  "x_centroid": np.arange(4, dtype=float),
                  "y_centroid": np.arange(4, dtype=float)}).to_csv(
        tmp_path / "cells.csv", index=False)

    obj = load_xenium(tmp_path)
    assert set(obj.feature_names()) == {"KIT", "TPSAB1"}   # controls gone

    obj_all = load_xenium(tmp_path, keep_controls=True)
    assert set(obj_all.feature_names()) == {n for n, _ in rows}


def test_load_xenium_csv_fallback_without_parquet_engine(tmp_path, monkeypatch):
    # A cells.parquet is preferred, but if no parquet engine is available the
    # loader must transparently fall back to cells.csv[.gz].
    genes = ["KIT", "TPSAB1"]
    cells = ["c0", "c1", "c2"]
    mat = sp.csc_matrix(np.array([[1, 2, 0], [0, 1, 3]], dtype=float))
    mtx_dir = tmp_path / "cell_feature_matrix"
    mtx_dir.mkdir()
    with gzip.open(mtx_dir / "matrix.mtx.gz", "wb") as fh:
        scipy.io.mmwrite(fh, mat)
    with gzip.open(mtx_dir / "features.tsv.gz", "wt") as fh:
        for g in genes:
            fh.write(f"{g}\t{g}\tGene Expression\n")
    with gzip.open(mtx_dir / "barcodes.tsv.gz", "wt") as fh:
        fh.write("\n".join(cells) + "\n")
    coords = pd.DataFrame({"cell_id": cells,
                           "x_centroid": [0.0, 1.0, 2.0],
                           "y_centroid": [2.0, 1.0, 0.0]})
    coords.to_csv(tmp_path / "cells.csv.gz", index=False)
    (tmp_path / "cells.parquet").write_bytes(b"not-a-real-parquet")  # unreadable

    def _boom(*a, **k):
        raise ImportError("no parquet engine")
    monkeypatch.setattr(pd, "read_parquet", _boom)

    obj = load_xenium(tmp_path)          # must not raise; uses the csv
    assert len(obj.cell_names()) == 3
    assert len(obj.get_tissue_coordinates()) == 3
