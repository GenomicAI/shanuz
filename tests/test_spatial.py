import numpy as np
import pandas as pd
import pytest
from shanuz import (
    Centroids,
    FOV,
    Molecules,
    Segmentation,
    create_centroids,
    create_fov,
    create_molecules,
    create_segmentation,
)


@pytest.fixture
def centroid_coords(cell_names):
    rng = np.random.default_rng(7)
    df = pd.DataFrame({
        "x": rng.uniform(0, 100, 20),
        "y": rng.uniform(0, 100, 20),
        "cell": cell_names,
    })
    return df


@pytest.fixture
def seg_coords(cell_names):
    rows = []
    rng = np.random.default_rng(8)
    for c in cell_names[:5]:
        for _ in range(4):  # 4 polygon vertices per cell
            rows.append({"x": rng.uniform(0, 10), "y": rng.uniform(0, 10), "cell": c})
    return pd.DataFrame(rows)


@pytest.fixture
def mol_coords(cell_names):
    rng = np.random.default_rng(9)
    n = 100
    df = pd.DataFrame({
        "x": rng.uniform(0, 100, n),
        "y": rng.uniform(0, 100, n),
        "gene": [f"gene_{i % 10}" for i in range(n)],
        "cell": [cell_names[i % 20] for i in range(n)],
    })
    return df


# ---- Centroids ----

def test_centroids_cells(centroid_coords, cell_names):
    c = create_centroids(centroid_coords)
    assert c.cells() == cell_names


def test_centroids_coords(centroid_coords):
    c = create_centroids(centroid_coords)
    df = c.get_tissue_coordinates()
    assert set(df.columns) == {"x", "y"}
    assert len(df) == 20


def test_centroids_subset(centroid_coords, cell_names):
    c = create_centroids(centroid_coords)
    sub = c.subset(cell_names[:5])
    assert len(sub.cells()) == 5


def test_centroids_rename(centroid_coords):
    c = create_centroids(centroid_coords)
    new_names = [f"new_{i}" for i in range(20)]
    renamed = c.rename_cells(new_names)
    assert renamed.cells() == new_names


# ---- Segmentation ----

def test_seg_cells(seg_coords):
    s = create_segmentation(seg_coords)
    assert len(s.cells()) == 5


def test_seg_subset(seg_coords):
    s = create_segmentation(seg_coords)
    cells = s.cells()[:2]
    sub = s.subset(cells)
    assert set(sub.cells()) == set(cells)


def test_seg_simplify(seg_coords):
    s = create_segmentation(seg_coords)
    simplified = s.simplify(tol=0.1)
    assert isinstance(simplified, Segmentation)


# ---- Molecules ----

def test_mol_cells(mol_coords, cell_names):
    m = create_molecules(mol_coords)
    assert len(m.cells()) > 0


def test_mol_features(mol_coords):
    m = create_molecules(mol_coords)
    assert len(m.features()) == 10


def test_mol_subset(mol_coords, cell_names):
    m = create_molecules(mol_coords)
    sub = m.subset(cell_names[:5])
    sub_cells = set(sub.cells())
    assert sub_cells <= set(cell_names[:5])


# ---- FOV ----

def test_fov_cells(centroid_coords, cell_names):
    fov = create_fov(centroid_coords, type_="centroids")
    assert fov.cells() == cell_names


def test_fov_get_boundary(centroid_coords):
    fov = create_fov(centroid_coords, type_="centroids")
    b = fov.get_boundaries("centroids")
    assert isinstance(b, Centroids)


def test_fov_subset(centroid_coords, cell_names):
    fov = create_fov(centroid_coords, type_="centroids")
    sub = fov.subset(cell_names[:5])
    assert len(sub.cells()) == 5


def test_fov_repr(centroid_coords):
    fov = create_fov(centroid_coords, type_="centroids")
    assert "FOV" in repr(fov)
