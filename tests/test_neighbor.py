import numpy as np
import pytest
from shanuz import Neighbor
from shanuz.graph import Graph


@pytest.fixture
def nn(cell_names):
    rng = np.random.default_rng(1)
    idx = rng.integers(0, 20, size=(20, 5))
    dist = rng.uniform(0, 1, size=(20, 5))
    return Neighbor(nn_idx=idx, nn_dist=dist, cell_names=cell_names)


def test_cells(nn, cell_names):
    assert nn.cells() == cell_names


def test_dim(nn):
    assert nn.dim() == (20, 5)


def test_indices_shape(nn):
    assert nn.indices().shape == (20, 5)


def test_distances_shape(nn):
    assert nn.distances().shape == (20, 5)


def test_rename_cells(nn):
    new_names = [f"x_{i}" for i in range(20)]
    renamed = nn.rename_cells(new_names=new_names)
    assert renamed.cells() == new_names


def test_as_graph(nn):
    g = nn.as_graph()
    assert isinstance(g, Graph)
    assert g.shape == (20, 20)


def test_shape_mismatch():
    with pytest.raises(ValueError):
        Neighbor(
            nn_idx=np.zeros((10, 5), int),
            nn_dist=np.zeros((10, 5)),
            cell_names=[f"c_{i}" for i in range(5)],  # wrong length
        )


def test_repr(nn):
    assert "Neighbor" in repr(nn)
