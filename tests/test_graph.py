import numpy as np
import pytest
import scipy.sparse as sp
from shanuz import Graph, Neighbor
from shanuz.graph import as_graph


@pytest.fixture
def simple_graph(cell_names):
    mat = sp.eye(20, format="csc")
    return Graph(matrix=mat, cell_names=cell_names)


def test_cells(simple_graph, cell_names):
    assert simple_graph.cells() == cell_names


def test_shape(simple_graph):
    assert simple_graph.shape == (20, 20)


def test_assay_used_default(simple_graph):
    assert simple_graph.default_assay() is None


def test_set_assay(simple_graph):
    simple_graph.set_default_assay("RNA")
    assert simple_graph.default_assay() == "RNA"


def test_as_graph_from_matrix(cell_names):
    mat = sp.eye(20, format="csc")
    g = as_graph(mat, cell_names=cell_names)
    assert isinstance(g, Graph)
    assert g.shape == (20, 20)


def test_as_neighbor(simple_graph):
    n = simple_graph.as_neighbor()
    assert isinstance(n, Neighbor)


def test_repr(simple_graph):
    assert "Graph" in repr(simple_graph)
