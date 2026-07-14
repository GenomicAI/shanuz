"""Tests for load_merscope (Vizgen MERSCOPE / Seurat's LoadVizgen)."""
import numpy as np
import pandas as pd
import pytest

from shanuz import load_merscope


def _write_merscope(tmp_path, cell_id_col="cell", with_blanks=True, with_fov=True):
    """Write a minimal MERSCOPE bundle: cell_by_gene.csv + cell_metadata.csv."""
    cells = [f"cell-{i}" for i in range(6)]
    expr = {
        cell_id_col: cells,
        "Gad1": [1, 0, 2, 0, 1, 3],
        "Slc17a7": [0, 1, 1, 0, 2, 0],
        "Sox9": [2, 2, 0, 1, 0, 1],
    }
    if with_blanks:
        expr["Blank-1"] = [9, 9, 9, 9, 9, 9]
        expr["Blank-42"] = [7, 7, 7, 7, 7, 7]
    pd.DataFrame(expr).to_csv(tmp_path / "cell_by_gene.csv", index=False)

    meta = {
        cell_id_col: cells,
        "center_x": np.arange(6, dtype=float),
        "center_y": np.arange(6, dtype=float)[::-1],
        "volume": np.full(6, 100.0),
    }
    if with_fov:
        meta["fov"] = ["A", "A", "A", "B", "B", "B"]
    pd.DataFrame(meta).to_csv(tmp_path / "cell_metadata.csv", index=False)
    return cells


def test_load_merscope(tmp_path):
    cells = _write_merscope(tmp_path)

    obj = load_merscope(tmp_path)

    assert obj.active_assay == "Vizgen"
    assert obj.project_name == "MERSCOPE"
    assert set(obj.feature_names()) == {"Gad1", "Slc17a7", "Sox9"}  # blanks dropped
    assert obj.cell_names() == cells
    # One image per FOV, with every cell placed.
    assert obj.image_names() == ["A", "B"]
    coords = obj.get_tissue_coordinates()
    assert len(coords) == 6
    np.testing.assert_allclose(
        coords.sort_values("cell")["x"].to_numpy(), np.arange(6, dtype=float)
    )
    # Non-coordinate metadata carried over.
    assert "volume" in obj.meta_data.columns


def test_load_merscope_keeps_blanks_when_asked(tmp_path):
    _write_merscope(tmp_path)
    obj = load_merscope(tmp_path, keep_controls=True)
    assert {"Blank-1", "Blank-42"} <= set(obj.feature_names())


def test_load_merscope_single_image_without_fov(tmp_path):
    _write_merscope(tmp_path, with_fov=False)
    obj = load_merscope(tmp_path)
    assert len(obj.image_names()) == 1
    assert len(obj.get_tissue_coordinates()) == 6


def test_load_merscope_unnamed_cell_id_column(tmp_path):
    """Some Vizgen exports leave the leading cell-id column unnamed."""
    _write_merscope(tmp_path, cell_id_col="Unnamed: 0")
    obj = load_merscope(tmp_path)
    assert len(obj.cell_names()) == 6
    assert set(obj.feature_names()) == {"Gad1", "Slc17a7", "Sox9"}


def test_load_merscope_missing_files_raise(tmp_path):
    with pytest.raises(FileNotFoundError, match="cell_by_gene"):
        load_merscope(tmp_path)


def test_load_merscope_missing_centers_raise(tmp_path):
    _write_merscope(tmp_path)
    # Drop the center_x/center_y columns -> loader must complain clearly.
    meta = pd.read_csv(tmp_path / "cell_metadata.csv").drop(
        columns=["center_x", "center_y"])
    meta.to_csv(tmp_path / "cell_metadata.csv", index=False)
    with pytest.raises(ValueError, match="center_x"):
        load_merscope(tmp_path)
