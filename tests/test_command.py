from datetime import datetime
from shanuz import ShanuzCommand
from shanuz.command import log_shanuz_command


def test_create():
    cmd = ShanuzCommand(name="normalize", call_string="normalize(object, scale=1e4)")
    assert cmd.name == "normalize"
    assert isinstance(cmd.time_stamp, datetime)


def test_param_access():
    cmd = ShanuzCommand(name="f", params={"scale": 1e4, "method": "LogNormalize"})
    assert cmd.scale == 1e4
    assert cmd.method == "LogNormalize"


def test_getitem():
    cmd = ShanuzCommand(name="f", params={"x": 1})
    assert cmd["name"] == "f"
    assert cmd["x"] == 1


def test_as_list():
    cmd = ShanuzCommand(name="f", params={"a": 1}, call_string="f(a=1)")
    d = cmd.as_list()
    assert "a" in d
    assert "call_string" in d


def test_log_command(small_seurat):
    cmd = log_shanuz_command(small_seurat, "test_fn", params={"resolution": 0.5})
    assert cmd.name == "test_fn"
    assert len(small_seurat.commands) == 1


def test_repr():
    cmd = ShanuzCommand(name="f", call_string="f()", params={"x": 1})
    assert "ShanuzCommand" in repr(cmd)


# ----------------------------------------------------------------------
# The command log — R fidelity
#
# Found by the object-internals tutorial (T-obj): `log_shanuz_command` was a
# public export with zero call sites, so `obj.commands` was always empty
# where Seurat logs one entry per pipeline step.
# ----------------------------------------------------------------------


def test_command_key_matches_seurats_naming():
    cmd = log_shanuz_command(
        _Recorder(), "FindNeighbors", assay="RNA", reduction="pca")
    # Seurat indexes obj@commands by name.assay[.reduction]; scripts look
    # entries up by exactly this string.
    assert cmd.key == "FindNeighbors.RNA.pca"
    assert cmd.name == "FindNeighbors"


def test_command_key_without_a_reduction():
    cmd = log_shanuz_command(_Recorder(), "NormalizeData", assay="RNA")
    assert cmd.key == "NormalizeData.RNA"


class _Recorder:
    """Minimal stand-in for the object's command list."""

    def __init__(self):
        self.commands = []


def test_standard_pipeline_logs_the_same_steps_as_seurat(small_seurat):
    # Measured against Seurat 5.5.1 on pbmc3k, which logs exactly these five
    # keys for this sequence of calls.
    from shanuz import (find_neighbors, find_variable_features, normalize_data,
                        run_pca, scale_data)

    normalize_data(small_seurat)
    find_variable_features(small_seurat, nfeatures=20)
    scale_data(small_seurat)
    run_pca(small_seurat, n_pcs=5)
    find_neighbors(small_seurat, dims=list(range(5)), k_param=5)

    assert [c.key for c in small_seurat.commands] == [
        "NormalizeData.RNA",
        "FindVariableFeatures.RNA",
        "ScaleData.RNA",
        "RunPCA.RNA",
        "FindNeighbors.RNA.pca",
    ]
