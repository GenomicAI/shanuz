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
