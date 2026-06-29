import numpy as np
import pytest
from shanuz import LogMap


def test_setget():
    lm = LogMap()
    lm["a"] = [True, False, True]
    assert np.array_equal(lm["a"], [True, False, True])


def test_contains():
    lm = LogMap({"x": np.array([True])})
    assert "x" in lm
    assert "y" not in lm


def test_delete():
    lm = LogMap({"a": np.array([True])})
    del lm["a"]
    assert "a" not in lm


def test_len():
    lm = LogMap({"a": np.ones(3, bool), "b": np.zeros(3, bool)})
    assert len(lm) == 2


def test_iter():
    lm = LogMap({"a": np.ones(2, bool), "b": np.zeros(2, bool)})
    assert set(lm) == {"a", "b"}


def test_bitwise_and():
    lm1 = LogMap({"k": np.array([True, False, True])})
    lm2 = LogMap({"k": np.array([True, True, False])})
    result = lm1 & lm2
    assert np.array_equal(result["k"], [True, False, False])


def test_bitwise_or():
    lm1 = LogMap({"k": np.array([True, False, False])})
    lm2 = LogMap({"k": np.array([False, True, False])})
    result = lm1 | lm2
    assert np.array_equal(result["k"], [True, True, False])


def test_invert():
    lm = LogMap({"k": np.array([True, False])})
    result = ~lm
    assert np.array_equal(result["k"], [False, True])


def test_copy():
    lm = LogMap({"k": np.array([True])})
    lm2 = lm.copy()
    lm2["k"][0] = False
    assert lm["k"][0] is np.bool_(True) or lm["k"][0] == True


def test_repr():
    lm = LogMap({"a": np.ones(1, bool)})
    assert "LogMap" in repr(lm)
