"""Tests for MULTI-seq demultiplexing (v0.9.0 Specialized Assays).

  * multiseq_demux  (multiseq.py)

Like the HTODemux tests, a synthetic barcode matrix with *known* ground truth —
clean singlets per barcode, doublets spanning barcode pairs, empty-ish negatives —
is built and ``multiseq_demux`` is asked to recover it. The Seurat metadata
columns (``MULTI_ID`` / ``MULTI_classification``), the identity, and the
singlet / doublet / negative accuracy are checked against ground truth, for both
the single-pass and ``autothresh`` paths. Network-free and deterministic.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz.multiseq import multiseq_demux  # noqa: E402
from shanuz.preprocessing import normalize_data  # noqa: E402
from shanuz.shanuz import create_shanuz_object  # noqa: E402

TAGS = ["HTO-A", "HTO-B", "HTO-C", "HTO-D"]
DOUBLET_PAIRS = [(0, 1), (1, 2), (2, 3), (0, 2), (1, 3)]


def _hashing_counts(seed=0):
    """A barcode count matrix (4 × N) plus per-cell ground-truth labels."""
    rng = np.random.default_rng(seed)
    H = len(TAGS)
    cols, truth = [], []

    bg = 1.0
    for h in range(H):                       # 30 singlets per barcode
        for _ in range(30):
            col = rng.poisson(bg, size=H).astype(float)
            col[h] = rng.poisson(100.0)
            cols.append(col)
            truth.append(TAGS[h])

    for _ in range(4):                       # 20 doublets across barcode pairs
        for a, b in DOUBLET_PAIRS:
            col = rng.poisson(bg, size=H).astype(float)
            col[a] = rng.poisson(70.0)
            col[b] = rng.poisson(70.0)
            cols.append(col)
            truth.append("Doublet")

    for _ in range(15):                      # 15 negatives (only background)
        cols.append(rng.poisson(bg, size=H).astype(float))
        truth.append("Negative")

    counts = np.asarray(cols).T              # H × N
    return counts, np.asarray(truth)


def _hashing_object(seed=0):
    counts, truth = _hashing_counts(seed)
    cells = [f"c{i}" for i in range(counts.shape[1])]
    obj = create_shanuz_object(
        counts=sp.csc_matrix(counts), assay="HTO",
        feature_names=TAGS, cell_names=cells,
        meta_data=pd.DataFrame(index=cells),
    )
    return obj, truth


# ----------------------------------------------------------------------
# metadata columns & identity
# ----------------------------------------------------------------------


def test_writes_multi_metadata_columns():
    obj, _ = _hashing_object()
    multiseq_demux(obj)
    for col in ("MULTI_ID", "MULTI_classification"):
        assert col in obj.meta_data.columns
    assert list(obj.meta_data["MULTI_ID"]) == list(obj.meta_data["MULTI_classification"])


def test_multi_id_domain():
    obj, _ = _hashing_object()
    multiseq_demux(obj)
    vals = set(obj.meta_data["MULTI_ID"])
    assert vals <= set(TAGS) | {"Doublet", "Negative"}


def test_idents_set_to_multi_id():
    obj, _ = _hashing_object()
    multiseq_demux(obj)
    assert len(obj.idents) == len(obj.cell_names())
    assert list(obj.idents) == list(obj.meta_data["MULTI_ID"])


# ----------------------------------------------------------------------
# recovery against ground truth
# ----------------------------------------------------------------------


def test_singlets_recovered():
    obj, truth = _hashing_object()
    multiseq_demux(obj)
    call = obj.meta_data["MULTI_ID"].to_numpy()
    singlet = np.isin(truth, TAGS)
    # a called singlet is neither Doublet nor Negative, and names its true barcode
    assert (call[singlet] == truth[singlet]).mean() >= 0.9


def test_doublets_recovered():
    obj, truth = _hashing_object()
    multiseq_demux(obj)
    call = obj.meta_data["MULTI_ID"].to_numpy()
    doublet = truth == "Doublet"
    assert (call[doublet] == "Doublet").mean() >= 0.7


def test_negatives_recovered():
    obj, truth = _hashing_object()
    multiseq_demux(obj)
    call = obj.meta_data["MULTI_ID"].to_numpy()
    negative = truth == "Negative"
    assert (call[negative] == "Negative").mean() >= 0.7


# ----------------------------------------------------------------------
# options, autothresh, bookkeeping & guards
# ----------------------------------------------------------------------


def test_thresholds_stored_in_misc():
    obj, _ = _hashing_object()
    multiseq_demux(obj)
    info = obj.misc["multiseq_demux"]["HTO"]
    assert set(info["thresholds"]) == set(TAGS)
    assert info["quantile"] == 0.7
    assert info["autothresh"] is False
    # well-separated barcodes each learn a finite threshold
    assert all(np.isfinite(t) for t in info["thresholds"].values())


def test_quantile_shifts_threshold():
    obj_lo, _ = _hashing_object()
    obj_hi, _ = _hashing_object()
    multiseq_demux(obj_lo, quantile=0.3)
    multiseq_demux(obj_hi, quantile=0.9)
    lo = obj_lo.misc["multiseq_demux"]["HTO"]["thresholds"]
    hi = obj_hi.misc["multiseq_demux"]["HTO"]["thresholds"]
    # a larger quantile pushes each cutoff toward the positive peak
    for tag in TAGS:
        assert hi[tag] >= lo[tag] - 1e-9


def test_autothresh_recovers_singlets():
    obj, truth = _hashing_object()
    multiseq_demux(obj, autothresh=True)
    info = obj.misc["multiseq_demux"]["HTO"]
    assert info["autothresh"] is True
    call = obj.meta_data["MULTI_ID"].to_numpy()
    singlet = np.isin(truth, TAGS)
    assert (call[singlet] == truth[singlet]).mean() >= 0.9


def test_normalize_false_uses_data_layer():
    obj, truth = _hashing_object()
    normalize_data(obj, normalization_method="CLR", margin=2, assay="HTO")
    multiseq_demux(obj, normalize=False)
    call = obj.meta_data["MULTI_ID"].to_numpy()
    singlet = np.isin(truth, TAGS)
    assert (call[singlet] == truth[singlet]).mean() >= 0.9


def test_deterministic():
    obj_a, _ = _hashing_object()
    obj_b, _ = _hashing_object()
    multiseq_demux(obj_a)
    multiseq_demux(obj_b)
    assert list(obj_a.meta_data["MULTI_ID"]) == list(obj_b.meta_data["MULTI_ID"])


def test_requires_two_barcodes():
    rng = np.random.default_rng(0)
    counts = sp.csc_matrix(rng.poisson(5.0, size=(1, 20)).astype(float))
    cells = [f"c{i}" for i in range(20)]
    obj = create_shanuz_object(
        counts=counts, assay="HTO",
        feature_names=["HTO-A"], cell_names=cells,
        meta_data=pd.DataFrame(index=cells),
    )
    with pytest.raises(ValueError):
        multiseq_demux(obj)
