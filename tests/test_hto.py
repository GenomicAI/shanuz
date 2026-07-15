"""Tests for cell-hashing demultiplexing (v0.9.0 Specialized Assays).

  * hto_demux  (hto.py)

A synthetic hashtag matrix is built with *known* ground truth — clean singlets
for each of four hashtags, doublets spanning tag pairs, and empty-ish negatives —
and ``hto_demux`` is asked to recover it. The Seurat metadata columns, the
``hash.ID`` identity, and the singlet/doublet/negative accuracy are all checked
against that ground truth. Network-free and deterministic (fixed RNG + seed).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz.hto import hto_demux  # noqa: E402
from shanuz.preprocessing import normalize_data  # noqa: E402
from shanuz.shanuz import create_shanuz_object  # noqa: E402

TAGS = ["HTO-A", "HTO-B", "HTO-C", "HTO-D"]
DOUBLET_PAIRS = [(0, 1), (1, 2), (2, 3), (0, 2), (1, 3)]


def _hashing_counts(seed=0):
    """A hashtag count matrix (4 × N) plus per-cell ground-truth labels."""
    rng = np.random.default_rng(seed)
    H = len(TAGS)
    cols, truth = [], []

    # Background staining is ~uniform per antibody across all droplets; a cell is
    # "positive" for a tag only when its on-target signal towers over that floor.
    bg = 1.0
    for h in range(H):                       # 30 singlets per hashtag
        for _ in range(30):
            col = rng.poisson(bg, size=H).astype(float)
            col[h] = rng.poisson(100.0)
            cols.append(col)
            truth.append(TAGS[h])

    for _ in range(4):                       # 20 doublets across tag pairs
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


def test_writes_seurat_metadata_columns():
    obj, _ = _hashing_object()
    hto_demux(obj)
    for col in (
        "HTO_maxID", "HTO_secondID", "HTO_margin",
        "HTO_classification", "HTO_classification.global", "hash.ID",
    ):
        assert col in obj.meta_data.columns


def test_global_classification_domain():
    obj, _ = _hashing_object()
    hto_demux(obj)
    vals = set(obj.meta_data["HTO_classification.global"])
    assert vals <= {"Negative", "Singlet", "Doublet"}


def test_idents_set_to_hash_id():
    obj, _ = _hashing_object()
    hto_demux(obj)
    assert len(obj.idents) == len(obj.cell_names())
    assert list(obj.idents) == list(obj.meta_data["hash.ID"])


def test_margin_nonnegative_and_ids_distinct():
    obj, _ = _hashing_object()
    hto_demux(obj)
    margin = obj.meta_data["HTO_margin"].to_numpy()
    assert (margin >= -1e-9).all()
    max_id = obj.meta_data["HTO_maxID"].to_numpy()
    second_id = obj.meta_data["HTO_secondID"].to_numpy()
    assert (max_id != second_id).all()


# ----------------------------------------------------------------------
# recovery against ground truth
# ----------------------------------------------------------------------


def test_singlets_recovered():
    obj, truth = _hashing_object()
    hto_demux(obj)
    gclass = obj.meta_data["HTO_classification.global"].to_numpy()
    hash_id = obj.meta_data["hash.ID"].to_numpy()

    singlet = np.isin(truth, TAGS)
    assert (gclass[singlet] == "Singlet").mean() >= 0.9
    # correct singlets are assigned to their true hashtag
    assert (hash_id[singlet] == truth[singlet]).mean() >= 0.9


def test_doublets_recovered():
    obj, truth = _hashing_object()
    hto_demux(obj)
    gclass = obj.meta_data["HTO_classification.global"].to_numpy()
    doublet = truth == "Doublet"
    assert (gclass[doublet] == "Doublet").mean() >= 0.75
    # a called doublet's classification names its two tags, sorted + joined
    hto_class = obj.meta_data["HTO_classification"].to_numpy()
    for lab in hto_class[doublet & (gclass == "Doublet")]:
        assert "_" in lab and lab == "_".join(sorted(lab.split("_")))


def test_negatives_recovered():
    obj, truth = _hashing_object()
    hto_demux(obj)
    gclass = obj.meta_data["HTO_classification.global"].to_numpy()
    negative = truth == "Negative"
    assert (gclass[negative] == "Negative").mean() >= 0.8


# ----------------------------------------------------------------------
# bookkeeping, options & guards
# ----------------------------------------------------------------------


def test_cutoffs_stored_in_misc():
    obj, _ = _hashing_object()
    hto_demux(obj)
    info = obj.misc["hto_demux"]["HTO"]
    assert set(info["cutoffs"]) == set(TAGS)
    assert info["ncenters"] == len(TAGS) + 1
    assert all(np.isfinite(c) for c in info["cutoffs"].values())


def test_init_overrides_center_count():
    obj, _ = _hashing_object()
    hto_demux(obj, init=6)
    assert obj.misc["hto_demux"]["HTO"]["ncenters"] == 6


def test_normalize_false_uses_data_layer():
    obj, truth = _hashing_object()
    normalize_data(obj, normalization_method="CLR", margin=2, assay="HTO")
    hto_demux(obj, normalize=False)
    gclass = obj.meta_data["HTO_classification.global"].to_numpy()
    singlet = np.isin(truth, TAGS)
    assert (gclass[singlet] == "Singlet").mean() >= 0.9


def test_deterministic():
    obj_a, _ = _hashing_object()
    obj_b, _ = _hashing_object()
    hto_demux(obj_a, seed=7)
    hto_demux(obj_b, seed=7)
    assert list(obj_a.meta_data["hash.ID"]) == list(obj_b.meta_data["hash.ID"])


def test_requires_two_hashtags():
    rng = np.random.default_rng(0)
    counts = sp.csc_matrix(rng.poisson(5.0, size=(1, 20)).astype(float))
    cells = [f"c{i}" for i in range(20)]
    obj = create_shanuz_object(
        counts=counts, assay="HTO",
        feature_names=["HTO-A"], cell_names=cells,
        meta_data=pd.DataFrame(index=cells),
    )
    with pytest.raises(ValueError):
        hto_demux(obj)


def test_unsupported_kfunc_raises():
    obj, _ = _hashing_object()
    with pytest.raises(NotImplementedError):
        hto_demux(obj, kfunc="clara")
