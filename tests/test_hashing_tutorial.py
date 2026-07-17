"""Tests for the cell-hashing tutorial's helpers and demux plumbing.

Covers the pure species/concordance helpers exactly, and exercises
load_hashing_object → run_demux → summarize → the species figure on a tiny
*synthetic* barnyard object (the real GSE108313 download is monkeypatched away),
so the whole thing runs in CI without network. The demux algorithms themselves
are covered by tests/test_hto.py and tests/test_multiseq.py; here we check the
tutorial wires them up and recovers cleanly-separated singlets.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tutorials.pbmc_hashing_tutorial import (  # noqa: E402
    call_concordance, load_hashing_object, run_demux, short_tags,
    species_labels, species_mask, summarize,
)


# --------------------------------------------------------------------------- #
# Pure helpers
# --------------------------------------------------------------------------- #

def test_short_tags_strips_the_oligo_sequence():
    assert short_tags(["BatchA-AGGACCATCCAA", "BatchH-TATCACATCGGT", "Plain"]) == [
        "BatchA", "BatchH", "Plain"
    ]


def test_species_mask_upper_case_is_human():
    mask = species_mask(["MALAT1", "CD3E", "A1BG-AS1", "Xkr4", "mt-Nd1", "0610007N19Rik"])
    assert list(mask) == [True, True, True, False, False, False]


def test_species_labels_covers_every_branch():
    genes = ["HUMANA", "HUMANB", "Mousea", "Mouseb"]
    counts = sp.csc_matrix(np.array(
        [[10, 0, 5, 0],    # HUMANA  (human)
         [10, 0, 5, 0],    # HUMANB  (human)
         [0, 10, 5, 0],    # Mousea  (mouse)
         [0, 10, 5, 0]],   # Mouseb  (mouse)
        dtype=float,
    ))
    frac, label = species_labels(counts, genes)
    assert list(label) == ["human", "mouse", "mixed", "empty"]
    assert np.isclose(frac[0], 1.0)
    assert np.isclose(frac[1], 0.0)
    assert np.isclose(frac[2], 0.5)
    assert np.isnan(frac[3])


def test_call_concordance_fraction_and_length_guard():
    assert call_concordance(["a", "b", "c", "d"], ["a", "x", "c", "d"]) == 0.75
    with pytest.raises(ValueError):
        call_concordance(["a"], ["a", "b"])
    with pytest.raises(ValueError):
        call_concordance([], [])


# --------------------------------------------------------------------------- #
# load → demux → summarize on a synthetic barnyard (no network)
# --------------------------------------------------------------------------- #

_TAGS = ["BatchA", "BatchB", "BatchC", "BatchD"]


def _synthetic_raw(n_singlet=25, n_neg=20, seed=0):
    """Raw arrays shaped like ``pbmc_hashing()`` returns: a clean 4-hashtag panel.

    Each singlet is strongly positive for exactly one tag; negatives are flat
    background. RNA is human-dominant (a couple of mouse genes present so the
    species columns are populated). Returns the loader's 5-tuple.
    """
    rng = np.random.default_rng(seed)
    hto_names = [f"{t}-OLIGO{i}" for i, t in enumerate(_TAGS)]
    genes = ["GENEA", "GENEB", "GENEC", "Mousea", "Mouseb"]  # 3 human, 2 mouse

    cells, hto_cols, rna_cols, design = [], [], [], []
    ntag = len(_TAGS)
    for t in range(ntag):
        for k in range(n_singlet):
            col = np.full(ntag, 3.0)
            col[t] = 400.0 + rng.integers(0, 50)
            hto_cols.append(col)
            design.append(_TAGS[t])
            cells.append(f"s{t}_{k}")
    for k in range(n_neg):
        hto_cols.append(np.full(ntag, 3.0))
        design.append("Negative")
        cells.append(f"neg_{k}")

    n = len(cells)
    for _ in range(n):
        rna_cols.append([rng.poisson(60), rng.poisson(60), rng.poisson(60),
                         rng.poisson(1), rng.poisson(1)])

    hto = sp.csc_matrix(np.array(hto_cols).T)
    rna = sp.csc_matrix(np.array(rna_cols, dtype=float).T)
    return rna, genes, hto, hto_names, cells, np.array(design)


def test_load_hashing_object_wires_assay_and_species_metadata(monkeypatch):
    rna, genes, hto, hto_names, cells, _ = _synthetic_raw()
    monkeypatch.setattr(
        "tutorials.pbmc_hashing_tutorial.pbmc_hashing",
        lambda data_dir=None, force_download=False: (rna, genes, hto, hto_names, cells),
    )
    obj, tags = load_hashing_object(min_cells=0)

    assert tags == _TAGS                       # oligo suffix stripped
    assert "HTO" in obj.assays
    assert obj.assays["HTO"]._all_feature_names == _TAGS
    for col in ("nCount_HTO", "nCount_human", "nCount_mouse", "human_frac", "species"):
        assert col in obj.meta_data.columns
    # human-dominant synthetic RNA → every cell called human
    assert set(obj.meta_data["species"]) == {"human"}


def test_run_demux_recovers_clean_singlets_and_reports(monkeypatch):
    rna, genes, hto, hto_names, cells, design = _synthetic_raw()
    monkeypatch.setattr(
        "tutorials.pbmc_hashing_tutorial.pbmc_hashing",
        lambda data_dir=None, force_download=False: (rna, genes, hto, hto_names, cells),
    )
    obj, tags = load_hashing_object(min_cells=0)
    run_demux(obj, positive_quantile=0.99)

    for col in ("HTO_classification.global", "hash.ID", "HTO_maxID", "MULTI_ID"):
        assert col in obj.meta_data.columns

    summ = summarize(obj, tags, verbose=False)
    assert set(summ["hto_global"]) <= {"Singlet", "Doublet", "Negative"}
    assert set(summ["multi_global"]) <= {"Singlet", "Doublet", "Negative"}

    # Cleanly-separated singlets should be recovered with the right sample tag.
    hash_id = np.asarray(obj.meta_data["hash.ID"].values)
    singlet = design != "Negative"
    correct = np.mean(hash_id[singlet] == design[singlet])
    assert correct >= 0.7, f"only {correct:.0%} of clean singlets recovered"


def test_species_scatter_figure_has_three_classes(monkeypatch):
    import matplotlib.pyplot as plt

    from tutorials.generate_hashing_plots import species_scatter

    rna, genes, hto, hto_names, cells, _ = _synthetic_raw()
    monkeypatch.setattr(
        "tutorials.pbmc_hashing_tutorial.pbmc_hashing",
        lambda data_dir=None, force_download=False: (rna, genes, hto, hto_names, cells),
    )
    obj, tags = load_hashing_object(min_cells=0)
    run_demux(obj, positive_quantile=0.99)

    fig = species_scatter(obj)
    ax = fig.axes[0]
    assert ax.get_xlabel().startswith("log10 human")
    assert ax.get_ylabel().startswith("log10 mouse")
    # one scatter collection per global class present
    assert len(ax.collections) >= 1
    plt.close(fig)
