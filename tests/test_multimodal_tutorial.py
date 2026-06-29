"""Tests for the multimodal CITE-seq tutorial's cell annotation.

Builds a tiny synthetic two-assay (RNA + ADT) object with cleanly-separated
signal and checks the combined protein-priority / RNA-fallback gating in
annotate_cells(). Network-free.
"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz.shanuz import create_shanuz_object  # noqa: E402
from shanuz.assay5 import create_assay5_object  # noqa: E402
from shanuz.preprocessing import normalize_data  # noqa: E402
from tutorials.cbmc_citeseq_tutorial import annotate_cells  # noqa: E402


def _two_assay_object(rna_levels, adt_levels, n=6):
    clusters = list(rna_levels)
    rna_genes = sorted({g for lv in rna_levels.values() for g in lv} | {"MALAT1"})
    adt_prots = sorted({p for lv in adt_levels.values() for p in lv})
    rg = {g: i for i, g in enumerate(rna_genes)}
    ag = {p: i for i, p in enumerate(adt_prots)}

    cells, idents = [], []
    rmat = np.zeros((len(rna_genes), len(clusters) * n))
    amat = np.zeros((len(adt_prots), len(clusters) * n))
    col = 0
    for cl in clusters:
        for _ in range(n):
            rmat[rg["MALAT1"], col] = 50
            for g, v in rna_levels[cl].items():
                rmat[rg[g], col] = v
            for p, v in adt_levels[cl].items():
                amat[ag[p], col] = v
            cells.append(f"cell_{col}")
            idents.append(cl)
            col += 1

    obj = create_shanuz_object(
        counts=sp.csc_matrix(rmat), assay="RNA",
        feature_names=rna_genes, cell_names=cells,
    )
    normalize_data(obj)
    obj.assays["ADT"] = create_assay5_object(
        counts=sp.csc_matrix(amat), feature_names=adt_prots,
        cell_names=cells, key="adt_",
    )
    normalize_data(obj, assay="ADT", normalization_method="CLR", margin=2)
    obj.idents = idents
    return obj


def test_annotate_cells_protein_then_rna_fallback():
    obj = _two_assay_object(
        rna_levels={
            "t": {}, "b": {}, "nk": {}, "mono": {},
            "plt": {"PPBP": 200, "PF4": 200},     # protein-panel-less → RNA
        },
        adt_levels={
            "t":   {"CD3": 200, "CD4": 200},
            "b":   {"CD19": 200},
            "nk":  {"CD16": 200, "CD56": 200},
            "mono": {"CD14": 200},
            "plt": {},
        },
    )
    anno = annotate_cells(obj)
    assert anno["t"] == "CD4 T"
    assert anno["b"] == "B"
    assert anno["nk"] == "NK"
    assert anno["mono"] == "CD14+ Mono"
    assert anno["plt"] == "Platelet"      # resolved by RNA, not protein


def test_annotate_cells_cd8_split():
    # A non-T cluster ("b") is needed so CD3 (margin=2 CLR) is actually
    # discriminative — otherwise a uniformly-high protein centers to ~0.
    obj = _two_assay_object(
        rna_levels={"cd4": {}, "cd8": {}, "b": {}},
        adt_levels={
            "cd4": {"CD3": 200, "CD4": 200, "CD8": 1},
            "cd8": {"CD3": 200, "CD8": 200, "CD4": 1},
            "b":   {"CD19": 200},
        },
    )
    anno = annotate_cells(obj)
    assert anno["cd4"] == "CD4 T"
    assert anno["cd8"] == "CD8 T"
    assert anno["b"] == "B"
