"""Tests for pooled CRISPR-screen analysis (v0.9.0 Specialized Assays).

  * calc_perturb_sig  (mixscape.py)
  * run_mixscape      (mixscape.py)
  * mixscape_lda      (mixscape.py)

A synthetic ECCITE-like screen is built with *known* ground truth: non-targeting
(NT) controls, plus two target genes whose cells are a mixture of true knockouts
(KO — a set of response genes shifted up/down) and non-perturbed escapers (NP —
indistinguishable from NT). The perturbation signature is computed and mixscape
is asked to recover, per cell, whether it is KO / NP / NT. The Seurat metadata
columns (``mixscape_class`` / ``.global`` / ``_p_ko``), the identity, the misc
bookkeeping, and the KO/NP/NT recovery are all checked. Network-free and
deterministic (fixed RNG + seeds).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shanuz.mixscape import calc_perturb_sig, mixscape_lda, run_mixscape  # noqa: E402
from shanuz.preprocessing import normalize_data, scale_data  # noqa: E402
from shanuz.reduction import run_pca  # noqa: E402
from shanuz.shanuz import create_shanuz_object  # noqa: E402

N_GENES = 80
N_RESP = 12          # first N_RESP genes respond to perturbation
GENES = ["G1", "G2"]


def _screen_counts(seed=0):
    """A gene × cell count matrix plus per-cell guide class and KO/NP/NT truth."""
    rng = np.random.default_rng(seed)
    base = rng.uniform(4.0, 12.0, size=N_GENES)
    up = np.arange(0, N_RESP // 2)             # response genes driven up in KO
    down = np.arange(N_RESP // 2, N_RESP)      # response genes driven down in KO

    cols, guide, truth = [], [], []

    def emit(lam):
        depth = rng.uniform(0.8, 1.25)         # per-cell depth, removed by PRTB
        cols.append(rng.poisson(np.clip(lam * depth, 0.0, None)))

    for _ in range(80):                        # non-targeting controls
        emit(base.copy())
        guide.append("NT")
        truth.append("NT")

    # Two guides with opposite response-gene directions (overlapping DE gene set).
    for gene, up_set, down_set in [("G1", up, down), ("G2", down, up)]:
        for i in range(50):
            lam = base.copy()
            is_ko = i < 35                     # 35 KO / 15 NP escapers per gene
            if is_ko:
                lam[up_set] *= 4.0
                lam[down_set] *= 0.2
            emit(lam)
            guide.append(gene)
            truth.append("KO" if is_ko else "NP")

    counts = np.asarray(cols).T.astype(float)  # genes × cells
    return counts, np.array(guide), np.array(truth)


def _screen_object(seed=0, split=False):
    counts, guide, truth = _screen_counts(seed)
    genes = [f"g{i}" for i in range(counts.shape[0])]
    cells = [f"c{i}" for i in range(counts.shape[1])]
    meta = pd.DataFrame({"gene": guide}, index=cells)
    if split:
        # interleaved replicate label, so every replicate carries NT + guide cells
        meta["rep"] = [f"rep{i % 2}" for i in range(len(cells))]
    obj = create_shanuz_object(
        counts=sp.csc_matrix(counts), assay="RNA",
        feature_names=genes, cell_names=cells, meta_data=meta,
    )
    normalize_data(obj)
    scale_data(obj)
    run_pca(obj, n_pcs=15)
    return obj, guide, truth


# ----------------------------------------------------------------------
# calc_perturb_sig
# ----------------------------------------------------------------------


def test_creates_prtb_assay():
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    assert "PRTB" in obj.assays
    prtb = obj.assays["PRTB"]
    assert prtb.features() == obj.assays["RNA"].features()
    assert prtb.cells() == obj.cell_names()


def test_perturb_sig_zeroes_controls():
    # An NT cell minus the mean of nearby NT cells should be near zero, whereas a
    # true KO cell keeps a large residual on the response genes.
    obj, guide, truth = _screen_object()
    calc_perturb_sig(obj)
    sig = obj.assays["PRTB"].layer_data("data")
    sig = np.asarray(sig.todense() if sp.issparse(sig) else sig)
    resp = np.arange(N_RESP)
    nt_mag = np.abs(sig[np.ix_(resp, np.where(truth == "NT")[0])]).mean()
    ko_mag = np.abs(sig[np.ix_(resp, np.where(truth == "KO")[0])]).mean()
    assert ko_mag > 3 * nt_mag


def test_calc_perturb_sig_requires_nt():
    obj, _, _ = _screen_object()
    obj.meta_data["gene"] = "G1"               # no NT cells left
    with pytest.raises(ValueError):
        calc_perturb_sig(obj)


# ----------------------------------------------------------------------
# run_mixscape — metadata, identity, domain
# ----------------------------------------------------------------------


def test_writes_mixscape_columns():
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    run_mixscape(obj)
    for col in ("mixscape_class", "mixscape_class.global", "mixscape_class_p_ko"):
        assert col in obj.meta_data.columns


def test_mixscape_class_domain():
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    run_mixscape(obj)
    allowed = {"NT"} | {f"{g} KO" for g in GENES} | {f"{g} NP" for g in GENES}
    assert set(obj.meta_data["mixscape_class"]) <= allowed
    assert set(obj.meta_data["mixscape_class.global"]) <= {"KO", "NP", "NT"}


def test_idents_set_to_mixscape_class():
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    run_mixscape(obj)
    assert list(obj.idents) == list(obj.meta_data["mixscape_class"])


# ----------------------------------------------------------------------
# run_mixscape — recovery against ground truth
# ----------------------------------------------------------------------


def test_nt_cells_stay_nt():
    obj, _, truth = _screen_object()
    calc_perturb_sig(obj)
    run_mixscape(obj)
    gclass = obj.meta_data["mixscape_class.global"].to_numpy()
    assert (gclass[truth == "NT"] == "NT").all()


def test_knockouts_recovered():
    obj, _, truth = _screen_object()
    calc_perturb_sig(obj)
    run_mixscape(obj)
    gclass = obj.meta_data["mixscape_class.global"].to_numpy()
    ko = truth == "KO"
    assert (gclass[ko] == "KO").mean() >= 0.8


def test_escapers_recovered():
    obj, _, truth = _screen_object()
    calc_perturb_sig(obj)
    run_mixscape(obj)
    gclass = obj.meta_data["mixscape_class.global"].to_numpy()
    npc = truth == "NP"
    assert (gclass[npc] == "NP").mean() >= 0.6


def test_recovery_with_split_by():
    obj, _, truth = _screen_object(split=True)
    calc_perturb_sig(obj, split_by="rep")
    run_mixscape(obj)
    gclass = obj.meta_data["mixscape_class.global"].to_numpy()
    assert (gclass[truth == "KO"] == "KO").mean() >= 0.8
    assert (gclass[truth == "NT"] == "NT").all()


# ----------------------------------------------------------------------
# run_mixscape — posterior, bookkeeping, options & guards
# ----------------------------------------------------------------------


def test_posterior_column():
    obj, _, truth = _screen_object()
    calc_perturb_sig(obj)
    run_mixscape(obj)
    p = obj.meta_data["mixscape_class_p_ko"].to_numpy()
    assert np.isnan(p[truth == "NT"]).all()          # NT cells get no posterior
    gene_p = p[truth != "NT"]
    assert np.all((gene_p >= -1e-9) & (gene_p <= 1 + 1e-9))


def test_bookkeeping_in_misc():
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    run_mixscape(obj)
    info = obj.misc["mixscape"]["PRTB"]
    assert info["nt_class"] == "NT"
    assert info["prtb_type"] == "KO"
    assert set(info["genes"]) == set(GENES)
    for g in GENES:
        assert info["genes"][g]["n_de"] >= 5           # response genes are found
        assert info["genes"][g]["n_ko"] > 0


def test_min_de_genes_forces_np():
    obj, _, truth = _screen_object()
    calc_perturb_sig(obj)
    run_mixscape(obj, min_de_genes=1000)               # unreachable → all NP
    gclass = obj.meta_data["mixscape_class.global"].to_numpy()
    assert (gclass[truth != "NT"] == "NP").all()
    assert not (gclass == "KO").any()


def test_min_cells_forces_np():
    obj, _, truth = _screen_object()
    calc_perturb_sig(obj)
    run_mixscape(obj, min_cells=1000)                  # every guide too small → NP
    gclass = obj.meta_data["mixscape_class.global"].to_numpy()
    assert not (gclass == "KO").any()


def test_prtb_type_labels():
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    run_mixscape(obj, prtb_type="KD")                  # knock-down screen naming
    assert "mixscape_class_p_kd" in obj.meta_data.columns
    assert set(obj.meta_data["mixscape_class.global"]) <= {"KD", "NP", "NT"}


def test_deterministic():
    obj_a, _, _ = _screen_object()
    obj_b, _, _ = _screen_object()
    calc_perturb_sig(obj_a)
    calc_perturb_sig(obj_b)
    run_mixscape(obj_a)
    run_mixscape(obj_b)
    assert list(obj_a.meta_data["mixscape_class"]) == list(obj_b.meta_data["mixscape_class"])


def test_run_mixscape_requires_nt():
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    obj.meta_data["gene"] = "G1"                        # no NT cells
    with pytest.raises(ValueError):
        run_mixscape(obj)


# ----------------------------------------------------------------------
# mixscape_lda
# ----------------------------------------------------------------------
#
# NPCS is held below the N_RESP response genes: a guide only contributes a block
# if it clears npcs + 1 DE genes, and the synthetic screen has ~12 to give.

NPCS = 5


def test_lda_creates_reduction():
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    mixscape_lda(obj, npcs=NPCS)
    assert "lda" in obj.reductions
    lda = obj.reductions["lda"]
    assert lda.cells() == obj.cell_names()
    # 3 guide classes (NT, G1, G2) → 2 discriminant directions
    assert lda.cell_embeddings.shape == (len(obj.cell_names()), 2)


def test_lda_only_needs_calc_perturb_sig():
    # MixscapeLDA groups by the raw guide label, so run_mixscape is not required.
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    mixscape_lda(obj, npcs=NPCS)                        # no run_mixscape call
    assert "mixscape_class" not in obj.meta_data.columns
    assert "lda" in obj.reductions


def test_lda_block_features():
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    mixscape_lda(obj, npcs=NPCS)
    lda = obj.reductions["lda"]
    # one npcs-wide block per contributing guide, named "<gene>_PC_<k>"
    assert lda.features() == [f"{g}_PC_{k + 1}" for g in GENES for k in range(NPCS)]
    assert lda.misc["genes_used"] == GENES
    assert lda.feature_loadings.shape == (len(GENES) * NPCS, 2)


def test_lda_separates_guide_classes():
    obj, _, truth = _screen_object()
    calc_perturb_sig(obj)
    mixscape_lda(obj, npcs=NPCS)
    guessed = obj.meta_data["lda_assignments"].to_numpy()
    guide = obj.meta_data["gene"].to_numpy()
    # NP escapers genuinely look like NT, so score only the separable cells.
    separable = truth != "NP"
    assert (guessed[separable] == guide[separable]).mean() >= 0.8


def test_lda_posterior_columns():
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    mixscape_lda(obj, npcs=NPCS)
    cols = [f"LDAP_{c}" for c in ["G1", "G2", "NT"]]
    for col in cols:
        assert col in obj.meta_data.columns
    post = obj.meta_data[cols].to_numpy()
    assert np.allclose(post.sum(axis=1), 1.0)
    assert ((post >= -1e-9) & (post <= 1 + 1e-9)).all()


def test_lda_assignment_domain():
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    mixscape_lda(obj, npcs=NPCS)
    assert set(obj.meta_data["lda_assignments"]) <= {"NT", *GENES}


def test_lda_preserves_idents():
    # run_mixscape leaves mixscape_class as the identity; LDA must hand it back.
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    run_mixscape(obj)
    before = list(obj.idents)
    mixscape_lda(obj, npcs=NPCS)
    assert list(obj.idents) == before


def test_lda_reduction_key():
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    mixscape_lda(obj, npcs=NPCS, reduction_name="prtb_lda", reduction_key="PLDA_")
    assert "prtb_lda" in obj.reductions
    assert obj.reductions["prtb_lda"].key == "PLDA_"


def test_lda_insufficient_de_raises():
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    with pytest.raises(ValueError):
        mixscape_lda(obj, npcs=1000)                    # no guide clears npcs + 1


def test_lda_requires_nt():
    obj, _, _ = _screen_object()
    calc_perturb_sig(obj)
    obj.meta_data["gene"] = "G1"                        # no NT cells
    with pytest.raises(ValueError):
        mixscape_lda(obj, npcs=NPCS)


def test_lda_deterministic():
    obj_a, _, _ = _screen_object()
    obj_b, _, _ = _screen_object()
    calc_perturb_sig(obj_a)
    calc_perturb_sig(obj_b)
    mixscape_lda(obj_a, npcs=NPCS)
    mixscape_lda(obj_b, npcs=NPCS)
    assert np.allclose(
        obj_a.reductions["lda"].cell_embeddings,
        obj_b.reductions["lda"].cell_embeddings,
    )
