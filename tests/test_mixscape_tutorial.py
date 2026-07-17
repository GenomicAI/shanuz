"""Network-free tests for the Mixscape tutorial (thp1_mixscape_tutorial.py).

Two layers, no download required:

  * the pure helpers (perturbation_table / responsive_genes / call_concordance /
    variable_feature_list) on hand-built inputs, and
  * the full tutorial pipeline — load_screen_object → prep_reduction →
    run_signature → run_classification → run_lda → summarize — driven over a
    *synthetic* ECCITE-like screen with known knockout truth, by monkeypatching
    the dataset loader. This exercises the tutorial's wiring (assays, metadata
    columns, the report dict) without touching the network; the numerical
    fidelity of Mixscape itself lives in tests/test_mixscape.py.
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

matplotlib.use("Agg")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import tutorials.thp1_mixscape_tutorial as tut  # noqa: E402
from tutorials.thp1_mixscape_tutorial import (  # noqa: E402
    call_concordance,
    perturbation_table,
    responsive_genes,
    variable_feature_list,
)

N_GENES = 200
N_RESP = 24               # first N_RESP genes respond to a knockout
GUIDES = ["G1", "G2"]


# ----------------------------------------------------------------------
# Pure helpers
# ----------------------------------------------------------------------


def test_perturbation_table_counts_and_order():
    meta = pd.DataFrame({
        "gene": ["NT"] * 4 + ["G1"] * 5 + ["G2"] * 2,
        "mixscape_class.global": (
            ["NT"] * 4 + ["KO", "KO", "KO", "NP", "NP"] + ["NP", "NP"]
        ),
    })
    tbl = perturbation_table(meta)
    assert list(tbl.columns) == ["gene", "n_cells", "n_ko", "n_np", "ko_rate"]
    assert list(tbl["gene"]) == ["G1", "G2"]          # NT dropped, sorted by rate
    row = tbl.set_index("gene")
    assert row.loc["G1", "n_cells"] == 5
    assert row.loc["G1", "n_ko"] == 3
    assert row.loc["G1", "n_np"] == 2
    assert row.loc["G1", "ko_rate"] == pytest.approx(0.6)
    assert row.loc["G2", "ko_rate"] == pytest.approx(0.0)


def test_perturbation_table_respects_prtb_type():
    meta = pd.DataFrame({
        "gene": ["NT", "G1", "G1"],
        "mixscape_class.global": ["NT", "KD", "NP"],
    })
    tbl = perturbation_table(meta, prtb_type="KD")
    assert tbl.set_index("gene").loc["G1", "n_ko"] == 1
    assert tbl.set_index("gene").loc["G1", "ko_rate"] == pytest.approx(0.5)


def test_responsive_genes_threshold():
    tbl = pd.DataFrame({"gene": ["A", "B", "C"], "ko_rate": [0.9, 0.5, 0.2]})
    assert responsive_genes(tbl, 0.5) == ["A", "B"]
    assert responsive_genes(tbl, 0.95) == []


def test_call_concordance_and_errors():
    assert call_concordance(["a", "b", "c"], ["a", "b", "c"]) == 1.0
    assert call_concordance(["a", "b"], ["a", "x"]) == pytest.approx(0.5)
    with pytest.raises(ValueError):
        call_concordance([], [])
    with pytest.raises(ValueError):
        call_concordance(["a", "b"], ["a"])


# ----------------------------------------------------------------------
# Synthetic screen driven through the tutorial pipeline (no network)
# ----------------------------------------------------------------------


def _synthetic_screen(seed=0):
    """Return the loader's 6-tuple for a small ECCITE-like screen with KO truth.

    NT controls plus two guides whose knockout cells shift a shared block of
    response genes; a minority of each guide's cells escape (NP) and look like NT.
    """
    rng = np.random.default_rng(seed)
    base = rng.uniform(4.0, 12.0, size=N_GENES)
    up = np.arange(0, N_RESP // 2)
    down = np.arange(N_RESP // 2, N_RESP)

    cols, gene, guide_id, rep, truth = [], [], [], [], []

    def emit(lam, g, gid, r):
        depth = rng.uniform(0.85, 1.2)
        cols.append(rng.poisson(np.clip(lam * depth, 0.0, None)))
        gene.append(g)
        guide_id.append(gid)
        rep.append(r)

    for i in range(150):                               # non-targeting controls
        emit(base.copy(), "NT", f"NTg{i % 4}", f"rep{i % 2 + 1}")
        truth.append("NT")
    for gi, (g, up_set, down_set) in enumerate(
        [("G1", up, down), ("G2", down, up)]
    ):
        for i in range(80):
            lam = base.copy()
            is_ko = i < 55                             # 55 KO / 25 escapers
            if is_ko:
                lam[up_set] *= 4.0
                lam[down_set] *= 0.2
            emit(lam, g, f"{g}g{i % 3}", f"rep{i % 2 + 1}")
            truth.append("KO" if is_ko else "NP")

    counts = np.asarray(cols).T.astype(float)          # genes × cells
    n_cells = counts.shape[1]
    genes = [f"g{i}" for i in range(N_GENES)]
    cells = [f"c{i}" for i in range(n_cells)]
    meta = pd.DataFrame(
        {
            "gene": gene,
            "guide_ID": guide_id,
            "replicate": rep,
            "crispr": ["NT" if g == "NT" else "Perturbed" for g in gene],
            "Phase": "G1",
            "S.Score": rng.normal(size=n_cells),
            "G2M.Score": rng.normal(size=n_cells),
            "percent.mito": rng.uniform(0, 5, size=n_cells),
        },
        index=cells,
    )
    adt = sp.csc_matrix(rng.poisson(50, size=(4, n_cells)).astype(float))
    adt_names = ["CD86", "PDL1", "PDL2", "CD366"]
    return (sp.csc_matrix(counts), genes, adt, adt_names, meta,
            cells), np.array(truth)


@pytest.fixture(scope="module")
def screen():
    from _pytest.monkeypatch import MonkeyPatch

    mp = MonkeyPatch()
    data, truth = _synthetic_screen()
    mp.setattr(tut, "thp1_eccite", lambda data_dir=None: data)
    obj = tut.load_screen_object()
    hvg = tut.prep_reduction(obj)
    tut.run_signature(obj, hvg)
    tut.run_classification(obj)
    tut.run_lda(obj, npcs=5)
    yield obj, truth
    mp.undo()


def test_load_builds_rna_and_adt(monkeypatch):
    data, _ = _synthetic_screen()
    monkeypatch.setattr(tut, "thp1_eccite", lambda data_dir=None: data)
    obj = tut.load_screen_object()
    assert "RNA" in obj.assays and "ADT" in obj.assays
    assert obj.assays["ADT"].features() == ["CD86", "PDL1", "PDL2", "CD366"]
    assert "nCount_ADT" in obj.meta_data.columns
    assert set(obj.meta_data["gene"]) == {"NT", *GUIDES}


def test_prep_sets_variable_features(monkeypatch):
    data, _ = _synthetic_screen()
    monkeypatch.setattr(tut, "thp1_eccite", lambda data_dir=None: data)
    obj = tut.load_screen_object()
    hvg = tut.prep_reduction(obj)
    assert hvg == variable_feature_list(obj, "RNA")
    assert len(hvg) > 0
    assert obj.reductions["pca"].cell_embeddings.shape[1] >= tut.NDIMS


def test_pipeline_writes_mixscape_columns(screen):
    obj, _ = screen
    for col in ("mixscape_class", "mixscape_class.global", "mixscape_class_p_ko"):
        assert col in obj.meta_data.columns
    assert "PRTB" in obj.assays
    assert "lda" in obj.reductions


def test_pipeline_recovers_knockouts(screen):
    obj, truth = screen
    gclass = obj.meta_data["mixscape_class.global"].to_numpy()
    assert (gclass[truth == "NT"] == "NT").all()            # controls stay control
    assert (gclass == "KO").sum() > 0                        # some KO recovered
    # the strong synthetic effect should be caught for the majority of true KOs
    assert (gclass[truth == "KO"] == "KO").mean() >= 0.6


def test_summarize_report_dict(screen):
    obj, _ = screen
    out = tut.summarize(obj, verbose=False)
    assert out["n_cells"] == len(obj.cell_names())
    assert set(out["global"]) <= {"KO", "NP", "NT"}
    assert out["n_genes_tested"] == len(GUIDES)
    assert isinstance(out["table"], pd.DataFrame)
    assert list(out["table"]["gene"]) and set(out["table"]["gene"]) <= set(GUIDES)


def test_report_concordance_reads_r_calls(screen, tmp_path):
    obj, _ = screen
    cells = obj.cell_names()
    # Fabricate an R-calls CSV that agrees with Python everywhere → concordance 1.
    r = pd.DataFrame({
        "cell": cells,
        "R_mixscape_global": obj.meta_data["mixscape_class.global"].to_numpy(),
        "R_mixscape_class": obj.meta_data["mixscape_class"].to_numpy(),
    })
    path = tmp_path / "r_calls.csv"
    r.to_csv(path, index=False)
    agree = tut.report_concordance(obj, r_calls_path=path, verbose=False)
    assert agree["R_mixscape_global"] == pytest.approx(1.0)
    assert agree["R_mixscape_class"] == pytest.approx(1.0)


def test_report_concordance_absent_returns_none(screen, tmp_path):
    obj, _ = screen
    assert tut.report_concordance(obj, r_calls_path=tmp_path / "nope.csv",
                                  verbose=False) is None
