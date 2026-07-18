"""Network-free tests for the cell-cycle tutorial (tutorials/thp1_cellcycle_tutorial.py).

Covers the pure metric helpers (phase concordance, score correlation, phase
distribution, scoreboard, gene resolution) directly, and drives the whole
pipeline on a small synthetic dataset with planted S-phase and G2/M-phase
populations — so cell_cycle_scoring has real phases to recover — never touching
the network or the real THP-1 download.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import tutorials.thp1_cellcycle_tutorial as tut  # noqa: E402
from tutorials.thp1_cellcycle_tutorial import (  # noqa: E402
    phase_concordance,
    score_correlation,
    phase_distribution,
    build_scoreboard,
)
from shanuz.module_score import CC_GENES  # noqa: E402


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

def test_phase_concordance_perfect_and_partial():
    p = np.array(["G1", "S", "G2M", "G1"])
    assert phase_concordance(p, p) == pytest.approx(1.0)
    assert phase_concordance(p, np.array(["G1", "S", "G1", "S"])) == pytest.approx(0.5)


def test_phase_concordance_rejects_bad_shapes():
    with pytest.raises(ValueError):
        phase_concordance(np.array(["G1"]), np.array(["G1", "S"]))
    with pytest.raises(ValueError):
        phase_concordance(np.array([]), np.array([]))


def test_score_correlation_tracks_and_returns_both():
    rng = np.random.default_rng(0)
    a = rng.normal(size=200)
    b = a + rng.normal(scale=0.01, size=200)     # nearly identical -> ~1.0
    out = score_correlation(a, b)
    assert set(out) == {"pearson", "spearman"}
    assert out["pearson"] > 0.99
    assert out["spearman"] > 0.99
    # An unrelated vector correlates far lower.
    assert score_correlation(a, rng.normal(size=200))["pearson"] < 0.5


def test_score_correlation_rejects_bad_shapes():
    with pytest.raises(ValueError):
        score_correlation(np.array([1.0]), np.array([1.0, 2.0]))


def test_phase_distribution_orders_and_counts_all_phases():
    phases = np.array(["G1", "G1", "S", "G2M", "S", "G1"])
    dist = phase_distribution(phases)
    assert list(dist["phase"]) == ["G1", "S", "G2M"]
    assert list(dist["n"]) == [3, 2, 1]
    assert dist["fraction"].sum() == pytest.approx(1.0)
    # A phase absent from the data still appears, at zero.
    only_g1 = phase_distribution(np.array(["G1", "G1"]))
    assert list(only_g1["n"]) == [2, 0, 0]


def test_build_scoreboard_column_order():
    rows = [{"metric": "S.Score", "pearson": 0.99, "spearman": 0.98}]
    board = build_scoreboard(rows)
    assert list(board.columns) == ["metric", "pearson", "spearman", "concordance"] or \
        list(board.columns) == ["metric", "pearson", "spearman"]
    assert list(board["metric"]) == ["S.Score"]


# ---------------------------------------------------------------------------
# Synthetic end-to-end pipeline (no network)
# ---------------------------------------------------------------------------

def _synthetic_thp1(seed=0):
    """A tiny THP-1 stand-in with planted S and G2/M cycling populations.

    Returns the loader's 6-tuple (rna, genes, adt, adt_names, meta, cells) so it
    can substitute for shanuz.datasets.thp1_eccite. Cells fall in three groups:
    S-cyclers (S-phase genes up), G2/M-cyclers (G2/M genes up), and resting — so
    cell_cycle_scoring has all three phases to recover. A third of cells also
    carry the interferon program, to exercise add_module_score.
    """
    rng = np.random.default_rng(seed)
    from tutorials.thp1_cellcycle_tutorial import IFN_PROGRAM

    s_genes = CC_GENES["s_genes"]
    g2m_genes = CC_GENES["g2m_genes"]
    filler = [f"g{i}" for i in range(400)]                 # control-bin pool
    genes = list(dict.fromkeys(s_genes + g2m_genes + IFN_PROGRAM + filler))
    G = len(genes)
    gidx = {g: i for i, g in enumerate(genes)}

    groups = {"Scyc": 70, "G2Mcyc": 70, "rest": 90}
    cols, cells = [], []
    c = 0
    for grp, n in groups.items():
        for _ in range(n):
            base = rng.gamma(0.3, size=G) + 0.05
            block = s_genes if grp == "Scyc" else g2m_genes if grp == "G2Mcyc" else []
            for g in block:
                base[gidx[g]] += 6.0
            if rng.random() < 0.3:
                for g in IFN_PROGRAM:
                    base[gidx[g]] += 4.0
            cols.append(rng.poisson(base * 4000.0 / base.sum()))
            cells.append(f"{grp}_{c}")
            c += 1
    rna = sp.csc_matrix(np.asarray(cols).T)                # G x n_cells
    adt = sp.csc_matrix((2, len(cells)))                   # dummy ADT
    meta = pd.DataFrame({"Phase": ["G1"] * len(cells)}, index=cells)  # dummy published
    return rna, genes, adt, ["adt1", "adt2"], meta, cells


@pytest.fixture(scope="module")
def scored():
    """Run the full pipeline once on the synthetic dataset, gene lists in a tmp dir."""
    data = _synthetic_thp1()
    mp = pytest.MonkeyPatch()
    mp.setattr(tut, "thp1_eccite", lambda data_dir=None: data)
    tmp = Path(tempfile.mkdtemp())
    mp.setattr(tut, "FIGURES", tmp)
    obj, summary = tut.run_full(verbose=False)
    yield tut, obj, summary, tmp
    mp.undo()


def test_scoring_writes_all_columns(scored):
    _tut, obj, _summary, _tmp = scored
    for col in (tut.PHASE_COL, tut.S_COL, tut.G2M_COL, tut.IFN_NAME):
        assert col in obj.meta_data.columns
    assert set(np.asarray(obj.meta_data[tut.PHASE_COL])).issubset(set(tut.PHASES))


def test_scoring_recovers_planted_phases(scored):
    _tut, _obj, summary, _tmp = scored
    dist = summary["phase_distribution"].set_index("phase")
    # The planted S and G2/M cyclers must surface as non-empty phases.
    assert dist.loc["S", "n"] > 0
    assert dist.loc["G2M", "n"] > 0
    assert summary["n_cells"] == 230


def test_writes_resolved_gene_lists(scored):
    _tut, obj, _summary, tmp = scored
    for fname, requested in (("s_genes.txt", CC_GENES["s_genes"]),
                             ("g2m_genes.txt", CC_GENES["g2m_genes"])):
        assert (tmp / fname).exists()
        written = (tmp / fname).read_text().split()
        present = set(obj.assays["RNA"].features())
        # Every written gene is one that was requested AND is in the assay.
        assert written == [g for g in requested if g in present]


def test_report_concordance_reads_fabricated_r_calls(scored):
    _tut, obj, _summary, tmp = scored
    meta = obj.meta_data
    cells = obj.cell_names()
    # Fabricate R calls equal to Python's -> concordance 1.0, correlations 1.0.
    df = pd.DataFrame({
        "cell": cells,
        "R_Phase": np.asarray(meta[tut.PHASE_COL]),
        "R_S_Score": np.asarray(meta[tut.S_COL]),
        "R_G2M_Score": np.asarray(meta[tut.G2M_COL]),
        "R_IFN": np.asarray(meta[tut.IFN_NAME]),
    })
    df.to_csv(tmp / "r_calls.csv", index=False)
    out = tut.report_concordance(obj, verbose=False)
    assert out is not None
    assert out["phase_concordance"] == pytest.approx(1.0)
    assert out["s_score"]["pearson"] == pytest.approx(1.0, abs=1e-9)
    assert out["g2m_score"]["pearson"] == pytest.approx(1.0, abs=1e-9)


def test_report_concordance_absent_returns_none(scored):
    _tut, obj, _summary, tmp = scored
    missing = tmp / "nope.csv"
    assert tut.report_concordance(obj, r_calls_path=missing, verbose=False) is None
