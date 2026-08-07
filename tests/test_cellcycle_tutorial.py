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
from truecell.module_score import CC_GENES  # noqa: E402


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
    can substitute for truecell.datasets.thp1_eccite. Cells fall in three groups:
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


# ---------------------------------------------------------------------------
# The deterministic regime and the multi-program call
# ---------------------------------------------------------------------------

def test_the_deterministic_regime_is_actually_deterministic(scored):
    """`nbin=1` + `ctrl=pool` must remove the seed from the answer.

    This is the premise the exact comparison rests on: with one bin and a draw
    that exhausts it, `sample(n, n)` is a permutation, so the control *set* is
    forced and only its summation order is random. If that stopped being true —
    a different default, a change to the pool — the exact comparison would
    quietly become a loose one, and nothing else here would notice.
    """
    from truecell import add_module_score

    _tut, obj, _summary, _tmp = scored
    genes = list(obj.assays["RNA"].features())
    program = [g for g in tut.IFN_PROGRAM if g in genes][:5]
    if len(program) < 2:
        pytest.skip("synthetic fixture carries too few IFN genes")

    seen = []
    for seed in (1, 4242):
        add_module_score(obj, features={"Probe": program},
                         nbin=1, ctrl=len(genes), seed=seed)
        seen.append(np.asarray(obj.meta_data["Probe"], dtype=float).copy())
    assert np.max(np.abs(seen[0] - seen[1])) <= tut.EXACT_TOLERANCE


def test_the_default_regime_is_not_deterministic(scored):
    """The control: at default settings the seed *does* move the score.

    Without this, a bug that ignored `seed` entirely would make the test above
    pass for the wrong reason.
    """
    from truecell import add_module_score

    _tut, obj, _summary, _tmp = scored
    genes = list(obj.assays["RNA"].features())
    program = [g for g in tut.IFN_PROGRAM if g in genes][:5]
    if len(program) < 2:
        pytest.skip("synthetic fixture carries too few IFN genes")

    seen = []
    for seed in (1, 4242):
        add_module_score(obj, features={"Probe2": program}, nbin=4, ctrl=3,
                         seed=seed)
        seen.append(np.asarray(obj.meta_data["Probe2"], dtype=float).copy())
    assert np.max(np.abs(seen[0] - seen[1])) > tut.EXACT_TOLERANCE


def test_scoring_writes_the_exact_and_multi_columns(scored):
    _tut, obj, _summary, _tmp = scored
    assert tut.EXACT_NAME in obj.meta_data.columns
    for i in (1, 2, 3):
        assert f"{tut.MULTI_NAME}{i}" in obj.meta_data.columns


def test_exact_comparison_flags_a_difference_above_tolerance(scored):
    """The exact check must be able to fail, not just to pass."""
    _tut, obj, _summary, tmp = scored
    meta, cells = obj.meta_data, obj.cell_names()
    base = dict(cell=cells, R_Phase=np.asarray(meta[tut.PHASE_COL]),
                R_S_Score=np.asarray(meta[tut.S_COL]),
                R_G2M_Score=np.asarray(meta[tut.G2M_COL]),
                R_IFN=np.asarray(meta[tut.IFN_NAME]))

    exact = np.asarray(meta[tut.EXACT_NAME], dtype=float)
    pd.DataFrame({**base, "R_IFN_exact": exact}).to_csv(tmp / "r_calls.csv",
                                                        index=False)
    assert tut.report_concordance(obj, verbose=False)["exact"]["within_tolerance"]

    nudged = exact.copy()
    nudged[0] += 1e-6            # far above EXACT_TOLERANCE, far below any score
    pd.DataFrame({**base, "R_IFN_exact": nudged}).to_csv(tmp / "r_calls.csv",
                                                          index=False)
    out = tut.report_concordance(obj, verbose=False)
    assert not out["exact"]["within_tolerance"]
    assert out["exact"]["max_abs_diff"] == pytest.approx(1e-6, rel=1e-6)


def test_transposed_multi_columns_are_detected(scored):
    """Two programs swapped must be caught.

    A transposition leaves every per-position correlation *high* whenever the
    programs score similarly, so `pearson > 0.9` on each column proves nothing.
    What catches it is that each Python column has to match its own R column
    better than it matches either of the others.
    """
    _tut, obj, _summary, tmp = scored
    meta, cells = obj.meta_data, obj.cell_names()
    cols = [np.asarray(meta[f"{tut.MULTI_NAME}{i}"], dtype=float) for i in (1, 2, 3)]
    base = dict(cell=cells, R_Phase=np.asarray(meta[tut.PHASE_COL]),
                R_S_Score=np.asarray(meta[tut.S_COL]),
                R_G2M_Score=np.asarray(meta[tut.G2M_COL]),
                R_IFN=np.asarray(meta[tut.IFN_NAME]),
                R_IFN_exact=np.asarray(meta[tut.EXACT_NAME], dtype=float))

    pd.DataFrame({**base, "R_Multi1": cols[0], "R_Multi2": cols[1],
                  "R_Multi3": cols[2]}).to_csv(tmp / "r_calls.csv", index=False)
    assert tut.report_concordance(obj, verbose=False)["multi_not_transposed"]

    # Programs 1 and 2 swapped on the R side.
    pd.DataFrame({**base, "R_Multi1": cols[1], "R_Multi2": cols[0],
                  "R_Multi3": cols[2]}).to_csv(tmp / "r_calls.csv", index=False)
    assert not tut.report_concordance(obj, verbose=False)["multi_not_transposed"]


def test_run_scoring_actually_uses_the_deterministic_settings(scored):
    """`run_scoring` must compute the EXACT column in the deterministic regime.

    Mutation testing caught this gap: dropping `nbin=1, ctrl=n_features` from
    `run_scoring` left every other test in this file green, because they all
    fabricate the R side *from Python's own column* — so a column computed the
    wrong way still matched itself. The property has to be asserted on the
    pipeline's output, not on a hand-made call with the right arguments.

    Seed-invariance is the property that *defines* the regime, but it cannot be
    the assertion here: this fixture carries few enough genes that ctrl=100
    exhausts every bin at the defaults too, so both columns come out
    seed-invariant (2.2e-15) and the check would not discriminate. Asserting
    that the EXACT column differs from a default-settings recomputation does
    discriminate, because nbin=1 draws a different control set.
    """
    from truecell import add_module_score

    _tut, obj, _summary, _tmp = scored
    tut.run_scoring(obj, seed=1)
    exact = np.asarray(obj.meta_data[tut.EXACT_NAME], dtype=float).copy()

    genes = set(obj.assays["RNA"].features())
    program = [g for g in tut.IFN_PROGRAM if g in genes]
    add_module_score(obj, features={"DfltExact": program}, seed=1)
    at_defaults = np.asarray(obj.meta_data["DfltExact"], dtype=float)

    assert np.max(np.abs(exact - at_defaults)) > tut.EXACT_TOLERANCE


def test_run_scoring_uses_non_default_settings_for_the_multi_call(scored):
    """The multi call must run at nbin=12/ctrl=40, not at the defaults.

    Also from mutation testing: reverting those arguments changed nothing any
    test could see. Recomputing the same three programs at the defaults and
    requiring a difference pins the settings themselves.
    """
    from truecell import add_module_score

    _tut, obj, _summary, _tmp = scored
    tut.run_scoring(obj, seed=1)
    stored = [np.asarray(obj.meta_data[f"{tut.MULTI_NAME}{i}"], dtype=float).copy()
              for i in (1, 2, 3)]

    genes = set(obj.assays["RNA"].features())
    programs = [[g for g in CC_GENES["s_genes"] if g in genes],
                [g for g in CC_GENES["g2m_genes"] if g in genes],
                [g for g in tut.IFN_PROGRAM if g in genes]]
    add_module_score(obj, features=programs, name="Dflt", seed=1)
    at_defaults = [np.asarray(obj.meta_data[f"Dflt{i}"], dtype=float)
                   for i in (1, 2, 3)]

    assert any(np.max(np.abs(s - d)) > tut.EXACT_TOLERANCE
               for s, d in zip(stored, at_defaults))
