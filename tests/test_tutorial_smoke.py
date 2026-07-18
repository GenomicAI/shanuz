"""End-to-end tutorial runs.

**These do not run in CI, and are not meant to be mistaken for coverage.** They
need the cached datasets (~200 MB under ``~/.shanuz_data``) and take minutes, so
they are gated behind an explicit opt-in:

    SHANUZ_TUTORIAL_SMOKE=1 pytest tests/test_tutorial_smoke.py -v

The gate is an env var rather than a bare "skip if the data is missing" check on
purpose. A test that silently skips wherever the data happens to be absent reads
as green in CI while proving nothing — which is the failure mode that let
``pbmc3k_tutorial.py`` ship broken in the first place. Requiring the opt-in means
a skip here always means "nobody asked for this", never "this passed".

What they cover that ``test_tutorial_marker_tables.py`` cannot: the actual script
path, top to bottom, against real data. The pandas 3 regression lived in the
middle of ``run_tutorial()`` — importable helpers and unit fixtures would not
have reached it. Worth running before cutting a release.
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = Path(os.path.expanduser("~/.shanuz_data"))

OPT_IN = os.environ.get("SHANUZ_TUTORIAL_SMOKE") == "1"

# script -> the dataset directory it needs cached
TUTORIALS = [
    ("pbmc3k_tutorial.py", "pbmc3k"),
    ("pbmc8k_subclustering_tutorial.py", "pbmc8k"),
    ("cbmc_citeseq_tutorial.py", "cbmc"),
    ("pbmc3k_sctransform_tutorial.py", "pbmc3k"),
    ("pbmc_hashing_tutorial.py", "pbmc_hashing"),
    ("thp1_mixscape_tutorial.py", "thp1_eccite"),
    ("ifnb_integration_tutorial.py", "ifnb"),
    ("panc8_reference_mapping_tutorial.py", "panc8"),
    ("thp1_cellcycle_tutorial.py", "thp1_eccite"),
    ("pbmc3k_dimreduc_tutorial.py", "pbmc3k"),
    ("ifnb_sketch_tutorial.py", "ifnb"),
    ("pbmc3k_objects_tutorial.py", "pbmc3k"),
]

pytestmark = pytest.mark.skipif(
    not OPT_IN,
    reason="set SHANUZ_TUTORIAL_SMOKE=1 to run the tutorials end-to-end "
           "(needs ~200MB of cached datasets, takes minutes)",
)


@pytest.mark.parametrize("script,dataset", TUTORIALS)
def test_tutorial_runs_to_completion(script, dataset):
    """The script exits 0 and reaches its final section.

    Deliberately asserts on the exit code and the tail of stdout rather than on
    any number the tutorial prints: cluster counts and marker rankings are
    library-version dependent (see the R-comparison notes), and pinning them here
    would turn a dependency bump into a test failure that says nothing about the
    tutorial. What this pins is that the script *runs*.
    """
    if not (DATA_ROOT / dataset).is_dir():
        pytest.skip(f"dataset {dataset!r} not cached under {DATA_ROOT}")

    proc = subprocess.run(
        [sys.executable, f"tutorials/{script}"],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=1800,
    )
    assert proc.returncode == 0, (
        f"tutorials/{script} exited {proc.returncode}\n"
        f"--- last 40 lines of stdout ---\n{os.linesep.join(proc.stdout.splitlines()[-40:])}\n"
        f"--- stderr ---\n{proc.stderr[-3000:]}"
    )
    assert proc.stdout.strip(), f"tutorials/{script} produced no output"


def test_ifnb_rpca_matches_seurat_batch_mixing():
    """Bug 2 regression, on real data: RPCA must actually integrate ifnb.

    Before the reciprocal-PCA fix, shanuz RPCA reached batch-mixing entropy 0.22
    on ifnb against Seurat's 0.91 — it barely mixed the CTRL/STIM batches and
    recovered cell type *below* the uncorrected baseline. The fix (Seurat's
    per-object scaling + the stacked-embedding SD/L2 normalisation) lifts it to
    ~0.87. This is the only test that reproduces the *emergent* under-integration:
    it needs the real dataset's many overlapping cell types and pervasive IFN
    response, which no synthetic fixture stands in for (both were checked and the
    pre-fix code integrated them fine). Asserts RPCA clears a floor the pre-fix
    code could not, that CCA (the working baseline) stays high, and that RPCA now
    recovers cell type better than doing nothing.
    """
    if not (DATA_ROOT / "ifnb").is_dir():
        pytest.skip("dataset 'ifnb' not cached")

    sys.path.insert(0, str(REPO_ROOT))
    from tutorials.ifnb_integration_tutorial import run_full

    _, summary = run_full(methods=("cca", "rpca"), do_umap=False, verbose=False)
    board = summary["scoreboard"].set_index("method")
    assert board.loc["rpca", "batch_entropy"] >= 0.7            # pre-fix was 0.22
    assert board.loc["cca", "batch_entropy"] >= 0.95            # baseline holds
    assert (board.loc["rpca", "ari_celltype"]
            > board.loc["uncorrected (PCA)", "ari_celltype"])   # pre-fix was below


def test_panc8_reference_mapping_recovers_celltype():
    """Reference mapping on real data: label transfer must recover the query's types.

    On panc8, transferring ``celltype`` from the CEL-seq2 reference to the
    SMART-seq2 query reaches ~0.98 accuracy against ground truth (Seurat ~0.99),
    with the abundant cell types recovered near-perfectly and high-confidence
    calls. Asserts the transfer clears a high floor, that it found a healthy
    anchor set, and that the mean prediction score stays high — the observable
    signatures a broken projection or anchor step would sink.
    """
    if not (DATA_ROOT / "panc8").is_dir():
        pytest.skip("dataset 'panc8' not cached")

    sys.path.insert(0, str(REPO_ROOT))
    from tutorials.panc8_reference_mapping_tutorial import run_full

    _ref, _query, anchors, _pred, summary = run_full(do_umap=False, verbose=False)
    assert summary["accuracy"] >= 0.9                        # observed 0.985
    assert len(anchors.anchors) > 500                        # observed ~3550
    assert summary["scoreboard"]["mean_score"].iloc[0] >= 0.9


def test_thp1_cellcycle_recovers_phases():
    """Cell-cycle scoring on real data: THP-1 must show real S/G2M populations.

    THP-1 is a proliferating line, so cell_cycle_scoring should place a
    substantial minority of cells in S and G2/M — not collapse everything to G1
    (which a broken module score, e.g. all-zero, would do since G1 is the
    both-scores-<=0 default). Also checks agreement with Papalexi's published
    phase, a signal that the scoring is producing sensible values. Observed:
    G1 0.72 / S 0.15 / G2M 0.13, published concordance 0.88.
    """
    if not (DATA_ROOT / "thp1_eccite").is_dir():
        pytest.skip("dataset 'thp1_eccite' not cached")

    sys.path.insert(0, str(REPO_ROOT))
    from tutorials.thp1_cellcycle_tutorial import run_full

    _obj, summary = run_full(verbose=False)
    dist = summary["phase_distribution"].set_index("phase")
    assert dist.loc["S", "fraction"] >= 0.08            # observed 0.15
    assert dist.loc["G2M", "fraction"] >= 0.06          # observed 0.13
    assert dist.loc["G1", "fraction"] <= 0.90           # not collapsed to all-G1
    assert summary["published_concordance"] >= 0.75     # observed 0.88


def test_pbmc3k_dimreduc_extras_hold_up():
    """JackStraw / ICA / t-SNE on real data: the observable signatures must survive.

    Deliberately *not* an assertion on ``score_jackstraw``'s per-PC scores. On
    real data those saturate at 0.0 for every PC — the finding this tutorial
    reports — and pinning the current behaviour here would cement the defect
    into the suite, which is exactly how the CLR and SCT bugs stayed green for
    so long. What is asserted instead holds both before and after any fix: the
    per-feature p-values are valid probabilities, the leading PCs carry more
    signal than the trailing ones, and t-SNE preserves its input's structure.
    """
    if not (DATA_ROOT / "pbmc3k").is_dir():
        pytest.skip("dataset 'pbmc3k' not cached")

    sys.path.insert(0, str(REPO_ROOT))
    from tutorials.pbmc3k_dimreduc_tutorial import run_full

    obj, summary = run_full(verbose=False)
    emp = obj.reductions["pca"].jackstraw.empirical_p_values
    assert np.isfinite(emp).all()
    assert emp.min() >= 0.0 and emp.max() <= 1.0

    # PC 1-5 carry the real cell-type structure; PC 16-20 are largely noise.
    n_sig = (emp <= 1e-5).sum(axis=0)
    assert n_sig[:5].mean() > n_sig[-5:].mean()

    assert summary["n_ics"] == 20
    assert summary["tsne_knn_vs_pca"] >= 0.3     # observed ~0.47; random is ~0.01


def test_ifnb_leverage_tracks_rarity():
    """The defect regression, on the data where the method's purpose is visible.

    Leverage sampling exists to keep rare cell states. On ifnb the annotations
    span 4,362 cells (CD14 Mono) down to 55 (Eryth), so "does leverage rise as a
    population gets rarer" is directly measurable — Spearman between a type's
    mean leverage and its size, which is -0.93 in both tools.

    This is the check no synthetic fixture can stand in for. Poisson clusters
    around a shared baseline do *not* reproduce it (several were tried, and R
    agrees with shanuz on those to 1e-5 while showing no enrichment either): real
    rare types are transcriptionally extreme, not merely scarce. Before the fix
    the full-rank scores were nearly flat — max/median 1.3 against R's 6.5 — so
    sampling by them was uniform sampling, and this correlation is what that
    destroys.
    """
    if not (DATA_ROOT / "ifnb").is_dir():
        pytest.skip("dataset 'ifnb' not cached")

    sys.path.insert(0, str(REPO_ROOT))
    from tutorials.ifnb_sketch_tutorial import run_full

    _, summary = run_full(verbose=False, do_projection=False, do_lazy=False)

    exact = summary["leverage"]["exact"]
    assert np.isfinite(exact).all() and exact.min() >= 0
    # Seurat truncates at 50 components, so the exact regime's scores sum to 50.
    assert exact.sum() == pytest.approx(50.0, abs=1e-3)
    # Real spread, not the near-uniform weights the old implementation produced.
    assert exact.max() / np.median(exact) > 2.5

    # The payoff: rarer types score higher. Observed -0.93 in both tools.
    assert summary["rarity_spearman"]["exact"] < -0.5

    # And the sketch acts on it — the rarest types are over-represented against
    # a same-size uniform draw, which is the control that makes this meaningful.
    lev = summary["sketch_composition"]["LeverageScore"]
    rarest = lev.tail(3).index
    assert lev.loc[rarest, "fold"].mean() > 1.0


@pytest.mark.parametrize("script,dataset", [TUTORIALS[0]])
def test_pbmc3k_prints_the_marker_table(script, dataset):
    """The exact block that the pandas 3 regression killed.

    ``run_tutorial`` crashed with KeyError: 'cluster' immediately after printing
    the "Top 3 markers per cluster:" header, so a run that reaches the header but
    dies right after still looks half-right in a log. Assert a cluster line
    actually follows it.
    """
    if not (DATA_ROOT / dataset).is_dir():
        pytest.skip(f"dataset {dataset!r} not cached under {DATA_ROOT}")

    proc = subprocess.run(
        [sys.executable, f"tutorials/{script}"],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=1800,
    )
    assert proc.returncode == 0, proc.stderr[-3000:]
    lines = proc.stdout.splitlines()
    header = next(
        (i for i, ln in enumerate(lines) if "Top 3 markers per cluster" in ln), None
    )
    assert header is not None, "tutorial never reached the marker table"
    cluster_lines = [
        ln for ln in lines[header + 1: header + 12] if ln.strip().startswith("Cluster ")
    ]
    assert len(cluster_lines) >= 8, (
        f"expected a line per cluster after the header, got {cluster_lines!r}"
    )
    # "Cluster 0: GENE, GENE, GENE" — genes present, not an empty table
    assert all(":" in ln and ln.split(":", 1)[1].strip() for ln in cluster_lines)


def test_pbmc3k_object_model_round_trips():
    """The object layer on real data, where the anchors are exact.

    Nothing here is stochastic, so these are equalities rather than thresholds.
    The split/join round trip is the one that matters: before this tutorial it
    returned a differently-named layer holding the right numbers in the wrong
    order, with the assay's own cell vector still in the original order — a
    silent misalignment of every column against the metadata indexing it.
    """
    if not (DATA_ROOT / "pbmc3k").is_dir():
        pytest.skip("dataset 'pbmc3k' not cached")

    sys.path.insert(0, str(REPO_ROOT))
    from tutorials.pbmc3k_objects_tutorial import run_full

    _, anchors = run_full(verbose=False)

    round_trip = anchors["split_join"]
    assert round_trip["layers_after_split"] == ["counts.batch1", "counts.batch2"]
    assert round_trip["layer_name_restored"]
    assert round_trip["cell_order_restored"]
    assert round_trip["matrix_restored"]

    # The no-argument join on a prepared assay — the call every script makes.
    assert anchors["join_all_layers"]["error"] is None

    # FetchData reaches all three kinds of variable, and returns numbers.
    fetch = anchors["fetch"]
    assert fetch["columns"] == ["nCount_RNA", "nFeature_RNA", "CD3E", "PC_1"]
    assert fetch["n_rows"] == anchors["shape"]["n_cells"]
    assert isinstance(fetch["CD3E_sum"], float)
    assert fetch["pc1_matches_embeddings"]

    # The command log, keyed the way Seurat keys it.
    assert anchors["commands"] == [
        "NormalizeData.RNA",
        "FindVariableFeatures.RNA",
        "ScaleData.RNA",
        "RunPCA.RNA",
        "FindNeighbors.RNA.pca",
    ]
