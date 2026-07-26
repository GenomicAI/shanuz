"""The differential-expression test suite — shanuz against Seurat 5.5.1.

Wave 3's first tutorial, and the last big untested surface in the library:
``find_markers`` offers eight statistical tests and **none of them had ever been
compared to R**. Their unit tests assert self-consistency on synthetic fixtures,
which is exactly the shape of coverage that let the CLR and SCTransform defects
survive.

The comparison runs both tools on the **same cells** — Python clusters pbmc3k and
writes the cluster-0-vs-cluster-1 assignment to ``figures_de/groups.csv``, which
the R side reads — so nothing here is measuring a clustering difference. What is
left is the test.

What it found
-------------
Two defects, both fixed here.

1. **``avg_log2FC`` put the pseudocount in the wrong place.** Seurat 5's
   ``log1pdata.mean.fxn`` is ``log2((sum(expm1(x)) + 1) / n)`` — one pseudocount
   added to the group's *total*, worth ``1/n`` on the mean scale. shanuz computed
   ``log2(mean(expm1(x)) + 1)``, adding a whole count to the *mean*. That floors
   every fold change near zero: a gene detected in 0 % of one cluster and 24 % of
   the other read **-1.26** where Seurat reads **-9.92**.

   This is the most-read column in any DE table, and it is worse than a display
   problem: ``logfc_threshold`` **filters on it**, so the error changed which
   genes came back at all. At Seurat's own default of 0.1, shanuz returned 4,903
   genes where Seurat returned 13,009 (Jaccard 0.377); at 0.25, 2,298 against
   11,931 (Jaccard **0.193**). Fewer than one gene in five agreed.

   Telling: where both groups actually express the gene the two formulas nearly
   agree (Spearman 0.990) — the error is concentrated in sparse, marker-like
   genes, which is precisely what differential expression is looking for. After
   the fix, **7.1e-15 across all 13,712 genes**.

2. **``negbinom`` was a different test.** Seurat's ``GLMDETest`` fits
   ``MASS::glm.nb`` — dispersion estimated by **maximum likelihood** — and reads
   the **Wald** p-value off the group coefficient. shanuz used a fixed
   method-of-moments dispersion and a **likelihood-ratio** test: a different
   estimator *and* a different statistic. It read HLA-DRA at 5.5e-128 against
   R's 1.1e-321. After the fix the p-values agree **exactly** (median |log10
   ratio| 0.000) for every gene detected above 5 %; what disagreement remains
   sits below that, where the negative-binomial GLM is fitting almost-empty rows
   and Seurat's own ``min.cells.feature`` drops the genes anyway.

Differences left standing, and why
----------------------------------
* **``deseq2`` is not Seurat's DESeq2.** Seurat's ``DESeq2DETest`` builds a
  ``DESeqDataSet`` with **one column per cell** and tests cells as replicates.
  shanuz sums counts per sample and tests at the sample level. Treating cells as
  replicates is the practice Squair et al. (2021) showed inflates false
  positives, so the pseudobulk route is the better statistics — and because it
  **requires ``sample_col``**, it cannot silently be mistaken for the per-cell
  test; it raises instead. Reported here rather than "fixed" in either direction.
* **``mast`` is a hand-rolled hurdle model**, not a call to the MAST package —
  which is not installable as a Python dependency. Spearman 0.947 on p-values and
  the same top 50 genes.
* **Seurat rounds ``myAUC`` to three decimals** inside ``DifferentialAUC``, so
  the ROC comparison is exact only to 5e-4. That is R's rounding, not a
  divergence — worth stating, because it looks like one.
* **R's ``wilcox`` returns ``NaN``** for genes with no expression in either
  group; shanuz returns ``p = 1``. A test that cannot be run has no evidence
  against the null, so 1 is the more useful answer, and R's NaN set is a subset
  of shanuz's p=1 set.

Usage
-----
    python tutorials/pbmc3k_de_tutorial.py
    Rscript tutorials/pbmc3k_de_verify.R      # writes figures_de/r_<test>.csv
    python tutorials/pbmc3k_de_tutorial.py --report
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from shanuz import create_shanuz_object, percentage_feature_set
from shanuz.clustering import find_clusters
from shanuz.datasets import pbmc3k
from shanuz.markers import find_markers
from shanuz.neighbors import find_neighbors
from shanuz.preprocessing import find_variable_features, normalize_data, scale_data
from shanuz.reduction import run_pca
from tutorials.bands import (
    Band, check_bands, check_shared_groups, render_verdicts,
)

FIGURES = Path(__file__).parent / "figures_de"
IDENT_1, IDENT_2 = "0", "1"
TOP_N = 50

# shanuz test -> Seurat's spelling of the same test.
TEST_MAP = {
    "wilcox": "wilcox",
    "t": "t",
    "bimod": "bimod",
    "LR": "LR",
    "negbinom": "negbinom",
    "roc": "roc",
    "mast": "MAST",
    "deseq2": "DESeq2",
}

# Seurat rounds myAUC to 3 dp in DifferentialAUC, so the ROC comparison cannot be
# tighter than that however correct both sides are. Everything else is compared
# on its own terms; no blanket tolerance.
AUC_TOLERANCE = 5e-4
# avg_log2FC is pure arithmetic on the shared matrix — it should be exact.
LOG2FC_TOLERANCE = 1e-12

# --- Declared bands --------------------------------------------------------
# The parity table below used to live only in de_vignette.md, so nothing failed
# when a number moved: `deseq2`'s top-50 overlap was described as "25/50" in
# prose and had drifted to 22 without anyone noticing, and the same silence
# would have covered a real regression in any other row.
#
# Seven of the eight tests are ports of the same statistic over the same cells,
# so their bands are exact. `deseq2` is the one measurement of a *deliberate*
# divergence — shanuz aggregates counts per sample and tests samples, Seurat's
# `DESeq2DETest` tests cells as replicates — so its band comes from measurement:
# resampling the pseudo-replicate split 20 times moves the overlap over 20-26
# (median 22), and the previous cluster assignment gave 25.
RANKED_TESTS = ("wilcox", "t", "bimod", "LR", "negbinom", "mast")

_PARITY_TOP50 = (
    "The same statistic on the same cells: the 50 most significant genes must "
    "be the same 50. A single dropped gene here is a regression, not scatter — "
    "this has read 50/50 on two different cluster assignments.")

BANDS: dict[str, Band] = {
    **{f"{t} top50": Band(TOP_N, TOP_N, _PARITY_TOP50, fmt=".0f")
       for t in RANKED_TESTS},
    "deseq2 top50": Band(
        15, 32,
        "A divergence measurement, not a parity target: shanuz tests pseudobulk "
        "samples, Seurat's DESeq2DETest tests cells. 20-26 over 20 resampled "
        "replicate splits and 25 on the previous cluster assignment. Bounded "
        "well below 50 on purpose — reaching parity would mean the pseudobulk "
        "aggregation had stopped happening and `sample_col` was being ignored.",
        fmt=".0f"),
    **{f"{t} rho>5%": Band(low, 1.0, why) for t, low, why in (
        ("wilcox", 0.9999,
         "Identical rank-sum statistic; measured exactly 1.0."),
        ("t", 0.9999, "Identical Welch t; measured exactly 1.0."),
        ("bimod", 0.9999, "Identical likelihood-ratio test; measured exactly 1.0."),
        ("LR", 0.9999, "Identical logistic-regression LRT; measured exactly 1.0."),
        ("negbinom", 0.88,
         "Both fit a negative-binomial GLM but not with the same optimiser, so "
         "the agreement is high rather than exact: 0.9165 here, 0.9194 before."),
        ("mast", 0.99,
         "shanuz's hurdle model is hand-rolled rather than a call to the MAST "
         "package, so this is the closest a reimplementation gets: 0.9979."),
    )},
    "deseq2 rho>5%": Band(
        0.12, 0.30,
        "Pseudobulk against per-cell, so a low correlation is the expected "
        "result and a high one would be the surprise. 0.1902-0.2021 over the 20 "
        "resampled splits.", ),
    "max |dlog2FC| (parity tests)": Band(
        0, LOG2FC_TOLERANCE,
        "avg_log2FC is arithmetic on the shared matrix with no statistics in "
        "it, so every test that shares Seurat's cell-level definition must "
        "agree to floating point: 6.4e-15 measured.", fmt=".2e"),
    "roc max |dAUC|": Band(
        0, AUC_TOLERANCE,
        "Seurat rounds myAUC to three decimals inside DifferentialAUC, so this "
        "cannot be tighter than half a unit in the last place however correct "
        "both sides are.", fmt=".2e"),
}


def build(data_dir=None):
    """pbmc3k through the standard pipeline, clustered."""
    counts, genes, cells = pbmc3k(data_dir=data_dir)
    obj = create_shanuz_object(counts=counts, assay="RNA", min_cells=3,
                               min_features=200, project="pbmc3k_de",
                               feature_names=genes, cell_names=cells)
    percentage_feature_set(obj, pattern=r"^MT-", col_name="percent.mt")
    md = obj.meta_data
    keep = (md["nFeature_RNA"] > 200) & (md["nFeature_RNA"] < 2500) & (md["percent.mt"] < 5)
    obj = obj.subset(cells=list(md.index[keep]))

    normalize_data(obj)
    find_variable_features(obj, selection_method="vst", nfeatures=2000)
    scale_data(obj, features=obj.assays["RNA"]._all_feature_names)
    run_pca(obj, n_pcs=50, features=obj.assays["RNA"].variable_features,
            reduction_name="pca")
    find_neighbors(obj, dims=range(10), k_param=20)
    find_clusters(obj, resolution=0.5, algorithm=1, random_seed=0)
    return obj


def shared_groups(obj) -> pd.Series:
    """The two clusters both tools will test, written out for the R side.

    Exported rather than re-derived because Louvain numbering is not guaranteed
    to agree across implementations, and a clustering difference would masquerade
    as a DE difference.
    """
    idents = pd.Series([str(i) for i in obj.idents],
                       index=list(obj.assays["RNA"].cells()), name="group")
    return idents[idents.isin([IDENT_1, IDENT_2])]


def run_tests(obj, groups: pd.Series, tests=None) -> dict[str, pd.DataFrame]:
    sub = obj.subset(cells=list(groups.index))
    sub.idents = list(groups.values)
    # A replicate label for the pseudobulk path. Three pseudo-replicates per
    # group: enough for DESeq2 to estimate dispersion, and deterministic.
    sub.meta_data["rep"] = [f"{g}_r{i % 3}" for i, g in enumerate(groups.values)]

    out = {}
    for test in (tests or TEST_MAP):
        kwargs = {"sample_col": "rep"} if test == "deseq2" else {}
        out[test] = find_markers(sub, IDENT_1, IDENT_2, test_use=test,
                                 logfc_threshold=0, min_pct=0, **kwargs)
    return out


def compare(py: pd.DataFrame, r: pd.DataFrame, test: str) -> dict:
    """One test's agreement with Seurat, on the genes both scored."""
    from scipy.stats import spearmanr

    shared = py.index.intersection(r.index)
    res: dict = {"n_python": len(py), "n_r": len(r), "n_shared": len(shared)}

    if "avg_log2FC" in py and "avg_log2FC" in r:
        d = np.abs(py.loc[shared, "avg_log2FC"] - r.loc[shared, "avg_log2FC"])
        res["log2fc_max_abs_diff"] = float(d.max())
        res["log2fc_exact"] = bool(d.max() <= LOG2FC_TOLERANCE)

    if test == "roc":
        d = np.abs(py.loc[shared, "myAUC"] - r.loc[shared, "myAUC"])
        res["auc_max_abs_diff"] = float(d.max())
        res["auc_within_seurat_rounding"] = bool(d.max() <= AUC_TOLERANCE)
        return res

    pp, pr = py.loc[shared, "p_val"], r.loc[shared, "p_val"]
    # R writes NaN where a test could not be run and 0 where it underflowed;
    # neither carries a rank, so both are excluded rather than imputed.
    ok = pp.notna() & pr.notna() & (pp > 0) & (pr > 0)
    res["n_compared"] = int(ok.sum())
    res["r_nan"] = int(pr.isna().sum())
    res["r_underflow_zero"] = int((pr == 0).sum())
    if ok.sum() > 2:
        res["p_spearman"] = float(spearmanr(pp[ok], pr[ok])[0])
        res["p_max_log10_ratio"] = float(np.abs(np.log10(pp[ok] / pr[ok])).max())
    # Ranked on the same `ok` genes the correlation uses. Including R's underflow
    # zeros here would rank 168 genes that are all exactly 0 against each other,
    # and the arbitrary tie-break made a perfectly-agreeing wilcox read 0/50.
    top_py = set(pp[ok].sort_values().head(TOP_N).index)
    top_r = set(pr[ok].sort_values().head(TOP_N).index)
    res[f"top{TOP_N}_overlap"] = len(top_py & top_r)

    # Where the gene is actually expressed, agreement should be strong; the tail
    # is near-empty rows where every NB/hurdle fit is ill-conditioned.
    expressed = shared[(py.loc[shared, "pct.1"] > 0.05) | (py.loc[shared, "pct.2"] > 0.05)]
    e = expressed[ok.reindex(expressed, fill_value=False)]
    if len(e) > 2:
        res["p_spearman_expressed"] = float(spearmanr(pp[e], pr[e])[0])
        res["n_expressed"] = len(e)
    return res


def run_full(data_dir=None, verbose=True):
    FIGURES.mkdir(exist_ok=True)
    obj = build(data_dir)
    groups = shared_groups(obj)
    groups.to_csv(FIGURES / "groups.csv", header=True)

    results = run_tests(obj, groups)
    for test, df in results.items():
        df.to_csv(FIGURES / f"py_{test}.csv")

    summary = {"n_cells": len(groups),
               "group_sizes": groups.value_counts().to_dict(),
               "n_genes": int(len(next(iter(results.values()))))}
    (FIGURES / "py_summary.json").write_text(json.dumps(summary, indent=2))
    if verbose:
        print(f"\n  {summary['n_cells']} cells "
              f"(cluster {IDENT_1}: {summary['group_sizes'].get(IDENT_1)}, "
              f"cluster {IDENT_2}: {summary['group_sizes'].get(IDENT_2)}), "
              f"{summary['n_genes']} genes")
        for test, df in results.items():
            print(f"    {test:9s} {len(df):6d} genes")
    return obj, results


def report_concordance():
    if not (FIGURES / "py_wilcox.csv").exists():
        raise SystemExit("Run the tutorial first — figures_de/py_*.csv are missing.")
    missing = [t for t in TEST_MAP if not (FIGURES / f"r_{t.lower()}.csv").exists()]
    if missing:
        raise SystemExit(
            f"Missing R output for: {', '.join(missing)}\n"
            "Run: Rscript tutorials/pbmc3k_de_verify.R"
        )
    rows = []
    for test in TEST_MAP:
        py = pd.read_csv(FIGURES / f"py_{test}.csv", index_col=0)
        r = pd.read_csv(FIGURES / f"r_{test.lower()}.csv", index_col=0)
        # Before anything is compared: was this R table computed on the groups
        # the Python side just wrote? Nothing else here can tell the difference
        # between a port that regressed and an R run left over from a previous
        # clustering, and the second reads exactly like the first.
        check_shared_groups(py, r, source=f"figures_de/r_{test.lower()}.csv")
        rows.append({"test": test, **compare(py, r, test)})
    table = pd.DataFrame(rows).set_index("test")
    _print_report(table)
    verdicts = check_bands(BANDS, measure_bands(table))
    table.attrs["band_verdicts"] = verdicts
    table.attrs["bands_hold"] = render_verdicts(
        verdicts, "Declared bands — every row of the parity table.")
    return table


def measure_bands(table: pd.DataFrame) -> dict[str, float]:
    """Pull the numbers :data:`BANDS` judges out of the concordance table.

    Kept separate from :func:`compare` so the bands can be checked against a
    table read from anywhere — including a deliberately damaged one in the
    tests, which is the only way to know the assertions can fail.
    """
    measured: dict[str, float] = {}
    for test in table.index:
        if f"{test} top50" in BANDS:
            measured[f"{test} top50"] = table.loc[test].get(
                f"top{TOP_N}_overlap", float("nan"))
        if f"{test} rho>5%" in BANDS:
            measured[f"{test} rho>5%"] = table.loc[test].get(
                "p_spearman_expressed", float("nan"))
    # One band over every test that shares Seurat's cell-level fold change.
    # deseq2 is excluded by name rather than by threshold: its fold change is
    # computed on summed counts, so ~3.5 is correct there and would silently
    # set the maximum for everyone else.
    cell_level = [t for t in table.index if t != "deseq2"]
    measured["max |dlog2FC| (parity tests)"] = float(
        table.loc[cell_level, "log2fc_max_abs_diff"].max())
    if "auc_max_abs_diff" in table.columns:
        measured["roc max |dAUC|"] = float(table.loc["roc", "auc_max_abs_diff"])
    return measured


def _print_report(table: pd.DataFrame) -> None:
    print(f"\n  {'test':10s} {'shared':>7s} {'max|dFC|':>10s} "
          f"{'p spearman':>11s} {'expressed':>10s} {'top' + str(TOP_N):>7s}")
    print("  " + "-" * 62)
    for test, row in table.iterrows():
        fc = row.get("log2fc_max_abs_diff", float("nan"))
        sp = row.get("p_spearman", float("nan"))
        ex = row.get("p_spearman_expressed", float("nan"))
        top = row.get(f"top{TOP_N}_overlap", float("nan"))
        if test == "roc":
            print(f"  {test:10s} {int(row['n_shared']):7d} {fc:10.2e} "
                  f"{'AUC ' + format(row['auc_max_abs_diff'], '.1e'):>11s} "
                  f"{'(Seurat rounds to 3dp)':>19s}")
            continue
        print(f"  {test:10s} {int(row['n_shared']):7d} {fc:10.2e} "
              f"{sp:11.6f} {ex:10.4f} {int(top):4d}/{TOP_N}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", default=None)
    parser.add_argument("--report", action="store_true",
                        help="compare against figures_de/r_<test>.csv")
    args = parser.parse_args()
    if args.report:
        table = report_concordance()
        # Exit non-zero so `--report` can be used as a check rather than read as
        # a table. A band that only prints is the situation this replaced.
        if not table.attrs.get("bands_hold", True):
            raise SystemExit(1)
        return
    run_full(data_dir=args.data_dir)
    print(f"\n  Wrote {FIGURES}")
    print("  Next: Rscript tutorials/pbmc3k_de_verify.R")


if __name__ == "__main__":
    main()
