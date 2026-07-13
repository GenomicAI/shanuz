"""Differential expression / marker gene detection.

Mirrors Seurat's FindMarkers() and FindAllMarkers().
"""
from __future__ import annotations

import warnings
from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu


def _roc_auc(x1: np.ndarray, x2: np.ndarray) -> tuple[float, float]:
    """Return (AUC, power) for classifying group 1 vs group 2 by expression.

    AUC is the Mann-Whitney U statistic normalised to [0, 1]; power = |2·AUC−1|
    (Seurat's classification power), 0 = random, 1 = perfect separation.
    """
    n1, n2 = len(x1), len(x2)
    if n1 == 0 or n2 == 0:
        return 0.5, 0.0
    try:
        u, _ = mannwhitneyu(x1, x2, alternative="two-sided", method="asymptotic")
        auc = u / (n1 * n2)
    except ValueError:
        auc = 0.5
    return float(auc), float(abs(2.0 * auc - 1.0))


def _lr_pvalue(expr: np.ndarray, group: np.ndarray, latent: Optional[np.ndarray]) -> float:
    """Logistic-regression likelihood-ratio test (Seurat's 'LR').

    Fits group ~ expr (+ latent) vs the reduced group ~ (latent); the LRT on the
    dropped expression term is χ²(df=1).
    """
    import statsmodels.api as sm
    from scipy.stats import chi2

    n = len(group)
    full_cols = [np.ones(n), expr]
    red_cols = [np.ones(n)]
    if latent is not None and latent.size:
        full_cols.append(latent)
        red_cols.append(latent)
    X_full = np.column_stack(full_cols)
    X_red = np.column_stack(red_cols)
    try:
        with warnings.catch_warnings():
            # Marker genes often (near-)perfectly separate the groups; that is
            # the signal, not an error — silence statsmodels' separation noise.
            warnings.simplefilter("ignore")
            full = sm.GLM(group, X_full, family=sm.families.Binomial()).fit()
            red = sm.GLM(group, X_red, family=sm.families.Binomial()).fit()
        stat = red.deviance - full.deviance
        return float(chi2.sf(max(stat, 0.0), df=1))
    except Exception:
        return 1.0


def _negbinom_pvalue(counts: np.ndarray, group: np.ndarray, latent: Optional[np.ndarray]) -> float:
    """Negative-binomial GLM likelihood-ratio test on counts (Seurat's 'negbinom').

    Fits counts ~ group (+ latent) vs counts ~ (latent) with a fixed
    moment-estimated dispersion, LRT on the group term is χ²(df=1).
    """
    import statsmodels.api as sm
    from scipy.stats import chi2

    y = counts.astype(float)
    m = y.mean()
    if m <= 0:
        return 1.0
    v = y.var()
    alpha = max((v - m) / (m * m), 1e-6) if v > m else 1e-6

    n = len(y)
    full_cols = [np.ones(n), group]
    red_cols = [np.ones(n)]
    if latent is not None and latent.size:
        full_cols.append(latent)
        red_cols.append(latent)
    X_full = np.column_stack(full_cols)
    X_red = np.column_stack(red_cols)
    try:
        fam = sm.families.NegativeBinomial(alpha=alpha)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            full = sm.GLM(y, X_full, family=fam).fit()
            red = sm.GLM(y, X_red, family=fam).fit()
        stat = red.deviance - full.deviance
        return float(chi2.sf(max(stat, 0.0), df=1))
    except Exception:
        return 1.0


def _mast_pvalue(expr: np.ndarray, group: np.ndarray, latent: Optional[np.ndarray]) -> float:
    """MAST two-part hurdle likelihood-ratio test (Finak 2015; Seurat's 'MAST').

    Fits a *discrete* logistic model of detection (``expr > 0``) and a
    *continuous* Gaussian model of the log-expression among detected cells, each
    as ``~ group (+ latent)``. Because the hurdle likelihood factorises into its
    detection and magnitude parts, the combined LR statistic is the sum of the
    two components' statistics tested on the sum of their degrees of freedom.
    Components that carry no information (constant detection, or magnitude seen in
    only one group) contribute 0 df and are dropped.
    """
    import statsmodels.api as sm
    from scipy.stats import chi2

    n = len(expr)
    detect = (expr > 0).astype(float)

    def _design(mask: np.ndarray, include_group: bool) -> np.ndarray:
        cols = [np.ones(int(mask.sum()))]
        if include_group:
            cols.append(group[mask])
        if latent is not None and latent.size:
            cols.append(latent[mask])
        return np.column_stack(cols)

    stat = 0.0
    df = 0

    # ---- discrete component: logistic LRT on the group term ------------------
    if 0.0 < detect.sum() < n:  # detection varies → group term is estimable
        allmask = np.ones(n, dtype=bool)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                full = sm.GLM(detect, _design(allmask, True),
                              family=sm.families.Binomial()).fit()
                red = sm.GLM(detect, _design(allmask, False),
                             family=sm.families.Binomial()).fit()
            d = red.deviance - full.deviance
            if np.isfinite(d) and d > 0:
                stat += d
                df += 1
        except Exception:
            pass

    # ---- continuous component: Gaussian LRT among detected cells -------------
    pos = expr > 0
    if pos.sum() >= 3 and np.unique(group[pos]).size > 1 and np.ptp(expr[pos]) > 0:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                full = sm.OLS(expr[pos], _design(pos, True)).fit()
                red = sm.OLS(expr[pos], _design(pos, False)).fit()
            d = 2.0 * (full.llf - red.llf)
            if np.isfinite(d) and d > 0:
                stat += d
                df += 1
        except Exception:
            pass

    if df == 0:
        return 1.0
    return float(chi2.sf(stat, df=df))


def find_markers(
    seurat,
    ident_1: Union[str, list[str]],
    ident_2: Optional[Union[str, list[str]]] = None,
    assay: Optional[str] = None,
    layer: Optional[str] = None,
    test_use: str = "wilcox",
    only_pos: bool = False,
    min_pct: float = 0.1,
    logfc_threshold: float = 0.25,
    features: Optional[list[str]] = None,
    latent_vars: Optional[list[str]] = None,
    sample_col: Optional[str] = None,
    max_cells_per_ident: Optional[int] = None,
    random_seed: int = 1,
) -> pd.DataFrame:
    """Find differentially expressed marker genes.

    Mirrors R's FindMarkers(pbmc, ident.1 = 2).

    Parameters
    ----------
    ident_1         : cluster label(s) for group 1
    ident_2         : cluster label(s) for group 2 (None = all others)
    test_use        : statistical test — 'wilcox' (default), 't', 'LR'
                      (logistic-regression LRT), 'negbinom' (negative-binomial
                      GLM LRT on counts), 'mast' (MAST two-part hurdle LRT on
                      log-normalized data), 'deseq2' (pseudobulk DESeq2 — sums
                      counts per sample then tests sample-level, requires
                      ``sample_col``; needs ``pip install shanuz[deseq2]``), or
                      'roc' (AUC classifier power).
    only_pos        : only return positive markers
    min_pct         : minimum fraction cells expressing gene in either group
    logfc_threshold : minimum log2 fold-change filter
    features        : restrict to these genes (default: all)
    latent_vars     : metadata columns to regress out as covariates in the
                      'LR', 'negbinom', and 'mast' models (Seurat's latent.vars).
                      For 'mast', pass the cellular detection rate here to match
                      Seurat's default CDR covariate.
    sample_col      : metadata column identifying pseudobulk replicates (donor /
                      sample); required for ``test_use='deseq2'``, ignored
                      otherwise.
    max_cells_per_ident : downsample each group to this many cells

    Returns
    -------
    For 'wilcox' / 't' / 'LR' / 'negbinom': DataFrame with columns
    p_val, avg_log2FC, pct.1, pct.2, p_val_adj (sorted by p_val).
    For 'roc': columns myAUC, avg_diff, power, avg_log2FC, pct.1, pct.2
    (sorted by power), with no p-value — matching Seurat.
    """
    assay_name = assay or seurat.active_assay
    assay_obj = seurat.assays[assay_name]
    cells = seurat.cell_names()
    idents = list(seurat.idents)

    # Resolve ident strings
    ident_1_set = {str(ident_1)} if isinstance(ident_1, str) else {str(i) for i in ident_1}
    if ident_2 is None:
        ident_2_set = {str(i) for i in set(idents) if str(i) not in ident_1_set}
    else:
        ident_2_set = {str(ident_2)} if isinstance(ident_2, str) else {str(i) for i in ident_2}

    # Cell indices for each group
    cells_1 = [c for c, i in zip(cells, idents) if str(i) in ident_1_set]
    cells_2 = [c for c, i in zip(cells, idents) if str(i) in ident_2_set]

    if not cells_1:
        raise ValueError(f"No cells found with ident {ident_1}.")
    if not cells_2:
        raise ValueError(f"No cells found for comparison group.")

    # Optional downsampling
    if max_cells_per_ident is not None:
        rng = np.random.default_rng(random_seed)
        if len(cells_1) > max_cells_per_ident:
            cells_1 = list(rng.choice(cells_1, max_cells_per_ident, replace=False))
        if len(cells_2) > max_cells_per_ident:
            cells_2 = list(rng.choice(cells_2, max_cells_per_ident, replace=False))

    # Get expression matrix for all genes (features × cells)
    data, feature_names = _get_expression_matrix(assay_obj, layer)

    cell_idx_map = {c: i for i, c in enumerate(cells)}
    idx_1 = [cell_idx_map[c] for c in cells_1]
    idx_2 = [cell_idx_map[c] for c in cells_2]

    if sp.issparse(data):
        mat1 = data[:, idx_1].toarray().astype(float)  # (features × n1)
        mat2 = data[:, idx_2].toarray().astype(float)  # (features × n2)
    else:
        mat1 = np.asarray(data)[:, idx_1].astype(float)
        mat2 = np.asarray(data)[:, idx_2].astype(float)

    # Restrict features
    if features is not None:
        feat_set = set(features)
        feat_mask = np.array([f in feat_set for f in feature_names])
    else:
        feat_mask = np.ones(len(feature_names), dtype=bool)

    # Percent cells expressing (> 0)
    pct1 = (mat1 > 0).mean(axis=1)
    pct2 = (mat2 > 0).mean(axis=1)

    # Pre-filter: gene must be expressed in at least min_pct of either group
    pct_mask = (pct1 >= min_pct) | (pct2 >= min_pct)
    combined_mask = feat_mask & pct_mask

    # Log2 fold change. Data is log1p-normalized, so to match Seurat's FoldChange()
    # we un-log each cell (expm1), average within the group, then re-log:
    #   avg_log2FC = log2(mean(expm1(x1)) + 1) - log2(mean(expm1(x2)) + 1)
    # NOTE: the mean must be taken AFTER expm1, not before — expm1(mean(x)) is the
    # geometric-style mean and systematically compresses fold-changes (Jensen).
    group1_mean = np.expm1(mat1).mean(axis=1)
    group2_mean = np.expm1(mat2).mean(axis=1)
    avg_log2fc = np.log2(group1_mean + 1) - np.log2(group2_mean + 1)

    # Pre-filter by logfc_threshold
    if logfc_threshold > 0:
        fc_mask_arr = np.abs(avg_log2fc) >= logfc_threshold
        combined_mask = combined_mask & fc_mask_arr

    test_indices = np.where(combined_mask)[0]

    # Per-cell covariates for the regression-based tests (LR / negbinom).
    latent = None
    if latent_vars and test_use in ("LR", "negbinom", "mast"):
        lat1 = seurat.meta_data.loc[cells_1, latent_vars].to_numpy(dtype=float)
        lat2 = seurat.meta_data.loc[cells_2, latent_vars].to_numpy(dtype=float)
        latent = np.vstack([lat1, lat2])

    # ---- ROC test: returns AUC / power, no p-value (matches Seurat) ----------
    if test_use == "roc":
        if len(test_indices) == 0:
            return pd.DataFrame(
                columns=["myAUC", "avg_diff", "power", "avg_log2FC", "pct.1", "pct.2"]
            )
        aucs = np.empty(len(test_indices))
        powers = np.empty(len(test_indices))
        avg_diff = np.empty(len(test_indices))
        for i, fi in enumerate(test_indices):
            auc, power = _roc_auc(mat1[fi, :], mat2[fi, :])
            aucs[i] = auc
            powers[i] = power
            avg_diff[i] = mat1[fi, :].mean() - mat2[fi, :].mean()
        roc_res = pd.DataFrame(
            {
                "myAUC": aucs,
                "avg_diff": avg_diff,
                "power": powers,
                "avg_log2FC": avg_log2fc[test_indices],
                "pct.1": pct1[test_indices],
                "pct.2": pct2[test_indices],
            },
            index=[feature_names[i] for i in test_indices],
        )
        if only_pos:
            roc_res = roc_res[roc_res["avg_log2FC"] > 0]
        return roc_res.sort_values("power", ascending=False)

    if len(test_indices) == 0:
        return pd.DataFrame(
            columns=["p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]
        )

    # ---- pseudobulk DESeq2: sample-level test, not per-cell -------------------
    if test_use == "deseq2":
        return _deseq2_pseudobulk(
            seurat, assay_obj, cells_1, cells_2, sample_col,
            feature_names, test_indices, pct1, pct2, only_pos,
        )

    # ---- p-value-based tests -------------------------------------------------
    p_vals = np.ones(len(test_indices))

    if test_use == "wilcox":
        for i, fi in enumerate(test_indices):
            x1 = mat1[fi, :]
            x2 = mat2[fi, :]
            if x1.sum() == 0 and x2.sum() == 0:
                p_vals[i] = 1.0
            else:
                # mannwhitneyu (asymptotic) applies the tie correction and
                # continuity correction that base-R wilcox.test / presto use —
                # essential for scRNA data, which is dominated by zero ties.
                # scipy.stats.ranksums does NOT correct for ties.
                try:
                    _, p = mannwhitneyu(
                        x1, x2, alternative="two-sided",
                        use_continuity=True, method="asymptotic",
                    )
                except ValueError:
                    # Raised only when every value in both groups is identical.
                    p = 1.0
                p_vals[i] = p if not np.isnan(p) else 1.0
    elif test_use == "t":
        from scipy.stats import ttest_ind
        for i, fi in enumerate(test_indices):
            x1 = mat1[fi, :]
            x2 = mat2[fi, :]
            _, p = ttest_ind(x1, x2, equal_var=False)
            p_vals[i] = p if not np.isnan(p) else 1.0
    elif test_use == "LR":
        n1, n2 = mat1.shape[1], mat2.shape[1]
        group = np.concatenate([np.ones(n1), np.zeros(n2)])
        for i, fi in enumerate(test_indices):
            expr = np.concatenate([mat1[fi, :], mat2[fi, :]])
            p_vals[i] = _lr_pvalue(expr, group, latent)
    elif test_use == "mast":
        n1, n2 = mat1.shape[1], mat2.shape[1]
        group = np.concatenate([np.ones(n1), np.zeros(n2)])
        for i, fi in enumerate(test_indices):
            expr = np.concatenate([mat1[fi, :], mat2[fi, :]])
            p_vals[i] = _mast_pvalue(expr, group, latent)
    elif test_use == "negbinom":
        counts_mat, _ = _get_expression_matrix(assay_obj, "counts")
        if sp.issparse(counts_mat):
            c1 = counts_mat[:, idx_1].toarray()
            c2 = counts_mat[:, idx_2].toarray()
        else:
            c1 = np.asarray(counts_mat)[:, idx_1]
            c2 = np.asarray(counts_mat)[:, idx_2]
        n1, n2 = c1.shape[1], c2.shape[1]
        group = np.concatenate([np.ones(n1), np.zeros(n2)])
        for i, fi in enumerate(test_indices):
            cnts = np.concatenate([c1[fi, :], c2[fi, :]])
            p_vals[i] = _negbinom_pvalue(cnts, group, latent)
    else:
        raise ValueError(
            f"Unsupported test_use: {test_use!r}. "
            "Use 'wilcox', 't', 'LR', 'negbinom', 'mast', 'deseq2', or 'roc'."
        )

    # Bonferroni correction (Seurat default: multiply by total gene count)
    n_total = len(feature_names)
    p_val_adj = np.minimum(p_vals * n_total, 1.0)

    results = pd.DataFrame(
        {
            "p_val": p_vals,
            "avg_log2FC": avg_log2fc[test_indices],
            "pct.1": pct1[test_indices],
            "pct.2": pct2[test_indices],
            "p_val_adj": p_val_adj,
        },
        index=[feature_names[i] for i in test_indices],
    )

    if only_pos:
        results = results[results["avg_log2FC"] > 0]

    return results.sort_values("p_val")


def find_all_markers(
    seurat,
    assay: Optional[str] = None,
    layer: Optional[str] = None,
    test_use: str = "wilcox",
    only_pos: bool = False,
    min_pct: float = 0.1,
    logfc_threshold: float = 0.25,
    sample_col: Optional[str] = None,
    max_cells_per_ident: Optional[int] = None,
    random_seed: int = 1,
) -> pd.DataFrame:
    """Find marker genes for each cluster vs all others.

    Mirrors R's FindAllMarkers(pbmc, only.pos = TRUE).

    Returns a single DataFrame with an extra 'cluster' column.
    """
    clusters = sorted(set(str(i) for i in seurat.idents))
    all_results = []

    for cluster in clusters:
        try:
            df = find_markers(
                seurat,
                ident_1=cluster,
                ident_2=None,
                assay=assay,
                layer=layer,
                test_use=test_use,
                only_pos=only_pos,
                min_pct=min_pct,
                logfc_threshold=logfc_threshold,
                sample_col=sample_col,
                max_cells_per_ident=max_cells_per_ident,
                random_seed=random_seed,
            )
            if len(df) > 0:
                df = df.copy()
                df["cluster"] = cluster
                df["gene"] = df.index
                all_results.append(df)
        except Exception as e:
            print(f"Warning: cluster {cluster} marker finding failed: {e}")

    if not all_results:
        return pd.DataFrame(
            columns=["p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene"]
        )

    combined = pd.concat(all_results, axis=0)
    combined = combined[["cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]]
    return combined.sort_values(["cluster", "p_val"])


def find_conserved_markers(
    seurat,
    ident_1: Union[str, list[str]],
    grouping_var: str,
    ident_2: Optional[Union[str, list[str]]] = None,
    assay: Optional[str] = None,
    layer: Optional[str] = None,
    test_use: str = "wilcox",
    only_pos: bool = False,
    min_pct: float = 0.1,
    logfc_threshold: float = 0.25,
    features: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Find markers conserved across the levels of a grouping variable.

    Mirrors R's ``FindConservedMarkers(obj, ident.1, grouping.var = "stim")``:
    runs :func:`find_markers` for ``ident_1`` vs ``ident_2`` independently within
    each level of ``grouping_var``, keeps only genes detected as markers in
    *every* level, and combines their per-level p-values with Fisher's method
    (:func:`scipy.stats.combine_pvalues`).

    Parameters
    ----------
    ident_1      : cluster label(s) for group 1.
    grouping_var : metadata column whose levels define the independent
                   comparisons (e.g. condition, batch, donor).
    ident_2      : cluster label(s) for group 2 (None = all other cells).
    (remaining args are forwarded to :func:`find_markers`.)

    Returns
    -------
    DataFrame indexed by gene with, for each level ``g``, the columns
    ``{g}_p_val, {g}_avg_log2FC, {g}_pct.1, {g}_pct.2, {g}_p_val_adj`` plus
    ``max_pval`` (worst per-level p-value) and ``combined_p_val`` (Fisher-combined
    across levels), sorted by ``combined_p_val``. Only genes that are markers in
    all levels are returned.
    """
    from scipy.stats import combine_pvalues

    if grouping_var not in seurat.meta_data.columns:
        raise KeyError(
            f"grouping_var {grouping_var!r} not found in meta_data "
            f"(columns: {list(seurat.meta_data.columns)})."
        )

    cells = seurat.cell_names()
    group_of = seurat.meta_data.loc[cells, grouping_var].astype(str)
    levels = sorted(group_of.unique())

    per_level: dict[str, pd.DataFrame] = {}
    for level in levels:
        level_cells = [c for c, g in zip(cells, group_of) if g == level]
        sub = seurat.subset(cells=level_cells)
        try:
            df = find_markers(
                sub,
                ident_1=ident_1,
                ident_2=ident_2,
                assay=assay,
                layer=layer,
                test_use=test_use,
                only_pos=only_pos,
                min_pct=min_pct,
                logfc_threshold=logfc_threshold,
                features=features,
            )
        except ValueError as e:
            warnings.warn(
                f"Skipping {grouping_var}={level!r}: {e}", RuntimeWarning, stacklevel=2
            )
            continue
        if len(df) > 0:
            per_level[level] = df

    if not per_level:
        raise ValueError(
            f"No level of {grouping_var!r} yielded markers for the requested comparison."
        )

    # Genes must be markers in every retained level.
    common = set.intersection(*(set(df.index) for df in per_level.values()))

    used = list(per_level)
    cols: dict[str, pd.Series] = {}
    for level in used:
        df = per_level[level].loc[list(common)]
        for c in df.columns:
            cols[f"{level}_{c}"] = df[c]
    result = pd.DataFrame(cols, index=list(common))

    if test_use == "roc":
        # ROC has no p-value; conservation is summarised by the min power.
        power_cols = [f"{level}_power" for level in used]
        result["min_power"] = result[power_cols].min(axis=1)
        return result.sort_values("min_power", ascending=False)

    pval_cols = [f"{level}_p_val" for level in used]
    result["max_pval"] = result[pval_cols].max(axis=1)
    if len(used) == 1:
        result["combined_p_val"] = result[pval_cols[0]]
    else:
        result["combined_p_val"] = [
            combine_pvalues(result.loc[g, pval_cols].to_numpy(dtype=float),
                            method="fisher").pvalue
            for g in result.index
        ]
    return result.sort_values("combined_p_val")


def _deseq2_pseudobulk(
    seurat,
    assay_obj,
    cells_1: list[str],
    cells_2: list[str],
    sample_col: Optional[str],
    feature_names: list[str],
    test_indices: np.ndarray,
    pct1: np.ndarray,
    pct2: np.ndarray,
    only_pos: bool,
) -> pd.DataFrame:
    """Pseudobulk DESeq2 test (Seurat's ``test.use = "DESeq2"``).

    Sums raw counts to one pseudobulk profile per (group, ``sample_col``) — the
    same aggregation as :func:`shanuz.aggregate.aggregate_expression` — then fits
    a DESeq2 model with design ``~condition`` and contrasts group 1 vs group 2.
    A positive ``avg_log2FC`` (DESeq2's ``log2FoldChange``) means up in group 1.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as e:  # pragma: no cover - exercised only without the dep
        raise ImportError(
            "test_use='deseq2' requires pydeseq2. Install with "
            "`pip install shanuz[deseq2]`."
        ) from e

    if sample_col is None:
        raise ValueError(
            "test_use='deseq2' is a pseudobulk test and requires sample_col — the "
            "replicate/donor column to aggregate cells into per-sample profiles."
        )
    if sample_col not in seurat.meta_data.columns:
        raise KeyError(f"sample_col {sample_col!r} not found in meta_data.")

    counts_mat, _ = _get_expression_matrix(assay_obj, "counts")
    counts_sub = counts_mat[test_indices, :]  # tested genes × all cells
    cell_pos = {c: i for i, c in enumerate(seurat.cell_names())}
    samples = seurat.meta_data[sample_col].astype(str)
    gene_names = [feature_names[i] for i in test_indices]

    def _pseudobulk(group_cells: list[str], cond: str) -> dict[str, np.ndarray]:
        by_sample: dict[str, list[int]] = {}
        for c in group_cells:
            by_sample.setdefault(samples[c], []).append(cell_pos[c])
        cols: dict[str, np.ndarray] = {}
        for samp, idxs in by_sample.items():
            summed = counts_sub[:, idxs].sum(axis=1)
            cols[f"{cond}::{samp}"] = np.asarray(summed).ravel()
        return cols

    pb = {**_pseudobulk(cells_1, "group1"), **_pseudobulk(cells_2, "group2")}
    sample_names = list(pb)
    condition = ["group1" if s.startswith("group1::") else "group2" for s in sample_names]

    n1, n2 = condition.count("group1"), condition.count("group2")
    if n1 < 2 or n2 < 2:
        warnings.warn(
            f"DESeq2 pseudobulk has {n1} vs {n2} replicate(s) in {sample_col!r}; "
            "dispersion estimates are unreliable without ≥2 replicates per group.",
            RuntimeWarning,
            stacklevel=2,
        )

    # pydeseq2 wants samples × genes, integer counts.
    counts_df = pd.DataFrame(
        np.column_stack([pb[s] for s in sample_names]).T.astype(int),
        index=sample_names,
        columns=gene_names,
    )
    metadata = pd.DataFrame({"condition": condition}, index=sample_names)

    dds = DeseqDataSet(counts=counts_df, metadata=metadata, design="~condition", quiet=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["condition", "group1", "group2"], quiet=True)
    stat.summary()
    res = stat.results_df.reindex(gene_names)

    out = pd.DataFrame(
        {
            "p_val": res["pvalue"].fillna(1.0).to_numpy(),
            "avg_log2FC": res["log2FoldChange"].fillna(0.0).to_numpy(),
            "pct.1": pct1[test_indices],
            "pct.2": pct2[test_indices],
            "p_val_adj": res["padj"].fillna(1.0).to_numpy(),
        },
        index=gene_names,
    )
    if only_pos:
        out = out[out["avg_log2FC"] > 0]
    return out.sort_values("p_val")


# ------------------------------------------------------------------
# Helper
# ------------------------------------------------------------------

def _get_expression_matrix(assay_obj, layer: Optional[str]):
    """Return (matrix features×cells, feature_names) using best available layer."""
    from .assay5 import Assay5
    from .assay import Assay

    if isinstance(assay_obj, Assay5):
        feature_names = assay_obj._all_feature_names
        if layer is not None and layer in assay_obj.layers:
            return assay_obj.layers[layer], feature_names
        for candidate in ("data", "counts"):
            if candidate in assay_obj.layers:
                return assay_obj.layers[candidate], feature_names
        raise ValueError("No expression data layer found in Assay5.")
    else:
        feature_names = assay_obj._feature_names
        from ._sparse import is_matrix_empty
        if layer == "counts":
            return assay_obj.counts, feature_names
        if layer == "scale_data" or layer == "scale.data":
            return assay_obj.scale_data, feature_names
        # Prefer log-normalized data
        if not is_matrix_empty(assay_obj.data):
            return assay_obj.data, feature_names
        return assay_obj.counts, feature_names
