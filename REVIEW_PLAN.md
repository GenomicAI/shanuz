# Package work from the domain-expert review

An expert reviewer read the Truecell package and manuscript and returned eleven
comments. The manuscript items are parked (`truecell_paper/community_paper/REVIEW_ROUND2_PLAN.md`
holds those). This document is the subset that lands in **this repository**: the
API gaps and validation gaps the review exposed.

Every finding below was checked against the source and against **live R Seurat
5.5.1**, not inferred from the review text.

The review is a useful instrument here precisely because it is *not* a software
review. It is a working biologist describing what she does at the keyboard, and
two of her routine habits turn out to be things this package cannot do.

---

## Verified gaps

### P1 — `find_clusters` cannot do what she does every day  ✅ DONE

> *"I typically run 3 resolutions at the same time and compare them before
> deciding which one to use for downstream analysis."* (comment 4)

This is not a preference, it is the standard Seurat idiom, and `FindClusters`
supports it directly. Truecell does not.

Verified against Seurat 5.5.1 — `FindClusters(obj, resolution = c(0.4, 0.8, 1.2))`:

```
COLUMNS: "RNA_snn_res.0.4" "RNA_snn_res.0.8" "RNA_snn_res.1.2"
Idents            == RNA_snn_res.1.2   (the LAST resolution)
seurat_clusters   == RNA_snn_res.1.2   (the LAST resolution)
```

Truecell's `find_clusters` (`truecell/clustering.py`) has **three** gaps against
this, in increasing order of severity:

| # | Gap | Consequence |
|---|---|---|
| a | `resolution: float` — scalar only | The three-resolution idiom cannot be expressed; the user writes a loop and hand-rolls the column names |
| b | No `cluster_name` parameter | Seurat's `cluster.name` override is unavailable |
| c | **Never writes a `{graph}_res.{r}` column at all** | Even at one resolution. A ported script that reads `obj[["RNA_snn_res.0.5"]]` — which is what every Seurat script does — **fails** |

(c) is the real defect and it is independent of the review: it is a plain
API-fidelity break that exists today at the default settings, and `grep` finds no
`_res.` anywhere in `clustering.py`. Only `seurat_clusters` is written.

**Fix.** Accept `float | Sequence[float]`; loop the clustering call; write one
`{graph_name}_res.{r}` column per resolution; set `seurat_clusters` **and**
`idents` from the **last** resolution in the sequence. Keep the scalar path
byte-identical to today's behaviour for `seurat_clusters` so nothing downstream
moves — the new column is additive.

Two details worth pinning from R rather than guessing, since this is exactly the
class of convention the project has been bitten by before:

- **Number formatting in the column name.** `0.4` → `res.0.4`, but what does
  `1.0` render as, or `0.50`? R's `paste0` on a numeric drops trailing zeros;
  Python's `str(1.0)` gives `"1.0"`. Pin it against R for a spread of values
  including integers, trailing zeros and more than one decimal place.
- **Seed behaviour across the sequence.** Does Seurat re-seed per resolution or
  let the RNG run on? This determines whether resolution *k*'s partition depends
  on how many resolutions preceded it. Test it; do not assume.

### P2 — `average_expression` does not exist  ✅ DONE

Seurat ships **two** group-summary functions and Truecell has only one:

| Seurat | Truecell | Semantics (verified) |
|---|---|---|
| `AggregateExpression` | `aggregate_expression` ✓ | Sum of raw counts per group |
| `AverageExpression` | **missing** | `rowMeans(expm1(data))` per group |

`AverageExpression` is *not* a count mean, and the difference is not subtle —
on a 30 × 12 test object gene 1 gave **332.84** from `AverageExpression` against
a count mean of **3.17**. It back-transforms the log1p-normalised `data` layer
through `expm1`, then means. Confirmed by elimination:

```
AverageExpression : 332.8377 260.9408 276.3951 196.5893
mean(expm1(data)) : 332.8377 260.9408 276.3951 196.5893   <- match
expm1(mean(data)) : 257.1996 113.7119 211.9232 169.9232
```

**This is cheap to implement correctly**, because the hard part already exists.
`markers.py::_row_expm1_sum` (added in #92) computes row-wise `sum(expm1(m))`
while preserving sparsity — `expm1(0) == 0`, so the transform touches stored
values only. `average_expression` is that helper divided by the group size.
Reuse it rather than writing a second copy; a divergence between the fold-change
path and the average-expression path would be invisible.

Scope check before building: `AverageExpression` also takes `group.by`, `assays`,
`features`, `layer`, `return.seurat`, and its behaviour on the `counts` vs `data`
layer differs. Pin each against R.

---

## Validation gaps

### T1 — `p_val_adj` is never compared. Anywhere.  ✅ DONE (with T2)

> *"maybe we should show whether identical adj p-values occur (sometimes, people
> base their manual annotation on adj p-values)"* (comment 7)

`grep` for `p_val_adj|padj|adjust` across `tutorials/pbmc3k_de_tutorial.py`
returns **nothing**. The DE evaluation — the one the reviewer called the best
section of the paper, carrying the headline "seven of eight tests reproduce the
top 50 exactly" — compares raw `p_val` and never touches the adjusted column.

That column is what a biologist actually reads. It is in every marker table the
package emits, and no test or tutorial has ever checked it against R.

**Fix.** Extend `compare()` in the DE tutorial with, per test:
- fraction of shared genes whose `p_val_adj` is bit-identical
- fraction agreeing to 1 / 3 / 6 significant figures
- **agreement on the decision at 0.05 and 0.01** — the one that determines what
  gets annotated
- count of genes where the two disagree on significance, listed rather than
  summarised

**Expect one honest complication.** Both tools Bonferroni-correct over the gene
universe, so identical adjusted p needs identical raw p *and* identical *n*.
Where the universes differ the adjusted values differ **for a reason that is not
a defect** — report the cause, not just a lower number.

**Reuse, don't rewrite:** `pbmc3k_de_tutorial.py:243–256` already handles R
writing `NaN` for un-runnable tests and `0` for underflow (168 genes on PBMC 3k),
and already records that naively including the underflow zeros made a
perfectly-agreeing Wilcoxon score 0/50. That logic applies unchanged here.

### T2 — logFC is compared by max-abs-diff, never by rank  ✅ DONE

The tutorial computes `log2fc_max_abs_diff` (a worst-case bound) and Spearman on
raw p, but no rank correlation on `avg_log2FC`. Comments 3 and 11.1.

Max-abs-diff and rank correlation fail differently: a uniform scale error leaves
ranks perfect and blows up the max, while a handful of swapped mid-table genes
leaves the max tiny and moves the ranks. Reporting both is strictly more
informative than either, and the second is what a reader ranking markers cares
about.

**Fix.** Add Spearman ρ and Kendall τ on `avg_log2FC`, plus top-*k* overlap by
logFC at k = 10/25/50 alongside the existing top-50-by-p.

### T3 — pseudobulk has no R-verified tutorial  ◐ RISK CLOSED BY TESTS

`AggregateExpression` appears in **zero** R verify scripts. `aggregate_expression`
has nine tests in `tests/test_pseudobulk_conserved.py` and every one of them
checks Python against Python — `agg["d1"] == dense[:, ::2].sum(axis=1)`.

This is the CLR risk profile exactly: a kernel that is formulaically correct,
verified against a Python-side re-derivation of its own formula, and therefore
unable to detect a convention mismatch. That is the defect this project's whole
validation strategy was designed around after §defect.

**It was a defect, not just a coverage gap.** Pinning against R found that
`return_object=True` left the raw sums in the `data` layer where Seurat writes
`log1p(sums / colSums × 10000)` — 14 against Seurat's 6.98. Fixed, with
`normalization_method` / `scale_factor` matching Seurat's arguments. The
sum-of-counts default and multi-column `group_by` were already right.

Also learned: **Seurat's `AggregateExpression` has no `layer` argument** and
always sums `counts`. Truecell's is a superset with a matching default — not a
fidelity break, now documented as deliberate.

`tests/test_pseudobulk_vs_r.py` pins six behaviours against R. The full tutorial
below is still worth building for the DE half, but the *silent-defect* risk is
closed.

**Fix.** A paired pseudobulk tutorial (#19), ifnb CTRL-vs-STIM per cell type,
with declared anchors: the aggregated matrices agree to floating point, then the
DESeq2 result agrees on logFC, adjusted p and the significant-gene set.

⚠️ **State the caveat in the tutorial itself:** ifnb has no true biological
replicates. Pseudo-replicates split from one sample give a result that is
*comparable between implementations* but is **not** valid biological inference.
The tutorial measures implementation concordance on a shared computation. Anyone
who knows pseudobulk will check this, and getting it wrong costs more than the
tutorial gains.

### T4 — `add_module_score` is R-verified only through cell-cycle  ◐ PART 1 WAS ALREADY DONE

> *"I breathe and eat gene module scores in my analyses."* (comment 11.3)

`AddModuleScore` appears in exactly one R script, `thp1_cellcycle_verify.R`, via
`cell_cycle_scoring` on the Tirosh S/G2M sets. That path is well measured
(Pearson 0.9982, 96.6% phase concordance on 20,729 cells).

But it is one gene set, at default `nbin=24` / `ctrl=100`, through a wrapper.
Direct `add_module_score` on an arbitrary program — the reviewer's actual use —
is unverified against R, as are `pool`, `search`, and any non-default binning.

**Fix, two parts.**
1. **Free:** the cell-cycle vignette should say that it *is* a module-score
   result. A reader looking for `AddModuleScore` fidelity currently has no way to
   find the 0.998 that already answers them.
2. Extend the cell-cycle tutorial, or add a small paired one, scoring 2–3
   arbitrary programs (a T-cell activation signature, an interferon-response set
   on ifnb, a myeloid program) at non-default `nbin`/`ctrl`. Report per-cell
   Pearson/Spearman **and** threshold-call concordance — the derived call is what
   gets used.

   Control for the RNG honestly: `AddModuleScore` samples control genes from
   expression-matched bins and NumPy's RNG is not R's, so exact agreement is not
   available. Report correlation at default `ctrl` and at a raised `ctrl`; if the
   residual shrinks, it is sampling noise, and if it does not, that is a finding.

### T5 — no resolution-stability tutorial

Follows P1. Once `find_clusters` takes a sequence, the guided-clustering tutorial
should show the idiom and compare against R at each resolution, so the fidelity
claim covers the parameter range users actually scan rather than one point in it.

### T6 — UMAP fidelity framing lives in one vignette  ✅ DONE

> *"UMAP is not mathematically unique, so comparing coordinates may not be so
> relevant… emphasize that UMAP was evaluated as implementation fidelity only."*
> (comment 6)

`tutorials/pbmc3k_tutorial.md:597–608` has a correct and well-written note saying
the layout difference is expected because `uwot` and `umap-learn` differ in
initialisation and optimisation. Other tutorials that produce UMAPs do not carry
it, and no tutorial states the *positive* form: UMAP is validated as a functional
correspondence, and is deliberately excluded from the numeric agreement metrics
because the embedding is not identified up to anything stronger than topology.

**Fix.** Promote that note to a shared statement in `tutorials/README.md` and
`docs/fidelity.md`, and cross-reference from each vignette producing a UMAP. Cheap
and it converts an absence into a stated design decision.

---

## Not package gaps (recorded so they are not re-investigated)

| Review item | Why it needs no code |
|---|---|
| 2 — which HVGs differ | `meta_features` already exposes Seurat's `HVFInfo()` columns (`mean`, `variance.expected`, `variance.standardized`, and the `mvp.*` set under the dispersion method). Inspecting the selection boundary needs no new API |
| 5 — PCA eigenvector sign/order | Already handled correctly where it matters: `tutorials/pbmc3k_tutorial.py:574` matches components on `abs(corrcoef)` with Hungarian assignment, i.e. sign- and order-invariant. Seurat has no such helper either, so this is a docs point, not an API one |
| 11.2 — pathway enrichment | Seurat has no enrichment function; nothing to port. Belongs in a tutorial as a downstream consumer of `find_all_markers` output, if anywhere |
| 8 — software engineering | No action |

---

## Sequence

| # | Work | Depends on | Size |
|---|---|---|---|
| 1 | ~~**P1** — multi-resolution `find_clusters` + `{graph}_res.{r}` column + `cluster_name`~~ **done** | — | Medium |
| 2 | ~~**T1 + T2** — adjusted-p and logFC-rank comparison in the DE tutorial~~ **done** | — | Small |
| 3 | ~~**P2** — `average_expression`, reusing `_row_expm1_sum`~~ **done** | — | Small |
| 4 | ~~**T6** — UMAP fidelity statement promoted to shared docs~~ **done** | — | Small |
| 5 | ~~**T4.1** — label the cell-cycle result as a module-score result~~ **already true in the package** | — | — |
| 6 | **T5** — resolution-stability coverage in the guided tutorial | 1 | Small |
| 7 | **T4.2** — arbitrary-program module-score verification | Nothing | Medium |
| 8 | **T3** — paired pseudobulk tutorial (#19) | 3 (if it uses `average_expression`) | Medium |

Items 1–5 are independent of each other and of any environment work.

**Cross-cutting rule, from this project's own history:** every new comparison
gets an anchor and a declared tolerance, and every new branch gets mutation-tested
rather than trusted because the suite is green. Marker code has now hidden a real
defect behind a green suite twice, and T1 exists because a validated-looking DE
comparison silently omitted a whole column.

---

## Found while doing the work: `_get_expression_matrix` had two defects  ✅ FIXED

Not in the review. Surfaced while building `average_expression`, which needed to
read the `scale.data` layer.

1. **`layer="scale_data"` silently returned the `data` layer.** The Assay5 layer
   dict is keyed `scale.data`; only that spelling was matched, so the underscore
   form missed and fell through to the `data` fallback. Right shape, no warning,
   normalized values where scaled ones were asked for.
2. **`layer="scale.data"` returned a 10-row matrix labelled with 30 names.**
   `scale.data` holds only the scaled subset, so every row read as a different
   gene. Identical in kind to the defect fixed in `reduction.py` under #66 — the
   same one had survived in this function.

Reachable from `find_markers(layer=...)`, `aggregate_expression` and
`average_expression`. **`reduction.py` and therefore PCA were never affected** —
it has its own accessor, which is where #66 was fixed.

The full suite passed both before and after the fix, which is the point: nothing
depended on the broken behaviour, and nothing tested the correct behaviour
either. Nine mutants now cover it.

---

## Found while doing the work: the R↔Python CSV handoff was lossy on both sides

Not in the review, and it affects **every** tutorial that compares numbers
through a CSV, not just the DE one. Fixed in the DE tutorial; the others are
untouched and should be checked.

| Direction | Default behaviour | Round-trips float64? |
|---|---|---|
| R `write.csv` | 15 significant digits | **No** |
| R `sprintf("%.17g")` | not correctly rounded — emits digits denoting a *different* double | **No** |
| R `sprintf("%a")` (C99 hex) | transcribes the IEEE-754 bits | **Yes** |
| pandas `to_csv` | shortest round-trippable repr | Yes |
| pandas `read_csv` (default) | misparses ~⅓ of random doubles by an ULP | **No** |
| pandas `read_csv(float_precision="round_trip")` | correctly rounded | **Yes** |

The R formatter result is the surprising one: raising the digit count does not
help, because R's `sprintf` is not correctly rounded at high precision.
`0.1234567890123456789` (shortest exact form `…568`) comes back as `…571`, a
different double, and Python parses R's own string to that different double too.

**Consequence.** Any claim of the form "these agree exactly" made by reading two
CSVs is, before this fix, partly a measurement of the two languages' text
formatters. Tolerance-based claims are unaffected as long as the tolerance is
well above ~1e-15 relative, which most of the declared bands are.

**How to apply.** Read every CSV with `float_precision="round_trip"`, and where
bit-identity is actually the question, have R write a `%a` side table.

**Scope of the remaining exposure:** `grep -rn "read_csv" tutorials/*.py | grep -v
float_precision` finds **42 call sites across 17 files**. The tutorials making
the strongest exactness claims are the ones to check first — the object model
("91 of 91 anchors, no tolerance"), the Visium container ("24 of 24, coordinates
to max|Δx| = 0"), out-of-core ("bit-identical" on-disk vs in-memory) and spatial
statistics (Moran's I to 1.6e-14). Filed as a background task; not started.
