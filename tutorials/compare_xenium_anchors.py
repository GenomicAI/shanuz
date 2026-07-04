#!/usr/bin/env python
"""Compare shanuz (Python) vs R Seurat deterministic anchors for the Xenium
spatial tutorial, and print a MATCH/DIFF table with an overall verdict.

Reads both JSONs from tutorials/figures_spatial/:
  * anchors.json      — written by generate_spatial_plots.py
  * r_reference.json  — written by xenium_spatial_verify.R
"""
import json
from pathlib import Path

FIG = Path(__file__).parent / "figures_spatial"
py = json.loads((FIG / "anchors.json").read_text())
r = json.loads((FIG / "r_reference.json").read_text())


def rel(a, b):
    return abs(a - b) / max(1.0, abs(b))


rows, ok_all = [], True


def check(name, a, b, tol=1e-6, exact=False):
    global ok_all
    ok = (a == b) if exact else (rel(float(a), float(b)) <= tol)
    ok_all = ok_all and ok
    rows.append((name, a, b, "MATCH" if ok else "DIFF"))


check("n_cells_raw", py["n_cells_raw"], r["n_cells_raw"], exact=True)
check("n_genes", py["n_genes"], r["n_genes"], exact=True)
check("n_cells_qc", py["n_cells_qc"], r["n_cells_qc"], exact=True)
for ct in sorted(set(py["celltype_counts"]) | set(r["celltype_counts"])):
    check(f"celltype[{ct}]", py["celltype_counts"].get(ct),
          r["celltype_counts"].get(ct), exact=True)
check("n_focal", py["n_focal"], r["n_focal"], exact=True)
check("focal_nn_median", py["focal_nn_median"], r["focal_nn_median"])
check("focal_nn_mean", py["focal_nn_mean"], r["focal_nn_mean"])
check("focal_local_density_mean", py["focal_local_density_mean"],
      r["focal_local_density_mean"])
check("region_ymed", py["region_ymed"], r["region_ymed"])
check("composition_chisq_p", py["composition_chisq_p"],
      r["composition_chisq_p"], tol=1e-4)
for g in sorted(py["composition"]):
    check(f"comp[{g}].log2_ratio", py["composition"][g]["log2_ratio"],
          r["composition"][g]["log2_ratio"])
    check(f"comp[{g}].p", py["composition"][g]["p"],
          r["composition"][g]["p"], tol=1e-4)
    check(f"comp[{g}].padj", py["composition"][g]["padj"],
          r["composition"][g]["padj"], tol=1e-4)

w = max(len(n) for n, *_ in rows)
print(f"{'anchor':<{w}}  {'shanuz':>16}  {'R Seurat':>16}  verdict")
print("-" * (w + 44))
for n, a, b, v in rows:
    fa = f"{a:.8g}" if isinstance(a, float) else str(a)
    fb = f"{b:.8g}" if isinstance(b, float) else str(b)
    print(f"{n:<{w}}  {fa:>16}  {fb:>16}  {v}")

print("\nOdds-ratio (definitional difference: scipy sample OR vs R conditional MLE):")
for g in sorted(py["composition"]):
    print(f"  {g:<16} shanuz={py['composition'][g]['odds_ratio']:.4f}  "
          f"R={r['composition'][g]['odds_ratio']:.4f}")

print("\nStructural (stochastic, reported not matched): "
      f"clusters shanuz={py['n_clusters']} R={r['n_clusters']}; "
      f"niches shanuz={py['n_niches']} R={r['n_niches']}")
print("\n" + ("ALL DETERMINISTIC ANCHORS MATCH" if ok_all else "SOME ANCHORS DIFFER"))
