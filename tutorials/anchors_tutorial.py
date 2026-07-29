#!/usr/bin/env python
"""Anchor internals — ``find_integration_anchors`` / ``integrate_data`` vs Seurat v4.

The integration tutorial (``integration_vignette.md``) compares *clusterings*:
it asks whether the two tools recover the same structure after correction. That
is the right end-to-end check, and it is also lossy — a partition can agree
while the anchors underneath it do not.

This tutorial compares the anchors themselves. For both reductions it asks:

  * which mutual-nearest-neighbour **pairs** does each tool call an anchor,
  * what **score** does each give them, and
  * what does the **corrected expression** of the query half come out as.

Run order (each step writes what the next one reads)::

    python  tutorials/anchors_tutorial.py       # writes the cell list + HVGs
    Rscript tutorials/anchors_verify.R          # writes the Seurat anchors
    python  tutorials/anchors_tutorial.py --report

Data: the ifnb counts exported once by ``tutorials/export_seuratdata.R ifnb``.
A fixed 2,400-cell subsample (1,200 CTRL / 1,200 STIM) keeps the R side to a
couple of minutes; the anchor algorithm does not care about scale, and the
full-dataset behaviour is what ``integration_vignette.md`` already covers.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from truecell.anchors import find_integration_anchors, integrate_data  # noqa: E402
from truecell.preprocessing import (  # noqa: E402
    find_variable_features,
    normalize_data,
    scale_data,
)
from truecell.truecell import create_truecell_object  # noqa: E402

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures_anchors"
DATA = Path.home() / ".truecell_data" / "ifnb"

N_PER_GROUP = 1200
N_FEATURES = 2000
DIMS = 30
SUBSAMPLE_SEED = 7

# Anchor agreement is a set comparison, so it is reported rather than asserted
# against a tolerance; these are the thresholds the vignette's claims rest on.
MIN_RECALL = {"cca": 0.95, "rpca": 0.95}


# ---------------------------------------------------------------- data
def _read_counts():
    with gzip.open(DATA / "matrix.mtx.gz") as fh:
        mat = sio.mmread(fh).tocsc()
    with gzip.open(DATA / "features.tsv.gz", "rt") as fh:
        features = [ln.split("\t")[0].strip() for ln in fh]
    with gzip.open(DATA / "barcodes.tsv.gz", "rt") as fh:
        barcodes = [ln.strip() for ln in fh]
    meta = pd.read_csv(DATA / "metadata.csv", index_col=0).loc[barcodes]
    return mat, features, barcodes, meta


def _subsample(barcodes, meta):
    """A fixed 1,200-per-condition draw, written out so R uses the same cells."""
    rng = np.random.default_rng(SUBSAMPLE_SEED)
    out = {}
    for group in ("CTRL", "STIM"):
        pool = [b for b in barcodes if meta.loc[b, "stim"] == group]
        out[group] = sorted(rng.choice(pool, size=N_PER_GROUP, replace=False).tolist())
    return out


def build_pair():
    """Two per-condition objects, normalized and scaled on shared anchor features."""
    mat, features, barcodes, meta = _read_counts()
    # CreateSeuratObject(min.cells = 3) runs on the whole object before the split,
    # so the gene filter has to be computed on the full matrix, not per batch.
    keep = np.asarray((mat > 0).sum(axis=1)).ravel() >= 3
    kept_features = [f for f, k in zip(features, keep) if k]
    sub = mat[keep, :]

    cells = _subsample(barcodes, meta)
    FIG.mkdir(exist_ok=True)
    for group, ids in cells.items():
        (FIG / f"cells_{group}.txt").write_text("\n".join(ids) + "\n")

    objects = {}
    for group, ids in cells.items():
        idx = [barcodes.index(c) for c in ids]
        obj = create_truecell_object(
            counts=sub[:, idx], assay="RNA", feature_names=kept_features,
            cell_names=ids, meta_data=meta.loc[ids],
        )
        normalize_data(obj)
        find_variable_features(obj, selection_method="vst", nfeatures=N_FEATURES)
        objects[group] = obj
    return objects


def shared_features(objects):
    """The anchor features, written out so R integrates on the same basis."""
    per = [set(o.get_assay().variable_features) for o in objects.values()]
    ranked = objects["CTRL"].get_assay().variable_features
    common = [f for f in ranked if all(f in s for s in per)]
    if len(common) < N_FEATURES:
        # Fill from the other object's ranking, the way SelectIntegrationFeatures
        # does, so both tools still integrate on N_FEATURES genes.
        for f in objects["STIM"].get_assay().variable_features:
            if f not in common:
                common.append(f)
            if len(common) == N_FEATURES:
                break
    feats = common[:N_FEATURES]
    (FIG / "anchor_features.txt").write_text("\n".join(feats) + "\n")
    return feats


# ---------------------------------------------------------------- the run
def run(reduction, objects, feats):
    pair = [objects["CTRL"], objects["STIM"]]
    anchors = find_integration_anchors(
        pair, anchor_features=feats, reduction=reduction, dims=DIMS
    )
    merged = integrate_data(anchors)
    return anchors, merged


def _dense(obj, feats):
    m = obj.get_assay().layer_data("data", features=list(feats))
    return np.asarray(m.todense()) if sp.issparse(m) else np.asarray(m)


def summarise(reduction, anchors, merged, objects, feats):
    c1 = objects["CTRL"].cell_names()
    c2 = objects["STIM"].cell_names()
    pairs = [(c1[int(i)], c2[int(j)]) for i, j in
             zip(anchors.anchors["cell1"], anchors.anchors["cell2"])]
    frame = pd.DataFrame(pairs, columns=["cell1", "cell2"])
    frame["score"] = anchors.anchors["score"].to_numpy()
    frame.to_csv(FIG / f"py_anchors_{reduction}.csv", index=False)

    raw = _dense(objects["STIM"], feats)
    got = np.asarray(merged.get_assay("integrated").layer_data("data", features=feats))
    names = merged.cell_names()
    delta = got[:, [names.index(c) for c in c2]] - raw
    stats = {
        "n_anchors": int(len(frame)),
        "score_mean": float(frame["score"].mean()),
        "delta_mean_abs": float(np.abs(delta).mean()),
        "delta_max_abs": float(np.abs(delta).max()),
        "delta_frac_nonzero": float(np.mean(delta != 0)),
    }
    (FIG / f"py_summary_{reduction}.json").write_text(json.dumps(stats, indent=2))
    return stats


# ---------------------------------------------------------------- comparison
def compare(reduction):
    """truecell's anchors against Seurat's, for one reduction."""
    r_path = FIG / f"r_anchors_{reduction}.csv"
    if not r_path.exists():
        return None
    r = pd.read_csv(r_path)
    p = pd.read_csv(FIG / f"py_anchors_{reduction}.csv")
    r_pairs = set(zip(r["cell1"], r["cell2"]))
    p_pairs = set(zip(p["cell1"], p["cell2"]))
    shared = r_pairs & p_pairs

    rs = r.set_index(["cell1", "cell2"])["score"]
    ps = p.set_index(["cell1", "cell2"])["score"]
    order = list(shared)
    x = np.array([ps.loc[k] for k in order]) if order else np.array([])
    y = np.array([rs.loc[k] for k in order]) if order else np.array([])

    out = {
        "reduction": reduction,
        "n_seurat": len(r_pairs),
        "n_truecell": len(p_pairs),
        "n_shared": len(shared),
        "recall": len(shared) / max(len(r_pairs), 1),
        "precision": len(shared) / max(len(p_pairs), 1),
        "score_corr": float(np.corrcoef(x, y)[0, 1]) if len(order) > 1 else float("nan"),
        "score_max_abs_diff": float(np.abs(x - y).max()) if len(order) else float("nan"),
        "score_identical": float(np.mean(np.isclose(x, y, atol=1e-9))) if len(order) else 0.0,
        "score_mean_truecell": float(x.mean()) if len(order) else float("nan"),
        "score_mean_seurat": float(y.mean()) if len(order) else float("nan"),
    }
    r_sum = FIG / f"r_summary_{reduction}.json"
    p_sum = FIG / f"py_summary_{reduction}.json"
    if r_sum.exists() and p_sum.exists():
        rj, pj = json.loads(r_sum.read_text()), json.loads(p_sum.read_text())
        out["delta_mean_abs"] = (pj["delta_mean_abs"], rj["delta_mean_abs"])
        out["delta_frac_nonzero"] = (pj["delta_frac_nonzero"], rj["delta_frac_nonzero"])
    return out


def report():
    rows = [c for c in (compare("cca"), compare("rpca")) if c]
    if not rows:
        print("No Seurat anchors found. Run `Rscript tutorials/anchors_verify.R` first.")
        return 1
    print("=" * 78)
    print("Anchor agreement — truecell vs Seurat 5.5.1")
    print("=" * 78)
    print(f"{'reduction':<10}{'Seurat':>9}{'truecell':>9}{'shared':>9}"
          f"{'recall':>9}{'prec':>8}{'score r':>10}{'same':>8}")
    ok = True
    for c in rows:
        print(f"{c['reduction']:<10}{c['n_seurat']:>9}{c['n_truecell']:>9}{c['n_shared']:>9}"
              f"{100*c['recall']:>8.1f}%{100*c['precision']:>7.1f}%"
              f"{c['score_corr']:>10.5f}{100*c['score_identical']:>7.1f}%")
        ok &= c["recall"] >= MIN_RECALL[c["reduction"]]
    print()
    for c in rows:
        if "delta_mean_abs" in c:
            py, rr = c["delta_mean_abs"]
            pyf, rrf = c["delta_frac_nonzero"]
            print(f"  {c['reduction']}: correction over the query half — "
                  f"mean|d| {py:.6f} vs Seurat {rr:.6f}, "
                  f"frac nonzero {pyf:.4f} vs {rrf:.4f}")
    print("\n  recall    = fraction of Seurat's anchors truecell also found")
    print("  precision = fraction of truecell's anchors Seurat also found")
    print("  score r   = correlation of the anchor scores on the shared pairs")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", action="store_true",
                    help="compare against the Seurat anchors and exit")
    args = ap.parse_args()
    if args.report:
        return report()

    print("Loading ifnb and building the CTRL/STIM pair ...")
    objects = build_pair()
    feats = shared_features(objects)
    for obj in objects.values():
        scale_data(obj, features=feats)
    print(f"  {N_PER_GROUP} CTRL / {N_PER_GROUP} STIM cells, "
          f"{len(feats)} anchor features -> {FIG.name}/")

    for reduction in ("cca", "rpca"):
        anchors, merged = run(reduction, objects, feats)
        s = summarise(reduction, anchors, merged, objects, feats)
        print(f"  {reduction:<5} {s['n_anchors']:>5} anchors  "
              f"score mean {s['score_mean']:.4f}  "
              f"correction mean|d| {s['delta_mean_abs']:.6f}")
    print("\nNow run:  Rscript tutorials/anchors_verify.R")
    print("then:     python tutorials/anchors_tutorial.py --report")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
