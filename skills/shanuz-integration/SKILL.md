---
name: shanuz-integration
description: Use when combining single-cell datasets in shanuz — batch correction with Harmony, CCA or RPCA via integrate_layers or the v4 anchor path, and reference mapping / label transfer with find_transfer_anchors, transfer_data, map_query and project_umap. Includes how to score whether the integration actually worked.
---

# Integration and reference mapping

Load the `shanuz` skill first. Two distinct jobs:

- **Integration** — several datasets that should be one, batch effect removed.
  You keep all the cells.
- **Reference mapping** — one annotated reference, one unannotated query. The
  reference is not modified; the query gets labels and coordinates.

## Integration — the v5 path (use this)

```python
import shanuz

obj = ref.merge(query)                       # or build one object holding both batches
obj.meta_data["batch"] = ...                 # the batch column

shanuz.normalize_data(obj)
shanuz.find_variable_features(obj)
shanuz.scale_data(obj)
shanuz.run_pca(obj, n_pcs=50)

shanuz.integrate_layers(obj, method="harmony", orig_reduction="pca", group_by="batch")
# → obj.reductions["harmony"]

shanuz.find_neighbors(obj, reduction="harmony", dims=range(30))
shanuz.find_clusters(obj, resolution=0.5)
shanuz.run_umap(obj, reduction="harmony", dims=range(30))
```

`method`: `"harmony"` · `"cca"` · `"rpca"`. Every method corrects
`orig_reduction` and writes a **new reduction of the same shape** — which is what
makes them interchangeable behind one call.

- `group_by=` is **required for every method**, not just Harmony.
- `new_reduction` defaults to the method name (`"harmony"`, `"cca"`, `"rpca"`).
- Downstream steps must be pointed at it: `find_neighbors(reduction="harmony")`
  and `run_umap(reduction="harmony")`. Forgetting that is the single most common
  way to "integrate" and see no change.

Harmony directly:

```python
shanuz.run_harmony(obj, group_by="batch", reduction="pca", dims=range(30),
                   reduction_name="harmony", theta=None, max_iter_harmony=10)
```

Needs `pip install shanuz[integration]` (harmonypy).

### Choosing a method

| Method | When |
|---|---|
| `harmony` | Default. Fast, scales, robust across many batches. |
| `cca` | Batches with genuinely shared structure and modest cell numbers. Aggressive — can over-correct real biology. |
| `rpca` | Large datasets, or batches with limited cell-type overlap. Conservative. |

On the ifnb IFN-β benchmark all three reach batch mixing **0.991** — above
Seurat's own 0.917. Method choice is about the data's shape, not about accuracy
here.

## Integration — the v4 anchor path

Use when you want the corrected **expression matrix** rather than a corrected
embedding (`IntegrateData` semantics), or when reproducing a v4 pipeline.

```python
anchors = shanuz.find_integration_anchors(
    [obj_a, obj_b], anchor_features=None, reduction="cca", dims=30,
    k_anchor=5, k_filter=200, k_score=30, reference=0,
)
combined = shanuz.integrate_data(anchors, new_assay="integrated", k_weight=100)
# integrate_data returns a NEW object — rebind it

# Or correct an existing embedding instead of the expression:
dr = shanuz.integrate_embeddings(anchors, reduction=obj.reductions["pca"],
                                 new_reduction="integrated_dr")
```

`reduction="cca"` or `"rpca"`. Each object must be normalized, have variable
features and be scaled before anchoring. `dims=30` here is a **count**, not a
0-based range — one of the few places in the API where that is true.

RPCA agrees with **100 % of Seurat's v4 anchors** (649/649) and 30/30 PCs on the
v5 embedding; CCA anchors agree 99.9 %.

## Reference mapping

```python
anchors = shanuz.find_transfer_anchors(
    reference=ref, query=qry, reduction="pcaproject", dims=30,
)

# Labels
pred = shanuz.transfer_data(anchors, refdata="cell_type", k_weight=50)
qry.meta_data["predicted.id"] = pred["predicted.id"]
qry.meta_data["prediction.score.max"] = pred["prediction.score.max"]

# Labels + the reference's own UMAP, in one call
pred = shanuz.map_query(anchors, refdata="cell_type",
                        reference_reduction="pca", reduction_model="umap")
# → qry.reductions["ref.umap"]

# UMAP projection alone
shanuz.project_umap(qry, ref, reduction="pca", umap_reduction="umap")
```

- `reduction="pcaproject"` (default) projects the query into the reference's PCA
  — the usual choice, and the cheap one. `"cca"` is the alternative when the two
  datasets share little.
- `refdata` is a metadata column name on the reference, or an array/Series
  aligned to its cells. Pass `refdata_features=` and an expression matrix to
  **impute** reference expression onto the query instead of labels.
- `transfer_data` returns a DataFrame — assign the columns you want yourself.
- `map_query` needs the reference to carry a fitted UMAP model
  (`run_umap` on the reference, before mapping).

Label transfer is **98.71 % per-cell concordant** with R Seurat on the panc8
cross-technology benchmark, both ~98.5 % accurate against held-out ground truth.

**Always keep the prediction score.** `prediction.score.max` is how you find the
query cells whose type is absent from the reference — they get a confident-looking
label and a low score. A mapping result reported without its score distribution
is not a finished analysis.

## Did it work? Score it

Never judge integration by looking at a UMAP.

```python
import numpy as np
from sklearn.metrics import adjusted_rand_score, silhouette_score

emb = obj.embeddings("harmony", dims=list(range(30)))

# Batch mixing: 1 - normalized batch silhouette. Higher is better mixed.
batch_sil = silhouette_score(emb, obj.meta_data["batch"])
batch_mix = 1 - (batch_sil + 1) / 2

# Biology kept: clusters should still recover known cell types
ari = adjusted_rand_score(obj.meta_data["cell_type"], obj.meta_data["seurat_clusters"])
```

Both numbers together, always. Batch mixing alone is maximised by destroying the
biology; cell-type ARI alone is maximised by not integrating at all.

## Traps

| Symptom | Cause |
|---|---|
| Integration "does nothing" | Downstream `find_neighbors` / `run_umap` still on `reduction="pca"`. |
| `TypeError` / missing batch | `group_by=` omitted — required for all three `integrate_layers` methods. |
| Anchor finding crashes or is empty | Objects not scaled, or no shared variable features. Anchor input is `layer="scale.data"`. |
| Cell types merge that shouldn't | CCA over-correcting; try `rpca`, or reduce `k_anchor`. |
| Query cells confidently mislabelled | A type missing from the reference. Check `prediction.score.max`. |
| Cluster count differs from Seurat | Expected. On ifnb RPCA, Seurat's deeper modularity search buys 0.17 % modularity by splitting CD14 Mono **on batch**; shanuz's coarser partition scores ARI 0.92 to the annotations against Seurat's 0.74. |

## Reference

- [Batch integration](https://genomicai.github.io/shanuz/tutorials/integration_vignette/) — Harmony/CCA/RPCA on ifnb, with the scoring code.
- [Reference mapping](https://genomicai.github.io/shanuz/tutorials/refmap_vignette/) — panc8 celseq2 → smartseq2.
- [Anchor internals](https://genomicai.github.io/shanuz/tutorials/anchors_vignette/) — the anchors and embeddings themselves, against Seurat.
