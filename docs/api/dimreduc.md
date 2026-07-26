# Dimensional reduction

Linear first, then the embeddings you look at. Every one of these writes a
[`DimReduc`](objects.md#shanuz.dimreduc.DimReduc) into `obj.reductions` under a
key, and records the features it actually used — not the features it was asked
for, which are not always the same set.

`jack_straw` is the permutation test for how many PCs to keep. It is worth
reading its docstring before trusting the number: R's `JackRandom` seeds each
replicate from its loop index and is therefore deterministic, while this one
seeds from its `seed` argument and moves. Across 60 seeds on PBMC 3k it keeps
12–15 PCs, mode 13, which is R's answer. That spread is asserted as a band, not
described in prose — see [Fidelity](../fidelity.md#bands).

## Linear

::: shanuz.reduction.run_pca

::: shanuz.reduction.run_spca

::: shanuz.reduction.run_ica

::: shanuz.glmpca.glm_pca

## Non-linear embeddings

::: shanuz.umap.run_umap

::: shanuz.reduction.run_tsne

## How many components to keep

::: shanuz.jackstraw.jack_straw

::: shanuz.jackstraw.score_jackstraw
