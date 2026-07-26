# Preprocessing

Counts to something you can do statistics on. Two routes, both of Seurat's:
log-normalize → find variable features → scale, or `sctransform` in one call.

`find_variable_features` has two selectors and they honour different arguments —
`"vst"` and `"disp"` take `nfeatures`, `"mvp"` takes the mean and dispersion
cutoffs. That is Seurat's behaviour, not a quirk of the port; the docstring says
which is which.

Verified against Seurat in [PBMC 3k](../tutorials/pbmc3k_tutorial.md) and, for the
regularized-NB route, per fitted gene in
[SCTransform](../tutorials/sctransform_vignette.md).

## Log-normalize workflow

::: shanuz.preprocessing.normalize_data

::: shanuz.preprocessing.find_variable_features

::: shanuz.preprocessing.scale_data

::: shanuz.preprocessing.percentage_feature_set

## Regularized negative binomial

::: shanuz.sctransform.sctransform
