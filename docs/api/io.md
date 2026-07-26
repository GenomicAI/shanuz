# Loading data

## From 10x output

::: shanuz.io.read_10x

## Bundled datasets

`shanuz.datasets` stands in for R's `SeuratData`. Each loader downloads on first
call into `~/.shanuz_data/` and returns a ready `Shanuz` object; the whole set is
roughly 770 MB cached. These are the datasets the [tutorials](../tutorials/README.md)
run on, which is what makes each tutorial reproducible from a clean machine.

::: shanuz.datasets.pbmc3k

::: shanuz.datasets.pbmc8k

::: shanuz.datasets.cbmc_citeseq

::: shanuz.datasets.pbmc_hashing

::: shanuz.datasets.ifnb

::: shanuz.datasets.panc8

::: shanuz.datasets.thp1_eccite

::: shanuz.datasets.xenium_mouse_brain

::: shanuz.datasets.visium_mouse_brain

## AnnData interoperability

::: shanuz.compat.anndata.as_anndata

::: shanuz.compat.anndata.from_anndata
