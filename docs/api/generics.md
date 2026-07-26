# Generics

R's Seurat is built on S4 generics — `Cells(x)`, `Idents(x)`, `LayerData(x, ...)`
— that dispatch on whatever object you hand them. `shanuz.generics` is the
same surface as free functions, so code ported from R reads the way it did in R.

Everything here also exists as a method or property on the objects themselves.
Use whichever fits; they are the same code path.

```python
from shanuz import generics as g

g.cells(pbmc)            # Cells(pbmc)
g.idents(pbmc)           # Idents(pbmc)
g.layer_data(pbmc, "counts")   # LayerData(pbmc, layer = "counts")
```

::: shanuz.generics
    options:
      show_root_heading: false
      show_root_toc_entry: false
      members_order: alphabetical
      show_source: false
      heading_level: 2
