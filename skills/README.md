# truecell agent skills

Ten skills that teach an LLM agent to use `truecell` correctly — the API contracts
that break code silently, the decisions each analysis step forces, and the
places truecell and R Seurat genuinely differ.

They are plain Markdown with YAML frontmatter, so they work as
[Claude Agent Skills](https://docs.claude.com/en/docs/agents-and-tools/agent-skills/overview)
and as context for any other model.

## The set

| Skill | Load it for |
|---|---|
| [`truecell`](truecell/SKILL.md) | **Start here.** Install, the six API contracts, the canonical pipeline, routing. Bundles the full [API map](truecell/reference/api-map.md) and [object model](truecell/reference/object-model.md). |
| [`truecell-workflow`](truecell-workflow/SKILL.md) | A standard scRNA-seq run: QC thresholds, LogNormalize vs SCTransform, how many PCs, resolution, annotation. |
| [`truecell-differential-expression`](truecell-differential-expression/SKILL.md) | Marker genes, the eight `test_use` options, pseudobulk, conserved markers. |
| [`truecell-integration`](truecell-integration/SKILL.md) | Batch correction (Harmony/CCA/RPCA), label transfer, reference mapping, and scoring whether it worked. |
| [`truecell-multimodal`](truecell-multimodal/SKILL.md) | CITE-seq + WNN, cell hashing, pooled CRISPR (Mixscape). Includes the CLR `margin` rule. |
| [`truecell-spatial`](truecell-spatial/SKILL.md) | Xenium / Visium / CosMx / MERSCOPE, niches, spatially variable features, spatial plots. |
| [`truecell-at-scale`](truecell-at-scale/SKILL.md) | Leverage sketching and on-disk `LazyMatrix`, for data that doesn't fit in RAM. |
| [`truecell-plotting`](truecell-plotting/SKILL.md) | All 17 plotting functions, their Seurat equivalents, and headless saving. |
| [`truecell-from-seurat`](truecell-from-seurat/SKILL.md) | Porting R Seurat code, and comparing the two tools' numbers honestly. |
| [`truecell-dev`](truecell-dev/SKILL.md) | Contributing to truecell itself: tests, lint, docs, the fidelity apparatus, release conventions. |

Each `SKILL.md` stands alone. The router skill points at the others but does not
depend on them being loaded.

## Using them

### Claude Code

`.claude/skills` in this repo is a symlink to this directory, so the skills are
discovered automatically when Claude Code runs here. To use them from another
project:

```bash
ln -s /path/to/truecell/skills ~/.claude/skills/truecell
```

Or copy individual skill directories into `.claude/skills/`.

### Claude.ai / Projects

Upload the `SKILL.md` files (and `truecell/reference/*.md`) as project knowledge.
Names and descriptions in the frontmatter are what make the right one surface.

### Any other model

Concatenate what the task needs — the router plus one domain skill is usually
enough, and the whole set is small:

```bash
cat skills/truecell/SKILL.md skills/truecell-workflow/SKILL.md
cat skills/truecell/reference/api-map.md          # when parameter names matter
```

## Keeping them true

Every signature, default and measured number in these files came from the
package or from a recorded R comparison, not from memory. When the API changes,
the places to re-derive are:

```bash
python -c "import truecell, inspect; print([n for n in truecell.__all__])"
python -c "import truecell, inspect; print(inspect.signature(truecell.find_markers))"
```

and, for the fidelity claims, <https://genomicai.github.io/truecell/fidelity/>.
