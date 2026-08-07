---
name: truecell-dev
description: Use when contributing to the truecell package itself rather than using it — running the test suite, ruff and mypy, building the docs site, writing or re-running a tutorial's R-comparison, adding a declared band, and the repo's release and changelog conventions. Carries the project's verification standard, which is stricter than "the tests pass".
---

# Developing truecell

For working **on** the package. To use it, load `truecell` instead.

Repo: <https://github.com/GenomicAI/truecell> · Python 3.12+ · hatchling · MIT.

## Setup

```bash
git clone https://github.com/GenomicAI/truecell.git
cd truecell
uv sync --all-extras --locked
source .venv/bin/activate
```

**`--locked`, not `uv pip install -e ".[all]"`.** The lock is what CI installs
and what the committed tutorial figures were drawn with; resolving fresh gives
you a different scientific stack and the figures move for reasons no anchor
reports. Regenerate figures only from a locked environment.

`[all]` = analysis + anndata + integration + deseq2 + dev + docs. The docs
toolchain is in it **on purpose**: `tests/test_docs.py` builds the site with
`--strict`, and that check can only run if mkdocs is in the same environment CI
tests in.

## Layout

| Path | What |
|---|---|
| `truecell/` | The package. 53 source modules; `spatial/`, `compat/`, `mixins/` are subpackages. |
| `tests/` | ~70 test modules. |
| `tutorials/` | 18 vignettes + their Python and R scripts + `figures_*/`. Also the R-comparison apparatus. |
| `docs/` | MkDocs site. `docs/tutorials` is a **symlink** to `../tutorials`. |
| `tools/` | Build-time helpers: `griffe_sphinx_roles.py`, `mkdocs_html_relpaths.py`. |
| `.github/workflows/` | `ci.yml`, `docs.yml`. |

## The checks

```bash
pytest tests/ -q          # 972 passed, 25 skipped at 8a2b523
ruff check truecell         # clean
ruff check .              # clean (tests/ and tutorials/ ignore E402 by scope)
mypy                      # clean, 53 source files
mkdocs build --strict     # the docs check CI runs
```

Scope is pinned in `pyproject.toml` (`[tool.mypy] files = ["truecell"]`), so a bare
`mypy` checks exactly what CI checks. `ruff` takes scope from the command line —
`ruff check truecell` and `ruff check .` are both correct and easy to confuse when
comparing counts.

`E402` is ignored for `tests/*` and `tutorials/*` by scope: each file bootstraps
`sys.path` before importing truecell so it runs straight from a clone, which puts
every later import below a statement. That is one deliberate line per file, not
121 scattered `# noqa`s.

## CI

`ci.yml` — matrix `3.12` / `3.13`. 3.12 is the floor `requires-python` declares,
and the bottom leg is the only thing that proves the floor still works. 3.14 is
deliberately absent: `harmonypy` publishes manylinux wheels only through cp313,
and letting the resolver avoid a source build makes it backtrack to harmonypy
0.2.0, which drags in torch.

The floor tracks [SPEC 0](https://scientific-python.org/specs/spec-0000/) —
numpy/scipy/pandas/scikit-learn's own three-year window — not CPython's longer
EOL calendar.

ruff and mypy run **advisory** in CI (`|| true`). Both are clean now; that is
maintained by running them locally, not by the gate.

A separate `build` job builds the sdist and wheel, runs `twine check`, then
installs the wheel **clean into `/tmp`** and asserts `truecell.__version__` matches
`pyproject.toml`. That catches what the editable test job cannot: whether the
wheel is complete, whether a plain `pip install truecell` works without the
optional scientific stack, and whether the metadata version is right.

`docs.yml` builds with `--strict` and deploys to Pages.

`main` is protected. A PR needs four passing checks: `test (3.12)`,
`test (3.13)`, `Build and verify the distributions`, `Build the site`.

## The fidelity apparatus

This is what the repo is actually built around. Each tutorial is **three files**:

| File | Role |
|---|---|
| `tutorials/<name>_tutorial.py` | Runs the analysis in Python, writes a numeric handoff |
| `tutorials/<name>_verify.R` | Runs the same analysis under real Seurat, writes its own |
| `tutorials/<name>_vignette.md` | The write-up, both sides' code, the comparison |

```bash
python tutorials/pbmc3k_de_tutorial.py       # writes the handoff
Rscript tutorials/pbmc3k_de_verify.R         # reads it, writes R's answers
python tutorials/pbmc3k_de_tutorial.py --report   # compares; exits non-zero if it should
```

**Order matters and the tools enforce it.** Python first — it writes the handoff
R reads. References were taken on **Seurat 5.5.1**. Data caches to
`~/.truecell_data/` (~770 MB for the full set).

Two datasets (`ifnb`, `panc8`) exist only as curated SeuratData `.rda` objects
and need a one-time bridge so both languages read byte-identical counts:

```bash
Rscript tutorials/export_seuratdata.R ifnb     # ~394 MB, first run only
Rscript tutorials/export_seuratdata.R panc8    # ~117 MB
```

### Anchors, not vibes

Each comparison declares specific quantities. An **anchor** is one named number
that must match, at a stated tolerance or exactly. "91 of 91 anchors match, no
tolerance" is checkable line by line; "the object model works" is not.

Two comparison designs, chosen deliberately: some tutorials pin nothing (both
sides run their own pipeline end to end — what a user experiences); others pin
the cells and features and hand them to R, so a difference in an earlier stage
cannot pose as a difference in the one under test.

### Bands

`tutorials/bands.py` holds the numbers that legitimately move, as `Band` objects
carrying a range **and the reason for it**. `--report` prints a verdict per band
and exits non-zero outside one.

Two rules, both load-bearing:

- **Every band came from a sweep, not one run.** The JackStraw band is 60 seeds;
  the DESeq2 band is 20 resampled replicate splits.
- **A missing measurement fails.** `Band.holds(nan)` is `False` — a measurement
  that quietly vanished is exactly how a stale reference goes unnoticed.

`bands.py` also carries `StaleReferenceError` and the handoff guards. Those exist
because a Python run from 25 July was compared against R files from 19 July,
taken on a *different* cluster assignment, and printed a full parity table
showing `wilcox` at 48 of 50. Against a matching R run it is 50 of 50. **The
regression was in the reference.** Re-run Python and then R before trusting any
report.

### Tutorial tests

```bash
pytest tests/test_tutorial_marker_tables.py                      # fast, no data, in CI
TRUECELL_TUTORIAL_SMOKE=1 pytest tests/test_tutorial_smoke.py -v   # real runs, needs cached data
```

The second is **opt-in rather than skip-when-missing**, so a skip always means
nobody asked, never that it passed. Run it before cutting a release — a green
unit suite says nothing about whether the tutorials still work end to end.

The `tutorials` CI job covers the **PBMC 3k slice only** (11 tests, dataset
cached, a skip counts as a failure). Everything needing one of the other eight
datasets — ~200 MB in total — still runs nowhere but a developer's machine, so
the pre-release run above is not optional. The rule exists because it was
skipped: 1.0.0 shipped with the platelet cluster captioned "DC" in the pbmc3k
headline figure, under a guard that was correct and had never executed.

## The verification standard

Higher than "the tests pass", and it is why the port's claims hold up.

- **Test the integration, not the component.** Green suites here have repeatedly
  contained checks that were doing nothing.
- **Mutation-test every new guard.** Break the thing on purpose and confirm the
  test goes red. A guard that has never failed is a hypothesis.
- **A documented caveat is where a defect hides.** Twice a genuine bug sat behind
  a comment explaining it as an expected language difference — the SCTransform
  model was one, JackStraw's null another. A "known difference" nobody has
  measured is not a finding.
- **A fix can make the headline number worse and still be the fix.** Sketching's
  `project_data` scored *above* Seurat while broken. A divergence that flatters
  the port is a reason to look harder.
- **Prove RNG differences distribution-against-distribution** over matched seeds.
  Single-run pairs were actively misleading on the sketch composition.
- **Read CI logs for counts, not for green ticks.**

The tutorials are the debugging apparatus, not a showcase: across the suite they
have found and fixed 55+ defects, including one in R Seurat itself (the Visium
radius).

## Docs

```bash
mkdocs serve            # watches docs/, truecell/ and tools/
mkdocs build --strict
```

- `docs/tutorials` is a **symlink** to `../tutorials`, so one copy of each
  vignette serves both the site and the repo page. That is also why
  `exclude_docs` in `mkdocs.yml` filters the tutorial working directory's
  intermediates (CSV/JSON/npy/stores).
- The `toc` slugify is set to **GitHub's**, not python-markdown's, so
  cross-vignette heading links work in both places.
- `tools/mkdocs_html_relpaths.py` is a build hook that re-anchors raw-HTML
  `<img src>` paths the way MkDocs already re-anchors Markdown-syntax ones. The
  vignettes put the R and Python figures side by side in HTML tables, which
  MkDocs passes through untouched. Keep figure paths **relative to the source
  file** — a `../` that reaches out of `tutorials/` fixes the site and breaks
  GitHub, and `tests/test_docs.py` fails on it.
- mkdocstrings reads annotations **statically** through griffe, so the 19
  `if TYPE_CHECKING` annotations (matplotlib's `Figure`, `Neighbor`, `Truecell`)
  render and cross-link even though `typing.get_type_hints()` raises `NameError`
  on all of them.

## Conventions

- **Docstrings are NumPy style** (`Parameters ----------`), and they are the
  primary API documentation — the site renders them. Record the *decision*: which
  of Seurat's code paths a function follows, where a default was chosen to agree
  with `Seurat::` rather than the wider Python ecosystem, and where the two
  genuinely differ.
- **Comments explain why, not what.** The existing ones are long where the reason
  is non-obvious; match that.
- **`CHANGELOG.md` is the authority on what shipped**, `ROADMAP.md` on what is
  planned. They are not the same thing — a milestone landing on `main` is not a
  release. `tests/test_packaging.py` cross-checks release tags against the
  changelog, which is why CI checks out with `fetch-depth: 0`.
- One PR per coherent change, with the changelog entry in it.
