"""The documentation site cannot silently rot.

Four things go stale on their own and none of them break anything at import
time, which is exactly why they need asserting:

* a new export lands in ``shanuz.__all__`` and nobody adds it to an API page;
* a function is renamed and the ``:::`` directive pointing at it keeps its old
  name — mkdocstrings then renders an empty page section;
* a nav entry outlives the file it names;
* a vignette's figure is moved or regenerated under a new name.

The strict build at the end catches all four *and* every broken cross-reference,
but it needs the docs toolchain installed. These tests do not, so the checks
that matter most still run on a plain ``pip install -e ".[all]"``.
"""
import importlib
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import shanuz  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
API = DOCS / "api"

_DIRECTIVE = re.compile(r"^:::\s+(?P<target>[\w.]+)\s*$", re.MULTILINE)
# Markdown images, minus the absolute and remote ones.
_IMAGE = re.compile(r"!\[[^\]]*\]\((?P<src>(?!https?:|/)[^)\s]+)")


def api_targets() -> dict[str, Path]:
    """Every ``::: shanuz.x.y`` on the API pages, mapped to the page it is on."""
    found: dict[str, Path] = {}
    for page in sorted(API.glob("*.md")):
        for match in _DIRECTIVE.finditer(page.read_text()):
            found[match.group("target")] = page
    return found


def resolve(dotted: str):
    """Import-and-getattr a dotted path, module or attribute."""
    try:
        return importlib.import_module(dotted)
    except ImportError:
        module, _, attr = dotted.rpartition(".")
        return getattr(importlib.import_module(module), attr)


# ---------------------------------------------------------------------------
# The API pages against the package
# ---------------------------------------------------------------------------

def defining_module(name: str, obj: object) -> str:
    """Where an export really lives.

    `__module__` answers for functions and classes. Data constants like
    `CC_GENES` have none, so find the submodule that holds the same object.
    """
    module = getattr(obj, "__module__", None)
    if module:
        return module
    if isinstance(obj, type(shanuz)):  # a re-exported module
        return obj.__name__
    import pkgutil

    for info in pkgutil.walk_packages(shanuz.__path__, "shanuz."):
        try:
            candidate = importlib.import_module(info.name)
        except ImportError:
            continue
        if getattr(candidate, name, None) is obj:
            return info.name
    return ""


def test_every_public_export_is_on_an_api_page():
    documented = set(api_targets())
    missing = []
    for name in shanuz.__all__:
        if name == "__version__":
            continue
        obj = getattr(shanuz, name)
        module = defining_module(name, obj)
        if isinstance(obj, type(shanuz)):
            # A re-exported module — `generics` is rendered whole, `plotting`
            # one function at a time. Either counts.
            covered = module in documented or any(t.startswith(f"{module}.") for t in documented)
            canonical = f"::: {module}"
        else:
            canonical = f"::: {module}.{name}"
            covered = f"{module}.{name}" in documented
        if not covered:
            missing.append(f"{name} (would be `{canonical}`)")
    assert not missing, (
        "public exports with no API page — add a `::: ` directive under docs/api/:\n  "
        + "\n  ".join(sorted(missing))
    )


def test_every_directive_points_at_something_that_exists():
    broken = []
    for target, page in sorted(api_targets().items()):
        try:
            resolve(target)
        except (ImportError, AttributeError) as exc:
            broken.append(f"{page.name}: ::: {target} ({exc.__class__.__name__})")
    assert not broken, "API directives naming objects that no longer exist:\n  " + "\n  ".join(broken)


def test_no_public_export_is_documented_twice():
    """Two pages rendering the same symbol makes one of the two anchors dead."""
    seen: dict[str, list[str]] = {}
    for page in sorted(API.glob("*.md")):
        for match in _DIRECTIVE.finditer(page.read_text()):
            seen.setdefault(match.group("target"), []).append(page.name)
    duplicated = {t: pages for t, pages in seen.items() if len(pages) > 1}
    assert not duplicated, f"symbols rendered on more than one page: {duplicated}"


# ---------------------------------------------------------------------------
# The site's own files
# ---------------------------------------------------------------------------

def nav_entries() -> list[str]:
    """The markdown paths named in mkdocs.yml's nav, without parsing the YAML.

    `mkdocs.yml` carries a `!!python/object/apply:` tag for the slugifier, which
    `yaml.safe_load` refuses and `yaml.unsafe_load` would execute. The nav is a
    flat list of `.md` paths either way, so match those directly.
    """
    text = (ROOT / "mkdocs.yml").read_text()
    nav = text[text.index("\nnav:"):]
    return re.findall(r"^\s*-\s+(?:[^:\n]+:\s*)?([\w./-]+\.md)\s*$", nav, re.MULTILINE)


def test_every_nav_entry_exists():
    entries = nav_entries()
    assert len(entries) > 25, f"nav parsing found only {len(entries)} pages; the regex is wrong"
    missing = [e for e in entries if not (DOCS / e).exists()]
    assert not missing, f"nav names files that are not in docs/: {missing}"


def test_every_api_page_is_in_the_nav():
    entries = set(nav_entries())
    orphans = [
        f"api/{p.name}" for p in sorted(API.glob("*.md")) if f"api/{p.name}" not in entries
    ]
    assert not orphans, f"API pages that no reader can navigate to: {orphans}"


def test_every_vignette_is_in_the_nav():
    entries = set(nav_entries())
    vignettes = {
        f"tutorials/{p.name}"
        for p in sorted((ROOT / "tutorials").glob("*.md"))
    }
    missing = sorted(vignettes - entries)
    assert not missing, f"vignettes missing from the site nav: {missing}"


def test_every_markdown_image_resolves():
    """Figures are committed, and a vignette linking a missing one is a 404."""
    broken = []
    pages = list(DOCS.glob("*.md")) + list(API.glob("*.md")) + list((ROOT / "tutorials").glob("*.md"))
    for page in pages:
        for match in _IMAGE.finditer(page.read_text()):
            src = match.group("src").split("#")[0]
            if not (page.parent / src).resolve().exists():
                broken.append(f"{page.relative_to(ROOT)} -> {src}")
    assert not broken, "markdown images with no file behind them:\n  " + "\n  ".join(broken)


def test_the_tutorials_symlink_still_points_at_the_tutorials():
    link = DOCS / "tutorials"
    assert link.is_symlink(), "docs/tutorials must stay a symlink, not become a copy"
    assert link.resolve() == (ROOT / "tutorials").resolve()


# ---------------------------------------------------------------------------
# The griffe extension
# ---------------------------------------------------------------------------

sys.path.insert(0, str(ROOT / "tools"))


@pytest.fixture(scope="module")
def roles():
    pytest.importorskip("griffe", reason="needs the docs toolchain")
    import griffe_sphinx_roles

    return griffe_sphinx_roles


def test_a_role_pointing_at_a_documented_symbol_becomes_a_link(roles):
    ext = roles.SphinxRolesExtension()
    out = ext._substitute(roles._ROLE.search(":func:`run_pca`"))
    assert out == "[`run_pca`][shanuz.reduction.run_pca]"


def test_a_role_pointing_at_a_private_helper_degrades_to_code(roles):
    """Never emit a link the site cannot resolve — --strict would fail the build."""
    ext = roles.SphinxRolesExtension()
    out = ext._substitute(roles._ROLE.search(":func:`_get_scaled_data`"))
    assert out == "`_get_scaled_data`"


def test_a_role_pointing_outside_the_package_degrades_to_code(roles):
    ext = roles.SphinxRolesExtension()
    out = ext._substitute(roles._ROLE.search(":class:`pandas.DataFrame`"))
    assert out == "`pandas.DataFrame`"


def test_the_tilde_prefix_shortens_the_link_text(roles):
    ext = roles.SphinxRolesExtension()
    out = ext._substitute(roles._ROLE.search(":class:`~shanuz.spatial.visium.VisiumV2`"))
    assert out == "[`VisiumV2`][shanuz.spatial.visium.VisiumV2]"


def test_parameter_prose_moves_out_of_the_type_slot(roles):
    doc = (
        "Do a thing.\n\n"
        "Parameters\n"
        "----------\n"
        "min_features : min features a cell must have to be kept\n"
        "assay        : assay name, wrapping onto\n"
        "               a second line\n\n"
        "Returns\n"
        "-------\n"
        "``seurat``, modified.\n"
    )
    out = roles._reflow_parameters(doc)
    assert "min_features\n    min features a cell must have to be kept" in out
    assert "assay\n    assay name, wrapping onto\n    a second line" in out
    # The section after Parameters must survive untouched.
    assert out.endswith("Returns\n-------\n``seurat``, modified.\n")


def test_the_reflow_never_drops_a_word_from_any_docstring(roles):
    """A transform that silently ate prose would be invisible on the rendered page."""
    import inspect
    import pkgutil

    checked = 0
    for info in pkgutil.walk_packages(shanuz.__path__, "shanuz."):
        try:
            module = importlib.import_module(info.name)
        except ImportError:
            continue
        for name in dir(module):
            obj = getattr(module, name, None)
            doc = getattr(obj, "__doc__", None)
            if not doc or getattr(obj, "__module__", None) != info.name:
                continue
            before = inspect.cleandoc(doc)
            after = roles._reflow_parameters(before)
            checked += 1
            assert "".join(before.split()).replace(":", "") == \
                   "".join(after.split()).replace(":", ""), f"{info.name}.{name} lost content"
    assert checked > 100, f"only inspected {checked} docstrings; the walk is broken"


# ---------------------------------------------------------------------------
# The build itself
# ---------------------------------------------------------------------------

def test_the_site_builds_strict(tmp_path):
    """--strict turns every broken link, anchor and docstring warning into a failure."""
    pytest.importorskip("mkdocs", reason="needs the docs toolchain")
    import subprocess

    result = subprocess.run(
        [sys.executable, "-m", "mkdocs", "build", "--strict", "--site-dir", str(tmp_path / "site")],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
