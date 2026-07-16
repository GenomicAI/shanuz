"""Packaging and release metadata (v0.10.0).

Nothing here tests behaviour. These guard the seams where one fact is written
down twice — pyproject.toml, the installed distribution's metadata, the git
tags, and CHANGELOG.md — and can drift apart without a single behavioural test
noticing.
"""

from __future__ import annotations

import re
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest
from packaging.version import InvalidVersion, Version

import shanuz

try:  # stdlib from 3.11; `tomli` backfills 3.10 via the [dev] extra
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - only on Python 3.10
    import tomli as tomllib

REPO_ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = REPO_ROOT / "pyproject.toml"
CHANGELOG = REPO_ROOT / "CHANGELOG.md"

FALLBACK_VERSION = "0.0.0+unknown"

# "## [0.2.0] - 2026-07-05" or "## [Unreleased]" -> ("0.2.0", "2026-07-05") / ("Unreleased", "")
_HEADING = re.compile(r"^## \[([^\]]+)\](?:\s+-\s+(\d{4}-\d{2}-\d{2}))?\s*$", re.M)


def _changelog_headings() -> list[tuple[str, str]]:
    return _HEADING.findall(CHANGELOG.read_text())


# ----------------------------------------------------------------------
# __version__
# ----------------------------------------------------------------------


def test_version_agrees_with_the_installed_distribution():
    """``__version__`` and the distribution metadata report the same string.

    Note what this does *not* pin, despite being the obvious place to look for
    it: hard-code ``__version__ = "0.2.0"`` again and this still passes, because
    the literal and the metadata agree. The *mechanism* is pinned by
    ``test_version_falls_back_when_the_distribution_is_missing`` — a literal
    cannot fall back. Kept for the case where the two genuinely disagree.
    """
    from importlib.metadata import version

    assert shanuz.__version__ == version("shanuz")


def test_version_is_pep440():
    assert Version(shanuz.__version__)


def test_version_matches_pyproject():
    """Catches a stale editable install — a hazard the metadata lookup introduced.

    ``__version__`` resolves through the *installed* dist-info, which is a
    snapshot written at install time: the directory is literally named
    ``shanuz-<version>.dist-info``. Editing pyproject.toml does not rewrite it,
    so an editable install keeps reporting the old number until someone
    reinstalls (verified: pyproject at 9.9.9 still imported as 0.2.0). The
    hard-coded string this replaced could not drift this way, so the single
    source of truth is only worth having if something checks it is still single.
    """
    declared = tomllib.loads(PYPROJECT.read_text())["project"]["version"]
    assert shanuz.__version__ == declared, (
        f"__version__ is {shanuz.__version__!r} but pyproject.toml declares "
        f"{declared!r}. The installed metadata is stale — reinstall with "
        f"`uv pip install -e .`"
    )


def test_version_falls_back_when_the_distribution_is_missing():
    """A source tree on sys.path with nothing installed must still import.

    Load-bearing well beyond the fallback: this is the *only* test in the suite
    that fails if ``__version__`` reverts to a hard-coded literal, since a
    literal cannot fall back (verified by mutation). Don't drop it as a slow
    edge case — it is what holds pyproject.toml as the single source of truth.

    Subprocessed because the lookup runs once, at import. The in-process
    alternative is ``importlib.reload`` on a 50-module package mid-suite, which
    rebinds every class object and breaks ``isinstance`` in whatever runs next.
    """
    code = textwrap.dedent(
        f"""
        import importlib.metadata

        def _missing(name):
            raise importlib.metadata.PackageNotFoundError(name)

        importlib.metadata.version = _missing

        import shanuz
        assert shanuz.__version__ == {FALLBACK_VERSION!r}, shanuz.__version__
        print("ok")
        """
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "ok"


def test_the_fallback_version_is_itself_pep440():
    """So ``Version(shanuz.__version__)`` parses on that path too, and sorts low."""
    assert Version(FALLBACK_VERSION) < Version("0.1.0")


# ----------------------------------------------------------------------
# CHANGELOG.md
# ----------------------------------------------------------------------


def test_changelog_exists_and_leads_with_unreleased():
    assert CHANGELOG.is_file(), "CHANGELOG.md is missing from the repo root"
    headings = [name for name, _ in _changelog_headings()]
    assert headings, "CHANGELOG.md has no `## [version]` headings"
    assert headings[0] == "Unreleased", "`## [Unreleased]` must come first"


def test_changelog_documents_the_current_version():
    """Whatever ``pip install shanuz`` reports has to be findable in the changelog."""
    headings = [name for name, _ in _changelog_headings()]
    assert shanuz.__version__ in headings, (
        f"__version__ is {shanuz.__version__!r} with no `## [{shanuz.__version__}]` "
        f"section; the changelog has {headings}"
    )


def test_changelog_releases_are_valid_versions_in_descending_order():
    versions = []
    for name, _ in _changelog_headings():
        if name == "Unreleased":
            continue
        try:
            versions.append(Version(name))
        except InvalidVersion:
            pytest.fail(f"`## [{name}]` is not a PEP 440 version")
    assert versions == sorted(versions, reverse=True), "releases are out of order"


def test_changelog_releases_are_dated():
    for name, date in _changelog_headings():
        if name != "Unreleased":
            assert date, f"`## [{name}]` is missing its `- YYYY-MM-DD` date"


def test_every_tagged_release_is_in_the_changelog():
    """A tag with no entry is a release nobody wrote down."""
    try:
        proc = subprocess.run(
            ["git", "tag", "--list", "v*"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
    except OSError:  # pragma: no cover - git absent
        pytest.skip("git is not available")
    if proc.returncode != 0:  # pragma: no cover - not a checkout
        pytest.skip("not a git checkout")

    tagged = {tag[1:] for tag in proc.stdout.split() if tag.startswith("v")}
    if not tagged:  # pragma: no cover - shallow clone with no tags
        pytest.skip("no release tags visible")

    documented = {name for name, _ in _changelog_headings()}
    missing = sorted(tagged - documented)
    assert not missing, f"tagged but undocumented in CHANGELOG.md: {missing}"
