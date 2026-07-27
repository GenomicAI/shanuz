"""Make raw-HTML ``<img>`` paths survive the same rewrite Markdown images get.

MkDocs resolves relative links written in Markdown syntax against the *source*
file and rewrites them to point at the *built* page, which under
``use_directory_urls`` sits one directory deeper than its source:
``tutorials/svf_vignette.md`` builds to ``tutorials/svf_vignette/index.html``,
so ``![](figures_svf/x.png)`` correctly becomes ``../figures_svf/x.png``.

Raw HTML in a Markdown document is passed through untouched, attributes and
all. That is normally harmless, but the vignettes use ``<img>`` inside HTML
tables to put the R figure and the Python figure side by side — most of the
figure references on this site are written that way. Untouched, every one of
them resolves a directory too deep and 404s, which is what happened to the ten
vignettes built entirely out of those tables.

Rewriting the source to ``../figures_svf/x.png`` would fix the site and break
GitHub, where the vignettes are also read and where the path *is* relative to
the source file. So the rewrite happens here, at build time, against the same
target the Markdown treeprocessor uses: the site keeps working and the repo
page keeps working, from one correct path in the source.

A reference with no file behind it is warned about rather than silently
rewritten, so ``mkdocs build --strict`` in CI fails on it.
"""
from __future__ import annotations

import logging
import posixpath
import re

log = logging.getLogger("mkdocs.hooks.html_relpaths")

# `src` on an <img>, captured with its quotes so the replacement can put them
# back untouched. Deliberately not a general HTML parser: the vignettes use
# exactly one raw-HTML tag with a path in it, and a regex that only matches
# that tag cannot quietly start rewriting something else.
_IMG_SRC = re.compile(r'(?P<head><img\b[^>]*?\bsrc=)(?P<quote>["\'])(?P<src>[^"\']+)(?P=quote)')

# Anything already absolute, remote, inline, or a bare fragment is left alone.
_ABSOLUTE = re.compile(r"^(?:[a-zA-Z][\w+.-]*:|//|/|#)")


def target_of(src: str, src_uri: str) -> str:
    """Where a source-relative path lands, as a path from the docs root."""
    src_dir = posixpath.dirname(src_uri)
    return posixpath.normpath(posixpath.join(src_dir, src) if src_dir else src)


def rewrite(src: str, src_uri: str, url: str) -> str:
    """Re-anchor one path from the source file's directory onto the page's URL.

    ``url`` is ``page.file.url``: ``tutorials/x/`` with directory URLs and
    ``tutorials/x.html`` without. ``dirname`` gives the directory the browser
    resolves against in either case.
    """
    return posixpath.relpath(target_of(src, src_uri), posixpath.dirname(url) or ".")


def on_page_markdown(markdown: str, page, config, files, **kwargs) -> str:
    """Rewrite before rendering, not after.

    By ``on_page_content`` the Markdown images are already ``<img>`` tags and
    already re-anchored, and a rewrite there would move them a second time.
    Here the only ``<img>`` in the document is the raw HTML one — python-markdown
    stashes those verbatim, which is why the treeprocessor never sees them and
    why editing them at this stage survives rendering untouched.
    """
    def replace(match: re.Match) -> str:
        src = match.group("src")
        if _ABSOLUTE.match(src):
            return match.group(0)
        target = target_of(src, page.file.src_uri)
        if files.get_file_from_path(target) is None:
            log.warning(
                "%s: <img src=%r> has no file behind it (looked for %r)",
                page.file.src_uri, src, target,
            )
            return match.group(0)
        new = rewrite(src, page.file.src_uri, page.file.url)
        return f'{match.group("head")}{match.group("quote")}{new}{match.group("quote")}'

    return _IMG_SRC.sub(replace, markdown)
