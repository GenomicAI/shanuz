"""Make the package's two docstring conventions legible to mkdocstrings.

Both fixes are build-time translations rather than source rewrites, for the same
reason: the source form is the one a reader of the code and a caller of `help()`
sees, it is consistent across the package, and it is not wrong — it just is not
what a NumPy-format parser expects.

**Sphinx roles.** The docstrings cross-reference with ``:func:`run_pca``` and
friends, 129 times. mkdocstrings does not understand that syntax, so left alone
they render as the literal text ``:func:`` followed by a code span.

**Parameter descriptions.** Every one of the 370 documented parameters is
written ``name : what it does`` — prose after the colon, never a type, because
the types are in the signatures where they belong. NumPy format reads that slot
as the type, so without this the site prints "min features a cell must have to
be kept" in the Type column and leaves the description empty.

Roles are rewritten to mkdocs-autorefs links **only when the target is one this
site actually renders**. Anything else — private helpers filtered out of the
API pages, `pandas.DataFrame`, `scipy.sparse.load_npz` — degrades to a plain
code span. That rule is the point of the extension rather than an edge case:
autorefs treats an unresolvable target as a warning, CI builds with ``--strict``,
so a link this file is not sure of would fail the build. Better to render it as
code than to guess.
"""

from __future__ import annotations

import re
from typing import Any

from griffe import Extension, Object

# `:role:`target`` — optionally `~target`, meaning "display the last component
# only", and optionally `text <target>`, Sphinx's explicit-title form.
_ROLE = re.compile(
    r":(?:py:)?(?P<role>func|meth|class|obj|attr|mod|data|exc):"
    r"`(?P<body>[^`]+)`"
)
_EXPLICIT_TITLE = re.compile(r"^(?P<title>.*?)\s*<(?P<target>[^>]+)>$")

# A `Parameters` / `Other Parameters` header and its underline.
_PARAM_HEADER = re.compile(r"^(?P<name>Parameters|Other Parameters)[ \t]*$")
_UNDERLINE = re.compile(r"^-{3,}[ \t]*$")
# `name : prose`, at the section's own indentation. `*args` and `**kwargs` count.
_PARAM_ENTRY = re.compile(r"^(?P<name>\*{0,2}\w[\w ,]*?)[ \t]*:[ \t]?(?P<desc>\S.*)$")


def _reflow_parameters(text: str) -> str:
    """Move each parameter's prose off the type slot and onto its own line.

    ``name : what it does`` becomes ``name`` followed by an indented
    description, which is the shape NumPy format expects. Leaving the type slot
    empty is what makes griffe fall back to the signature annotation — the one
    place in this package where the type is actually declared.
    """
    lines = text.splitlines()
    out: list[str] = []
    i = 0
    while i < len(lines):
        header = _PARAM_HEADER.match(lines[i])
        if not (header and i + 1 < len(lines) and _UNDERLINE.match(lines[i + 1])):
            out.append(lines[i])
            i += 1
            continue
        out.append(lines[i])
        out.append(lines[i + 1])
        i += 2
        while i < len(lines):
            line = lines[i]
            # A blank line ends the section only when what follows is not an
            # indented continuation — some sections have a gap mid-list.
            if not line.strip():
                nxt = next((s for s in lines[i + 1:] if s.strip()), "")
                if not nxt or not nxt.startswith((" ", "\t")) and not _PARAM_ENTRY.match(nxt):
                    break
                out.append(line)
                i += 1
                continue
            if line[:1].strip() == "" :  # already indented: a continuation line
                out.append("    " + line.strip())
                i += 1
                continue
            entry = _PARAM_ENTRY.match(line)
            if entry is None:
                break
            out.append(entry.group("name").strip())
            out.append("    " + entry.group("desc").strip())
            i += 1
    rebuilt = "\n".join(out)
    return rebuilt + "\n" if text.endswith("\n") else rebuilt


def _render_targets() -> dict[str, str]:
    """Map every name a docstring might use to the anchor this site publishes.

    The anchors mkdocstrings emits are canonical paths — `truecell.reduction.run_pca`,
    not `truecell.run_pca` — because that is what the `:::` directives on the API
    pages name. So the index is built from the canonical path of each public
    export and of each public method on the exported classes, and a bare name
    resolves only when exactly one object answers to it.
    """
    import inspect

    import truecell

    canonical: dict[str, str] = {}
    ambiguous: set[str] = set()

    def offer(alias: str, path: str) -> None:
        if alias in canonical and canonical[alias] != path:
            ambiguous.add(alias)
        canonical.setdefault(alias, path)

    for name in truecell.__all__:
        obj = getattr(truecell, name, None)
        module = getattr(obj, "__module__", None)
        if obj is None or not module or not module.startswith("truecell"):
            continue
        path = f"{module}.{name}"
        offer(name, path)
        offer(path, path)
        offer(f"truecell.{name}", path)
        if inspect.isclass(obj):
            for attr, member in vars(obj).items():
                if attr.startswith("_") or not callable(member):
                    continue
                offer(attr, f"{path}.{attr}")
                offer(f"{name}.{attr}", f"{path}.{attr}")

    for alias in ambiguous:
        del canonical[alias]
    return canonical


class SphinxRolesExtension(Extension):
    """Translate Sphinx roles in every docstring griffe collects."""

    def __init__(self) -> None:
        self._targets: dict[str, str] | None = None

    @property
    def targets(self) -> dict[str, str]:
        # Built on first use, not in __init__: importing truecell while griffe is
        # still loading modules is a cycle waiting to happen.
        if self._targets is None:
            self._targets = _render_targets()
        return self._targets

    def _substitute(self, match: re.Match[str]) -> str:
        body = match.group("body").strip()
        title: str | None = None
        if (explicit := _EXPLICIT_TITLE.match(body)) is not None:
            title, body = explicit.group("title"), explicit.group("target").strip()

        abbreviate = body.startswith("~")
        body = body.lstrip("~")
        if title is None:
            title = body.rsplit(".", 1)[-1] if abbreviate else body

        target = self.targets.get(body)
        if target is None:
            return f"`{title}`"
        return f"[`{title}`][{target}]"

    def on_object(self, *, obj: Object, **kwargs: Any) -> None:
        docstring = obj.docstring
        if docstring is None:
            return
        value = docstring.value
        if ":" in value:
            value = _ROLE.sub(self._substitute, value)
        docstring.value = _reflow_parameters(value)
