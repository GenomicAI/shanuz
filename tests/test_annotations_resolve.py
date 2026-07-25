"""Every annotation in `shanuz/` must name something the module actually binds.

`from __future__ import annotations` turns every signature into a string that is
never evaluated, so an annotation may name a symbol that exists nowhere and
nothing complains — the module imports, the function runs, the tests pass. What
breaks is anything that *reads* the types: a checker, an IDE, or a documentation
generator resolving `-> "plt.Figure"` when no `plt` exists at module scope.

That is not hypothetical. Twenty-one annotations across four modules named
symbols bound only inside function bodies (`plt = _mpl()`) or imported only at
call time to dodge a circular import (`Graph`/`Neighbor`/`Shanuz`). The fix is a
`if TYPE_CHECKING:` block, which type checkers and static doc tools read and the
interpreter never executes — and which is easy to delete later while every other
test stays green. This module is the thing that would notice.
"""
import ast
import builtins
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
PKG = REPO_ROOT / "shanuz"

# Modules whose annotations depend on a `TYPE_CHECKING` block, and the symbol
# each one defers. Listed explicitly so deleting a block is a failure with a
# name attached rather than a silent loss of coverage.
DEFERRED = {
    "plotting.py": "Figure",
    "graph.py": "Neighbor",
    "neighbor.py": "Graph",
    "compat/anndata.py": "Shanuz",
}


def _module_scope_names(tree: ast.Module) -> set:
    """Names bound at module scope — the only scope an annotation can see.

    Deliberately does *not* recurse into function or class bodies. An earlier
    draft of this check walked every scope, which counted the function-local
    `plt = _mpl()` as a binding for the module-level annotation `plt.Figure`
    and so reported a clean run against the very tree that had all 21 defects.
    `if` and `try` bodies *are* followed, since that is where both
    `TYPE_CHECKING` blocks and optional-dependency imports live.
    """
    names = set(dir(builtins))

    def visit(body):
        for node in body:
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                for alias in node.names:
                    names.add((alias.asname or alias.name).split(".")[0])
            elif isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                names.add(node.name)
            elif isinstance(node, ast.Assign):
                for target in node.targets:
                    for sub in ast.walk(target):
                        if isinstance(sub, ast.Name):
                            names.add(sub.id)
            elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
                names.add(node.target.id)
            elif isinstance(node, ast.If):
                visit(node.body)
                visit(node.orelse)
            elif isinstance(node, ast.Try):
                visit(node.body)
                visit(node.orelse)
                visit(node.finalbody)
                for handler in node.handlers:
                    visit(handler.body)

    visit(tree.body)
    return names


def _roots(expr: ast.expr) -> list:
    """Identifiers an annotation depends on, unwrapping strings and generics.

    `Optional[Graph]` depends on both `Optional` and `Graph`; `plt.Figure`
    depends on `plt`, not on `Figure`; `"Shanuz"` is parsed and then treated
    like the expression it spells.
    """
    if isinstance(expr, ast.Constant):
        if isinstance(expr.value, str):
            try:
                return _roots(ast.parse(expr.value, mode="eval").body)
            except SyntaxError:
                return []
        return []                                   # None, literal ints in Literal[...]
    if isinstance(expr, ast.Name):
        return [expr.id]
    if isinstance(expr, ast.Attribute):
        return _roots(expr.value)                   # np.ndarray -> np
    if isinstance(expr, ast.Subscript):
        return _roots(expr.value) + _roots(expr.slice)
    if isinstance(expr, ast.BinOp):                 # int | None
        return _roots(expr.left) + _roots(expr.right)
    if isinstance(expr, (ast.Tuple, ast.List)):
        return [r for elt in expr.elts for r in _roots(elt)]
    return []


def _annotations(tree: ast.Module):
    """Yield (root_name, owner, lineno) for every annotation in the module."""
    for node in ast.walk(tree):
        found = []
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.returns is not None:
                found.append(node.returns)
            args = node.args
            for arg in [*args.posonlyargs, *args.args, *args.kwonlyargs,
                        args.vararg, args.kwarg]:
                if arg is not None and arg.annotation is not None:
                    found.append(arg.annotation)
            owner = node.name
        elif isinstance(node, ast.AnnAssign) and node.annotation is not None:
            found.append(node.annotation)
            owner = ast.unparse(node.target)
        else:
            continue
        for expr in found:
            for root in _roots(expr):
                yield root, owner, expr.lineno


def _sources():
    # utf-8-sig rather than utf-8: no source carries a BOM any more (see
    # `test_no_source_file_starts_with_a_byte_order_mark`), but decoding one
    # away here keeps a stray BOM from surfacing as a confusing SyntaxError
    # inside an annotation test. The dedicated guard below is what names it.
    return [(p, ast.parse(p.read_text(encoding="utf-8-sig"), filename=str(p)))
            for p in sorted(PKG.rglob("*.py"))]


def test_no_source_file_starts_with_a_byte_order_mark():
    """A UTF-8 BOM is legal Python and breaks every tool that reads the text.

    CPython's import machinery decodes source as `utf-8-sig`, so a BOM never
    shows up as a runtime failure — `import`, `inspect.getsource`, ruff and mypy
    all handle `shanuz/command.py` without complaint, which is why one survived
    in it unnoticed. What breaks is the ordinary idiom for reading source back:

        ast.parse(Path(mod.__file__).read_text())

    `Path.read_text()` defaults to plain utf-8 and keeps the U+FEFF, and
    `ast.parse` then rejects the file with `invalid non-printable character`.
    That is not hypothetical — the AST walk in this very module had to decode
    around it, and any future contributor writing a codemod, a doc generator or
    a source-level lint hits the same wall on a file that looks fine.

    Scans the whole repo rather than just the package: the reason the BOM is
    worth removing has nothing to do with which directory it lives in.
    """
    offenders = [
        p.relative_to(REPO_ROOT)
        for p in sorted(REPO_ROOT.rglob("*.py"))
        if ".venv" not in p.parts and "build" not in p.parts
        and p.read_bytes().startswith(b"\xef\xbb\xbf")
    ]
    assert not offenders, (
        "Python sources beginning with a UTF-8 BOM — these import fine but "
        "cannot be read with Path.read_text() and parsed:\n  "
        + "\n  ".join(str(p) for p in offenders)
    )


def test_every_annotation_names_a_module_scope_binding():
    """The whole package, not a sampled subset — the defect was spread over four files."""
    sources = _sources()
    assert len(sources) > 40, f"only {len(sources)} modules found; is PKG right?"

    unresolved = [
        f"{path.relative_to(PKG.parent)}:{line} in {owner}() -> {root!r}"
        for path, tree in sources
        for root, owner, line in _annotations(tree)
        if root not in _module_scope_names(tree)
    ]
    assert not unresolved, (
        "annotations naming symbols bound nowhere at module scope — a type "
        "checker and any doc generator both resolve these to nothing:\n  "
        + "\n  ".join(unresolved)
    )


@pytest.mark.parametrize("relpath,symbol", sorted(DEFERRED.items()))
def test_deferred_symbols_stay_inside_type_checking(relpath, symbol):
    """The import must be in a `TYPE_CHECKING` block, not merely present.

    Hoisting it to module scope silences the checker just as well, so this
    asserts the placement rather than the resolution. For `plotting.py` the
    difference is a hard matplotlib dependency; for the other three it is an
    import cycle at interpreter start.
    """
    tree = ast.parse((PKG / relpath).read_text(encoding="utf-8-sig"))

    guarded, unguarded = False, False
    for node in tree.body:
        if isinstance(node, ast.If) and "TYPE_CHECKING" in ast.unparse(node.test):
            for sub in ast.walk(node):
                if isinstance(sub, (ast.Import, ast.ImportFrom)):
                    if any(a.asname == symbol or a.name == symbol for a in sub.names):
                        guarded = True
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            if any(a.asname == symbol or a.name == symbol for a in node.names):
                unguarded = True

    assert guarded, f"{relpath} no longer defers {symbol!r} under TYPE_CHECKING"
    assert not unguarded, (
        f"{relpath} imports {symbol!r} at module scope; the TYPE_CHECKING block "
        f"exists precisely so this import never runs"
    )


def test_importing_shanuz_still_does_not_import_matplotlib():
    """`shanuz/__init__.py` imports `plotting` eagerly, so this is load-bearing.

    matplotlib is optional — `_mpl()` raises a pip-install message when it is
    absent. Moving `from matplotlib.figure import Figure` out of the
    `TYPE_CHECKING` block would turn `import shanuz` into a hard requirement for
    it on every install, and no other test in the suite would notice, because
    CI has matplotlib installed.

    Run in a subprocess: by the time this file executes, the rest of the suite
    has long since put matplotlib in `sys.modules`.
    """
    code = (
        "import sys; import shanuz; "
        "assert 'shanuz.plotting' in sys.modules, 'plotting no longer imported eagerly'; "
        "sys.exit(1 if 'matplotlib' in sys.modules else 0)"
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert proc.returncode == 0, (
        "importing shanuz now pulls in matplotlib, making an optional "
        f"dependency mandatory.\n{proc.stdout}{proc.stderr}"
    )


@pytest.mark.parametrize("relpath,symbol", sorted(DEFERRED.items()))
def test_deferred_symbols_are_absent_at_runtime(relpath, symbol):
    """The complement of the AST check: confirm the block really does not execute.

    If `TYPE_CHECKING` were ever truthy — or the block were replaced with a
    plain `try: import ... except ImportError:` — these names would appear in
    the module namespace and the runtime deferral would be gone.
    """
    import importlib

    name = "shanuz." + relpath[:-3].replace("/", ".")
    module = importlib.import_module(name)
    assert not hasattr(module, symbol), (
        f"{name}.{symbol} exists at runtime; the TYPE_CHECKING block is executing"
    )
