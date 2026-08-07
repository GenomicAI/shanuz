"""Every tutorial CSV read that feeds a number must round-trip float64.

pandas' default CSV *reader* is not correctly rounded — it misparses about a
third of random doubles by an ULP. `to_csv` was never the problem; it already
writes the shortest round-trippable form. So a tutorial that writes a value,
hands it to R, reads both back and reports "these agree exactly" is partly
measuring the parser.

This is a lint, not a numeric check: it pins the *convention* across all
eighteen tutorials so the next `pd.read_csv` does not quietly reintroduce it.
The numeric consequence was small but real — fixing 37 call sites moved six
reported max-difference figures, most of them downward (PBMC 8k's `percent.mt`
5.773e-15 -> 5.329e-15, PBMC 3k's VST mean relative difference 1.548e-14 ->
4.973e-15). No declared band moved.

The R side of the same problem is worse and is not lintable from here: R's
`write.csv` renders 15 significant digits, and raising it does not help because
R's own `sprintf("%.17g")` is not correctly rounded either. Where bit-identity
is actually the question, the R script writes a C99 hex-float side table — see
`write_exact` in `tutorials/pbmc3k_de_verify.R`.
"""
import re
from pathlib import Path

import pytest

TUTORIALS = Path(__file__).parent.parent / "tutorials"

# Reads that carry no float into a comparison. Each is exempt for a stated
# reason; adding to this set should require one.
EXEMPT = {
    ("anchors_tutorial.py", "metadata.csv"): "cell labels, no floats compared",
    ("visium_tutorial.py", "header=header"): "wrapped in len() — a row count",
    ("pbmc3k_de_tutorial.py", "dtype=str"): "the hex reader; floats via float.fromhex",
}

READ_CSV = re.compile(r"read_csv\(")


def _call_text(source: str, start: int) -> str:
    """The text of the call beginning at `start`, to its matching paren."""
    depth, i, quote = 0, start, None
    while i < len(source):
        ch = source[i]
        if quote:
            if ch == "\\":
                i += 2
                continue
            if ch == quote:
                quote = None
        elif ch in "\"'":
            quote = ch
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return source[start:i + 1]
        i += 1
    return source[start:]


def _offending_calls(path: Path) -> list[str]:
    source = path.read_text()
    bad = []
    for m in READ_CSV.finditer(source):
        call = _call_text(source, m.end() - 1)
        if "float_precision" in call:
            continue
        if any(name == path.name and marker in call
               for (name, marker) in EXEMPT):
            continue
        line = source[:m.start()].count("\n") + 1
        bad.append(f"{path.name}:{line}  {' '.join(call.split())[:90]}")
    return bad


@pytest.mark.parametrize(
    "path", sorted(TUTORIALS.glob("*.py")), ids=lambda p: p.name)
def test_read_csv_round_trips_float64(path):
    bad = _offending_calls(path)
    assert not bad, (
        "pd.read_csv without float_precision=\"round_trip\" misparses ~1/3 of "
        "doubles by an ULP, so any exactness claim built on these values is "
        "partly measuring the parser:\n  " + "\n  ".join(bad))


def test_the_exemptions_still_exist():
    """A stale exemption would silently widen this lint.

    If one of the exempt calls is edited away, the entry stops matching anything
    and quietly permits a real offender in the same file.
    """
    for (name, marker), reason in EXEMPT.items():
        path = TUTORIALS / name
        assert path.exists(), f"{name} is gone; drop its exemption ({reason})"
        assert marker in path.read_text(), (
            f"exemption {name}:{marker!r} ({reason}) no longer matches anything "
            "— remove it or update the marker")
