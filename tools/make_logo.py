#!/usr/bin/env python3
"""Generate the shanuz logo into ``docs/assets/logo/``.

The mark is a point cloud whose density traces a lowercase ``s``. Every dot is a
cell, and pointillism is the nod to Georges Seurat -- the painter the R package
this ports is named for. The ``s`` is two lobes meeting at the waist, each
carrying a different colour: two implementations, one shape.

Three decisions here were arrived at by looking at the result, and are worth not
re-litigating:

* **The spine is a spline, not swept arcs.** Two arcs -- the obvious
  construction for an ``s`` -- close into a ring and read as a "no entry" sign.
  It is a centripetal Catmull-Rom through a hand-placed skeleton instead, with
  exact 180-degree rotational symmetry about the waist. The symmetry is load
  bearing: it puts the waist at precisely half the arc length, so the two-lobe
  colour split needs no fudge factor.
* **The wordmark has no font.** It is monoline geometric lowercase built from
  the same circles as the mark and emitted as sampled polylines, so there is no
  font to license and none to fail to load. Its ``s`` reuses the mark's
  skeleton, split into the same two colours at the same waist.
* **The favicon is solid, not dotted.** Rendered at 16px and compared
  side by side, the dot pattern silts up into mush; the silhouette survives. The
  favicon and the small-size glyph therefore draw the same skeleton as one
  stroke, which keeps the family resemblance at sizes the dots cannot hold.

Colours are the documentation site's, from ``docs/stylesheets/extra.css``.

Regenerate with::

    python tools/make_logo.py

Add ``--png`` to also write the raster exports, which needs Chrome (for
rendering) and Pillow (for downsampling); both are optional and the SVGs are
the actual deliverable.
"""

from __future__ import annotations

import argparse
import math
import pathlib
import random
import shutil
import subprocess
import sys
import tempfile

# --- palette -----------------------------------------------------------------
PINE = "#17423a"        # --md-primary-fg-color
AMBER = "#b2560d"       # --md-accent-fg-color, light ground
AMBER_SOFT = "#eb9d5c"  # --md-accent-fg-color, dark ground
CREAM = "#f3f6f4"       # --md-primary-bg-color

# --- the `s` skeleton --------------------------------------------------------
# Normalised to the unit box, v = 0 at the top. Only the top half is written
# down; the bottom half is its 180-degree rotation about the waist at (.5, .5).
_S_TOP = [
    (0.88, 0.090),   # top terminal, upper right
    (0.60, 0.005),   # apex
    (0.24, 0.045),
    (0.03, 0.210),   # leftmost
    (0.17, 0.395),
    (0.50, 0.500),   # waist
]
S_SKELETON = _S_TOP + [(1 - u, 1 - v) for u, v in reversed(_S_TOP[:-1])]


def catmull_rom(pts, n_per_seg: int = 60, alpha: float = 0.5):
    """Centripetal Catmull-Rom through `pts`, with the endpoints duplicated."""
    p = [pts[0]] + list(pts) + [pts[-1]]
    out = []
    for i in range(len(p) - 3):
        p0, p1, p2, p3 = p[i:i + 4]

        def knot(ti, a, b):
            d = math.dist(a, b)
            return ti + (d ** alpha if d else 1e-6)

        t0 = 0.0
        t1 = knot(t0, p0, p1)
        t2 = knot(t1, p1, p2)
        t3 = knot(t2, p2, p3)
        for k in range(n_per_seg):
            t = t1 + (t2 - t1) * k / n_per_seg

            def lerp(a, b, ta, tb):
                w = 0.0 if tb == ta else (t - ta) / (tb - ta)
                return (a[0] + (b[0] - a[0]) * w, a[1] + (b[1] - a[1]) * w)

            a1, a2, a3 = lerp(p0, p1, t0, t1), lerp(p1, p2, t1, t2), lerp(p2, p3, t2, t3)
            b1, b2 = lerp(a1, a2, t0, t2), lerp(a2, a3, t1, t3)
            out.append(lerp(b1, b2, t1, t2))
    out.append(pts[-1])
    return out


UNIT_S = catmull_rom(S_SKELETON)


def spine(x0: float, y0: float, w: float, h: float):
    """The skeleton mapped into a box, as (x, y, t) with t the arc fraction."""
    pts = [(x0 + u * w, y0 + v * h) for u, v in UNIT_S]
    acc = [0.0]
    for i in range(1, len(pts)):
        acc.append(acc[-1] + math.dist(pts[i], pts[i - 1]))
    return [(x, y, s / acc[-1]) for (x, y), s in zip(pts, acc)]


# --- the mark ----------------------------------------------------------------
SPINE = spine(28.0, 20.0, 44.0, 60.0)      # the `s`, inside the 100-unit box
CROSS_LO, CROSS_HI = 0.43, 0.57            # the waist is t = 0.5, by symmetry


def _cloud(n: int, seed: int, sigma: float, r_lo: float, r_hi: float):
    """Scatter dots about the spine: perpendicular jitter dominates, so the
    stroke stays a stroke, and dots near the spine are the large ones."""
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        i = min(int(rng.random() * (len(SPINE) - 1)), len(SPINE) - 2)
        x, y, t = SPINE[i]
        tx, ty = SPINE[i + 1][0] - x, SPINE[i + 1][1] - y
        norm = math.hypot(tx, ty) or 1.0
        perp, along = rng.gauss(0, sigma), rng.gauss(0, sigma * 0.45)
        dx = (-ty / norm) * perp + (tx / norm) * along
        dy = (tx / norm) * perp + (ty / norm) * along
        edge = min(1.0, abs(perp) / (sigma * 2.1))
        r = r_lo + (r_hi - r_lo) * ((1 - edge) ** 0.6) * (0.82 + 0.36 * rng.random())
        out.append((x + dx, y + dy, r, t, edge))
    return out


def _field(n: int, seed: int, keep_out: float = 10.5):
    """A few ambient cells -- the rest of the dataset the `s` cluster sits in."""
    rng = random.Random(seed)
    coarse = SPINE[::6]
    out, tries = [], 0
    while len(out) < n and tries < n * 90:
        tries += 1
        x, y = rng.uniform(12, 88), rng.uniform(12, 88)
        if min(math.hypot(x - sx, y - sy) for sx, sy, _ in coarse) >= keep_out:
            out.append((x, y, 0.85 + rng.random() * 0.45))
    return out


def _lobe(t: float, edge: float, on_dark: bool):
    """`on_dark`: drawn for a dark page, so the tile is cream and the ink pine."""
    hi = PINE if on_dark else CREAM
    lo = AMBER if on_dark else AMBER_SOFT
    if t < CROSS_LO:
        fill = hi
    elif t > CROSS_HI:
        fill = lo
    else:
        fill = lo if (t - CROSS_LO) / (CROSS_HI - CROSS_LO) > 0.5 else hi
    return fill, round(max(0.3, 1.0 - 0.62 * edge ** 1.4), 3)


HEX = " ".join(
    f"{50 + 48 * math.cos(math.radians(a)):.2f},{50 + 48 * math.sin(math.radians(a)):.2f}"
    for a in range(-90, 270, 60)
)


def _ground(shape: str, on_dark: bool) -> str:
    fill = CREAM if on_dark else PINE
    rim = AMBER if on_dark else AMBER_SOFT
    if shape == "tile":
        return f'<rect width="100" height="100" rx="23" fill="{fill}"/>'
    if shape == "hex":
        return (f'<polygon points="{HEX}" fill="{fill}"/>'
                f'<polygon points="{HEX}" fill="none" stroke="{rim}" stroke-width="2.8"/>')
    if shape == "circle":
        return (f'<circle cx="50" cy="50" r="50" fill="{fill}"/>'
                f'<circle cx="50" cy="50" r="45.5" fill="none" stroke="{rim}"'
                f' stroke-width="1.3" opacity="0.6"/>')
    return ""


def mark_parts(shape="tile", on_dark=False, n=190, seed=5, ambient=12, sigma=2.5):
    amb_fill = PINE if on_dark else CREAM
    body = [f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{amb_fill}" opacity="0.24"/>'
            for x, y, r in (_field(ambient, seed + 3) if ambient else [])]
    for x, y, r, t, e in _cloud(n, seed, sigma, 0.8, 2.35):
        fill, op = _lobe(t, e, on_dark)
        body.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}"'
                    f' fill="{fill}" opacity="{op}"/>')
    return _ground(shape, on_dark), "\n    ".join(body)


def _svg(view: str, body: str, label: str = "shanuz") -> str:
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{view}" role="img"'
            f' aria-label="{label}">\n  {body}\n</svg>\n')


def mark(shape="tile", on_dark=False, **kw) -> str:
    bg, body = mark_parts(shape, on_dark, **kw)
    return _svg("0 0 100 100", f'{bg}\n  <g>\n    {body}\n  </g>')


def glyph(on_dark=False, shape="tile", stroke=13.0) -> str:
    """The same skeleton as one stroke. This is what survives at 16px."""
    pts = spine(26.0, 17.0, 48.0, 66.0)
    half = len(pts) // 2
    hi = PINE if on_dark else CREAM
    lo = AMBER if on_dark else AMBER_SOFT
    out = [_ground(shape, on_dark)]
    for seg, col in ((pts[:half + 1], hi), (pts[half:], lo)):
        d = " ".join(("M" if i == 0 else "L") + f"{x:.1f},{y:.1f}"
                     for i, (x, y, _) in enumerate(seg))
        out.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{stroke}"'
                   f' stroke-linecap="round" stroke-linejoin="round"/>')
    return _svg("0 0 100 100", "\n  ".join(out))


# --- the wordmark ------------------------------------------------------------
# Monoline geometric lowercase, y up from the baseline. Every letter's *outer*
# silhouette -- the path plus half the stroke either side -- spans the x-height,
# so the paths are inset by half a stroke top and bottom. Getting that inset
# onto the `s` but not the other letters left the `s` 15 units short.
SW = 15.0                              # stroke weight
BASE, XTOP = SW / 2, 100 - SW / 2      # baseline 7.5, x-height line 92.5
ATOP = 140 - SW / 2                    # ascender line 132.5
RAD = (XTOP - BASE) / 2                # 42.5: the radius of the round letters
ROUND_W = 2 * RAD                      # 85: their path width
S_W = 63.0                             # the `s` is narrower, as an `s` is
TRACK = 30.0
KERN = {"s": -9.0}                     # the `s`'s open lower terminal needs less


def _semi(cx, cy, r, a0, a1, n=48):
    return [(cx + r * math.cos(math.radians(a)), cy + r * math.sin(math.radians(a)))
            for a in (a0 + (a1 - a0) * i / (n - 1) for i in range(n))]


def letter(ch: str):
    """Returns (strokes, advance); a stroke is (points, role) in y-up coords."""
    mid = BASE + RAD                                   # 50: the shoulder line
    if ch in "hn":
        return ([([(0, BASE), (0, ATOP if ch == "h" else XTOP)], ""),
                 (_semi(RAD, mid, RAD, 180, 0) + [(ROUND_W, BASE)], "")], ROUND_W)
    if ch == "u":
        return ([([(0, XTOP), (0, mid)] + _semi(RAD, mid, RAD, 180, 360), ""),
                 ([(ROUND_W, XTOP), (ROUND_W, BASE)], "")], ROUND_W)
    if ch == "a":
        return ([(_semi(RAD, mid, RAD, 90, 450), ""),
                 ([(ROUND_W, XTOP), (ROUND_W, BASE)], "")], ROUND_W)
    if ch == "z":
        w = 70.0
        return ([([(0, XTOP), (w, XTOP), (0, BASE), (w, BASE)], "")], w)
    if ch == "s":
        pts = [(u * S_W, XTOP - v * (XTOP - BASE)) for u, v in UNIT_S]
        half = len(pts) // 2
        return ([(pts[:half + 1], "hi"), (pts[half:], "lo")], S_W)
    raise ValueError(f"no letterform for {ch!r}")


def wordmark_parts(ink: str, lo: str, word: str = "shanuz"):
    x, parts = 0.0, []
    for ch in word:
        strokes, adv = letter(ch)
        for pts, role in strokes:
            col = lo if role == "lo" else ink
            d = " ".join(("M" if j == 0 else "L") + f"{x + px:.1f},{140 - py:.1f}"
                         for j, (px, py) in enumerate(pts))
            parts.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{SW}"'
                         f' stroke-linecap="round" stroke-linejoin="round"/>')
        x += adv + TRACK + KERN.get(ch, 0.0)
    return "\n    ".join(parts), x - TRACK - KERN.get(word[-1], 0.0)


def wordmark(on_dark=False, word="shanuz") -> str:
    ink = CREAM if on_dark else PINE
    lo = AMBER_SOFT if on_dark else AMBER
    inner, w = wordmark_parts(ink, lo, word)
    return _svg(f"{-SW / 2} 0 {w + SW:.1f} 140", inner, word)


def lockup(shape="tile", on_dark=False, word="shanuz", solid=False) -> str:
    ink = CREAM if on_dark else PINE
    lo = AMBER_SOFT if on_dark else AMBER
    inner, w = wordmark_parts(ink, lo, word)
    if solid:
        art = glyph(on_dark, shape).split("\n", 1)[1].rsplit("</svg>", 1)[0].strip()
    else:
        bg, body = mark_parts(shape, on_dark)
        art = f"{bg}\n  <g>\n    {body}\n  </g>"

    box, pad, xheight = 100.0, 26.0, 40.0     # wordmark x-height: 40 of the 100 box
    sc = xheight / 100.0
    total = box + pad + (w + SW) * sc
    ty = box / 2 + xheight / 2 - 1.5 - 140 * sc
    return _svg(
        f"0 0 {total:.1f} 100",
        f'{art}\n  <g transform="translate({box + pad + SW / 2 * sc:.2f},{ty:.2f})'
        f' scale({sc:.4f})">\n    {inner}\n  </g>',
        word,
    )


# --- outputs -----------------------------------------------------------------
SVGS = {
    # primary: the pointillist mark
    "shanuz-mark.svg":             lambda: mark("tile"),
    "shanuz-mark-inverse.svg":     lambda: mark("tile", on_dark=True),
    # horizontal lockup, mark plus wordmark
    "shanuz-lockup.svg":           lambda: lockup("tile"),
    "shanuz-lockup-inverse.svg":   lambda: lockup("tile", on_dark=True),
    # the wordmark alone
    "shanuz-wordmark.svg":         lambda: wordmark(),
    "shanuz-wordmark-inverse.svg": lambda: wordmark(on_dark=True),
    # the R-ecosystem sticker shape, since this is a port of an R package
    "shanuz-hex.svg":              lambda: mark("hex"),
    "shanuz-hex-inverse.svg":      lambda: mark("hex", on_dark=True),
    # simplified: what reads below roughly 48px
    "shanuz-glyph.svg":            lambda: glyph(),
    "shanuz-glyph-inverse.svg":    lambda: glyph(on_dark=True),
    # the glyph with no tile of its own, for dropping onto a ground that is
    # already coloured -- the documentation header, which is pine in both the
    # light and the dark palette, so it is the cream-ink file that goes there.
    "shanuz-glyph-open.svg":         lambda: glyph(shape="none"),
    "shanuz-glyph-open-inverse.svg": lambda: glyph(shape="none", on_dark=True),
    "favicon.svg":                 lambda: glyph(),
}

# name -> (source svg, width, height); height None means square
PNGS = {
    "shanuz-mark-512.png":   ("shanuz-mark.svg", 512, None),
    "shanuz-mark-256.png":   ("shanuz-mark.svg", 256, None),
    "shanuz-mark-128.png":   ("shanuz-mark.svg", 128, None),
    "shanuz-hex-512.png":    ("shanuz-hex.svg", 512, None),
    "apple-touch-icon.png":  ("favicon.svg", 180, None),
    "favicon-48.png":        ("favicon.svg", 48, None),
    "favicon-32.png":        ("favicon.svg", 32, None),
    "favicon-16.png":        ("favicon.svg", 16, None),
    # Both grounds, because the README is read on GitHub and GitHub has a dark
    # mode; the `<picture>` there switches between these two.
    "shanuz-lockup-1200.png": ("shanuz-lockup.svg", 1200, None),
    "shanuz-lockup-inverse-1200.png": ("shanuz-lockup-inverse.svg", 1200, None),
}

CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"


def _render_png(svg_path: pathlib.Path, out: pathlib.Path, width: int) -> None:
    """Render at a whole multiple of the target size, then downsample.

    Chrome refuses window sizes below about 50px, and downsampling a large
    render is what a favicon pipeline does anyway -- it antialiases far better
    than asking the renderer for 16 pixels.

    Two things here were got wrong once, and both showed up in the exported
    lockup rather than in the square assets, which is what made them easy to
    miss:

    * **The window is an exact multiple of the output grid, and the SVG is
      given those same pixel dimensions outright** instead of being sized
      against the viewport in ``vw``/``vh``. A 1024-wide window for the lockup
      wants to be 271.33 tall and can only be 271, and laid out that way the
      art was placed nine pixels down inside it -- so the bottom of the tile
      fell outside the frame and its rounded corners were cut mid-curve.
      Sizing the element and the window to the same integers leaves nothing to
      disagree about. It was *not* the inline-element baseline gap, which is
      the obvious suspect and makes no difference here.
    * **The render is never smaller than the output.** A fixed 1024 meant the
      1200px lockup was being scaled *up*, which is what softened its edges.
    """
    from PIL import Image

    svg = svg_path.read_text()
    vb = svg.split('viewBox="', 1)[1].split('"', 1)[0].split()
    ratio = float(vb[3]) / float(vb[2])
    out_w, out_h = width, max(1, round(width * ratio))
    k = max(2, math.ceil(1024 / out_w))     # at least 2x, and at least ~1024px
    big_w, big_h = out_w * k, out_h * k

    head, rest = svg.split("<svg ", 1)
    sized = (f'{head}<svg style="position:absolute;top:0;left:0;'
             f'width:{big_w}px;height:{big_h}px" preserveAspectRatio="none" {rest}')
    with tempfile.TemporaryDirectory() as td:
        page = pathlib.Path(td) / "p.html"
        shot = pathlib.Path(td) / "p.png"
        page.write_text(f"<style>html,body{{margin:0;padding:0}}</style>{sized}")
        subprocess.run(
            [CHROME, "--headless", "--disable-gpu", f"--screenshot={shot}",
             f"--window-size={big_w},{big_h}", "--hide-scrollbars",
             "--default-background-color=00000000", str(page)],
            check=True, capture_output=True,
        )
        im = Image.open(shot).convert("RGBA")
    im.resize((out_w, out_h), Image.LANCZOS).save(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out", nargs="?", default=None,
                    help="output directory (default: docs/assets/logo)")
    ap.add_argument("--png", action="store_true", help="also write the raster exports")
    args = ap.parse_args()

    out = pathlib.Path(args.out) if args.out else \
        pathlib.Path(__file__).resolve().parent.parent / "docs" / "assets" / "logo"
    out.mkdir(parents=True, exist_ok=True)

    for name, fn in SVGS.items():
        (out / name).write_text(fn())
    print(f"wrote {len(SVGS)} SVGs to {out}")

    if not args.png:
        print("(re-run with --png for the raster exports)")
        return 0

    if not pathlib.Path(CHROME).exists() and not shutil.which("chrome"):
        print("PNG export needs Chrome; skipping", file=sys.stderr)
        return 1
    try:
        import PIL  # noqa: F401
    except ImportError:
        print("PNG export needs Pillow; skipping", file=sys.stderr)
        return 1

    for name, (src, width, _) in PNGS.items():
        _render_png(out / src, out / name, width)
    print(f"wrote {len(PNGS)} PNGs to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
