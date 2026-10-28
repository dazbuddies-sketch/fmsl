"""The open case of Đinh, Guo, Nguyễn and Phạm, "Limits of real bivariate rational functions"
(arXiv:2202.04889): the denominator vanishes along curves through the origin, not only at the
origin. Their paper computes the range of limit values when the zero is isolated (the range is
then a closed interval) and leaves this case open, noting the range may be disconnected in
[-oo, +oo].

For each example this script prints
  * FMSL's proven set of limit values at (0, 0), and
  * an independent exact check: the exact range of f on the circle of radius 1e-12, computed
    with a rational parametrization of the circle and exact one-variable range finding (no
    Newton-Puiseux, no zooming). As the radius shrinks these ranges converge to the set of limit
    values, so the two must agree.

    pip install fmsl
    python open_case.py
"""
import math, time
from fractions import Fraction as Fr

from fmsl.verify import limit_at_origin
from fmsl.space import parse3
from fmsl.radial import _polar
from fmsl.plane import image_1d
from fmsl.algebraic import pgcd, pdivmod, ptrim
from fmsl.cluster import _Inf

EXAMPLES = [
    # (formula, zero set of the denominator through the origin)
    ("x^2/(x^2 - y^2)", "the lines y = x and y = -x"),
    ("(x^2 + y^2)/(x*y)", "the two axes"),
    ("(x^2 + x*y + y^2)/(x^2 - 2*y^2)", "the lines x = sqrt(2) y and x = -sqrt(2) y"),
    ("(x^2 + y^2)^2/(x^3*y)", "the two axes"),
    ("(x^4 + y^4)/(x^2*y^2 - x^5)", "the y-axis and the cusp y^2 = x^3"),
    ("(x^2 - y^2)/(x^2 - 2*y^4)", "the parabolas x = sqrt(2) y^2 and x = -sqrt(2) y^2"),
    ("y^2/(x - y)^2", "the line y = x (a double zero)"),
    ("x*y/(x - y)", "the line y = x"),
]
RHO = Fr(1, 10 ** 12)


def flt(v):
    return (math.inf if v.sgn > 0 else -math.inf) if isinstance(v, _Inf) else float(v)


def circle_range(text, rho=RHO):
    """exact range of f on the circle of radius rho (as float intervals)"""
    N, D = parse3(text.replace(" ", ""))
    m = max(i + j for (i, j, k) in list(N) + list(D))
    out = []
    for chart in (1, 2):
        Nc, Dc = _polar(N, m, chart), _polar(D, m, chart)
        def at(P):
            q = [Fr(0)] * (max((b for _, b in P), default=0) + 1)
            for (a, b), v in P.items(): q[b] += v * rho ** a
            return ptrim(q)
        n, d = at(Nc), at(Dc)
        if not d: continue
        if not n: out.append((0.0, 0.0)); continue
        if chart == 2:                                  # chart 2 adds only the point (-rho, 0)
            if d[0] != 0: v = float(n[0] / d[0]) if n else 0.0; out.append((v, v))
            continue
        g = pgcd(n, d); n, d = ptrim(pdivmod(n, g)[0]), ptrim(pdivmod(d, g)[0])
        if len(n) <= 1 and len(d) <= 1: v = float(n[0] / d[0]); out.append((v, v))
        else: out += [(flt(a), flt(b)) for a, b in image_1d(n, d)]
    return out


def on_circle(ivs, big=1e9):
    """intervals as arcs of the circle [-oo, +oo] (atan scale), merged; huge values count as oo"""
    def a(v): return math.copysign(math.pi / 2, v) if abs(v) > big else math.atan(v)
    arcs = sorted((a(lo), a(hi)) for lo, hi in ivs)
    merged = []
    for lo, hi in arcs:
        if merged and lo <= merged[-1][1] + 1e-9: merged[-1][1] = max(merged[-1][1], hi)
        else: merged.append([lo, hi])
    return merged


def same(A, B, tol=1e-6):
    return len(A) == len(B) and all(abs(p - q) < tol for x, y in zip(A, B) for p, q in zip(x, y))


def show(arcs):
    def v(t):
        if abs(abs(t) - math.pi / 2) < 1e-9: return "+oo" if t > 0 else "-oo"
        return f"{math.tan(t):.6g}"
    return " u ".join("{" + v(lo) + "}" if abs(hi - lo) < 1e-12 else f"[{v(lo)}, {v(hi)}]" for lo, hi in arcs)


def main():
    passed = 0
    t0 = time.time()
    for f, zero in EXAMPLES:
        r = limit_at_origin(f, ["x", "y"])
        eng = on_circle([tuple(c) for c in r["components"]])
        chk = on_circle(circle_range(f))
        ok = r["proven"] and same(eng, chk)
        passed += ok
        print(f"{f}\n  denominator zero along {zero}\n"
              f"  FMSL (proven={r['proven']}): {r['set']}\n"
              f"  exact range on the circle of radius 1e-12: {show(chk)}\n"
              f"  {'AGREE' if ok else 'DISAGREE'}\n")
    print(f"{passed} of {len(EXAMPLES)} agree ({time.time() - t0:.1f} s)")
    return passed == len(EXAMPLES)


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
