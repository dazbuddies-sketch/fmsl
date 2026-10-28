"""Independent check of fmsl.radial (python benchmarks/radial_crosscheck.py [n]).

fmsl.radial computes the fluid state of f(x, y, r), r = sqrt(x^2 + y^2), from the polar blow-up
and the two-variable engine at special directions. This check uses neither the special directions
nor the engine: on the circle of radius rho (rational), the rational parametrization makes f a
rational function of one variable t, whose exact range over t in R (plus the point t = oo) is
computed by the exact one-variable range algorithm. As rho -> 0 these circle images converge to the
fluid state; they must match at rho = 10^-12 (coverage of the extended line, as in
loop_crosscheck.py, and every isolated value present in both).
"""
import os, random, sys, time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src")); sys.path.insert(0, HERE)
from fmsl.radial import radial_limit, _norm_text, _polar
from fmsl.space import parse3
from fmsl.plane import image_1d
from fmsl.algebraic import pgcd, pdivmod, ptrim
from fmsl.cluster import _Inf
import loop_crosscheck as lc

KNOWN = ["x*y/r", "(r-x)/y", "x^2/(r^2+r^3)", "x/r", "x^2*y/r^3", "(r-x)/r", "x*y/(r*(x+y))",
         "1/r", "x/r^2", "(x^2-y^2)/r^2", "(r^2-x^2)/y^2", "x^3/(r^3+y^3)", "(x+r)/(y+r)", "x*y/(x^2+r^4)"]


def circle_image(text, rho):
    N, D = parse3(_norm_text(text))
    m = max(i + j for (i, j, k) in list(N) + list(D))
    out = []
    for chart in (1, 2):
        Nc, Dc = _polar(N, m, chart), _polar(D, m, chart)
        def at(P):
            deg = max((b for _, b in P), default=0)
            q = [Fr(0)] * (deg + 1)
            for (a, b), v in P.items(): q[b] += v * rho ** a
            return ptrim(q)
        n, d = at(Nc), at(Dc)
        if not d: continue
        if not n: out.append((Fr(0), Fr(0))); continue
        if chart == 2:                                  # only the point t = oo, i.e. u = 0
            if d[0] != 0: v = n[0] / d[0] if n else Fr(0); out.append((v, v))
            continue
        g = pgcd(n, d); n, d = ptrim(pdivmod(n, g)[0]), ptrim(pdivmod(d, g)[0])
        if len(n) <= 1 and len(d) <= 1: v = n[0] / d[0]; out.append((v, v))
        else: out += image_1d(n, d)
    f = lambda v: (float("inf") if v.sgn > 0 else float("-inf")) if isinstance(v, _Inf) else float(v)
    return [(f(a), f(b)) for a, b in out]


def points_agree(A, B, tol=1e-6, big=1e9):
    """every isolated point of A lies in B (within tol), and vice versa; the coverage test alone
       cannot see a single point moving. Values beyond `big` count as +-oo: at radius 1e-12 any blow-up is at least 1/r = 1e12."""
    INF = float("inf")
    def cl(v): return INF if v > big else (-INF if v < -big else v)
    def inside(v, S):
        v = cl(v)
        for lo, hi in S:
            lo, hi = cl(lo), cl(hi)
            if abs(v) == INF:
                if v in (lo, hi): return True
                continue
            pad = tol * max(1, abs(v))
            if lo - pad <= v <= hi + pad: return True
        return False
    pts = lambda S: [a for a, b in S if cl(a) == cl(b)]
    return all(inside(a, B) for a in pts(A)) and all(inside(a, A) for a in pts(B))


def agrees(img, eng):
    return lc.mismatch_line(lc.coverage_line(img), lc.coverage_line(eng)) < 0.002 and points_agree(img, eng)


def rand_text(gen):
    def mono():
        i, j, k = gen.randint(0, 2), gen.randint(0, 2), gen.randint(0, 2)
        return f"({gen.choice([-2, -1, 1, 2])})*x^{i}*y^{j}*r^{k}"
    num = " + ".join(mono() for _ in range(gen.randint(1, 2)))
    den = " + ".join(mono() for _ in range(gen.randint(1, 3)))
    return f"({num})/({den})"


def main(n=60):
    t0 = time.time(); gen = random.Random(17)
    funcs = KNOWN + [rand_text(gen) for _ in range(n)]
    agree = bad = skipped = 0
    for f in funcs:
        try:
            r = radial_limit(f)
        except (ZeroDivisionError, ValueError) as e:
            skipped += 1; continue
        if not r["proven"]: skipped += 1; continue
        eng = [((float("inf") if a.sgn > 0 else float("-inf")) if isinstance(a, _Inf) else float(a),
                (float("inf") if b.sgn > 0 else float("-inf")) if isinstance(b, _Inf) else float(b))
               for a, b in r["components"]]
        img = circle_image(f, Fr(1, 10 ** 12))
        if agrees(img, eng): agree += 1
        else: bad += 1; print("DISAGREE", f, r["set"], [(round(a, 4), round(b, 4)) for a, b in img])
    print(f"{agree + bad} functions compared ({len(KNOWN)} hand-picked): {agree} agree with the exact circle "
          f"image at radius 1e-12, {bad} disagree; {skipped} skipped (not proven or undefined); {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 60)
