"""Independent cross-check of the three-variable engine (fmsl.space).

The fluid state of f = N/D at the origin is the limit, as r -> 0, of the image of f on any small
closed surface around the origin.  Here: the surface of the cube |x|, |y|, |z| <= r.  On each of
its six faces f is a two-variable rational function over a square, and its exact image is
computed by fmsl.plane.image_2d (critical values inside the face, exact one-variable images
on the edges, corners, then one exact test value per gap).  No Newton polyhedron, no weights,
no zooming -- only the shared low-level algebra of fmsl.plane, which is itself tested against
dense sampling.

For r = 10^-3 and 10^-6 the script measures how much of the extended line [-oo, +oo]
(position atan(value)) differs between the cube image and the engine's exact signed fluid
state.  Run:  python benchmarks/cube_crosscheck.py [n_random]
"""
import math, random, sys, time
from fractions import Fraction as Fr

sys.path.insert(0, "src"); sys.path.insert(0, "benchmarks")
from fmsl.space import profile3_polys, parse3
from fmsl.plane import image_2d
from fmsl.cluster import _Inf
from loop_crosscheck import coverage_line, mismatch_line

def face_poly(P, axis, sign, r):
    """P restricted to the face {coordinate[axis] = sign*r}, the other two = r*(c, d)"""
    out = {}
    others = [i for i in range(3) if i != axis]
    for m, v in P.items():
        coef = v * (sign * r) ** m[axis] * r ** (m[others[0]] + m[others[1]])
        key = (m[others[0]], m[others[1]])
        out[key] = out.get(key, 0) + coef
    return {k: v for k, v in out.items() if v != 0}

def _shift(P, c0, d0):
    """P(c0 + u, d0 + v) as a dict in (u, v)"""
    from math import comb
    out = {}
    for (i, j), v in P.items():
        for a in range(i + 1):
            for b in range(j + 1):
                k = (a, b)
                out[k] = out.get(k, 0) + v * comb(i, a) * comb(j, b) * Fr(c0) ** (i - a) * Fr(d0) ** (j - b)
    return {k: v for k, v in out.items() if v != 0}

def base_point_values(P, Q):
    """Points inside the face where P = Q = 0 (after cancelling common factors): the face
       function is 0/0 there, and every value of its two-variable fluid state at that point is
       a limit of values on the face.  image_2d only sees values AT points, so these are added
       here, computed by the two-variable engine (itself checked by benchmarks/loop_crosscheck).
       Returns (intervals, exact); irrational or boundary base points make the result inexact."""
    from fmsl.plane import resultant_d, at_c, bgcd, bdivexact, _poly_roots_in
    from fmsl.multivar import profile_polys
    from fmsl.algebraic import A
    from fmsl.cluster import cmp
    g = bgcd(P, Q)
    if len(g) > 1 or (g and next(iter(g)) != (0, 0)):
        P, Q = bdivexact(P, g), bdivexact(Q, g)
    res = resultant_d(P, Q)
    if not res or all(v == 0 for v in res): return [], False
    ivs, exact = [], True
    for c0 in _poly_roots_in(res, Fr(-1), Fr(1)) if len(res) > 1 else []:
        p, q = at_c(P, c0), at_c(Q, c0)
        if isinstance(c0, A):
            # check whether an irrational c0 carries a common zero at all
            exact = False; continue
        from fmsl.cluster import _real_roots
        for d0 in (_poly_roots_in(q, Fr(-1), Fr(1)) if len(q) > 1 else []):
            if isinstance(d0, A):
                exact = False; continue
            from fmsl.algebraic import peval
            if peval(p, d0) != 0: continue
            if abs(c0) == 1 or abs(d0) == 1:          # on the face boundary: not handled
                exact = False; continue
            cl = profile_polys(_shift(P, c0, d0), _shift(Q, c0, d0))["cluster"]
            if not cl.get("exact", False) or "signed_components" not in cl:
                exact = False; continue
            ivs += [tuple(iv) for iv in cl["signed_components"]]
    return ivs, exact

def cube_image(N, D, r):
    ivs, exact = [], True
    for axis in range(3):
        for sign in (1, -1):
            P, Q = face_poly(N, axis, sign, r), face_poly(D, axis, sign, r)
            if not Q: return None, False
            iv, ex = image_2d(P, Q, box=(Fr(-1), Fr(1)))
            ivs += iv; exact = exact and ex
            iv, ex = base_point_values(P, Q)
            ivs += iv; exact = exact and ex
    return ivs, exact

def to_float(ivs):
    f = lambda v: (math.inf if v.sgn > 0 else -math.inf) if isinstance(v, _Inf) else float(v)
    return [(f(a), f(b)) for a, b in ivs]

CLASSIC = ['x*y*z/(x^2+y^2+z^2)', 'x*y/(x^2+y^2+z^2)', 'x^2*y/(x^4+y^2+z^2)',
           '(x^2+y^2-z^2)/(x^2+y^2+z^2)', 'x^4/((x-y)^2+x^4+z^2)', 'x*y*z/(x^3+y^3+z^3)',
           '1/(x^2+y^2+z^2)', '(x*y+y*z+z*x)/(x^2+y^2+z^2)', 'x^2*y^2*z^2/(x^6+y^6+z^6)',
           'x*y*z/(x^4+y^4+z^2)', '(y^2-x*z)/(x^2+z^2+y^4)']

def _job(f):
    try:
        N, D = parse3(f)
    except ZeroDivisionError:
        return f, None
    t = time.time()
    r = profile3_polys(N, D)
    te = time.time() - t
    if not r["complete"] and r["verdict"] == "UNPROVEN":
        return f, ("UNPROVEN", sorted({e['why'] for e in r['entries'] if e['kind'] == 'unresolved'}))
    ce = coverage_line(to_float(r["signed_components"]))
    ms = []
    for k in KS:
        im, ex = cube_image(N, D, Fr(1, 10 ** k))
        if im is None: return f, None
        ms.append(mismatch_line(coverage_line(to_float(im)), ce))
    return f, (r["signed_text"], ms, te, r["complete"])

KS = [3, 6]

if __name__ == "__main__":
    import multiprocessing as mp
    nrand = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    gen = random.Random(31)
    def rp(n):
        P = {}
        for _ in range(n):
            m = (gen.randint(0, 3), gen.randint(0, 3), gen.randint(0, 3))
            P[m] = P.get(m, 0) + gen.choice([-3, -2, -1, 1, 2, 3])
        return " + ".join(f"({v})*x^{a}*y^{b}*z^{c}" for (a, b, c), v in P.items() if v) or "1"
    funcs = CLASSIC + [f"({rp(gen.randint(1, 3))})/({rp(gen.randint(2, 3))})" for _ in range(nrand)]
    agree = n = unres = 0; T = time.time()
    with mp.Pool(2) as pool:
        for f, res in pool.imap_unordered(_job, funcs):
            if res is None: continue
            if res[0] == "UNPROVEN":
                unres += 1; print(f"{f[:60]:60s} UNPROVEN ({', '.join(res[1])})", flush=True); continue
            text, ms, te, complete = res
            n += 1; good = ms[-1] < 0.002; agree += good
            tag = "" if good else "   <-- CHECK"
            part = "" if complete else "  (engine partial)"
            print(f"{f[:60]:60s} {text[:40]:40s} " + " ".join(f"{m:.4f}" for m in ms) + f"  [{te:.1f}s]{part}{tag}", flush=True)
    print(f"\n{n} functions compared ({unres} left UNPROVEN by the engine), {time.time() - T:.0f}s")
    print(f"cube image matches the exact fluid state at r = 1e-{KS[-1]}: {agree}/{n}")
