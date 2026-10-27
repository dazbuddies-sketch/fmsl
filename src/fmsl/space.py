"""
fmsl.space -- fluid states of rational functions of THREE variables at the origin.

Every value f = N/D approaches near 0 is reached along a Puiseux arc (curve selection lemma).
An arc with x not identically 0 has the form
        x = sigma s^q,   y = c s^p + ...,   z = d s^r + ...      (s -> 0+)
i.e. weights (1, a, b) with a = p/q, b = r/q.  For the weight (1, a, b) the leading behaviour is
        f  ~  s^delta * fN(c, d) / fD(c, d)
where fN, fD are the initial forms (the terms of least weighted order).  The (a, b) quarter-plane
is cut by the lines where two terms tie (inside N, inside D, between N and D) into
    faces     -- the forms are monomials: f -> 0, oo (with signs) or a constant;
    edges     -- the forms depend on one monomial u = c^v2 d^v3: a one-variable problem;
    vertices  -- the forms are genuine polynomials in (c, d): the range of fN/fD over the plane
                 is computed exactly (fmsl.plane).
Roots of the forms are settled by the implicit function / intermediate value theorems when
they are simple, and by zooming (a change of variables, then the same analysis) otherwise.
Arcs inside the coordinate planes are handled by the two-variable engine.

The cluster set is then the closure of the union of all these values; it is connected on the
circle R u {oo} (the same theorem as in two variables), which is checked.
Version 1 limits (reported as UNPROVEN, never guessed): zooming along a whole curve of roots
that is not a coordinate line, zooming at an irrational point, depth > MAX_DEPTH.
"""
from fractions import Fraction as Fr
from math import gcd

from .compiler import Parser
from .gauss import G
from .algebraic import A, ptrim, pmul, peval, sign_of, squarefree
from .cluster import (cmp, POS, NEG, _Inf, _real_roots, _rational_between, _bounds, set_unique,
                      _union_line, _describe_line, union_on_circle, describe_arc)
from .plane import (bclean, bdc, bdd, in_d, at_c, has_real_points, image_2d, image_1d, bgcd,
                    bdivexact, resultant_d, beval, _poly_roots_in, _samples, _gcd_any)
from .multivar import profile_polys

MAX_DEPTH = 4
BUDGET = 4000

# ------------------------------------------------------------------ parsing (x, y, z)
ONE3 = {(0, 0, 0): Fr(1)}
def _add3(P, Q):
    R = dict(P)
    for k, v in Q.items():
        R[k] = R.get(k, 0) + v
        if R[k] == 0: del R[k]
    return R
def _mul3(P, Q):
    R = {}
    for a, u in P.items():
        for b, v in Q.items():
            k = (a[0] + b[0], a[1] + b[1], a[2] + b[2]); R[k] = R.get(k, 0) + u * v
    return {k: v for k, v in R.items() if v != 0}
def _pow3(P, n):
    R = ONE3
    for _ in range(n): R = _mul3(R, P)
    return R

def parse3(text):
    ast = Parser(text).parse()
    def ev(e):
        t = e[0]
        if t == "num":
            v = e[1]
            if isinstance(v, G):
                if not v.is_real(): raise ValueError("real coefficients only")
                v = v.re
            return ({(0, 0, 0): Fr(v)} if v else {}, ONE3)
        if t == "var":
            k = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1)}.get(e[1])
            if k is None: raise NameError(f"unknown variable {e[1]!r} (use x, y, z)")
            return ({k: Fr(1)}, ONE3)
        if t == "neg":
            n, d = ev(e[1]); return ({k: -v for k, v in n.items()}, d)
        if t == "pow":
            n, d = ev(e[1]); k = e[2]
            return (_pow3(n, k), _pow3(d, k)) if k >= 0 else (_pow3(d, -k), _pow3(n, -k))
        if t == "call": raise ValueError("rational functions only")
        (n1, d1), (n2, d2) = ev(e[1]), ev(e[2])
        if t == "+": return (_add3(_mul3(n1, d2), _mul3(n2, d1)), _mul3(d1, d2))
        if t == "-": return (_add3(_mul3(n1, d2), {k: -v for k, v in _mul3(n2, d1).items()}), _mul3(d1, d2))
        if t == "*": return (_mul3(n1, n2), _mul3(d1, d2))
        if t == "/":
            if not n2: raise ZeroDivisionError("division by the zero polynomial")
            return (_mul3(n1, d2), _mul3(d1, n2))
    return ev(ast)

# ------------------------------------------------------------------ the (a, b) arrangement
def _norm_line(al, be, ga):
    g = gcd(gcd(abs(al), abs(be)), abs(ga)) or 1
    al, be, ga = al // g, be // g, ga // g
    if al < 0 or (al == 0 and be < 0): al, be, ga = -al, -be, -ga
    return (al, be, ga)

def _meets_quadrant(al, be, ga):
    if al == 0: return be != 0 and ga * be > 0
    if be == 0: return ga * al > 0
    if al * be < 0: return True
    return ga * al > 0

def _lines(N, D):
    mons = [(m, 0) for m in N] + [(m, 1) for m in D]
    out = set()
    for i in range(len(mons)):
        for j in range(i + 1, len(mons)):
            (m1, _), (m2, _) = mons[i], mons[j]
            al, be, ga = m1[1] - m2[1], m1[2] - m2[2], m2[0] - m1[0]
            if al == 0 and be == 0: continue
            L = _norm_line(al, be, ga)
            if _meets_quadrant(*L): out.add(L)
    return sorted(out)

def _on_line(L, a, b): return L[0] * a + L[1] * b - L[2]

def _line_points(L):
    """parametrisation of the line inside the open quadrant: (param range lo, hi, to_ab)"""
    al, be, ga = L
    if be != 0:
        # b = (ga - al a)/be ; a > 0 and b > 0
        lo, hi = Fr(0), None
        if al != 0:
            t = Fr(ga) / al                                  # where b = 0
            if (be > 0 and al > 0) or (be < 0 and al < 0): hi = t if hi is None else min(hi, t)
            else: lo = max(lo, t)
        to_ab = lambda a, al=al, be=be, ga=ga: (a, Fr(ga - al * a) / be) if not isinstance(a, Fr) else (a, (ga - al * a) / be)
        return lo, hi, (lambda a: (Fr(a), (Fr(ga) - al * Fr(a)) / be)), (lambda a, b: a)
    return Fr(0), None, (lambda b: (Fr(ga) / al, Fr(b))), (lambda a, b: b)

def _in_region(region, a, b):
    if not region: return True
    for val, key in ((a, "a"), (b, "b")):
        lo, hi = region.get(key, (None, None))
        if lo is not None and not val > lo: return False
        if hi is not None and not val < hi: return False
    return True

def _cells(N, D, region=None):
    """sample weights: vertices, one point per edge piece, and points in the faces next to them,
       restricted to the open region (bounds on a and b; their boundary lines take part in the
       arrangement so that faces next to them are sampled).  Returns list of (a, b, kind)."""
    out = _cells_raw(N, D, _region_lines(region))
    return [c for c in out if _in_region(region, c[0], c[1])]

def _region_lines(region):
    extra = []
    for key, (al, be) in (("a", (1, 0)), ("b", (0, 1))):
        for t in (region or {}).get(key, (None, None)):
            if t is not None and t > 0: extra.append((Fr(al), Fr(be), Fr(t)))
    return extra

def _cells_raw(N, D, extra_lines=()):
    Ls = _lines(N, D) + [L for L in extra_lines]
    if not Ls: return [(Fr(1), Fr(1), "face")]
    verts = set()
    for i in range(len(Ls)):
        for j in range(i + 1, len(Ls)):
            (a1, b1, g1), (a2, b2, g2) = Ls[i], Ls[j]
            det = a1 * b2 - a2 * b1
            if det == 0: continue
            a = Fr(g1 * b2 - g2 * b1) / det; b = Fr(a1 * g2 - a2 * g1) / det
            if a > 0 and b > 0: verts.add((a, b))
    out = [(a, b, "vertex") for a, b in sorted(verts)]
    for L in Ls:
        lo, hi, to_ab, par = _line_points(L)
        ps = sorted({par(a, b) for (a, b) in verts if _on_line(L, a, b) == 0})
        pts = []
        if not ps:
            pts.append((lo + hi) / 2 if hi is not None else lo + 1)
        else:
            pts.append((lo + ps[0]) / 2)
            for u, v in zip(ps, ps[1:]): pts.append((u + v) / 2)
            pts.append((ps[-1] + hi) / 2 if hi is not None else ps[-1] + 1)
        for t in pts:
            a, b = to_ab(t)
            if not (a > 0 and b > 0): continue
            out.append((a, b, "edge"))
            # the two faces on either side
            n = (Fr(L[0]), Fr(L[1]))
            eps = None
            for L2 in Ls:
                if L2 == L: continue
                v0 = _on_line(L2, a, b); dv = L2[0] * n[0] + L2[1] * n[1]
                if dv != 0:
                    e = abs(v0 / dv) / 3
                    eps = e if eps is None else min(eps, e)
            for comp, nc in ((a, n[0]), (b, n[1])):
                if nc != 0:
                    e = abs(comp / nc) / 3
                    eps = e if eps is None else min(eps, e)
            eps = eps or Fr(1, 3)
            for sgn in (1, -1):
                a2, b2 = a + sgn * eps * n[0], b + sgn * eps * n[1]
                if a2 > 0 and b2 > 0: out.append((a2, b2, "face"))
    return out

# ------------------------------------------------------------------ initial forms
def _form(P, a, b, sigma):
    if not P: return None, {}
    w = {m: m[0] + a * m[1] + b * m[2] for m in P}
    mm = min(w.values())
    f = {}
    for m, v in P.items():
        if w[m] == mm:
            f[(m[1], m[2])] = f.get((m[1], m[2]), 0) + v * sigma ** m[0]
    return mm, bclean(f)

def _structure(fN, fD):
    exps = list(fN) + list(fD)
    if len(fN) == 1 and len(fD) == 1: return ("mono", None)
    e0 = exps[0]
    dirs = [(e[0] - e0[0], e[1] - e0[1]) for e in exps if e != e0]
    v = next((d for d in dirs if d != (0, 0)), None)
    if v is None: return ("mono", None)
    if all(d[0] * v[1] - d[1] * v[0] == 0 for d in dirs):
        g = gcd(abs(v[0]), abs(v[1])); v = (v[0] // g, v[1] // g)
        if v[0] < 0 or (v[0] == 0 and v[1] < 0): v = (-v[0], -v[1])
        return ("line", v)
    return ("plane", None)

def _to_u(f, v):
    """f(c, d) = c^j0 d^k0 * g(u) with u = c^v0 d^v1 (all exponents of f on one line of
       direction v).  Returns ((j0, k0), g as an ascending coefficient list)."""
    ts = {(j, k): ((j // v[0]) if v[0] else (k // v[1])) for (j, k) in f}
    tmin = min(ts.values())
    e_min = next(e for e, t in ts.items() if t == tmin)
    g = [Fr(0)] * (max(ts.values()) - tmin + 1)
    for e, val in f.items(): g[ts[e] - tmin] += val
    return e_min, ptrim(g)

def _halfline_signs(g, h, positive):
    """signs of g(u)/h(u) for u in (0, oo) (positive=True) or (-oo, 0)"""
    lo, hi = (Fr(0), None) if positive else (None, Fr(0))
    bps = _poly_roots_in(squarefree(pmul(g, h)), lo, hi) if len(ptrim(pmul(g, h))) > 1 else []
    out = set()
    for t in _samples(bps, lo, hi):
        a, b = peval(g, t), peval(h, t)
        if a != 0 and b != 0: out.add(sign_of(a) * sign_of(b))
    return out

def _sign_set_2d(fN, fD):
    """signs of fN/fD over the torus c, d != 0 (exact cell sampling)"""
    F = bclean(_mulb(fN, fD))
    out = set()
    L = in_d(F)
    bps = [Fr(0)]
    if len(L) > 1:
        bps += _poly_roots_in(L[-1])
        disc = resultant_d(F, bdd(F))
        if disc: bps += _poly_roots_in(disc)
    for c0 in _samples(bps):
        if c0 == 0: continue
        g = at_c(F, c0)
        dbps = [Fr(0)] + (_poly_roots_in(g) if len(g) > 1 else [])
        for d0 in _samples(dbps):
            if d0 == 0: continue
            a, b = beval(fN, c0, d0), beval(fD, c0, d0)
            if a != 0 and b != 0: out.add(sign_of(a) * sign_of(b))
    return out

def _mulb(P, Q):
    R = {}
    for (a, b), u in P.items():
        for (c, d), v in Q.items():
            R[(a + c, b + d)] = R.get((a + c, b + d), 0) + u * v
    return R

def _regular_crossing(F, Gc):
    """True if F(c, d) = 0 has a real point with c, d != 0 where F changes sign (simple root in d)
       and Gc != 0 there: then by the implicit function theorem F = 0 is a smooth real curve
       crossing the torus, and Gc keeps its sign nearby."""
    L = in_d(F)
    if len(L) <= 1:
        for c0 in _poly_roots_in(L[0] if L else []):
            if cmp(c0, 0) == 0: continue
            # vertical line c = c0: F changes sign across it if c0 is a simple root
            from .algebraic import pderiv
            if sign_of(peval_any(pderiv(L[0]), c0)) != 0:
                for d0 in (Fr(1), Fr(-1), Fr(2)):
                    if sign_of(beval(Gc, c0, d0)) != 0: return True
        return False
    bps = [Fr(0)] + _poly_roots_in(L[-1])
    disc = resultant_d(F, bdd(F))
    if disc: bps += _poly_roots_in(disc)
    r = resultant_d(F, Gc) if Gc else []
    if r: bps += _poly_roots_in(r)
    for c0 in _samples(bps):
        if c0 == 0: continue
        p = at_c(F, c0)
        g = squarefree(p)
        if len(g) <= 1: continue
        from .algebraic import pderiv
        dp = pderiv(p)
        for d0 in _poly_roots_in(g):
            if cmp(d0, 0) == 0: continue
            # only a SIMPLE root is a sign-changing smooth crossing; a repeated root (e.g. of a
            # square such as (y^2 z - x)^2) touches zero without crossing
            if sign_of(peval_any(dp, d0)) == 0: continue
            if sign_of(beval(Gc, c0, d0)) != 0: return True
    return False

def peval_any(p, x):
    acc = Fr(0)
    for v in reversed(p): acc = acc * x + v
    return acc

def _common_points(A_, B_):
    """real points with c, d != 0 where A_ = B_ = 0, if they are isolated.
       Returns (rational points, flag: irrational or non-isolated common zeros exist)."""
    pts, other = [], False
    if not A_ or not B_: return pts, other
    g = bgcd(A_, B_)
    if any(k != (0, 0) for k in g):
        if has_real_points(g, avoid_zero=True): other = True
        A_, B_ = bdivexact(A_, g), bdivexact(B_, g)
    degA, degB = max(k for _, k in A_), max(k for _, k in B_)
    if degA == 0 or degB == 0:
        U, V = (A_, B_) if degA == 0 else (B_, A_)
        cs = _poly_roots_in(in_d(U)[0])
        for c0 in cs:
            if cmp(c0, 0) == 0: continue
            if isinstance(c0, A): other = True; continue
            gv = at_c(V, c0)
            if not gv: other = True; continue                 # V vanishes on the whole line
            for d0 in _poly_roots_in(gv):
                if cmp(d0, 0) == 0: continue
                if isinstance(d0, A): other = True
                else: pts.append((c0, d0))
        return pts, other
    res = resultant_d(A_, B_)
    if not res: return pts, True
    for c0 in _poly_roots_in(res):
        if cmp(c0, 0) == 0: continue
        gg = _gcd_any(at_c(A_, c0), at_c(B_, c0))
        if len(gg) <= 1: continue
        if isinstance(c0, A) or any(isinstance(v, A) for v in gg): other = True; continue
        for d0 in _poly_roots_in(gg):
            if cmp(d0, 0) == 0: continue
            if isinstance(d0, A): other = True
            else: pts.append((c0, d0))
    return pts, other

def _special_points(F, Gc):
    """real points (c, d != 0) where F = 0 is singular, or meets Gc = 0"""
    pts, other = [], False
    Fc, Fd = bdc(F), bdd(F)
    if Fc and Fd:
        p1, o1 = _common_points(Fc, Fd)
        pts += [q for q in p1 if beval(F, *q) == 0]; other |= o1 and _sing_irrational(F, Fc, Fd)
    elif not Fc or not Fd:
        H = Fd if not Fc else Fc
        if H:
            p1, o1 = _common_points(F, H); pts += p1; other |= o1
        else:
            other = True
    if Gc:
        p2, o2 = _common_points(F, Gc); pts += [q for q in p2 if q not in pts]; other |= o2
    return pts, other

def _sing_irrational(F, Fc, Fd):
    """are there irrational / non-isolated SINGULAR points of F = 0 (not just critical points
       of F off the curve)?  Conservative: True unless every flagged common zero is off F = 0."""
    return True

# ------------------------------------------------------------------ zooming
def _zoom(P, sigma, q, shifts):
    """x = sigma s^q;  for each of y, z: (p, c0) -> s^p (c0 + Y)   or None -> unchanged.
       Returns a polynomial in (s, Y, Z)."""
    out = {}
    from math import comb
    for (i, j, k), v in P.items():
        base = {(q * i, 0, 0): v * sigma ** i}
        for idx, (e, sh) in enumerate(((j, shifts[0]), (k, shifts[1]))):
            if sh is None:
                base = {(m[0], m[1] + (e if idx == 0 else 0), m[2] + (e if idx == 1 else 0)): w for m, w in base.items()}
                continue
            p, c0 = sh
            nb = {}
            for m, w in base.items():
                for t in range(e + 1):
                    coef = comb(e, t) * c0 ** (e - t)
                    if coef == 0: continue
                    key = (m[0] + p * e, m[1] + (t if idx == 0 else 0), m[2] + (t if idx == 1 else 0))
                    nb[key] = nb.get(key, 0) + w * coef
            base = nb
        for m, w in base.items():
            out[m] = out.get(m, 0) + w
    return {m: w for m, w in out.items() if w != 0}

# ------------------------------------------------------------------ the analysis
def _analyze(N, D, depth, sigmas, label, out, budget, region=None, planes=True):
    vertex_jobs = []
    ctx = {"region": region, "planes": planes, "sigmas": sigmas}
    for (a, b, kind) in _cells(N, D, region):
        for sigma in sigmas:
            if budget[0] <= 0:
                out.append({"kind": "unresolved", "why": "work budget", "where": label}); return vertex_jobs
            budget[0] -= 1
            mN, fN = _form(N, a, b, sigma)
            mD, fD = _form(D, a, b, sigma)
            if not fD: continue
            where = f"{label}weight (1, {a}, {b}) x{'>' if sigma > 0 else '<'}0"
            if not fN:
                out.append({"kind": "value", "value": Fr(0), "where": where}); continue
            delta = mN - mD
            st, v = _structure(fN, fD)
            if st == "mono":
                (jn, kn), n = next(iter(fN.items())); (jd, kd), dd = next(iter(fD.items()))
                r = n / dd
                if delta > 0: out.append({"kind": "value", "value": Fr(0), "where": where})
                elif delta < 0:
                    odd = (jn - jd) % 2 or (kn - kd) % 2
                    out.append({"kind": "inf", "signs": {1, -1} if odd else {sign_of(r)}, "where": where})
                elif (jn, kn) == (jd, kd):
                    out.append({"kind": "value", "value": r, "where": where})
                else:
                    odd = (jn - jd) % 2 or (kn - kd) % 2
                    iv = (NEG, POS) if odd else ((Fr(0), POS) if r > 0 else (NEG, Fr(0)))
                    out.append({"kind": "interval", "iv": [iv], "where": where + " (monomial ratio)"})
                continue
            if st == "line":
                _edge(N, D, fN, fD, v, a, b, sigma, delta, depth, label, where, out, budget, ctx)
                continue
            # vertex: genuine two-variable forms
            if delta > 0:
                out.append({"kind": "value", "value": Fr(0), "where": where})
                if _regular_crossing(fD, fN):
                    out.append({"kind": "all", "where": where, "why": "a smooth real pole surface crosses here while f is near 0"})
                    continue
                if _binomial_zero_set_zoom(N, D, fD, a, b, sigma, depth, label, where, out, budget, ctx):
                    continue
                _vertex_zoom(N, D, fD, fN, a, b, sigma, depth, label, where, out, budget, zero_of="D", ctx=ctx)
            elif delta < 0:
                sg = _sign_set_2d(fN, fD)
                out.append({"kind": "inf", "signs": sg or {1, -1}, "where": where})
                if _regular_crossing(fN, fD):
                    out.append({"kind": "all", "where": where, "why": "f changes sign through 0 here while it is unbounded"})
                    continue
                _vertex_zoom(N, D, fN, fD, a, b, sigma, depth, label, where, out, budget, zero_of="N", ctx=ctx)
            else:
                vertex_jobs.append((fN, fD, where))
                _vertex_zoom(N, D, fD, fN, a, b, sigma, depth, label, where, out, budget, zero_of="D", only_special=True, ctx=ctx)
    return vertex_jobs

def _edge(N, D, fN, fD, v, a, b, sigma, delta, depth, label, where, out, budget, ctx=None):
    eN, gN = _to_u(fN, v); eD, gD = _to_u(fD, v)
    # domain of u = c^v0 d^v1 over c, d != 0
    pos_only = (v[0] % 2 == 0) and (v[1] % 2 == 0)
    halves = (True,) if pos_only else (True, False)
    pj, pk = eN[0] - eD[0], eN[1] - eD[1]
    if delta > 0:
        out.append({"kind": "value", "value": Fr(0), "where": where})
    elif delta < 0:
        sg = set()
        for (sc, sd) in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
            upos = (sc ** v[0]) * (sd ** v[1]) > 0
            pre = (sc ** (pj % 2)) * (sd ** (pk % 2))
            sg |= {pre * t for t in _halfline_signs(gN, gD, upos)}
        out.append({"kind": "inf", "signs": sg or {1, -1}, "where": where})
    else:
        # f -> u^m gN(u)/gD(u)
        if v[0]: m = Fr(pj, v[0])
        else: m = Fr(pk, v[1])
        if m.denominator != 1 or (pj, pk) != (int(m) * v[0], int(m) * v[1]):
            out.append({"kind": "unresolved", "why": "edge prefactor", "where": where}); return
        m = int(m)
        n_, d_ = list(gN), list(gD)
        if m > 0: n_ = [Fr(0)] * m + n_
        elif m < 0: d_ = [Fr(0)] * (-m) + d_
        ivs = []
        for upos in halves:
            ivs += image_1d(n_, d_, Fr(0), None) if upos else image_1d(n_, d_, None, Fr(0))
        out.append({"kind": "interval", "iv": ivs, "where": where + f" (depends on u = c^{v[0]} d^{v[1]})"})
    # roots of the forms along the edge
    roots = []
    for g, other, is_den in ((gD, gN, True), (gN, gD, False)):
        if not is_den and delta >= 0: continue
        for u0 in _poly_roots_in(g):
            if cmp(u0, 0) == 0: continue
            if pos_only and cmp(u0, 0) < 0: continue
            from .algebraic import pderiv
            simple = sign_of(peval_any(pderiv(g), u0)) != 0
            other_nz = sign_of(peval_any(other, u0)) != 0
            if simple and other_nz:
                if is_den:
                    out.append({"kind": "inf", "signs": {1, -1}, "where": where + " pole surface (implicit function theorem)"})
                    if delta > 0: out.append({"kind": "all", "where": where, "why": "f crosses through infinity here while it is near 0"})
                else:
                    out.append({"kind": "all", "where": where, "why": "f changes sign through 0 here while it is unbounded"})
                continue
            roots.append(u0)
    for u0 in set_unique(roots):
        if v[1] == 0 and v[0] == 1:
            _do_zoom(N, D, a, b, sigma, ("y", u0), depth, label, where, out, budget, ctx)
        elif v[0] == 0 and v[1] == 1:
            _do_zoom(N, D, a, b, sigma, ("z", u0), depth, label, where, out, budget, ctx)
        elif abs(v[0]) == 1 or abs(v[1]) == 1:
            _curve_zoom(N, D, a, b, sigma, v, u0, depth, label, where, out, budget, ctx)
        else:
            out.append({"kind": "unresolved", "why": f"zoom along the curve c^{v[0]} d^{v[1]} = {u0}", "where": where})

def _vertex_zoom(N, D, F, Gc, a, b, sigma, depth, label, where, out, budget, zero_of, only_special=False, ctx=None):
    """zoom at the special real points of F = 0 (singular, or shared with Gc)"""
    if not has_real_points(F, avoid_zero=True) and not only_special:
        return
    pts, other = _special_points(F, Gc)
    # a non-isolated real zero set without regular crossings (e.g. a double curve) needs a
    # family zoom, which version 1 does not do
    if not _regular_crossing(F, {(0, 0): Fr(1)}) and has_real_points(F, avoid_zero=True) and not pts:
        other = True
    if other:
        out.append({"kind": "unresolved", "why": "zoom at irrational or non-isolated special points", "where": where})
    for (c0, d0) in pts:
        _do_zoom(N, D, a, b, sigma, ("yz", (c0, d0)), depth, label, where, out, budget, ctx)

def _do_zoom(N, D, a, b, sigma, what, depth, label, where, out, budget, ctx=None):
    if depth >= MAX_DEPTH:
        out.append({"kind": "unresolved", "why": "depth limit", "where": where}); return
    q = a.denominator * b.denominator // gcd(a.denominator, b.denominator)
    p, r = int(a * q), int(b * q)
    if what[0] == "y":
        if isinstance(what[1], A):
            out.append({"kind": "unresolved", "why": "irrational zoom point", "where": where}); return
        sh = ((p, what[1]), None); txt = f"y = {what[1]}*|x|^{a}+..."
    elif what[0] == "z":
        if isinstance(what[1], A):
            out.append({"kind": "unresolved", "why": "irrational zoom point", "where": where}); return
        sh = (None, (r, what[1])); txt = f"z = {what[1]}*|x|^{b}+..."
    else:
        c0, d0 = what[1]
        sh = ((p, c0), (r, d0)); txt = f"y = {c0}*|x|^{a}+..., z = {d0}*|x|^{b}+..."
    N2, D2 = _zoom(N, sigma, q, sh), _zoom(D, sigma, q, sh)
    lab = f"{label}[{txt} x{'>' if sigma > 0 else '<'}0] "
    # the new first coordinate is s with x = sigma s^q: weights of a kept coordinate scale by q
    region = (ctx or {}).get("region") or {}
    new_region = {}
    if what[0] == "y" and "b" in region:
        new_region["b"] = tuple(None if t is None else t * q for t in region["b"])
    if what[0] == "z" and "a" in region:
        new_region["a"] = tuple(None if t is None else t * q for t in region["a"])
    planes = (ctx or {}).get("planes", True)
    jobs = _analyze(N2, D2, depth + 1, (1,), lab, out, budget, new_region or None, planes)
    if planes: _planes(N2, D2, depth + 1, lab, out)
    _finish_vertices(jobs, out)

def _binomial_zero_set_zoom(N, D, fD, a, b, sigma, depth, label, where, out, budget, ctx):
    """Vertex with f -> 0 away from fD = 0, where fD alone depends only on u = c^v0 d^v1
       (e.g. fD = (c^2 d - 1)^2).  Then fD = 0 on the torus is the union of the curves
       u = u0 over the real roots u0 of gD, and zooming along each one covers every path that
       stays near fD = 0, exactly as for an edge.  Returns False if this does not apply."""
    if len(fD) < 2: return False
    st, v = _structure({}, fD)
    if st != "line": return False
    if not (abs(v[0]) == 1 or abs(v[1]) == 1): return False
    _, gD = _to_u(fD, v)
    pos_only = (v[0] % 2 == 0) and (v[1] % 2 == 0)
    roots = [u0 for u0 in _poly_roots_in(gD) if cmp(u0, 0) != 0 and not (pos_only and cmp(u0, 0) < 0)]
    for u0 in set_unique(roots):
        if v == (1, 0):
            _do_zoom(N, D, a, b, sigma, ("y", u0), depth, label, where, out, budget, ctx)
        elif v == (0, 1):
            _do_zoom(N, D, a, b, sigma, ("z", u0), depth, label, where, out, budget, ctx)
        else:
            _curve_zoom(N, D, a, b, sigma, v, u0, depth, label, where, out, budget, ctx)
    return True

def _curve_zoom(N, D, a, b, sigma, v, u0, depth, label, where, out, budget, ctx):
    """zoom along the curve c^v0 d^v1 = u0 of an edge, where one exponent is +-1.
       W = x^k y^v0 z^v1 (k = -(v0 a + v1 b)) has weight 0 along the edge and tends to
       w0 = sigma^k u0.  If v0 = +-1, y = (W x^-k z^-v1)^v0 exactly (x, z != 0): the new
       coordinates (x, W - w0, z) turn the curve into a coordinate point.  A path in the new
       coordinates is a genuine path into the origin only if y -> 0 too, i.e. the old weight
       of y, -(k + v1 b')/v0 ... is positive: this is passed on as a bound on the weight of z."""
    if depth >= MAX_DEPTH:
        out.append({"kind": "unresolved", "why": "depth limit", "where": where}); return
    if isinstance(u0, A):
        out.append({"kind": "unresolved", "why": "irrational zoom point", "where": where}); return
    k = -(v[0] * a + v[1] * b)
    if k.denominator != 1:
        out.append({"kind": "unresolved", "why": "curve zoom: non-integral exponent", "where": where}); return
    k = int(k)
    w0 = Fr(sigma) ** k * u0
    # solve for the coordinate with exponent +-1
    if abs(v[0]) == 1: sol, oth, e_sol, e_oth, kept_key = 1, 2, v[0], v[1], "b"
    else:              sol, oth, e_sol, e_oth, kept_key = 2, 1, v[1], v[0], "a"
    # y_sol = (W x^-k y_oth^-e_oth)^(1/e_sol) with e_sol = +-1  ->  exponents of W, x, y_oth
    sgn = e_sol                                    # 1/e_sol = e_sol for +-1
    def sub(P):
        R = {}
        for m, val in P.items():
            j = m[sol]                              # power of the eliminated coordinate
            ex = [m[0] - k * sgn * j, 0, 0]
            ex[oth] = m[oth] - e_oth * sgn * j
            wpow = sgn * j                          # power of W (may be negative)
            key = (ex[0], wpow, ex[oth])
            R[key] = R.get(key, 0) + val
        # clear negative powers of x, W, the other coordinate (same factor for N and D)
        return R
    Ns, Ds = sub(N), sub(D)
    shift = [max(0, -min(key[i] for key in list(Ns) + list(Ds))) for i in range(3)]
    def finish(R):
        out_ = {}
        from math import comb
        for (ix, iw, io), val in R.items():
            ix += shift[0]; iw += shift[1]; io += shift[2]
            # W = w0 + Y
            for t in range(iw + 1):
                c = comb(iw, t) * w0 ** (iw - t)
                if c == 0: continue
                key = [ix, 0, 0]; key[sol] = t; key[oth] = io
                key = tuple(key)
                out_[key] = out_.get(key, 0) + val * c
        return {kk: vv for kk, vv in out_.items() if vv != 0}
    N2, D2 = finish(Ns), finish(Ds)
    if not D2:
        out.append({"kind": "unresolved", "why": "curve zoom: denominator vanishes on the curve", "where": where}); return
    # genuine-path condition: old weight of the eliminated coordinate > 0
    #   weight(y_sol) = sgn * (weight(W) - k - e_oth * w_oth) with weight(W) = 0 along W -> w0
    #   => sgn * (-k - e_oth * w) > 0, a bound on w = weight of the other coordinate
    old_region = (ctx or {}).get("region") or {}
    elim_key = "a" if sol == 1 else "b"
    lo, hi = old_region.get(kept_key, (None, None))
    # old weight of the eliminated coordinate as a function of w (weight of the kept one):
    #   L(w) = sgn * (-k - e_oth * w);  it must lie in (max(0, lo_e), hi_e)
    lo_e, hi_e = old_region.get(elim_key, (None, None))
    lo_e = Fr(0) if lo_e is None else max(Fr(0), lo_e)
    cons = [(-sgn * e_oth, -sgn * k - lo_e)]                 # L(w) - lo_e > 0
    if hi_e is not None: cons.append((sgn * e_oth, hi_e + sgn * k))   # hi_e - L(w) > 0
    for coef, const in cons:                                # coef*w + const > 0
        if coef > 0:
            t = Fr(-const) / coef; lo = t if lo is None else max(lo, t)
        elif coef < 0:
            t = Fr(-const) / coef; hi = t if hi is None else min(hi, t)
        elif const <= 0:
            return                                          # no genuine paths on this curve
    if lo is not None and hi is not None and lo >= hi: return
    if hi is not None and hi <= 0: return
    region = {kept_key: (lo, hi)}
    lab = f"{label}[c^{v[0]} d^{v[1]} = {u0} x{'>' if sigma > 0 else '<'}0] "
    jobs = _analyze(N2, D2, depth + 1, (ctx or {}).get("sigmas", (1, -1)), lab, out, budget, region, planes=False)
    _finish_vertices(jobs, out)

def _planes(N, D, depth, label, out):
    """arcs inside coordinate planes, by the two-variable engine"""
    specs = [((1, 2), "x = 0"), ((0, 2), "y = 0"), ((0, 1), "z = 0")] if depth == 0 else [((0, 2), "Y = 0"), ((0, 1), "Z = 0")]
    for keep, name in specs:
        drop = ({0, 1, 2} - set(keep)).pop()
        N2 = {(m[keep[0]], m[keep[1]]): v for m, v in N.items() if m[drop] == 0}
        D2 = {(m[keep[0]], m[keep[1]]): v for m, v in D.items() if m[drop] == 0}
        if not D2: continue
        if not N2:
            out.append({"kind": "value", "value": Fr(0), "where": label + name}); continue
        sides = (1, -1) if depth == 0 else (1,)
        pr = profile_polys(N2, D2, sides=sides, y_axis=(depth == 0))
        cl = pr["cluster"]
        if not cl.get("complete", False) or "signed_components" not in cl:
            out.append({"kind": "unresolved", "why": "two-variable plane not resolved", "where": label + name})
            if "signed_components" not in cl: continue
        out.append({"kind": "interval", "iv": cl["signed_components"], "where": label + name + " (two-variable engine)",
                    "exact": cl.get("exact", True)})

def _finish_vertices(jobs, out):
    """vertices with delta = 0: the exact range of fN/fD over the plane; the values already
       found elsewhere are candidates for its endpoints (asymptotic directions)"""
    known = []
    for e in out:
        if e["kind"] == "value": known.append(e["value"])
        elif e["kind"] == "interval":
            for lo, hi in e["iv"]:
                for x in (lo, hi):
                    if not isinstance(x, _Inf): known.append(x)
    for fN, fD, where in jobs:
        try:
            ivs, ex = image_2d(fN, fD, extra=set_unique(known))
            out.append({"kind": "interval", "iv": ivs, "where": where + " (two-variable form, exact range)", "exact": ex})
        except (ArithmeticError, NotImplementedError) as exc:
            out.append({"kind": "unresolved", "why": f"vertex range: {exc}", "where": where})

def profile3(text):
    N, D = parse3(text)
    return profile3_polys(N, D)

def _cancel_monomial_content(N, D):
    """divide N and D by their common monomial factor x^a y^b z^c (f is unchanged where defined;
       a common factor that vanishes on a coordinate plane would otherwise look like a whole
       surface of 0/0 points)"""
    if not N or not D: return N, D
    common = tuple(min(m[i] for m in list(N) + list(D)) for i in range(3))
    if common == (0, 0, 0): return N, D
    sub = lambda P: {tuple(m[i] - common[i] for i in range(3)): v for m, v in P.items()}
    return sub(N), sub(D)

def profile3_polys(N, D):
    N, D = _cancel_monomial_content(N, D)
    out, budget = [], [BUDGET]
    jobs = _analyze(N, D, 0, (1, -1), "", out, budget)
    _planes(N, D, 0, "", out)
    _finish_vertices(jobs, out)
    return {"N": N, "D": D, "entries": out, **_summarise(out)}

def _summarise(out):
    items, complete, exact = [], True, True
    for e in out:
        k = e["kind"]
        if k == "value": items.append((e["value"], e["value"]))
        elif k == "inf":
            for t in e["signs"]: items.append((POS, POS) if t > 0 else (NEG, NEG))
        elif k == "all": items.append("all")
        elif k == "interval":
            items += list(e["iv"]); exact = exact and e.get("exact", True)
        elif k == "unresolved": complete = False
    signed = _union_line(items) if items else []
    def to_circle(it):
        if it != "all" and isinstance(it[0], _Inf) and it[0].sgn < 0 and isinstance(it[1], _Inf) and it[1].sgn < 0:
            return (POS, POS)
        return it
    comps, connected, whole = union_on_circle([to_circle(i) for i in items]) if items else ([], True, False)
    # verdict
    pts = [c for c in signed]
    finite_pts = [c for c in pts if not isinstance(c[0], _Inf) and cmp(c[0], c[1]) == 0]
    nondeg = [c for c in pts if cmp(c[0], c[1]) != 0]
    if nondeg or len(pts) > 1 and not all(isinstance(c[0], _Inf) for c in pts):
        verdict = "NO LIMIT (proven)"
    elif not complete:
        verdict = "UNPROVEN"
    elif len(pts) == 1 and finite_pts:
        verdict = f"LIMIT = {finite_pts[0][0]} (proven)"
    elif pts and all(isinstance(c[0], _Inf) for c in pts):
        verdict = "DIVERGES to infinity (proven)"
    else:
        verdict = "UNPROVEN"
    return {"verdict": verdict, "signed_components": signed, "signed_text": _describe_line(signed),
            "components": comps, "connected": connected, "whole": whole,
            "text": describe_arc(comps, connected, whole), "complete": complete, "exact": exact}

def describe3(text):
    r = profile3(text)
    lines = [f"f(x, y, z) = {text}", f"  verdict: {r['verdict']}",
             f"  fluid state (signed): {r['signed_text']}",
             f"  on the circle R u {{oo}}: {r['text']}" + ("" if r["connected"] else "  (NOT connected: a path family is missing)")]
    if not r["complete"]:
        why = sorted({e['why'] for e in r['entries'] if e['kind'] == 'unresolved'})
        lines.append(f"  partial: {', '.join(why)}")
    return "\n".join(lines)
