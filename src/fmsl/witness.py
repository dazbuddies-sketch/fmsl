"""Proven limit values in any number of variables, each with an explicit curve (an inner bound).

For f = P/Q in n variables, every weight vector w (positive integers) and every real point c give
the curve x(t) = (c1 t^w1, ..., cn t^wn), t -> 0+. Writing in_w P for the terms of least weighted
degree d_P (and the same for Q),

    P(x(t)) = t^d_P (in_w P(c) + O(t)),     Q(x(t)) = t^d_Q (in_w Q(c) + O(t)),

so when in_w Q(c) != 0 the limit of f along the curve is exactly
    in_w P(c) / in_w Q(c)          if d_P = d_Q,
    0                               if d_P > d_Q,
    +oo or -oo (sign of in_w P(c)/in_w Q(c))   if d_P < d_Q and in_w P(c) != 0.
Every such value is a limit value of f at the origin, proven by its curve.

Whole intervals: along a segment c(s) = (1 - s) a + s b, s in [0, 1], the map s -> [in_w P : in_w Q]
on the circle R u {oo} is continuous when in_w P and in_w Q have no common root on the segment
(decided exactly by a gcd and Sturm's theorem). Between consecutive real roots of in_w Q it is a
continuous real function, and at a root it tends to +oo or -oo with a sign read off exactly. Each
piece therefore covers a whole closed interval of values (intermediate value theorem), and each of
those values is a limit value of f (the closure of limit values is the set of limit values).

This is the n-variable form of the weighted initial forms that fmsl.space uses exactly in three
variables. Here it is only an inner bound: every value reported is proven, but values may be
missing (weights beyond the search bound, curves needing more than one term, zooming). It is the
exact answer whenever it covers the whole circle, and it proves "the limit does not exist" whenever
it contains two different values. No other tool computes either in four or more variables.

    from fmsl.witness import inner_bound
    r = inner_bound("x*w/(x^2+y^2+z^2+w^2)", ["x", "y", "z", "w"])
    r["set"], r["whole_circle"], r["limit_exists"]   # '[-1/2, 1/2]', False, False
"""
import random
from fractions import Fraction as Fr
from itertools import product
from math import gcd, inf, isinf

from .compiler import Parser
from .gauss import G

__all__ = ["inner_bound", "parse_n", "curve_point"]


# ------------------------------------------------------------------ polynomials in n variables
def _add(P, Q):
    R = dict(P)
    for k, v in Q.items():
        R[k] = R.get(k, 0) + v
        if R[k] == 0: del R[k]
    return R


def _mul(P, Q):
    R = {}
    for a, u in P.items():
        for b, v in Q.items():
            k = tuple(i + j for i, j in zip(a, b)); R[k] = R.get(k, 0) + u * v
    return {k: v for k, v in R.items() if v != 0}


def parse_n(text, names):
    """(P, Q) as dicts {exponent tuple: Fraction} with f = P/Q"""
    n = len(names)
    one = {(0,) * n: Fr(1)}
    def pw(P, k):
        R = one
        for _ in range(k): R = _mul(R, P)
        return R
    def ev(e):
        t = e[0]
        if t == "num":
            v = e[1]
            if isinstance(v, G):
                if not v.is_real(): raise ValueError("real coefficients only")
                v = v.re
            return ({(0,) * n: Fr(v)} if v else {}, one)
        if t == "var":
            if e[1] not in names: raise NameError(f"unknown variable {e[1]!r}")
            i = names.index(e[1]); return ({tuple(int(j == i) for j in range(n)): Fr(1)}, one)
        if t == "neg":
            a, b = ev(e[1]); return ({k: -v for k, v in a.items()}, b)
        if t == "pow":
            a, b = ev(e[1]); k = e[2]
            if not isinstance(k, int): raise ValueError("whole-number powers only")
            return (pw(a, k), pw(b, k)) if k >= 0 else (pw(b, -k), pw(a, -k))
        if t == "call": raise ValueError("rational functions only")
        (a1, b1), (a2, b2) = ev(e[1]), ev(e[2])
        if t == "+": return (_add(_mul(a1, b2), _mul(a2, b1)), _mul(b1, b2))
        if t == "-": return (_add(_mul(a1, b2), {k: -v for k, v in _mul(a2, b1).items()}), _mul(b1, b2))
        if t == "*": return (_mul(a1, a2), _mul(b1, b2))
        if t == "/":
            if not a2: raise ZeroDivisionError("division by the zero polynomial")
            return (_mul(a1, b2), _mul(b1, a2))
    return ev(Parser(text).parse())


def _initial(P, w):
    d = min(sum(i * j for i, j in zip(k, w)) for k in P)
    return d, {k: v for k, v in P.items() if sum(i * j for i, j in zip(k, w)) == d}


def _ev(P, c):
    s = Fr(0)
    for k, v in P.items():
        m = v
        for ci, e in zip(c, k):
            if e: m *= ci ** e
        s += m
    return s


# ------------------------------------------------------------------ one variable, exact
def _utrim(p):
    p = list(p)
    while p and p[-1] == 0: p.pop()
    return p


def _uadd(p, q):
    n = max(len(p), len(q))
    return _utrim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)])


def _umul(p, q):
    if not p or not q: return []
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): r[i + j] += a * b
    return _utrim(r)


def _urem(p, q):
    p = list(p)
    while len(p) >= len(q) and p:
        f = p[-1] / q[-1]; s = len(p) - len(q)
        for i, b in enumerate(q): p[s + i] -= f * b
        p = _utrim(p)
    return p


def _ugcd(p, q):
    while q: p, q = q, _urem(p, q)
    return [a / p[-1] for a in p] if p else p


def _uderiv(p): return _utrim([i * a for i, a in enumerate(p)][1:])


def _ueval(p, x):
    s = Fr(0)
    for a in reversed(p): s = s * x + a
    return s


def _sturm(p):
    seq = [p, _uderiv(p)]
    while seq[-1]:
        r = _urem(seq[-2], seq[-1])
        if not r: break
        seq.append([-a for a in r])
    return seq


def _changes(seq, x):
    vals = [v for v in (_ueval(q, x) for q in seq) if v != 0]
    return sum(1 for a, b in zip(vals, vals[1:]) if (a > 0) != (b > 0))


def _roots(p, lo, hi):
    """disjoint rational intervals (a, b], one per distinct real root of p in (lo, hi), sorted"""
    p = _utrim(p)
    if len(p) <= 1: return []
    g = _ugcd(p, _uderiv(p))
    if len(g) > 1:                                  # squarefree part
        q, r = p, []
        # p / g by long division
        q = list(p); out = [Fr(0)] * (len(p) - len(g) + 1)
        while len(q) >= len(g) and q:
            f = q[-1] / g[-1]; s = len(q) - len(g); out[s] = f
            for i, b in enumerate(g): q[s + i] -= f * b
            q = _utrim(q)
        p = _utrim(out)
    seq = _sturm(p)
    out, stack = [], [(lo, hi)]
    while stack:
        a, b = stack.pop()
        k = _changes(seq, a) - _changes(seq, b)
        if k == 0: continue
        if k == 1 and b - a < Fr(1, 10 ** 6): out.append((a, b)); continue
        m = (a + b) / 2
        if _ueval(p, m) == 0:
            out.append((m - Fr(1, 10 ** 9) * (b - a), m)); stack += [(a, m - Fr(1, 10 ** 9) * (b - a)), (m, b)]
            continue
        stack += [(a, m), (m, b)]
    return sorted(out)


def _along(P, a, b):
    """P(c(s)) for c(s) = a + s (b - a), as a polynomial in s"""
    lin = [[ai, bi - ai] for ai, bi in zip(a, b)]
    s = []
    for k, v in P.items():
        m = [v]
        for L, e in zip(lin, k):
            for _ in range(e): m = _umul(m, _utrim(L))
        s = _uadd(s, m)
    return s


# ------------------------------------------------------------------ the inner bound
def _segment(pP, pQ):
    """the values of [pP : pQ] covered on s in [0, 1], as closed intervals of the extended line,
       or None if pP and pQ share a root there (the path is not continuous on the circle)"""
    if len(_utrim(pQ)) == 0: return None
    g = _ugcd(pP, pQ) if pP else pQ
    if len(g) > 1 and (_roots(g, Fr(-1, 10 ** 9), Fr(1)) or _ueval(g, 0) == 0): return None
    roots = _roots(pQ, Fr(0), Fr(1))
    if _ueval(pQ, 0) == 0 or _ueval(pQ, 1) == 0: return None
    cuts = [Fr(0)] + [x for iv in roots for x in iv] + [Fr(1)]
    val = lambda s: _ueval(pP, s) / _ueval(pQ, s)
    out = []
    for i in range(0, len(cuts), 2):                # pieces [cut_i, cut_{i+1}] free of roots of pQ
        lo, hi = cuts[i], cuts[i + 1]
        left = val(lo) if i == 0 else None
        right = val(hi) if i + 2 == len(cuts) else None
        mid = val((lo + hi) / 2)
        sgnQ = 1 if _ueval(pQ, (lo + hi) / 2) > 0 else -1
        def at_root(iv):                            # the signed infinity approached next to a root
            a, b = iv
            while _roots(pP, a, b) or _ueval(pP, b) == 0:
                m = (a + b) / 2; a, b = (a, m) if _roots(pQ, a, m) else (m, b)
            return inf if (1 if _ueval(pP, b) > 0 else -1) * sgnQ > 0 else -inf
        if left is None: left = at_root((cuts[i - 1], cuts[i]))
        if right is None: right = at_root((cuts[i + 1], cuts[i + 2]))
        for u, v in ((left, mid), (mid, right)):
            out.append((min(u, v), max(u, v)))
    return out


def _climb(iP, iQ, c, sgn, gen, steps=300):
    """a point where sgn * in_w P / in_w Q is (locally) large: a float search, returned rational.
       It only chooses where to look; the value there is then computed exactly."""
    fl = lambda P: [(k, float(v)) for k, v in P.items()]
    fP, fQ = fl(iP), fl(iQ)
    def ev(P, x):
        s = 0.0
        for k, v in P:
            m = v
            for xi, e in zip(x, k): m *= xi ** e
            s += m
        return s
    def val(x):
        q = ev(fQ, x)
        return sgn * ev(fP, x) / q if q else -1e300
    x = [float(ci) for ci in c]; best = val(x); h = 0.5
    for _ in range(steps):
        y = [xi + h * gen.uniform(-1, 1) for xi in x]
        if max(abs(t) for t in y) > 1e3: continue
        v = val(y)
        if v > best and abs(v) < 1e12: x, best = y, v
        else: h *= 0.98
    return tuple(Fr(t).limit_denominator(10 ** 6) for t in x)


def _union(ivs):
    out = []
    for lo, hi in sorted(ivs):
        if out and lo <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], hi))
        else: out.append((lo, hi))
    return out


def _covers_circle(U):
    return len(U) == 1 and U[0] == (-inf, inf) or (len(U) == 2 and U[0][0] == -inf and U[1][1] == inf
                                                       and U[0][1] >= U[1][0])


def _weights(n, k):
    for w in product(range(1, k + 1), repeat=n):
        g = 0
        for x in w: g = gcd(g, x)
        if g == 1: yield w


def _sub(P, w, c0, j, n):
    """P after x_j = c0_j t^w_j and x_i = t^w_i (c0_i + y_i) (i != j); new variables (t, y_i...)"""
    def var(k): return {tuple(int(i == k) for i in range(n)): Fr(1)}
    others = [i for i in range(n) if i != j]
    X = []
    for i in range(n):
        tw = {tuple([w[i]] + [0] * (n - 1)): Fr(1)}
        if i == j: X.append({k: v * c0[i] for k, v in tw.items()})
        else: X.append(_mul(tw, _add({(0,) * n: c0[i]} if c0[i] else {}, var(1 + others.index(i)))))
    cache = {}
    def pw(i, e):
        if (i, e) not in cache: cache[(i, e)] = {(0,) * n: Fr(1)} if e == 0 else _mul(pw(i, e - 1), X[i])
        return cache[(i, e)]
    out = {}
    for k, v in P.items():
        m = {(0,) * n: v}
        for i, e in enumerate(k):
            if e: m = _mul(m, pw(i, e))
        out = _add(out, m)
    return out


def curve_point(witness, s):
    """the point at parameter s > 0 on the curve that proves a witnessed value (exact)"""
    _, w, c, chain = witness if len(witness) == 4 else (*witness, [])
    x = [ci * s ** wi for ci, wi in zip(c, w)]
    for w0, c0, j in reversed(chain):               # undo the zooms, innermost first
        n = len(c0); others = [i for i in range(n) if i != j]
        t, ys = x[0], x[1:]
        x = [c0[i] * t ** w0[i] if i == j else t ** w0[i] * (c0[i] + ys[others.index(i)]) for i in range(n)]
    return tuple(x)


def _zeros(iQ, n, positive, gen, tries=40):
    """rational points c with in_w Q(c) = 0 (and c_p > 0 for p in positive): grid search and
       solving for a variable of degree 1"""
    grid = [Fr(k) for k in (-3, -2, -1, 1, 2, 3)] + [Fr(1, 2), Fr(-1, 2)]
    pick = lambda i: abs(gen.choice(grid)) if i in positive else gen.choice(grid)
    found = []
    for _ in range(tries):
        c = tuple(pick(i) for i in range(n))
        if _ev(iQ, c) == 0: found.append(c)
    for k in range(n):
        if max(m[k] for m in iQ) != 1: continue
        A = {m: v for m, v in iQ.items() if m[k] == 1}; B = {m: v for m, v in iQ.items() if m[k] == 0}
        for _ in range(4):
            c = [pick(i) for i in range(n)]; c[k] = Fr(1)
            a = _ev(A, c)
            if a == 0: continue
            c[k] = -_ev(B, c) / a
            if c[k] == 0 or (k in positive and c[k] < 0): continue
            found.append(tuple(c))
    return list(dict.fromkeys(found))


def _scan(P, Q, n, positive, gen, max_weight, samples):
    grid = [Fr(k) for k in (-3, -2, -1, 1, 2, 3)] + [Fr(1, 2), Fr(-1, 2), Fr(1, 3), Fr(-1, 3)]
    pick = lambda i: abs(gen.choice(grid)) if i in positive else gen.choice(grid)
    ivs, wit, zooms = [], [], []
    for w in _weights(n, 2 * max_weight):
        dP, iP = _initial(P, w)
        dQ, iQ = _initial(Q, w)
        if max(w) > max_weight and dP != dQ: continue      # beyond the box: only ties, where values vary
        for c0 in _zeros(iQ, n, positive, gen) if len(iQ) > 1 else []:
            if _ev(iP, c0) != 0 or dP != dQ: zooms.append((w, c0))
        pts = [tuple(pick(i) for i in range(n)) for _ in range(samples)]
        pts = [c for c in pts if _ev(iQ, c) != 0]
        if not pts: continue
        if dP > dQ:
            ivs.append((0, 0)); wit.append(("0", w, pts[0])); continue
        if dP < dQ:
            for c in pts:
                p = _ev(iP, c)
                if p != 0:
                    v = inf if (p > 0) == (_ev(iQ, c) > 0) else -inf
                    ivs.append((v, v)); wit.append((str(v), w, c))
            continue
        zp = [c for c in (_zeros(iP, n, positive, gen, tries=10) if len(iP) > 1 else []) if _ev(iQ, c) != 0][:3]
        pts += zp                                   # points with value exactly 0: segments from them reach 0
        for c in pts:
            v = _ev(iP, c) / _ev(iQ, c); ivs.append((v, v)); wit.append((str(v), w, c))
        cl = [_climb(iP, iQ, c, sgn, gen) for c in pts[:2] for sgn in (1, -1)]
        cl = [c for c in cl if _ev(iQ, c) != 0 and all(c[i] > 0 for i in positive)]
        pts += cl
        for c in cl:
            v = _ev(iP, c) / _ev(iQ, c); ivs.append((v, v)); wit.append((str(v), w, c))
        for a, b in list(zip(pts, pts[1:])) + [(pts[0], c) for c in cl] + [(z, c) for z in zp for c in cl]:
            seg = _segment(_along(iP, a, b), _along(iQ, a, b))
            if seg:
                ivs += seg; wit.append(("segment", w, (a, b)))
    return ivs, wit, zooms


def inner_bound(expr, names, max_weight=4, samples=6, seed=0, depth=1, max_zooms=8):
    """proven limit values of the rational function expr at the origin, any number of variables.
       Returns {"set", "components", "whole_circle", "limit_exists" (False when two values are
       proven, else None: not decided), "proven", "witnesses": [(value or "segment", weight,
       point(s)[, zoom chain])]}; curve_point(witness, s) gives the proving curve.

       Zooming (depth >= 1): at a rational point c0 where in_w Q vanishes, the curves
       x_j = c0_j t^w_j, x_i = t^w_i (c0_i + y_i), with t -> 0+ and y -> 0, still go to the
       origin, so the limit values of the new function of (t, y) at 0 with t > 0 are limit values
       of f. This reaches curves with two or more terms, such as x = y sqrt(-z) (through the
       rational point (1, 1, -1) of x^2 + z y^2)."""
    P, Q = parse_n(expr, names)
    n = len(names)
    gen = random.Random(seed)
    ivs, wit = [], []
    todo = [(P, Q, (), [])]
    for level in range(depth + 1):
        nxt = []
        for P1, Q1, pos, chain in todo:
            a, b, zooms = _scan(P1, Q1, n, pos, gen, max_weight, samples)
            ivs += a
            wit += [x + (chain,) if chain else x for x in b]
            if level == depth: continue
            seen = set()
            for w, c0 in zooms:
                if len(nxt) >= max_zooms: break
                j = next(i for i in range(n) if c0[i] != 0)
                key = (w, tuple(c0[i] / c0[j] ** Fr(w[i], w[j]) if w[i] % w[j] == 0 else c0[i] for i in range(n)))
                if key in seen: continue
                seen.add(key)
                nxt.append((_sub(P1, w, c0, j, n), _sub(Q1, w, c0, j, n), (0,), chain + [(w, c0, j)]))
        todo = nxt
        if _covers_circle(_union([(float(a) if not isinf(a) else a, float(b) if not isinf(b) else b)
                                  for a, b in ivs])): break
    U = _union([(float(a) if not isinf(a) else a, float(b) if not isinf(b) else b) for a, b in ivs])
    whole = _covers_circle(U)
    if whole: U = [(-inf, inf)]
    distinct = len(U) > 1 or (U and U[0][0] != U[0][1])
    fmt = lambda v: "+oo" if v == inf else "-oo" if v == -inf else f"{v:.6g}"
    text = " u ".join("{" + fmt(a) + "}" if a == b else f"[{fmt(a)}, {fmt(b)}]" for a, b in U)
    # Part 1 of docs/theory/n-variables.md: the answer is connected on the circle, so of the gaps
    # between the proven pieces, at most one is outside it; every other gap is filled
    gaps = [] if whole else [(U[i][1], U[i + 1][0]) for i in range(len(U) - 1)]
    if not whole and not (U[0][0] == -inf and U[-1][1] == inf): gaps.append((U[-1][1], U[0][0]))
    note = None
    if len(gaps) >= 2:
        note = ("the answer is connected, so it contains every gap between these pieces except at most "
                "one of: " + ", ".join(f"({fmt(a)}, {fmt(b)})" for a, b in gaps))
    return {"set": text, "components": U, "whole_circle": whole, "limit_exists": False if distinct else None,
            "proven": "inner bound" if not whole else "exact (the whole circle)", "witnesses": wit,
            "gaps": gaps, "connected_note": note}
