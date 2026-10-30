"""Complete answers in any number of variables, for denominators with an isolated zero.

Suppose the denominator Q has, for some weight vector w, a positive definite initial form:
Q_w(u) > 0 for u != 0, where Q = Q_w + (terms of higher weighted degree d). Weighted scaling
x = lambda . u (x_i = lambda^w_i u_i) maps every point near 0 to lambda -> 0+ and u on the
"weighted cube" { max_i |u_i| = 1 } (a compact set: the 2n faces where one |u_i| = 1). There
Q(lambda . u) = lambda^d (Q_w(u) + O(lambda)) and P(lambda . u) = lambda^e (P_w(u) + O(lambda)),
uniformly, with Q_w >= m > 0. Hence, with g = P_w / Q_w on the faces:
    e > d                      the limit exists and is 0;
    e = d                      the set of limit values is exactly the range [min g, max g] (the
                               faces form a connected set for n >= 2, and g is continuous there);
    e < d, P_w takes both signs  every value: on each sphere f runs from large positive to large
                               negative, so it takes every real value (intermediate values);
    e < d, P_w > 0 (or < 0) on the faces   the limit is +oo (or -oo).
All of the conditions are proven by rigorous interval arithmetic (outward-rounded floating
point) with branch and bound over the faces; min g and max g are enclosed to within tol, with
the lower ends of the enclosures attained by explicit points (so they are proven values as well).
This is the classical weighted polar-coordinates argument, made automatic and rigorous; it covers
the usual textbook problems in any number of variables. When the hypothesis cannot be proven,
the result falls back to fmsl.witness.inner_bound (proven values, possibly incomplete).

    from fmsl.nvar import limit_n
    limit_n("x*w/(x^2+y^2+z^2+w^2)", ["x", "y", "z", "w"])["set"]     # '[-1/2, 1/2] (+- 1e-06)'
"""
import heapq, math
from fractions import Fraction as Fr
from itertools import product
from math import gcd, inf, nextafter

from .witness import parse_n, _initial, inner_bound

__all__ = ["limit_n"]


# ------------------------------------------------------------------ outward-rounded intervals
def _dn(x): return nextafter(x, -inf)
def _up(x): return nextafter(x, inf)


def _imul(a, b):
    p = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return (_dn(min(p)), _up(max(p)))


def _ipow(a, k):
    if k == 0: return (1.0, 1.0)
    lo, hi = a
    if k % 2 == 0 and lo <= 0 <= hi:
        m = max(-lo, hi)
        return (0.0, _up(m ** k) * (1 + 1e-15))
    c = (lo ** k, hi ** k)
    return (_dn(min(c)) - abs(min(c)) * 1e-15, _up(max(c)) + abs(max(c)) * 1e-15)


def _ieval(terms, box):
    lo = hi = 0.0
    for k, (clo, chi) in terms:
        m = (clo, chi)
        for iv, e in zip(box, k):
            if e: m = _imul(m, _ipow(iv, e))
        lo, hi = _dn(lo + m[0]), _up(hi + m[1])
    return lo, hi


def _fterms(P):
    out = []
    for k, v in P.items():
        f = float(v)
        out.append((k, (_dn(f), _up(f)) if Fr(f) != v else (f, f)))
    return out


def _pev(P, x):
    s = Fr(0)
    for k, v in P.items():
        m = v
        for xi, e in zip(x, k):
            if e: m *= xi ** e
        s += m
    return s


def _faces(n):
    for i in range(n):
        for s in (-1, 1):
            yield i, s


def _box(face, n):
    i, s = face
    return tuple((float(s), float(s)) if j == i else (-1.0, 1.0) for j in range(n))


def _split(box, fixed):
    widths = [(b - a, j) for j, (a, b) in enumerate(box) if j != fixed]
    _, j = max(widths)
    a, b = box[j]; m = (a + b) / 2
    return (box[:j] + ((a, m),) + box[j + 1:], box[:j] + ((m, b),) + box[j + 1:])


def _centre(box):
    return tuple(Fr((a + b) / 2) for a, b in box)


def _positive(Qw, n, budget=20000):
    """True if Qw > 0 on every face is proven, False if a point with Qw <= 0 is found, else None"""
    import random
    from collections import deque
    gen = random.Random(0)
    T = _fterms(Qw)
    for i, sg in _faces(n):                          # a quick look for a point where Qw <= 0
        for _ in range(50):
            x = [Fr(gen.randint(-8, 8), 8) for _ in range(n)]; x[i] = Fr(sg)
            if _pev(Qw, x) <= 0: return False
    for face in _faces(n):
        stack = deque([_box(face, n)])
        while stack:
            budget -= 1
            if budget < 0: return None
            b = stack.popleft()
            lo, hi = _ieval(T, b)
            if lo > 0: continue
            if _pev(Qw, _centre(b)) <= 0: return False
            stack += _split(b, face[0])
    return True


def _deriv(P, j):
    out = {}
    for k, v in P.items():
        if k[j]:
            kk = tuple(e - (i == j) for i, e in enumerate(k)); out[kk] = out.get(kk, 0) + v * k[j]
    return out


def _negative(H, n, budget=200000):
    """True if H < 0 on every face is proven (mean-value interval form, branch and bound);
       otherwise a rational point where H >= 0 (or None when the budget runs out)"""
    T = _fterms(H)
    G = [_fterms(_deriv(H, j)) for j in range(n)]
    for face in _faces(n):
        stack = [_box(face, n)]
        while stack:
            budget -= 1
            if budget < 0: return None
            b = stack.pop()
            ub = _ieval(T, b)[1]
            if ub < 0: continue
            c = tuple((lo + hi) / 2 for lo, hi in b)
            hc = _ieval(T, tuple((x, x) for x in c))[1]
            mv = hc
            for j in range(n):
                if j == face[0]: continue
                g = _ieval(G[j], b)
                r = (b[j][1] - b[j][0]) / 2
                mv = _up(mv + _up(max(abs(g[0]), abs(g[1])) * r))
            if mv < 0: continue
            cf = _centre(b)
            if _pev(H, cf) >= 0: return cf
            stack += _split(b, face[0])
    return True


def _climb_faces(Pw, Qw, n, sign, starts=12, steps=400):
    import random
    gen = random.Random(1)
    fP = [(k, float(v)) for k, v in Pw.items()]; fQ = [(k, float(v)) for k, v in Qw.items()]
    def ev(T, x):
        s = 0.0
        for k, v in T:
            m = v
            for xi, e in zip(x, k):
                if e: m *= xi ** e
            s += m
        return s
    best, arg = -inf, None
    for i, sg in _faces(n):
        for _ in range(starts):
            x = [gen.uniform(-1, 1) for _ in range(n)]; x[i] = float(sg)
            q = ev(fQ, x); v = sign * ev(fP, x) / q; h = 0.3
            for _ in range(steps):
                y = [min(1.0, max(-1.0, xj + h * gen.uniform(-1, 1))) for xj in x]; y[i] = float(sg)
                q = ev(fQ, y)
                if q > 0 and sign * ev(fP, y) / q > v: x, v = y, sign * ev(fP, y) / q
                else: h *= 0.985
            if v > best: best, arg = v, x
    return tuple(Fr(t).limit_denominator(10 ** 9) for t in arg)


def _range(Pw, Qw, n, sign, tol):
    """(attained, bound, point): max over the faces of sign * Pw/Qw is in [attained, bound], with
       attained reached exactly at the rational point and bound - attained <= tol; or None"""
    arg = _climb_faces(Pw, Qw, n, sign)
    for _ in range(6):
        best = sign * _pev(Pw, arg) / _pev(Qw, arg)
        beta = best + Fr(tol)
        H = {k: sign * v for k, v in Pw.items()}
        for k, v in Qw.items():
            H[k] = H.get(k, 0) - beta * v
        H = {k: v for k, v in H.items() if v != 0}
        r = _negative(H, n)
        if r is True: return best, float(beta), arg
        if r is None: return None
        arg = r                                     # a better point: try again from it
    return None


def _weights(n, k):
    for w in product(range(1, k + 1), repeat=n):
        g = 0
        for x in w: g = gcd(g, x)
        if g == 1: yield w


def _fmt(v):
    if v == inf: return "+oo"
    if v == -inf: return "-oo"
    f = Fr(v).limit_denominator(1000)
    return str(f) if abs(float(f) - float(v)) < 1e-9 else f"{float(v):.10g}"


def limit_n(expr, names, tol=1e-6, max_weight=6):
    """the set of limit values of a rational function at the origin, n variables. Returns
       {"set", "components", "limit_exists", "limit", "proven", "method", "explanation"}.
       proven is True for a complete answer (ends enclosed within tol), else "inner bound"."""
    P, Q = parse_n(expr, names)
    n = len(names)
    for w in sorted(_weights(n, max_weight), key=lambda w: (sum(w), w)):
        d, Qw = _initial(Q, w)
        if len(Qw) < n: continue                     # a positive definite form involves every variable
        if _positive(Qw, n) is not True: continue
        e, Pw = _initial(P, w)
        base = {"method": "weighted scaling", "weight": w, "proven": True}
        if e > d:
            return {**base, "set": "{0}", "components": [(0.0, 0.0)], "limit_exists": True, "limit": "0",
                    "explanation": f"Q's initial form for weight {w} is positive definite (proven), and P has "
                                   f"higher weighted order ({e} > {d}), so f -> 0."}
        if e < d:
            pos = _positive(Pw, n); neg = _positive({k: -v for k, v in Pw.items()}, n)
            if pos is True or neg is True:
                s = "+oo" if pos is True else "-oo"
                return {**base, "set": "{" + s + "}", "components": [(inf, inf) if pos else (-inf, -inf)],
                        "limit_exists": True, "limit": s,
                        "explanation": f"weight {w}: Q's initial form is positive definite, P has lower order "
                                       f"({e} < {d}) and its initial form has one sign, so f -> {s}."}
            if pos is False and neg is False:
                return {**base, "set": "[-oo, +oo]", "components": [(-inf, inf)], "limit_exists": False,
                        "limit": None,
                        "explanation": f"weight {w}: Q's initial form is positive definite, P has lower order "
                                       f"({e} < {d}) and its initial form takes both signs, so on every small "
                                       f"sphere f takes every real value."}
            break
        # e == d: the range of Pw/Qw on the faces
        c0 = next(iter(Pw.items()))
        for t in sorted({tol, max(tol, 1e-4), max(tol, 1e-3)}):   # a coarser enclosure when a fine one is too costly
            hi = _range(Pw, Qw, n, 1, t); lo = _range(Pw, Qw, n, -1, t)
            if hi is not None and lo is not None: tol = t; break
        if hi is None or lo is None: break
        a_att, a_bound = -lo[0], -lo[1]; b_att, b_bound = hi[0], hi[1]
        if a_att == b_att and all(_pev(Pw, x) == a_att * _pev(Qw, x) for x in (lo[2], hi[2])) and \
                not {k: v for k, v in _sub_lin(Pw, Qw, a_att).items()}:
            s = _fmt(a_att)
            return {**base, "set": "{" + s + "}", "components": [(float(a_att), float(a_att))],
                    "limit_exists": True, "limit": s,
                    "explanation": f"weight {w}: P's initial form is exactly {s} times Q's, which is positive "
                                   f"definite, so f -> {s}."}
        A = _fmt(a_att) if abs(float(a_att) - a_bound) < 1e-12 else _fmt(a_att)
        B = _fmt(b_att)
        return {**base, "set": f"[{A}, {B}] (+- {tol:g})", "components": [(float(a_att), float(b_att))],
                "enclosure": [(a_bound, float(a_att)), (float(b_att), b_bound)],
                "limit_exists": False, "limit": None,
                "explanation": f"weight {w}: Q's initial form is positive definite (proven) and P's has the "
                               f"same order {d}, so the limit values are exactly the range of their ratio on "
                               f"the weighted cube: its minimum lies in [{a_bound:.10g}, {float(a_att):.10g}] "
                               f"and its maximum in [{float(b_att):.10g}, {b_bound:.10g}] (rigorous interval "
                               f"bounds; the inner ends are attained)."}
    r = inner_bound(expr, names)
    return {"set": r["set"], "components": r["components"], "limit_exists": r["limit_exists"], "limit": None,
            "proven": True if r["whole_circle"] else "inner bound", "method": "witness curves",
            "explanation": "no weight with a positive definite initial denominator was proven; the values "
                           "shown are proven by explicit curves" + ("" if r["whole_circle"] else
                           ", but others may be missing")}


def _sub_lin(Pw, Qw, c):
    out = dict(Pw)
    for k, v in Qw.items():
        out[k] = out.get(k, 0) - c * v
        if out[k] == 0: del out[k]
    return out
