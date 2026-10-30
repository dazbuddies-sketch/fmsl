"""Indeterminacy curves in three variables (docs/theory/n-variables.md, Part 2).

For f = P/Q in three variables, P and Q can both vanish along a whole curve through the origin,
where f is 0/0 at every point. This never happens in two variables.

Theorem (Part 2 of the note). If there are points p != 0 of {P = Q = 0}, arbitrarily close to the
origin, where dP(p) and dQ(p) are linearly independent, then every value, +oo and -oo included, is
a limit value of f at the origin.

Proof: near such p, P and Q are two of a set of local coordinates, so f looks like u/v, which takes
every value, and tends to +oo and to -oo, in every neighbourhood of p; letting p -> 0 gives the
same at the origin.

This module looks for such points on lines through the origin with small integer directions v: the
line {t v} lies in {P = Q = 0} exactly when P(t v) and Q(t v) are the zero polynomial in t, and dP
and dQ are independent at t v for all small t != 0 exactly when some 2x2 minor of the matrix with
rows grad P(t v), grad Q(t v) is a nonzero polynomial in t. Both tests are exact. Finding such a
line proves the answer; not finding one proves nothing, and the engine then runs as before.
"""
from fractions import Fraction as Fr
from itertools import product
from math import gcd

__all__ = ["clean_line"]


def _on_line(P, v):
    """P(t v) as a polynomial in t: {degree: coefficient}"""
    out = {}
    for (i, j, k), c in P.items():
        d = i + j + k
        out[d] = out.get(d, 0) + c * Fr(v[0]) ** i * Fr(v[1]) ** j * Fr(v[2]) ** k
    return {d: c for d, c in out.items() if c != 0}


def _grad(P):
    g = [{}, {}, {}]
    for (i, j, k), c in P.items():
        for axis, e in enumerate((i, j, k)):
            if e:
                key = tuple(x - (1 if a == axis else 0) for a, x in enumerate((i, j, k)))
                g[axis][key] = g[axis].get(key, 0) + c * e
    return g


def _pmul(a, b):
    out = {}
    for d1, c1 in a.items():
        for d2, c2 in b.items():
            out[d1 + d2] = out.get(d1 + d2, 0) + c1 * c2
    return {d: c for d, c in out.items() if c != 0}


def _psub(a, b):
    out = dict(a)
    for d, c in b.items(): out[d] = out.get(d, 0) - c
    return {d: c for d, c in out.items() if c != 0}


def directions(limit=3):
    seen = []
    for v in product(range(-limit, limit + 1), repeat=3):
        if not any(v): continue
        g = 0
        for x in v: g = gcd(g, abs(x))
        if g != 1: continue
        if next(x for x in v if x) < 0: continue
        seen.append(v)
    return sorted(seen, key=lambda v: (sum(abs(x) for x in v), v))


def clean_line(P, Q, limit=3):
    """a direction v (integer triple) of a line through the origin inside {P = Q = 0} along which
       dP and dQ are independent near (not at) the origin, or None. P, Q: dicts {(i, j, k): coeff}"""
    if not P or not Q: return None
    gP, gQ = _grad(P), _grad(Q)
    for v in directions(limit):
        if _on_line(P, v) or _on_line(Q, v): continue          # the line is not inside {P = Q = 0}
        a = [_on_line(g, v) for g in gP]
        b = [_on_line(g, v) for g in gQ]
        for i, j in ((0, 1), (0, 2), (1, 2)):
            if _psub(_pmul(a[i], b[j]), _pmul(a[j], b[i])):   # a nonzero minor, as a polynomial in t
                return v
    return None
