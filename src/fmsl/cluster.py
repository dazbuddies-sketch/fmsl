"""
fmsl.cluster -- the fluid state of a two-variable 0/0 point as an exact arc.

Theorem (connectedness).  Let f = N/D with N, D real polynomials without a common factor
through the origin.  Then f is a continuous map from a punctured disc around 0 into the
circle  R u {oo}  (the real projective line), and its CLUSTER SET -- all values f approaches
near the origin -- is a nonempty, closed, CONNECTED subset of that circle:
        a point (a limit, possibly oo),  a closed arc,  or the whole circle.
Proof: C = intersection over r of closure(f(punctured disc of radius r)); each set is the
closure of a continuous image of a connected set, hence compact and connected, and they
decrease; a decreasing intersection of compact connected sets is connected.

Every value in C is reached along some Puiseux arc (curve selection lemma), and the profile
classifies every such arc, so C is the closure of the union of
    * the exact values along fixed paths,
    * the images R(R) of the path-dependence functions R(c) = fN(c)/fD(c),
    * oo where f is unbounded, and everything where an "every value" entry was proven.
Computing that union exactly and checking that it IS connected is an internal consistency
test: a disconnected answer would contradict the theorem, so some family of paths would be
missing.
"""
from fractions import Fraction as Fr
from functools import cmp_to_key

from .algebraic import (A, AlgField, isolate_real_roots, squarefree, pmul, psub, pderiv,
                        ptrim, peval, sign_of, same_value)


class _Inf:
    def __init__(s, sgn): s.sgn = sgn
    def __repr__(s): return "+oo" if s.sgn > 0 else "-oo"

POS, NEG = _Inf(1), _Inf(-1)
_INEXACT = [False]


def _fields_cmp(u, v):
    """exact comparison of two algebraic numbers from different fields (u != v assumed)"""
    from .algebraic import _ieval
    for _ in range(200):
        a1, b1 = _ieval(u.F.reduce(u.c), u.F.lo, u.F.hi)
        a2, b2 = _ieval(v.F.reduce(v.c), v.F.lo, v.F.hi)
        if b1 < a2: return -1
        if b2 < a1: return 1
        u.F.refine(); v.F.refine()
    if _exact_equal(u, v): return 0
    for _ in range(2000):                              # not equal: they must separate
        a1, b1 = _ieval(u.F.reduce(u.c), u.F.lo, u.F.hi)
        a2, b2 = _ieval(v.F.reduce(v.c), v.F.lo, v.F.hi)
        if b1 < a2: return -1
        if b2 < a1: return 1
        u.F.refine(); v.F.refine()
    raise ArithmeticError("could not separate two algebraic numbers")


def _charpoly(u):
    """characteristic polynomial (over Q) of multiplication by u on Q[t]/p: u is a root of it"""
    from .algebraic import pdivmod
    P = u.F.p; n = len(P) - 1
    cols = []
    for k in range(n):
        r = u.F.reduce(pmul(u.c, [Fr(0)] * k + [Fr(1)]))
        cols.append([r[i] if i < len(r) else Fr(0) for i in range(n)])
    M = [[cols[j][i] for j in range(n)] for i in range(n)]
    # Faddeev-LeVerrier, exact
    coeffs = [Fr(1)]; Mk = [[Fr(0)] * n for _ in range(n)]; I = [[Fr(int(i == j)) for j in range(n)] for i in range(n)]
    c_prev = Fr(1)
    for k in range(1, n + 1):
        Mk = [[sum(M[i][l] * Mk[l][j] for l in range(n)) + c_prev * I[i][j] for j in range(n)] for i in range(n)]
        AM = [[sum(M[i][l] * Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        c_prev = -sum(AM[i][i] for i in range(n)) / k
        coeffs.append(c_prev)
    return list(reversed(coeffs))                      # ascending powers


def _exact_equal(u, v):
    """exact equality of algebraic numbers from different fields"""
    from .algebraic import _ieval, count_roots
    q = squarefree(_charpoly(u))
    qv = Fr(0)                                         # q(v), evaluated exactly in v's field
    for cf in reversed(q): qv = qv * v + cf
    if (qv != 0) if not isinstance(qv, A) else bool(qv): return False
    # both are roots of q: equal iff they lie in the same isolating interval of q
    for _ in range(2000):
        a1, b1 = _ieval(u.F.reduce(u.c), u.F.lo, u.F.hi)
        a2, b2 = _ieval(v.F.reduce(v.c), v.F.lo, v.F.hi)
        lo, hi = min(a1, a2), max(b1, b2)
        pad = (hi - lo) or Fr(1, 10 ** 200)              # shrinks with the enclosures
        if count_roots(q, lo - pad, hi) == 1: return True
        if b1 < a2 or b2 < a1: return False
        u.F.refine(); v.F.refine()
    raise ArithmeticError("could not decide equality")


def cmp(u, v):
    """total order on extended reals: -oo < finite < +oo, exact for Fractions and algebraics"""
    if isinstance(u, _Inf) or isinstance(v, _Inf):
        ku = u.sgn if isinstance(u, _Inf) else 0
        kv = v.sgn if isinstance(v, _Inf) else 0
        return (ku > kv) - (ku < kv)
    if isinstance(u, A) and isinstance(v, A) and u.F is not v.F:
        eq = same_value(u, v)
        if eq is True: return 0
        return _fields_cmp(u, v)
    return sign_of(u - v)


def _bounds(r):
    """rational (lo, hi) with lo <= r <= hi, refinable (by refining r's field).
       For an algebraic number g(alpha) this is an interval enclosure of g over alpha's
       isolating interval -- NOT alpha's interval itself."""
    if isinstance(r, A):
        from .algebraic import _ieval
        if r.c == [Fr(0), Fr(1)]: return r.F.lo, r.F.hi           # alpha itself
        return _ieval(r.F.reduce(r.c), r.F.lo, r.F.hi)
    return r, r


def _rational_between(r1, r2):
    """a rational strictly between r1 < r2 (exact numbers)"""
    if cmp(r1, r2) >= 0: raise ArithmeticError("rational_between: r1 >= r2")
    for _ in range(4000):
        lo1, hi1 = _bounds(r1); lo2, hi2 = _bounds(r2)
        if hi1 < lo2:
            from .multivar import _simplest_between
            return _simplest_between(hi1, lo2)                 # small denominators: fast later work
        if hi1 == lo2 and not isinstance(r1, A) and not isinstance(r2, A): return hi1
        if isinstance(r1, A): r1.F.refine()
        if isinstance(r2, A): r2.F.refine()
    raise ArithmeticError("could not separate")


def _real_roots(p, lo=None, hi=None):
    """distinct real roots of a rational polynomial, as exact numbers, sorted.
       With lo/hi (rationals), only roots in (lo, hi] are guaranteed (others may appear).
       Rational roots are found exactly and divided out first, so that an algebraic number
       is only used for a genuinely irrational root."""
    from .multivar import _rational_roots
    from .algebraic import pdivmod
    p = squarefree(ptrim([Fr(v) for v in p]))
    out = []
    if p and p[0] == 0:
        out.append(Fr(0)); p = p[1:]
    for r in _rational_roots(p, lo, hi):
        out.append(r); p = pdivmod(p, [-r, Fr(1)])[0]
    for iv in (isolate_real_roots(p, lo, hi) if len(p) > 1 else []):
        if iv[0] == "rational": out.append(iv[1])
        else: out.append(A.root(AlgField(p, iv[1], iv[2])))
    return sorted(out, key=cmp_to_key(cmp))


def image_of_rational_function(n, d):
    """closure of { n(c)/d(c) : c real, d(c) != 0 } as a list of closed extended intervals.
       n, d: coefficient lists over Q, coprime, ratio not constant.  Exact."""
    n, d = ptrim([Fr(v) for v in n]), ptrim([Fr(v) for v in d])
    R = lambda c: peval(n, c) / peval(d, c)
    g = psub(pmul(pderiv(n), d), pmul(n, pderiv(d)))          # critical points: R' = 0
    prod = squarefree(pmul(g, d)) if g else squarefree(d)
    bps = _real_roots(prod) if len(prod) > 1 else []
    is_pole = [sign_of(peval(d, b) if not isinstance(b, A) else _evalA(d, b)) == 0 for b in bps]
    dn, dd = len(n) - 1, len(d) - 1
    lim = Fr(0) if dn < dd else (n[-1] / d[-1] if dn == dd else None)   # None: infinite
    pieces = []
    edges = [None] + list(range(len(bps))) + [None]
    for k in range(len(bps) + 1):
        L = None if k == 0 else bps[k - 1]
        Rr = None if k == len(bps) else bps[k]
        if L is None and Rr is None: m1, m2 = Fr(0), Fr(1)
        elif L is None: lo = _bounds(Rr)[0] - 1; m1, m2 = lo - 1, lo
        elif Rr is None: hi = _bounds(L)[1] + 1; m1, m2 = hi, hi + 1
        else:
            m2 = _rational_between(L, Rr)
            m1 = _rational_between(L, m2)
        v1, v2 = R(m1), R(m2)
        inc = v2 > v1
        def end(b_idx, is_left):
            if b_idx is None:                                  # c -> -oo (left) or +oo (right)
                if lim is not None: return lim
            elif not is_pole[b_idx]:
                return _value(n, d, bps[b_idx])
            # infinite end: direction of monotone motion decides the sign
            if is_left: return NEG if inc else POS
            return POS if inc else NEG
        eL = end(None if L is None else k - 1, True)
        eR = end(None if Rr is None else k, False)
        lo_, hi_ = (eL, eR) if cmp(eL, eR) <= 0 else (eR, eL)
        pieces.append((lo_, hi_))
    return pieces


def _evalA(p, b):
    acc = Fr(0)
    for v in reversed(p): acc = acc * b + v
    return acc


def _value(n, d, b):
    return _evalA(n, b) / _evalA(d, b)


def union_on_circle(items):
    """items: closed intervals (lo, hi) of the extended line, or the string 'all'.
       Returns (components, connected, whole) where +oo and -oo are the same point oo."""
    if any(it == "all" for it in items): return [(NEG, POS)], True, True
    iv = []
    for lo, hi in items:
        if isinstance(lo, _Inf) and lo.sgn > 0: lo = hi = POS     # the point oo
        iv.append((lo, hi))
    iv.sort(key=cmp_to_key(lambda s, t: cmp(s[0], t[0])))
    comps = []
    for lo, hi in iv:
        if comps and cmp(lo, comps[-1][1]) <= 0:
            if cmp(hi, comps[-1][1]) > 0: comps[-1] = (comps[-1][0], hi)
        else:
            comps.append((lo, hi))
    # glue through oo
    touches_neg = [c for c in comps if isinstance(c[0], _Inf) and c[0].sgn < 0]
    touches_pos = [c for c in comps if isinstance(c[1], _Inf) and c[1].sgn > 0]
    whole = len(comps) == 1 and isinstance(comps[0][0], _Inf) and comps[0][0].sgn < 0 \
        and isinstance(comps[0][1], _Inf) and comps[0][1].sgn > 0
    n = len(comps)
    if touches_neg and touches_pos and touches_neg[0] is not touches_pos[-1]:
        n -= 1                                                  # one arc through oo
        if n == 1 and cmp(touches_pos[-1][0], touches_neg[0][1]) <= 0:
            whole = True                                        # the two ends overlap
    return comps, n == 1, whole


def _fmt(v):
    if isinstance(v, _Inf): return "oo"
    if isinstance(v, A): return f"{float(v):.6g} (algebraic)"
    if isinstance(v, Fr) and v.denominator > 10 ** 12: return f"{float(v):.6g} (approx.)"
    return str(v)


def describe_arc(comps, connected, whole):
    if whole: return "every value: the whole circle R u {oo}"
    neg = [c for c in comps if isinstance(c[0], _Inf) and c[0].sgn < 0]
    pos = [c for c in comps if isinstance(c[1], _Inf) and c[1].sgn > 0]
    if connected and neg and pos and len(comps) == 2:
        a, b = pos[-1][0], neg[0][1]
        if isinstance(a, _Inf): return f"the arc from oo to {_fmt(b)}  ([-oo, {_fmt(b)}], oo included)"
        if isinstance(b, _Inf): return f"the arc from {_fmt(a)} to oo  ([{_fmt(a)}, +oo], oo included)"
        return f"the arc from {_fmt(a)} up through oo to {_fmt(b)}  ([{_fmt(a)}, oo] u [-oo, {_fmt(b)}])"
    parts = []
    for lo, hi in comps:
        if cmp(lo, hi) == 0: parts.append("{" + _fmt(lo) + "}")
        else: parts.append(f"[{'-oo' if isinstance(lo, _Inf) else _fmt(lo)}, {'+oo' if isinstance(hi, _Inf) else _fmt(hi)}]")
    return " u ".join(parts)


def _inf_signs(fN, fD):
    """signs taken by fN(c)/fD(c) for real c != 0 (fN, fD rational): the signs of f -> oo in a
       sector where f is unbounded.  None if the coefficients are algebraic (not decided)."""
    if any(isinstance(v, A) for v in list(fN) + list(fD)):
        nz = lambda p: [(j, v) for j, v in enumerate(p) if v != 0]
        a, b = nz(fN), nz(fD)
        if len(a) == 1 and len(b) == 1:                     # monomials n c^j / d c^k
            sg = sign_of(a[0][1]) * sign_of(b[0][1])
            return {sg, -sg} if (a[0][0] - b[0][0]) % 2 else {sg}
        return None
    bps = _real_roots(squarefree(pmul(fN, fD)) + [])
    bps = sorted(set_unique(bps + [Fr(0)]), key=cmp_to_key(cmp))
    tests = []
    if bps:
        tests.append(_bounds(bps[0])[0] - 1); tests.append(_bounds(bps[-1])[1] + 1)
        for u, v in zip(bps, bps[1:]): tests.append(_rational_between(u, v))
    else:
        tests = [Fr(-1), Fr(1)]
    out = set()
    for t in tests:
        n, d = peval(fN, t), peval(fD, t)
        if n != 0 and d != 0: out.add(sign_of(n) * sign_of(d))
    return out


def set_unique(xs):
    out = []
    for x in xs:
        if not any(cmp(x, y) == 0 for y in out): out.append(x)
    return out


def _union_line(items):
    """union on the extended line [-oo, +oo] (+oo and -oo distinct); items as in union_on_circle"""
    if any(it == "all" for it in items): return [(NEG, POS)]
    iv = sorted(items, key=cmp_to_key(lambda s, t: cmp(s[0], t[0])))
    comps = []
    for lo, hi in iv:
        if comps and cmp(lo, comps[-1][1]) <= 0:
            if cmp(hi, comps[-1][1]) > 0: comps[-1] = (comps[-1][0], hi)
        else:
            comps.append((lo, hi))
    return comps


def _describe_line(comps):
    parts = []
    for lo, hi in comps:
        f = lambda v: ("+oo" if v.sgn > 0 else "-oo") if isinstance(v, _Inf) else _fmt(v)
        parts.append("{" + f(lo) + "}" if cmp(lo, hi) == 0 else f"[{f(lo)}, {f(hi)}]")
    return " u ".join(parts)


def cluster_set(entries):
    """exact cluster set from profile entries.
       Returns dict: text (on the circle R u {oo}), signed_text / signed_components (on the
       extended line, +oo and -oo distinct), components, connected, whole, exact, complete."""
    items, exact, complete = [], True, True
    signs_exact = True
    _INEXACT[0] = False
    for e in entries:
        k = e["kind"]
        if k == "value":
            items.append((e["value"], e["value"]))
        elif k == "inf":
            sg = e.get("signs")
            if sg is None and "fN" in e: sg = _inf_signs(e["fN"], e["fD"])
            if sg is None: sg = {1, -1}; signs_exact = False
            for t in sg: items.append((POS, POS) if t > 0 else (NEG, NEG))
        elif k == "all":
            items.append("all")
        elif k == "unresolved":
            complete = False
        elif k == "depends":
            fN, fD = e["R"]
            if not all(not isinstance(v, A) for v in list(fN) + list(fD)):
                # coefficients in a number field: approximate them (40 digits); flag as inexact
                fN = [Fr(v.F.approx(v.c)) if isinstance(v, A) else Fr(v) for v in fN]
                fD = [Fr(v.F.approx(v.c)) if isinstance(v, A) else Fr(v) for v in fD]
                exact = False
            from .algebraic import pgcd, pdivmod
            g = pgcd(fN, fD)
            n, d = pdivmod(fN, g)[0], pdivmod(fD, g)[0]
            if len(ptrim(n)) <= 1 and len(ptrim(d)) <= 1:
                v = (ptrim(n) or [Fr(0)])[0] / d[0]; items.append((v, v)); continue
            items += image_of_rational_function(n, d)
    signed = _union_line(items)
    def to_circle(it):                       # on the circle the point -oo is the point oo
        if it != "all" and isinstance(it[0], _Inf) and it[0].sgn < 0 and isinstance(it[1], _Inf) and it[1].sgn < 0:
            return (POS, POS)
        return it
    circ = [to_circle(it) for it in items]
    comps, connected, whole = union_on_circle(circ)
    exact = exact and not _INEXACT[0]
    return {"components": comps, "connected": connected, "whole": whole, "exact": exact,
            "complete": complete, "text": describe_arc(comps, connected, whole),
            "signed_components": signed, "signed_text": _describe_line(signed), "signs_exact": signs_exact}
