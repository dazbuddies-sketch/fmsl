"""
fmsl.plane -- exact real algebra in two variables (c, d), for the three-variable engine.

  * resultants (Sylvester matrix, fraction-free Bareiss elimination over Q[c])
  * bivariate gcd (primitive pseudo-remainder sequence over Q[c][d])
  * does a real plane curve F(c, d) = 0 have real points (in a box, avoiding a finite set)?
    -- one-level cylindrical decomposition: breakpoints in c from the leading coefficient,
       the discriminant and extra polynomials; one rational sample per cell; Sturm counts in d
  * the exact range of a two-variable rational function R = P/Q over the plane or a box:
    candidate endpoints (critical values, boundary values, supplied values), then one rational
    test value per gap decides whether the gap belongs to the range.
All decisions are exact (rationals, real algebraic numbers, Sturm sequences).
"""
from fractions import Fraction as Fr
from functools import cmp_to_key

from .algebraic import (A, AlgField, ptrim, padd, psub, pmul, pdivmod, pgcd, pderiv, peval,
                        squarefree, count_roots, root_bound, isolate_real_roots, sign_of)
from .cluster import (cmp, POS, NEG, _Inf, _rational_between, _bounds, _real_roots,
                      set_unique)

# ----------------------------- bivariate polynomials {(j, k): coeff} = sum coeff c^j d^k
def bclean(P): return {k: v for k, v in P.items() if v != 0}

def badd(P, Q):
    R = dict(P)
    for k, v in Q.items(): R[k] = R.get(k, 0) + v
    return bclean(R)

def bscale(P, s): return bclean({k: v * s for k, v in P.items()})

def bmul(P, Q):
    R = {}
    for (a, b), u in P.items():
        for (c, d), v in Q.items():
            R[(a + c, b + d)] = R.get((a + c, b + d), 0) + u * v
    return bclean(R)

def bdc(P): return bclean({(j - 1, k): v * j for (j, k), v in P.items() if j > 0})
def bdd(P): return bclean({(j, k - 1): v * k for (j, k), v in P.items() if k > 0})

def in_d(P):
    """P as a polynomial in d: list (index = power of d) of univariate polynomials in c"""
    if not P: return []
    n = max(k for _, k in P)
    out = [[] for _ in range(n + 1)]
    for (j, k), v in P.items():
        row = out[k]
        while len(row) <= j: row.append(Fr(0))
        row[j] += v
    return [ptrim(r) for r in out]

def from_d(L):
    return bclean({(j, k): v for k, r in enumerate(L) for j, v in enumerate(r)})

def swap(P): return {(k, j): v for (j, k), v in P.items()}

def at_c(P, c0):
    """univariate polynomial in d of P(c0, d); c0 rational or algebraic"""
    L = in_d(P)
    return ptrim([peval_any(r, c0) for r in L])

def peval_any(p, x):
    acc = Fr(0)
    for v in reversed(p): acc = acc * x + v
    return acc

def beval(P, c0, d0):
    return sum((v * c0 ** j * d0 ** k for (j, k), v in P.items()), Fr(0))

# ----------------------------- determinants and resultants over Q[c]
def _exact_div(a, b):
    q, r = pdivmod(a, b)
    if r: raise ArithmeticError("Bareiss: inexact division")
    return q

def det_poly(M):
    """determinant of a square matrix with entries in Q[c] (fraction-free Bareiss)"""
    n = len(M)
    if n == 0: return [Fr(1)]
    M = [[ptrim(list(e)) for e in row] for row in M]
    sign, prev = 1, [Fr(1)]
    for k in range(n - 1):
        if not M[k][k]:
            sw = next((i for i in range(k + 1, n) if M[i][k]), None)
            if sw is None: return []
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = _exact_div(psub(pmul(M[k][k], M[i][j]), pmul(M[i][k], M[k][j])), prev)
            M[i][k] = []
        prev = M[k][k]
    d = M[n - 1][n - 1]
    return [v * sign for v in d] if d else []

def resultant_d(P, Q):
    """Res_d(P, Q) as a univariate polynomial in c"""
    A_, B_ = in_d(P), in_d(Q)
    m, n = len(A_) - 1, len(B_) - 1
    if m < 0 or n < 0: return []
    if m == 0: return _ppow(A_[0], n)
    if n == 0: return _ppow(B_[0], m)
    N = m + n
    rows = []
    for i in range(n):
        row = [[] for _ in range(N)]
        for k in range(m + 1): row[i + k] = A_[m - k]
        rows.append(row)
    for i in range(m):
        row = [[] for _ in range(N)]
        for k in range(n + 1): row[i + k] = B_[n - k]
        rows.append(row)
    return det_poly(rows)

def _ppow(p, e):
    r = [Fr(1)]
    for _ in range(e): r = pmul(r, p)
    return r

# ----------------------------- bivariate gcd (primitive PRS over Q[c][d])
def _content(L):
    g = []
    for r in L:
        if r: g = pgcd(g, r) if g else [v / r[-1] for v in r]
    return g or [Fr(1)]

def _prim(L):
    g = _content(L)
    return [ptrim(pdivmod(r, g)[0]) if r else [] for r in L], g

def _prem(A_, B_):
    """pseudo-remainder of A by B in d (coefficients in Q[c])"""
    A_ = [list(r) for r in A_]
    m, n = len(A_) - 1, len(B_) - 1
    lb = B_[-1]
    while len(A_) - 1 >= n and any(A_):
        k = len(A_) - 1 - n
        la = A_[-1]
        A_ = [pmul(r, lb) for r in A_]
        for i, b in enumerate(B_):
            A_[i + k] = psub(A_[i + k], pmul(la, b))
        while A_ and not A_[-1]: A_.pop()
    return A_

def bgcd(P, Q):
    """gcd of two bivariate polynomials (up to a rational factor)"""
    if not P: return dict(Q)
    if not Q: return dict(P)
    A_, ca = _prim(in_d(P)); B_, cb = _prim(in_d(Q))
    cont = pgcd(ca, cb) if len(ca) > 1 and len(cb) > 1 else [Fr(1)]
    if len(A_) < len(B_): A_, B_ = B_, A_
    while B_ and len(B_) > 1:
        R = _prem(A_, B_)
        if not R: break
        A_, B_ = B_, _prim(R)[0]
    G = B_ if (B_ and len(B_) > 1) else [[Fr(1)]]
    if len(B_) == 1 and B_:      # gcd has degree 0 in d
        G = [[Fr(1)]]
    return bmul(from_d(G), {(j, 0): v for j, v in enumerate(cont)})

def bdivexact(P, G):
    """P / G exactly (both bivariate), via division in d over Q(c) with exact remainders"""
    A_, B_ = in_d(P), in_d(G)
    n = len(B_) - 1; lb = B_[-1]
    Qt = [[] for _ in range(max(len(A_) - n, 1))]
    A_ = [list(r) for r in A_]
    while len(A_) - 1 >= n and any(A_):
        k = len(A_) - 1 - n
        t = _exact_div(A_[-1], lb)
        Qt[k] = t
        for i, b in enumerate(B_): A_[i + k] = psub(A_[i + k], pmul(t, b))
        while A_ and not A_[-1]: A_.pop()
    if any(A_): raise ArithmeticError("not an exact division")
    return from_d(Qt)

# ----------------------------- real points of a curve
def _poly_roots_in(p, lo=None, hi=None):
    """exact real roots of a rational univariate polynomial, optionally inside (lo, hi)"""
    p = ptrim(p)
    if len(p) <= 1: return []
    rl = lo if isinstance(lo, Fr) else None
    rh = hi if isinstance(hi, Fr) else None
    out = _real_roots(p, rl, rh)
    if lo is not None: out = [r for r in out if cmp(r, lo) > 0]
    if hi is not None: out = [r for r in out if cmp(r, hi) < 0]
    return out

def _samples(bps, lo=None, hi=None):
    """one rational point in each open cell of (lo, hi) cut at the sorted exact breakpoints"""
    pts = sorted(set_unique(bps), key=cmp_to_key(cmp))
    if not pts:
        if lo is None and hi is None: return [Fr(0)]
        if lo is None: return [hi - 1]
        if hi is None: return [lo + 1]
        return [(lo + hi) / 2]
    out = [_bounds(pts[0])[0] - 1 if lo is None else _rational_between(lo, pts[0])]
    for u, v in zip(pts, pts[1:]): out.append(_rational_between(u, v))
    out.append(_bounds(pts[-1])[1] + 1 if hi is None else _rational_between(pts[-1], hi))
    return out

def _count_real_d(g, dlo=None, dhi=None, avoid_zero=False):
    g = squarefree(ptrim(g))
    if len(g) <= 1: return 0
    B = root_bound(g)
    lo = -B - 1 if dlo is None else dlo
    hi = B + 1 if dhi is None else dhi
    n = count_roots(g, lo, hi)                    # roots in (lo, hi]
    if dhi is not None and peval(g, dhi) == 0: n -= 1
    if avoid_zero and peval(g, Fr(0)) == 0 and lo < 0 < hi: n -= 1
    return n

def has_real_points(F, box=None, extra=(), avoid_zero=False):
    """True if the curve F(c, d) = 0 has a real point with (c, d) in the open box
       (box = (lo, hi) for both coordinates, or None for the plane), off the finitely many
       points where the polynomials in `extra` vanish together with F.
       Complete for curve components of dimension 1 (isolated real points are not searched;
       callers only ask at regular values, where none exist)."""
    F = bclean(F)
    if not F: return True
    lo, hi = (None, None) if box is None else box
    L = in_d(F)
    if len(L) == 1:                                     # F = f(c): vertical lines
        return bool(_poly_roots_in(L[0], lo, hi)) and not (avoid_zero and all(cmp(r, 0) == 0 for r in _poly_roots_in(L[0], lo, hi)))
    Lp, cont = _prim(L)
    if len(cont) > 1:
        rts = [r for r in _poly_roots_in(cont, lo, hi) if not (avoid_zero and cmp(r, 0) == 0)]
        if rts: return True                             # a vertical line c = root
    Fp = from_d(Lp)
    bps = []
    bps += _poly_roots_in(Lp[-1], lo, hi)
    disc = resultant_d(Fp, bdd(Fp))
    if disc: bps += _poly_roots_in(disc, lo, hi)
    for G in extra:
        r = resultant_d(Fp, G)
        if r: bps += _poly_roots_in(r, lo, hi)
    if box is not None:
        for e in (lo, hi):
            bps += _poly_roots_in(at_d(Fp, e), lo, hi)
    if avoid_zero:
        bps += [Fr(0)] if (lo is None or lo < 0) and (hi is None or hi > 0) else []
        bps += _poly_roots_in(at_d(Fp, Fr(0)), lo, hi)
    for c0 in _samples(bps, lo, hi):
        if avoid_zero and c0 == 0: continue
        if _count_real_d(at_c(Fp, c0), lo, hi, avoid_zero) > 0: return True
    return False

def at_d(P, d0):
    """univariate polynomial in c of P(c, d0)"""
    return at_c(swap(P), d0)

# ----------------------------- one-variable images (exact), on R or an interval
def image_1d(n, d, lo=None, hi=None):
    """closure of { n(t)/d(t) : t in (lo, hi), d(t) != 0 } as signed closed intervals.
       n, d rational univariate polynomials, coprime."""
    from .algebraic import coef
    n, d = ptrim([coef(v) for v in n]), ptrim([coef(v) for v in d])
    if not d: return []
    if not n: return [(Fr(0), Fr(0))]
    g = pgcd(n, d)
    if len(g) > 1: n, d = pdivmod(n, g)[0], pdivmod(d, g)[0]
    if len(n) == 1 and len(d) == 1: v = n[0] / d[0]; return [(v, v)]
    R = lambda t: peval(n, t) / peval(d, t)
    crit = psub(pmul(pderiv(n), d), pmul(n, pderiv(d)))
    bps = set_unique(_poly_roots_in(crit, lo, hi) + _poly_roots_in(d, lo, hi))
    bps = sorted(bps, key=cmp_to_key(cmp))
    ends = ([lo] if lo is not None else [None]) + bps + ([hi] if hi is not None else [None])
    is_pole = lambda b: b is not None and sign_of(peval_any(d, b)) == 0
    dn, dd = len(n) - 1, len(d) - 1
    lim = Fr(0) if dn < dd else (n[-1] / d[-1] if dn == dd else None)
    out = []
    for L, Rr in zip(ends, ends[1:]):
        if L is None and Rr is None: m1, m2 = Fr(0), Fr(1)
        elif L is None: b = _bounds(Rr)[0] - 1; m1, m2 = b - 1, b
        elif Rr is None: b = _bounds(L)[1] + 1; m1, m2 = b, b + 1
        else:
            m2 = _rational_between(L, Rr); m1 = _rational_between(L, m2)
        inc = R(m2) > R(m1)
        def end(b, left):
            if b is None:
                if lim is not None: return lim
            elif not is_pole(b):
                return peval_any(n, b) / peval_any(d, b)
            if left: return NEG if inc else POS
            return POS if inc else NEG
        e1, e2 = end(L, True), end(Rr, False)
        out.append((e1, e2) if cmp(e1, e2) <= 0 else (e2, e1))
    return out

# ----------------------------- critical values of R = P/Q in the plane
_RECORD = []          # when a list is pushed here, critical points are recorded as (value, c, d)

def _add(vals, v, c0, d0):
    vals.append(v)
    if _RECORD: _RECORD[-1].append((v, c0, d0))

def critical_values(P, Q, box=None):
    """values R(c, d) at the real critical points (grad R = 0, Q != 0), exactly when the
       critical points are isolated with a linear fibre gcd; otherwise sampled exactly on the
       critical curve.  Returns (values, exact_flag)."""
    G1 = badd(bmul(bdc(P), Q), bscale(bmul(P, bdc(Q)), -1))
    G2 = badd(bmul(bdd(P), Q), bscale(bmul(P, bdd(Q)), -1))
    if not G1 and not G2: return [], True
    vals, exact = [], True
    lo, hi = (None, None) if box is None else box
    h = bgcd(G1, G2) if G1 and G2 else (G1 or G2)
    if h and max((k for _, k in h), default=0) + max((j for j, _ in h), default=0) > 0:
        # a curve of critical points (R constant along each of its components): sample it
        hd = in_d(h)
        bps = []
        if len(hd) > 1:
            bps += _poly_roots_in(hd[-1], lo, hi)
            disc = resultant_d(h, bdd(h))
            if disc: bps += _poly_roots_in(disc, lo, hi)
            bps += _poly_roots_in(resultant_d(h, Q), lo, hi) if resultant_d(h, Q) else []
            for c0 in _samples(bps, lo, hi):
                for d0 in _poly_roots_in(at_c(h, c0), lo, hi):
                    q = beval(Q, c0, d0)
                    if sign_of(q) != 0: _add(vals, beval(P, c0, d0) / q, c0, d0)
        else:
            for c0 in _poly_roots_in(hd[0], lo, hi):          # vertical critical lines
                for d0 in _samples([], lo, hi):
                    q = beval(Q, c0, d0)
                    if sign_of(q) != 0: _add(vals, beval(P, c0, d0) / q, c0, d0)
        G1 = bdivexact(G1, h) if G1 else G1
        G2 = bdivexact(G2, h) if G2 else G2
    if G1 and G2 and max(k for _, k in G1) > 0 and max(k for _, k in G2) > 0:
        res = resultant_d(G1, G2)
        if res:
            for c0 in _poly_roots_in(res, lo, hi):
                g = _gcd_any(at_c(G1, c0), at_c(G2, c0))
                g = ptrim(g)
                if len(g) == 2:
                    d0 = -g[0] / g[1]
                    if box is not None and not (cmp(d0, lo) > 0 and cmp(d0, hi) < 0): continue
                    q = beval(Q, c0, d0)
                    if sign_of(q) != 0: _add(vals, beval(P, c0, d0) / q, c0, d0)
                elif len(g) > 2:
                    try:
                        vals += _exact_fibre_values(P, Q, g, c0, lo, hi)
                    except (ArithmeticError, NotImplementedError):
                        exact = False                      # number field too large: 40 digits, flagged
                        vals += _approx_fibre_values(P, Q, g, c0, lo, hi)
    elif G1 or G2:
        # one gradient component vanishes identically: critical points = zeros of the other
        G = G1 or G2
        if max(k for _, k in G) == 0:
            for c0 in _poly_roots_in(in_d(G)[0], lo, hi):
                for d0 in _samples([], lo, hi):
                    q = beval(Q, c0, d0)
                    if sign_of(q) != 0: _add(vals, beval(P, c0, d0) / q, c0, d0)
    return vals, exact

def _gcd_any(a, b):
    """univariate gcd with rational or (single-field) algebraic coefficients"""
    a, b = ptrim(a), ptrim(b)
    while b:
        q = [Fr(0)] * max(len(a) - len(b) + 1, 1)
        a = list(a)
        while len(a) >= len(b) and a:
            k = len(a) - len(b); f = a[-1] / b[-1]
            for i, bv in enumerate(b): a[i + k] = a[i + k] - f * bv
            while a and (a[-1] == 0): a.pop()
        a, b = b, a
    return [v / a[-1] for v in a] if a else a

def _exact_fibre_values(P, Q, g, c0, lo, hi):
    """critical points over c0 with a nonlinear fibre g(d) = 0, exactly: the real roots of g
       (coefficients in Q(c0)) are real algebraic numbers, and R is evaluated in the compositum
       of the two number fields (roadmap #3)"""
    out = []
    for d0 in _real_roots(g):
        if lo is not None and not (cmp(d0, lo) > 0 and cmp(d0, hi) < 0): continue
        q = beval(Q, c0, d0)
        if sign_of(q) != 0: _add(out, beval(P, c0, d0) / q, c0, d0)
    return out

def _approx_fibre_values(P, Q, g, c0, lo, hi):
    """critical points over c0 with a nonlinear fibre: approximate d by bisection on a rational
       polynomial obtained by approximating the coefficients (flagged inexact by the caller)"""
    gr = [Fr(v.F.approx(v.c)) if isinstance(v, A) else Fr(v) for v in g]
    c0r = Fr(c0.F.approx(c0.c)) if isinstance(c0, A) else c0
    out = []
    for d0 in _poly_roots_in(gr, lo, hi):
        d0r = Fr(d0.F.approx(d0.c)) if isinstance(d0, A) else d0
        q = beval(Q, c0r, d0r)
        if q != 0: out.append(beval(P, c0r, d0r) / q)
    return out

# ----------------------------- the range of R = P/Q
def image_2d(P, Q, box=None, extra=(), torus=False):
    """closure of R(region) for R = P/Q as signed closed intervals of the extended line.
       region: the plane (box=None) or the open square box=(lo, hi)^2.
       extra: values known to lie on the boundary of what matters (they become candidates).
       Returns (intervals, exact_flag)."""
    P, Q = bclean(P), bclean(Q)
    if not Q: return [], True
    if not P: return [(Fr(0), Fr(0))], True
    g = bgcd(P, Q)
    if len(g) > 1 or any(k != (0, 0) for k in g):
        P, Q = bdivexact(P, g), bdivexact(Q, g)
    # constant?
    k0 = next(iter(Q)); r = P.get(k0, 0) / Q[k0]
    if bclean(badd(P, bscale(Q, -r))) == {}: return [(r, r)], True
    cand, exact = critical_values(P, Q, box)
    cand = list(cand) + [v for v in extra if not isinstance(v, _Inf)]
    lo, hi = (None, None) if box is None else box
    # boundary / asymptotic directions along the coordinate lines
    Pd, Qd = in_d(P), in_d(Q)
    at_inf = []
    if box is None:
        # d -> +-oo with c fixed, and c -> +-oo with d fixed: limits are rational functions
        for (PP, QQ) in ((P, Q), (swap(P), swap(Q))):
            a, b = in_d(PP), in_d(QQ)
            if len(a) == len(b):
                for (l, h) in image_1d(a[-1], b[-1]): cand += [x for x in (l, h) if not isinstance(x, _Inf)]
        # escaping to infinity in other ways: the fluid state of R at the points at infinity,
        # in the charts (1/u, 1/v), (1/u, v), (u, 1/v), by the two-variable engine
        at_inf, ok = _values_at_infinity(P, Q)
        exact = exact and ok
        for (l, h) in at_inf: cand += [x for x in (l, h) if not isinstance(x, _Inf)]
    else:
        for e in (lo, hi):
            for (PP, QQ) in ((at_d(P, e), at_d(Q, e)), (at_c(P, e), at_c(Q, e))):
                for (l, h) in image_1d(PP, QQ, lo, hi): cand += [x for x in (l, h) if not isinstance(x, _Inf)]
            for e2 in (lo, hi):
                q = beval(Q, e, e2)
                if q: cand.append(beval(P, e, e2) / q)
    cand = sorted(set_unique(cand), key=cmp_to_key(cmp))
    # test each open gap with one rational value
    gaps = []
    edges = [NEG] + cand + [POS]
    out = []
    for u, v in zip(edges, edges[1:]):
        if isinstance(u, _Inf) and isinstance(v, _Inf): w = Fr(0)
        elif isinstance(u, _Inf): w = _bounds(v)[0] - 1
        elif isinstance(v, _Inf): w = _bounds(u)[1] + 1
        else: w = _rational_between(u, v)
        F = badd(P, bscale(Q, -w))
        if has_real_points(F, box, extra=(Q,), avoid_zero=torus):
            out.append((u, v))
    return out + list(at_inf), exact

def _values_at_infinity(P, Q):
    """cluster values of R = P/Q as (c, d) escapes to infinity (charts at infinity).
       These are limits of values of R, so they belong to the closure of its range."""
    from .multivar import profile_polys
    degc = max(j for j, _ in list(P) + list(Q)); degd = max(k for _, k in list(P) + list(Q))
    out, ok = [], True
    for inv_c, inv_d in ((True, True), (True, False), (False, True)):
        def tr(M):
            return {((degc - j) if inv_c else j, (degd - k) if inv_d else k): v for (j, k), v in M.items()}
        Pn, Qn = tr(P), tr(Q)
        # cancel common powers of the chart coordinates
        mj = min(j for j, _ in list(Pn) + list(Qn)); mk = min(k for _, k in list(Pn) + list(Qn))
        Pn = {(j - mj, k - mk): v for (j, k), v in Pn.items()}; Qn = {(j - mj, k - mk): v for (j, k), v in Qn.items()}
        try:
            cl = profile_polys(Pn, Qn)["cluster"]
            out += list(cl.get("signed_components", []))
            ok = ok and cl.get("exact", False) and cl.get("complete", False)
        except (ArithmeticError, NotImplementedError, ZeroDivisionError):
            ok = False
    return out, ok
