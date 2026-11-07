"""
fmsl.celedon -- Celedon points hidden inside the zoom (parametric Newton-Puiseux).

fmsl.moving finds candidate Celedon points from the TOP-LEVEL Newton data of a family
f_t = N(x, y, t) / D(x, y, t).  A change of character that happens only after the engine
zooms in along a curve y = c0 |x|^a + ... is found there by sampling, and reported as an
approximate "change near ...".

This module runs the zoom with t kept as a SYMBOL.  At every tie weight a and side of x, the
zoom points are the repeated or shared real roots c0 of the initial forms.  When such a root
is a rational function of t, c0(t) = P(t)/Q(t), it is found exactly (power series at a
generic time, Pade approximant, then an exact symbolic check), and the substitution

        x = sigma * s^q,   y = s^p * (P(t) + w) / Q(t)       (a = p/q, s > 0)

gives a new polynomial family in (s, w, t) whose own top-level candidates are computed by
fmsl.moving.candidate_times -- and then zoomed again, to a fixed depth.

What is proven (see CELEDON_DEEP.md for the statement and proof):
  * every time returned is an exact real algebraic number;
  * the returned list contains every time at which the engine's zoom tree changes, along all
    zoom branches whose roots are rational functions of t, down to the given depth.
What is NOT covered, and is reported in 'unsupported':
  * zoom roots that are irrational algebraic functions of t (e.g. c0 = sqrt(t)) -- those
    branches fall back to fmsl.moving's sampling, exactly as before;
  * zoom trees deeper than the depth limit.
Extra candidates are harmless: fmsl.moving keeps a candidate only where the type of the
snapshot actually changes.
"""
from fractions import Fraction as Fr
from math import comb

from .algebraic import ptrim, pmul, padd, pgcd, pdivmod, squarefree
from .cluster import cmp, _real_roots
from .plane import bgcd, bdivexact, resultant_d
from .moving import (_coeff_polys, _tie_weights, _initial_form, _dc, _bmul, _cdeg,
                     candidate_times)


# ----------------------------------------------------------------- small helpers
def _in_c(P):
    """bivariate dict {(t power, c power)} -> list over c of polynomials in t"""
    n = _cdeg(P)
    L = [[] for _ in range(n + 1)]
    for (k, j), v in P.items():
        r = L[j]
        while len(r) <= k: r.append(Fr(0))
        r[k] += v
    return [ptrim(r) for r in L]


def _from_c(L):
    return {(k, j): v for j, r in enumerate(L) for k, v in enumerate(r) if v}


def _peval(p, x):
    acc = Fr(0)
    for v in reversed(p): acc = acc * x + v
    return acc


def _pshift(p, t0):
    """p(u + t0) as a polynomial in u"""
    out = []
    for k, v in enumerate(p):
        if not v: continue
        for e in range(k + 1):
            c = v * comb(k, e) * (t0 ** (k - e) if k - e else 1)
            while len(out) <= e: out.append(Fr(0))
            out[e] += c
    return ptrim(out)


def _ppow(p, e):
    r = [Fr(1)]
    for _ in range(e): r = pmul(r, p)
    return r


def _sqfree_c(g):
    """squarefree part of g in c over Q(t)"""
    if _cdeg(g) <= 0: return g
    G = bgcd(g, _dc(g))
    if _cdeg(G) <= 0: return g
    return bdivexact(g, G)


def _drop_c_factor(g):
    """remove the factor c (the root c0 = 0 belongs to other weights)"""
    m = min((j for (_, j) in g), default=0)
    return {(k, j - m): v for (k, j), v in g.items()} if m else g


def _spec(L, t0):
    return ptrim([_peval(r, t0) for r in L])


# ------------------------------------------------- roots of g(c, t) in Q(t), exactly
def _series_root(Lu, r0, order):
    """power-series root c(u) = r0 + r1 u + ... of sum_j Lu[j](u) c^j (simple root at u=0)"""
    def series_val(cs):
        # value of g at c = cs (truncated series), truncated to `order` terms
        acc = [Fr(0)] * order
        for j in reversed(range(len(Lu))):
            # acc = acc * cs + Lu[j]
            new = [Fr(0)] * order
            for a, x in enumerate(acc):
                if x:
                    for b, y in enumerate(cs):
                        if a + b < order and y: new[a + b] += x * y
            for a, x in enumerate(Lu[j][:order]): new[a] += x
            acc = new
        return acc
    dg = sum(j * _peval(Lu[j], 0) * r0 ** (j - 1) for j in range(1, len(Lu)))
    if dg == 0: return None
    cs = [r0]
    for n in range(1, order):
        val = series_val(cs + [Fr(0)])
        cs.append(-val[n] / dg)
    return cs


def _solve(M, b):
    """exact Gaussian elimination; None if singular"""
    n = len(M)
    A = [row[:] + [b[i]] for i, row in enumerate(M)]
    for col in range(n):
        piv = next((r for r in range(col, n) if A[r][col] != 0), None)
        if piv is None: return None
        A[col], A[piv] = A[piv], A[col]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col] / A[col][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[col])]
    return [A[i][n] / A[i][i] for i in range(n)]


def _pade(cs, L):
    """[L/L] Pade approximant P/Q of the series cs (needs 2L+1 terms), Q(0) = 1"""
    if L == 0: return [cs[0]], [Fr(1)]
    M = [[cs[n - j] if n - j >= 0 else Fr(0) for j in range(1, L + 1)] for n in range(L + 1, 2 * L + 1)]
    b = [-cs[n] for n in range(L + 1, 2 * L + 1)]
    q = _solve(M, b)
    if q is None: return None
    Q = [Fr(1)] + q
    P = [sum((Q[j] * cs[n - j] for j in range(0, min(n, L) + 1)), Fr(0)) for n in range(L + 1)]
    return ptrim(P), ptrim(Q)


def _vanishes(Lt, P, Q):
    """is sum_j Lt[j](t) (P/Q)^j identically zero?"""
    n = len(Lt) - 1
    acc = []
    for j, r in enumerate(Lt):
        if r: acc = padd(acc, pmul(r, pmul(_ppow(P, j), _ppow(Q, n - j))))
    return not ptrim(acc)


def rational_roots_in_t(g, tries=(0, 1, -1, 2, -2, Fr(1, 2), 3, -3, Fr(1, 3), 5, 7, -5)):
    """real roots c0(t) of g(c, t) that are rational functions of t, as (P, Q) polynomial
       pairs with Q != 0; and whether every real root branch was accounted for (False if some
       branch is not a rational function of t)."""
    g = _sqfree_c(_drop_c_factor(g))
    if _cdeg(g) <= 0: return [], True
    Lt = _in_c(g)
    dt = max(len(r) for r in Lt) - 1
    for t0 in tries:
        t0 = Fr(t0)
        s = _spec(Lt, t0)
        if len(s) - 1 != len(Lt) - 1: continue                    # leading coefficient vanishes
        if len(ptrim(pgcd(s, ptrim([i * s[i] for i in range(1, len(s))])))) > 1: continue
        Lu = [_pshift(r, t0) for r in Lt]
        found, complete = [], True
        for r0 in _real_roots(s):
            if not isinstance(r0, Fr):
                complete = False                                   # irrational at a rational time
                continue
            cs = _series_root(Lu, r0, 2 * dt + 2)
            ok = False
            for L in range(0, dt + 1):
                pq = _pade(cs, L)
                if pq is None: continue
                Pu, Qu = pq
                P, Q = _pshift(Pu, -t0), _pshift(Qu, -t0)
                if _vanishes(Lt, P, Q):
                    found.append((P, Q)); ok = True; break
            if not ok: complete = False
        # Every root in Q(t) is rational at t0, so it is among the real roots above.  The roots
        # NOT found (complex at t0, or real but not a rational function of t) may still be real
        # zoom points at other times, so the branch is complete only if all deg_c roots were found.
        if len(found) != len(Lt) - 1: complete = False
        return found, complete
    return [], False


# ------------------------------------------------------------------ the substitution
def zoom_substitute(N, D, a, sigma, P, Q):
    """N, D in x^i y^j t^k  ->  N1, D1 in s^i w^j t^k under x = sigma s^q,
       y = s^p (P(t) + w) / Q(t), both multiplied by Q^K (same K, so the ratio is unchanged)"""
    p, q = a.numerator, a.denominator
    K = max([j for (_, j, _) in N] + [j for (_, j, _) in D])
    Ppows = [_ppow(P, e) for e in range(K + 1)]
    Qpows = [_ppow(Q, e) for e in range(K + 1)]

    def sub(F):
        out = {}
        for (i, j, k), v in F.items():
            si = q * i + p * j
            sg = sigma ** i
            for e in range(j + 1):                                 # (P + w)^j = sum C(j,e) P^(j-e) w^e
                tp = pmul(Ppows[j - e], Qpows[K - j])
                c = v * sg * comb(j, e)
                for kk, u in enumerate(tp):
                    if u:
                        key = (si, e, k + kk)
                        out[key] = out.get(key, 0) + c * u
        return {kk: vv for kk, vv in out.items() if vv != 0}
    return sub(N), sub(D)


# ------------------------------------------------------------------ the main search
def _collision_polys(fN, fD):
    """polynomials in t whose roots are where ANY two roots of fN*fD collide, escape to
       c = 0 or infinity -- the generic-validity conditions of the symbolic zoom"""
    out = []
    h = _sqfree_c(_drop_c_factor(_bmul(fN, fD)))
    if _cdeg(h) > 0:
        r = ptrim(resultant_d(h, _dc(h)))
        if len(r) > 1: out.append(r)
        Lc = _in_c(h)
        for r in (Lc[-1], Lc[0]):
            if len(r) > 1: out.append(r)
    return out


def deep_candidate_times(N, D, depth=3, _path=(), _collect=None):
    """exact candidate Celedon points, including those that arise inside the zoom.
       Returns (candidates, unsupported): candidates as in fmsl.moving.candidate_times
       [(time, [reasons])]; unsupported lists zoom branches not covered (irrational in t,
       or beyond the depth limit), so the caller knows where sampling is still the only check."""
    polys, unsupported = [], []
    CN, CD = _coeff_polys(N), _coeff_polys(D)
    if not CN or not CD: return [], unsupported
    where = " > ".join(_path) if _path else "top level"
    for a in _tie_weights(set(CN) | set(CD)):
        for sigma in (1, -1):
            fN, fD = _initial_form(CN, a, sigma), _initial_form(CD, a, sigma)
            if not fN or not fD: continue
            side = "x>0" if sigma > 0 else "x<0"
            for r in _collision_polys(fN, fD):
                polys.append((r, f"zoom roots collide (weight {a}, {side}, {where})"))
            zoomers = []
            for P_, why in ((fN, "repeated numerator root"), (fD, "repeated denominator root")):
                if _cdeg(P_) > 0:
                    G = bgcd(P_, _dc(P_))
                    if _cdeg(G) > 0: zoomers.append(G)
            G = bgcd(fN, fD)
            if _cdeg(G) > 0: zoomers.append(G)
            seen = []
            for G in zoomers:
                roots, complete = rational_roots_in_t(G)
                if not complete:
                    unsupported.append(f"zoom root not rational in t (weight {a}, {side}, {where})")
                for P, Q in roots:
                    if not P: continue                              # c0 = 0: other weights
                    key = (tuple(P), tuple(Q))
                    if key in seen: continue
                    seen.append(key)
                    if len(Q) > 1:
                        polys.append((Q, f"zoom root has a pole (weight {a}, {side}, {where})"))
                    label = f"y = ({_ptxt(P)})/({_ptxt(Q)})*|x|^{a} ({side})" if len(Q) > 1 \
                        else f"y = ({_ptxt(P)})*|x|^{a} ({side})"
                    if depth <= 0:
                        unsupported.append(f"depth limit at {label}")
                        continue
                    N1, D1 = zoom_substitute(N, D, a, sigma, P, Q)
                    if _collect is not None: _collect.append((N1, D1))
                    for t, why in candidate_times(N1, D1):
                        for w in why:
                            polys.append(("root", t, f"{w} [inside zoom {label}]"))
                    sub, uns = deep_candidate_times(N1, D1, depth - 1, _path + (label,), _collect)
                    for t, why in sub:
                        for w in why: polys.append(("root", t, w))
                    unsupported.extend(uns)
    found = []
    def add(r, why):
        for item in found:
            if cmp(item[0], r) == 0: item[1].add(why); return
        found.append([r, {why}])
    for item in polys:
        if item[0] == "root":
            add(item[1], item[2])
        else:
            L, why = item
            sq = squarefree(L)
            if len(sq) <= 1: continue
            for r in _real_roots(sq): add(r, why)
    from functools import cmp_to_key
    found.sort(key=cmp_to_key(lambda u, v: cmp(u[0], v[0])))
    return [(r, sorted(w)) for r, w in found], sorted(set(unsupported))


def all_candidate_times(N, D, depth=3, zoomed=None):
    """top-level candidates (fmsl.moving) merged with the deep ones.  If `zoomed` is a list,
       the zoomed families (N1, D1) are appended to it (same ratio as N/D, new coordinates)."""
    top = candidate_times(N, D)
    deep, unsupported = deep_candidate_times(N, D, depth, _collect=zoomed)
    found = [[r, set(w)] for r, w in top]
    for r, w in deep:
        for item in found:
            if cmp(item[0], r) == 0: item[1].update(w); break
        else:
            found.append([r, set(w)])
    from functools import cmp_to_key
    found.sort(key=cmp_to_key(lambda u, v: cmp(u[0], v[0])))
    return [(r, sorted(w)) for r, w in found], unsupported


def _ptxt(p):
    terms = []
    for k, v in enumerate(p):
        if not v: continue
        terms.append(f"{v}" if k == 0 else (f"{v}*t" if k == 1 else f"{v}*t^{k}"))
    return " + ".join(terms) or "0"
