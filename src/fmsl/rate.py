"""The rate: how fast a limit is approached (roadmap #4).

If f(x, y) -> L at the origin, the rate exponent is the best theta in

        |f(x, y) - L|  <=  C * r^theta     near the origin,   r = sqrt(x^2 + y^2).

It is rational, it is attained along an arc (so it cannot be improved), and it equals the
smallest value, over all arcs gamma(t) -> 0, of

        ord_t (f - L)(gamma(t))  /  ord_t |gamma(t)|.

That is a Lojasiewicz-type exponent. FMSL computes it exactly by Newton-Puiseux zooming:

  * along arcs y = c |x|^b (x > 0 and x < 0, every b > 0) and along the axes, the order of
    g = f - L = (N - L D)/D is ord N' - ord D for all but finitely many c: a piecewise linear
    function of b, so its minimum is at a corner of the Newton polygons (or at b = 1);
  * the exceptional c are the real roots of the initial form of D (where D is smaller than
    usual, so g is larger): there FMSL zooms, y = |x|^b (c + y1), and repeats. Roots of the
    initial form of N only make g smaller, so they cannot lower the minimum.
Since f has a limit, D has no real zero curve through the origin, so the zooming stops; a
depth cap makes it report UNPROVEN instead of looping on inputs it does not expect.

    from fmsl.rate import rate
    rate("x^2*y/(x^2+y^2)")["exponent"]       # Fraction(1, 1): |f| <= C r
    rate("x^4*y/(x^6+y^2)")["exponent"]       # 1, attained along y = c |x|^3

One variable: the rate is the order of the first nonzero term of f - L (also computed here).
"""
from fractions import Fraction as Fr
from math import gcd

from .multivar import parse_rational, padd, pmul, profile

MAX_DEPTH = 12


def _limit_value(text):
    from .cluster import _Inf
    cl = profile(text)["cluster"]
    comps = cl.get("signed_components") or []
    if not (cl.get("exact") and cl.get("complete", True)) or len(comps) != 1: return None
    a, b = comps[0]
    if isinstance(a, _Inf) or a != b: return None
    return a


def _w(P, b):
    return min(i + b * j for (i, j) in P)


def _edge_slopes(P):
    """slopes b > 0 at which the minimum of i + b j over the terms of P is attained twice"""
    out = set()
    pts = list(P)
    for (i1, j1) in pts:
        for (i2, j2) in pts:
            if j2 > j1 and i1 > i2:
                b = Fr(i1 - i2, j2 - j1)
                m = _w(P, b)
                if sum(1 for (i, j) in pts if i + b * j == m) >= 2: out.add(b)
    return out


def _initial_poly(P, b):
    """the initial form of P for weight (1, b), as a coefficient list in c (low to high)"""
    m = _w(P, b)
    terms = {j: v for (i, j), v in P.items() if i + b * j == m}
    return [terms.get(j, 0) for j in range(max(terms) + 1)]


def _zoom(P, p, q, c):
    """substitute u = v^q, w = v^p (c + w'), and divide out nothing (orders are tracked)"""
    from math import comb
    R = {}
    for (i, j), v in P.items():
        e = q * i + p * j
        ck = 1
        powers = [1]
        for _ in range(j): powers.append(powers[-1] * c)
        for k in range(j + 1):                          # (c + w')^j = sum C(j,k) c^(j-k) w'^k
            coef = v * comb(j, k) * powers[j - k]
            if coef != 0:
                R[(e, k)] = R.get((e, k), 0) + coef
    return {k: v for k, v in R.items() if v != 0}


def _ctext(c):
    if isinstance(c, Fr): return str(c.numerator) if c.denominator == 1 else f"{c.numerator}/{c.denominator}"
    return f"{float(c):.6g}(algebraic)"


def _etext(e):
    return str(e.numerator) if e.denominator == 1 else f"({e.numerator}/{e.denominator})"


class _Unproven(Exception): pass


def _explore(P, Q, S, depth, sigma, terms, E, Qtot, best, top):
    """P = numerator of f - L, Q = denominator, in (t, w) with t > 0; the arc so far is
       y = sum(terms) + |x|^E * w, |x| = t^Qtot; S = ord_t |gamma| (None at the top: min(1, b))."""
    from .cluster import _real_roots
    if depth > MAX_DEPTH: raise _Unproven("zoom deeper than the cap")
    side = "x>0" if sigma > 0 else "x<0"
    def arc(extra):
        parts = [f"{_ctext(c)}*|x|^{_etext(e)}" for c, e in terms] + ([extra] if extra else [])
        return "y = " + (" + ".join(parts) if parts else "0") + f"  ({side})"
    def offer(val, where, generic=None):
        if best["value"] is None or val < best["value"]:
            best["value"], best["arc"] = val, where
            best["arc_data"] = {"sigma": sigma, "terms": list(terms), "generic_exponent": generic}
    slopes = _edge_slopes(P) | _edge_slopes(Q) | ({Fr(1)} if top else set())
    for b in sorted(slopes):
        s = min(Fr(1), b) if top else S
        offer((_w(P, b) - _w(Q, b)) / s, arc(f"c*|x|^{_etext(E + b / Qtot)}") + "  for all but finitely many c",
              E + b / Qtot)
    # the exact arc w = 0 (b -> infinity)
    P0 = {i: v for (i, j), v in P.items() if j == 0}; Q0 = {i: v for (i, j), v in Q.items() if j == 0}
    if Q0 and P0:
        offer(Fr(min(P0) - min(Q0)) / (1 if top else S), arc(""))
    elif not Q0:
        raise _Unproven("the denominator vanishes along a curve")
    # zoom where the initial form of the denominator has a real root c != 0
    for b in sorted(_edge_slopes(Q)):
        ini = _initial_poly(Q, b)
        try:
            roots = [c for c in _real_roots(ini) if c != 0]
        except ArithmeticError as exc:
            raise _Unproven(f"algebraic numbers too large ({exc})")
        p, q = b.numerator, b.denominator
        for c in roots:
            S2 = (min(p, q) if top else S * q)
            _explore(_zoom(P, p, q, c), _zoom(Q, p, q, c), S2, depth + 1, sigma,
                     terms + [(c, E + b / Qtot)], E + b / Qtot, Qtot * q, best, False)


def rate(text):
    """The exact rate exponent theta of a two-variable rational function with a limit L at the
       origin: |f - L| <= C r^theta near 0, attained along an arc (so theta is optimal)."""
    N, D = parse_rational(text)
    L = _limit_value(text)
    if L is None:
        return {"limit": None, "exponent": None, "proven": True,
                "text": "no limit at the origin (or not proven): the rate is defined when a limit exists"}
    G = padd(N, {k: -L * v for k, v in D.items()})
    if not G:
        return {"limit": L, "exponent": None, "proven": True, "arc": None,
                "text": f"f = {L} identically near the origin (where defined)"}
    best = {"value": None, "arc": None}
    try:
        for sigma in (1, -1):
            P = {(i, j): v * (sigma ** i) for (i, j), v in G.items()}
            Q = {(i, j): v * (sigma ** i) for (i, j), v in D.items()}
            _explore(P, Q, None, 0, sigma, [], Fr(0), 1, best, True)
        for tau in (1, -1):                             # the y-axis: x = 0, y = tau t
            P0 = {j: v * tau ** j for (i, j), v in G.items() if i == 0}
            Q0 = {j: v * tau ** j for (i, j), v in D.items() if i == 0}
            if Q0 and P0 and (best["value"] is None or Fr(min(P0) - min(Q0)) < best["value"]):
                best["value"], best["arc"] = Fr(min(P0) - min(Q0)), f"x = 0  ({'y>0' if tau > 0 else 'y<0'})"
                best["arc_data"] = {"y_axis": tau}
    except _Unproven as exc:
        return {"limit": L, "exponent": None, "proven": False, "text": f"rate not determined: {exc} [UNPROVEN]"}
    th = best["value"]
    Ltxt = _ctext(L) if isinstance(L, Fr) else str(L)
    return {"limit": L, "exponent": th, "proven": True, "arc": best["arc"], "arc_data": best.get("arc_data"),
            "text": f"|f - {Ltxt}| <= C * r^{_etext(th)} near the origin, and no larger exponent works: "
                    f"along {best['arc']} the difference is of exact order r^{_etext(th)} [proven]"}


def rate_1d(expr, var="x"):
    """One variable: the order of the first nonzero term of f - L (same idea, one direction)."""
    from .compiler import Parser, build
    from .verify import limit_at_origin
    r = limit_at_origin(expr, [var])
    if not r["limit_exists"]:
        return {"limit": None, "exponent": None, "proven": True, "text": "no limit at 0"}
    from .engine import Essential
    try:
        vv = build(Parser(f"({expr}) - ({r['limit']})").parse(), var).val()
    except Essential:
        return {"limit": r["limit"], "exponent": None, "proven": False,
                "text": "the rate is not decided for this expression (an essential factor) [UNPROVEN]"}
    if vv[0] == "zero":
        return {"limit": r["limit"], "exponent": None, "proven": bool(vv[1]), "text": "f - L is identically zero"}
    k = Fr(vv[1])
    return {"limit": r["limit"], "exponent": k, "proven": bool(vv[2]),
            "text": f"|f - {r['limit']}| ~ C*|{var}|^{_etext(k)} as {var} -> 0"}
