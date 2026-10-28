"""Limits with the distance r = sqrt(x^2 + y^2): exact fluid states of f(x, y, r) at the origin.

Distances appear everywhere (potentials 1/r, norms, volatilities), and r is not a rational
function of x and y. FMSL handles it by the polar blow-up, done with rational functions only:

        x = r (1 - t^2)/(1 + t^2),   y = r 2t/(1 + t^2),   r > 0,   t in R u {oo}

maps (0, eps) x circle one-to-one onto the punctured disc of radius eps, so the fluid state of f
at the origin is exactly the set of values F(r, t) = f(x, y, r) approaches as r -> 0+, with t
anywhere on the circle (a compact set). Writing F = r^(a-b) N1/D1:

  * away from finitely many special t, F tends to the value of r^(a-b) N1(0,t)/D1(0,t): the closure
    of these values is an exact one-variable image (or {0}, or the signs of +-oo);
  * at each special t0 (a real root of N1(0, t) D1(0, t), rational or algebraic, and t = oo in the
    second chart t = 1/u), the local fluid state of F at (r, t) = (0, t0) with r > 0 comes from
    the two-variable engine.
The union of these finitely many closed sets is the fluid state. The point need not be the origin:
`radial_limit(text, point=(a, b))` uses r = sqrt((x-a)^2 + (y-b)^2).

    from fmsl.radial import radial_limit
    radial_limit("x*y/r")["set"]               # {0}
    radial_limit("(r - x)/y")["set"]           # [-oo, +oo]  (it is tan of half the angle)
    radial_limit("x^2/(r^2 + r^3)")["set"]     # [0, 1]
In expressions, write r or sqrt(x^2+y^2).
"""
import re
from fractions import Fraction as Fr
from math import comb


def _norm_text(text):
    s = text.replace(" ", "")
    for pat in ("sqrt(x^2+y^2)", "sqrt(y^2+x^2)", "(x^2+y^2)^(1/2)", "(y^2+x^2)^(1/2)"):
        s = s.replace(pat, "r")
    if "sqrt" in s or "^(1/2)" in s:
        raise ValueError("only the distance sqrt(x^2+y^2) is supported as a square root here")
    return re.sub(r"\br\b", "z", s)


def _ppow(p, n):
    out = [Fr(1)]
    for _ in range(n):
        q = [Fr(0)] * (len(out) + len(p) - 1)
        for i, a in enumerate(out):
            for j, b in enumerate(p):
                q[i + j] += a * b
        out = q
    return out


def _polar(P, m, chart):
    """substitute x, y, z(=r) in P (dict over (i, j, k)) and clear the denominator (1+t^2)^m;
       returns a dict over (r power, t power)"""
    if chart == 1:   one_minus, two_t, one_plus = [Fr(1), 0, Fr(-1)], [0, Fr(2)], [Fr(1), 0, Fr(1)]
    else:            one_minus, two_t, one_plus = [Fr(-1), 0, Fr(1)], [0, Fr(2)], [Fr(1), 0, Fr(1)]
    R = {}
    for (i, j, k), v in P.items():
        poly = _ppow(one_minus, i)
        for p in (_ppow(two_t, j), _ppow(one_plus, m - i - j)):
            q = [Fr(0)] * (len(poly) + len(p) - 1)
            for a, x in enumerate(poly):
                for b, y in enumerate(p):
                    q[a + b] += x * y
            poly = q
        for b, c in enumerate(poly):
            if c: R[(i + j + k, b)] = R.get((i + j + k, b), 0) + v * c
    return {key: c for key, c in R.items() if c != 0}


def _lowest_r(P):
    a = min(i for i, _ in P)
    return a, [P.get((a, b), Fr(0)) for b in range(max(b for i, b in P if i == a) + 1)]


def _shift_t(P, t0):
    """P(r, t0 + tau) as a dict over (r power, tau power)"""
    R = {}
    for (i, b), v in P.items():
        for k in range(b + 1):
            c = v * comb(b, k) * (t0 ** (b - k))
            if c != 0: R[(i, k)] = R.get((i, k), 0) + c
    return {key: c for key, c in R.items() if c != 0}


def radial_limit(text, point=(0, 0)):
    """the complete set of values f(x, y, r), r = sqrt((x-a)^2 + (y-b)^2), approaches as (x, y) -> (a, b)"""
    try:
        return _radial(text, point)
    except ArithmeticError as exc:
        return {"expression": text, "components": None, "set": None, "limit_exists": None, "limit": None,
                "proven": False, "notes": [f"not decided: {exc} [UNPROVEN]"]}


def _radial(text, point):
    from .space import parse3
    from .multivar import profile_polys
    from .cluster import _union_line, _real_roots, image_of_rational_function, POS, NEG, _describe_line
    from .algebraic import pgcd, pdivmod, ptrim
    s = _norm_text(text)
    a0, b0 = Fr(point[0]), Fr(point[1])
    if a0 or b0:
        s = re.sub(r"\bx\b", f"(x+({a0}))", s); s = re.sub(r"\by\b", f"(y+({b0}))", s)
    N, D = parse3(s)
    m = max(i + j for (i, j, k) in list(N) + list(D))
    items, exact, complete, notes = [], True, True, []
    Nc, Dc = _polar(N, m, 1), _polar(D, m, 1)
    if not Nc: return _result([(Fr(0), Fr(0))], True, True, text, ["f = 0 identically"])
    aN, N1 = _lowest_r(Nc); aD, D1 = _lowest_r(Dc)
    # generic directions
    g = pgcd(ptrim(N1), ptrim(D1))
    n1, d1 = ptrim(pdivmod(N1, g)[0]), ptrim(pdivmod(D1, g)[0])
    if aN > aD:
        items.append((Fr(0), Fr(0)))
    elif len(n1) <= 1 and len(d1) <= 1:
        val = n1[0] / d1[0]
        if aN == aD: items.append((val, val))
        else: items.append((POS, POS) if val > 0 else (NEG, NEG))
    else:
        img = image_of_rational_function(n1, d1)
        if aN == aD:
            items += img
        else:                                          # r^(negative) * values: only the signs survive
            from .cluster import cmp
            if any(cmp(hi, 0) > 0 for lo, hi in img): items.append((POS, POS))
            if any(cmp(lo, 0) < 0 for lo, hi in img): items.append((NEG, NEG))
    # special directions: roots of N1 * D1 in chart 1, and t = oo (u = 0) in chart 2
    from .cluster import cmp
    special = []
    for p in (ptrim(N1), ptrim(D1)):
        if len(p) > 1:
            for t0 in _real_roots(p):
                try:
                    dup = any(cmp(t0, u) == 0 for u in special)
                except ArithmeticError:
                    dup = False                        # different fields, too large to compare: keep both
                if not dup: special.append(t0)
    def local(P, Q, t0):
        nonlocal exact, complete
        Pn, Qn = _shift_t(P, t0), _shift_t(Q, t0)
        try:
            r = profile_polys(Pn, Qn, sides=(1,), y_axis=False)
        except ArithmeticError:                        # number field too large: say so, claim nothing
            complete = False
            return []
        cl = r["cluster"]
        exact = exact and bool(cl.get("exact")); complete = complete and bool(cl.get("complete", True))
        return list(cl.get("signed_components", []))
    for t0 in special:
        items += local(Nc, Dc, t0)
    N2, D2 = _polar(N, m, 2), _polar(D, m, 2)
    items += local(N2, D2, Fr(0))
    comps = _union_line(items)
    return _result(comps, exact, complete, text, [f"special directions: {len(special)} + t = oo"])


def _result(comps, exact, complete, text, notes):
    from .cluster import _describe_line, _Inf
    single = len(comps) == 1 and comps[0][0] == comps[0][1] and not isinstance(comps[0][0], _Inf)
    return {"expression": text, "components": comps, "set": _describe_line(comps),
            "limit_exists": single, "limit": comps[0][0] if single else None,
            "proven": exact and complete, "notes": notes}
