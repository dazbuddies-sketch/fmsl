"""Exact global minimum and maximum of a function of two variables over the whole plane, and
whether it is attained.

    from fmsl.optimize import global_min
    global_min("x^2 + (x*y - 1)^2")
        -> infimum 0, attained False: approached only as (x, y) -> infinity (e.g. x = 1/y, y -> oo)
    global_min("x^4 + y^4 - 2*x*y")
        -> infimum -1/2, attained True at (x, y) = (-0.707107, -0.707107) and (0.707107, 0.707107)

For f = P/Q with P, Q polynomials and Q without real zeros (in particular every polynomial):

  * the infimum is the lower end of the exact range of f over the plane (fmsl.plane.image_2d:
    critical values, values approached at infinity, and one exact test value per gap);
  * if the infimum m is attained at a point p, then p is a critical point (grad f = 0), so m is a
    critical value; conversely a critical value equal to m is attained there. So:
        attained  <=>  m is one of the exact critical values;
  * otherwise m is approached only along sequences that leave every bounded set (on a bounded
    closed set f attains its minimum), that is, at infinity.
Everything is exact: rationals and real algebraic numbers, compared by Sturm sequences.
"""
from functools import cmp_to_key

from .cluster import cmp, _Inf
from .multivar import parse_rational
from . import plane


def _fmt(v):
    from fractions import Fraction as Fr
    if isinstance(v, _Inf): return "+oo" if v.sgn > 0 else "-oo"
    if not isinstance(v, Fr):
        r = Fr(float(v)).limit_denominator(10 ** 6)          # is it rational? decided exactly
        if cmp(v, r) != 0: return f"{float(v):.6g} (algebraic)"
        v = r
    return str(v.numerator) if v.denominator == 1 else f"{v.numerator}/{v.denominator}"


def _extreme(text, sense):
    N, D = parse_rational(text)
    if sense < 0: N = {k: -v for k, v in N.items()}                 # max f = -min(-f)
    Pb, Qb = plane.bclean(N), plane.bclean(D)
    if not Qb: raise ZeroDivisionError("denominator is identically zero")
    if len(Qb) > 1 or (0, 0) not in Qb:
        if plane.has_real_points(Qb, None):
            raise ValueError("the denominator vanishes somewhere in the plane; only functions defined "
                             "everywhere are handled (every polynomial is)")
    plane._RECORD.append([])
    try:
        rng, exact = plane.image_2d(Pb, Qb)
        crit = plane._RECORD[-1]
    finally:
        plane._RECORD.pop()
    lows = [lo for lo, hi in rng]
    m = min(lows, key=cmp_to_key(cmp))
    sgn = -1 if sense < 0 else 1
    neg = lambda v: (_Inf(-v.sgn) if isinstance(v, _Inf) else -v) if sense < 0 else v
    word = "minimum" if sense > 0 else "maximum"
    inf_word = "infimum" if sense > 0 else "supremum"
    if isinstance(m, _Inf):
        return {"value": _fmt(neg(m)), "attained": False, "points": [], "proven": exact,
                "text": f"unbounded: the {inf_word} is {_fmt(neg(m))}" + ("" if exact else " [not fully exact]")}
    pts = [(c, d) for v, c, d in crit if cmp(v, m) == 0]
    if pts:
        uniq = []
        for p in pts:
            if not any(cmp(p[0], q[0]) == 0 and cmp(p[1], q[1]) == 0 for q in uniq): uniq.append(p)
        where = " and ".join(f"({_fmt(c)}, {_fmt(d)})" for c, d in uniq[:6]) + (" ..." if len(uniq) > 6 else "")
        return {"value": _fmt(neg(m)), "value_exact": neg(m), "attained": True,
                "points": [(_fmt(c), _fmt(d)) for c, d in uniq], "proven": exact,
                "text": f"the global {word} is {_fmt(neg(m))}, attained, for example, at (x, y) = {where}"
                        + (" [proven]" if exact else " [not fully exact]")}
    return {"value": _fmt(neg(m)), "value_exact": neg(m), "attained": False, "points": [], "proven": exact,
            "text": f"the {inf_word} is {_fmt(neg(m))} but it is NOT attained: it is approached only as "
                    f"(x, y) -> infinity, and no point of the plane reaches it"
                    + (" [proven]" if exact else " [not fully exact]")}


def global_min(text):
    """exact infimum of f(x, y) over the plane, and whether (and where) it is attained"""
    return _extreme(text, 1)


def global_max(text):
    """exact supremum of f(x, y) over the plane, and whether (and where) it is attained"""
    return _extreme(text, -1)
