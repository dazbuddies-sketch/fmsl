"""Limits at infinity and limits from one side, in two and three variables, by exact changes of
variables that bring them back to a limit at the origin.

At infinity. The inversion q -> q/|q|^2 is a bijection of the punctured plane (or space) onto
itself with |q/|q|^2| = 1/|q|. So (x, y) -> oo (|(x, y)| -> oo, in any manner) exactly when the
new point tends to 0, and the set of values f approaches at infinity equals the set of values of

    g(u, v) = f(u/(u^2+v^2), v/(u^2+v^2))         at (u, v) -> 0,

a rational function when f is. A distance sqrt(Q(x, y)), Q a quadratic form, becomes
sqrt(Q(u, v))/(u^2+v^2), because Q(q/|q|^2) = Q(q)/|q|^4; the distance machinery then removes the
remaining root.

    limit_at("x*y/(x^2+y^2)", "oo")["set"]          # [-1/2, 1/2]
    limit_at("(x^2+y)/(x^2+y^2)", "oo")["set"]      # [0, 1]

From one side. x -> a from above (x > a) is x = a + u^2 with u -> 0, u != 0: the map u -> a + u^2
is onto (a, a + eps) and u -> 0 exactly when x -> a. From below, x = a - u^2. Several variables
may be restricted at once (a quadrant: x > 0 and y > 0).

    limit_at("x*y/(x^2+y^2)", [0, 0], side={"x": "+", "y": "+"})["set"]   # [0, 1/2]
"""
import re
from fractions import Fraction as Fr

__all__ = ["at_infinity", "one_sided", "RegionError"]


class RegionError(ValueError):
    pass


def _calls(expr):
    return set(re.findall(r"([A-Za-z_]\w*)\s*\(", expr))


def at_infinity(expr, vs):
    """the expression g, in the same variable names, whose fluid state at the origin is the fluid
       state of expr as |(variables)| -> oo"""
    from .distance import _sqrt_spans, _half_powers, _poly
    s = _half_powers(expr.replace(" ", ""))
    other = _calls(s) - {"sqrt"}
    if other:
        raise RegionError(f"at infinity FMSL handles rational functions and distances; {sorted(other)[0]} "
                          "has no limit set FMSL can compute there")
    r2 = "(" + "+".join(f"{v}^2" for v in vs) + ")"
    spans = _sqrt_spans(s)
    out, last, roots = [], 0, []
    for k, (a, b, arg) in enumerate(spans):
        P = _poly(arg, vs)
        if any(sum(e) != 2 for e in P):
            raise RegionError("at infinity a square root must be of a quadratic form, such as sqrt(x^2 + 2*y^2)")
        out.append(s[last:a]); out.append(f"__R{k}__"); last = b; roots.append(arg)
    out.append(s[last:]); t = "".join(out)
    tmp = {v: f"__V{k}__" for k, v in enumerate(vs)}
    for v, tok in tmp.items(): t = re.sub(rf"\b{re.escape(v)}\b", tok, t)
    for v, tok in tmp.items(): t = t.replace(tok, f"({v}/{r2})")
    for k, arg in enumerate(roots):                     # sqrt(Q(x)) -> sqrt(Q(u))/|u|^2
        t = t.replace(f"__R{k}__", f"(sqrt({arg})/{r2})")
    return t


def one_sided(expr, vs, point, side):
    """the expression g whose fluid state at the origin is the fluid state of expr at `point`
       with each variable in `side` restricted to one side ("+": above, "-": below)"""
    if "sqrt" in expr:
        from .distance import reduce_distance, DistanceError
        try:
            has_distance = reduce_distance(expr, vs, point) is not None
        except DistanceError:
            has_distance = True
        if has_distance: raise RegionError("one-sided limits are not combined with distances yet")
    for v, sd in side.items():
        if v not in vs: raise RegionError(f"unknown variable {v!r} in side")
        if sd not in ("+", "-"): raise RegionError("side values are '+' (from above) or '-' (from below)")
    t = expr
    tmp = {v: f"__S{k}__" for k, v in enumerate(vs)}
    for v, tok in tmp.items(): t = re.sub(rf"\b{re.escape(v)}\b", tok, t)
    for (v, tok), p in zip(tmp.items(), point):
        p = Fr(p)
        if v in side:
            core = f"{v}^2" if side[v] == "+" else f"(-{v}^2)"
        else:
            core = v
        t = t.replace(tok, f"({core})" if p == 0 else f"(({p})+{core})")
    return t
