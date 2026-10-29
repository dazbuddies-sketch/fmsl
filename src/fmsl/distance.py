"""Distances in two and three variables, reduced exactly to the rational engines.

A formula may contain a Euclidean distance d = sqrt((x-a)^2 + (y-b)^2 [+ (z-c)^2]) to the point
where the limit is taken. FMSL removes the square root with a polynomial change of variables that
squares the distance away (the real form of z -> z^2 and its three-dimensional cousin):

    two variables:    x = u^2 - v^2,          y = 2uv,                  d = u^2 + v^2
    three variables:  x = u^2 - v^2 - w^2,    y = 2uv,    z = 2uw,      d = u^2 + v^2 + w^2

Indeed (u^2 - v^2 - w^2)^2 + (2uv)^2 + (2uw)^2 = (u^2 + v^2 + w^2)^2. The map is onto (every point
near the origin has a preimage), and q -> 0 exactly when its image tends to 0, because |q|^2 is the
distance. So the set of values f approaches at the point equals the set of values the rational
function g(u, v, w) = f(x(q), y(q), z(q), d(q)) approaches at q = 0: limits along sequences
correspond one-to-one in both directions. The rational two- and three-variable engines then give
the answer, proven, with nothing estimated.

A distance in only some of the variables (sqrt(x^2 + y^2) with a third variable z, a cylinder) works
the same way; the other variables are left alone.

    from fmsl.verify import limit_at_origin
    limit_at_origin("x*y*z/sqrt(x^2+y^2+z^2)^3")["set"]     # [-1/(3 sqrt 3), 1/(3 sqrt 3)]
"""
import re
from fractions import Fraction as Fr

__all__ = ["reduce_distance", "DistanceError"]


class DistanceError(ValueError):
    pass


def _sqrt_spans(s):
    """(start, end, argument) of each sqrt( ... ) in s, outermost first"""
    out, i = [], 0
    while True:
        k = s.find("sqrt(", i)
        if k < 0: return out
        depth, j = 0, k + 4
        while j < len(s):
            if s[j] == "(": depth += 1
            elif s[j] == ")":
                depth -= 1
                if depth == 0: break
            j += 1
        if depth: raise DistanceError("unbalanced parentheses after sqrt")
        out.append((k, j + 1, s[k + 5:j]))
        i = j + 1


def _poly(text, vs):
    """the polynomial `text` in the variables vs as {exponent tuple: Fraction}"""
    from .space import parse3
    from .verify import _rename
    names = ["x", "y", "z"][:len(vs)]
    t = _rename(text, vs, names)
    if len(vs) == 2: t = f"({t})+0*z"
    N, D = parse3(t)
    if len(D) != 1 or (0, 0, 0) not in D:
        raise DistanceError("the square root must be of a polynomial")
    c = D[(0, 0, 0)]
    return {k[:len(vs)]: Fr(v) / c for k, v in N.items() if v}


def _as_distance(arg, vs, point):
    """the set S of variables if arg = sum over S of (v - point_v)^2, else raise"""
    P = _poly(arg, vs)
    n = len(vs)
    e = lambda i, p: tuple(p if k == i else 0 for k in range(n))
    S = [i for i in range(n) if P.get(e(i, 2), 0) == 1]
    expect = {}
    for i in S:
        expect[e(i, 2)] = Fr(1)
        if point[i]: expect[e(i, 1)] = -2 * point[i]
    const = sum((point[i] ** 2 for i in S), Fr(0))
    if const: expect[(0,) * n] = const
    if len(S) < 2 or P != expect:
        where = ", ".join(f"{v} = {p}" for v, p in zip(vs, point))
        raise DistanceError(f"only a Euclidean distance to the limit point ({where}) is supported "
                            f"under a square root, such as sqrt({'+'.join(f'{v}^2' for v in vs)})")
    return S


def _prod(it):
    out = Fr(1)
    for v in it: out *= v
    return out


def _half_powers(s):
    """(P)^(k/2) -> sqrt(P)^k, so that a distance written as a power is recognized"""
    while True:
        m = re.search(r"\)\^\((-?\d+)/2\)", s)
        if not m: return s
        depth, j = 0, m.start()
        while j >= 0:
            if s[j] == ")": depth += 1
            elif s[j] == "(":
                depth -= 1
                if depth == 0: break
            j -= 1
        if j < 0: return s
        k = int(m.group(1))
        s = s[:j] + "sqrt" + s[j:m.start() + 1] + f"^({k})" + s[m.end():]


def reduce_distance(expr, variables, point=None, lead=0):
    """the rational expression g in the same variable names, to be studied at the origin, whose fluid
       state there equals the fluid state of expr at `point`; None if expr has no distance under a
       square root (other square roots, such as sqrt(1 + x), are left for fmsl.analytic).
       lead picks which variable of the distance plays u (any choice gives the same answer)."""
    if "sqrt" not in expr and "/2)" not in expr.replace(" ", ""): return None
    vs = list(variables)
    if len(vs) not in (2, 3): raise DistanceError("distances are supported in two or three variables")
    point = [Fr(p) for p in (point or [0] * len(vs))]
    s = _half_powers(expr.replace(" ", ""))
    spans = []
    for a, b, arg in _sqrt_spans(s):
        if "sqrt" in arg: continue                       # nested roots: left to fmsl.analytic
        try:
            P = _poly(arg, vs)
        except Exception:
            continue                                    # not a polynomial: left to fmsl.analytic
        at = sum((v * _prod(p ** k for p, k in zip(point, e)) for e, v in P.items()), Fr(0))
        if at != 0: continue                            # sqrt(c + ...), c != 0: analytic at the point
        spans.append((a, b, arg))
    if not spans: return None
    sets = [tuple(_as_distance(a, vs, point)) for _, _, a in spans]
    if len(set(sets)) != 1: raise DistanceError("only one distance (one square root) per formula")
    S = list(sets[0])
    S = S[lead % len(S):] + S[:lead % len(S)]
    # placeholders: the distance, then each variable
    out, last = [], 0
    for a, b, _ in spans:
        out.append(s[last:a]); out.append("__D__"); last = b
    out.append(s[last:]); t = "".join(out)
    tmp = {v: f"__V{k}__" for k, v in enumerate(vs)}
    for v, tok in tmp.items(): t = re.sub(rf"\b{re.escape(v)}\b", tok, t)
    sub = {}
    u = vs[S[0]]
    for i, v in enumerate(vs):
        if i == S[0]:
            core = f"{u}^2" + "".join(f"-{vs[j]}^2" for j in S[1:])
        elif i in S:
            core = f"2*{u}*{v}"
        else:
            core = v
        c = point[i]
        sub[tmp[v]] = f"({core})" if c == 0 else f"(({c})+{core})"
    for tok, rep in sub.items(): t = t.replace(tok, rep)
    t = t.replace("__D__", "(" + "+".join(f"{vs[j]}^2" for j in S) + ")")
    return t
