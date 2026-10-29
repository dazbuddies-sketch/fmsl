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

Weighted distances. More generally, sqrt(Q) for a positive-definite quadratic form Q in the
distances to the point (sqrt(x^2 + 2y^2), sqrt(2x^2 + 2xy + y^2 + 3z^2), ...). If Q takes a rational
square value s^2 at some rational vector v0, a rational change of basis starting with v0 and
orthogonal for Q writes Q = s^2 (X1^2 + D2 X2^2 + D3 X3^2) with rational D2, D3 > 0, and

    X1 = u^2 - D2 v^2 - D3 w^2,   X2 = 2uv,   X3 = 2uw   give   sqrt(Q) = s (u^2 + D2 v^2 + D3 w^2),

since (u^2 - D2 v^2 - D3 w^2)^2 + D2 (2uv)^2 + D3 (2uw)^2 = (u^2 + D2 v^2 + D3 w^2)^2. The map is
onto (u^2 = (s^-1 sqrt(Q) + X1)/2, v = X2/(2u), w = X3/(2u); when u = 0 a whole ellipse maps to the
point), and q -> 0 exactly when its image does, because u^2 + D2 v^2 + D3 w^2 is comparable to |q|^2
and sqrt(Q) to the distance. So, as before, the fluid states agree. When Q takes no rational square
value (for 3x^2 + 3y^2 none exists, by Legendre's theorem), the root would bring in an irrational
constant, and FMSL refuses.

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


def _shift(P, point):
    """P(point + h) as a polynomial in h"""
    from math import comb
    out = {}
    for e, v in P.items():
        parts = [{(k,): Fr(comb(n, k)) * point[i] ** (n - k) for k in range(n + 1)} for i, n in enumerate(e)]
        terms = [((), v)]
        for part in parts:
            terms = [(ex + (k,), c * w) for ex, c in terms for (k,), w in part.items() if w]
        for ex, c in terms:
            out[ex] = out.get(ex, 0) + c
    return {e: c for e, c in out.items() if c}


def _rsqrt(c):
    from math import isqrt
    if c <= 0: return None
    p, q = c.numerator, c.denominator
    a, b = isqrt(p), isqrt(q)
    return Fr(a, b) if a * a == p and b * b == q else None


def _as_form(arg, vs, point):
    """(S, M): the variables involved and the symmetric matrix of the quadratic form Q with
       arg = Q(v - point), which must be positive definite on S; else raise"""
    P = _shift(_poly(arg, vs), point)
    n = len(vs)
    where = ", ".join(f"{v} = {p}" for v, p in zip(vs, point))
    if any(sum(e) != 2 for e in P):
        raise DistanceError(f"under a square root FMSL supports a positive-definite quadratic form in the "
                            f"distances to the limit point ({where}), such as sqrt(x^2 + 2*y^2)")
    S = sorted({i for e in P for i in range(n) if e[i]})
    M = [[Fr(0)] * len(S) for _ in S]
    for e, c in P.items():
        idx = [S.index(i) for i in range(n) for _ in range(e[i])]
        a, b = idx
        if a == b: M[a][a] += c
        else: M[a][b] += c / 2; M[b][a] += c / 2
    if len(S) < 2:
        raise DistanceError("a square root of a single square, such as sqrt(x^2) = |x|, is not supported here")
    for k in range(1, len(S) + 1):                      # Sylvester's criterion
        if _det([row[:k] for row in M[:k]]) <= 0:
            raise DistanceError(f"the quadratic form under the square root must be positive definite near ({where})")
    return S, M


def _det(A):
    A = [list(r) for r in A]; n = len(A); d = Fr(1)
    for i in range(n):
        p = next((r for r in range(i, n) if A[r][i] != 0), None)
        if p is None: return Fr(0)
        if p != i: A[i], A[p] = A[p], A[i]; d = -d
        d *= A[i][i]
        for r in range(i + 1, n):
            f = A[r][i] / A[i][i]
            A[r] = [x - f * y for x, y in zip(A[r], A[i])]
    return d


def _bil(M, a, b):
    return sum(a[i] * M[i][j] * b[j] for i in range(len(a)) for j in range(len(b)))


def _square_vectors(M, limit=12):
    """integer vectors v (primitive, up to sign) with Q(v) a nonzero rational square, small first"""
    from itertools import product
    from math import gcd
    n = len(M); seen = []
    for size in range(1, limit + 1):
        vecs = sorted(product(range(-size, size + 1), repeat=n), key=lambda v: (sum(1 for x in v if x), [-x for x in v]))
        for v in vecs:
            if max(abs(x) for x in v) != size: continue
            g = 0
            for x in v: g = gcd(g, abs(x))
            if g != 1: continue
            if next(x for x in v if x) < 0: continue
            if _rsqrt(_bil(M, v, v)) is not None:
                seen.append([Fr(x) for x in v])
                if len(seen) >= 6: return seen
    return seen


def _basis(M, v0):
    """b_1 = v0, then a Q-orthogonal completion (Gram-Schmidt in the bilinear form of Q)"""
    n = len(M); B = [v0]
    for j in range(n):
        e = [Fr(1 if k == j else 0) for k in range(n)]
        for b in B:
            e = [x - _bil(M, e, b) / _bil(M, b, b) * y for x, y in zip(e, b)]
        if any(e) and _bil(M, e, e) != 0: B.append(e)
        if len(B) == n: break
    return B


def _prod(it):
    out = Fr(1)
    for v in it: out *= v
    return out


def has_half_power(expr):
    """a power (P)^(k/2), the other way to write a square root"""
    return re.search(r"\)\^\(-?\d+/2\)", expr.replace(" ", "")) is not None


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
    if "sqrt" not in expr and not has_half_power(expr): return None
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
    forms = [_as_form(a, vs, point) for _, _, a in spans]
    if any(f != forms[0] for f in forms): raise DistanceError("only one square root (one distance) per formula")
    S, M = forms[0]
    cands = _square_vectors(M)
    if not cands:
        raise DistanceError("this quadratic form takes no rational square value at small integer points, "
                            "so the square root cannot be removed with rational coefficients")
    v0 = cands[lead % len(cands)]
    B = _basis(M, v0)
    d = [_bil(M, b, b) for b in B]
    s0 = _rsqrt(d[0])
    D = [di / d[0] for di in d]
    names = [vs[i] for i in S]                          # new coordinates (u, v_2, ...) reuse these names
    u = names[0]
    X = [f"({u}^2" + "".join(f"-({D[k]})*{names[k]}^2" for k in range(1, len(S))) + ")"] + \
        [f"(2*{u}*{names[k]})" for k in range(1, len(S))]
    # placeholders: the root, then each variable
    out, last = [], 0
    for a, b, _ in spans:
        out.append(s[last:a]); out.append("__D__"); last = b
    out.append(s[last:]); t = "".join(out)
    tmp = {v: f"__V{k}__" for k, v in enumerate(vs)}
    for v, tok in tmp.items(): t = re.sub(rf"\b{re.escape(v)}\b", tok, t)
    for i, v in enumerate(vs):
        if i in S:
            j = S.index(i)
            core = "+".join(f"({B[k][j]})*{X[k]}" for k in range(len(S)) if B[k][j] != 0) or "0"
        else:
            core = v
        c = point[i]
        t = t.replace(tmp[v], f"({core})" if c == 0 else f"(({c})+{core})")
    t = t.replace("__D__", f"(({s0})*(" + "+".join([f"{u}^2"] + [f"({D[k]})*{names[k]}^2" for k in range(1, len(S))]) + "))")
    return t
