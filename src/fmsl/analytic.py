"""sin, cos, exp, log, sqrt and friends in two and three variables, reduced to rational functions
with a proof that nothing is lost.

Write f = A/B, where A and B are analytic at the origin (built from polynomials with +, -, *, /,
integer powers, and sin, cos, tan, exp, sinh, cosh, tanh, atan, asin of arguments that vanish at
the origin, log(1 + u) and sqrt(c^2 + u) with u vanishing there). Let A_K, B_K be their Taylor
polynomials of total degree K, computed exactly (all coefficients are rational).

Theorem (the Lojasiewicz route, as in Strzebonski, arXiv:2102.01242). Let a, b be the orders of A
and B at the origin (b > 0). Suppose that, for an integer L, the rational function
(x^2 + y^2 [+ z^2])^L / B_Kb^2 is bounded near the origin, that is |B_Kb| >= c r^L, where
Kb >= max(L, 2L - a), and let Ka >= 2L - b. Then A/B and A_Ka/B_Kb have the same set of limit values
at the origin.

Proof. A - A_Ka = O(r^(Ka+1)) and B - B_Kb = O(r^(Kb+1)) because A and B are analytic, while
A_Ka = O(r^a) and B_Kb = O(r^b). Since Kb + 1 > L, |B| >= |B_Kb| - C r^(Kb+1) >= (c/2) r^L near 0,
and
    |A/B - A_Ka/B_Kb| = |(A - A_Ka) B_Kb - A_Ka (B - B_Kb)| / |B B_Kb|
                      <= C' (r^(Ka+1+b) + r^(a+Kb+1)) / r^(2L) -> 0.
Along every sequence tending to the origin the two quotients differ by something tending to 0,
so they have the same limit values (including +-oo).

Both conditions are checked by FMSL's rational engines: the boundedness of r^(2L)/B_K^2 is itself
a proven fluid state, and the fluid state of A_K/B_K is the answer. When no such L exists (the
denominator vanishes along curves through the origin), the method does not apply and FMSL says
"not proven" rather than guess.
"""
from fractions import Fraction as Fr
from math import factorial

__all__ = ["analytic_limit", "taylor_quotient", "NotAnalytic"]


class NotAnalytic(ValueError):
    pass


# ---------------------------------------------------------------- truncated power series
class S:
    """a power series in n variables, truncated at total degree K: {exponent tuple: Fraction}"""
    __slots__ = ("c", "n", "K")

    def __init__(s, c, n, K):
        s.c = {e: v for e, v in c.items() if v != 0 and sum(e) <= K}; s.n = n; s.K = K

    @staticmethod
    def const(v, n, K): return S({(0,) * n: Fr(v)}, n, K)

    @staticmethod
    def var(i, n, K): return S({tuple(1 if j == i else 0 for j in range(n)): Fr(1)}, n, K)

    def c0(s): return s.c.get((0,) * s.n, Fr(0))

    def order(s): return min((sum(e) for e in s.c), default=None)

    def __add__(s, o):
        c = dict(s.c)
        for e, v in o.c.items(): c[e] = c.get(e, 0) + v
        return S(c, s.n, s.K)

    def __neg__(s): return S({e: -v for e, v in s.c.items()}, s.n, s.K)

    def __sub__(s, o): return s + (-o)

    def scale(s, k): return S({e: v * k for e, v in s.c.items()}, s.n, s.K)

    def __mul__(s, o):
        c = {}
        for e1, v1 in s.c.items():
            d1 = sum(e1)
            for e2, v2 in o.c.items():
                if d1 + sum(e2) > s.K: continue
                e = tuple(a + b for a, b in zip(e1, e2))
                c[e] = c.get(e, 0) + v1 * v2
        return S(c, s.n, s.K)

    def compose(s, coeffs):
        """sum_k coeffs[k] * s^k for s with zero constant term (coeffs: list of Fractions)"""
        out, p = S.const(coeffs[0], s.n, s.K), S.const(1, s.n, s.K)
        for k in range(1, len(coeffs)):
            p = p * s
            if not p.c: break
            if coeffs[k]: out = out + p.scale(coeffs[k])
        return out

    def inverse(s):
        c = s.c0()
        if c == 0: raise ZeroDivisionError
        u = (s - S.const(c, s.n, s.K)).scale(Fr(1) / c)          # s = c (1 + u)
        return u.compose([Fr((-1) ** k) for k in range(s.K + 1)]).scale(Fr(1) / c)


def _binom(a, k):
    out = Fr(1)
    for j in range(k): out = out * (a - j) / (j + 1)
    return out


def _series(name, K):
    """Maclaurin coefficients up to degree K (all rational)"""
    r = range(K + 1)
    if name == "sin": return [Fr((-1) ** ((k - 1) // 2), factorial(k)) if k % 2 else Fr(0) for k in r]
    if name == "cos": return [Fr((-1) ** (k // 2), factorial(k)) if k % 2 == 0 else Fr(0) for k in r]
    if name == "exp": return [Fr(1, factorial(k)) for k in r]
    if name == "sinh": return [Fr(1, factorial(k)) if k % 2 else Fr(0) for k in r]
    if name == "cosh": return [Fr(1, factorial(k)) if k % 2 == 0 else Fr(0) for k in r]
    if name == "log1p": return [Fr(0)] + [Fr((-1) ** (k + 1), k) for k in r if k > 0]
    if name == "atan": return [Fr((-1) ** ((k - 1) // 2), k) if k % 2 else Fr(0) for k in r]
    if name == "asin":                                   # sum binom(2m, m) x^(2m+1) / (4^m (2m+1))
        return [Fr(factorial(2 * (k // 2)), 4 ** (k // 2) * factorial(k // 2) ** 2 * k) if k % 2 else Fr(0) for k in r]
    raise KeyError(name)


def _rational_sqrt(c):
    """sqrt of a positive Fraction if rational, else None"""
    from math import isqrt
    p, q = c.numerator, c.denominator
    if p <= 0: return None
    a, b = isqrt(p), isqrt(q)
    return Fr(a, b) if a * a == p and b * b == q else None


# ---------------------------------------------------------------- expression -> (A, B) series
def _eval(e, names, K):
    """(N, D) truncated series with f = N/D as analytic germs at the origin"""
    n = len(names)
    one = S.const(1, n, K)
    t = e[0]
    if t == "num":
        v = e[1]
        if not isinstance(v, Fr): raise NotAnalytic("complex numbers are not supported here")
        return S.const(v, n, K), one
    if t == "var":
        if e[1] not in names: raise NotAnalytic(f"unknown name {e[1]!r}")
        return S.var(names.index(e[1]), n, K), one
    if t == "neg":
        N, D = _eval(e[1], names, K); return -N, D
    if t in ("+", "-"):
        (N1, D1), (N2, D2) = _eval(e[1], names, K), _eval(e[2], names, K)
        if D1.c == one.c and D2.c == one.c: return (N1 + N2 if t == "+" else N1 - N2), one
        return (N1 * D2 + N2 * D1 if t == "+" else N1 * D2 - N2 * D1), D1 * D2
    if t == "*":
        (N1, D1), (N2, D2) = _eval(e[1], names, K), _eval(e[2], names, K)
        return N1 * N2, D1 * D2
    if t == "/":
        (N1, D1), (N2, D2) = _eval(e[1], names, K), _eval(e[2], names, K)
        return N1 * D2, D1 * N2
    if t == "pow":
        N, D = _eval(e[1], names, K)
        p = e[2]
        if isinstance(p, int):
            if p < 0: N, D, p = D, N, -p
            RN, RD = one, one
            for _ in range(p): RN, RD = RN * N, RD * D
            return RN, RD
        s = _germ(N, D, "a fractional power")             # s^(p/q) with s(0) = c > 0
        c = s.c0()
        root = _rational_power(c, p)
        if root is None: raise NotAnalytic("a fractional power needs a base whose value at the point has a rational power")
        u = (s - S.const(c, n, K)).scale(Fr(1) / c)
        return u.compose([_binom(p, k) for k in range(K + 1)]).scale(root), one
    if t == "call":
        name = e[1]
        N, D = _eval(e[2], names, K)
        s = _germ(N, D, name)
        c = s.c0()
        if name in ("tan", "tanh"):
            if c != 0: _bad(name)
            a, b = ("sin", "cos") if name == "tan" else ("sinh", "cosh")
            return s.compose(_series(a, K)), s.compose(_series(b, K))
        if name in ("sin", "cos", "exp", "sinh", "cosh", "atan", "asin"):
            if c != 0: _bad(name)
            return s.compose(_series(name, K)), one
        if name == "log":
            if c != 1: raise NotAnalytic("log needs an argument equal to 1 at the point, as in log(1 + u)")
            return (s - one).compose(_series("log1p", K)), one
        if name == "sqrt":
            root = _rational_sqrt(c) if c > 0 else None
            if root is None:
                raise NotAnalytic("sqrt needs an argument whose value at the point is a positive rational square; "
                                  "a distance sqrt(x^2+y^2) is handled separately")
            u = (s - S.const(c, n, K)).scale(Fr(1) / c)
            return u.compose([_binom(Fr(1, 2), k) for k in range(K + 1)]).scale(root), one
        raise NotAnalytic(f"{name} is not supported in several variables")
    raise NotAnalytic(f"unsupported expression {t!r}")


def _bad(name):
    raise NotAnalytic(f"{name} needs an argument equal to 0 at the point (e.g. {name}(x*y))")


def _germ(N, D, what):
    if D.c0() == 0:
        raise NotAnalytic(f"the argument of {what} must be defined at the point (no division by something that is 0 there)")
    return N * D.inverse()


def _rational_power(c, p):
    from math import isqrt
    if c <= 0: return None
    a, b = p.numerator, p.denominator
    def root(m):
        r = round(m ** (1.0 / b))
        for k in (r - 1, r, r + 1):
            if k >= 0 and k ** b == m: return k
        return None
    num, den = root(c.numerator), root(c.denominator)
    if num is None or den is None: return None
    return Fr(num, den) ** a


def _text(P, names):
    terms = []
    for e, v in sorted(P.c.items()):
        mono = "*".join(f"{x}^{k}" for x, k in zip(names, e) if k)
        terms.append(f"({v})" + ("*" + mono if mono else ""))
    return " + ".join(terms) if terms else "0"


def taylor_quotient(expr, names, K):
    """(A_K, B_K) as series for expr = A/B (names: the variable names, e.g. ["x", "y"])"""
    from .compiler import Parser
    return _eval(Parser(expr).parse(), list(names), K)


def is_rational(expr):
    import re
    return not re.search(r"[a-z_]\w*\s*\(", expr)


# ---------------------------------------------------------------- the proof
def _trunc(P, K):
    return S(P.c, P.n, K)


def analytic_limit(expr, names, engine, k_max=24, extra_L=6):
    """(components, proven, explanation) for expr in the variables `names` (["x","y"] or
       ["x","y","z"]); engine(text) -> (components, proven, explanation, witnesses) is FMSL's
       rational engine for that number of variables.

       With a = order of A, b = order of B, and |B_Kb| >= c r^L, the bound
           |A/B - A_Ka/B_Kb| <= C (r^(Ka+1+b) + r^(Kb+1+a)) / r^(2L)
       tends to 0 when Ka >= 2L - b, Kb >= 2L - a and Kb >= L, which are the degrees used."""
    import math
    r2 = "(" + "+".join(f"{x}^2" for x in names) + ")"
    K = 2
    A, B = taylor_quotient(expr, names, K)
    while not B.c and K < k_max:                       # find the order of the denominator
        K += 2; A, B = taylor_quotient(expr, names, K)
    if not B.c:
        raise NotAnalytic(f"the denominator vanishes to order > {k_max} at the point (or identically)")
    if B.c0() != 0:                                     # f is analytic at the point: its value is the limit
        v = A.c0() / B.c0()
        return [(v, v)], True, f"analytic at the point; value {v}"
    b = B.order()
    tried = []
    for L in range(b, b + extra_L + 1):
        a = A.order() if A.c else K + 1                 # a lower bound for the order of A is enough
        Ka, Kb = max(0, 2 * L - b), max(L, 2 * L - a)
        need = max(Ka, Kb)
        if need > k_max: break
        if need > K:
            K = need; A, B = taylor_quotient(expr, names, K)
            a = A.order() if A.c else K + 1
            Ka, Kb = max(0, 2 * L - b), max(L, 2 * L - a)
        Bk = _trunc(B, Kb)
        comps, proven, _, _ = engine(f"{r2}^{L}/({_text(Bk, names)})^2")
        tried.append(L)
        if comps is None or not proven or any(math.isinf(u) or math.isinf(v) for u, v in comps):
            continue
        Ak = _trunc(A, Ka)
        head = (f"Taylor polynomials A of degree {Ka} and B of degree {Kb}; |B| >= c r^{L} proven "
                f"(r^{2 * L}/B^2 bounded), so the remainders change no limit value")
        if not Ak.c:
            return [(0.0, 0.0)], True, head + "; the numerator's polynomial is 0, so the limit is 0"
        comps, proven, expl, _ = engine(f"({_text(Ak, names)})/({_text(Bk, names)})")
        return comps, proven, head + "; " + expl
    return None, False, ("not proven: the denominator's zero at the point is not isolated, or no bound "
                         f"|denominator| >= c r^L was found for L in {tried} (degree <= {k_max})")
