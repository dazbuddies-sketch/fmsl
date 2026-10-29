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


# ---------------------------------------------------------------- expression -> (A, B)
class _SeriesRing:
    """truncated power series at the origin"""
    def __init__(s, names, K): s.names, s.n, s.K, s.N = list(names), len(names), K, K
    def const(s, v): return S.const(v, s.n, s.K)
    def var(s, name): return S.var(s.names.index(name), s.n, s.K)
    def is_one(s, P): return P.c == {(0,) * s.n: 1}


class Undecided(Exception):
    """an exact restriction needs a value FMSL cannot represent exactly (e.g. sin of a non-constant)"""


class R:
    """an element of Q[t][w]/(w^J): the restriction of a germ to a line (t along it, w across it),
       kept to order J across; t-coefficients are exact polynomials"""
    __slots__ = ("c", "J")

    def __init__(s, c, J): s.c = {e: v for e, v in c.items() if v != 0 and e[1] < J}; s.J = J

    def c0(s):
        w0 = {a: v for (a, b), v in s.c.items() if b == 0}
        if any(a > 0 for a in w0): raise Undecided()
        return w0.get(0, Fr(0))

    def __add__(s, o):
        if isinstance(o, O): return O(min(s.w_order(), o.lb))
        c = dict(s.c)
        for e, v in o.c.items(): c[e] = c.get(e, 0) + v
        return R(c, s.J)

    def __neg__(s): return R({e: -v for e, v in s.c.items()}, s.J)
    def __sub__(s, o): return s + (-o)
    def scale(s, k): return R({e: v * k for e, v in s.c.items()}, s.J)

    def __mul__(s, o):
        if isinstance(o, O): return O(s.w_order() + o.lb)
        c = {}
        for (a1, b1), v1 in s.c.items():
            for (a2, b2), v2 in o.c.items():
                if b1 + b2 >= s.J: continue
                e = (a1 + a2, b1 + b2); c[e] = c.get(e, 0) + v1 * v2
        return R(c, s.J)

    def compose(s, coeffs):
        out, p = R({(0, 0): Fr(coeffs[0])}, s.J), R({(0, 0): Fr(1)}, s.J)
        for k in range(1, len(coeffs)):
            p = p * s
            if not p.c: break
            if coeffs[k]: out = out + p.scale(coeffs[k])
        return out

    def inverse(s):
        c = s.c0()
        if c == 0: raise Undecided()
        u = (s - R({(0, 0): c}, s.J)).scale(Fr(1) / c)
        return u.compose([Fr((-1) ** k) for k in range(s.J + 1)]).scale(Fr(1) / c)

    def w_order(s):
        return min((b for (a, b) in s.c), default=s.J)


class O:
    """an analytic germ restricted to the line whose expansion across it is not known exactly,
       only a lower bound lb for its order in w"""
    __slots__ = ("lb",)
    def __init__(s, lb): s.lb = lb
    def __add__(s, o): return O(min(s.lb, _lb(o)))
    __radd__ = __add__
    def __neg__(s): return s
    def __sub__(s, o): return O(min(s.lb, _lb(o)))
    def __rsub__(s, o): return O(min(s.lb, _lb(o)))
    def scale(s, k): return s if k != 0 else R({}, 1)
    def __mul__(s, o): return O(s.lb + _lb(o))
    __rmul__ = __mul__
    def w_order(s): return s.lb


def _lb(x):
    return x.lb if isinstance(x, O) else x.w_order()


def _opaque_call(name, s):
    """F(s) when s is not exactly known: F(0) = 0 for these functions, so F(s) = s G(s)"""
    if name in ("sin", "sinh", "atan", "asin", "tan", "tanh"): return O(_lb(s))
    return O(0)


class _LineRing:
    """germs restricted to a line through the origin: each variable becomes an element of R"""
    def __init__(s, subs, J): s.subs, s.J, s.N = subs, J, J + 1
    def const(s, v): return R({(0, 0): Fr(v)}, s.J)
    def var(s, name): return s.subs[name]
    def is_one(s, P): return isinstance(P, R) and P.c == {(0, 0): 1}


def _eval(e, ring):
    """(N, D) with f = N/D as analytic germs at the origin, in the given ring"""
    one = ring.const(1)
    N_ = ring.N
    t = e[0]
    if t == "num":
        v = e[1]
        if not isinstance(v, Fr): raise NotAnalytic("complex numbers are not supported here")
        return ring.const(v), one
    if t == "var":
        if e[1] not in getattr(ring, "names", e[1]) and not isinstance(ring, _LineRing):
            raise NotAnalytic(f"unknown name {e[1]!r}")
        return ring.var(e[1]), one
    if t == "neg":
        N, D = _eval(e[1], ring); return -N, D
    if t in ("+", "-"):
        (N1, D1), (N2, D2) = _eval(e[1], ring), _eval(e[2], ring)
        if ring.is_one(D1) and ring.is_one(D2): return (N1 + N2 if t == "+" else N1 - N2), one
        return (N1 * D2 + N2 * D1 if t == "+" else N1 * D2 - N2 * D1), D1 * D2
    if t == "*":
        (N1, D1), (N2, D2) = _eval(e[1], ring), _eval(e[2], ring)
        return N1 * N2, D1 * D2
    if t == "/":
        (N1, D1), (N2, D2) = _eval(e[1], ring), _eval(e[2], ring)
        return N1 * D2, D1 * N2
    if t == "pow":
        N, D = _eval(e[1], ring)
        p = e[2]
        if isinstance(p, int):
            if p < 0: N, D, p = D, N, -p
            RN, RD = one, one
            for _ in range(p): RN, RD = RN * N, RD * D
            return RN, RD
        s = _germ(N, D, "a fractional power")             # s^(p/q) with s(0) = c > 0
        if isinstance(ring, _LineRing):
            try:
                if isinstance(s, O): raise Undecided()
                c = s.c0()
            except Undecided:
                return O(0), one
        else:
            c = s.c0()
        root = _rational_power(c, p)
        if root is None: raise NotAnalytic("a fractional power needs a base whose value at the point has a rational power")
        u = (s - ring.const(c)).scale(Fr(1) / c)
        return u.compose([_binom(p, k) for k in range(N_ + 1)]).scale(root), one
    if t == "call":
        name = e[1]
        N, D = _eval(e[2], ring)
        s = _germ(N, D, name)
        if isinstance(ring, _LineRing):
            try:
                if isinstance(s, O): raise Undecided()
                c = s.c0()
            except Undecided:
                if name in ("tan", "tanh"): return _opaque_call(name, s), O(0)
                return _opaque_call(name, s), one
        else:
            c = s.c0()
        if name in ("tan", "tanh"):
            if c != 0: _bad(name)
            a, b = ("sin", "cos") if name == "tan" else ("sinh", "cosh")
            return s.compose(_series(a, N_)), s.compose(_series(b, N_))
        if name in ("sin", "cos", "exp", "sinh", "cosh", "atan", "asin"):
            if c != 0: _bad(name)
            return s.compose(_series(name, N_)), one
        if name == "log":
            if c != 1: raise NotAnalytic("log needs an argument equal to 1 at the point, as in log(1 + u)")
            return (s - one).compose(_series("log1p", N_)), one
        if name == "sqrt":
            root = _rational_sqrt(c) if c > 0 else None
            if root is None:
                raise NotAnalytic("sqrt needs an argument whose value at the point is a positive rational square; "
                                  "a distance sqrt(x^2+y^2) is handled separately")
            u = (s - ring.const(c)).scale(Fr(1) / c)
            return u.compose([_binom(Fr(1, 2), k) for k in range(N_ + 1)]).scale(root), one
        raise NotAnalytic(f"{name} is not supported in several variables")
    raise NotAnalytic(f"unsupported expression {t!r}")


def _bad(name):
    raise NotAnalytic(f"{name} needs an argument equal to 0 at the point (e.g. {name}(x*y))")


def _germ(N, D, what):
    if isinstance(D, O) or isinstance(N, O):
        return O(_lb(N))                                # 1/D is analytic, so the order is at least N's
    if isinstance(D, R):
        try:
            if D.c0() == 0: raise Undecided()
            return N * D.inverse()
        except Undecided:
            return O(_lb(N))
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
    return _eval(Parser(expr).parse(), _SeriesRing(names, K))


def is_rational(expr):
    import re
    return not re.search(r"[a-z_]\w*\s*\(", expr)


# ---------------------------------------------------------------- the proof
def _trunc(P, K):
    return S(P.c, P.n, K)


# ---------------------------------------------------------------- common factors along lines
def _lines(B):
    """rational lines through the origin inside the zero set of B's lowest-degree part:
       (alpha, beta, subs) with l = alpha x + beta y, and the substitution putting t along the line
       and w = l across it"""
    from .cluster import _real_roots
    b = B.order()
    h = {e: v for e, v in B.c.items() if sum(e) == b}
    out = []
    poly = [h.get((b - j, j), Fr(0)) for j in range(b + 1)]          # h(1, s)
    while len(poly) > 1 and poly[-1] == 0: poly.pop()
    if len(poly) > 1:
        for s0 in _real_roots(poly):
            if isinstance(s0, Fr):
                out.append((-s0, Fr(1), {"x": R({(1, 0): Fr(1)}, 99), "y": R({(1, 0): s0, (0, 1): Fr(1)}, 99)}))
    if h.get((0, b), Fr(0)) == 0:                                   # the line x = 0
        out.append((Fr(1), Fr(0), {"x": R({(0, 1): Fr(1)}, 99), "y": R({(1, 0): Fr(1)}, 99)}))
    return out


def _multiplicities(e, subs, J):
    """(order of N, order of D) across the line, each capped at J; None if not decidable exactly"""
    subs = {k: R(v.c, J) for k, v in subs.items()}
    try:
        N, D = _eval(e, _LineRing(subs, J))
    except (Undecided, NotAnalytic, ZeroDivisionError):
        return None
    return min(_lb(N), J), min(_lb(D), J)


def _div_line(P, alpha, beta, m, K):
    """P / (alpha x + beta y)^m for a two-variable series known to be divisible, truncated at K"""
    for _ in range(m):
        out = {}
        for d in range(0, P.K + 1):
            hs = [P.c.get((d - j, j), Fr(0)) for j in range(d + 1)]
            if not any(hs): continue
            if beta == 0:                                            # divide by alpha x
                if hs[d] != 0: raise ArithmeticError("not divisible")
                for j in range(d): out[(d - 1 - j, j)] = hs[j] / alpha
                continue
            q = [Fr(0)] * d                                          # h(s) / (beta s + alpha)
            rem = list(hs)
            for j in range(d, 0, -1):
                q[j - 1] = rem[j] / beta
                rem[j] -= q[j - 1] * beta; rem[j - 1] -= q[j - 1] * alpha
            if rem[0] != 0: raise ArithmeticError("not divisible")
            for j in range(d): out[(d - 1 - j, j)] = q[j]
        P = S(out, P.n, P.K - 1)
    return S(P.c, P.n, K)


def _quotient_series(expr, names):
    """series(K) -> (A, B) for f = A/B after cancelling every rational line through the origin on
       which both the numerator and the denominator vanish (two variables), and a note"""
    from .compiler import Parser
    e = Parser(expr).parse()
    plain = lambda K: _eval(e, _SeriesRing(names, K))
    if len(names) != 2: return plain, ""
    K = 2
    A, B = plain(K)
    while not B.c and K < 24: K += 2; A, B = plain(K)
    if not B.c or B.c0() != 0: return plain, ""
    cancel, notes = [], []
    for alpha, beta, subs in _lines(B):
        subs = {names[0]: subs["x"], names[1]: subs["y"]}
        mk = _multiplicities(e, subs, B.order() + 1)
        if mk is None: continue
        j, k = mk
        m = min(j, k)
        if m > 0:
            cancel.append((alpha, beta, m))
            lin = f"{alpha}*{names[0]} + {beta}*{names[1]}".replace("1*", "")
            notes.append(f"({lin})^{m}")
    if not cancel: return plain, ""
    M = sum(m for _, _, m in cancel)
    def series(K):
        A, B = plain(K + M)
        for alpha, beta, m in cancel:
            A = _div_line(A, alpha, beta, m, A.K - m) if A.c else S({}, 2, A.K - m)
            B = _div_line(B, alpha, beta, m, B.K - m)
        return S(A.c, 2, K), S(B.c, 2, K)
    return series, ("both numerator and denominator vanish on " + ", ".join(notes) +
                    " = 0 (proven by exact restriction to each line); cancelled exactly")


def analytic_limit(expr, names, engine, k_max=24, extra_L=6):
    """(components, proven, explanation) for expr in the variables `names` (["x","y"] or
       ["x","y","z"]); engine(text) -> (components, proven, explanation, witnesses) is FMSL's
       rational engine for that number of variables.

       With a = order of A, b = order of B, and |B_Kb| >= c r^L, the bound
           |A/B - A_Ka/B_Kb| <= C (r^(Ka+1+b) + r^(Kb+1+a)) / r^(2L)
       tends to 0 when Ka >= 2L - b, Kb >= 2L - a and Kb >= L, which are the degrees used.

       In two variables, a rational line l = 0 through the origin on which both A and B vanish is
       first cancelled: A and B are divided by l^m exactly, where m is the smaller of their orders
       across the line. Those orders are bounded from below exactly, by restricting the expression
       to the line (t along it, w = l across it) in Q[t][w]/(w^J); restriction commutes with every
       operation used. A part that cannot be restricted exactly (sin(t) along the line) is kept
       only through a lower bound on its order (F(s) = s G(s) when F(0) = 0; products add orders,
       sums take the minimum), so the cancelled power never exceeds the true one."""
    import math
    series, cnote = _quotient_series(expr, names)
    pre = (cnote + "; ") if cnote else ""
    r2 = "(" + "+".join(f"{x}^2" for x in names) + ")"
    K = 2
    A, B = series(K)
    while not B.c and K < k_max:                       # find the order of the denominator
        K += 2; A, B = series(K)
    if not B.c:
        raise NotAnalytic(f"the denominator vanishes to order > {k_max} at the point (or identically)")
    if B.c0() != 0:                                     # f is analytic at the point: its value is the limit
        v = A.c0() / B.c0()
        return [(v, v)], True, pre + f"analytic at the point; value {v}"
    b = B.order()
    tried = []
    for L in range(b, b + extra_L + 1):
        a = A.order() if A.c else K + 1                 # a lower bound for the order of A is enough
        Ka, Kb = max(0, 2 * L - b), max(L, 2 * L - a)
        need = max(Ka, Kb)
        if need > k_max: break
        if need > K:
            K = need; A, B = series(K)
            a = A.order() if A.c else K + 1
            Ka, Kb = max(0, 2 * L - b), max(L, 2 * L - a)
        Bk = _trunc(B, Kb)
        comps, proven, _, _ = engine(f"{r2}^{L}/({_text(Bk, names)})^2")
        tried.append(L)
        if comps is None or not proven or any(math.isinf(u) or math.isinf(v) for u, v in comps):
            continue
        Ak = _trunc(A, Ka)
        head = (pre + f"Taylor polynomials A of degree {Ka} and B of degree {Kb}; |B| >= c r^{L} proven "
                f"(r^{2 * L}/B^2 bounded), so the remainders change no limit value")
        if not Ak.c:
            return [(0.0, 0.0)], True, head + "; the numerator's polynomial is 0, so the limit is 0"
        comps, proven, expl, _ = engine(f"({_text(Ak, names)})/({_text(Bk, names)})")
        return comps, proven, head + "; " + expl
    return None, False, (pre + "not proven: the denominator's zero at the point is not isolated, or no bound "
                         f"|denominator| >= c r^L was found for L in {tried} (degree <= {k_max})")
