"""
fmsl.algebraic -- exact real algebraic numbers.

A real algebraic number alpha is represented by a square-free polynomial p with rational
coefficients and a rational isolating interval (lo, hi] containing exactly one real root
of p. Elements of Q(alpha) are polynomials in alpha reduced modulo p.

  * zero test:  g(alpha) = 0  <=>  h = gcd(g, p) has its root in (lo, hi]   (Sturm count)
                when g(alpha) != 0 but gcd(g, p) is nontrivial, p is replaced by p / h,
                which still has alpha as a root (so p need not be irreducible up front)
  * inverse:    extended Euclid once gcd(g, p) = 1
  * sign:       interval evaluation of g on (lo, hi], refining by bisection
All decisions are exact; floating point is used only for display.
"""
from fractions import Fraction as Fr
from numbers import Rational

_NUM = (int, Rational)

# ------------------------------ univariate polynomials over Q ------------------------------
def ptrim(p):
    p = list(p)
    while p and p[-1] == 0: p.pop()
    return p

def padd(a, b):
    n = max(len(a), len(b))
    return ptrim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])

def psub(a, b): return padd(a, [-v for v in b])

def pmul(a, b):
    if not a or not b: return []
    r = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i + j] += x * y
    return ptrim(r)

def pdivmod(a, b):
    a = ptrim(a); b = ptrim(b)
    if not b: raise ZeroDivisionError("polynomial division by zero")
    q = [Fr(0)] * max(len(a) - len(b) + 1, 1)
    while len(a) >= len(b) and a:
        k = len(a) - len(b); f = Fr(a[-1]) / b[-1]; q[k] = f
        for i, bv in enumerate(b): a[i + k] -= f * bv
        a = ptrim(a)
    return ptrim(q), a

def pgcd(a, b):
    a, b = ptrim(a), ptrim(b)
    while b: a, b = b, pdivmod(a, b)[1]
    return [Fr(v) / a[-1] for v in a] if a else a

def pderiv(p): return ptrim([i * p[i] for i in range(1, len(p))])

def peval(p, x):
    acc = Fr(0)
    for v in reversed(p): acc = acc * x + v
    return acc

def squarefree(p):
    p = ptrim(p)
    return pdivmod(p, pgcd(p, pderiv(p)))[0] if len(p) > 1 else p

def _primitive_positive(p):
    """p scaled by a POSITIVE rational to integer coefficients with gcd 1 (signs unchanged)"""
    from math import gcd
    den = 1
    for v in p:
        v = Fr(v); den = den * v.denominator // gcd(den, v.denominator)
    ints = [int(Fr(v) * den) for v in p]
    g = 0
    for v in ints: g = gcd(g, abs(v))
    return [Fr(v // g) for v in ints] if g > 1 else [Fr(v) for v in ints]

def sturm_sequence(p):
    """Sturm sequence, each member rescaled by a positive constant (keeps the coefficients
       small without changing any sign)"""
    seq = [_primitive_positive(ptrim(p)), _primitive_positive(pderiv(p))]
    while len(seq[-1]) > 1:
        r = pdivmod(seq[-2], seq[-1])[1]
        if not r: break
        seq.append(_primitive_positive([-v for v in r]))
    return [s for s in seq if s]

def _changes(vals):
    vals = [v for v in vals if v != 0]
    return sum(1 for a, b in zip(vals, vals[1:]) if (a > 0) != (b > 0))

_STURM_CACHE = {}

def _sturm_of(p):
    """(square-free part, Sturm sequence as integer coefficient lists), cached"""
    key = tuple(Fr(v) for v in ptrim(p))
    hit = _STURM_CACHE.get(key)
    if hit is None:
        sq = squarefree(list(key))
        seq = [[int(v) for v in s] for s in sturm_sequence(sq)] if len(sq) > 1 else []
        hit = (sq, seq)
        if len(_STURM_CACHE) > 2000: _STURM_CACHE.clear()
        _STURM_CACHE[key] = hit
    return hit

def _sign_at(coeffs, x):
    """sign of an integer polynomial at the rational x = a/b (b > 0), integer arithmetic only:
       p(a/b) * b^n = sum_i c_i a^i b^(n-i)"""
    x = Fr(x); a, b = x.numerator, x.denominator
    n = len(coeffs) - 1
    tot = 0
    apow, bpow = 1, b ** n
    for i in range(n + 1):
        tot += coeffs[i] * apow * bpow
        apow *= a
        if i < n: bpow //= b
    return (tot > 0) - (tot < 0)

def count_roots(p, lo, hi):
    """distinct real roots of p in (lo, hi]"""
    sq, seq = _sturm_of(p)
    if len(sq) <= 1: return 0
    def ch(x):
        sg = [t for t in (_sign_at(s, x) for s in seq) if t != 0]
        return sum(1 for u, v in zip(sg, sg[1:]) if u != v)
    return ch(lo) - ch(hi)

def root_bound(p):
    p = ptrim(p)
    return 1 + max(abs(Fr(v) / p[-1]) for v in p[:-1]) if len(p) > 1 else Fr(1)

def isolate_real_roots(p, lo=None, hi=None):
    """isolating intervals (lo, hi] for the distinct real roots of p (rational p), optionally
       only those in (lo, hi].  Returns list of ('rational', r) or ('interval', lo, hi)."""
    p = squarefree(p)
    if len(p) <= 1: return []
    B = root_bound(p)
    L = -B if lo is None else max(-B, Fr(lo))
    H = B if hi is None else min(B, Fr(hi))
    if L >= H: return []
    out, stack = [], [(L, H)]
    while stack:
        lo, hi = stack.pop()
        n = count_roots(p, lo, hi)
        if n == 0: continue
        if n == 1:
            if peval(p, hi) == 0: out.append(("rational", hi))
            else: out.append(("interval", lo, hi))
            continue
        mid = (lo + hi) / 2
        stack += [(lo, mid), (mid, hi)]
    return sorted(out, key=lambda r: r[1])

def _ieval(p, lo, hi):
    """interval enclosure of p over [lo, hi]"""
    a = b = Fr(0)
    for v in reversed(p):
        cands = (a * lo, a * hi, b * lo, b * hi)
        a, b = min(cands) + v, max(cands) + v
    return a, b

# ------------------------------ the number field Q(alpha) ------------------------------
class AlgField:
    """alpha = the unique root of the square-free polynomial p in (lo, hi]"""
    def __init__(s, p, lo, hi):
        s.p, s.lo, s.hi = squarefree(p), Fr(lo), Fr(hi)
        assert count_roots(s.p, s.lo, s.hi) == 1, "interval must isolate exactly one root"

    def refine(s):
        mid = (s.lo + s.hi) / 2
        if count_roots(s.p, s.lo, mid) == 1: s.hi = mid
        else: s.lo = mid

    def reduce(s, g):
        return pdivmod(g, s.p)[1] if len(ptrim(g)) >= len(s.p) else ptrim(g)

    def is_zero(s, g):
        r = s.reduce(g)
        if not r: return True
        h = pgcd(r, s.p)
        if len(h) <= 1: return False
        if count_roots(h, s.lo, s.hi) == 1: return True
        s.p = pdivmod(s.p, h)[0]                    # alpha is a root of the cofactor
        return False

    def sign(s, g):
        if s.is_zero(g): return 0
        r = s.reduce(g)
        while True:
            a, b = _ieval(r, s.lo, s.hi)
            if a > 0: return 1
            if b < 0: return -1
            s.refine()

    def approx(s, g, tol=Fr(1, 10 ** 40)):
        while s.hi - s.lo > tol: s.refine()
        return peval(s.reduce(g), (s.lo + s.hi) / 2)

    def inverse(s, g):
        if s.is_zero(g): raise ZeroDivisionError("division by an algebraic number equal to 0")
        a, m = s.reduce(g), s.p
        # extended Euclid: find u with u*a = 1 mod m (gcd(a, m) = 1 after is_zero)
        r0, r1, u0, u1 = m, a, [], [Fr(1)]
        while len(r1) > 1:
            q, r = pdivmod(r0, r1)
            r0, r1 = r1, r
            u0, u1 = u1, psub(u0, pmul(q, u1))
        if not r1: raise ArithmeticError("unexpected common factor")
        return s.reduce([v / r1[0] for v in u1])


class A:
    """element of Q(alpha) for a fixed real algebraic alpha"""
    __slots__ = ("F", "c")

    def __init__(s, F, c):
        s.F, s.c = F, F.reduce([Fr(v) for v in c])

    @staticmethod
    def root(F): return A(F, [0, 1])

    def _co(s, o):
        if isinstance(o, A):
            if o.F is not s.F: raise NotImplementedError("arithmetic across different algebraic fields")
            return o.c
        if isinstance(o, _NUM): return [Fr(o)] if o else []
        return None

    def _wrap(s, c):
        r = A(s.F, c)
        return r.rational() if len(r.c) <= 1 else r

    def rational(s):
        return s.c[0] if s.c else Fr(0)

    def __add__(s, o):
        c = s._co(o)
        return NotImplemented if c is None else s._wrap(padd(s.c, c))
    __radd__ = __add__
    def __sub__(s, o):
        c = s._co(o)
        return NotImplemented if c is None else s._wrap(psub(s.c, c))
    def __rsub__(s, o):
        c = s._co(o)
        return NotImplemented if c is None else s._wrap(psub(c, s.c))
    def __neg__(s): return s._wrap([-v for v in s.c])
    def __mul__(s, o):
        c = s._co(o)
        return NotImplemented if c is None else s._wrap(pmul(s.c, c))
    __rmul__ = __mul__
    def __truediv__(s, o):
        c = s._co(o)
        if c is None: return NotImplemented
        if isinstance(o, _NUM):
            if o == 0: raise ZeroDivisionError
            return s._wrap([v / Fr(o) for v in s.c])
        return s._wrap(pmul(s.c, s.F.inverse(c)))
    def __rtruediv__(s, o):
        c = s._co(o)
        return NotImplemented if c is None else s._wrap(pmul(c, s.F.inverse(s.c)))
    def __pow__(s, n):
        r = A(s.F, [1])
        b = s
        if n < 0: b, n = 1 / s, -n
        while n:
            if n & 1: r = r * b
            b = b * b; n >>= 1
        return r if isinstance(r, A) and len(r.c) > 1 else (r.rational() if isinstance(r, A) else r)

    def sign(s): return s.F.sign(s.c)
    def __bool__(s): return not s.F.is_zero(s.c)
    def __eq__(s, o):
        c = s._co(o)
        return NotImplemented if c is None else s.F.is_zero(psub(s.c, c))
    def __ne__(s, o):
        r = s.__eq__(o)
        return r if r is NotImplemented else not r
    __hash__ = None
    def __lt__(s, o): return (s - o).sign() < 0 if isinstance(s - o, A) else (s - o) < 0
    def __gt__(s, o): return (s - o).sign() > 0 if isinstance(s - o, A) else (s - o) > 0
    def __le__(s, o): return not s.__gt__(o)
    def __ge__(s, o): return not s.__lt__(o)
    def __abs__(s): return -s if s.sign() < 0 else s
    def __float__(s): return float(s.F.approx(s.c))
    def __repr__(s): return f"{float(s):.10g}"
    __str__ = __repr__


def sign_of(v):
    """sign of a Fraction or algebraic number"""
    if isinstance(v, A): return v.sign()
    return (v > 0) - (v < 0)

def same_value(u, v):
    """exact equality across Fractions and algebraic numbers. Returns True/False, or None
       if the two numbers live in different algebraic fields and could not be separated."""
    if isinstance(u, A) and isinstance(v, A) and u.F is not v.F:
        for _ in range(200):
            a1, b1 = _ieval(u.F.reduce(u.c), u.F.lo, u.F.hi)
            a2, b2 = _ieval(v.F.reduce(v.c), v.F.lo, v.F.hi)
            if b1 < a2 or b2 < a1: return False
            u.F.refine(); v.F.refine()
        return None
    return u == v
