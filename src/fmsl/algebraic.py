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

def _div(u, v):
    """exact quotient; rationals stay Fractions, algebraic numbers stay algebraic"""
    if isinstance(u, A) or isinstance(v, A): return u / v
    return Fr(u) / v

def is_alg_poly(p):
    """True if some coefficient is a (non-rational) algebraic number"""
    return any(isinstance(v, A) for v in p)

def pdivmod(a, b):
    a = ptrim(a); b = ptrim(b)
    if not b: raise ZeroDivisionError("polynomial division by zero")
    q = [Fr(0)] * max(len(a) - len(b) + 1, 1)
    while len(a) >= len(b) and a:
        k = len(a) - len(b); f = _div(a[-1], b[-1]); q[k] = f
        for i, bv in enumerate(b): a[i + k] -= f * bv
        a = ptrim(a)
    return ptrim(q), a

def pgcd(a, b):
    a, b = ptrim(a), ptrim(b)
    while b: a, b = b, pdivmod(a, b)[1]
    return [_div(v, a[-1]) for v in a] if a else a

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

def coef(v):
    """a coefficient as an exact number: algebraic numbers stay algebraic, the rest become Fractions"""
    return v if isinstance(v, A) else Fr(v)

def count_roots(p, lo, hi):
    """distinct real roots of p in (lo, hi]"""
    if is_alg_poly(p):
        return _count_generic(_sturm_generic(_lift_poly(ptrim(p), _coeff_field(p)))[1], lo, hi)
    sq, seq = _sturm_of(p)
    if len(sq) <= 1: return 0
    def ch(x):
        sg = [t for t in (_sign_at(s, x) for s in seq) if t != 0]
        return sum(1 for u, v in zip(sg, sg[1:]) if u != v)
    return ch(lo) - ch(hi)

def root_bound(p):
    p = ptrim(p)
    if is_alg_poly(p):
        return 1 + max(_abs_upper(v) for v in p[:-1]) / _abs_lower(p[-1]) if len(p) > 1 else Fr(1)
    return 1 + max(abs(Fr(v) / p[-1]) for v in p[:-1]) if len(p) > 1 else Fr(1)

def isolate_real_roots(p, lo=None, hi=None):
    """isolating intervals (lo, hi] for the distinct real roots of p (rational p), optionally
       only those in (lo, hi].  Returns list of ('rational', r) or ('interval', lo, hi)."""
    if is_alg_poly(p): return isolate_real_roots_alg(p, lo, hi)[0]
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
        # alpha is the only root of the square-free p in (lo, hi], and it is simple, so p
        # changes sign exactly there: one sign evaluation halves the interval
        mid = (s.lo + s.hi) / 2
        ip = _sturm_of(s.p)[1]
        if not ip:
            if count_roots(s.p, s.lo, mid) == 1: s.hi = mid
            else: s.lo = mid
            return
        P = ip[0]                                  # p scaled by a positive constant, integers
        sm = _sign_at(P, mid)
        if sm == 0:                                # alpha = mid exactly (rational)
            s.lo, s.hi = (s.lo + mid) / 2, mid
            return
        sh = _sign_at(P, s.hi)
        if sh == 0 or sh != sm: s.lo = mid         # root in (mid, hi]
        else: s.hi = mid                           # root in (lo, mid]

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
            if o.F is not s.F: return None          # handled by _cross
            return o.c
        if isinstance(o, _NUM): return [Fr(o)] if o else []
        return None

    def _cross(s, o):
        """(s, o) as elements of one common field, or None if o is not a number"""
        if isinstance(o, A) and o.F is not s.F:
            L = unify_fields(s.F, o.F)
            return embed(s, L), embed(o, L)
        return None

    def _wrap(s, c):
        r = A(s.F, c)
        return r.rational() if len(r.c) <= 1 else r

    def rational(s):
        return s.c[0] if s.c else Fr(0)

    def __add__(s, o):
        x = s._cross(o)
        if x: return x[0] + x[1]
        c = s._co(o)
        return NotImplemented if c is None else s._wrap(padd(s.c, c))
    __radd__ = __add__
    def __sub__(s, o):
        x = s._cross(o)
        if x: return x[0] - x[1]
        c = s._co(o)
        return NotImplemented if c is None else s._wrap(psub(s.c, c))
    def __rsub__(s, o):
        x = s._cross(o)
        if x: return x[1] - x[0]
        c = s._co(o)
        return NotImplemented if c is None else s._wrap(psub(c, s.c))
    def __neg__(s): return s._wrap([-v for v in s.c])
    def __mul__(s, o):
        x = s._cross(o)
        if x: return x[0] * x[1]
        c = s._co(o)
        return NotImplemented if c is None else s._wrap(pmul(s.c, c))
    __rmul__ = __mul__
    def __truediv__(s, o):
        x = s._cross(o)
        if x: return x[0] / x[1]
        c = s._co(o)
        if c is None: return NotImplemented
        if isinstance(o, _NUM):
            if o == 0: raise ZeroDivisionError
            return s._wrap([v / Fr(o) for v in s.c])
        return s._wrap(pmul(s.c, s.F.inverse(c)))
    def __rtruediv__(s, o):
        x = s._cross(o)
        if x: return x[1] / x[0]
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
        x = s._cross(o)
        if x: return x[0] == x[1]
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
        for _ in range(10):
            a1, b1 = _ieval(u.F.reduce(u.c), u.F.lo, u.F.hi)
            a2, b2 = _ieval(v.F.reduce(v.c), v.F.lo, v.F.hi)
            if b1 < a2 or b2 < a1: return False
            u.F.refine(); v.F.refine()
        try:
            from .cluster import _exact_equal            # exact: characteristic-polynomial test
            return _exact_equal(u, v)
        except ArithmeticError:
            return None
    return u == v


# ---------------- polynomials with algebraic coefficients, and several fields at once ----------------
# A polynomial whose coefficients lie in a field Q(alpha) is handled exactly:
#   * Sturm sequences work over any ordered field, so roots are counted and isolated by rational
#     intervals using exact field arithmetic (signs of field elements are exact);
#   * a root rho is represented through the NORM  N(t) = Res_y(m(y), p(y, t))  (m = defining
#     polynomial of alpha), a rational polynomial having every root of p among its roots;
#     the isolating interval is refined until it isolates rho among the roots of N as well;
#   * numbers from two different fields Q(alpha), Q(beta) are combined in Q(gamma),
#     gamma = beta + k alpha (primitive element), with alpha, beta expressed in Q(gamma)
#     through a gcd over Q(gamma).  This is done automatically by the arithmetic of A.

def _coeff_field(p):
    F = None
    for v in p:
        if isinstance(v, A):
            if F is None: F = v.F
            elif v.F is not F: F = unify_fields(F, v.F)
    return F

def _lift_poly(p, F):
    return [embed(v, F) if isinstance(v, A) else Fr(v) for v in p]

def _abs_upper(v):
    if not isinstance(v, A): return abs(Fr(v))
    a, b = _ieval(v.F.reduce(v.c), v.F.lo, v.F.hi)
    return max(abs(a), abs(b))

def _abs_lower(v):
    """a positive lower bound for |v|, v != 0"""
    if not isinstance(v, A): return abs(Fr(v))
    while True:
        a, b = _ieval(v.F.reduce(v.c), v.F.lo, v.F.hi)
        if a > 0: return a
        if b < 0: return -b
        v.F.refine()

def _sturm_generic(p):
    p = squarefree(ptrim(p))
    seq = [p, pderiv(p)]
    while len(seq[-1]) > 1:
        r = pdivmod(seq[-2], seq[-1])[1]
        if not r: break
        seq.append([-v for v in r])
    return p, [t for t in seq if t]

def _count_generic(seq, lo, hi):
    def ch(x):
        sg = [t for t in (sign_of(peval(s, x)) for s in seq) if t != 0]
        return sum(1 for u, v in zip(sg, sg[1:]) if u != v)
    return ch(lo) - ch(hi)

def isolate_real_roots_alg(p, lo=None, hi=None):
    """like isolate_real_roots, for a polynomial with coefficients in one algebraic field"""
    F = _coeff_field(p)
    p = _lift_poly(ptrim(p), F)
    sq, seq = _sturm_generic(p)
    if len(sq) <= 1: return [], sq
    B = 1 + max(_abs_upper(v) for v in sq[:-1]) / _abs_lower(sq[-1])
    L = -B if lo is None else max(-B, Fr(lo))
    H = B if hi is None else min(B, Fr(hi))
    if L >= H: return [], sq
    out, stack = [], [(L, H)]
    while stack:
        a, b = stack.pop()
        n = _count_generic(seq, a, b)
        if n == 0: continue
        if n == 1:
            out.append(("rational", b) if sign_of(peval(sq, b)) == 0 else ("interval", a, b)); continue
        mid = (a + b) / 2
        stack += [(a, mid), (mid, b)]
    return sorted(out, key=lambda r: r[1]), sq

def norm_poly(p):
    """N(t) = Res_y(m(y), p(y, t)): rational, and every root of p is a root of N"""
    from .plane import resultant_d
    F = _coeff_field(p)
    p = _lift_poly(p, F)
    P = {}
    for i, v in enumerate(p):
        cs = v.c if isinstance(v, A) else ([Fr(v)] if v else [])
        for j, c in enumerate(cs):
            if c: P[(i, j)] = P.get((i, j), 0) + c
    M = {(0, j): Fr(c) for j, c in enumerate(F.p) if c}
    return ptrim(resultant_d(P, M))

def roots_alg(p, lo=None, hi=None):
    """exact real roots (Fractions or A) of a polynomial with algebraic coefficients, in (lo, hi]"""
    from .multivar import _rational_roots
    if _coeff_field(p) is None:
        from .cluster import _real_roots
        return _real_roots(p, lo, hi)
    ivs, sq = isolate_real_roots_alg(p, lo, hi)
    if not ivs: return []
    if len(sq) == 2:                                      # linear: the root lies in the same field
        r = -sq[0] / sq[1]
        return [r]
    N = squarefree(norm_poly(sq))
    out = []
    seq = _sturm_generic(sq)[1]
    for iv in ivs:
        if iv[0] == "rational": out.append(iv[1]); continue
        a, b = iv[1], iv[2]
        while count_roots(N, a, b) != 1:                  # refine with p's own Sturm sequence
            mid = (a + b) / 2
            if sign_of(peval(sq, mid)) == 0: a = b = None; out.append(mid); break
            if _count_generic(seq, a, mid) == 1: b = mid
            else: a = mid
        if a is None: continue
        rr = [r for r in _rational_roots(N, a, b) if a < r <= b]
        out.append(rr[0] if rr else A.root(AlgField(N, a, b)))
    return out

# ---- several fields: embeddings into a common field ----
_EMB = {}       # (id(E), id(L)) -> (E, L, image of E's generator in L)

def _record(E, L, img):
    _EMB[(id(E), id(L))] = (E, L, img)

def _image(E, L):
    hit = _EMB.get((id(E), id(L)))
    return hit[2] if hit and hit[0] is E and hit[1] is L else None

def embed(x, L):
    """the element x (Fraction or A) as an element of the field L"""
    if not isinstance(x, A) or x.F is L: return x
    img = _image(x.F, L)
    if img is None: raise NotImplementedError("no embedding known")
    acc = Fr(0)
    for c in reversed(x.F.reduce(x.c)): acc = acc * img + c
    return acc

def unify_fields(F, K):
    """a field containing both F and K (records the embeddings)"""
    if F is K: return F
    if _image(F, K) is not None: return K
    if _image(K, F) is not None: return F
    for (iE, iL), (E, L, img) in list(_EMB.items()):
        if E is F and _image(K, L) is not None: return L
    return _compositum(F, K)

MAX_FIELD_DEGREE = 16

class FieldTooLarge(ArithmeticError):
    pass

def _compositum(F, K):
    from .plane import resultant_d
    if (len(F.p) - 1) * (len(K.p) - 1) > MAX_FIELD_DEGREE:
        raise FieldTooLarge(f"common field of degree {(len(F.p) - 1) * (len(K.p) - 1)} > {MAX_FIELD_DEGREE}")
    for k in range(1, 40):
        # R(t) = Res_y(F.p(y), K.p(t - k y)); gamma = beta + k alpha is a root
        P = {}
        for j, c in enumerate(K.p):
            if not c: continue
            from math import comb
            for e in range(j + 1):          # (t - k y)^j = sum comb(j,e) t^(j-e) (-k y)^e
                key = (j - e, e)
                P[key] = P.get(key, 0) + Fr(c) * comb(j, e) * (-k) ** e
        M = {(0, j): Fr(c) for j, c in enumerate(F.p) if c}
        R = ptrim(resultant_d(P, M))
        if len(R) <= 1 or len(pgcd(R, pderiv(R))) > 1: continue
        # isolate gamma inside K.interval + k * F.interval
        for _ in range(400):
            lo, hi = K.lo + k * F.lo, K.hi + k * F.hi
            if count_roots(R, lo, hi) == 1: break
            F.refine(); K.refine()
        else:
            continue
        L = AlgField(R, lo, hi)
        g = A.root(L)
        # alpha = the common root of F.p(y) and K.p(g - k y) (a gcd over L)
        G2 = []
        for j, c in enumerate(K.p):
            if not c: continue
            term = [Fr(c)]
            for _ in range(j): term = pmul(term, [g, Fr(-k)])
            G2 = padd(G2, term)
        G = pgcd([Fr(c) for c in F.p], G2)
        if len(G) != 2: continue
        alpha = -G[0] / G[1]
        beta = g - k * alpha
        _record(F, L, alpha); _record(K, L, beta)
        # anything already embedded in F or K is embedded in L too
        for (iE, iX), (E, X, img) in list(_EMB.items()):
            if X is F or X is K:
                if _image(E, L) is None:
                    _record(E, L, embed(img, L) if isinstance(img, A) else img)
        return L
    raise ArithmeticError("no primitive element found")
