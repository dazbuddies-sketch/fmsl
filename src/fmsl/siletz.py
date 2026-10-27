"""
fmsl.siletz -- Siletz bounds: PROVEN error bounds for the emitted series branch.

Named after the Siletz Precision System in MJ Celedon's FMS ("structured precision
between defined values"). This is its rigorous counterpart: every quantity is an
exact rational enclosure with outward rounding, so the final number is a proof.

Method
  1. Circle bound B(h, r): an upper bound for |h(z)| on the complex circle |z| = r,
     computed recursively over the expression graph:
        |a+b| <= B(a)+B(b)   |a*b| <= B(a)B(b)   |e^w| <= e^|w|   |sin w| <= sinh|w|
        |cos w| <= cosh|w|   |a/b| <= B(a)/L(b)
     L(b, r) is a lower bound for |b| from its dominant term:
        |b(z)| >= |z|^v (|e0| - sum_{j>=1} |e_j| r^j)   for all 0 < |z| <= r,
     which also PROVES b has no zeros in the punctured disk (so h is analytic there).
     The tail of that sum is itself bounded by Cauchy's estimate at a larger radius.
  2. Cauchy's estimate: with g = h * x^(-lo) analytic in |z| <= r,
        |d_j| <= B(g, r) / r^j   for every coefficient not yet computed.
     So the truncation error of K terms on |x| <= h0 < r is at most
        B(g, r) (h0/r)^K / (1 - h0/r).
  3. Floating point: coefficient rounding (u per coefficient) plus Higham's Horner
     bound gamma_{2K} * sum |d_j| |x|^j, plus gamma_m for the x^lo multiplications.
Everything is relative to a proven lower bound of |g| on |x| <= h0.

Scope: the bound covers the series branch |x| < h0 (the dangerous region). The direct
formula used for |x| >= h0 calls the platform math library, which gives no formal
accuracy guarantee, so that branch remains measured, not proven.
"""
from fractions import Fraction as Fr
from math import ceil, floor
from .gauss import G

U = Fr(1, 2 ** 53)            # unit roundoff of IEEE double
NTERMS = 24                   # exact terms used in dominant-term lower bounds


def _grid(q):
    """grid spacing ~ 2^-40 RELATIVE to |q| (about 12 significant digits kept)"""
    e = q.numerator.bit_length() - q.denominator.bit_length()
    return Fr(2) ** (e - 40)

def up(q):
    """round a nonnegative rational UP to ~40 significant bits (keeps sizes small, stays rigorous)"""
    if q <= 0: return Fr(0) if q == 0 else -down(-q)
    g = _grid(q); return ceil(q / g) * g

def down(q):
    if q <= 0: return Fr(0) if q == 0 else -up(-q)
    g = _grid(q); return floor(q / g) * g

def gamma(n, u):
    return n * u / (1 - n * u)

def _series_tail_factor(y, n):
    # sum_{k>=n} y^k/k! <= y^n/n! * 1/(1 - y/(n+1)),  valid for y < n+1
    term = Fr(1)
    for k in range(1, n + 1): term = term * y / k
    return term * (n + 1) / (n + 1 - y)

def _exp_family_up(y, parity):
    """rigorous upper bound of e^y (parity None), cosh y (0), sinh y (1) for rational y >= 0"""
    y = up(y)
    n = max(20, 2 * ceil(y) + 10)
    s, term = Fr(0), Fr(1)
    for k in range(n):
        if parity is None or k % 2 == parity: s += term
        term = term * y / (k + 1)
    return up(s + _series_tail_factor(y, n))

def exp_up(y):  return _exp_family_up(y, None)
def cosh_up(y): return _exp_family_up(y, 0)
def sinh_up(y): return _exp_family_up(y, 1)

def _absu(c): return G(c).abs_upper()
def _absl(c): return G(c).abs_lower()


class Siletz:
    def __init__(s):
        s.memo = {}

    # ---------------- circle bound ----------------
    def B(s, v, r):
        """upper bound of |v(z)| on |z| = r; None if analyticity cannot be proven"""
        key = (id(v), r)
        if key in s.memo: return s.memo[key]
        s.memo[key] = None                      # guards against pathological recursion
        res = s._B(v, r)
        s.memo[key] = None if res is None else up(res)
        return s.memo[key]

    def _B(s, v, r):
        K = v.kind
        if K == "rat":
            N, d, vd = v.data
            num = sum((_absu(c) * r ** k for k, c in enumerate(N)), Fr(0))
            low = _absl(d[0]) - sum((_absu(c) * r ** k for k, c in enumerate(d) if k >= 1), Fr(0))
            if low <= 0: return None            # cannot prove the denominator zero-free
            return num / (r ** vd * low)
        if K == "exp": return exp_up(_absu(v.data) * r)
        if K == "sin": return sinh_up(_absu(v.data) * r)
        if K == "cos": return cosh_up(_absu(v.data) * r)
        if K == "neg": return s.B(v.a, r)
        if K in ("+", "*"):
            a, b = s.B(v.a, r), s.B(v.b, r)
            if a is None or b is None: return None
            return a + b if K == "+" else a * b
        if K == "/":
            a, lb = s.B(v.a, r), s.L(v.b, r)
            if a is None or lb is None: return None
            return a / lb
        if K in ("cexp", "csin", "ccos"):
            f = v.data["f"]
            bf = s.B(f, r)
            if bf is None: return None
            return {"cexp": exp_up, "csin": sinh_up, "ccos": cosh_up}[K](bf)
        return None

    # ---------------- dominant-term lower bound ----------------
    def L(s, b, r):
        """lower bound of |b(z)| on |z| = r, valid with proof that b != 0 on 0 < |z| <= r"""
        val = b.val()
        if val[0] != "ord": return None
        v = val[1]
        e0 = _absl(b.coeff(v))
        head = sum((_absu(b.coeff(v + j)) * r ** j for j in range(1, NTERMS + 1)), Fr(0))
        tail = None
        for R in (2 * r, 3 * r, Fr(3, 2) * r):
            BR = s.B(b, R)
            if BR is None: continue
            Mg = BR / R ** v if v >= 0 else BR * R ** (-v)          # bound for b * z^(-v) on |z| = R
            q = r / R
            tail = Mg * q ** (NTERMS + 1) / (1 - q)
            break
        if tail is None: return None
        low = e0 - head - tail
        if low <= 0: return None
        return down(low) * (r ** v if v >= 0 else 1 / r ** (-v))


def _sqrt_down(q):
    """rational lower bound of sqrt(q), q >= 0"""
    from math import isqrt
    if q <= 0: return Fr(0)
    sc = 2 ** 60
    return Fr(isqrt(q.numerator * sc * sc // q.denominator), sc)

def _poly_at(d, K, x):
    acc = G(0)
    for c in reversed(d[:K]): acc = acc * x + c
    return acc


def proven_series_error(v, lo, coeffs, K, h0, complex_mode, pow_mults, max_depth=12):
    """Proven bound on the relative error of the emitted series branch for real 0 < |x| < h0.
       The lower bound on |g(x)| is proven on the REAL interval by subdivision: on [a, b],
       |p(x)| >= |p(m)| - sup|p'| (b - a)/2, with p evaluated exactly at the midpoint m.
       pow_mults: floating multiplications/divisions used for x^lo.
       Returns (bound as Fraction, radius used) or (None, None)."""
    h0 = Fr(h0)
    u = 2 * U if complex_mode else U
    d = coeffs
    absd = [_absu(c) for c in d[:K]]
    sil = Siletz()
    e_pow = gamma(pow_mults, u) if pow_mults else Fr(0)
    g2K = gamma(2 * K, u)
    best = (None, None)
    for f in (Fr(5, 4), Fr(3, 2), Fr(2), Fr(3), Fr(4), Fr(6), Fr(8), Fr(12), Fr(16)):
        r = up(h0 * f)
        Bh = sil.B(v, r)
        if Bh is None: continue
        Mg = Bh / r ** lo if lo >= 0 else Bh * r ** (-lo)

        def piece(a, b, depth):
            """worst relative error bound on the real segment [a, b] (|a| <= |b|), or None"""
            bm = max(abs(a), abs(b))
            q = bm / r
            T = Mg * q ** K / (1 - q)
            S = sum((absd[j] * bm ** j for j in range(K)), Fr(0))
            D = sum((j * absd[j] * bm ** (j - 1) for j in range(1, K)), Fr(0))
            m = (a + b) / 2
            pm = _poly_at(d, K, m)
            L = _sqrt_down(pm.re * pm.re + pm.im * pm.im) - D * abs(b - a) / 2 - T
            if L <= 0:
                if depth >= max_depth: return None
                l1 = piece(a, m, depth + 1)
                if l1 is None: return None
                l2 = piece(m, b, depth + 1)
                if l2 is None: return None
                return max(l1, l2)
            e_trunc = T / L
            e_round = (u + g2K * (1 + u)) * S / L
            return up((1 + e_trunc) * (1 + e_round) * (1 + e_pow) - 1)

        worst = Fr(0)
        n0 = 8
        ok = True
        for side in (1, -1):
            for i in range(n0):
                a, b = side * h0 * i / n0, side * h0 * (i + 1) / n0
                e = piece(a, b, 0)
                if e is None: ok = False; break
                worst = max(worst, e)
            if not ok: break
        if not ok: continue
        if best[0] is None or worst < best[0]: best = (worst, r)
    return best
