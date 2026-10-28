"""
fmsl.engine -- the Fluid Transition Engine.

Every value is a node in an expression DAG. Its behavior at the transition point
x = 0 is a LAZY Laurent series with exact coefficients in Q(i) (Gaussian rationals),
computed on demand.

Certificate (proof that "no nonzero term exists"):
  A certified value is E/Q with Q a polynomial (degree <= dQ) and E in
      span{ x^j * e^(q(x)) : q in spectrum, j < multiplicity(q) }
  where each q is a polynomial exponent with q(0) = 0 (sin and cos use +-i*q).
  These basis functions are linearly independent, so the map "E -> its first M+1
  Taylor coefficients" is injective for a finite M:
    * linear exponents (e^(lam x), sin, cos, polynomials): M = N - 1, because E
      solves a constant-coefficient linear ODE of order N;
    * polynomial exponents (e^(x^2), sin(x^3), ...): smallest M found exactly by
      incremental Gaussian elimination over Q(i).
  Hence E != 0  =>  ord(E) <= M, and scanning a finite window either finds the
  exact leading term or proves the value is identically zero (C_e).

Outside the certified class (e.g. exp(sin x), division by sin x) coefficients are
still exact, but "is it identically zero?" is reported as UNPROVEN -- never guessed.
(Richardson's theorem: zero-testing is undecidable for general expressions.)

Essential singularities (exp, sin, cos of a pole) have no Laurent series; they are
detected exactly and reported with their real-line behavior.
"""
from fractions import Fraction as Fr
from math import factorial, atan2, degrees
from .gauss import G

DEPTH = 60            # scan depth for uncertified values
_stats = {"coeffs": 0}

class Essential(Exception):
    pass

class V:
    __slots__ = ("kind", "a", "b", "lo", "_r", "spec", "pdeg", "dQ", "memo", "_val", "data")

    def __init__(s, kind, a=None, b=None, lo=0, spec=None, pdeg=None, dQ=0, data=None):
        s.kind, s.a, s.b, s.lo, s.spec, s.pdeg, s.dQ = kind, a, b, lo, spec, pdeg, dQ
        s.data, s.memo, s._val = data, {}, "?"
        s._r = "?"

    @property
    def r(s):
        """certificate: if E != 0 then ord(E) <= r - 1   (None = uncertified)"""
        if s._r == "?":
            if s.spec is None: s._r = None
            else:
                N = sum(s.spec.values())
                if all(len(q) <= 1 for q in s.spec): s._r = N          # linear: ODE bound
                else: s._r = _rank_bound(s.spec) + 1                   # polynomial exponents
                if s.pdeg is not None: s._r = min(s._r, s.pdeg + 1)
        return s._r

    # ----------------------------- constructors -----------------------------
    @staticmethod
    def rat(N, D=(1,)):
        """N/D with polynomial coefficient lists (index = power of x)."""
        N = [G(c) for c in N]; D = [G(c) for c in D]
        while N and N[-1] == 0: N.pop()
        while D and D[-1] == 0: D.pop()
        if not D: raise ZeroDivisionError("denominator identically zero")
        vd = next(i for i, c in enumerate(D) if c != 0)
        vn = next((i for i, c in enumerate(N) if c != 0), None)
        lo = (vn - vd) if vn is not None else 0
        return V("rat", lo=lo, spec=({ZERO: len(N)} if N else {}), pdeg=len(N) - 1, dQ=len(D) - 1,
                 data=(N, D[vd:], vd))

    @staticmethod
    def exp(lam=1):
        lam = G(lam); return V("exp", lo=0, spec={(lam.pair(),): 1}, data=lam)
    @staticmethod
    def sin(lam=1):  return V("sin", lo=1, spec=_trig_spec([G(lam)]), data=G(lam))
    @staticmethod
    def cos(lam=1):  return V("cos", lo=0, spec=_trig_spec([G(lam)]), data=G(lam))

    # ----------------------------- arithmetic -------------------------------
    def __add__(s, o):
        if isinstance(o, EssFactor): return NotImplemented
        o = _lift(o)
        cert = s.spec is not None and o.spec is not None
        spec = _union(_shift(s.spec, o.dQ), _shift(o.spec, s.dQ)) if cert else None
        pdeg = max(s.pdeg + o.dQ, o.pdeg + s.dQ) if cert and s.pdeg is not None and o.pdeg is not None else None
        return V("+", s, o, min(s.lo, o.lo), spec, pdeg, s.dQ + o.dQ)
    __radd__ = __add__
    def __neg__(s):    return V("neg", s, None, s.lo, s.spec, s.pdeg, s.dQ)
    def __sub__(s, o):
        if isinstance(o, EssFactor): return NotImplemented
        return s + (-_lift(o))
    def __rsub__(s, o): return _lift(o) - s
    def __mul__(s, o):
        if isinstance(o, EssFactor): return NotImplemented
        o = _lift(o)
        cert = s.spec is not None and o.spec is not None
        pdeg = s.pdeg + o.pdeg if cert and s.pdeg is not None and o.pdeg is not None else None
        spec = _prod(s.spec, o.spec) if cert else None
        return V("*", s, o, s.lo + o.lo, spec, pdeg, s.dQ + o.dQ)
    __rmul__ = __mul__
    def __truediv__(s, o):
        if isinstance(o, EssFactor): return NotImplemented
        o = _lift(o)
        vb = o.val()
        if vb[0] == "zero":
            raise ZeroDivisionError("divisor is exactly C_e" if vb[1] else
                                    f"divisor indistinguishable from 0 through order {vb[2]}")
        single = (o.spec is not None and len(o.spec) == 1 and list(o.spec.values())[0] == 1
                  and list(o.spec)[0] != ZERO)
        if s.spec is not None and o.spec is not None and o.pdeg is not None:
            spec = _shift(s.spec, o.dQ)                      # divide by polynomial numerator
            pdeg = s.pdeg + o.dQ if s.pdeg is not None else None
            return V("/", s, o, s.lo - vb[1], spec, pdeg, s.dQ + o.pdeg, data=vb[1])
        if s.spec is not None and single:                    # divide by c * e^{q}: multiply by e^{-q}
            q = list(o.spec)[0]
            spec = _prod(_shift(s.spec, o.dQ), {_kneg(q): 1})
            return V("/", s, o, s.lo - vb[1], spec, None, s.dQ, data=vb[1])
        return V("/", s, o, s.lo - vb[1], None, None, s.dQ, data=vb[1])
    def __rtruediv__(s, o):
        if isinstance(o, EssFactor): return NotImplemented
        return _lift(o) / s
    def __pow__(s, n):
        out = V.rat([1])
        for _ in range(n): out = out * s
        return out

    # ------------------------ composition (uncertified) ---------------------
    @staticmethod
    def compose(fn, f, spec=None, keep=False):
        if isinstance(f, EssFactor): raise Essential(f._undecided(f"inside {fn}(...)"))
        f = _lift(f)
        neg = [(k, f.coeff(k)) for k in range(f.lo, 0) if f.coeff(k) != 0]
        if neg:
            m, c = -neg[0][0], neg[0][1]          # inner pole of order m, coeff c
            if keep and m == -f.lo and all(f.coeff(k) == 0 for k in range(f.lo + 1, 0)) and f.coeff(0) == 0:
                # f = c x^-m + (vanishing terms): keep fn(c x^-m + ...) as an essential factor,
                # so that a product with a regular factor can still be decided (squeeze theorem)
                return EssFactor(fn, m, c, f)
            raise Essential(_ess_text(fn, m, c))
        if f.coeff(0) != 0:
            raise ValueError(f"{fn}(f) with f(0) = {f.coeff(0)} != 0 has transcendental "
                             "coefficients (e.g. e^c); outside exact rational arithmetic")
        core = {"f": f, "s": {}, "c": {}, "e": {}}
        kind = {"exp": "cexp", "sin": "csin", "cos": "ccos"}[fn]
        return V(kind, lo=0, spec=spec, data=core)

    # --------------------------- lazy coefficients --------------------------
    def coeff(s, k):
        if k < s.lo: return G(0)
        if k in s.memo: return s.memo[k]
        _stats["coeffs"] += 1
        K = s.kind
        if K == "rat":
            N, d, vd = s.data
            j = k + vd                                   # coefficient j of N/d
            acc = (N[j] if j < len(N) else 0) - sum(
                d[i] * s.coeff(k - i) for i in range(1, min(j, len(d) - 1) + 1))
            r = acc / d[0]
        elif K == "exp": r = s.data ** k / factorial(k)
        elif K == "sin": r = (s.data ** k / factorial(k)) * (0, 1, 0, -1)[k % 4]
        elif K == "cos": r = (s.data ** k / factorial(k)) * (1, 0, -1, 0)[k % 4]
        elif K == "neg": r = -s.a.coeff(k)
        elif K == "+":   r = s.a.coeff(k) + s.b.coeff(k)
        elif K == "*":
            r = sum(s.a.coeff(i) * s.b.coeff(k - i) for i in range(s.a.lo, k - s.b.lo + 1))
        elif K == "/":
            vb = s.data; lead = s.b.coeff(vb)
            r = (s.a.coeff(k + vb) - sum(s.b.coeff(vb + i) * s.coeff(k - i)
                                         for i in range(1, k - s.lo + 1))) / lead
        elif K in ("cexp", "csin", "ccos"):
            r = _comp_coeff(s.data, K[1:], k)
        elif K == "clog":                                # L = log f, f(0) = 1:  L' f = f'
            f, L = s.data["f"], s.data["L"]
            for n in range(1, k + 1):
                if n not in L:
                    L[n] = (n * f.coeff(n) - sum(j * L[j] * f.coeff(n - j) for j in range(1, n))) / n
            r = G(0) if k == 0 else L[k]
        elif K == "cpow":                                # g = f^a, f(0) = c:  g' f = a f' g
            f, g, a, c = s.data["f"], s.data["g"], s.data["a"], s.data["c"]
            for n in range(1, k + 1):
                if n not in g:
                    g[n] = sum((a * j - (n - j)) * f.coeff(j) * g[n - j] for j in range(1, n + 1)) / (n * c)
            r = g[k]
        s.memo[k] = r
        return r

    # ------------------------------ queries ---------------------------------
    def val(s):
        """('ord', k, certified) | ('zero', certified, scanned_to)"""
        if s._val != "?": return s._val
        lo = max(s.lo, -s.dQ) if s.r is not None else s.lo
        hi = s.r - 1 if s.r is not None else s.lo + DEPTH
        res = ("zero", s.r is not None, hi)
        for k in range(lo, hi + 1):
            if s.coeff(k) != 0:
                res = ("ord", k, s.r is not None); break
        s._val = res
        return res

    def describe(s):
        v = s.val()
        if v[0] == "zero":
            return ("C_e  (certified: identically zero)" if v[1]
                    else f"zero through order {v[2]}  (UNCERTIFIED - no proof it is exactly 0)")
        k, c = v[1], s.coeff(v[1])
        tag = "" if v[2] else "  [uncertified class, but this term is exact]"
        if k > 0:  return f"removable, value 0 (vanishes to order {k}, lead {c}){tag}"
        if k == 0: return f"removable, value {c}{tag}"
        return f"pole of order {-k}: {infinity_flow(k, c)} (lead {c}){tag}"

    def certified(s): return s.r is not None
    def window(s):    return None if s.r is None else s.r - 1 - max(s.lo, -s.dQ) + 1


class EssFactor:
    """R(x) * fn(f(x)), where f has a pole of order m with leading coefficient c at 0 and R is an
       ordinary (Laurent) value.  Only products and quotients with ordinary values are kept;
       anything else (sums, two essential factors, coefficients) raises Essential as before.
       limits() decides the one-sided limits rigorously when that is possible:
         sin / cos of a pole with real lead: bounded by 1 and oscillating; times R vanishing to
           order k > 0 -> 0 on both sides (squeeze theorem); otherwise no limit;
         exp of a pole with real lead: on a side where the exponent -> -oo it is flat (beats every
           power), so the product -> 0; where it -> +oo it beats every power, so the product
           -> +-oo with the sign of R on that side."""
    def __init__(s, fn, m, c, f, R=None):
        s.fn, s.m, s.c, s.f = fn, m, G(c), f
        s.R = R if R is not None else V.rat([1])
    def _text(s): return _ess_text(s.fn, s.m, s.c)
    def _undecided(s, how):
        return (f"ESSENTIAL singularity: {s.fn} of a pole of order {s.m} (lead {s.c}) {how}; "
                f"the limit is not decided [UNPROVEN]")
    def _fail(s, *a): raise Essential(s._undecided("inside a sum"))
    def val(s): raise Essential(s._text())
    def coeff(s, k): raise Essential(s._text())
    def __mul__(s, o):
        if isinstance(o, EssFactor): raise Essential(s._undecided("times another essential factor"))
        return EssFactor(s.fn, s.m, s.c, s.f, s.R * _lift(o))
    __rmul__ = __mul__
    def __truediv__(s, o):
        if isinstance(o, EssFactor): raise Essential(s._undecided("divided by another essential factor"))
        return EssFactor(s.fn, s.m, s.c, s.f, s.R / _lift(o))
    def __rtruediv__(s, o): raise Essential(s._undecided("in a denominator"))
    def __neg__(s): return EssFactor(s.fn, s.m, s.c, s.f, -s.R)
    __add__ = __radd__ = __sub__ = __rsub__ = _fail
    def __pow__(s, n):
        if n == 1: return s
        raise Essential(s._undecided(f"raised to the power {n}"))

    def limits(s):
        """(left, right, text) with left/right in {Fraction 0, '+oo', '-oo', 'osc'}, or None"""
        if not s.c.is_real(): return None
        rv = s.R.val()
        if rv[0] == "zero":
            if not rv[1]: return None
            return (Fr(0), Fr(0), "identically 0 (the factor in front is C_e) [proven]")
        k, a = rv[1], s.R.coeff(rv[1])
        if not G(a).is_real(): return None
        a = G(a).re
        if s.fn in ("sin", "cos"):
            if k > 0:
                return (Fr(0), Fr(0), f"limit 0 on both sides: |{s.fn}| <= 1 while the factor in front "
                                      f"vanishes to order {k} (squeeze theorem) [proven]")
            return ("osc", "osc", f"{s.fn} of a pole of order {s.m} oscillates"
                                  + (" with growing amplitude" if k < 0 else "") + "; no limit on either side [proven]")
        # exp
        out = []
        for sg in (-1, 1):
            expo = s.c.re * (sg ** s.m)                  # sign of the exponent c x^-m on this side
            if expo < 0: out.append(Fr(0))
            else:
                side = a * (sg ** k) if k >= 0 else a * (sg ** (-k))
                out.append("+oo" if side > 0 else "-oo")
        txt = (f"exp of a pole of order {s.m} (lead {s.c}) times a factor of order {k}: "
               f"x->0-: {out[0]}; x->0+: {out[1]} (exp beats every power) [proven]")
        return (out[0], out[1], txt)


ZERO = ()          # spectrum key of the zero exponent (plain polynomials)

def _kstrip(k):
    k = list(k)
    while k and k[-1] == (0, 0): k.pop()
    return tuple(k)

def _kadd(p, q):
    n = max(len(p), len(q)); z = (Fr(0), Fr(0))
    return _kstrip(tuple((a[0] + b[0], a[1] + b[1]) for a, b in
                         zip(p + (z,) * (n - len(p)), q + (z,) * (n - len(q)))))

def _kneg(q): return tuple((-a, -b) for a, b in q)

def _trig_spec(coeffs):
    """sin/cos of q(x) = sum coeffs[k-1] x^k : exponents +i q and -i q"""
    ip = _kstrip(tuple((-G(c).im, G(c).re) for c in coeffs))
    return {ip: 1, _kneg(ip): 1}

# ---- exact complex-rational arithmetic for the rank certificate ----
def _cm(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def _cinv(a):
    d = a[0] * a[0] + a[1] * a[1]; return (a[0] / d, -a[1] / d)

def _exp_series(q, n):
    """Taylor coefficients 0..n of e^{q(x)}, q complex-rational, q(0)=0: k g_k = sum j q_j g_{k-j}"""
    g = [(Fr(1), Fr(0))]
    for k in range(1, n + 1):
        acc = (Fr(0), Fr(0))
        for j in range(1, min(k, len(q)) + 1):
            t = _cm(q[j - 1], g[k - j]); acc = (acc[0] + j * t[0], acc[1] + j * t[1])
        g.append((acc[0] / k, acc[1] / k))
    return g

_rank_cache = {}
def _rank_bound(spec):
    """Smallest M such that the Taylor coefficients 0..M of the basis
       {x^j e^{q}: q in spec, j < mult} have full rank N. Exact, terminates by independence."""
    key = tuple(sorted(spec.items()))
    if key in _rank_cache: return _rank_cache[key]
    basis = [(q, j) for q, m in spec.items() for j in range(m)]
    N = len(basis)
    if N == 0: _rank_cache[key] = -1; return -1
    series = {}
    def col(q, j, n):
        if q not in series or len(series[q]) <= n:
            cap = max(2 * len(series.get(q, [])), n + 1, 16)   # per-exponent growth
            series[q] = _exp_series(q, cap)
        return series[q][n - j] if n >= j else (Fr(0), Fr(0))
    pivots = []                     # list of (pivot column, row vector) in echelon form
    n = -1
    while len(pivots) < N:
        n += 1
        row = [col(q, j, n) for q, j in basis]
        for pc, pr in pivots:
            c = row[pc]
            if c != (0, 0):
                row = [(a[0] - t[0], a[1] - t[1]) for a, t in zip(row, (_cm(c, b) for b in pr))]
        nz = next((i for i, a in enumerate(row) if a != (0, 0)), None)
        if nz is not None:
            inv = _cinv(row[nz]); row = [_cm(inv, a) for a in row]
            pivots.append((nz, row))
    _rank_cache[key] = n
    return n

def _shift(spec, d):
    """multiply by a polynomial of degree d: every root gains multiplicity d"""
    return {k: m + d for k, m in spec.items()}

def _union(a, b):
    out = dict(a)
    for k, m in b.items(): out[k] = max(out.get(k, 0), m)
    return out

def _prod(a, b):
    """x^i e^{p x} * x^j e^{q x} = x^(i+j) e^{(p+q) x}: roots add, multiplicities m1+m2-1"""
    out = {}
    for p, m1 in a.items():
        for q, m2 in b.items():
            k = _kadd(p, q)
            out[k] = max(out.get(k, 0), m1 + m2 - 1)
    return out

def _lift(o):
    return o if isinstance(o, V) else V.rat([o])

def _comp_coeff(core, which, k):
    """Coefficients of exp(f), sin(f), cos(f) for f(0)=0, via the ODE recurrences
       g' = f' g ;  s' = f' c ;  c' = -f' s."""
    f, S, C, E = core["f"], core["s"], core["c"], core["e"]
    def e(n):
        if n not in E:
            E[n] = Fr(1) if n == 0 else sum(j * f.coeff(j) * e(n - j) for j in range(1, n + 1)) / n
        return E[n]
    def sn(n):
        if n not in S:
            S[n] = Fr(0) if n == 0 else sum(j * f.coeff(j) * cs(n - j) for j in range(1, n + 1)) / n
        return S[n]
    def cs(n):
        if n not in C:
            C[n] = Fr(1) if n == 0 else -sum(j * f.coeff(j) * sn(n - j) for j in range(1, n + 1)) / n
        return C[n]
    for n in range(k):            # fill iteratively to avoid deep recursion
        e(n) if which == "exp" else (sn(n), cs(n))
    return e(k) if which == "exp" else (sn(k) if which == "sin" else cs(k))

def infinity_flow(k, c):
    """Infinity Flow (after Celedon's operator I): the exact direction in which a pole
       runs to infinity from each side. Near 0, h ~ c x^k with k < 0."""
    c = G(c)
    if c.is_real():
        left = "+inf" if (c > 0) == ((-k) % 2 == 0) else "-inf"
        right = "+inf" if c > 0 else "-inf"
        return f"left {left}, right {right}"
    ang = degrees(atan2(float(c.im), float(c.re)))
    angl = (ang + (180 if k % 2 else 0) + 180) % 360 - 180
    return (f"|value| -> inf like |lead|/|x|^{-k}, heading in direction {ang:.1f} deg "
            f"from the right and {angl:.1f} deg from the left")

def _ess_text(fn, m, c):
    c = G(c)
    if not c.is_real():
        if fn == "exp" and c.re == 0:
            return (f"ESSENTIAL singularity: exp of a purely imaginary pole of order {m} (lead {c}); "
                    f"|value| stays 1 on the real line while the phase spins infinitely fast; "
                    f"no Laurent series exists")
        return (f"ESSENTIAL singularity: {fn} of a complex pole of order {m} (lead {c}); "
                f"no limit on the real line; no Laurent series exists")
    if fn == "exp":
        right = "+inf" if c > 0 else "0 (flat: every derivative -> 0)"
        lsign = c if m % 2 == 0 else -c
        left = "+inf" if lsign > 0 else "0 (flat: every derivative -> 0)"
        return (f"ESSENTIAL singularity: exp of a pole of order {m} (lead {c}); "
                f"x->0-: {left}; x->0+: {right}; no Laurent series exists")
    return (f"ESSENTIAL singularity: {fn} of a pole of order {m} (lead {c}); "
            f"oscillates without limit on both sides; no Laurent series exists")

# convenient atoms
x = V.rat([0, 1])
# keep=True: an essential singularity inside is kept as an EssFactor (used by the FMSL language,
# so that products such as x*sin(1/x) can be decided); keep=False raises Essential at once.
def Exp(f, keep=False):  return V.exp(f) if not isinstance(f, (V, EssFactor)) else (_atom_or_comp("exp", f, keep))
def Sin(f, keep=False):  return V.sin(f) if not isinstance(f, (V, EssFactor)) else (_atom_or_comp("sin", f, keep))
def Cos(f, keep=False):  return V.cos(f) if not isinstance(f, (V, EssFactor)) else (_atom_or_comp("cos", f, keep))

def _rational_power(c, a):
    """c^a for rational c > 0 and rational a, if it is rational; else None"""
    c, a = Fr(c), Fr(a)
    def root(n, q):
        if n < 0: return None
        r = round(n ** (1.0 / q)) if n else 0
        for t in (r - 1, r, r + 1):
            if t >= 0 and t ** q == n: return t
        return None
    p, q = a.numerator, a.denominator
    num, den = root(c.numerator, q), root(c.denominator, q)
    if num is None or den is None: return None
    return Fr(num, den) ** p


def _lead_real(f):
    """f(0) as a real rational, after checking f has no pole and real f(0)"""
    if isinstance(f, EssFactor): raise Essential(f._undecided("inside a power or log"))
    neg = [(k, f.coeff(k)) for k in range(f.lo, 0) if f.coeff(k) != 0]
    if neg: raise Essential(f"the argument has a pole of order {-neg[0][0]}: not a power series")
    c0 = f.coeff(0)
    if hasattr(c0, "is_real") and not c0.is_real(): raise ValueError("complex values are not supported here")
    return Fr(c0.re) if hasattr(c0, "re") else Fr(c0)


def Log(f):
    """log f, for f tending to 1 (exact power series). f -> c != 1 would need the transcendental
       constant log c; f -> 0 is a logarithmic singularity, not a power series."""
    f = _lift(f)
    c = _lead_real(f)
    if c == 0: raise Essential("logarithmic singularity: log of a quantity tending to 0 (not a power series)")
    if c != 1:
        raise ValueError(f"log(f) with f(0) = {c} != 1 has the transcendental constant log({c}); "
                         "write it as log(f/c) + log(c) with f/c -> 1, or move the point")
    return V("clog", lo=0, data={"f": f, "L": {}})


def PowR(f, a):
    """f^a for a rational exponent a, for f tending to c > 0 with c^a rational (e.g. c = 1)."""
    a = Fr(a)
    if a.denominator == 1:
        n = int(a); b = _lift(f)
        return b ** n if n >= 0 else 1 / (b ** -n)
    f = _lift(f)
    c = _lead_real(f)
    if c == 0: raise Essential("branch point: a fractional power of a quantity tending to 0 (not a power series)")
    if c < 0: raise ValueError("a fractional power of a negative quantity is not real")
    g0 = _rational_power(c, a)
    if g0 is None:
        raise ValueError(f"{c}^{a} is irrational; outside exact rational arithmetic")
    return V("cpow", lo=0, data={"f": f, "g": {0: G(g0)}, "a": a, "c": c})


def _atom_or_comp(fn, f, keep=False):
    if isinstance(f, EssFactor): raise Essential(f._undecided(f"inside {fn}(...)"))
    """Certified atom when f is exactly lam*x; certified composition when f is a
       polynomial with f(0)=0; otherwise uncertified composition."""
    # certified polynomial (pdeg known, no denominator): read its exact coefficients
    if f.r is not None and f.pdeg is not None and f.dQ == 0:
        c = [f.coeff(k) for k in range(0, max(f.pdeg, 1) + 1)]
        if all(v == 0 for v in c):                       # argument is exactly 0
            return V.rat([0 if fn == "sin" else 1])
        if c[0] == 0 and all(v == 0 for v in c[2:]) and c[1] != 0:
            return {"exp": V.exp, "sin": V.sin, "cos": V.cos}[fn](c[1])
        if c[0] == 0:                                     # polynomial exponent q(x)
            spec = ({_kstrip(tuple(G(v).pair() for v in c[1:])): 1} if fn == "exp"
                    else _trig_spec(c[1:]))
            return V.compose(fn, f, spec, keep=keep)
    return V.compose(fn, f, keep=keep)
