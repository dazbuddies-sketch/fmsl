"""
fmsl.multivar -- fluid states of two-variable rational functions at the origin.

In one variable a 0/0 point has ONE state. In two variables the state depends on the
path of approach. This module computes the full state PROFILE, exactly:

  * Every real analytic path into the origin is (Puiseux) y = c x^a + higher terms,
    the axis x = 0, or a refinement of these. Writing a = p/q and substituting
    x = sigma s^q, y = c s^p (s -> 0+, sigma = +-1) turns a path into a ONE-variable
    problem that the 1-D engine solves exactly.
  * Weighted orders: monomial x^i y^j has weight i + a j. The lowest-weight terms of
    N and D (their "initial forms") decide the state along y = c x^a:
        delta(a) = ord_a(N) - ord_a(D)
        delta > 0 -> 0,   delta < 0 -> infinity,   delta = 0 -> in_N(c) / in_D(c)
  * The state only changes at CRITICAL curvatures a, where two monomials tie in weight
    (the slopes of the Newton polygons): the balance points. Between them the state is
    determined by a single monomial on each side.
  * Where an initial form has a real root c0, the leading terms cancel along curves
    tangent to y = c0 x^a; the engine zooms in (Newton-Puiseux step: y = x^a (c0 + w))
    and repeats, up to a depth limit. Rational roots are refined exactly; irrational
    real roots (detected with Sturm sequences) are reported as unresolved.

Verdict:
  * two paths with different exact states  -> NO LIMIT (proven)
  * all paths agree, nothing unresolved     -> LIMIT = L (proven, for rational functions)
  * otherwise                               -> reported as unproven, with the open curves.
"""
from fractions import Fraction as Fr
from math import gcd, isqrt
from itertools import product

from .compiler import Parser
from .gauss import G
from .algebraic import A, AlgField, isolate_real_roots, sign_of, same_value

MAX_DEPTH = 6
WORK_BUDGET = 1500      # max path records per function; beyond it, regions are reported UNPROVEN

# ------------------------------ bivariate polynomials ------------------------------
def padd(P, Q):
    R = dict(P)
    for k, v in Q.items():
        R[k] = R.get(k, 0) + v
        if R[k] == 0: del R[k]
    return R

def pneg(P): return {k: -v for k, v in P.items()}

def pmul(P, Q):
    R = {}
    for (i1, j1), a in P.items():
        for (i2, j2), b in Q.items():
            k = (i1 + i2, j1 + j2)
            R[k] = R.get(k, 0) + a * b
    return {k: v for k, v in R.items() if v != 0}

ONE = {(0, 0): Fr(1)}

def ppow(P, n):
    R = ONE
    for _ in range(n): R = pmul(R, P)
    return R

def parse_rational(text):
    """Parse an expression in x, y with + - * / ^ into (N, D) polynomial dicts."""
    ast = Parser(text).parse()
    def ev(e):
        t = e[0]
        if t == "num":
            v = e[1]
            if isinstance(v, G):
                if not v.is_real(): raise ValueError("multivar: real coefficients only")
                v = v.re
            return ({(0, 0): Fr(v)} if v else {}, ONE)
        if t == "var":
            if e[1] == "x": return ({(1, 0): Fr(1)}, ONE)
            if e[1] == "y": return ({(0, 1): Fr(1)}, ONE)
            raise NameError(f"unknown variable {e[1]!r} (use x and y)")
        if t == "neg":
            n, d = ev(e[1]); return (pneg(n), d)
        if t == "pow":
            n, d = ev(e[1]); k = e[2]
            return (ppow(n, k), ppow(d, k)) if k >= 0 else (ppow(d, -k), ppow(n, -k))
        if t == "call":
            raise ValueError("multivar currently handles rational functions (no sin/cos/exp)")
        (n1, d1), (n2, d2) = ev(e[1]), ev(e[2])
        if t == "+": return (padd(pmul(n1, d2), pmul(n2, d1)), pmul(d1, d2))
        if t == "-": return (padd(pmul(n1, d2), pneg(pmul(n2, d1))), pmul(d1, d2))
        if t == "*": return (pmul(n1, n2), pmul(d1, d2))
        if t == "/":
            if not n2: raise ZeroDivisionError("division by the zero polynomial")
            return (pmul(n1, d2), pmul(d1, n2))
    return ev(ast)

# ------------------------------ univariate helpers (in c) ------------------------------
def _trim(p):
    p = list(p)
    while p and p[-1] == 0: p.pop()
    return p

def _ueval(p, c):
    acc = Fr(0)
    for v in reversed(p): acc = acc * c + v
    return acc

def _udiv(a, b):
    a = _trim(a); b = _trim(b)
    q = [Fr(0)] * max(len(a) - len(b) + 1, 1)
    while len(a) >= len(b) and a:
        k = len(a) - len(b); f = a[-1] / b[-1]; q[k] = f
        for i, bv in enumerate(b): a[i + k] -= f * bv
        a = _trim(a)
    return _trim(q), a

def _ugcd(a, b):
    a, b = _trim(a), _trim(b)
    while b:
        a, b = b, _udiv(a, b)[1]
    return [v / a[-1] for v in a] if a else a

def _uderiv(p): return _trim([i * p[i] for i in range(1, len(p))])

def _all_pole_curves(fN, fD):
    """True when every nonzero root of fD is simple and not a root of fN: then each real root is a
    pole curve (implicit function theorem), whatever field the coefficients live in."""
    f = _trim(fD)
    while f and f[0] == 0: f = f[1:]                      # c = 0 belongs to other curvatures
    if len(f) <= 1: return True
    return len(_ugcd(f, _uderiv(f))) <= 1 and len(_ugcd(f, _trim(fN))) <= 1


def _real_root_count(p):
    """distinct real NONZERO roots of p; coefficients may be rational or algebraic (exact signs)"""
    p = _trim(p)
    while p and p[0] == 0: p = p[1:]
    if len(p) <= 1: return 0
    g = _ugcd(p, _uderiv(p))
    if len(g) > 1: p = _udiv(p, g)[0]
    seq = [p, _uderiv(p)]
    while len(_trim(seq[-1])) > 1:
        r = _udiv(seq[-2], seq[-1])[1]
        if not r: break
        seq.append([-v for v in r])
    seq = [q for q in seq if q]
    def changes(vals):
        vals = [v for v in vals if v != 0]
        return sum(1 for a, b in zip(vals, vals[1:]) if (a > 0) != (b > 0))
    at_pinf = [sign_of(q[-1]) for q in seq]
    at_minf = [sign_of(q[-1]) * (-1) ** (len(q) - 1) for q in seq]
    return changes(at_minf) - changes(at_pinf)

def _is_rational_poly(p):
    return all(not isinstance(v, A) for v in p)

def _irrational_roots(p):
    """exact real algebraic numbers for the irrational real roots of a rational polynomial"""
    p = _trim(p)
    while p and p[0] == 0: p = p[1:]
    for r in _rational_roots(p):                    # divide out rational roots
        while _ueval(p, r) == 0: p = _udiv(p, [-r, Fr(1)])[0]
    out = []
    for root in isolate_real_roots(p):
        if root[0] == "interval":
            out.append(A.root(AlgField(p, root[1], root[2])))
    return out

def _sturm_count(p):
    """number of distinct real roots of p (rational coefficients)"""
    p = _trim(p)
    if len(p) <= 1: return 0
    p = _udiv(p, _ugcd(p, _uderiv(p)))[0]            # square-free part
    seq = [p, _uderiv(p)]
    while len(_trim(seq[-1])) > 1:
        r = _udiv(seq[-2], seq[-1])[1]
        if not r: break
        seq.append([-v for v in r])
    def changes(vals):
        vals = [v for v in vals if v != 0]
        return sum(1 for a, b in zip(vals, vals[1:]) if (a > 0) != (b > 0))
    lead = lambda q: q[-1]
    at_pinf = [lead(q) for q in seq if q]
    at_minf = [lead(q) * (-1) ** (len(q) - 1) for q in seq if q]
    return changes(at_minf) - changes(at_pinf)

def _divisors(n):
    n = abs(n); out = set()
    for i in range(1, isqrt(n) + 1):
        if n % i == 0: out |= {i, n // i}
    return out

def _simplest_between(lo, hi):
    """the rational with the smallest denominator in the open interval (lo, hi), lo < hi.
       Continued-fraction recursion: logarithmic in the size of the numbers."""
    from math import floor
    lo, hi = Fr(lo), Fr(hi)
    if lo < 0 < hi: return Fr(0)
    if hi <= 0: return -_simplest_between(-hi, -lo)
    n = floor(lo)
    if n + 1 < hi: return Fr(n + 1)
    a, b = lo - n, hi - n                         # 0 <= a < b <= 1
    if a == 0: return n + Fr(1, floor(1 / b) + 1)
    return n + 1 / _simplest_between(1 / b, 1 / a)

def _rational_roots(p, lo=None, hi=None):
    """distinct nonzero rational roots of a rational polynomial p.
       Small coefficients: rational root theorem. Large coefficients (deep zooms): isolate the
       real roots exactly and test the simplest rational in each shrinking interval. A rational
       root missed this way is still handled exactly later, as an algebraic number."""
    p = _trim(p)
    while p and p[0] == 0: p = p[1:]
    if len(p) <= 1: return []
    den = 1
    for v in p: den = den * v.denominator // gcd(den, v.denominator)
    ints = [int(v * den) for v in p]
    if max(abs(ints[0]), abs(ints[-1])) < 10 ** 12:
        roots = []
        for num in _divisors(ints[0]):
            for dd in _divisors(ints[-1]):
                for sgn in (1, -1):
                    r = Fr(sgn * num, dd)
                    if r not in roots and _ueval(p, r) == 0: roots.append(r)
        return roots
    from .algebraic import isolate_real_roots, count_roots
    roots = []
    for iv in isolate_real_roots(p, lo, hi):
        if iv[0] == "rational":
            if iv[1] != 0: roots.append(iv[1])
            continue
        l, h = iv[1], iv[2]
        for _ in range(60):
            c = _simplest_between(l, h)
            if _ueval(p, c) == 0:
                roots.append(c); break
            mid = (l + h) / 2
            if count_roots(p, l, mid) == 1: h = mid
            else: l = mid
    return roots

def _initial_form(P, a, sigma):
    """(weighted order, coefficient list in c) of P along x = sigma s^q, y = c s^p"""
    if not P: return None, []
    w = {k: k[0] + a * k[1] for k in P}
    m = min(w.values())
    form = {}
    for (i, j), v in P.items():
        if w[(i, j)] == m:
            form[j] = form.get(j, 0) + v * sigma ** i
    n = max(form) + 1
    return m, _trim([form.get(j, Fr(0)) for j in range(n)])

def _critical_as(P, Q=None):
    """curvatures a > 0 where two monomials tie in weight i + a*j.
       With Q given: ties between a monomial of P and a monomial of Q (where delta(a) can
       change sign). Without Q: ties inside P (where its initial form changes)."""
    ks1 = list(P); ks2 = list(Q) if Q is not None else ks1
    out = set()
    for (i1, j1) in ks1:
        for (i2, j2) in ks2:
            if j1 != j2:
                a = Fr(i2 - i1, j1 - j2)
                if a > 0: out.add(a)
    return out

# ------------------------------ path evaluation (exact, via 1-D engine) ------------------------------
def _along(P, sigma, c, p, q):
    """univariate coefficient list (in s) of P(sigma s^q, c s^p); c rational or algebraic"""
    c = c if isinstance(c, A) else Fr(c)
    out = {}
    for (i, j), v in P.items():
        e = q * i + p * j
        out[e] = out.get(e, 0) + v * sigma ** i * c ** j
    out = {k: v for k, v in out.items() if v != 0}
    if not out: return []
    return [out.get(e, Fr(0)) for e in range(max(out) + 1)]

def _state_1d(Nc, Dc):
    """state of the polynomial ratio N(s)/D(s) as s -> 0+ (exact over Q or Q(alpha)):
       ('value', v) | ('inf', sign) | ('undef',) | ('undef0',)"""
    Nc, Dc = _trim(Nc), _trim(Dc)
    if not Dc:
        return ("undef0",) if not Nc else ("undef",)
    if not Nc: return ("value", Fr(0))
    vn = next(i for i, v in enumerate(Nc) if v != 0)
    vd = next(i for i, v in enumerate(Dc) if v != 0)
    lead = Nc[vn] / Dc[vd]
    if vn > vd: return ("value", Fr(0))
    if vn == vd: return ("value", lead)
    return ("inf", sign_of(lead))

# ------------------------------ profile ------------------------------
def _ctext(c):
    return f"{float(c):.6g}(algebraic)" if isinstance(c, A) else str(c)

def _fmt_exp(e):
    return str(e) if e.denominator == 1 else f"({e})"

def _rf_text(n, d):
    def poly(p):
        terms = []
        for j, v in enumerate(p):
            if v == 0: continue
            mono = "" if j == 0 else ("c" if j == 1 else f"c^{j}")
            coef = str(v) if (v not in (1, -1) or not mono) else ("" if v == 1 else "-")
            terms.append(coef + mono if coef not in ("", "-") or mono else coef)
        return " + ".join(terms).replace("+ -", "- ") or "0"
    return f"({poly(n)}) / ({poly(d)})"


def _analyze(N, D, depth, sigmas, curve, out, to_orig=lambda u, v: (u, v), expo=lambda b: b):
    """curve(a, c_text) -> human-readable path description in ORIGINAL coordinates;
       to_orig(u, v) maps current coordinates to the original (x, y), so every reported
       path can be checked by direct evaluation of the original function."""
    # balance points: ties inside N, inside D, and BETWEEN N and D (where delta changes sign)
    crit = sorted(_critical_as(N) | _critical_as(D) | _critical_as(N, D))
    probes = set(crit)
    pts = [Fr(0)] + crit
    for lo_, hi_ in zip(pts, pts[1:]): probes.add((lo_ + hi_) / 2)
    probes.add(crit[-1] * 2 if crit else Fr(1))
    if crit: probes.add(crit[0] / 2)
    for a in sorted(probes):
        p, q = a.numerator, a.denominator
        for sigma in sigmas:
            mN, fN = _initial_form(N, a, sigma)
            mD, fD = _initial_form(D, a, sigma)
            if not D:
                out.append({"kind": "undef", "path": curve(a, "c", sigma), "a": a, "sigma": sigma}); continue
            if not N:
                out.append({"kind": "value", "value": Fr(0), "path": curve(a, "c", sigma), "a": a, "sigma": sigma}); continue
            delta = mN - mD
            entry = {"a": a, "sigma": sigma, "delta": delta, "critical": a in crit and depth == 0,
                     "path": curve(a, "c", sigma), "fN": fN, "fD": fD, "depth": depth,
                     "point": (lambda c, t, a=a, sigma=sigma: to_orig(sigma * t ** a.denominator, c * t ** a.numerator))}
            if delta > 0:
                entry.update(kind="value", value=Fr(0))
            elif delta < 0:
                entry.update(kind="inf")
            else:
                g = _ugcd(fN, fD)
                rn, rd = _udiv(fN, g)[0], _udiv(fD, g)[0]
                if len(rn) <= 1 and len(rd) <= 1:
                    entry.update(kind="value", value=rn[0] / rd[0] if rn else Fr(0))
                else:
                    entry.update(kind="depends", R=(fN, fD), text=_rf_text(rn, rd))
            out.append(entry)
            # roots of initial forms that need zooming in
            # Denominator roots ALWAYS need zooming (D may be far smaller along those curves);
            # numerator roots only when delta < 0 (N may then be small enough to stay finite).
            need, open_roots = [], 0
            for k, f in enumerate([fD] + ([fN] if delta < 0 else [])):
                if _is_rational_poly(f):
                    need += _rational_roots(f) + _irrational_roots(f)     # exact, incl. sqrt(2)-type
                else:
                    g = _trim(f)
                    while g and g[0] == 0: g = g[1:]
                    if len(g) == 2:                                      # linear: exact root in the same field
                        need.append(-g[0] / g[1]); continue
                    cnt = _real_root_count(f)                            # roots in a field tower
                    if cnt and k == 0 and _all_pole_curves(fN, fD):
                        here = curve(a, "c", sigma) + f"  (c = any of the {cnt} real roots of the denominator's initial form)"
                        out.append({"kind": "inf", "a": a, "sigma": sigma, "pole_curve": True, "signs": {1, -1},
                                    "path": here + ": the denominator vanishes on real curves here; implicit function theorem"})
                        if delta > 0:
                            out.append({"kind": "all", "a": a, "sigma": sigma, "path": here,
                                        "why": "f crosses through infinity here while it is near 0 on both sides"})
                    elif cnt and k == 1 and _all_pole_curves(fD, fN):
                        out.append({"kind": "all", "a": a, "sigma": sigma,
                                    "path": curve(a, "c", sigma) + f"  (c = any of the {cnt} real roots of the numerator's initial form)",
                                    "why": "f changes sign through 0 here while it is unbounded on both sides"})
                    else:
                        open_roots += cnt
            uniq = []
            for c0 in need:
                if not any(same_value(c0, u) is True for u in uniq): uniq.append(c0)
            need = uniq
            if open_roots > 0:
                out.append({"kind": "unresolved", "a": a, "sigma": sigma,
                            "path": curve(a, "c", sigma) + "  (c = root of an initial form with algebraic coefficients)"})
            for c0 in need:
                # ---- resolution stopping rule: a SIMPLE root that is not shared needs no zoom ----
                vD, vN = sign_of(_ueval(fD, c0)), sign_of(_ueval(fN, c0))
                simple_pole = vD == 0 and vN != 0 and sign_of(_ueval(_uderiv(fD), c0)) != 0
                simple_zero = (delta < 0 and vN == 0 and vD != 0
                               and sign_of(_ueval(_uderiv(fN), c0)) != 0)
                if simple_pole or simple_zero:
                    here = curve(a, _ctext(c0), sigma) + " + higher terms"
                    if simple_pole:
                        # implicit function theorem: D = 0 on a real curve here, N != 0 beside it
                        out.append({"kind": "inf", "a": a, "sigma": sigma, "pole_curve": True, "signs": {1, -1},
                                    "path": here + " (the denominator vanishes on a real curve here;"
                                                   " implicit function theorem)"})
                    if simple_pole and delta > 0 or simple_zero:
                        # intermediate value theorem across the curve: on each small arc crossing it
                        # f runs from about 0 (or -big) to +-infinity (or +big) -> every real value
                        out.append({"kind": "all", "a": a, "sigma": sigma, "path": here,
                                    "why": ("f crosses through infinity here while it is near 0 on both sides"
                                            if simple_pole else
                                            "f changes sign through 0 here while it is unbounded on both sides")})
                    continue
                if depth >= MAX_DEPTH or len(out) >= WORK_BUDGET:
                    why = "depth limit" if depth >= MAX_DEPTH else "work budget"
                    # A simple root of fD that is not a root of fN: by the implicit function theorem
                    # D = 0 along a real curve inside this sector, and N stays away from 0 there, so
                    # f is unbounded in every neighbourhood of the origin.  Proven without zooming.
                    if (sign_of(_ueval(_uderiv(fD), c0)) != 0 and sign_of(_ueval(fN, c0)) != 0):
                        out.append({"kind": "inf", "a": a, "sigma": sigma, "pole_curve": True, "signs": {1, -1},
                                    "path": curve(a, _ctext(c0), sigma)
                                    + " (+ higher terms: a curve where the denominator vanishes and the"
                                      " numerator does not; implicit function theorem)"})
                        continue
                    out.append({"kind": "unresolved", "a": a, "sigma": sigma,
                                "path": curve(a, _ctext(c0), sigma) + f" (+ higher terms; {why})"})
                    continue
                # exactly on the curve y = c0 x^a
                st = _state_1d(_along(N, sigma, c0, p, q), _along(D, sigma, c0, p, q))
                base = curve(a, _ctext(c0), sigma)
                rec = {"a": a, "sigma": sigma, "path": base + " exactly", "depth": depth, "fixed_c": c0,
                       "point": (lambda c, t, a=a, sigma=sigma, c0=c0: to_orig(sigma * t ** a.denominator, c0 * t ** a.numerator))}
                if st[0] == "value": rec.update(kind="value", value=st[1])
                elif st[0] == "inf": rec.update(kind="inf", signs={st[1]})
                elif st[0] == "undef0": rec.update(kind="undef", note="0/0 along the whole curve (common factor)")
                else: rec.update(kind="undef", note="denominator is zero along the whole curve")
                out.append(rec)
                # zoom: y = s^p (c0 + w), x = sigma s^q  -> new polynomials in (s, w)
                def zoom(P):
                    R = {}
                    for (i, j), v in P.items():
                        base_e = q * i + p * j
                        # (c0 + w)^j expansion
                        coef = Fr(1)
                        for k in range(j + 1):
                            ck = coef * c0 ** (j - k)
                            if ck:
                                key = (base_e, k)
                                R[key] = R.get(key, 0) + v * sigma ** i * ck
                            coef = coef * (j - k) / (k + 1)
                    return {k: v for k, v in R.items() if v != 0}
                N2, D2 = zoom(N), zoom(D)
                # a term w = c s^b at the next level is the term c |x|^expo2(b) of y
                def expo2(b, expo=expo, q=q, p=p): return expo((p + b) / q)
                def curve2(a2, ctext, s2, a=a, c0=c0, sigma=sigma, expo2=expo2):
                    return curve(a, _ctext(c0), sigma) + f" + {ctext}*|x|^{_fmt_exp(expo2(a2))}"
                def to_orig2(u, v, to_orig=to_orig, sigma=sigma, q=q, p=p, c0=c0):
                    return to_orig(sigma * u ** q, u ** p * (c0 + v))
                _analyze(N2, D2, depth + 1, (1,), curve2, out, to_orig2, expo2)


def profile(text):
    """Full fluid-state profile of a rational function of x, y at the origin."""
    N, D = parse_rational(text)
    return profile_polys(N, D)


def profile_polys(N, D, sides=(1, -1), y_axis=True):
    """profile from polynomial dicts.  sides=(1,) restricts to x > 0 (used by the
       three-variable engine on planes where the first coordinate is a positive parameter)."""
    out = []
    def curve(a, ctext, sigma):
        side = "x>0" if sigma > 0 else "x<0"
        return f"y = {ctext}*|x|^{_fmt_exp(a)}  ({side})"
    _analyze(N, D, 0, tuple(sides), curve, out)
    # the axes
    for sigma in sides:
        st = _state_1d(_along(N, sigma, 0, 1, 1), _along(D, sigma, 0, 1, 1))
        out.append(_axis_entry(st, f"y = 0 (x-axis, {'x>0' if sigma > 0 else 'x<0'})"))
    for tau in ((1, -1) if y_axis else ()):
        Ny = _along({(j, i): v for (i, j), v in N.items()}, tau, 0, 1, 1)
        Dy = _along({(j, i): v for (i, j), v in D.items()}, tau, 0, 1, 1)
        out.append(_axis_entry(_state_1d(Ny, Dy), f"x = 0 (y-axis, {'y>0' if tau > 0 else 'y<0'})"))
    from .cluster import cluster_set
    try:
        cl = cluster_set(out)
    except (ArithmeticError, NotImplementedError) as exc:
        cl = {"text": f"not computed ({exc})", "connected": None, "exact": False, "complete": False}
    return {"N": N, "D": D, "entries": out, "verdict": _verdict(out), "cluster": cl}

def _axis_entry(st, path):
    e = {"path": path, "a": None, "sigma": None}
    if st[0] == "value": e.update(kind="value", value=st[1])
    elif st[0] == "inf": e.update(kind="inf", signs={st[1]})
    else: e.update(kind="undef", note="function undefined along this axis")
    return e

def _verdict(entries):
    vals, ambiguous = [], False
    for e in entries:
        if e["kind"] != "value": continue
        v = e["value"]
        if isinstance(v, A) and len(v.c) <= 1: v = v.rational()
        cmp = [same_value(v, u) for u in vals]
        if any(c is True for c in cmp): continue
        if any(c is None for c in cmp): ambiguous = True; continue
        vals.append(v)
    dep = [e for e in entries if e["kind"] == "depends"]
    every = [e for e in entries if e["kind"] == "all"]
    inf = [e for e in entries if e["kind"] == "inf"]
    unres = [e for e in entries if e["kind"] == "unresolved"]
    if dep or every or len(vals) > 1 or (vals and inf):
        reasons = []
        if every: reasons.append(f"every real value is taken arbitrarily close to the origin beside "
                                 f"{every[0]['path']} ({every[0]['why']}; intermediate value theorem)")
        if dep: reasons.append(f"value depends on the path constant c along {dep[0]['path']}: {dep[0]['text']}")
        if len(vals) > 1: reasons.append(f"different exact values along different paths: {sorted(vals, key=float)}")
        if vals and inf: reasons.append(f"finite along some paths, infinite along {inf[0]['path']}")
        return ("NO LIMIT (proven)", reasons)
    if unres or ambiguous:
        return ("UNPROVEN", [f"open: {e['path']}" for e in unres] or ["two algebraic values could not be separated"])
    if inf and not vals:
        return ("DIVERGES to infinity along every path (proven)", [])
    (L,) = vals
    return (f"LIMIT = {L} (proven)", [])


def describe(text):
    r = profile(text)
    lines = [f"f(x, y) = {text}", f"  verdict: {r['verdict'][0]}"]
    for reason in r["verdict"][1]: lines.append(f"    because {reason}")
    cl = r["cluster"]
    note = "" if cl.get("exact") else "  (approximate: coefficients in a number field)"
    if not cl.get("complete"): note += "  (partial: open paths remain)"
    lines.append(f"  fluid state (all values f approaches): {cl['text']}{note}")
    if cl.get("connected") is False and cl.get("complete"):
        lines.append("  WARNING: not connected -- this contradicts the connectedness theorem, so a path family is missing")
    crit = sorted({e['a'] for e in r['entries'] if e.get('critical')})
    if crit: lines.append(f"  balance points (critical curvatures a): {', '.join(str(a) for a in crit)}")
    return "\n".join(lines)
