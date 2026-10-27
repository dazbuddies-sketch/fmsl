"""
fmsl.bifurcation -- bifurcations of iterated maps as balance points.

For a one-parameter family of polynomial maps x -> f(x, r):
  fixed points:      F1(x, r) = f(x) - x = 0
  period-n points:   Qn(x, r) = (f^n(x) - x) / (f(x) - x)    (exact polynomial division)
At a bifurcation point (x0, r0), numerator and denominator of Qn vanish together: a 0/0
transition point. The new orbits are born along the real branches of Qn = 0 (or F1 = 0)
through (x0, r0). With u = x - x0 and e = r - r0, every branch has the Puiseux form
      e = c |u|^a + ...      (x-side sigma = +-1)
and is found EXACTLY from the Newton polygon: a is a balance curvature (an edge slope),
c a real root of the edge's initial form (rational or real algebraic).
So the size of the new orbit near the bifurcation is |u| ~ (e/c)^(1/a): the exponent and
the coefficient are both exact.
"""
from fractions import Fraction as Fr
from math import comb

from .multivar import _critical_as, _initial_form, _rational_roots, _irrational_roots, _trim
from .algebraic import A

# bivariate polynomials in (x, r): {(i, j): coeff}  meaning  coeff * x^i * r^j
def _add(P, Q):
    R = dict(P)
    for k, v in Q.items():
        R[k] = R.get(k, 0) + v
        if R[k] == 0: del R[k]
    return R

def _scale(P, s): return {k: v * s for k, v in P.items()}

def _mul(P, Q):
    R = {}
    for (a, b), u in P.items():
        for (c, d), v in Q.items():
            k = (a + c, b + d); R[k] = R.get(k, 0) + u * v
    return {k: v for k, v in R.items() if v != 0}

def compose(f, g):
    """f(g(x, r), r)"""
    out, gp, cache = {}, {(0, 0): Fr(1)}, {0: {(0, 0): Fr(1)}}
    maxi = max(i for i, _ in f)
    for i in range(1, maxi + 1):
        gp = _mul(gp, g); cache[i] = gp
    for (i, j), v in f.items():
        out = _add(out, _mul(cache[i], {(0, j): v}))
    return out

def iterate(f, n):
    g = {(1, 0): Fr(1)}
    for _ in range(n): g = compose(f, g)
    return g

def _lead(P):
    """leading term in lex order (x-degree first, then r-degree)"""
    k = max(P, key=lambda t: (t[0], t[1]))
    return k, P[k]

def divide_exact(P, Q):
    """exact division P / Q in Q[x, r] (lex order); raises if not exact"""
    P, quo = dict(P), {}
    (qi, qj), qc = _lead(Q)
    while P:
        (pi, pj), pc = _lead(P)
        if pi < qi or pj < qj: raise ArithmeticError("not an exact division")
        t = {(pi - qi, pj - qj): pc / qc}
        quo = _add(quo, t)
        P = _add(P, _scale(_mul(t, Q), -1))
    return quo

def period_polynomial(f, n):
    """Q_n = (f^n(x) - x) / (f(x) - x)"""
    x = {(1, 0): Fr(1)}
    return divide_exact(_add(iterate(f, n), _scale(x, -1)), _add(f, _scale(x, -1)))

def fixed_polynomial(f):
    return _add(f, {(1, 0): Fr(-1)})

def _exact(v): return v if isinstance(v, A) else Fr(v)

def shift(P, x0, r0):
    """P(u + x0, e + r0) as a polynomial in (u, e); x0, r0 rational or real algebraic"""
    x0, r0 = _exact(x0), _exact(r0)
    xp, rp = {0: Fr(1)}, {0: Fr(1)}
    out = {}
    for (i, j), v in P.items():
        for a in range(i + 1):
            for b in range(j + 1):
                if i - a not in xp: xp[i - a] = x0 ** (i - a)
                if j - b not in rp: rp[j - b] = r0 ** (j - b)
                c = v * comb(i, a) * comb(j, b) * xp[i - a] * rp[j - b]
                out[(a, b)] = out.get((a, b), 0) + c
    return {k: v for k, v in out.items() if v != 0}

def branches(P):
    """Real branches of P(u, e) = 0 through (0, 0): list of dicts
       {'a': curvature, 'sigma': side of u, 'c': coefficient (exact)}  meaning  e = c |u|^a,
       plus axis branches ('u = 0' / 'e = 0') when P vanishes identically on an axis."""
    out = []
    if all(i > 0 for (i, j) in P): out.append({"axis": "u = 0"})
    if all(j > 0 for (i, j) in P): out.append({"axis": "e = 0"})
    for a in sorted(_critical_as(P)):
        for sigma in (1, -1):
            _, form = _initial_form(P, a, sigma)
            if all(not isinstance(v, A) for v in form):
                roots = _rational_roots(form) + _irrational_roots(form)
            else:
                f = _trim(form)
                while f and f[0] == 0: f = f[1:]                 # drop the c = 0 root
                if len(f) == 2: roots = [-f[0] / f[1]]           # linear: exact root in Q(alpha)
                elif len(f) <= 1: roots = []
                else: raise NotImplementedError("initial form of degree > 1 over an algebraic field")
            for c in roots:
                out.append({"a": a, "sigma": sigma, "c": c})
    return out

def describe_branches(P):
    lines = []
    for b in branches(P):
        if "axis" in b:
            lines.append(f"branch {b['axis']}")
            continue
        a, c, sg = b["a"], b["c"], b["sigma"]
        side = "u > 0" if sg > 0 else "u < 0"
        cval = f"{float(c):.6g} (algebraic)" if isinstance(c, A) else str(c)
        lines.append(f"branch e = {cval} * |u|^{a}  on {side}   =>  |u| = (e / {cval})^(1/{a})")
    return lines


def _dx(P): return {(i - 1, j): v * i for (i, j), v in P.items() if i > 0}


def branch_stability(f, n, x0, r0, branch):
    """Is the period-n orbit born along `branch` attracting?  The multiplier
       lambda = d(f^n)/dx along the branch is a Puiseux series lambda0 + k s^m + ... ;
       at the bifurcation |lambda0| = 1, and the sign of k decides (exactly).
       Returns (lambda0, k, m, 'stable' | 'unstable')."""
    from .multivar import _along
    from .algebraic import sign_of
    M = shift(_dx(iterate(f, n)), x0, r0)
    lam0 = M.get((0, 0), 0)
    a, sg, c = branch["a"], branch["sigma"], branch["c"]
    p, q = a.numerator, a.denominator
    ser = _along({k: v for k, v in M.items() if k != (0, 0)}, sg, c, p, q)   # u = sg s^q, e = c s^p
    m = next((i for i, v in enumerate(ser) if v != 0), None)
    if m is None: return lam0, 0, None, "neutral to all orders"
    k = ser[m]
    if lam0 == 1: verdict = "stable" if sign_of(k) < 0 else "unstable"
    elif lam0 == -1: verdict = "stable" if sign_of(k) > 0 else "unstable"
    else: verdict = "stable" if abs(float(lam0)) < 1 else "unstable"
    return lam0, k, m, verdict
