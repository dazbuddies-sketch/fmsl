"""The non-properness (Jelonek) set of a polynomial map F: R^n -> R^n, n >= 2 (written for n = 3).

S_F = { lim F(p_k) : |p_k| -> oo, F(p_k) converges }. F is proper exactly over the complement of
S_F. Jelonek's real Jacobian conjecture says that a polynomial local diffeomorphism of R^n whose
S_F has codimension >= 2 is a bijection (so S_F is empty). Jelonek proved it for codimension
>= 3 and for n = 2; Pinchuk's map shows codimension 1 is not enough. Codimension 2 is open
for every n >= 3. A local diffeomorphism of R^3 whose S_F is a nonempty curve would be a
counterexample.

This module gives two bounds, and says which is which.

Inner bound (proven points, each with an arc):
  Arcs to infinity x_i = c_i s^a_i (a_i integers, some negative, s -> 0+), with one level of
  zooming x_i = s^a_i (c0_i + y_i) at rational points c0 where the leading forms that blow up
  vanish. When every component of F has order >= 0 along an arc, F tends to the explicit point
  L(c), the order-0 parts evaluated at c. Each such point is in S_F, proven by its arc.

Outer bound (a set containing S_F, exact):
  For each coordinate x_i, an exact polynomial relation R(x_i, F(x)) = 0 is found by linear
  algebra (modulo a large prime, then rational reconstruction and an exact check over Q). Write
  R = a(u) t^D + lower powers of t. Along a sequence with x_i -> oo and F -> u, dividing by x_i^D
  gives a(u) = 0. Every sequence to infinity has some coordinate unbounded, so
        S_F  is contained in  V(a_1) u ... u V(a_n)   (real points).
  Using a second set of coordinates (a random basis of linear forms) and intersecting gives a
  smaller set that still contains S_F.

Dimension of the outer bound: a random line meets a real surface with positive probability but
a curve with probability 0. Each line is checked exactly (Sturm's theorem on the restriction).
If no line among many has a real point of the outer bound, that is strong evidence that S_F has
dimension <= 1, but it is evidence, not proof, and the result says so.

    from fmsl.jelonek_n import analyse
    r = analyse(["x", "x*y", "z"], ["x", "y", "z"])
    r["inner_points"][:3], r["outer"]["surface_lines"], r["jacobian"]
"""
import random
from fractions import Fraction as Fr
from itertools import product
from math import gcd

from .witness import parse_n, _add, _mul, _ev, _initial, _roots, _utrim, _umul, _uadd, _ugcd, _ueval

__all__ = ["analyse", "relation", "inner_points", "outer_bound", "jacobian"]

PRIME = (1 << 61) - 1


# ------------------------------------------------------------------ parsing and basic algebra
def _polys(texts, names):
    out = []
    for t in texts:
        P, Q = parse_n(t, names)
        if len(Q) != 1 or any(k != (0,) * len(names) for k in Q): raise ValueError("polynomials only")
        q = next(iter(Q.values()))
        out.append({k: v / q for k, v in P.items()})
    return out


def _deriv(P, j):
    out = {}
    for k, v in P.items():
        if k[j]:
            kk = tuple(e - (i == j) for i, e in enumerate(k)); out[kk] = out.get(kk, 0) + v * k[j]
    return {k: v for k, v in out.items() if v != 0}


def jacobian(F, n):
    """det of the Jacobian matrix, as a polynomial dict"""
    M = [[_deriv(F[i], j) for j in range(n)] for i in range(n)]
    def det(rows, cols):
        if len(rows) == 1: return M[rows[0]][cols[0]]
        out = {}
        for k, c in enumerate(cols):
            sub = det(rows[1:], cols[:k] + cols[k + 1:])
            term = _mul(M[rows[0]][c], sub)
            if k % 2: term = {a: -b for a, b in term.items()}
            out = _add(out, term)
        return out
    return det(list(range(n)), list(range(n)))


# ------------------------------------------------------------------ outer bound: relations
def _modinv(a): return pow(a % PRIME, PRIME - 2, PRIME)


def _ratrec(a):
    """rational reconstruction of a modulo PRIME, or None"""
    a %= PRIME
    r0, r1, s0, s1 = PRIME, a, 0, 1
    bound = int((PRIME // 2) ** 0.5)
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound: return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)


def relation(F, n, ell, max_total=12, max_cols=2500):
    """a nonzero R(t, u) with R(ell(x), F(x)) = 0, as {(k, beta): coeff}, exact; or None.
       ell: a linear form, as a dict like the polys. Columns t^k u^beta are added by increasing
       weighted degree k + deg(F) |beta| until a dependency mod PRIME appears; it is then
       reconstructed over Q and checked exactly."""
    dF = max(sum(k) for P in F for k in P)
    mods = lambda P: {k: (v.numerator * _modinv(v.denominator)) % PRIME for k, v in P.items()}
    def mmul(A, B):
        R = {}
        for a, u in A.items():
            for b, v in B.items():
                k = tuple(i + j for i, j in zip(a, b)); R[k] = (R.get(k, 0) + u * v) % PRIME
        return {k: v for k, v in R.items() if v}
    Fm, Lm = [mods(P) for P in F], mods(ell)
    one = {(0,) * n: 1}
    tpow, upow = {0: one}, {(0,) * n: one}
    def tp(k):
        if k not in tpow: tpow[k] = mmul(tp(k - 1), Lm)
        return tpow[k]
    def up(beta):
        if beta not in upow:
            i = next(j for j in range(n) if beta[j]); b2 = tuple(x - (j == i) for j, x in enumerate(beta))
            upow[beta] = mmul(up(b2), Fm[i])
        return upow[beta]
    basis = []                                       # (pivot monomial, vector, combination)
    cols = []
    for W in range(1, max_total * dF + 1):
        for k in range(0, W + 1):
            rest = W - k
            if rest % dF: continue
            m = rest // dF
            for beta in product(range(m + 1), repeat=n):
                if sum(beta) != m: continue
                cols.append((k, beta))
                if len(cols) > max_cols: return None
                v = mmul(tp(k), up(beta)); comb = {len(cols) - 1: 1}
                for piv, bv, bc in basis:
                    c = v.get(piv, 0)
                    if c:
                        for key, x in bv.items():
                            v[key] = (v.get(key, 0) - c * x) % PRIME
                            if not v[key]: del v[key]
                        for key, x in bc.items():
                            comb[key] = (comb.get(key, 0) - c * x) % PRIME
                            if not comb[key]: del comb[key]
                if v:
                    piv = min(v); inv = _modinv(v[piv])
                    basis.append((piv, {a: b * inv % PRIME for a, b in v.items()},
                                  {a: b * inv % PRIME for a, b in comb.items()}))
                    continue
                if not any(cols[j][0] > 0 for j in comb): continue      # a relation among F alone
                lead = max(comb, key=lambda j: (cols[j][0], j)); inv = _modinv(comb[lead])
                R = {}
                for j, x in comb.items():
                    q = _ratrec(x * inv)
                    if q is None: break
                    R[cols[j]] = q
                else:
                    if _check(R, F, ell, n): return R
    return None


def _check(R, F, ell, n):
    out, cacheF, cacheL = {}, {}, {0: {(0,) * n: Fr(1)}}
    def lp(k):
        if k not in cacheL: cacheL[k] = _mul(lp(k - 1), ell)
        return cacheL[k]
    def fp(beta):
        if beta not in cacheF:
            if not any(beta): cacheF[beta] = {(0,) * n: Fr(1)}
            else:
                i = next(j for j in range(n) if beta[j]); b2 = tuple(x - (j == i) for j, x in enumerate(beta))
                cacheF[beta] = _mul(fp(b2), F[i])
        return cacheF[beta]
    for (k, beta), c in R.items():
        out = _add(out, {a: b * c for a, b in _mul(lp(k), fp(beta)).items()})
    return not out


def _leading(R, n):
    D = max(k for k, _ in R)
    return {beta: c for (k, beta), c in R.items() if k == D}, D


def outer_bound(F, n, bases=2, seed=0):
    """the leading coefficients a(u) (polynomials in u = F-values) for each linear form of each
       basis; S_F lies in the intersection over bases of the union over that basis's forms"""
    gen = random.Random(seed)
    out = []
    for b in range(bases):
        if b == 0:
            forms = [{tuple(int(j == i) for j in range(n)): Fr(1)} for i in range(n)]
        else:
            while True:
                M = [[gen.randint(-3, 3) for _ in range(n)] for _ in range(n)]
                if _det_int(M): break
            forms = [{tuple(int(j == i) for j in range(n)): Fr(M[r][i]) for i in range(n) if M[r][i]} for r in range(n)]
        leads = []
        for ell in forms:
            R = relation(F, n, ell)
            if R is None: return None
            a, D = _leading(R, n)
            leads.append({"form": ell, "degree": D, "lead": a})
        out.append(leads)
    return out


def _det_int(M):
    n = len(M)
    if n == 1: return M[0][0]
    return sum((-1) ** c * M[0][c] * _det_int([r[:c] + r[c + 1:] for r in M[1:]]) for c in range(n))


def _restrict(P, p0, v):
    """P(p0 + t v) as a univariate polynomial (coefficient list)"""
    lin = [[Fr(a), Fr(b)] for a, b in zip(p0, v)]
    s = []
    for k, c in P.items():
        m = [c]
        for L, e in zip(lin, k):
            for _ in range(e): m = _umul(m, _utrim(L))
        s = _uadd(s, m)
    return s


def surface_test(outer, n, lines=200, seed=1, box=6):
    """count random lines that meet the outer bound in a real point (exact per line).
       A real surface inside the outer bound is met by a positive share of lines; a curve almost
       never. Returns (lines meeting it, lines tried)"""
    gen = random.Random(seed)
    hits = 0
    for _ in range(lines):
        p0 = [Fr(gen.randint(-box * 89, box * 89), 89) for _ in range(n)]
        v = [Fr(gen.randint(-97, 97), gen.randint(1, 97)) for _ in range(n)]
        if not any(v): continue
        common = None
        for leads in outer:
            prod = [Fr(1)]
            for L in leads: prod = _umul(prod, _restrict(L["lead"], p0, v))
            prod = _utrim(prod)
            common = prod if common is None else _ugcd(common, prod)
            if not common: common = prod
        if common and len(common) > 1 and _roots(common, Fr(-10 ** 6), Fr(10 ** 6)): hits += 1
        elif not common: hits += 1                   # identically zero on the line
    return hits, lines


# ------------------------------------------------------------------ inner bound: arcs
def _order(P, w):
    return min(sum(a * b for a, b in zip(k, w)) for k in P) if P else None


def _init(P, w):
    m = _order(P, w)
    return m, {k: v for k, v in P.items() if sum(a * b for a, b in zip(k, w)) == m}


def _sub(P, w, c0, j, n):
    """x_j = c0_j s^w_j, x_i = s^w_i (c0_i + y_i): new variables (s, y_i for i != j), Laurent in s"""
    others = [i for i in range(n) if i != j]
    def var(k): return {tuple(int(i == k) for i in range(n)): Fr(1)}
    X = []
    for i in range(n):
        sw = {tuple([w[i]] + [0] * (n - 1)): Fr(1)}
        if i == j: X.append({k: v * c0[i] for k, v in sw.items()})
        else: X.append(_mul(sw, _add({(0,) * n: c0[i]} if c0[i] else {}, var(1 + others.index(i)))))
    cache = {}
    def pw(i, e):
        if (i, e) not in cache: cache[(i, e)] = {(0,) * n: Fr(1)} if e == 0 else _mul(pw(i, e - 1), X[i])
        return cache[(i, e)]
    out = {}
    for k, v in P.items():
        m = {(0,) * n: v}
        for i, e in enumerate(k):
            if e: m = _mul(m, pw(i, e))
        out = _add(out, m)
    return out


def _zeros(polys, n, positive, gen, tries=60):
    grid = [Fr(k) for k in (-3, -2, -1, 1, 2, 3)] + [Fr(1, 2), Fr(-1, 2)]
    pick = lambda i: abs(gen.choice(grid)) if i in positive else gen.choice(grid)
    found = []
    for _ in range(tries):
        c = tuple(pick(i) for i in range(n))
        if all(_ev(P, c) == 0 for P in polys): found.append(c)
    P0 = polys[0]
    for k in range(n):                               # fix the others, solve for c_k exactly
        if max(m[k] for m in P0) < 1: continue
        for _ in range(8):
            c = [pick(i) for i in range(n)]
            uni = {}
            for m, v in P0.items():
                rest = v
                for i, e in enumerate(m):
                    if i != k and e: rest *= c[i] ** e
                uni[m[k]] = uni.get(m[k], 0) + rest
            coeffs = _utrim([uni.get(d, Fr(0)) for d in range(max(uni) + 1)])
            for r in _rational_roots(coeffs):
                if r == 0 or (k in positive and r < 0): continue
                c2 = list(c); c2[k] = r
                if all(_ev(P, c2) == 0 for P in polys): found.append(tuple(c2))
    return list(dict.fromkeys(found))


def _rational_roots(p):
    """the rational roots of a univariate polynomial with rational coefficients"""
    from math import lcm, isqrt
    p = _utrim(p)
    if len(p) <= 1: return []
    L = 1
    for c in p: L = lcm(L, c.denominator)
    a = [int(c * L) for c in p]
    while a and a[0] == 0: a = a[1:]
    if len(a) <= 1: return [Fr(0)] if len(p) > len(a) else []
    def divs(m):
        m = abs(m); out = set()
        if m > 10 ** 12: return {1}
        for d in range(1, isqrt(m) + 1):
            if m % d == 0: out |= {d, m // d}
        return out
    out = set()
    for P in divs(a[0]):
        for Q in divs(a[-1]):
            for r in (Fr(P, Q), Fr(-P, Q)):
                if _ueval(p, r) == 0: out.add(r)
    if len(p) > len(a): out.add(Fr(0))
    return sorted(out)


def inner_points(F, n, K=3, depth=1, samples=4, seed=0, max_zooms=40):
    """proven points of S_F with their arcs: [(point, arc)], where arc = (weights, c, zoom chain)"""
    gen = random.Random(seed)
    out = []
    todo = [(F, None, [])]
    zooms_done = 0
    for level in range(depth + 1):
        nxt = []
        for G, pos, chain in todo:
            positive = () if pos is None else (0,)
            rng = range(-K, K + 1) if pos is None else range(1, K + 1)
            for w in product(rng, repeat=n):
                if pos is None and min(w) >= 0: continue          # the arc must go to infinity
                if pos is not None and min(w) <= 0: continue
                g = 0
                for x in w: g = gcd(g, abs(x))
                if g != 1: continue
                ords = [_init(P, w) if P else (10 ** 9, {}) for P in G]
                blow = [I for m, I in ords if m < 0]
                if not blow:
                    for _ in range(samples):
                        c = tuple((abs(x) if i in positive else x) for i, x in
                                  enumerate(gen.choice([Fr(k) for k in (-3, -2, -1, 1, 2, 3)] + [Fr(1, 2)])
                                            for _ in range(n)))
                        if pos is None and not any(c[i] for i in range(n) if w[i] < 0): continue
                        pt = tuple((_ev(I, c) if m == 0 else Fr(0)) for m, I in ords)
                        out.append((pt, (w, c, chain)))
                    continue
                if level == depth or zooms_done >= max_zooms: continue
                for c0 in _zeros(blow, n, positive, gen)[:2]:
                    if pos is None and not any(c0[i] for i in range(n) if w[i] < 0): continue
                    j = next(i for i in range(n) if c0[i] != 0 and (pos is not None or w[i] < 0))
                    nxt.append(([_sub(P, w, c0, j, n) for P in G], 0, chain + [(w, c0, j)]))
                    zooms_done += 1
        todo = nxt
    return out


def arc_point(arc, s):
    """the point at parameter s > 0 on an inner-bound arc (exact)"""
    w, c, chain = arc
    x = [ci * Fr(s) ** wi for ci, wi in zip(c, w)]
    for w0, c0, j in reversed(chain):
        n = len(c0); others = [i for i in range(n) if i != j]
        t, ys = x[0], x[1:]
        x = [c0[i] * t ** w0[i] if i == j else t ** w0[i] * (c0[i] + ys[others.index(i)]) for i in range(n)]
    return tuple(x)


# ------------------------------------------------------------------ everything together
def analyse(texts, names, lines=200, K=3):
    F = _polys(texts, names)
    n = len(names)
    J = jacobian(F, n)
    pts = inner_points(F, n, K=K)
    uniq = sorted(set(p for p, _ in pts), key=lambda p: tuple(float(x) for x in p))
    outer = outer_bound(F, n)
    hits = surface_test(outer, n, lines) if outer else None
    return {"jacobian": J, "inner_points": uniq, "inner_arcs": pts, "outer": outer,
            "surface_lines": hits}
