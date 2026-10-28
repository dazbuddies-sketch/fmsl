"""Independent check of fmsl.jelonek (python benchmarks/jelonek_crosscheck.py).

fmsl.jelonek finds the non-properness set S_F by Newton-Puiseux arcs at infinity. This check
uses no arcs. On random vertical lines a = r (and horizontal lines b = r) it

  1. finds every candidate point by elimination: a point (r, b0) of S_F, even a complex one,
     makes the number of solutions of F(x, y) = (r, b) drop as b -> b0, so b0 is a root of the
     leading coefficient (in x, and in y) of the resultants Res_y(f - r, g - b), Res_x(f - r, g - b),
     computed as polynomials in b (fraction-free Sylvester determinants);
     The lines are generic (r a random fraction such as 17/11), so they cross S_F transversally.
  2. decides each real candidate b0 directly: it solves F(x, y) = (r, b0 +- eps) exactly for
     eps = 1e-3 and 1e-6 and asks whether some real solution escapes (its size grows like a
     negative power of eps: a factor 1000 smaller eps, at least 2.5 times larger). Escaping solutions mean (r, b0) is in S_F;
  3. compares with fmsl.jelonek: the points it says are in S_F on that line must be exactly the
     candidates with escaping solutions.
It also checks every arc FMSL reports as a witness: along it, |p| -> infinity and F(p) tends to
the claimed point (exact rational arithmetic at s = 10^-12 and 10^-18).
"""
import os, random, sys, time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from fmsl.jelonek import jelonek_set, _poly, _peval
from fmsl.plane import bmul, badd, bscale, bdivexact, bclean, resultant_d, at_c
from fmsl.cluster import _real_roots, cmp

EXAMPLES = [  # (f, g, described S_F, checked by hand)
    ("x", "x*y", "the line a = 0"),
    ("x*y", "x*y^2-y", "the lines b = 0 and a = 1"),
    ("x^2*y", "x*y", "the lines a = 0 and b = 0"),
    ("x+y^2", "y", "empty (a polynomial automorphism)"),
    ("x^2*y^2+x*y", "y", "the half-line b = 0, a >= -1/4"),
    ("x*y", "x^2*y^2+x", "the parabola b = a^2"),
]


# ------------------------------------------------------------------ elimination (no arcs)
def _sylvester_det(A, B):
    """Res_y(A, B) for A, B lists (coefficients in y, low to high) of bivariate polys in (x, b)"""
    m, n = len(A) - 1, len(B) - 1
    if m < 0 or n < 0: return {}
    if m == 0 or n == 0:                               # Res = A0^n or B0^m
        base, e = (A[0], n) if m == 0 else (B[0], m)
        out = {(0, 0): Fr(1)}
        for _ in range(e): out = bmul(out, base)
        return out
    N = m + n
    M = []
    for i in range(n): M.append([A[m - (k - i)] if 0 <= k - i <= m else {} for k in range(N)])
    for i in range(m): M.append([B[n - (k - i)] if 0 <= k - i <= n else {} for k in range(N)])
    sign, prev = 1, {(0, 0): Fr(1)}
    for k in range(N - 1):
        if not M[k][k]:
            piv = next((r for r in range(k + 1, N) if M[r][k]), None)
            if piv is None: return {}
            M[k], M[piv] = M[piv], M[k]; sign = -sign
        for i in range(k + 1, N):
            for j in range(k + 1, N):
                M[i][j] = bdivexact(badd(bmul(M[k][k], M[i][j]), bscale(bmul(M[i][k], M[k][j]), -1)), prev) \
                    if (M[k][k] and M[i][j]) or (M[i][k] and M[k][j]) else {}
            M[i][k] = {}
        prev = M[k][k]
    return bscale(M[N - 1][N - 1], sign)


def _in_y(P, extra_b=None):
    """P(x, y) (+ a term in b) as a list over y of polys in (x, b)"""
    deg = max(j for (_, j) in P)
    out = [{} for _ in range(deg + 1)]
    for (i, j), v in P.items(): out[j][(i, 0)] = out[j].get((i, 0), 0) + v
    if extra_b: out[0] = badd(out[0], extra_b)
    return [bclean(c) for c in out]


def candidates(f, g, r):
    """real b0 on the line a = r where the solution count can drop (superset of S_F there)"""
    out = []
    for swap in (False, True):
        F = {((j, i) if swap else (i, j)): v for (i, j), v in f.items()}
        G = {((j, i) if swap else (i, j)): v for (i, j), v in g.items()}
        A = _in_y(badd(F, {(0, 0): -r}))
        B = _in_y(G, {(0, 1): Fr(-1)})
        R = _sylvester_det(A, B)                       # polynomial in (x, b)
        if not R: return None
        top = max(i for (i, _) in R)
        lc = [Fr(0)] * (max(k for (i, k) in R if i == top) + 1)
        for (i, k), v in R.items():
            if i == top: lc[k] += v
        while lc and lc[-1] == 0: lc.pop()
        if len(lc) > 1: out += _real_roots(lc)
    uniq = []
    for c in out:
        if not any(cmp(c, u) == 0 for u in uniq): uniq.append(c)
    return uniq


def _rat(v, tol=Fr(1, 10 ** 40)):
    if isinstance(v, Fr): return v
    from fmsl.cluster import _bounds
    lo, hi = _bounds(v)
    while hi - lo > tol: v.F.refine(); lo, hi = _bounds(v)
    return (lo + hi) / 2


def max_preimage(f, g, a, b):
    """largest |coordinate| of the real solutions of F = (a, b): x from the exact real roots of
       Res_y, rounded to 40 digits; y from the real roots of f(x0, y) = a; kept if g(x0, y) = b to
       within 1e-9 (relative)"""
    Fa, Gb = badd(f, {(0, 0): -a}), badd(g, {(0, 0): -b})
    Rx = resultant_d(Fa, Gb)                           # eliminate y: polynomial in x
    if not Rx or len(Rx) <= 1: return 0.0
    best = 0.0
    for x0 in _real_roots(Rx):
        xr = _rat(x0)
        fy = at_c(Fa, xr)
        while fy and fy[-1] == 0: fy.pop()
        ys = _real_roots(fy) if len(fy) > 1 else []
        for y0 in ys:
            yr = _rat(y0)
            terms = [v * xr ** i * yr ** j for (i, j), v in g.items()]
            gv = sum(terms)
            if abs(float(gv - b)) <= 1e-9 * (1 + sum(abs(float(t)) for t in terms)):
                best = max(best, abs(float(xr)), abs(float(yr)))
    return best


def escapes(f, g, a, b0):
    if not isinstance(b0, Fr):
        from fmsl.cluster import _bounds
        lo, hi = _bounds(b0)
        while hi - lo > Fr(1, 10 ** 12): b0.F.refine(); lo, hi = _bounds(b0)
        b0 = (lo + hi) / 2
    sizes = []
    for eps in (Fr(1, 10 ** 3), Fr(1, 10 ** 6)):
        sizes.append(max(max_preimage(f, g, a, b0 + eps), max_preimage(f, g, a, b0 - eps)))
    # escaping solutions grow like a negative power of eps (often eps^-1/3 or eps^-1/2): a factor
    # 1000 in eps then gives a factor >= 3 or so; bounded solutions barely move
    return sizes[1] > 2.5 * max(sizes[0], 1.0), sizes


def in_set(res, a, b):
    from fmsl.cluster import _real_roots as rr
    for L in res["lines"]:
        v, w = (a, b) if L["fixed"] == "a" else (b, a)
        if cmp(v, L["value"]) == 0 and any(cmp(lo, w) <= 0 <= cmp(hi, w) for lo, hi in L["ranges"]): return True
    for C in res["curves"]:
        ra = list(C["A"]); ra[0] -= a
        while ra and ra[-1] == 0: ra.pop()
        if len(ra) > 1 and any(cmp(_peval(C["B"], t), b) == 0 for t in rr(ra)): return True
    return any(cmp(p[0], a) == 0 and cmp(p[1], b) == 0 for p in res["points"])


ESCAPING = [0]


def line_check(f_text, g_text, r, horizontal=False):
    f, g = _poly(f_text), _poly(g_text)
    res = jelonek_set(f_text, g_text)
    if horizontal:                                       # b = r: exchange the roles of f and g
        f, g = g, f
    cands = candidates(f, g, r)
    if cands is None: return None
    bad = []
    for b0 in cands:
        esc, sizes = escapes(f, g, r, b0)
        a_, b_ = (b0, r) if horizontal else (r, b0)
        if esc: ESCAPING[0] += 1
        if esc != in_set(res, a_, b_): bad.append((float(a_), float(b_), esc, sizes))
    return len(cands), bad


def arc_check(f, g, piece_arc, chart, A, B, t):
    """along the reported arc with free coefficient t: |p| grows and F(p) -> (A(t), B(t))"""
    from math import lcm
    exps, coefs, E, Q = [], [], Fr(0), 1
    for c, b in piece_arc:
        E = E + Fr(b) / Q
        exps.append(E); coefs.append(t if c == "c" else _rat(c, Fr(1, 10 ** 60)))
        Q *= Fr(b).denominator
    L = lcm(*[e.denominator for e in exps], 1)
    target = (_peval(A, t), _peval(B, t))
    errs, sizes = [], []
    for m in (8, 16):
        u = Fr(1, 10 ** m); sv = u ** L                  # s = u^L, so every exponent is an integer
        var = chart[1] * (1 / sv)
        other = sum(cf * u ** int(e * L) for cf, e in zip(coefs, exps))
        x, y = (var, other) if chart[0] == "x" else (other, var)
        fv = sum(v * x ** i * y ** j for (i, j), v in f.items()); gv = sum(v * x ** i * y ** j for (i, j), v in g.items())
        errs.append(max(abs(float(fv - target[0])), abs(float(gv - target[1]))))
        sizes.append(max(abs(float(x)), abs(float(y))))
    return sizes[1] > sizes[0] > 1 and errs[1] < max(errs[0], 1e-300) * 0.5 + 1e-30


def arcs_ok(f_text, g_text):
    f, g = _poly(f_text), _poly(g_text)
    res = jelonek_set(f_text, g_text)
    pieces = [(pc["arc"], pc["chart"], pc["A"], pc["B"]) for L in res["lines"] for pc in L["arcs"]]
    pieces += [(C["arc"], C["chart"], C["A"], C["B"]) for C in res["curves"]]
    good = bad = 0
    for arc, chart, A, B in pieces:
        for t in (Fr(1, 3), Fr(-5, 7)):
            if arc_check(f, g, arc, chart, A, B, t): good += 1
            else: bad += 1; print("ARC FAILS", f_text, "|", g_text, arc, chart, t)
    return good, bad


def dominant(f_text, g_text):
    """the Jacobian determinant is not identically zero"""
    f, g = _poly(f_text), _poly(g_text)
    d = lambda P, k: {(i - (k == 0), j - (k == 1)): v * (i if k == 0 else j) for (i, j), v in P.items() if (i if k == 0 else j)}
    J = badd(bmul(d(f, 0), d(g, 1)), bscale(bmul(d(f, 1), d(g, 0)), -1))
    return bool(bclean(J))


def random_map(gen):
    def rp():
        T = {}
        for _ in range(gen.randint(2, 3)):
            i, j = gen.randint(0, 2), gen.randint(0, 2)
            if i + j: T[(i, j)] = T.get((i, j), 0) + gen.choice([-2, -1, 1, 2])
        return " + ".join(f"({v})*x^{i}*y^{j}" for (i, j), v in T.items() if v) or "x"
    return rp(), rp()


def main(n=40):
    t0 = time.time()
    gen = random.Random(3)
    maps = [(f, g) for f, g, _ in EXAMPLES] + [random_map(gen) for _ in range(n)]
    lines = bad_total = cand_total = skipped = 0
    arcs_good = arcs_bad = 0
    in_sf = [0]
    nondominant = 0
    for f, g in maps:
        if not dominant(f, g): nondominant += 1; continue      # F(R^2) is a curve: no line test applies
        a_ok, a_bad = arcs_ok(f, g); arcs_good += a_ok; arcs_bad += a_bad
        for k in range(2):
            r = Fr(gen.choice([-1, 1]) * gen.randint(1, 40), gen.choice([3, 7, 11, 13]))   # a generic line
            while r.denominator == 1: r += Fr(1, 7)
            for hz in (False, True):
                out = line_check(f, g, r, hz)
                if out is None: skipped += 1; continue
                nc, bad = out
                lines += 1; cand_total += nc
                if bad: bad_total += len(bad); print("DISAGREE", f, "|", g, "line", "b" if hz else "a", "=", r, bad)
    print(f"{len(maps)} maps ({len(EXAMPLES)} with known answers; {nondominant} with zero Jacobian skipped), {lines} random lines, {cand_total} candidate points "
          f"decided directly ({ESCAPING[0]} in S_F, {cand_total - ESCAPING[0]} not): {bad_total} disagree with "
          f"fmsl.jelonek; {skipped} lines skipped (degenerate); witness arcs: {arcs_good} confirmed, "
          f"{arcs_bad} failed; {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 40)
