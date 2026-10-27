"""Independent cross-check of the exact fluid-state arcs (fmsl.cluster).

The fluid state (cluster set) of f = N/D at the origin is the limit, as r -> 0, of the image of
f on ANY small loop around the origin.  Here the loop is the square |x|, |y| = r.  On each side
f is a ratio of polynomials in one variable v in [-1, 1]; it is monotone between consecutive
breakpoints (zeros of the derivative along the side, zeros of D, corners), so the image of the
loop is an exact union of intervals.  No Newton polygon, no Puiseux series, no zooming: only
exact rational arithmetic and Sturm root isolation.

For r = 10^-4, 10^-8, 10^-16, 10^-32 the script measures how much of the circle R u {oo}
(in the angle 2*atan(value)) differs between the loop image and the engine's exact arc.
The difference must shrink to zero.  Run:  python benchmarks/loop_crosscheck.py
"""
import math, random, sys, time
from fractions import Fraction as Fr

sys.path.insert(0, "src")
from fmsl.multivar import parse_rational, profile
from fmsl.algebraic import isolate_real_roots, squarefree, count_roots, ptrim, pmul, psub, AlgField, A, _ieval, pgcd, pdivmod
from fmsl.cluster import _Inf

GRID = 7200


def fl(x):
    """float, saturating at +-1e300 (only the angle 2*atan(value) is used)"""
    try: return float(x)
    except OverflowError: return 1e300 if x > 0 else -1e300


def dx(P): return {(i - 1, j): v * i for (i, j), v in P.items() if i > 0}
def dy(P): return {(i, j - 1): v * j for (i, j), v in P.items() if j > 0}

def restrict(P, side, r):
    """coefficients in v of P on a side of the square"""
    kind, s = side
    out = {}
    for (i, j), v in P.items():
        if kind == "x":   # x = s r, y = r v
            c, e = v * (s * r) ** i * r ** j, j
        else:             # y = s r, x = r v
            c, e = v * (s * r) ** j * r ** i, i
        out[e] = out.get(e, 0) + c
    return ptrim([out.get(e, Fr(0)) for e in range(max(out) + 1)]) if out else []

def ev(p, v):
    acc = Fr(0)
    for c in reversed(p): acc = acc * v + c
    return acc


class Root:
    def __init__(s, p, iv):
        s.p = p
        if iv[0] == "rational": s.lo = s.hi = iv[1]
        else: s.lo, s.hi = iv[1], iv[2]
    def refine(s):
        if s.lo == s.hi: return
        m = (s.lo + s.hi) / 2
        if ev(s.p, m) == 0: s.lo = s.hi = m
        elif count_roots(s.p, s.lo, m) == 1: s.hi = m
        else: s.lo = m
    def value(s, f, Nv, Dv):
        """f at the root, exactly: the root is an algebraic number, Nv/Dv is evaluated in its
           field and enclosed by interval arithmetic to 1e-12 (no sampling: peaks can be
           narrower than any fixed grid)"""
        if s.lo == s.hi: return fl(f(s.lo))
        F = AlgField(s.p, s.lo, s.hi); u = A.root(F)
        num = ev(Nv, u); den = ev(Dv, u)
        val = num / den
        if not isinstance(val, A): return fl(val)
        g = F.reduce(val.c)
        while True:
            a, b = _ieval(g, F.lo, F.hi)
            if b - a <= Fr(1, 10 ** 12) * (1 + abs(a)): return fl((a + b) / 2)
            for _ in range(4): F.refine()


def roots_in(p, lo=Fr(-1), hi=Fr(1)):
    p = squarefree(p)
    if len(p) <= 1: return []
    out = []
    for iv in isolate_real_roots(p):
        R = Root(p, iv)
        # a root exactly at a corner (v = -1 or 1) is the corner itself: the corner handles it
        if any(R.lo < e <= R.hi and ev(R.p, e) == 0 for e in (lo, hi)): continue
        while True:                                                   # root lies in (R.lo, R.hi]
            if R.lo != R.hi and ev(R.p, R.hi) == 0: R.lo = R.hi       # root exactly at the right end
            if R.hi <= lo or R.lo >= hi: break                        # outside the side
            if R.lo == R.hi:
                if lo < R.lo < hi: out.append(R)                      # (a root at a corner: corner handles it)
                break
            if lo < R.lo and R.hi < hi: out.append(R); break          # strictly inside
            R.refine()
    return out


def side_image(Nv, Dv, Gv):
    """closed extended intervals covered by f = Nv/Dv on v in [-1, 1]; Gv = numerator of f'."""
    f = lambda v: ev(Nv, v) / ev(Dv, v)
    bps = [(Fr(-1), "corner")]
    G2 = squarefree(Gv)                          # critical points that are also poles: keep as poles
    while True:
        h = pgcd(G2, Dv)
        if len(h) <= 1: break
        G2 = pdivmod(G2, h)[0]
    rs = [(R, "pole") for R in roots_in(Dv)] + [(R, "crit") for R in roots_in(G2)]
    # separate and sort
    while True:                                  # refine overlapping neighbours until disjoint
        rs.sort(key=lambda t: t[0].lo)
        bad = [(a, b) for a, b in zip(rs, rs[1:])
               if not (a[0].hi < b[0].lo or a[0].lo == a[0].hi == b[0].lo == b[0].hi)]
        if not bad: break
        for a, b in bad:
            for _ in range(16): a[0].refine(); b[0].refine()
    merged = []
    for R, k in rs:
        if merged and merged[-1][0].lo == merged[-1][0].hi == R.lo == R.hi:
            if k == "pole": merged[-1] = (R, "pole")
            continue
        merged.append((R, k))
    bps += merged + [(Fr(1), "corner")]
    out = []
    for (A_, ka), (B_, kb) in zip(bps, bps[1:]):
        a_hi = A_ if ka == "corner" else A_.hi
        b_lo = B_ if kb == "corner" else B_.lo
        if a_hi >= b_lo: continue
        m1, m2 = a_hi + (b_lo - a_hi) / 3, a_hi + 2 * (b_lo - a_hi) / 3
        inc = f(m2) > f(m1)
        def end(X, k, left):
            if k == "corner":
                if ev(Dv, X) == 0: k = "pole"
                else: return fl(f(X))
            if k == "pole": return (-math.inf if inc else math.inf) if left else (math.inf if inc else -math.inf)
            return X.value(f, Nv, Dv)
        e1, e2 = end(A_, ka, True), end(B_, kb, False)
        out.append((min(e1, e2), max(e1, e2)))
    return out


def loop_image(N, D, r):
    Nx, Ny, Dx, Dy = dx(N), dy(N), dx(D), dy(D)
    ivs = []
    for side in (("x", 1), ("x", -1), ("y", 1), ("y", -1)):
        Nv, Dv = restrict(N, side, r), restrict(D, side, r)
        if not Dv: return None
        if Nv:                                   # cancel common factors (f is the reduced ratio)
            g = pgcd(Nv, Dv)
            if len(g) > 1: Nv, Dv = pdivmod(Nv, g)[0], pdivmod(Dv, g)[0]
        # derivative of Nv/Dv along v: Nv' Dv - Nv Dv'
        dN = ptrim([k * Nv[k] for k in range(1, len(Nv))]); dD = ptrim([k * Dv[k] for k in range(1, len(Dv))])
        Gv = psub(pmul(dN, Dv), pmul(Nv, dD))
        if not Nv: ivs.append((0.0, 0.0)); continue
        if not Gv:
            v = fl(ev(Nv, Fr(0)) / ev(Dv, Fr(0))) if ev(Dv, Fr(0)) else math.inf
            ivs.append((v, v)); continue
        ivs += side_image(Nv, Dv, Gv)
    return ivs


def coverage(ivs):
    """boolean coverage of the circle grid by closed extended intervals; angle = 2 atan(v)"""
    cov = [False] * GRID
    ang = lambda v: math.copysign(math.pi, v) if math.isinf(v) else 2 * math.atan(v)
    for lo, hi in ivs:
        a, b = ang(lo), ang(hi)
        i0 = int(math.floor((a + math.pi) / (2 * math.pi) * GRID))
        i1 = int(math.ceil((b + math.pi) / (2 * math.pi) * GRID))
        for i in range(max(i0, 0), min(i1, GRID - 1) + 1): cov[i % GRID] = True
        if math.isinf(lo) or math.isinf(hi): cov[0] = cov[GRID - 1] = True    # oo
    return cov


def coverage_line(ivs):
    """coverage of the extended line [-oo, +oo] (+oo and -oo distinct); position = atan(value)"""
    cov = [False] * GRID
    ang = lambda v: math.copysign(math.pi / 2, v) if math.isinf(v) else math.atan(v)
    for lo, hi in ivs:
        i0 = int(math.floor((ang(lo) + math.pi / 2) / math.pi * (GRID - 1)))
        i1 = int(math.ceil((ang(hi) + math.pi / 2) / math.pi * (GRID - 1)))
        for i in range(max(i0, 0), min(i1, GRID - 1) + 1): cov[i] = True
    return cov


def signed_ivs(cl):
    f = lambda v: (math.inf if v.sgn > 0 else -math.inf) if isinstance(v, _Inf) else float(v)
    return [(f(lo), f(hi)) for lo, hi in cl["signed_components"]]


def mismatch_line(c1, c2, slack=4):
    diff = 0
    for i in range(GRID):
        if c1[i] != c2[i]:
            if not any(c2[j] != c2[i] for j in range(max(0, i - slack), min(GRID, i + slack + 1))): diff += 1
    return diff / GRID


def exact_ivs(cl):
    out = []
    for lo, hi in cl["components"]:
        f = lambda v: (math.inf if v.sgn > 0 else -math.inf) if isinstance(v, _Inf) else float(v)
        out.append((f(lo), f(hi)))
    return out


def mismatch(c1, c2, slack=4):
    """fraction of the circle where the two coverages differ, ignoring `slack` grid cells
       at each boundary of the exact set (grid rounding)"""
    diff = 0
    for i in range(GRID):
        if c1[i] != c2[i]:
            near = any(c2[(i + k) % GRID] != c2[i] for k in range(-slack, slack + 1))
            if not near: diff += 1
    return diff / GRID


if __name__ == "__main__":
    SCR = sys.argv[1] if len(sys.argv) > 1 else None
    classic = ['x*y/(x^2+y^2)', 'x^2*y/(x^4+y^2)', 'x^3/(x^2+y^2)', 'x^4/((x-y)^2+x^4)',
               '(y-x^3)^2*x/((y-x^3)^2 + x^7)', 'x^6/((y-x-x^2)^2+x^6)', 'x^3/(y-x^2)',
               'x^4/((y^2-2*x^2)^2+x^6)', 'x^2/(x^2+2*y^2)+y^2/(3*x^2+y^2)', 'x*y/(x^2-3*x*y+y^2)',
               '(x^2-y^2)/(x^2+y^2)', '1/(x^2+y^2)', 'x/y']
    gen = random.Random(11)
    def rpoly(n):
        P = {}
        for _ in range(n):
            i, j = gen.randint(0, 5), gen.randint(0, 5)
            P[(i, j)] = P.get((i, j), 0) + gen.choice([-3, -2, -1, 1, 2, 3])
        return " + ".join(f"({v})*x^{i}*y^{j}" for (i, j), v in P.items() if v) or "1"
    funcs = classic + [f"({rpoly(gen.randint(1, 3))})/({rpoly(gen.randint(2, 3))})" for _ in range(int(sys.argv[2]) if len(sys.argv) > 2 else 100)]
    ks = [4, 8, 16, 32]
    worst_final, rows, skipped = 0.0, [], 0
    T = time.time()
    for f in funcs:
        try:
            pr = profile(f)
        except ZeroDivisionError:
            continue
        N, D = pr["N"], pr["D"]
        cl = pr["cluster"]
        ce = coverage(exact_ivs(cl))
        ms = []
        for k in ks:
            try:
                im = loop_image(N, D, Fr(1, 10 ** k))
            except ZeroDivisionError:
                im = None
                print("  skipped (0/0 on the loop):", f, flush=True)
            if im is None: ms = None; break
            ms.append(mismatch(coverage(im), ce))
        if ms is None: skipped += 1; continue
        rows.append((f, cl["text"], ms))
        worst_final = max(worst_final, ms[-1])
        flag = "" if ms[-1] < 0.002 else "   <-- CHECK"
        if True:
            print(f"{f:45s} {cl['text'][:44]:44s} " + "  ".join(f"{m:.4f}" for m in ms) + flag, flush=True)
    n = len(rows)
    conv = sum(1 for _, _, ms in rows if ms[-1] < 0.002)
    print(f"\n{n} functions ({skipped} skipped), {time.time() - T:.0f}s")
    print(f"mismatch at r = 1e-32 below 0.2% of the circle: {conv}/{n};  worst: {worst_final:.4f}")
    print("median mismatch by r:", "  ".join(f"1e-{k}: {sorted(ms[i] for _, _, ms in rows)[n // 2]:.4f}" for i, k in enumerate(ks)))
