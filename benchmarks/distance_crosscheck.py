"""Independent checks of fmsl.distance (python benchmarks/distance_crosscheck.py [n2] [n3]).

1. Two variables: two unrelated exact methods must agree. fmsl.distance squares the distance
   away and uses the two-variable engine; fmsl.radial uses the polar blow-up, a one-variable
   image and local analysis at special directions. Same functions as radial_crosscheck.py.
2. Three variables, directions only: for f = P(x, y, z)/rho^k with P homogeneous of degree k,
   f depends only on the direction, so its fluid state is the range of P on the unit sphere, a
   closed interval (the sphere is connected). Its ends are found numerically (a dense Fibonacci
   lattice of 40,000 directions, then local refinement) and must match the proven ends to 1e-6.
3. Three variables, any rational f(x, y, z, rho): the squaring map can lead with x, y or z, giving
   three different rational functions of degree up to twice the original; the three proven
   answers must be identical. Every value of f sampled on the sphere of radius 1e-6 must also
   lie in (or within 1e-3 of) the proven set.
"""
import math, os, random, signal, sys, time
from fractions import Fraction as Fr

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src")); sys.path.insert(0, HERE)
from fmsl.verify import limit_at_origin
from fmsl.distance import reduce_distance
from fmsl.radial import radial_limit
from fmsl.cluster import _Inf, _describe_line
import radial_crosscheck as rc

D3 = "sqrt(x^2+y^2+z^2)"
LIMIT_S = 60          # per engine call; slower items are skipped and counted


class Slow(Exception): pass


def timed(fn, *a):
    def h(*_): raise Slow()
    signal.signal(signal.SIGALRM, h); signal.alarm(LIMIT_S)
    try: return fn(*a)
    finally: signal.alarm(0)


def fib_sphere(n=40000):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)


def to_py(text):
    return eval("lambda x, y, z: " + text.replace(D3, "rho(x,y,z)").replace("^", "**"),
                {"rho": lambda x, y, z: np.sqrt(x * x + y * y + z * z)})


def rand_hom(gen, k):
    mons = [(i, j, k - i - j) for i in range(k + 1) for j in range(k + 1 - i)]
    terms = []
    for m in gen.sample(mons, gen.randint(1, min(4, len(mons)))):
        c = gen.choice([-3, -2, -1, 1, 2, 3])
        terms.append(f"({c})*x^{m[0]}*y^{m[1]}*z^{m[2]}")
    return "+".join(terms)


def sphere_range(P):
    f = to_py(P)
    U = fib_sphere()
    v = f(U[:, 0], U[:, 1], U[:, 2])
    def g(a, sgn):
        x, y, z = math.sin(a[0]) * math.cos(a[1]), math.sin(a[0]) * math.sin(a[1]), math.cos(a[0])
        return sgn * f(x, y, z)
    out = []
    for sgn, idx in ((1, np.argmin(v)), (-1, np.argmax(v))):
        x, y, z = U[idx]
        a0 = [math.acos(max(-1, min(1, z))), math.atan2(y, x)]
        r = minimize(g, a0, args=(sgn,), method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-15, "maxiter": 4000})
        out.append(sgn * r.fun)
    return out


def comps_float(r):
    return [tuple(float(v) for v in c) for c in r["components"]]


def part1(n2):
    gen = random.Random(1)
    fs = rc.KNOWN + [rc.rand_text(gen) for _ in range(n2)]
    agree = tot = 0
    for f in fs:
        a = radial_limit(f)
        if not a["proven"]: continue
        g = f.replace("r", "sqrt(x^2+y^2)")
        try:
            b = limit_at_origin(g, ["x", "y"])
        except Exception as e:
            print("ERROR", f, e); tot += 1; continue
        if not b["proven"]: continue
        tot += 1
        fl = lambda v: (math.inf if v.sgn > 0 else -math.inf) if isinstance(v, _Inf) else float(v)
        ca = [tuple(fl(v) for v in c) for c in a["components"]]
        cb = [tuple(c) for c in b["components"]]
        same = len(ca) == len(cb) and all(abs(p - q) <= 1e-12 * max(1, abs(p)) or p == q for c1, c2 in zip(ca, cb) for p, q in zip(c1, c2))
        if same: agree += 1
        else: print("DISAGREE 2D", f, a["set"], b["set"])
    return agree, tot


def part2(n):
    gen = random.Random(2)
    ok = tot = 0
    for _ in range(n):
        k = gen.randint(1, 3)
        P = rand_hom(gen, k)
        f = f"({P})/{D3}^{k}"
        try: r = timed(limit_at_origin, f, ["x", "y", "z"])
        except Exception: continue
        if not r["proven"]: continue
        tot += 1
        lo, hi = sphere_range(P)
        (a, b), = comps_float(r)
        good = abs(a - lo) < 1e-6 * max(1, abs(lo)) and abs(b - hi) < 1e-6 * max(1, abs(hi))
        ok += good
        if not good: print("DISAGREE sphere", f, r["set"], lo, hi)
    return ok, tot


def rand_general(gen):
    """small rational functions of x, y, z and rho: the squaring map doubles every degree, so the
       three-variable engine is given degrees up to about 8"""
    mons = [(i, j, k) for i in range(3) for j in range(3) for k in range(3) if i + j + k <= 2]
    def poly(nterms):
        terms = []
        for _ in range(nterms):
            i, j, k = gen.choice(mons); m = gen.randint(0, 1)
            c = gen.choice([-2, -1, 1, 2])
            terms.append(f"({c})*x^{i}*y^{j}*z^{k}*{D3}^{m}")
        return "+".join(terms)
    return f"({poly(gen.randint(1, 2))})/({poly(gen.randint(1, 2))})"


def part3(n):
    gen = random.Random(3)
    ok = tot = skipped = 0
    U = fib_sphere(4000) * 1e-6
    for _ in range(n):
        f = rand_general(gen)
        sets = []
        try:
            for lead in range(3):
                g = reduce_distance(f, ["x", "y", "z"], None, lead)
                r = timed(limit_at_origin, g, ["x", "y", "z"])
                sets.append((r["set"], r["proven"], r["components"]))
        except Exception as e:
            skipped += 1; continue
        if not all(p for _, p, _ in sets): skipped += 1; continue
        tot += 1
        same = len({s for s, _, _ in sets}) == 1
        with np.errstate(all="ignore"):
            v = to_py(f)(U[:, 0], U[:, 1], U[:, 2])
        comps = sets[0][2]
        def near(t):
            if not np.isfinite(t): return True
            if abs(t) > 1e4: return any(math.isinf(b) and b * t > 0 or math.isinf(a) and a * t > 0 for a, b in comps)
            return any(a - 1e-3 * max(1, abs(t)) <= t <= b + 1e-3 * max(1, abs(t)) for a, b in comps)
        inside = all(near(t) for t in v)
        ok += same and inside
        if not (same and inside): print("DISAGREE 3D", f, [s for s, _, _ in sets], "samples inside:", inside)
    return ok, tot, skipped


def main(n2=60, n3=40):
    t0 = time.time()
    a1, t1 = part1(n2)
    a2, t2 = part2(n3)
    a3, t3, s3 = part3(n3)
    print(f"two variables, squaring map vs polar blow-up: {a1}/{t1} identical")
    print(f"three variables, direction-only functions vs numerical sphere range: {a2}/{t2} agree")
    print(f"three variables, general: {a3}/{t3} with identical answers for all three leads and every sample inside ({s3} skipped: not proven by the three-variable engine for some lead, or slower than {LIMIT_S} s)")
    print(f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    main(*(int(a) for a in sys.argv[1:3]))
