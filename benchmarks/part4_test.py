"""Test the refined Part 4 conjecture (docs/theory/n-variables.md) on functions whose numerator and
denominator both vanish on the z-axis without the clean condition.

Two exact quantities are compared:
  J      the values f approaches tangentially to the axis, i.e. with |(x, y)| = o(|z|): the fluid
         state at 0 of g(u, v, z) = f(z u, z v, z), computed exactly by the three-variable engine;
  slices the fluid states at points (0, 0, t) on the axis, for t = +-10^-6 (approximating their
         limit as t -> 0), also exact.
Conjecture: J is the limit of the slice fluid states, i.e. intermediate scales |z|^a, 1 < a < oo,
add no values. A disagreement beyond rounding of the t = 1e-6 approximation is a counterexample.
Families: random (z-dependent coefficients) and adversarial (denominators and numerators with
structure at intermediate scales, such as x^2 + y^2 - z^3 or x^2 + z*y^2).
"""
import math, os, random, signal, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from fmsl.verify import limit_at, limit_at_origin


class Slow(Exception): pass


def run(fn, *a, cap=40):
    def h(*_): raise Slow()
    signal.signal(signal.SIGALRM, h); signal.alarm(cap)
    try:
        r = fn(*a); return r if r["proven"] else None
    except (Slow, Exception):
        return None
    finally:
        signal.alarm(0)


def arcs(comps, big=1e4):
    a = lambda v: math.copysign(math.pi / 2, v) if math.isinf(v) or abs(v) > big else math.atan(v)
    out = sorted((a(lo), a(hi)) for lo, hi in comps)
    merged = []
    for lo, hi in out:
        if merged and lo <= merged[-1][1] + 1e-6: merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else: merged.append((lo, hi))
    return merged


def close(A, B, tol=3e-3):
    def inside(x, S): return any(lo - tol <= x <= hi + tol for lo, hi in S)
    ends = lambda S: [x for iv in S for x in iv]
    return all(inside(x, B) for x in ends(A)) and all(inside(x, A) for x in ends(B))


def substitute(f):
    return f.replace("x", "(z*U)").replace("y", "(z*V)").replace("U", "x").replace("V", "y")


def rand(gen, adversarial):
    c = lambda: gen.choice([-2, -1, 1, 2])
    zc = lambda: gen.choice(["1", "z", "(1+z)", "(2-z)", "z^2"])
    def quad():
        return "+".join(f"({c()})*{m}*{zc()}" for m in gen.sample(["x^2", "x*y", "y^2"], gen.randint(1, 3)))
    if not adversarial:
        return f"({quad()})/({quad()})"
    extra = gen.choice(["z^3", "z^4", "z^2*x", "z*x^2", "z^3*x", "z^2*y"])
    P = quad() + gen.choice(["", f"+({c()})*{extra}*x", f"+({c()})*x*y*z"])
    Q = gen.choice([f"x^2+y^2-({extra})*{gen.choice(['1','x','y'])}", f"x^2+z*y^2", f"y^2+z^2*x^2+({c()})*x*y*z",
                    quad() + f"+({c()})*z^3*x"])
    return f"({P})/({Q})"


def main(n=60):
    gen = random.Random(31)
    tally = {"agree": 0, "COUNTEREXAMPLE": 0, "undecided": 0}
    for i in range(n):
        f = rand(gen, adversarial=i % 2 == 1)
        rg = run(limit_at_origin, substitute(f), ["x", "y", "z"])
        rs = [run(limit_at, f, [0, 0, t], ["x", "y", "z"]) for t in ("1/1000000", "-1/1000000")]
        if rg is None or any(r is None for r in rs):
            tally["undecided"] += 1; continue
        J = arcs([tuple(c) for c in rg["components"]])
        S = arcs([tuple(c) for r in rs for c in r["components"]])
        ok = close(J, S)
        tally["agree" if ok else "COUNTEREXAMPLE"] += 1
        print("ok " if ok else "XX ", f, "| J =", rg["set"], "| slices:", [r["set"] for r in rs], flush=True)
    print(tally)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 60)
