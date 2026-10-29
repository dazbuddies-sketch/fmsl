"""Independent check of limits at infinity and from one side (fmsl.regions)
(python benchmarks/region_crosscheck.py [n]).

The original formula is evaluated directly, in 320-digit arithmetic, at about 12,000 points:
  * at infinity: points with |q| > 1e36, at mixed scales and on curves (the small-point sampler of
    analytic_crosscheck.py, scaled by 1e-30 and pushed out by q -> q/|q|^2, so every direction and
    many curve shapes are represented);
  * from one side: the same small points, with each restricted coordinate given the required sign.
No change of variables is used. Every sampled value must lie in the proven set (within 2e-3 on the
circle [-oo, +oo], as atan), and the values should fill it; damaged answers must be rejected.
"""
import math, os, random, sys, time
from decimal import Decimal as Dec

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src")); sys.path.insert(0, HERE)
from fmsl.verify import limit_at
from fmsl.compiler import Parser
import analytic_crosscheck as ac

KNOWN_INF = ["x*y/(x^2+y^2)", "(x^2+y)/(x^2+y^2)", "x^3*y/(x^6+y^2)", "(x^2-y^2)/(x^2+y^2+1)",
             "x/(x^2+y^2)", "x*y/sqrt(x^2+y^2)", "sqrt(x^2+2*y^2)/(x+2*y^2)", "(x+y)/sqrt(x^2+y^2)",
             "x^2*y/(x^4+y^2+1)", "(x*y+1)/(x^2+y^2+1)", "x*y*z/(x^2+y^2+z^2)^(3/2)", "(x+y+z)^2/(x^2+y^2+z^2)"]
KNOWN_SIDE = [("x*y/(x^2+y^2)", {"x": "+", "y": "+"}), ("x*y/(x^2+y^2)", {"x": "+", "y": "-"}),
              ("x^3*y/(x^6+y^2)", {"y": "+"}), ("y/x", {"x": "+", "y": "+"}), ("(x-y)/(x+y)", {"x": "+", "y": "+"}),
              ("sin(x*y)/(x^2+y^2)", {"x": "+"}), ("x^2*y/(x^4+y^2)", {"y": "-"}), ("x*y*z/(x^2+y^2+z^2)^2", {"z": "+"})]


def rand_rational(gen, names):
    def mono():
        return "*".join(f"{v}^{gen.randint(0, 3)}" for v in names)
    num = "+".join(f"({gen.choice([-2, -1, 1, 2])})*{mono()}" for _ in range(gen.randint(1, 2)))
    den = "+".join([f"{v}^{2 * gen.randint(1, 3)}" for v in names] +
                   [f"({gen.choice([-1, 1])})*{mono()}" for _ in range(gen.randint(0, 1))])
    return f"({num})/({den})"


def values(f, names, pts):
    e = Parser(f).parse()
    out = []
    for p in pts:
        try:
            v = ac.ev(e, dict(zip(names, p)))
        except (ZeroDivisionError, ArithmeticError, ValueError):
            continue
        out.append(math.atan(float(v)) if abs(v) < Dec("1e300") else math.copysign(math.pi / 2, float(v)))
    return out


def verdict(vals, comps):
    A = ac.arcs(comps)
    if any(not any(lo - 2e-3 <= v <= hi + 2e-3 for lo, hi in A) for v in vals): return "UNSOUND"
    for lo, hi in A:
        inside = sorted(v for v in vals if lo - 2e-3 <= v <= hi + 2e-3)
        if not inside or inside[0] > lo + 0.02 or inside[-1] < hi - 0.02 or \
           any(b - a > 0.05 for a, b in zip(inside, inside[1:])):
            return "coverage"
    return "ok"


def far(pts, scale=Dec(10) ** -30):
    """push small points out: q = s p / |s p|^2 with s = 1e-30, so |q| > 1e36 (limits at infinity
       can converge slowly, like |q|^(-1/3))"""
    out = []
    for p in pts:
        p = [v * scale for v in p]
        n2 = sum(v * v for v in p)
        out.append([v / n2 for v in p])
    return out


def signed(pts, names, side):
    return [[abs(v) if side.get(n) == "+" else (-abs(v) if side.get(n) == "-" else v) for v, n in zip(p, names)]
            for p in pts]


def main(n=30):
    gen = random.Random(17)
    cases = [(f, "oo", None) for f in KNOWN_INF] + [(f, None, sd) for f, sd in KNOWN_SIDE]
    for _ in range(n):
        names = ["x", "y"]
        cases.append((rand_rational(gen, names), "oo", None))
    for _ in range(n // 2):
        names = ["x", "y"]
        sd = gen.choice([{"x": "+"}, {"y": "-"}, {"x": "+", "y": "+"}, {"x": "-", "y": "+"}])
        cases.append((rand_rational(gen, names), None, sd))
    t0 = time.time(); tally = {}; damaged = caught = 0
    for f, pt, sd in cases:
        names = ["x", "y", "z"] if "z" in f else ["x", "y"]
        try:
            r = limit_at(f, "oo", names) if pt == "oo" else limit_at(f, [0] * len(names), names, side=sd)
            if not r["proven"]:
                st = "unproven"
            else:
                pts = ac.points(len(names), gen, 12000)
                pts = far(pts) if pt == "oo" else signed(pts, names, sd)
                vals = values(f, names, pts)
                st = verdict(vals, r["components"])
                bad = ac.shrink(r["components"])
                if bad is not None:
                    damaged += 1; caught += verdict(vals, bad) != "ok"
            if st not in ("ok",): print(st.upper(), f, pt or sd, r["set"], flush=True)
        except Exception as ex:
            st = "error"; print("ERROR", f, pt or sd, type(ex).__name__, ex, flush=True)
        tally[st] = tally.get(st, 0) + 1
    print(f"{len(cases)} limits: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())) +
          f"; damaged answers rejected {caught} of {damaged}  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 30)
