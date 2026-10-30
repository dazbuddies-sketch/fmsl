"""Independent check of analytic formulas whose polynomial denominator vanishes along curves
(fmsl.analytic, the curve method) (python benchmarks/curve_crosscheck.py [n]).

The original formula is evaluated directly in 320-digit arithmetic at about 12,000 points within
1e-6 of the origin (the sampler of analytic_crosscheck.py: mixed scales and curves y = c x^k).
Every value must lie in the proven set (+oo and -oo are told apart), the values should fill it,
and deliberately damaged answers must be rejected.
"""
import os, random, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src")); sys.path.insert(0, HERE)
import analytic_crosscheck as ac

KNOWN = ["sin(x)/(x-y)", "sin(x-y)/(x-y)^2", "exp(x)/(x-y)", "(1-cos(x))/(x^2-y^2)", "sin(x)^2/(x^2-y^2)",
         "sin(x*y)/(x^2-y^3)", "(exp(x)-1)/(y-x^2)", "sin(x)/(y^2-x^3)", "atan(x)*y/(x^2-2*y^2)",
         "log(1+x^2)/(x^2-y^2)", "(x^2+sin(y)^2)/(x*y)", "sin(x)*sin(y)/(x^2-y^2)", "cos(x)/(x*y)",
         "(1-cos(x*y))/(x^2*y^2-x^5)", "sinh(x)^2/((x-y)^2)", "tan(x)^2/(x^2-y^4)"]
DENS = ["x-y", "x^2-y^2", "y-x^2", "y^2-x^3", "x^2-2*y^2", "x*y", "(x-y)^2", "x^2-y^4", "x*(x-y)", "y^3-x^2"]
FUN = ["sin({})", "(exp({})-1)", "(1-cos({}))", "log(1+{})", "atan({})", "sinh({})", "cos({})", "exp({})"]


def rand(gen):
    m = gen.choice(["x", "y", "x*y", "x^2", "y^2", "x+y", "x-2*y"])
    num = gen.choice(FUN).format(m)
    if gen.random() < .5: num += "*" + gen.choice(["x", "y", "x^2", "x*y", "y^2"])
    if gen.random() < .3: num += "+" + gen.choice(FUN).format(gen.choice(["x", "y", "x*y"]))
    return f"({num})/({gen.choice(DENS)})"


def damage(comps):
    """wrong answers of the kinds this method could plausibly produce: a missing piece (one
       component dropped, or the wrong sign of infinity), or a finite end moved"""
    import math
    comps = [list(c) for c in comps]
    if len(comps) > 1: return comps[:-1]
    (lo, hi), = comps
    if math.isinf(lo) and math.isinf(hi) and lo < 0 < hi: return [[lo, 0.0]]       # whole line -> half
    if lo == hi: return [[-lo, -hi]] if math.isinf(lo) else [[lo + 0.1, hi + 0.1]]
    return ac.shrink(comps)


def main(n=60):
    gen = random.Random(23)
    fs = KNOWN + [rand(gen) for _ in range(n)]
    t0 = time.time(); tally = {}; caught = tot = 0
    for f in fs:
        try:
            st, info = ac.check(f, ["x", "y"], gen)
            if st not in ("unproven",):
                d, _ = ac.check(f, ["x", "y"], gen, damage=damage)
                if d not in ("skip", "unproven"): tot += 1; caught += d != "ok"
        except Exception as ex:
            st, info = "error", f"{type(ex).__name__}: {ex}"
        tally[st] = tally.get(st, 0) + 1
        if st in ("ok", "coverage") and info != "[-oo, +oo]": tally["(not the whole line)"] = tally.get("(not the whole line)", 0) + 1
        if st != "ok": print(st.upper(), f, info, flush=True)
    print(f"{len(fs)} functions: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())) +
          f"; damaged answers rejected {caught} of {tot}  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 60)
