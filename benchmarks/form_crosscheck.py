"""Independent checks of weighted distances sqrt(Q), Q a positive-definite quadratic form
(python benchmarks/form_crosscheck.py [n]).

1. Different square vectors: the reduction starts from a vector v0 with Q(v0) a rational square;
   different choices give different rational functions, and their proven answers must agree.
2. Direct evaluation: the original formula, sqrt included, evaluated in 80-digit arithmetic at
   about 12,000 points within 1e-6 of the origin (the sampler of analytic_crosscheck.py); every
   value must lie in the proven set, and the values should fill it.
Forms that take no rational square value are refused by FMSL; they are counted, not checked.
"""
import os, random, signal, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src")); sys.path.insert(0, HERE)
from fmsl.verify import limit_at_origin
from fmsl.distance import reduce_distance, DistanceError
import analytic_crosscheck as ac

LIMIT_S = 90


class Slow(Exception): pass


def timed(fn, *a):
    def h(*_): raise Slow()
    signal.signal(signal.SIGALRM, h); signal.alarm(LIMIT_S)
    try: return fn(*a)
    finally: signal.alarm(0)


def rand_form2(gen):
    while True:
        a, c = gen.randint(1, 5), gen.randint(1, 5); b = gen.randint(-4, 4)
        if 4 * a * c - b * b > 0: return f"{a}*x^2+({b})*x*y+{c}*y^2"


def rand_form3(gen):
    return f"{gen.randint(1, 4)}*x^2+{gen.randint(1, 4)}*y^2+{gen.randint(1, 4)}*z^2"


def rand_formula(gen, Q, names):
    R = f"sqrt({Q})"
    mono = lambda: "*".join(f"{v}^{gen.randint(0, 2)}" for v in names)
    shapes = [f"{mono()}*x/{R}^{gen.randint(1, 3)}", f"({R}-{gen.randint(1, 3)}*x)/({R}+y)",
              f"x*y/({R}*({R}+x))", f"({R}-x)/y", f"sin(x*y)/{R}^2", f"x^2*y/({R}^3+y^2)"]
    return gen.choice(shapes)


KNOWN = [("x/sqrt(x^2+2*y^2)", ["x", "y"]), ("y/sqrt(x^2+2*y^2)", ["x", "y"]),
         ("x/sqrt(x^2+x*y+y^2)", ["x", "y"]), ("(sqrt(x^2+4*y^2)-x)/y", ["x", "y"]),
         ("sin(x*y)/sqrt(x^2+2*y^2)^2", ["x", "y"]), ("x/sqrt(x^2+2*y^2+3*z^2)", ["x", "y", "z"]),
         ("z/sqrt(x^2+y^2+2*z^2)", ["x", "y", "z"]), ("x*y/(x^2+2*y^2+3*z^2)^(1/2)", ["x", "y", "z"])]


def main(n=40):
    gen = random.Random(5)
    cases = list(KNOWN)
    for _ in range(n):
        if gen.random() < .8:
            Q = rand_form2(gen); names = ["x", "y"]
        else:
            Q = rand_form3(gen); names = ["x", "y", "z"]
        cases.append((rand_formula(gen, Q, names), names))
    t0 = time.time()
    tally = {}
    for f, names in cases:
        try:
            sets = set()
            for lead in range(3):
                g = reduce_distance(f, names, None, lead)
                r = timed(limit_at_origin, g, names)
                if not r["proven"]: raise Slow()
                sets.add(r["set"])
            if len(sets) != 1:
                status = "LEADS DISAGREE"; print(status, f, sets, flush=True)
            else:
                status, info = ac.check(f, names, gen)
                if status != "ok": print(status.upper(), f, info, flush=True)
        except DistanceError as e:
            status = "refused"
        except Slow:
            status = "not proven in time"
        except Exception as e:
            status = "error"; print("ERROR", f, type(e).__name__, e, flush=True)
        tally[status] = tally.get(status, 0) + 1
    print(f"{len(cases)} formulas: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())) +
          f"  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 40)
