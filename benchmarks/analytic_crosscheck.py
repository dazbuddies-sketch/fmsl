"""Independent check of fmsl.analytic (python benchmarks/analytic_crosscheck.py [n]).

fmsl.analytic replaces sin, cos, exp, log, ... by exact Taylor polynomials and proves, with a
Lojasiewicz bound, that the remainders change no limit value. This check uses none of that: it
evaluates the original formula directly, in 80-digit decimal arithmetic (its own sin, cos, exp,
atan, ... from series; exp, ln and sqrt from the decimal module), at about 12,000 points within 1e-6 of the origin: random points at
mixed scales, and points on curves y = c x^k (k from 1/4 to 4, c from 1e-6 to 1e7) and
z = c x^k in three variables, where limit values are typically attained.

  * soundness: every sampled value lies in the proven set (within 2e-3 on the circle
    [-oo, +oo], measured as atan);
  * coverage: the sampled values fill the proven set (no gap wider than 0.05 inside any
    component, both ends reached within 0.02, every isolated value reached).
Soundness failures would be errors in FMSL. Coverage misses can also mean the sampling did not
reach a value; they are reported separately.
"""
import math, os, random, sys, time
from decimal import Decimal as Dec, getcontext
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from fmsl.verify import limit_at_origin
from fmsl.compiler import Parser

getcontext().prec = 80
ONE = Dec(1)


# ---------------------------------------------------------------- 80-digit elementary functions
def d_exp(x): return x.exp()            # the decimal module's own exp and ln, correctly rounded


def d_log(x):
    if x <= 0: raise ValueError
    return x.ln()


def d_sin_cos(x):
    s, c, t, n = Dec(0), Dec(0), ONE, 0
    while True:
        if abs(t) < Dec(10) ** -85 and n > 2: break
        if n % 4 == 0: c += t
        elif n % 4 == 1: s += t
        elif n % 4 == 2: c -= t
        else: s -= t
        n += 1; t = t * x / n
    return s, c


def d_atan(x):
    if abs(x) > Dec("0.5"):
        return 2 * d_atan(x / (1 + (1 + x * x).sqrt()))
    s, t, n, sg = Dec(0), x, 1, 1
    while abs(t) > Dec(10) ** -85:
        s += sg * t / n; t *= x * x; n += 2; sg = -sg
    return s


def d_asin(x): return d_atan(x / (1 - x * x).sqrt())


FN = {
    "sin": lambda x: d_sin_cos(x)[0], "cos": lambda x: d_sin_cos(x)[1],
    "tan": lambda x: (lambda s, c: s / c)(*d_sin_cos(x)), "exp": d_exp, "log": d_log,
    "sqrt": lambda x: x.sqrt(), "atan": d_atan, "asin": d_asin,
    "sinh": lambda x: (d_exp(x) - d_exp(-x)) / 2, "cosh": lambda x: (d_exp(x) + d_exp(-x)) / 2,
    "tanh": lambda x: (d_exp(x) - d_exp(-x)) / (d_exp(x) + d_exp(-x)),
}


def ev(e, env):
    t = e[0]
    if t == "num": return Dec(e[1].numerator) / Dec(e[1].denominator)
    if t == "var": return env[e[1]]
    if t == "neg": return -ev(e[1], env)
    if t == "+": return ev(e[1], env) + ev(e[2], env)
    if t == "-": return ev(e[1], env) - ev(e[2], env)
    if t == "*": return ev(e[1], env) * ev(e[2], env)
    if t == "/": return ev(e[1], env) / ev(e[2], env)
    if t == "pow":
        b, p = ev(e[1], env), e[2]
        if isinstance(p, int): return b ** p
        return d_exp(d_log(b) * Dec(p.numerator) / Dec(p.denominator))
    if t == "call": return FN[e[1]](ev(e[2], env))
    raise ValueError(t)


# ---------------------------------------------------------------- sampling
def points(n_vars, gen, n=6000):
    pts = []
    rnd = lambda: Dec(str(gen.uniform(1, 10))) * (1 if gen.random() < .5 else -1)
    for _ in range(n // 2):                            # mixed scales
        pts.append([rnd() * Dec(10) ** -gen.randint(7, 16) for _ in range(n_vars)])
    ks = [Fr(1, 4), Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(1), Fr(3, 2), Fr(2), Fr(3), Fr(4)]
    for _ in range(n - n // 2):                        # curves y = c x^k (and z = c' x^k')
        x = Dec(str(gen.uniform(1, 10))) * Dec(10) ** -gen.randint(8, 16)
        if gen.random() < .5: x = -x
        p = [x]
        for _ in range(n_vars - 1):
            k = gen.choice(ks)
            c = Dec(str(gen.uniform(-3, 3))) if gen.random() < .6 else \
                Dec(10) ** gen.randint(-6, 6) * Dec(str(gen.uniform(1, 10))) * gen.choice([-1, 1])
            ax = abs(x)
            p.append(c * (ax.ln() * Dec(k.numerator) / Dec(k.denominator)).exp())
        gen.shuffle(p)
        pts.append(p)
    return [p for p in pts if max(abs(v) for v in p) < Dec("1e-6")]


def arcs(comps):
    a = lambda v: math.copysign(math.pi / 2, v) if math.isinf(v) else math.atan(v)
    return [(a(lo), a(hi)) for lo, hi in comps]


def check(f, names, gen, n=12000, damage=None):
    r = limit_at_origin(f, names)
    if not r["proven"]: return "unproven", r["set"]
    comps = r["components"] if damage is None else damage(r["components"])
    if comps is None: return "skip", None
    A = arcs(comps)
    e = Parser(f).parse()
    vals = []
    for p in points(len(names), gen, n):
        try:
            v = ev(e, dict(zip(names, p)))
        except (ZeroDivisionError, ArithmeticError, ValueError):
            continue
        vals.append(math.atan(float(v)) if abs(v) < Dec("1e300") else math.copysign(math.pi / 2, float(v)))
    bad = [v for v in vals if not any(lo - 2e-3 <= v <= hi + 2e-3 for lo, hi in A)]
    if bad: return "UNSOUND", (r["set"], math.tan(bad[0]))
    miss = []
    for lo, hi in A:
        inside = sorted(v for v in vals if lo - 2e-3 <= v <= hi + 2e-3)
        if not inside or inside[0] > lo + 0.02 or inside[-1] < hi - 0.02 or \
           any(b - a > 0.05 for a, b in zip(inside, inside[1:])):
            miss.append((lo, hi))
    return ("ok" if not miss else "coverage"), r["set"]


KNOWN2 = ["sin(x*y)/(x^2+y^2)", "(1-cos(x*y))/(x^2+y^2)^2", "(exp(x*y)-1)/(x^2+y^2)",
          "log(1+x^2*y)/(x^4+y^2)", "sin(x)*sin(y)/(x^2+y^2)", "(x-sin(x))/(x^2+y^2)",
          "sin(x*y)/sqrt(x^2+y^2)", "(x-sin(x))/sqrt(x^2+y^2)^3", "tan(x*y)/(x^2+y^2)",
          "(exp(x)-1-x)/(x^2+y^2)", "sin(x^2+y^2)/(x^2+y^2)", "(sqrt(1+x)-1-x/2)/(x^2+y^2)",
          "atan(x*y)/(x^2+y^4)", "(cos(x)-cos(y))/(x^2+y^2)", "sin(x)*y/(x^2+y^2)",
          "(sin(x)-x)*y/(x^6+y^2)", "sinh(x)*y^2/(x^2+y^4)", "log(1+x*y)/(sin(x)^2+sin(y)^2)",
          "(exp(x*y)-cos(x*y))/(x^2+x*y+y^2)", "asin(x^2*y)/(x^4+y^2)"]
KNOWN3 = ["sin(x*y*z)/(x^2+y^2+z^2)^(3/2)", "(exp(x*y*z)-1)/(x^2+y^2+z^2)", "sin(x)*y*z/(x^2+y^2+z^2)",
          "(1-cos(x))*y/(x^2+y^2+z^2)^(3/2)"]

FUN = [("sin({})", 1), ("(exp({})-1)", 1), ("(1-cos({}))", 2), ("log(1+{})", 1), ("atan({})", 1),
       ("tan({})", 1), ("sinh({})", 1)]
DEN2 = ["x^2+y^2", "x^4+y^2", "x^2+y^6", "(x^2+y^2)^2", "x^2+x*y+y^2", "sin(x)^2+sin(y)^2",
        "(y-x^2)^2+x^4", "exp(x^2+y^2)-1", "x^2+y^2+x^3", "sqrt(x^2+y^2)^3"]


def rand2(gen):
    def mono():
        i, j = gen.randint(0, 2), gen.randint(0, 2)
        if i + j == 0: i = 1
        return f"x^{i}*y^{j}"
    num = []
    for _ in range(gen.randint(1, 2)):
        fn, _ = gen.choice(FUN)
        term = fn.format(mono())
        if gen.random() < .5: term += "*" + mono()
        num.append(f"({gen.choice([-2, -1, 1, 2])})*{term}")
    return f"({'+'.join(num)})/({gen.choice(DEN2)})"


def shrink(comps):
    """a damaged answer: the widest finite component loses 10% at its upper end, or a single
       finite value moves by 0.1 (None if there is nothing finite to damage)"""
    fin = [(i, lo, hi) for i, (lo, hi) in enumerate(comps) if not (math.isinf(lo) or math.isinf(hi))]
    if not fin: return None
    i, lo, hi = max(fin, key=lambda t: t[2] - t[1])
    out = [list(c) for c in comps]
    out[i] = [lo, hi - 0.1 * (hi - lo)] if hi > lo else [lo + 0.1, hi + 0.1]
    return out


def main(n=40):
    gen = random.Random(11)
    t0 = time.time()
    tally = {}
    cases = [(f, ["x", "y"]) for f in KNOWN2] + [(f, ["x", "y", "z"]) for f in KNOWN3] + \
            [(rand2(gen), ["x", "y"]) for _ in range(n)]
    for f, names in cases:
        try:
            status, info = check(f, names, gen)
        except Exception as ex:
            status, info = "error", f"{type(ex).__name__}: {ex}"
        tally[status] = tally.get(status, 0) + 1
        if status not in ("ok",): print(status.upper(), f, info)
    print(f"{len(cases)} functions: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())) +
          f"  ({time.time() - t0:.0f} s)")
    caught = tot = 0                                   # the check must reject damaged answers
    for f, names in cases[:len(KNOWN2) + len(KNOWN3)]:
        st, _ = check(f, names, gen, damage=shrink)
        if st in ("skip", "unproven"): continue
        tot += 1; caught += st != "ok"
    print(f"damaged answers rejected: {caught} of {tot}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 40)
