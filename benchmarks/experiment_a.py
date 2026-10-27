"""Experiment A of MJ Celedon's FMS research brief: singularity classification, head to head.

The brief asks for at least 1,000 expressions covering ten kinds of local behaviour, compared
with established symbolic and numerical methods against independently established ground truth,
measuring correct and false classifications, unresolved cases and time.

Ground truth comes from how each expression is BUILT (known series, AM-GM bounds, eigenvalues of
quadratic forms, explicit constructions), never from FMSL or from a baseline.

Methods
  FMSL      fmsl (one variable: the compiler; several variables: fmsl.multivar / fmsl.space)
  SymPy     sympy.limit from each side (one variable).  For several variables SymPy has no
            limit routine; the baseline is the standard practice of limits along paths
            (y = 0, x = 0, y = x, y = -x, y = x^2, ...; in 3-D also the diagonals).
  float     naive double-precision probing (labelled as such; the brief warns that beating it
            alone proves nothing).

Scoring
  one variable:      the pair (left limit, right limit), each a number, +oo, -oo or 'osc'
                     (no limit).  Correct only if both sides are right.
  several variables: VERDICT correct (limit exists or not) and SET correct (the exact set of
                     limit values, the fluid state, on the extended line).
  unresolved:        the method gave no answer, said it could not decide, timed out (20 s),
                     or does not support the expression.  Never counted as correct.

Run:  python benchmarks/experiment_a.py [--quick] [--no-sympy]
Writes benchmarks/results/experiment_a.json and benchmarks/results/experiment_a.md
"""
import json, math, os, random, re, signal, sys, time
from fractions import Fraction as Fr

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
INF, NINF, OSC = "+oo", "-oo", "osc"
TIMEOUT = 20


class Timeout(Exception):
    pass


def _alarm(*a): raise Timeout()


# ============================================================ generators (ground truth by construction)
def sgn(v): return (v > 0) - (v < 0)


def gen_one_variable(gen, n_per):
    """list of (category, expression in x at x -> 0, (left, right) truth, fmsl_supported)"""
    out = []
    def poly(deg, c0):                      # random polynomial with constant term c0 != 0
        cs = [c0] + [gen.choice([-3, -2, -1, 0, 1, 2, 3]) for _ in range(deg)]
        return " + ".join(f"({c})*x^{i}" for i, c in enumerate(cs) if c)
    def expand_mul(p, m):                    # x^m * p  as an expanded sum
        terms = []
        for t in p.split(" + "):
            c, e = re.match(r"\((-?\d+)\)\*x\^(\d+)", t).groups()
            terms.append(f"({c})*x^{int(e) + m}")
        return " + ".join(terms)
    # 1. removable singularities (rational, hidden common factor; and classical series)
    for _ in range(n_per):
        if gen.random() < 0.5:
            a, b = gen.choice([-3, -2, -1, 1, 2, 3]), gen.choice([-3, -2, -1, 1, 2, 3])
            m = gen.randint(1, 3)
            N, D = expand_mul(poly(gen.randint(0, 2), a), m), expand_mul(poly(gen.randint(0, 2), b), m)
            v = Fr(a, b)
            out.append(("1 removable", f"({N})/({D})", (v, v), True))
        else:
            k = gen.choice([1, 2, 3, Fr(1, 2)])
            kind = gen.randint(0, 2)
            if kind == 0: e, v = f"(exp({k}*x)-1)/x", Fr(k)
            elif kind == 1: e, v = f"(1-cos({k}*x))/x^2", Fr(k) ** 2 / 2
            else: e, v = f"sin({k}*x)/x", Fr(k)
            out.append(("1 removable", e, (v, v), True))
    # 2. poles of different orders (even orders: same sign on both sides)
    for _ in range(n_per):
        a, b = gen.choice([-3, -2, -1, 1, 2, 3]), gen.choice([-3, -2, -1, 1, 2, 3])
        m, k = gen.randint(0, 2), gen.choice([2, 4])
        N, D = expand_mul(poly(gen.randint(0, 2), a), m), expand_mul(poly(gen.randint(0, 2), b), m + k)
        s = INF if a * b > 0 else NINF
        out.append(("2 pole (even order)", f"({N})/({D})", (s, s), True))
    # 3. one-sided divergence (odd orders: opposite signs)
    for _ in range(n_per):
        a, b = gen.choice([-3, -2, -1, 1, 2, 3]), gen.choice([-3, -2, -1, 1, 2, 3])
        m, k = gen.randint(0, 2), gen.choice([1, 3])
        if gen.random() < 0.6:
            N, D = expand_mul(poly(gen.randint(0, 2), a), m), expand_mul(poly(gen.randint(0, 2), b), m + k)
            e = f"({N})/({D})"; r = a * b
        else:
            e = f"({a})*cos(x)/sin({b}*x)^{k}"; r = a * (b ** k) * 1   # a / (b x)^k
            r = a * b
        right = INF if r > 0 else NINF
        left = right if k % 2 == 0 else (NINF if right == INF else INF)
        out.append(("3 one-sided divergence", e, (left, right), True))
    # 4. oscillatory: x^m sin(k/x)
    for _ in range(n_per):
        k, m = gen.choice([1, 2, 3]), gen.choice([-1, 0, 1, 2])
        e = f"x^{m}*sin({k}/x)" if m else f"sin({k}/x)"
        if m < 0: e = f"sin({k}/x)/x"
        truth = (Fr(0), Fr(0)) if m > 0 else (OSC, OSC)
        out.append(("4 oscillatory", e, truth, True))
    # 5. piecewise / absolute value
    for _ in range(n_per):
        k = gen.choice([1, 2, 3])
        kind = gen.randint(0, 2)
        if kind == 0: e, t = f"Abs({k}*x)/x", (Fr(-k), Fr(k))
        elif kind == 1: e, t = f"x/Abs(x)^{2}", (NINF, INF)
        else: e, t = f"(Abs(x)+{k}*x)/x", (Fr(k - 1), Fr(k + 1))
        out.append(("5 piecewise", e, t, False))
    # 6. catastrophic cancellation (Taylor remainders)
    for _ in range(n_per):
        k = gen.choice([1, 2, 3])
        kind = gen.randint(0, 2)
        if kind == 0: e, v = f"(exp({k}*x)-1-{k}*x-{k}^2*x^2/2)/x^3", Fr(k ** 3, 6)
        elif kind == 1: e, v = f"(sin({k}*x)-{k}*x)/x^3", Fr(-k ** 3, 6)
        else: e, v = f"(cos({k}*x)-1+{k}^2*x^2/2)/x^4", Fr(k ** 4, 24)
        out.append(("6 catastrophic cancellation", e, (v, v), True))
    # 7. numerically ill-conditioned (deep cancellation, high order)
    for _ in range(n_per):
        m = gen.choice([2, 3, 4])
        kind = gen.randint(0, 1)
        if kind == 0: e, v = f"(1-cos(x^{m}))/x^{2 * m}", Fr(1, 2)
        else: e, v = f"(exp(x^{m})-1-x^{m})/x^{2 * m}", Fr(1, 2)
        out.append(("7 ill-conditioned", e, (v, v), True))
    # 9. essential singularities
    for _ in range(n_per):
        k, n = gen.choice([1, 2, 3]), gen.choice([0, 1, 2, 3])
        kind = gen.randint(0, 1)
        if kind == 0:
            e = f"exp(-{k}/x^2)/x^{n}" if n else f"exp(-{k}/x^2)"; t = (Fr(0), Fr(0))
        else:
            e = f"exp({k}/x)"; t = (Fr(0), INF)
        out.append(("9 essential", e, t, True))
    return out


def _eig2(a, b, c):
    """eigenvalues of [[a, b/2], [b/2, c]] as floats (truth for (a x^2 + b xy + c y^2)/(x^2+y^2))"""
    tr, d = a + c, math.sqrt((a - c) ** 2 + b * b)
    return ((tr - d) / 2, (tr + d) / 2)


def gen_several_variables(gen, n_per):
    """list of (category, expression, truth as a list of (lo, hi) on the extended line, nvars)"""
    out = []
    # 10. multivariable, path-dependent
    for _ in range(n_per):
        kind = gen.randint(0, 5)
        if kind == 0:                                  # quadratic form ratio: [lambda_min, lambda_max]
            a, b, c = gen.randint(-3, 3), gen.choice([-3, -2, -1, 1, 2, 3]), gen.randint(-3, 3)
            lo, hi = _eig2(a, b, c)
            out.append(("10 multivariable", f"(({a})*x^2+({b})*x*y+({c})*y^2)/(x^2+y^2)", [(lo, hi)], 2))
        elif kind == 1:                                # u v/(u^2+v^2) with u = x^p, v = y^q (AM-GM)
            p, q, k = gen.choice([1, 2, 3]), gen.choice([1, 2, 3]), gen.choice([1, 2, 3])
            if p % 2 == 0 and q % 2 == 0: t = [(0.0, k / 2)]
            else: t = [(-k / 2, k / 2)]
            out.append(("10 multivariable", f"{k}*x^{p}*y^{q}/(x^{2 * p}+y^{2 * q})", t, 2))
        elif kind == 2:                                # removable: polynomial over nonvanishing
            a, b = gen.choice([1, 2, 3]), gen.choice([1, 2, 3])
            v = gen.choice([-2, -1, 1, 2, 3])
            out.append(("10 multivariable", f"({v}*{a}+x*y)/({a}+x^2+y^2)", [(float(v), float(v))], 2))
        elif kind == 3:                                # limit 0 with a 0/0 point
            a, b = gen.choice([1, 2]), gen.choice([1, 2])
            out.append(("10 multivariable", f"x^{2 * a}*y^{2 * b}/(x^2+y^2)", [(0.0, 0.0)], 2))
        elif kind == 4:                                # 3-D diagonal quadratic form
            a, b, c = gen.randint(-3, 3), gen.randint(-3, 3), gen.randint(-3, 3)
            if a == b == c: c += 1
            t = [(float(min(a, b, c)), float(max(a, b, c)))]
            out.append(("10 multivariable", f"(({a})*x^2+({b})*y^2+({c})*z^2)/(x^2+y^2+z^2)", t, 3))
        else:                                          # 3-D: xyz-type, limit 0
            k = gen.choice([1, 2, 3])
            out.append(("10 multivariable", f"{k}*x*y*z/(x^2+y^2+z^2)", [(0.0, 0.0)], 3))
    # 8. parameter-dependent behaviour: instances of families at chosen parameter values
    for _ in range(n_per):
        kind = gen.randint(0, 2)
        if kind == 0:                                  # (x^2 - t y^2)/(x^2+y^2): between -t and 1
            t = gen.choice([-3, -2, -1, Fr(-1, 2), 0, Fr(1, 2), 1, 2, 3])
            lo, hi = sorted((float(-t), 1.0))
            out.append(("8 parameter-dependent", f"(x^2-({t})*y^2)/(x^2+y^2)", [(lo, hi)], 2))
        elif kind == 1:                                # x^2 y/(x^4 + t y^2): +-1/(2 sqrt t), or all
            t = gen.choice([-2, -1, 0, Fr(1, 4), 1, 2, 4])
            tr = [(-1 / (2 * math.sqrt(t)), 1 / (2 * math.sqrt(t)))] if t > 0 else [(-math.inf, math.inf)]
            out.append(("8 parameter-dependent", f"x^2*y/(x^4+({t})*y^2)", tr, 2))
        else:                                          # (t x^2 + y^2 - xy)/(x^2+y^2): eigenvalues
            t = gen.choice([-2, -1, 0, Fr(1, 4), 1, 2])
            lo, hi = _eig2(float(t), -1, 1)
            out.append(("8 parameter-dependent", f"(({t})*x^2+y^2-x*y)/(x^2+y^2)", [(lo, hi)], 2))
    return out


# ============================================================ grading helpers
def same_value(a, b, tol=1e-9):
    if isinstance(a, str) or isinstance(b, str): return a == b
    a, b = float(a), float(b)
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def same_sets(A, B, tol=1e-7):
    if A is None or B is None or len(A) != len(B): return False
    def f(v): return v if isinstance(v, float) else float(v)
    for (a1, b1), (a2, b2) in zip(sorted(A, key=lambda p: f(p[0])), sorted(B, key=lambda p: f(p[0]))):
        for u, v in ((a1, a2), (b1, b2)):
            u, v = f(u), f(v)
            if math.isinf(u) or math.isinf(v):
                if u != v: return False
            elif abs(u - v) > tol * max(1.0, abs(u)): return False
    return True


def verdict_of(S):
    return "LIMIT" if S and len(S) == 1 and same_value(S[0][0], S[0][1]) and not math.isinf(float(S[0][0])) else "NO LIMIT"


# ============================================================ FMSL
def fmsl_one(expr):
    from fmsl.compiler import run
    lines = []
    prog = f"fluid f(x) = {expr}\nat x -> 0: f\n"
    run(prog, out=lines.append)
    txt = " ".join(lines)
    if "ERROR" in txt: return None
    if "ESSENTIAL" in txt:                                # check first: its text mentions poles
        if "UNPROVEN" in txt or "not decided" in txt: return None
        if "limit 0 on both sides" in txt: return (Fr(0), Fr(0))
        if "oscillates" in txt: return (OSC, OSC)
        def side(sd):
            mm = re.search(re.escape(sd) + r": ([^ ;(]+)", txt)
            if not mm: return None
            return {"0": Fr(0), "+inf": INF, "-inf": NINF, "+oo": INF, "-oo": NINF}.get(mm.group(1))
        l, r = side("x->0-"), side("x->0+")
        return None if l is None or r is None else (l, r)
    m = re.search(r"removable: value (\S+)", txt)
    if m and "UNPROVEN" not in txt:
        v = m.group(1)
        if "vanishes" in txt: v = "0"
        v = Fr(v) if re.fullmatch(r"-?\d+(/\d+)?", v) else None
        return None if v is None else (v, v)
    m = re.search(r"pole of order (\d+).*\(lead ([^)]+)\)", txt)
    if m:
        n, c = int(m.group(1)), m.group(2).strip()
        try: c = float(Fr(c))
        except ValueError: return None
        right = INF if c > 0 else NINF
        left = right if n % 2 == 0 else (NINF if right == INF else INF)
        return (left, right)
    return None


def fmsl_several(expr, nv):
    if nv == 2:
        from fmsl.multivar import profile
        r = profile(expr)
        cl = r["cluster"]
        if not cl.get("exact", True) or not cl.get("complete", True): return None
        comps = cl.get("signed_components")
    else:
        from fmsl.space import profile3
        r = profile3(expr)
        if not r["complete"] and not r["signed_text"] == "[-oo, +oo]": return None
        comps = r["signed_components"]
    from fmsl.cluster import _Inf
    f = lambda v: (math.inf if v.sgn > 0 else -math.inf) if isinstance(v, _Inf) else float(v)
    return [(f(a), f(b)) for a, b in comps]


# ============================================================ SymPy baselines
def sympy_expr(expr, nv=1):
    import sympy as sp
    x, y, z = sp.symbols("x y z", real=True)
    return sp.sympify(expr.replace("^", "**"), locals={"x": x, "y": y, "z": z, "Abs": sp.Abs}), (x, y, z)


def _sympy_val(L):
    import sympy as sp
    if L is sp.oo: return INF
    if L is sp.S.NegativeInfinity: return NINF
    if isinstance(L, sp.AccumBounds) or L is sp.nan or L.has(sp.AccumBounds): return OSC
    if L.is_number and L.is_real:
        try: return Fr(str(sp.nsimplify(L))) if L.is_rational else float(L)
        except (ValueError, TypeError): return float(L)
    return None


def sympy_one(expr):
    import sympy as sp
    e, (x, _, _) = sympy_expr(expr)
    return (_sympy_val(sp.limit(e, x, 0, "-")), _sympy_val(sp.limit(e, x, 0, "+")))


def sympy_paths(expr, nv):
    """limits along standard paths; LIMIT if all agree (finite), else the set of path values"""
    import sympy as sp
    e, (x, y, z) = sympy_expr(expr, nv)
    s = sp.symbols("s", positive=True)
    if nv == 2:
        paths = [(s, 0), (0, s), (s, s), (s, -s), (s, s ** 2), (-s, s ** 2), (s, 2 * s), (s ** 2, s)]
        vals = [_sympy_val(sp.limit(e.subs({x: a, y: b}), s, 0, "+")) for a, b in paths]
    else:
        paths = [(s, 0, 0), (0, s, 0), (0, 0, s), (s, s, s), (s, -s, s), (s, s, 0), (s, 0, s), (0, s, s)]
        vals = [_sympy_val(sp.limit(e.subs({x: a, y: b, z: c}), s, 0, "+")) for a, b, c in paths]
    if any(v is None for v in vals): return None
    pts = []
    for v in vals:
        if v == OSC: return [(-math.inf, math.inf)]          # cannot tell; reported as everything
        fv = math.inf if v == INF else (-math.inf if v == NINF else float(v))
        if not any(same_value(fv, p) if not math.isinf(fv) else fv == p for p in pts): pts.append(fv)
    return [(p, p) for p in sorted(pts)]


# ============================================================ naive floats
def float_one(expr):
    env = {"exp": math.exp, "sin": math.sin, "cos": math.cos, "Abs": abs}
    code = expr.replace("^", "**")
    def ev(v):
        try: return eval(code, {"__builtins__": {}}, dict(env, x=v))
        except (ZeroDivisionError, OverflowError, ValueError): return math.nan
    def side(sg):
        a, b = ev(sg * 1e-4), ev(sg * 1e-7)
        if math.isnan(a) or math.isnan(b): return None
        if abs(b) > 1e6: return INF if b > 0 else NINF
        if abs(a - b) <= 1e-3 * max(1, abs(b)): return b
        return OSC
    return (side(-1), side(1))


# ============================================================ main
def timed(fn, *a):
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(TIMEOUT)
    t = time.time()
    try:
        r = fn(*a)
    except Timeout:
        r = "timeout"
    except Exception as exc:                    # any failure of a method is "unresolved"
        r = None
    finally:
        signal.alarm(0)
    return r, time.time() - t


def main():
    quick = "--quick" in sys.argv
    use_sympy = "--no-sympy" not in sys.argv
    if use_sympy:
        try:
            import sympy  # noqa: F401
        except ImportError:
            use_sympy = False
    gen = random.Random(1729)
    n_per = 12 if quick else 112                 # 9 one-variable categories, 2 several-variable
    one = gen_one_variable(gen, n_per)
    sev = gen_several_variables(gen, n_per)
    methods = ["FMSL"] + (["SymPy"] if use_sympy else []) + ["float"]
    stats = {}                                   # (category, method) -> counters
    def add(cat, m, key, dt=0.0):
        d = stats.setdefault((cat, m), {"n": 0, "correct": 0, "wrong": 0, "unresolved": 0, "time": 0.0,
                                        "verdict_ok": 0})
        if key: d[key] += 1
        d["time"] += dt
    rows = []
    for cat, expr, truth, supported in one:
        for m in methods:
            if m == "FMSL":
                ans, dt = timed(fmsl_one, expr) if supported else (None, 0.0)
            elif m == "SymPy":
                ans, dt = timed(sympy_one, expr)
            else:
                ans, dt = timed(float_one, expr)
            stats.setdefault((cat, m), None)
            d = stats.get((cat, m)) or {"n": 0, "correct": 0, "wrong": 0, "unresolved": 0, "time": 0.0, "verdict_ok": 0}
            stats[(cat, m)] = d
            d["n"] += 1; d["time"] += dt
            if ans in (None, "timeout") or any(v is None for v in ans):
                d["unresolved"] += 1; res = "unresolved"
            elif all(same_value(a, t, 1e-6 if m == "float" else 1e-9) for a, t in zip(ans, truth)):
                d["correct"] += 1; res = "correct"
            else:
                d["wrong"] += 1; res = "wrong"
            rows.append({"cat": cat, "expr": expr, "method": m, "result": res,
                         "answer": str(ans), "truth": str(truth), "time": round(dt, 4)})
    for cat, expr, truth, nv in sev:
        for m in methods:
            if m == "float": continue
            if m == "FMSL": ans, dt = timed(fmsl_several, expr, nv)
            else: ans, dt = timed(sympy_paths, expr, nv)
            d = stats.get((cat, m)) or {"n": 0, "correct": 0, "wrong": 0, "unresolved": 0, "time": 0.0, "verdict_ok": 0}
            stats[(cat, m)] = d
            d["n"] += 1; d["time"] += dt
            if ans in (None, "timeout"):
                d["unresolved"] += 1; res = "unresolved"
            else:
                if verdict_of(ans) == verdict_of(truth): d["verdict_ok"] += 1
                if same_sets(ans, truth): d["correct"] += 1; res = "correct"
                else: d["wrong"] += 1; res = "wrong"
            rows.append({"cat": cat, "expr": expr, "method": m, "result": res,
                         "answer": str(ans), "truth": str(truth), "time": round(dt, 4)})
    # ---- report
    cats = sorted({c for c, _ in stats})
    out = ["# Experiment A: singularity classification, head to head", "",
           f"{len(one)} one-variable and {len(sev)} several-variable expressions "
           f"({len(one) + len(sev)} total); ground truth by construction; time limit {TIMEOUT} s each.",
           "Unresolved = no answer, undecided, timed out, or not supported. Several variables: "
           "'correct' means the exact set of limit values; 'verdict' means limit / no limit only.", "",
           "| category | method | n | correct | wrong | unresolved | verdict right | time (s) |",
           "|---|---|---|---|---|---|---|---|"]
    tot = {}
    for c in cats:
        for m in methods:
            d = stats.get((c, m))
            if not d: continue
            vr = f"{d['verdict_ok']}" if c.startswith(("8", "10")) else "-"
            out.append(f"| {c} | {m} | {d['n']} | {d['correct']} | {d['wrong']} | {d['unresolved']} | {vr} | {d['time']:.1f} |")
            t = tot.setdefault(m, [0, 0, 0, 0, 0.0])
            t[0] += d["n"]; t[1] += d["correct"]; t[2] += d["wrong"]; t[3] += d["unresolved"]; t[4] += d["time"]
    out += ["", "| method | n | correct | wrong | unresolved | time (s) |", "|---|---|---|---|---|---|"]
    for m, t in tot.items():
        out.append(f"| {m} | {t[0]} | {t[1]} | {t[2]} | {t[3]} | {t[4]:.1f} |")
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    base = os.path.join(os.path.dirname(__file__), "results", "experiment_a")
    with open(base + ".md", "w") as fh: fh.write("\n".join(out) + "\n")
    with open(base + ".json", "w") as fh: json.dump(rows, fh, indent=0)
    print("\n".join(out))


if __name__ == "__main__":
    main()
