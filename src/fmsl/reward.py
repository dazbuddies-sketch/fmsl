"""FMSL as a reward function and problem source for training and evaluating AI models.

Training a model to reason with verifiable rewards needs two things: an endless supply of problems,
and an automatic grader that cannot be fooled. Limit problems in several variables suit this well:
the answer is a set (an interval, a union of pieces, a single value), so it is hard to guess, and
FMSL proves it.

    from fmsl.reward import problems, score
    for p in problems(5, seed=1):
        print(p["prompt"])
        print(score(p, "ANSWER: SET [-1/2, 1/2]"))

score() returns
    reward   1.0 for an exactly right answer; otherwise 0.8 x the overlap (intersection over union
             of claimed and true sets, measured on the circle [-oo, +oo] through atan, where a
             single value counts as an arc of width 0.01), so partly right answers earn partial
             credit and an exact answer is always strictly best; 0 if unparsed
    exact    the set (or the limit) is right, endpoints within 1e-3 relative
    verdict  right about whether the limit exists
    missing / extra   what the claim leaves out / includes wrongly (from fmsl.verify.check_claim)
Problems FMSL cannot prove are never produced, so every reward rests on a proof.
"""
import math, random, re

__all__ = ["score", "problems", "FAMILIES"]

_W = 0.01                                               # width given to a single value, in atan units


def _arcs(comps):
    a = lambda v: math.copysign(math.pi / 2, v) if math.isinf(v) else math.atan(v)
    out = []
    for lo, hi in comps:
        A, B = a(lo), a(hi)
        if B - A < _W: m = (A + B) / 2; A, B = m - _W / 2, m + _W / 2
        out.append((A, B))
    out.sort()
    merged = []
    for A, B in out:
        if merged and A <= merged[-1][1]: merged[-1] = (merged[-1][0], max(merged[-1][1], B))
        else: merged.append((A, B))
    return merged


def _measure(arcs): return sum(B - A for A, B in arcs)


def _intersect(P, Q):
    out = []
    for a, b in P:
        for c, d in Q:
            lo, hi = max(a, c), min(b, d)
            if hi > lo: out.append((lo, hi))
    return out


def _claim_text(answer):
    """the claim from a model's reply: the last 'ANSWER: LIMIT x' / 'ANSWER: SET s' line, or the text"""
    m = re.findall(r"ANSWER:\s*(LIMIT|SET)\s*(.+)", answer)
    if not m: return answer.strip()
    kind, rest = m[-1]
    return rest.strip() if kind == "SET" else "LIMIT " + rest.strip()


def score(problem, answer):
    """grade a reply against a problem from problems() (or any dict with f, variables, point, side)"""
    from .verify import check_claim, parse_claim
    claim = _claim_text(answer)
    point = problem.get("point")
    point = None if point in (None, "origin") else point
    try:
        C = parse_claim(claim)
    except Exception:
        return {"reward": 0.0, "exact": False, "verdict": False, "parsed": False}
    r = check_claim(problem["f"], claim, problem.get("variables"), point=point, side=problem.get("side"))
    if r["correct"] is None:
        return {"reward": None, "exact": None, "verdict": None, "parsed": True, "note": r["verdict"]}
    truth = problem.get("_components") or _components(problem)
    T, A = _arcs(truth), _arcs([tuple(c) for c in C])
    inter = _measure(_intersect(T, A))
    union = _measure(T) + _measure(A) - inter
    iou = inter / union if union > 0 else 0.0
    says_limit = len(C) == 1 and C[0][0] == C[0][1] and not math.isinf(C[0][0])
    exact = bool(r["correct"])
    return {"reward": 1.0 if exact else round(0.8 * iou, 4), "exact": exact,
            "verdict": says_limit == bool(r["limit_exists"]), "parsed": True,
            "missing": r["missing"], "extra": r["extra"], "truth": r["truth"]}


def _components(problem):
    from .verify import limit_at, limit_at_origin
    vs, p, sd = problem.get("variables"), problem.get("point"), problem.get("side")
    if p == "oo": r = limit_at(problem["f"], "oo", vs)
    elif sd: r = limit_at(problem["f"], [0] * len(vs), vs, side=sd)
    else: r = limit_at_origin(problem["f"], vs)
    return [tuple(c) for c in r["components"]]


# ---------------------------------------------------------------- problems
_WDEN = [("x^2+y^2", 1, 1, 2), ("x^4+y^2", 1, 2, 4), ("x^2+y^4", 2, 1, 4), ("(x^2+y^2)^2", 1, 1, 4),
         ("x^2+x*y+y^2", 1, 1, 2), ("x^6+y^2", 1, 3, 6), ("x^2+2*y^2", 1, 1, 2), ("x^4+x^2*y+y^2", 1, 2, 4)]
_FUN = ["sin({})", "(exp({})-1)", "(1-cos({}))", "log(1+{})", "atan({})", "tan({})", "sinh({})"]


def _monos(wx, wy, d):
    return [f"x^{i}*y^{j}" for i in range(0, 7) for j in range(0, 7) if i * wx + j * wy == d and i + j > 0]


def _num(gen, wx, wy, d, functions):
    terms = []
    for _ in range(gen.randint(1, 2)):
        c = gen.choice([-2, -1, 1, 2, 3])
        if functions and gen.random() < .7:
            fn = gen.choice(_FUN); k = 2 if fn.startswith("(1-cos") else 1
            opts = [(a, _monos(wx, wy, a), _monos(wx, wy, d - k * a) if d - k * a > 0 else ["1"]) for a in range(1, d // k + 1)]
            opts = [o for o in opts if o[1] and o[2]]
            if opts:
                a, inner, outer = gen.choice(opts)
                t = fn.format(gen.choice(inner)); o = gen.choice(outer)
                terms.append(f"({c})*{t}" + ("" if o == "1" else "*" + o)); continue
        terms.append(f"({c})*{gen.choice(_monos(wx, wy, d) or ['x*y'])}")
    return "+".join(terms)


def _rational(gen):
    D, wx, wy, d = gen.choice(_WDEN)
    return {"f": f"({_num(gen, wx, wy, d, False)})/({D})", "variables": ["x", "y"]}


def _analytic(gen):
    D, wx, wy, d = gen.choice(_WDEN)
    return {"f": f"({_num(gen, wx, wy, d, True)})/({D})", "variables": ["x", "y"]}


def _distance(gen):
    Q = gen.choice(["x^2+y^2", "x^2+2*y^2", "x^2+x*y+y^2", "4*x^2+y^2", "x^2+3*y^2"])
    R = f"sqrt({Q})"
    return {"f": gen.choice([f"x/{R}", f"y/{R}", f"({R}-x)/{R}", f"x^2*y/{R}^3", f"({R}-x)/y",
                             f"x*y/({R}*({R}+x))", f"sin(x*y)/{R}^2", f"(x^2-y^2)/({R}*({R}+y))"]),
            "variables": ["x", "y"]}


def _infinity(gen):
    D, wx, wy, d = gen.choice(_WDEN)
    return {"f": f"({_num(gen, wx, wy, d, False)})/({D})", "variables": ["x", "y"], "point": "oo"}


def _one_side(gen):
    p = gen.choice([_rational, _analytic])(gen)
    p["side"] = gen.choice([{"x": "+"}, {"y": "-"}, {"x": "+", "y": "+"}, {"x": "+", "y": "-"}])
    return p


FAMILIES = {"rational": _rational, "analytic": _analytic, "distance": _distance,
            "infinity": _infinity, "one_side": _one_side}

_PROMPT = """Let f({vars}) = {f}.
Find the set of ALL limit values of f {where}: every L (a real number, +oo or -oo) such that f(p_n) -> L
for some sequence of points p_n {seq} at which f is defined{cond}. If this set is a single real
number, the limit exists.
End with exactly one line: ANSWER: LIMIT <value>  or  ANSWER: SET <set>, where <set> is a union of
closed intervals and points such as [-1/2, 1/2] or [-oo, -3] u [3, +oo] or {{0}} u [1, 2]."""


def _prompt(p):
    vs = ", ".join(p["variables"])
    zeros = ", ".join("0" for _ in p["variables"])
    if p.get("point") == "oo":
        where, seq, cond = f"as |({vs})| -> oo", "with |p_n| -> oo", ""
    else:
        where, seq = f"as ({vs}) -> ({zeros})", f"-> ({zeros})"
        cond = ""
        if p.get("side"):
            cond = ", with " + " and ".join(f"{k} {'>' if s == '+' else '<'} 0" for k, s in p["side"].items())
    return _PROMPT.format(vars=vs, f=p["f"], where=where, seq=seq, cond=cond)


def problems(n, families=None, seed=0, limit_share=0.3, timeout=30):
    """n problems with proven answers: dicts with f, variables, point, side, prompt, answer (the
       proven set), limit_exists, family. Families: rational, analytic, distance, infinity, one_side
       (default: all, in turn). At most a share limit_share of them have a limit (a single value)."""
    import signal
    from .verify import limit_at, limit_at_origin
    gen = random.Random(seed)
    fams = list(families or FAMILIES)
    out, seen, tries = [], set(), 0
    def alarm(*_): raise TimeoutError()
    while len(out) < n and tries < 50 * n:
        tries += 1
        fam = fams[len(out) % len(fams)]
        p = FAMILIES[fam](gen)
        key = (p["f"], p.get("point"), str(p.get("side")))
        if key in seen: continue
        seen.add(key)
        try:
            old = signal.signal(signal.SIGALRM, alarm); signal.alarm(timeout)
        except ValueError:                              # not in the main thread: no time limit
            old = None
        try:
            vs = p["variables"]
            if p.get("point") == "oo": r = limit_at(p["f"], "oo", vs)
            elif p.get("side"): r = limit_at(p["f"], [0] * len(vs), vs, side=p["side"])
            else: r = limit_at_origin(p["f"], vs)
        except Exception:
            continue
        finally:
            if old is not None: signal.alarm(0)
        if not r["proven"] or r["components"] is None: continue
        if r["limit_exists"] and gen.random() > limit_share: continue
        p.update(family=fam, answer=r["set"], limit_exists=r["limit_exists"], limit=r["limit"],
                 _components=[tuple(c) for c in r["components"]])
        p.setdefault("point", "origin"); p.setdefault("side", None)
        p["prompt"] = _prompt(p)
        out.append(p)
    return out
