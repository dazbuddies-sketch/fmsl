"""LimitBench: limits of functions of several variables, with proven answers.

Each item asks for the complete set of values f approaches as the point tends to the origin
(the fluid state): whether the limit exists, its value if it does, and otherwise the closed set
of values.  Ground truth comes from two independent sources that must agree:
  * construction (known closed forms: AM-GM bounds, eigenvalues of quadratic forms, ...), or
  * FMSL's exact engine AND the independent exact loop-image check (benchmarks/loop_crosscheck.py);
items where they disagree, or where FMSL is not complete, are dropped.
Each item also records what the common "try a few paths" method concludes, and whether that is a
TRAP (the method's verdict is wrong).

Run:  python benchmarks/limitbench/generate.py      -> benchmarks/limitbench/limitbench.jsonl
"""
import json, math, os, random, sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src")); sys.path.insert(0, os.path.join(HERE, ".."))
from fmsl.multivar import profile, parse_rational
from fmsl.space import profile3
from fmsl.cluster import _Inf

PATHS2 = [((1, 1), None), (None, (1, 1)), ((1, 1), (1, 1)), ((1, 1), (-1, 1)), ((1, 1), (2, 1)),
          ((1, 1), (1, 2)), ((1, 1), (-1, 2)), ((1, 2), (1, 1))]


def _lim(TN, TD):
    def lead(T):
        T = {e: c for e, c in T.items() if c != 0}
        return (min(T), T[min(T)]) if T else None
    a, b = lead(TN), lead(TD)
    if b is None: return None
    if a is None: return Fr(0)
    if a[0] > b[0]: return Fr(0)
    if a[0] == b[0]: return Fr(a[1]) / b[1]
    return "+oo" if a[1] / b[1] > 0 else "-oo"


def path_values(N, D):
    out = []
    for fx, fy in PATHS2:
        def sub(P):
            T = {}
            for (i, j), c in P.items():
                if (i and fx is None) or (j and fy is None): continue
                cx, ex = fx if i else (1, 0); cy, ey = fy if j else (1, 0)
                e = (ex * i if i else 0) + (ey * j if j else 0)
                T[e] = T.get(e, 0) + c * Fr(cx) ** i * Fr(cy) ** j
            return T
        out.append(_lim(sub(N), sub(D)))
    return out


def f2s(v):
    if isinstance(v, _Inf): return "+oo" if v.sgn > 0 else "-oo"
    if isinstance(v, str): return v
    if isinstance(v, Fr): return str(v.numerator) if v.denominator == 1 else f"{v.numerator}/{v.denominator}"
    if isinstance(v, (int,)): return str(v)
    return repr(float(v)) if not isinstance(v, float) else repr(v)


def as_float(v):
    if isinstance(v, _Inf): return math.inf if v.sgn > 0 else -math.inf
    if v in ("+oo", "-oo"): return math.inf if v == "+oo" else -math.inf
    return float(v)


def item_from_fmsl(f, category, tier, gen_note):
    N, D = parse_rational(f)
    r = profile(f)
    cl = r["cluster"]
    if not (cl.get("exact") and cl.get("complete", True)): return None
    comps = cl["signed_components"]
    # independent check: exact image on small squares must match
    import loop_crosscheck as lc
    im = lc.loop_image(N, D, Fr(1, 10 ** 16))
    if im is None: return None
    if lc.mismatch_line(lc.coverage_line(im), lc.coverage_line(lc.signed_ivs(cl))) >= 0.002: return None
    pv = path_values(N, D)
    if any(v is None for v in pv): return None
    pvs = sorted({f2s(v) for v in pv})
    single = len(comps) == 1 and comps[0][0] == comps[0][1] and not isinstance(comps[0][0], _Inf)
    path_says_limit = len(pvs) == 1 and pvs[0] not in ("+oo", "-oo")
    return {"f": f, "variables": ["x", "y"], "category": category, "tier": tier,
            "limit_exists": single, "limit": f2s(comps[0][0]) if single else None,
            "fluid_state": cl["signed_text"],
            "fluid_state_float": [[as_float(a), as_float(b)] for a, b in comps],
            "path_test": {"values": pvs, "says_limit_exists": path_says_limit,
                          "trap": path_says_limit != single or (single and pvs[0] != f2s(comps[0][0]))},
            "verified_by": "FMSL exact engine + independent exact loop image (r = 1e-16)", "note": gen_note}


def main():
    gen = random.Random(2027)
    items = []
    def add(it):
        if it is not None and not any(o["f"] == it["f"] for o in items): items.append(it)
    # tier 1: classic, by construction and verified
    for p, q in [(1, 1), (1, 2), (2, 1), (2, 3), (3, 1), (1, 3), (3, 2), (2, 2)]:
        for k in (1, 2, 3):
            add(item_from_fmsl(f"{k}*x^{p}*y^{q}/(x^{2*p}+y^{2*q})", "AM-GM ratio", 1, "u v/(u^2+v^2) with u = x^p, v = y^q"))
    for a in range(-2, 3):
        for b in (-2, -1, 1, 2):
            for c in range(-2, 3):
                add(item_from_fmsl(f"(({a})*x^2+({b})*x*y+({c})*y^2)/(x^2+y^2)", "quadratic form", 1, "range = eigenvalues"))
    # tier 1b: the limit EXISTS (so that "no limit" is not a winning default answer)
    for p_ in range(0, 4):
        for q_ in range(0, 4):
            if p_ + q_ < 3: continue
            for k in (1, 2, -1):
                add(item_from_fmsl(f"({k})*x^{p_}*y^{q_}/(x^2+y^2)", "limit exists", 1, "numerator of higher order: limit 0"))
    for c in (-3, -2, -1, 1, 2, 3):
        for extra in ("x^3", "x*y^2", "x^2*y", "y^3", "x^4"):
            add(item_from_fmsl(f"({c}*x^2+{c}*y^2+{extra})/(x^2+y^2)", "limit exists", 1, f"limit {c}"))
            add(item_from_fmsl(f"({c}*x^4+{c}*y^2+{extra}*y)/(x^4+y^2)", "limit exists", 2, f"limit {c} at scale y ~ x^2"))
    for a_ in (1, 2):
        for b_ in (1, 2):
            add(item_from_fmsl(f"x^{2*a_+1}*y^{2*b_}/(x^{2*a_}+y^{2*b_})", "limit exists", 2, "AM-GM bound times a vanishing factor"))
            add(item_from_fmsl(f"(x^2*y^2)/(x^{2*a_}+y^{2*b_})", "limit exists", 2, "limit 0 when the degree wins"))
    # tier 2: path traps (curves y = c x^k that standard paths miss)
    for k in (3, 4, 5):
        for c in (1, 2, 3):
            add(item_from_fmsl(f"{c}*x^{k}*y/(x^{2*k}+y^2)", "hidden curve", 2, f"extremes on y ~ x^{k}"))
            add(item_from_fmsl(f"x^{2*k}/((y-{c}*x^{k})^2+x^{4*k})", "hidden curve", 2, f"value 1 only along y = {c} x^{k}"))
    for k in (2, 3):
        add(item_from_fmsl(f"x^{k}*y^2/(x^{2*k}+y^4)", "hidden curve", 2, "u v/(u^2+v^2) type"))
    # tier 3: random rational functions (FMSL + loop check must agree)
    def rpoly(n):
        P = {}
        for _ in range(n):
            i, j = gen.randint(0, 4), gen.randint(0, 4)
            if i + j == 0: continue
            P[(i, j)] = P.get((i, j), 0) + gen.choice([-3, -2, -1, 1, 2, 3])
        return " + ".join(f"({v})*x^{i}*y^{j}" for (i, j), v in P.items() if v) or "x"
    tries = 0
    while sum(1 for it in items if it["tier"] == 3) < 250 and tries < 2500:
        tries += 1
        try:
            add(item_from_fmsl(f"({rpoly(gen.randint(1, 3))})/({rpoly(gen.randint(2, 3))})", "random rational", 3, ""))
        except (ZeroDivisionError, ValueError, ArithmeticError):
            pass
    # tier 4: three variables, by construction
    for a in range(-2, 3):
        for b in range(-2, 3):
            for c in range(-2, 3):
                if a == b == c: continue
                if gen.random() > 0.35: continue
                f = f"(({a})*x^2+({b})*y^2+({c})*z^2)/(x^2+y^2+z^2)"
                r = profile3(f)
                lo, hi = min(a, b, c), max(a, b, c)
                if r["signed_text"] != f"[{lo}, {hi}]": continue           # construction and engine must agree
                items.append({"f": f, "variables": ["x", "y", "z"], "category": "3-D quadratic form", "tier": 4,
                              "limit_exists": False, "limit": None, "fluid_state": r["signed_text"],
                              "fluid_state_float": [[float(lo), float(hi)]], "path_test": None,
                              "verified_by": "construction (diagonal form: range = [min, max] of coefficients) + FMSL", "note": ""})
    for k in (1, 2, 3):
        for f, st in [(f"{k}*x*y*z/(x^2+y^2+z^2)", "{0}"), (f"{k}*(x*y+y*z+z*x)/(x^2+y^2+z^2)", f"[{f2s(Fr(-k, 2))}, {k}]")]:
            r = profile3(f)
            if r["signed_text"] == st:
                lim = st == "{0}"
                items.append({"f": f, "variables": ["x", "y", "z"], "category": "3-D classic", "tier": 4,
                              "limit_exists": lim, "limit": "0" if lim else None, "fluid_state": st,
                              "fluid_state_float": [[0.0, 0.0]] if lim else [[-k / 2, float(k)]], "path_test": None,
                              "verified_by": "construction + FMSL", "note": ""})
    for n, it in enumerate(items): it["id"] = f"LB-{n + 1:04d}"
    out = os.path.join(HERE, "limitbench.jsonl")
    with open(out, "w") as fh:
        for it in items: fh.write(json.dumps(it) + "\n")
    from collections import Counter
    print(len(items), "items;", dict(Counter(it["tier"] for it in items)), "by tier;",
          sum(1 for it in items if it["path_test"] and it["path_test"]["trap"]), "path-test traps;",
          sum(1 for it in items if it["limit_exists"]), "with a limit")


if __name__ == "__main__":
    main()
