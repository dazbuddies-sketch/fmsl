"""LimitBench-Beyond: limits beyond fractions of polynomials.

Five families, matching what FMSL 0.20 to 0.25 added:
  B1  sin, cos, exp, log, ... in two variables      e.g. (1 - cos(x*y))/(x^2+y^2)^2
  B2  a common factor along a line                  e.g. sin(x+y)*x/((x+y)*(x^2+y^2))
  B3  distances, plain and weighted, 2 and 3 vars   e.g. y/sqrt(x^2+2*y^2)
  B4  at infinity                                   e.g. (x^2+y)/(x^2+y^2) as |(x, y)| -> oo
  B5  from one side (half-planes, quadrants)        e.g. x*y/(x^2+y^2) with x > 0, y > 0
Ground truth: FMSL's proven answer, kept only when an independent check agrees: direct 320-digit
evaluation of the original formula at about 12,000 points (near the point, beyond 1e36, or with
the required signs) never leaves the proven set and fills it (benchmarks/analytic_crosscheck.py
and region_crosscheck.py). Items FMSL cannot prove, or the check cannot confirm, are dropped.

Run:  python benchmarks/limitbench/generate_beyond.py [n_per_family]  -> limitbench_beyond.jsonl
"""
import json, math, os, random, signal, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..")); sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
from fmsl.verify import limit_at, limit_at_origin
import analytic_crosscheck as ac
import region_crosscheck as rc


class _Slow(Exception): pass


def _alarm(*_): raise _Slow()


FUN = ["sin({})", "(exp({})-1)", "(1-cos({}))", "log(1+{})", "atan({})", "tan({})", "sinh({})"]
DEN = ["x^2+y^2", "x^4+y^2", "x^2+y^4", "(x^2+y^2)^2", "x^2+x*y+y^2", "x^6+y^2", "x^2+2*y^2"]
MONO = ["x", "y", "x*y", "x^2", "y^2", "x^2*y", "x*y^2", "x^3", "y^3"]


WDEN = [("x^2+y^2", 1, 1, 2), ("x^4+y^2", 1, 2, 4), ("x^2+y^4", 2, 1, 4), ("(x^2+y^2)^2", 1, 1, 4),
        ("x^2+x*y+y^2", 1, 1, 2), ("x^6+y^2", 1, 3, 6), ("x^2+2*y^2", 1, 1, 2), ("x^4+x^2*y+y^2", 1, 2, 4)]


def _monos(wx, wy, d):
    return [f"x^{i}*y^{j}" for i in range(0, 7) for j in range(0, 7) if i * wx + j * wy == d and i + j > 0]


def b1(gen):
    """numerator of the same weighted degree as the denominator, so the answer is rarely just 0"""
    D, wx, wy, d = gen.choice(WDEN)
    for _ in range(50):
        fn = gen.choice(FUN); k = 2 if fn.startswith("(1-cos") else 1
        a = gen.randint(1, d // k)
        inner, outer = _monos(wx, wy, a), _monos(wx, wy, d - k * a) if d - k * a > 0 else ["1"]
        if inner and outer: break
    num = fn.format(gen.choice(inner))
    o = gen.choice(outer)
    if o != "1": num += "*" + o
    if gen.random() < .3: num += f"+({gen.choice([-1, 1, 2])})*" + gen.choice(_monos(wx, wy, d) or ["x*y"])
    return f"({num})/({D})", ["x", "y"], None, None


def b2(gen):
    a, b = gen.choice([(1, 1), (1, -1), (1, -2), (2, 1), (3, -1)])
    l = f"({a}*x+({b})*y)"
    fn = gen.choice(["sin({})", "(exp({})-1)", "log(1+{})", "atan({})", "tan({})"])
    g = gen.choice(["1", "x", "y", "x^2", "x*y"])
    h = gen.choice(["1", "x^2+y^2", "x^4+y^2", "x^2+x*y+y^2"])
    inner = gen.choice(["", "*x", "*y"])
    return f"{fn.format(l + inner)}*{g}/({l}*({h}))", ["x", "y"], None, None


def b3(gen):
    if gen.random() < .7:
        Q = gen.choice(["x^2+y^2", "x^2+2*y^2", "x^2+x*y+y^2", "4*x^2+y^2", "x^2+3*y^2"])
        R = f"sqrt({Q})"
        f = gen.choice([f"x/{R}", f"y/{R}", f"x*y/{R}^2", f"({R}-x)/{R}", f"x^2*y/{R}^3", f"({R}-x)/y",
                        f"x*y/({R}*({R}+x))", f"sin(x*y)/{R}^2", f"(x^2-y^2)/({R}*({R}+y))"])
        return f, ["x", "y"], None, None
    Q = gen.choice(["x^2+y^2+z^2", "x^2+2*y^2+z^2", "x^2+y^2+2*z^2"])
    R = f"sqrt({Q})"
    f = gen.choice([f"x/{R}", f"z/{R}", f"x*y/{R}^2", f"(x*y+y*z)/{R}^2", f"({R}-z)/{R}"])
    return f, ["x", "y", "z"], None, None


def b4(gen):
    def mono(): return f"x^{gen.randint(0, 3)}*y^{gen.randint(0, 3)}"
    num = "+".join(f"({gen.choice([-2, -1, 1, 2])})*{mono()}" for _ in range(gen.randint(1, 2)))
    den = gen.choice(["x^2+y^2", "x^4+y^2", "x^2+y^4", "x^2+y^2+1", "x^4+y^4", "x^6+y^2", "x^2+x*y+y^2"])
    return f"({num})/({den})", ["x", "y"], "oo", None


def b5(gen):
    f = b1(gen)[0] if gen.random() < .4 else gen.choice([
        "x*y/(x^2+y^2)", "y/x", "(x-y)/(x+y)", "x^2*y/(x^4+y^2)", "x*y^2/(x^2+y^4)", "(x^2-y^2)/(x^2+y^2)",
        "x*y/(x^2+x*y+y^2)", "y^2/(x^2+y^2)", "(x+y)/(x-y)", "x^3*y/(x^6+y^2)"])
    sd = gen.choice([{"x": "+"}, {"y": "-"}, {"x": "+", "y": "+"}, {"x": "+", "y": "-"}, {"x": "-", "y": "+"}])
    return f, ["x", "y"], None, sd


FAMILIES = [("beyond: sin, exp, log", b1, 3), ("beyond: common factor", b2, 4), ("beyond: distance", b3, 4),
            ("beyond: at infinity", b4, 5), ("beyond: one side", b5, 5)]


def confirm(f, names, point, side, gen):
    """FMSL's proven answer, if the independent direct evaluation agrees; else None"""
    if point == "oo":
        r = limit_at(f, "oo", names)
    elif side:
        r = limit_at(f, [0] * len(names), names, side=side)
    else:
        r = limit_at_origin(f, names)
    if not r["proven"] or r["components"] is None: return None
    pts = ac.points(len(names), gen, 12000)
    if point == "oo": pts = rc.far(pts)
    elif side: pts = rc.signed(pts, names, side)
    vals = rc.values(f, names, pts)
    if rc.verdict(vals, r["components"]) != "ok": return None
    return r


def main(n=40, only=None):
    """only: a family index (0-4) to generate just that family into limitbench_beyond.part<k>.jsonl,
       so the five families can run in parallel; then `merge` joins the parts"""
    gen = random.Random(2026 + (only or 0))
    items, seen = [], set()
    for k, (cat, make, tier) in enumerate(FAMILIES):
        if only is not None and k != only: continue
        got, tries = 0, 0
        while got < n and tries < 25 * n:
            tries += 1
            f, names, point, side = make(gen)
            key = (f, point, json.dumps(side, sort_keys=True))
            if key in seen: continue
            seen.add(key)
            signal.signal(signal.SIGALRM, _alarm); signal.alarm(60)
            try:
                r = confirm(f, names, point, side, gen)
            except Exception:
                r = None
            finally:
                signal.alarm(0)
            if r is None: continue
            if r["limit_exists"] and gen.random() > 0.3: continue      # keep "the limit exists" a minority
            comps = [[a, b] for a, b in r["components"]]
            items.append({"f": f, "variables": names, "point": point or "origin", "side": side,
                          "category": cat, "tier": tier, "limit_exists": r["limit_exists"], "limit": r["limit"],
                          "fluid_state": r["set"], "fluid_state_float": comps, "path_test": None,
                          "verified_by": "FMSL proven + direct 320-digit evaluation (fills the set, never leaves it)",
                          "id": f"LBB-{len(items) + 1:04d}"})
            got += 1
            print(cat, got, f, point or "", side or "", r["set"], flush=True)
    out = "limitbench_beyond.jsonl" if only is None else f"limitbench_beyond.part{only}.jsonl"
    with open(os.path.join(HERE, out), "w") as fh:
        for it in items: fh.write(json.dumps(it) + "\n")
    from collections import Counter
    print(len(items), "items;", dict(Counter(it["category"] for it in items)),
          "; limit exists in", sum(it["limit_exists"] for it in items))


def merge():
    items = []
    for k in range(len(FAMILIES)):
        p = os.path.join(HERE, f"limitbench_beyond.part{k}.jsonl")
        items += [json.loads(l) for l in open(p)]
        os.remove(p)
    for i, it in enumerate(items): it["id"] = f"LBB-{i + 1:04d}"
    with open(os.path.join(HERE, "limitbench_beyond.jsonl"), "w") as fh:
        for it in items: fh.write(json.dumps(it) + "\n")
    print(len(items), "items")


if __name__ == "__main__":
    if sys.argv[1:2] == ["merge"]: merge()
    else: main(int(sys.argv[1]) if len(sys.argv) > 1 else 40, int(sys.argv[2]) if len(sys.argv) > 2 else None)
