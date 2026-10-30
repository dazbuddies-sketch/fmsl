"""Check fmsl.witness.inner_bound.

Part A (three variables): the inner bound must lie inside the exact answer of FMSL's
three-variable engine (soundness), and we count how often it equals it (completeness).
Part B (four variables, where no exact answer exists): every witnessed value is re-checked
independently by evaluating f exactly (rational arithmetic) along its witness curve
x = c t^w at t = 10^-20, 10^-40, 10^-80, which must approach the value.
Usage: python benchmarks/witness_check.py [nA] [nB]
"""
import math, os, random, signal, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src")); sys.path.insert(0, HERE)
from fmsl.verify import limit_at_origin
from fmsl.witness import inner_bound, parse_n, _ev, curve_point
from part4_test import rand, arcs


class Slow(Exception): pass


def exact3(f, cap=40):
    def h(*_): raise Slow()
    signal.signal(signal.SIGALRM, h); signal.alarm(cap)
    try:
        r = limit_at_origin(f, ["x", "y", "z"]); return r if r["proven"] else None
    except (Slow, Exception): return None
    finally: signal.alarm(0)


def inside(inner, outer, tol=1e-6):
    A, B = arcs(inner, big=1e12), arcs(outer, big=1e12)
    return all(any(lo - tol <= a and b <= hi + tol for lo, hi in B) for a, b in A)


def rand4(gen):
    c = lambda: gen.choice([-2, -1, 1, 2, 3])
    V = ["x", "y", "z", "w"]
    def mono(d):
        e = [0] * 4
        for _ in range(d): e[gen.randrange(4)] += 1
        return "*".join(f"{v}^{k}" for v, k in zip(V, e) if k) or "1"
    num = "+".join(f"({c()})*{mono(gen.randint(2, 4))}" for _ in range(gen.randint(1, 3)))
    den = gen.choice(["x^2+y^2+z^2+w^2", "x^2+y^4+z^2+w^6", "x^2+y^2-z*w", "x*y+z*w",
                      "x^4+y^2+z^2+w^2", "x^2+2*y^2+z^4+w^4+x*y"]) + \
          gen.choice(["", f"+({c()})*{mono(3)}"])
    return f"({num})/({den})"


def main(nA=40, nB=40):
    gen = random.Random(5)
    t = {"inside": 0, "OUTSIDE": 0, "equal": 0, "equal up to tiny gaps": 0, "engine unfinished": 0}
    for i in range(nA):
        f = rand(gen, i % 2 == 1)
        r = exact3(f)
        if r is None: t["engine unfinished"] += 1; continue
        ib = inner_bound(f, ["x", "y", "z"])
        ok = inside(ib["components"], [tuple(c) for c in r["components"]])
        eq = ok and inside([tuple(c) for c in r["components"]], ib["components"], tol=1e-3)
        t["inside" if ok else "OUTSIDE"] += 1
        tiny = len(ib["gaps"]) >= 2 and all(b - a < 1e-6 for a, b in ib["gaps"][:-1] + ib["gaps"][-1:] if b - a < 1)
        t["equal up to tiny gaps" if eq and len(ib["gaps"]) >= 2 else "equal"] += eq
        print("A", "ok " if ok else "XX ", "=" if eq else " ", f, "| exact", r["set"], "| inner", ib["set"], flush=True)
    print("A:", t, flush=True)
    names = ["x", "y", "z", "w"]
    u = {"witnesses checked": 0, "BAD": 0, "whole circle": 0, "limit does not exist": 0}
    for i in range(nB):
        f = rand4(gen)
        ib = inner_bound(f, names)
        P, Q = parse_n(f, names)
        u["whole circle"] += ib["whole_circle"]; u["limit does not exist"] += ib["limit_exists"] is False
        for wt in ib["witnesses"]:
            val = wt[0]
            if val == "segment" or val in ("inf", "-inf"): continue
            v = Fr(val) if "/" in val or val.lstrip("-").isdigit() else Fr(val)
            errs = []
            for k in (20, 40, 80):
                tt = Fr(1, 10 ** k)
                x = curve_point(wt, tt)
                q = _ev(Q, x)
                errs.append(abs(_ev(P, x) / q - v) if q else Fr(10 ** 9))
            good = errs[-1] < errs[0] or errs[-1] == 0
            good = good and errs[-1] < Fr(1, 10 ** 4)
            u["witnesses checked"] += 1; u["BAD"] += not good
            if not good: print("B BAD", f, val, w, c, [float(e) for e in errs], flush=True)
        print("B", f, "->", ib["set"], ib["proven"], flush=True)
    print("B:", u, flush=True)


if __name__ == "__main__":
    a = [int(x) for x in sys.argv[1:3]]
    main(*a)
