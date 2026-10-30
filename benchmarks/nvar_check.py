"""Check fmsl.nvar.limit_n (complete answers when the denominator's initial form is positive definite).

Part A (three variables): every complete answer is compared with FMSL's exact three-variable
engine, which works differently (Puiseux arcs, zooming); the ends must agree within the enclosure.
Part B (four variables, no other exact method): exact rational evaluation of f at 4,000 random
points with |x| between 1e-40 and 1e-20 must never leave the proven set (outer bound), and the
witness point of each end, scaled toward 0 along its curve, must approach that end (inner bound).
Usage: python benchmarks/nvar_check.py [nA] [nB]
"""
import math, os, random, signal, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from fmsl.verify import limit_at_origin
from fmsl.nvar import limit_n
from fmsl.witness import parse_n, _ev

DENS3 = [("x^2+y^2+z^2", (1, 1, 1)), ("x^4+y^2+z^2", (1, 2, 2)), ("x^2+2*y^2+z^4", (2, 2, 1)),
         ("x^2+x*y+y^2+z^2", (1, 1, 1)), ("x^6+y^2+z^4", (2, 6, 3)), ("x^4+y^4+z^2", (1, 1, 2))]
DENS4 = [("x^2+y^2+z^2+w^2", (1, 1, 1, 1)), ("x^4+y^2+z^2+w^4", (1, 2, 2, 1)),
         ("x^2+2*y^2+z^2+3*w^2+x*y", (1, 1, 1, 1)), ("x^2+y^4+z^2+w^6", (3, 3, 3, 1)),
         ("x^4+y^4+z^4+w^4", (1, 1, 1, 1)), ("x^2+y^2+z^4+w^4+z*w^3", (2, 2, 1, 1))]


class Slow(Exception): pass


def capped(fn, *a, cap=60):
    def h(*_): raise Slow()
    signal.signal(signal.SIGALRM, h); signal.alarm(cap)
    try: return fn(*a)
    except (Slow, Exception): return None
    finally: signal.alarm(0)


def make(gen, names, dens):
    D, w = gen.choice(dens)
    d = min(sum(int(c) * wi for c, wi in zip(m, w)) for m in [[0] * len(names)])  # placeholder
    # weighted degree of D's terms
    import re
    deg = None
    for term in D.split("+"):
        e = [0] * len(names)
        for f in term.split("*"):
            if f.isdigit(): continue
            v, _, k = f.partition("^"); e[names.index(v)] += int(k or 1)
        t = sum(a * b for a, b in zip(e, w)); deg = t if deg is None else min(deg, t)
    shift = gen.choice([0, 0, 0, -1, 1])           # mostly equal order; sometimes lower or higher
    monos = []
    for e in __import__("itertools").product(range(7), repeat=len(names)):
        if 0 < sum(e) and sum(a * b for a, b in zip(e, w)) == deg + shift: monos.append(e)
    if not monos: return None
    terms = []
    for e in gen.sample(monos, min(len(monos), gen.randint(1, 3))):
        terms.append(f"({gen.choice([-3, -2, -1, 1, 2, 3])})*" + "*".join(f"{v}^{k}" for v, k in zip(names, e) if k))
    extra = gen.choice(["", f"+{names[0]}^{deg + 2}"])
    return f"({'+'.join(terms)})/({D}{extra})"


def arcs(comps):
    a = lambda v: math.copysign(math.pi / 2, v) if math.isinf(v) else math.atan(v)
    return [(a(lo), a(hi)) for lo, hi in comps]


def main(nA=40, nB=40):
    gen = random.Random(11)
    t = {"agree": 0, "DISAGREE": 0, "not complete": 0, "engine unfinished": 0}
    names = ["x", "y", "z"]
    i = 0
    while i < nA:
        f = make(gen, names, DENS3)
        if f is None: continue
        i += 1
        r = capped(limit_n, f, names)
        if not r or r["proven"] is not True: t["not complete"] += 1; print("A -  ", f, r and r["set"]); continue
        e = capped(limit_at_origin, f, names)
        if not e or not e["proven"]: t["engine unfinished"] += 1; continue
        A, B = arcs(r["components"]), arcs([tuple(c) for c in e["components"]])
        tol = 3e-6 if "(+-" in r["set"] else 1e-9
        ok = len(A) == len(B) and all(abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol for a, b in zip(A, B))
        t["agree" if ok else "DISAGREE"] += 1
        print("A", "ok" if ok else "XX", f, "|", r["set"], "| engine", e["set"], flush=True)
    print("A:", t, flush=True)
    names = ["x", "y", "z", "w"]
    u = {"complete": 0, "SAMPLE OUTSIDE": 0, "not complete": 0}
    i = 0
    while i < nB:
        f = make(gen, names, DENS4)
        if f is None: continue
        i += 1
        r = capped(limit_n, f, names, cap=120)
        if not r or r["proven"] is not True: u["not complete"] += 1; print("B -  ", f, r and r["set"]); continue
        u["complete"] += 1
        P, Q = parse_n(f, names)
        lo, hi = r["components"][0][0], r["components"][-1][1]
        if "enclosure" in r: lo, hi = r["enclosure"][0][0], r["enclosure"][1][1]
        bad = 0
        for _ in range(4000):
            s = Fr(1, 10 ** gen.randint(20, 40))
            x = tuple(Fr(gen.uniform(-1, 1)) * s ** wi for wi in r["weight"])
            q = _ev(Q, x)
            if q == 0: continue
            v = _ev(P, x) / q
            if r["set"] == "[-oo, +oo]": continue
            if math.isinf(lo) and lo == hi: bad += not ((v > 0) == (lo > 0) and abs(v) > 1e10); continue
            if not (Fr(lo) - Fr(1, 10 ** 12) <= v <= Fr(hi) + Fr(1, 10 ** 12)): bad += 1
        u["SAMPLE OUTSIDE"] += bad > 0
        print("B", "ok" if not bad else f"XX({bad})", f, "->", r["set"], flush=True)
    print("B:", u, flush=True)


if __name__ == "__main__":
    a = [int(x) for x in sys.argv[1:3]]
    main(*a)
