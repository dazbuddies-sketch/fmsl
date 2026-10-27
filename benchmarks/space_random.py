"""How often the three-variable engine decides a random rational function of x, y, z.

200 random functions N/D (seed 7), N with 1-3 terms and D with 2-3 terms, exponents 0..3,
coefficients in {-3..3}.  Counts verdicts (LIMIT / NO LIMIT / DIVERGES / UNPROVEN) and times.
Run:  python benchmarks/space_random.py [n]
"""
import random, statistics, sys, time
sys.path.insert(0, "src")
from fmsl.space import profile3_polys, parse3

def rp(gen, n):
    P = {}
    for _ in range(n):
        m = (gen.randint(0, 3), gen.randint(0, 3), gen.randint(0, 3))
        P[m] = P.get(m, 0) + gen.choice([-3, -2, -1, 1, 2, 3])
    return " + ".join(f"({v})*x^{a}*y^{b}*z^{c}" for (a, b, c), v in P.items() if v) or "1"

def functions(n=200, seed=7):
    gen = random.Random(seed)
    out = []
    while len(out) < n:
        f = f"({rp(gen, gen.randint(1, 3))})/({rp(gen, gen.randint(2, 3))})"
        try:
            parse3(f)
        except ZeroDivisionError:
            continue
        out.append(f)
    return out

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    counts, times, unproven = {}, [], []
    for f in functions(n):
        N, D = parse3(f)
        t = time.time(); r = profile3_polys(N, D); times.append(time.time() - t)
        v = r["verdict"].split(" (")[0].split(" =")[0]
        v = "DIVERGES" if v.startswith("DIVERGES") else v
        counts[v] = counts.get(v, 0) + 1
        if v == "UNPROVEN":
            unproven.append((f, sorted({e["why"] for e in r["entries"] if e["kind"] == "unresolved"})))
    for k in sorted(counts): print(f"{k:10s} {counts[k]}")
    print(f"decided {n - counts.get('UNPROVEN', 0)}/{n}; median time {statistics.median(times):.2f}s, max {max(times):.1f}s")
    for f, why in unproven: print("  UNPROVEN:", f, why)
