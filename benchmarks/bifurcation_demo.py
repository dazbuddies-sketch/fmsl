"""Part B: bifurcations of iterated maps as balance points. Exact Newton-polygon predictions
vs brute-force iteration.  Run:  python benchmarks/bifurcation_demo.py [outdir]"""
import json, math, sys, time
from fractions import Fraction as Fr

import numpy as np

from fmsl.bifurcation import (iterate, divide_exact, shift, branches, period_polynomial,
                              fixed_polynomial, _add)
from fmsl.algebraic import AlgField, A, isolate_real_roots

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
logistic = {(1, 1): Fr(1), (2, 1): Fr(-1)}

# ---------- exact algebraic constants for the second period doubling ----------
p = [Fr(1), 0, Fr(-10), 0, Fr(1)]                      # alpha = sqrt2 + sqrt3
iv = isolate_real_roots(p)[-1]
FLD = AlgField(p, iv[1], iv[2]); al = A.root(FLD)
s6 = (al ** 2 - 5) / 2; s2 = (al ** 3 - 9 * al) / 2
R4 = 1 + s6
X4 = (2 + s6 + s2) / (2 * (1 + s6))

def closed_form(c):
    """k0 + k1 a + k2 a^2 + k3 a^3 with a = sqrt2 + sqrt3  ->  q0 + q2 sqrt2 + q3 sqrt3 + q6 sqrt6"""
    k = (list(c.c) + [Fr(0)] * 4)[:4] if isinstance(c, A) else [Fr(c), 0, 0, 0]
    q0, q2, q3, q6 = k[0] + 5 * k[2], k[1] + 11 * k[3], k[1] + 9 * k[3], 2 * k[2]
    return f"{q0} + ({q2})*sqrt2 + ({q3})*sqrt3 + ({q6})*sqrt6"

CASES = [
    # name, polynomial of the orbit condition, (x0, r0), map (float), r0 float, period, which side
    ("transcritical (logistic, r0 = 1)", fixed_polynomial(logistic), (0, 1),
     lambda x, r: r * x * (1 - x), 1.0, 1),
    ("saddle-node (x + r - x^2, r0 = 0)", fixed_polynomial({(1, 0): Fr(1), (0, 1): Fr(1), (2, 0): Fr(-1)}), (0, 0),
     lambda x, r: x + r - x * x, 0.0, 1),
    ("pitchfork (x + r x - x^3, r0 = 0)", fixed_polynomial({(1, 0): Fr(1), (1, 1): Fr(1), (3, 0): Fr(-1)}), (0, 0),
     lambda x, r: x + r * x - x ** 3, 0.0, 1),
    ("period doubling 1 (logistic, r0 = 3)", period_polynomial(logistic, 2), (Fr(2, 3), 3),
     lambda x, r: r * x * (1 - x), 3.0, 2),
    ("period doubling 2 (logistic, r0 = 1+sqrt6)",
     divide_exact(_add(iterate(logistic, 4), {(1, 0): Fr(-1)}), _add(iterate(logistic, 2), {(1, 0): Fr(-1)})),
     (X4, R4), lambda x, r: r * x * (1 - x), float(R4), 4),
]

def measure(fmap, r, x0, period, start_offset=0.01):
    """iterate until the period-`period` orbit settles; return |u| (distance of the attractor from x0)"""
    x = x0 + start_offset
    prev = None
    for k in range(1, 400):
        for _ in range(20000):
            x = fmap(x, r)
        orbit = [x]
        for _ in range(period - 1):
            orbit.append(fmap(orbit[-1], r))
        near = sorted(orbit, key=lambda t: abs(t - x0))
        if period == 1:
            val = abs(near[0] - x0)
        else:                                          # the two orbit points that split off x0
            val = abs(near[0] - near[1]) / 2
        if prev is not None and abs(val - prev) <= 1e-13 * max(val, 1e-300):
            return val, k * 20000
        prev = val
    return val, 400 * 20000

results = {"cases": []}
eps_list = [1e-2, 1e-3, 1e-4, 1e-5]
print(f"{'bifurcation':44s} {'exact branch':34s} {'e':>7s} {'predicted |u|':>14s} {'measured |u|':>14s} {'ratio':>8s} {'iterations':>10s}")
for name, P, (x0, r0), fmap, r0f, period in CASES:
    t0 = time.time()
    Pu = shift(P, x0, r0)
    brs = [b for b in branches(Pu) if "a" in b and b["c"] > 0]      # branches existing for e > 0
    b = brs[0]
    a, c = b["a"], b["c"]
    cf = float(c)
    exact_txt = f"e = {c if not isinstance(c, A) else f'{cf:.6g}'}*|u|^{a}"
    rows = []
    for e in eps_list:
        pred = (e / cf) ** (1 / float(a))
        meas, iters = measure(fmap, r0f + e, float(x0), period)
        rows.append({"e": e, "pred": pred, "meas": meas, "iters": iters})
        print(f"{name:44s} {exact_txt:34s} {e:7.0e} {pred:14.6e} {meas:14.6e} {meas / pred:8.5f} {iters:10d}")
    slopes = [math.log(rows[i + 1]["meas"] / rows[i]["meas"]) / math.log(eps_list[i + 1] / eps_list[i])
              for i in range(len(rows) - 1)]
    print(f"{'':44s} measured exponent (log-log slope): " + ", ".join(f"{s:.4f}" for s in slopes)
          + f"   exact: {1 / float(a):.4f}   ({time.time() - t0:.0f}s)")
    results["cases"].append({"name": name, "a": str(a), "c": cf,
                             "c_exact": closed_form(c) if isinstance(c, A) else str(c),
                             "rows": rows, "slopes": slopes})
print("\nsecond period-doubling coefficient, exact:", results["cases"][-1]["c_exact"],
      f"= {results['cases'][-1]['c']:.12g}")

# ---------- data for the figure: logistic bifurcation diagram ----------
rs = np.linspace(2.85, 3.57, 1400)
x = np.full_like(rs, 0.5)
for _ in range(4000): x = rs * x * (1 - x)
pts_r, pts_x = [], []
for _ in range(64):
    x = rs * x * (1 - x); pts_r.append(rs.copy()); pts_x.append(x.copy())
results["diagram"] = {"r": np.concatenate(pts_r).tolist(), "x": np.concatenate(pts_x).tolist()}
results["x4"] = float(X4); results["r4"] = float(R4)
Q2s = shift(period_polynomial(logistic, 2), Fr(2, 3), 3)
results["newton_Q2"] = [[i, j, float(v)] for (i, j), v in Q2s.items()]
with open(f"{OUT}/bifurcation_results.json", "w") as fh:
    json.dump(results, fh)
