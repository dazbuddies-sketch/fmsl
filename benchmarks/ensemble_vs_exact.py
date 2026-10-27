"""Part A: repeated random sampling (Monte Carlo) vs exact analysis of two-variable fluid states.
Run:  python benchmarks/ensemble_vs_exact.py [outdir]"""
import json, random, sys
from fractions import Fraction as Fr

import numpy as np

from fmsl.multivar import profile

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
RNG = np.random.default_rng(2026)


def poly_eval(P, X, Y):
    return sum(float(v) * X ** i * Y ** j for (i, j), v in P.items())


def sample_paths(K, rng=RNG, straight_only=False):
    """random approach paths y = c * |x|^a, x = sigma * s: curvature a log-uniform on [1/4, 4]"""
    sigma = rng.choice([-1.0, 1.0], K)
    a = np.ones(K) if straight_only else np.exp(rng.uniform(np.log(0.25), np.log(4.0), K))
    c = rng.normal(0, 2.0, K)
    return sigma, a, c


def ensemble_values(N, D, paths, s):
    sigma, a, c = paths
    x = sigma * s
    y = c * np.abs(x) ** a
    with np.errstate(all="ignore"):
        return poly_eval(N, x, y) / poly_eval(D, x, y)


def ensemble_verdict(N, D, K=1000, s=1e-9, straight_only=False, rng=RNG):
    v = ensemble_values(N, D, sample_paths(K, rng, straight_only), s)
    finite = np.isfinite(v) & (np.abs(v) < 1e6)
    if not finite.any():
        return "DIVERGES", None, v
    if not finite.all():
        return "NO LIMIT", None, v
    if v.max() - v.min() < 1e-3:
        return "LIMIT", float(np.median(v)), v
    return "NO LIMIT", None, v


def exact_class(r):
    v = r["verdict"][0]
    if v.startswith("LIMIT"): return "LIMIT", Fr(v.split("=")[1].split("(")[0].strip())
    if v.startswith("NO LIMIT"): return "NO LIMIT", None
    if v.startswith("DIVERGES"): return "DIVERGES", None
    return "UNPROVEN", None


def agree(ex, en):
    if ex[0] != en[0]: return False
    if ex[0] == "LIMIT": return abs(float(ex[1]) - en[1]) < 1e-2
    return True


def rpoly(rng, nterms):
    P = {}
    for _ in range(nterms):
        i, j = rng.randint(0, 6), rng.randint(0, 6)
        P[(i, j)] = P.get((i, j), 0) + Fr(rng.choice([-3, -2, -1, 1, 2, 3]))
    return " + ".join(f"({v})*x^{i}*y^{j}" for (i, j), v in P.items() if v) or "1"


if __name__ == "__main__":
    results = {}
    # ---------- 1. head-to-head on random functions ----------
    rng = random.Random(7)
    conf = {}
    rows = []
    while len(rows) < 300:
        f = f"({rpoly(rng, rng.randint(1, 3))})/({rpoly(rng, rng.randint(2, 3))})"
        try:
            r = profile(f)
        except ZeroDivisionError:
            continue
        ex = exact_class(r)
        if ex[0] == "UNPROVEN": continue
        en = ensemble_verdict(r["N"], r["D"])
        conf[(ex[0], en[0])] = conf.get((ex[0], en[0]), 0) + 1
        rows.append((f, ex[0], en[0], agree(ex, en)))
    n_agree = sum(1 for r in rows if r[3])
    print(f"1. Random functions with a proven exact verdict: {len(rows)}")
    print(f"   ensemble (1000 random paths, s = 1e-9) agrees with exact: {n_agree}/{len(rows)}")
    for k in sorted(conf): print(f"   exact {k[0]:9s} -> ensemble {k[1]:9s}: {conf[k]}")
    wrong = [r for r in rows if not r[3]]
    results["random"] = {"n": len(rows), "agree": n_agree,
                         "confusion": {f"{a}|{b}": v for (a, b), v in conf.items()},
                         "wrong_examples": [w[:3] for w in wrong[:8]]}

    # ---------- 2. the named traps ----------
    traps = {
        "parabola trap  x^2*y/(x^4+y^2)": "x^2*y/(x^4+y^2)",
        "two-level trap  x^4/((x-y)^2+x^4)": "x^4/(x^2 + y^2 - 2*x*y + x^4)",
        "thin-curve trap  x(y-x^3)^2/((y-x^3)^2+x^7)": "(y-x^3)^2*x/((y-x^3)^2 + x^7)",
        "sqrt(2) trap  x^4/((y^2-2x^2)^2+x^6)": "x^4/((y^2-2*x^2)^2 + x^6)",
        "honest limit  x^3/(x^2+y^2)": "x^3/(x^2+y^2)",
    }
    print("\n2. Named traps (exact verdict vs ensembles)")
    results["traps"] = {}
    for name, f in traps.items():
        r = profile(f)
        lines = ensemble_verdict(r["N"], r["D"], straight_only=True, rng=np.random.default_rng(1))
        curved = ensemble_verdict(r["N"], r["D"], rng=np.random.default_rng(1))
        big = ensemble_verdict(r["N"], r["D"], K=100000, rng=np.random.default_rng(1))
        print(f"   {name:46s} exact: {r['verdict'][0]:28s} | lines x1000: {lines[0]:9s}"
              f" | curves x1000: {curved[0]:9s} | curves x100000: {big[0]}")
        results["traps"][name] = {"exact": r["verdict"][0], "lines": lines[0], "curves_1e3": curved[0],
                                  "curves_1e5": big[0]}

    # ---------- 3. detection probability vs number of repetitions ----------
    print("\n3. Chance that K random paths reveal the trap (500 independent trials per K)")
    Ks = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000]
    results["detect"] = {}
    for name in list(traps)[:4]:
        r = profile(traps[name])
        probs = []
        for K in Ks:
            hits = 0
            rr = np.random.default_rng(K)
            for _ in range(500):
                v, _, _ = ensemble_verdict(r["N"], r["D"], K=K, rng=rr)
                hits += v != "LIMIT"
            probs.append(hits / 500)
        results["detect"][name] = probs
        print(f"   {name:46s} " + "  ".join(f"K={K}:{p:.2f}" for K, p in zip(Ks, probs)))
    results["Ks"] = Ks

    # ---------- 4. what the repeated samples show for the parabola trap ----------
    r = profile(traps["parabola trap  x^2*y/(x^4+y^2)"])
    v = ensemble_values(r["N"], r["D"], sample_paths(20000, np.random.default_rng(5)), 1e-9)
    results["parabola_samples"] = [float(t) for t in v if np.isfinite(t)]
    print(f"\n4. parabola trap, 20000 samples: values span [{np.nanmin(v):.3f}, {np.nanmax(v):.3f}]"
          f" (exact profile range [-1/2, 1/2]); share of samples with |f| > 0.01: {np.mean(np.abs(v) > 0.01):.3f}")
    with open(f"{OUT}/ensemble_results.json", "w") as fh:
        json.dump(results, fh)
