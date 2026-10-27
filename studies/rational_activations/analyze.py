"""The pre-registered analysis (PREREGISTRATION.md) of the runs written by train.py.

For every step and both activations: the EXACT test (FMSL) of whether the activation has a
genuine pole inside the range of its inputs, and the two baselines; then the instability events
and the H1/H2 scores.  Writes results.json.
"""
import glob, json, os, sys
from fractions import Fraction as Fr
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from fmsl.algebraic import pgcd, pdivmod, ptrim, count_roots, peval, squarefree, root_bound

HERE = os.path.dirname(__file__)


def exact_poles(a, b, lo, hi):
    """(pole in [lo, hi], real pole anywhere): exact, with the floats converted exactly.
       A root of Q is a genuine pole iff it is a root of Q / gcd(P, Q)."""
    P = ptrim([Fr(float(c)) for c in a]); Q = ptrim([Fr(float(c)) for c in b])
    if len(Q) <= 1: return (False, False)                    # constant denominator: no pole
    g = pgcd(P, Q) if P else Q
    Qr = squarefree(pdivmod(Q, g)[0]) if len(g) > 1 else squarefree(Q)
    if len(Qr) <= 1: return (False, False)
    lo, hi = Fr(float(lo)), Fr(float(hi))
    in_range = count_roots(Qr, lo, hi) > 0 or peval(Qr, lo) == 0
    B = root_bound(Qr)
    anywhere = count_roots(Qr, -B - 1, B + 1) > 0
    return (in_range, anywhere)


def float_check(b, lo, hi):
    z = np.linspace(lo, hi, 1000)
    Q = b[0] + z * (b[1] + z * b[2])
    return bool(np.min(np.abs(Q)) < 1e-3)


def spikes(series, start=100, factor=10.0, window=50):
    s = np.asarray(series, dtype=float)
    for t in range(max(start, window), len(s)):
        prev = s[t - window:t]
        if not np.isfinite(s[t]) or s[t] > factor * np.median(prev):
            return t
    return None


def first(flags):
    for t, f in enumerate(flags):
        if f: return t
    return None


def main():
    rows = []
    for fn in sorted(glob.glob(os.path.join(HERE, "runs", "*.npz"))):
        d = np.load(fn)
        n = len(d["loss"])
        exact_in, exact_any, flt = [], [], []
        for t in range(n):
            ein = eany = fl = False
            for L in ("1", "2"):
                a, b, (lo, hi) = d["a" + L][t], d["b" + L][t], d["r" + L][t]
                i, an = exact_poles(a, b, lo, hi)
                ein |= i; eany |= an; fl |= float_check(b, lo, hi)
            exact_in.append(ein); exact_any.append(eany); flt.append(fl)
        rows.append({"run": os.path.basename(fn), "lr": float(d["lr"]), "seed": int(d["seed"]), "steps": n,
                     "instability": spikes(d["loss"]), "pole_in_range": first(exact_in),
                     "pole_anywhere": first(exact_any), "float_check": first(flt),
                     "grad_spike": spikes(d["gnorm"]),
                     "pole_in_range_steps": int(sum(exact_in)), "float_steps": int(sum(flt)),
                     "disagree_steps": int(sum(1 for u, v in zip(exact_in, flt) if u != v))})
        r = rows[-1]
        print(f"{r['run']:22s} unstable@{r['instability']}  pole_in_range@{r['pole_in_range']}  "
              f"real_pole@{r['pole_anywhere']}  float@{r['float_check']}  gradspike@{r['grad_spike']}  "
              f"exact/float disagree on {r['disagree_steps']} steps", flush=True)

    def score(key, before=25, after=50):
        unstable = [r for r in rows if r["instability"] is not None]
        flagged = [r for r in rows if r[key] is not None]
        h1 = [r for r in unstable if r[key] is not None and r["instability"] - before <= r[key] <= r["instability"]]
        h2 = [r for r in flagged if r["instability"] is not None and r[key] <= r["instability"] <= r[key] + after]
        return {"unstable_runs": len(unstable), "flagged_runs": len(flagged),
                "H1_hits": len(h1), "H2_hits": len(h2),
                "H1_rate": (len(h1) / len(unstable)) if unstable else None,
                "H2_rate": (len(h2) / len(flagged)) if flagged else None}

    summary = {k: score(k) for k in ("pole_in_range", "float_check", "grad_spike", "pole_anywhere")}
    s = summary["pole_in_range"]
    def verdict(s):
        if s["unstable_runs"] < 5 or s["H1_rate"] is None or s["H2_rate"] is None:
            return "NOT TESTABLE (fewer than 5 unstable runs, or no flagged runs)"
        if s["H1_rate"] >= 0.9 and s["H2_rate"] >= 0.9: return "SUPPORTED"
        if s["H1_rate"] < 0.6 or s["H2_rate"] < 0.6: return "REFUTED"
        return "PARTIAL"
    summary["verdict_H"] = verdict(s)
    json.dump({"runs": rows, "summary": summary}, open(os.path.join(HERE, "results.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
