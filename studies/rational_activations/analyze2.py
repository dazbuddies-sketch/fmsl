"""Pre-registered analysis of experiment 2 (PREREGISTRATION_2.md).  Same exact pole test as
analyze.py; instability is absolute: loss not finite, or above its step-0 value (t >= 100)."""
import glob, json, os
import numpy as np
from analyze import exact_poles, float_check, spikes, first

HERE = os.path.dirname(__file__)


def instability(loss, start=100):
    L = np.asarray(loss, dtype=float)
    for t in range(start, len(L)):
        if not np.isfinite(L[t]) or L[t] > L[0]: return t
    return None


def main():
    rows = []
    for fn in sorted(glob.glob(os.path.join(HERE, "runs2", "*.npz"))):
        d = np.load(fn); n = len(d["loss"])
        ein, ean, flt = [], [], []
        for t in range(n):
            e1 = e2 = f = False
            for L in ("1", "2"):
                a, b, (lo, hi) = d["a" + L][t], d["b" + L][t], d["r" + L][t]
                if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b)) and np.isfinite(lo) and np.isfinite(hi)):
                    continue
                i, an = exact_poles(a, b, lo, hi); e1 |= i; e2 |= an
                z = np.linspace(lo, hi, 1000); Q = np.polyval(b[::-1], z); f |= bool(np.min(np.abs(Q)) < 1e-3)
            ein.append(e1); ean.append(e2); flt.append(f)
        r = {"run": os.path.basename(fn), "lr": float(d["lr"]), "seed": int(d["seed"]), "steps": n,
             "instability": instability(d["loss"]), "pole_in_range": first(ein), "pole_anywhere": first(ean),
             "float_check": first(flt), "grad_spike": spikes(d["gnorm"]),
             "disagree_steps": int(sum(1 for u, v in zip(ein, flt) if u != v)),
             "final_loss": float(d["loss"][-1]), "max_loss_after_100": float(np.nanmax(d["loss"][100:])) if n > 100 else None}
        rows.append(r)
        print(f"{r['run']:22s} unstable@{r['instability']}  pole_in_range@{r['pole_in_range']}  real_pole@{r['pole_anywhere']}  "
              f"float@{r['float_check']}  gradspike@{r['grad_spike']}  disagree={r['disagree_steps']}", flush=True)

    def score(key, before=25, after=50):
        unstable = [r for r in rows if r["instability"] is not None]
        flagged = [r for r in rows if r[key] is not None]
        h1 = [r for r in unstable if r[key] is not None and r["instability"] - before <= r[key] <= r["instability"]]
        h2 = [r for r in flagged if r["instability"] is not None and r[key] <= r["instability"] <= r[key] + after]
        return {"unstable_runs": len(unstable), "flagged_runs": len(flagged), "H1_hits": len(h1), "H2_hits": len(h2),
                "H1_rate": len(h1) / len(unstable) if unstable else None,
                "H2_rate": len(h2) / len(flagged) if flagged else None}
    summ = {k: score(k) for k in ("pole_in_range", "float_check", "grad_spike", "pole_anywhere")}
    s = summ["pole_in_range"]
    if s["unstable_runs"] < 5 or s["flagged_runs"] < 5: v = "NOT TESTABLE"
    elif s["H1_rate"] >= 0.9 and s["H2_rate"] >= 0.9: v = "SUPPORTED"
    elif s["H1_rate"] < 0.6 or s["H2_rate"] < 0.6: v = "REFUTED"
    else: v = "PARTIAL"
    summ["verdict_H"] = v
    json.dump({"runs": rows, "summary": summ}, open(os.path.join(HERE, "results2.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
