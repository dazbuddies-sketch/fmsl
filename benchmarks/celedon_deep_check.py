"""
Cross-check of the deep Celedon point search (fmsl.celedon, moving_profile(deep=True)).

For each family (classic, random, and families built to hide a change inside the zoom):
  * run moving_profile in the default mode and in deep mode;
  * compare: Celedon points found exactly, changes left as "change near ..." (unexplained);
  * dense check of deep mode: 12 rational times in every reported interval; any change of
    snapshot type not at a reported Celedon point is counted as a miss.

Run:  python benchmarks/celedon_deep_check.py [n_random] [n_hidden] [seed]
"""
import random, sys, time
from fractions import Fraction as Fr

sys.path.insert(0, "src")
from fmsl.moving import moving_profile, parse_moving, snapshot
from fmsl.cluster import cmp, _bounds

CLASSIC = ['(x^2-t*y^2)/(x^2+y^2)', 'x^2*y/(x^4+t*y^2)', 'x*y/(x^2+t*y^2)',
           '(x^3+t*y^3)/(x^2+y^2)', '(t*x^2+y^2-x*y)/(x^2+y^2)', '(x^2+t*y)/(x^2+y^2)',
           'x^2*y/(x^4+y^2+t*x^2*y)', '(x*y+t*x^2)/(x^2+(t-1)*y^2)']


def random_family(gen):                     # same generator as benchmarks/moving_crosscheck.py
    def coef():
        c = gen.choice([-2, -1, 1, 2])
        return gen.choice([f"({c})", f"({c})*t", f"({c})*(t-{gen.choice([1, 2, -1])})"])
    def poly(n, deg):
        terms = []
        for _ in range(n):
            i, j = gen.randint(0, deg), gen.randint(0, deg)
            if i + j == 0: continue
            terms.append(f"{coef()}*x^{i}*y^{j}")
        return " + ".join(terms) or "x"
    return f"({poly(gen.randint(1, 3), 3)})/({poly(gen.randint(2, 3), 4)})"


def hidden_family(gen):
    """the denominator vanishes to second order along y = c0(t) x, and the next zoom level
       carries a t-dependent coefficient that is invisible in the top-level coefficients"""
    c0 = gen.choice(["1", "2", "-1", "t", "(t+1)", "(1-t)", "2*t"])
    r = gen.choice([0, 1, -1, 2, Fr(1, 2)])
    b = gen.choice([1, 2, 3])
    # along y = c0 x + u:  x*y^3 = c0^3 x^4 + (higher order), so the x^4 coefficient at the
    # second zoom level is (t - r), while the top-level coefficients are t - r + b c0^3 and -b
    den = f"(y-({c0})*x)^2 + (t-({r})+{b}*({c0})^3)*x^4 - {b}*x*y^3"
    num = gen.choice(["x^4", "x^3*y", "x^4+x^2*y^2", "x^5"])
    return f"({num})/({den})"


def _fl(t):
    return float(t)


def dense_misses(text, res, per=12):
    """sample `per` times in every reported interval; count type changes inside an interval"""
    N, D = parse_moving(text)
    misses = []
    for iv in res["intervals"]:
        lo, hi = iv["lo"], iv["hi"]
        flo = _fl(_bounds(lo)[1]) if lo is not None else _fl(_bounds(hi)[0]) - 6 if hi is not None else -6.0
        fhi = _fl(_bounds(hi)[0]) if hi is not None else flo + 6
        if lo is None and hi is None: flo, fhi = -6.0, 6.0
        ts = []
        for k in range(1, per + 1):
            t = Fr(flo + (fhi - flo) * k / (per + 1)).limit_denominator(10 ** 6)
            if (lo is None or cmp(t, lo) > 0) and (hi is None or cmp(t, hi) < 0): ts.append(t)
        sigs = []
        for t in ts:
            s = snapshot(N, D, t)
            sigs.append(s["sig"] if s else None)
        for u in range(len(sigs) - 1):
            if sigs[u] != sigs[u + 1]: misses.append((ts[u], ts[u + 1]))
    return misses


def run(text):
    t0 = time.time()
    a = moving_profile(text, joint=False)
    t1 = time.time()
    b = moving_profile(text, joint=False, deep=True)
    t2 = time.time()
    m = dense_misses(text, b)
    return {"text": text,
            "default_cp": [str(s["t"]) for s in a["celedon_points"]],
            "default_unexplained": len(a["unexplained"]),
            "deep_cp": [str(s["t"]) for s in b["celedon_points"]],
            "deep_unexplained": len(b["unexplained"]),
            "deep_unsupported": len(b.get("deep_unsupported") or []),
            "dense_misses": len(m), "t_default": t1 - t0, "t_deep": t2 - t1}


if __name__ == "__main__":
    n_random = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    n_hidden = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 7
    gen = random.Random(seed)
    fams = [("classic", f) for f in CLASSIC]
    fams += [("random", random_family(gen)) for _ in range(n_random)]
    fams += [("hidden", hidden_family(gen)) for _ in range(n_hidden)]
    rows = []
    for kind, f in fams:
        try:
            r = run(f)
        except Exception as exc:                       # report, never hide
            r = {"text": f, "error": f"{type(exc).__name__}: {exc}"}
        r["kind"] = kind
        rows.append(r)
        if "error" in r:
            print(f"[{kind}] ERROR {f}: {r['error']}")
        else:
            print(f"[{kind}] {f}\n   default CP {r['default_cp']} unexplained {r['default_unexplained']}"
                  f" | deep CP {r['deep_cp']} unexplained {r['deep_unexplained']}"
                  f" unsupported {r['deep_unsupported']} dense misses {r['dense_misses']}"
                  f" | {r['t_default']:.1f}s / {r['t_deep']:.1f}s", flush=True)
    print("\nSUMMARY")
    for kind in ("classic", "random", "hidden"):
        R = [r for r in rows if r["kind"] == kind and "error" not in r]
        E = [r for r in rows if r["kind"] == kind and "error" in r]
        if not R and not E: continue
        print(f"  {kind}: {len(R)} run, {len(E)} errors; "
              f"default left {sum(r['default_unexplained'] > 0 for r in R)} with unexplained changes, "
              f"deep left {sum(r['deep_unexplained'] > 0 for r in R)}; "
              f"deep dense-sampling misses in {sum(r['dense_misses'] > 0 for r in R)}; "
              f"deep had unsupported branches in {sum(r['deep_unsupported'] > 0 for r in R)}; "
              f"lost points (default exact CP missing from deep) "
              f"{sum(any(c not in r['deep_cp'] for c in r['default_cp']) for r in R)}")
