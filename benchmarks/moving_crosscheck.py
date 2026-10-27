"""Independent checks of fmsl.moving (moving targets).

For random families f(x, y, t) = N/D with coefficients polynomial in t:
  1. every snapshot reported by fmsl.moving (a sample time in each interval) is compared with
     the exact image of f_t on small squares |x|, |y| = r  (benchmarks/loop_crosscheck.py);
  2. the moving state at every rational Celedon point is compared with the exact image of
     f(x, y, t0 + z) on the surface of small cubes  (benchmarks/cube_crosscheck.py);
  3. completeness of the Celedon points: 12 extra sample times per interval; any change of
     type there that the module did not report is counted as MISSED.
Run:  python benchmarks/moving_crosscheck.py [n_random]
"""
import random, sys, time
from fractions import Fraction as Fr
sys.path.insert(0, "src"); sys.path.insert(0, "benchmarks")
from fmsl.moving import (parse_moving, moving_profile, snapshot, _at_time, _shift_time,
                         _samples_between)
from fmsl.multivar import profile_polys
from fmsl.algebraic import A
from fmsl.cluster import cmp
from loop_crosscheck import loop_image, coverage_line, mismatch_line, signed_ivs
import cube_crosscheck as cc

CLASSIC = ['(x^2-t*y^2)/(x^2+y^2)', 'x^2*y/(x^4+t*y^2)', 'x*y/(x^2+t*y^2)',
           '(x^3+t*y^3)/(x^2+y^2)', '(t*x^2+y^2-x*y)/(x^2+y^2)', '(x^2+t*y)/(x^2+y^2)',
           'x^2*y/(x^4+y^2+t*x^2*y)', '(x*y+t*x^2)/(x^2+(t-1)*y^2)']


def random_family(gen):
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


def snap_check(N, D, t0, r=Fr(1, 10 ** 16)):
    N2, D2 = _at_time(N, t0), _at_time(D, t0)
    if not N2 or not D2: return 0.0
    cl = profile_polys(N2, D2)["cluster"]
    im = loop_image(N2, D2, r)
    if im is None: return None
    return mismatch_line(coverage_line(im), coverage_line(signed_ivs(cl)))


def moving_check(N, D, t0, text):
    import math
    N3, D3 = _shift_time(N, t0), _shift_time(D, t0)
    im, ex = cc.cube_image(N3, D3, Fr(1, 10 ** 6))
    if im is None: return None
    f = lambda v: math.inf if v == "+oo" else (-math.inf if v == "-oo" else None)
    from fmsl.space import profile3_polys
    r = profile3_polys(N3, D3)
    return mismatch_line(coverage_line(cc.to_float(im)), coverage_line(cc.to_float(r["signed_components"])))


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    gen = random.Random(2026)
    funcs = CLASSIC + [random_family(gen) for _ in range(n)]
    T = time.time()
    snaps = snap_bad = movs = mov_bad = missed = fams = 0
    for f in funcs:
        try:
            N, D = parse_moving(f)
            r = moving_profile(f)
        except ZeroDivisionError:
            continue
        fams += 1
        bad = []
        for iv in r["intervals"]:
            m = snap_check(N, D, iv["sample"])
            if m is None: continue
            snaps += 1
            if m >= 0.002: snap_bad += 1; bad.append(f"snapshot t={iv['sample']} mismatch {m:.3f}")
        for s in r["special"]:
            if isinstance(s["t"], A): continue
            m = moving_check(N, D, s["t"], s.get("moving"))
            if m is None: continue
            movs += 1
            if m >= 0.002: mov_bad += 1; bad.append(f"moving t={s['t']} mismatch {m:.3f}")
        # completeness: dense samples inside each reported interval must all have one type
        for iv in r["intervals"]:
            lo, hi = iv["lo"], iv["hi"]
            pts = sorted(set(_samples_between(lo, hi)))
            extra = []
            for p, q in zip(pts, pts[1:]):
                extra += [p + (q - p) * k / 5 for k in range(1, 5)]
            sigs = set()
            for p in pts + extra:
                if any(cmp(p, s["t"]) == 0 for s in r["special"]): continue
                s = snapshot(N, D, p)
                sigs.add(s["sig"] if s else None)
            if len(sigs) > 1 and not r["unexplained"]:
                missed += 1; bad.append(f"type changes inside ({lo}, {hi})")
        sp = ", ".join(str(s["t"]) if not isinstance(s["t"], A) else f"{float(s['t']):.4g}" for s in r["special"])
        print(f"{f[:62]:62s} special: [{sp}]" + (f"  unexplained: {len(r['unexplained'])}" if r["unexplained"] else "")
              + ("".join("\n      <-- " + b for b in bad)), flush=True)
    print(f"\n{fams} families, {time.time() - T:.0f}s")
    print(f"snapshots checked against the square loop: {snaps - snap_bad}/{snaps} agree")
    print(f"moving states checked against the cube:    {movs - mov_bad}/{movs} agree")
    print(f"intervals with an unreported type change:  {missed}")
