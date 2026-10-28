"""Independent check of fmsl.optimize.global_min (python benchmarks/optimize_crosscheck.py [n]).

global_min finds the infimum m from the exact range (critical values, values at infinity) and
decides "attained" by whether m is a critical value. This check uses neither: it asks the exact
real-points test of fmsl.plane (one-level cylindrical decomposition, Sturm counts) directly:
  * below:    for a rational a < m, f = a has NO real solution     (no value below m)
  * above:    for a rational b > m, f = b HAS a real solution      (m is the infimum, not lower)
  * attained: f = m has a real solution  <=>  global_min says attained   (m rational only;
              isolated solutions included, by also searching the fibres over the breakpoints)
a and b are taken within 1e-6 of m (from an isolating interval when m is irrational).

Functions: classics with known answers, and random polynomials of the form
S1^2 + S2^2 + T, with S1, S2 random of degree <= 2 and T of degree <= 2 (bounded below; some attain, some not).
"""
import os, random, signal, sys, time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from fmsl.optimize import global_min
from fmsl.multivar import parse_rational
from fmsl.plane import (has_real_points, bclean, badd, bscale, in_d, from_d, _prim, _poly_roots_in,
                        resultant_d, bdd, at_c)
from fmsl.cluster import _bounds, _Inf

KNOWN = [  # (f, infimum, attained)
    ("x^2+(x*y-1)^2", Fr(0), False),
    ("(1-x*y)^2+x^2*y^4", Fr(0), False),
    ("x^4+y^4-2*x*y", Fr(-1, 2), True),
    ("x^2*y^4+x^4*y^2-3*x^2*y^2+1", Fr(0), True),        # Motzkin: >= 0, not a sum of squares
    ("1/(1+x^2+y^2)", Fr(0), False),
    ("(x^2-1)^2+(x^2*y-x-1)^2", Fr(0), True),
    ("x^2+y^2", Fr(0), True),
]


def rpoly(gen, deg, n):
    T = {}
    for _ in range(n):
        i, j = gen.randint(0, deg), gen.randint(0, deg)
        if i + j <= deg + 1: T[(i, j)] = T.get((i, j), 0) + gen.choice([-2, -1, 1, 2])
    return " + ".join(f"({v})*x^{i}*y^{j}" for (i, j), v in T.items() if v) or "1"


def randoms(n, seed=11):
    gen = random.Random(seed)
    for _ in range(n):
        yield f"({rpoly(gen, 1, 3)})^2 + ({rpoly(gen, 1, 3)})^2 + ({rpoly(gen, 1, 2)})"


def solvable(N, D, a):
    return has_real_points(bclean(badd(N, bscale(D, -a))), None)


def solvable_including_isolated(N, D, a):
    """f = a has a real solution, isolated points included: the curve test, plus every fibre over
       the breakpoints (leading coefficient, discriminant), where isolated real points can only lie"""
    F = bclean(badd(N, bscale(D, -a)))
    if has_real_points(F, None): return True
    L = in_d(F)
    if len(L) == 1: return bool(_poly_roots_in(L[0]))
    Lp, cont = _prim(L)
    if len(cont) > 1 and _poly_roots_in(cont): return True
    Fp = from_d(Lp)
    bps = list(_poly_roots_in(Lp[-1]))
    disc = resultant_d(Fp, bdd(Fp))
    if disc: bps += _poly_roots_in(disc)
    return any(_poly_roots_in(at_c(Fp, c0)) for c0 in bps)


def check(f):
    r = global_min(f)
    m = r.get("value_exact")
    N, D = parse_rational(f)
    if m is None:                                           # unbounded below: some value below -10^6 exists
        return r, solvable(N, D, Fr(-10 ** 6)), None
    lo, hi = (m, m) if isinstance(m, Fr) else _bounds(m)
    while hi - lo > Fr(1, 10 ** 6):
        m.F.refine(); lo, hi = _bounds(m)
    a, b = lo - Fr(1, 10 ** 6), hi + Fr(1, 10 ** 6)
    ok = (not solvable(N, D, a)) and solvable(N, D, b)
    att = solvable_including_isolated(N, D, m) if isinstance(m, Fr) else None
    return r, ok, att


def main(n=60):
    t0 = time.time(); bad = 0; total = 0; att_checked = 0
    for f, m0, a0 in KNOWN:
        r, ok, att = check(f); total += 1
        good = ok and r["attained"] == a0 and r["value"] == (str(m0.numerator) if m0.denominator == 1 else f"{m0.numerator}/{m0.denominator}")
        if att is not None: att_checked += 1; good = good and att == r["attained"]
        if not good: bad += 1; print("MISMATCH", f, r["text"], ok, att)
    n_att = n_not = slow = 0
    def alarm(*_): raise TimeoutError()
    signal.signal(signal.SIGALRM, alarm)
    for f in randoms(n):
        signal.alarm(60)
        try:
            r, ok, att = check(f)
        except ValueError:
            continue
        except TimeoutError:
            slow += 1; print("over 60 s (skipped):", f, flush=True); continue
        finally:
            signal.alarm(0)
        total += 1
        if att is not None:
            att_checked += 1
            if att != r["attained"]: ok = False
        n_att += r["attained"]; n_not += not r["attained"]
        if not ok: bad += 1; print("MISMATCH", f, r["text"], att)
    print(f"{total} functions ({len(KNOWN)} classics with known answers): {bad} disagree with the "
          f"independent real-points test; 'attained' checked directly on {att_checked}; random: "
          f"{n_att} attained, {n_not} not attained, {slow} skipped (over 60 s); {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 60)
