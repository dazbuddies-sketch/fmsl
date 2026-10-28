"""Independent check of fmsl.rate (python benchmarks/rate_crosscheck.py): exact rational sample points very close to the origin, on circles
and on curves y = c x^k (k rational), and exact log|f - L| / log r. The estimate must not go below
theta (that would refute |f-L| <= C r^theta) and should come close to it (attainment)."""
import json, math, os, random, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from fractions import Fraction as Fr
from fmsl.rate import rate
from fmsl.multivar import parse_rational
def ev(P,x,y): return sum(v*x**i*y**j for (i,j),v in P.items())
def lg(q): return math.log(q.numerator)-math.log(q.denominator)
def estimate(N,D,L,arc_ks):
    rnd=random.Random(5); best=1e9
    for _ in range(1500):
        mode=rnd.random(); sx=rnd.choice([-1,1])
        if mode<0.25:
            m=Fr(rnd.randint(-50,50),rnd.randint(1,50)); t=Fr(1,10**rnd.randint(10,16))
            x,y=t*(1-m*m)/(1+m*m), t*2*m/(1+m*m)
        else:
            k=rnd.choice(arc_ks); p,q=k.numerator,k.denominator
            mn=min(p,q); s=Fr(1,10**rnd.randint(-(-10//mn),-(-16//mn))); c=Fr(rnd.randint(-300,300),100)
            x,y=sx*s**q, c*s**p
            if rnd.random()<0.5: x,y=y,x
        r2=x*x+y*y; d=ev(D,x,y)
        if d==0 or r2==0: continue
        g=abs(ev(N,x,y)/d-L)
        if g>0: best=min(best, lg(g)/(0.5*lg(r2)))
    return best

KS = [Fr(1, 4), Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(1), Fr(3, 2), Fr(2), Fr(5, 2), Fr(3), Fr(4), Fr(5)]


def on_arc(N, D, L, data):
    """the slope of log|f - L| against log r between two tiny points of the reported arc
       (attainment): equals theta up to higher-order terms; the constant cancels"""
    from math import lcm
    if data is None: return None
    out = []
    for m in (24, 36):
        if "y_axis" in data:
            t = Fr(1, 10 ** m); x, y = Fr(0), data["y_axis"] * t
        else:
            ex = [e for _, e in data["terms"]] + ([data["generic_exponent"]] if data["generic_exponent"] is not None else [])
            Q = lcm(*[e.denominator for e in ex]) if ex else 1
            s = Fr(1, 10 ** max(2, m // Q)); ax = s ** Q
            y = Fr(0)
            for c, e in data["terms"]:
                cf = c if isinstance(c, Fr) else Fr(c.F.approx(c.c, Fr(1, 10 ** 80)))
                y += cf * s ** int(e * Q)
            if data["generic_exponent"] is not None: y += Fr(7, 3) * s ** int(data["generic_exponent"] * Q)
            x = data["sigma"] * ax
        d = ev(D, x, y)
        if d == 0: return None
        g = abs(ev(N, x, y) / d - L)
        if g == 0: return None
        out.append((lg(g), 0.5 * lg(x * x + y * y)))
    return (out[1][0] - out[0][0]) / (out[1][1] - out[0][1])


def cases():
    """LimitBench's two-variable functions with a limit, and a family that needs zooming, with a
    closed-form rate: x^m/((y - x^k)^2 + x^(2n)), m > 2n > 2k: theta = m - 2n, along y = x^k + c x^n."""
    for l in open(os.path.join(HERE, "limitbench", "limitbench.jsonl")):
        it = json.loads(l)
        if it["limit_exists"] and len(it["variables"]) == 2: yield it["f"], None
    for k in (1, 2, 3):
        for n in range(k + 1, k + 3):
            for m in (2 * n + 1, 2 * n + 2, 2 * n + 3):
                yield f"x^{m}/((y-x^{k})^2+x^{2 * n})", Fr(m - 2 * n)
                yield f"y^{m}/((x-y^{k})^2+y^{2 * n})", Fr(m - 2 * n)


def main():
    n = below = above = wrong = unproven = 0
    for f, known in cases():
        r = rate(f)
        if r["exponent"] is None:
            if not r["proven"]: unproven += 1; print("UNPROVEN", f)
            continue
        n += 1; th = r["exponent"]
        if known is not None and th != known: wrong += 1; print("WRONG closed form", f, th, known)
        N, D = parse_rational(f)
        est = estimate(N, D, r["limit"], KS)
        if est < float(th) - 0.05: below += 1; print("BELOW (refutes the bound)", f, th, round(est, 3))
        att = on_arc(N, D, r["limit"], r.get("arc_data"))
        if att is not None and att < float(th) - 0.02: below += 1; print("BELOW on its own arc", f, th, round(att, 3))
        elif att is None or att > float(th) + 0.02:
            above += 1; print("not attained on the reported arc", f, th, att, r["arc"])
    print(f"{n} rates checked: {wrong} differ from the closed form, {below} refuted by sampling, "
          f"{above} not attained on the reported arc, {unproven} unproven")


if __name__ == "__main__":
    main()
