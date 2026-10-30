"""Part 2 test: P = x a + y b, Q = x c + y d share the z-axis. If det[[a,b],[c,d]] is not identically 0
on the z-axis, the z-axis has clean points (dP, dQ independent) arbitrarily close to 0, and the
theorem predicts the fluid state is the whole circle. Compare with FMSL's 3D engine."""
import sys, random, signal, time
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from fractions import Fraction as Fr
from fmsl.space import profile3
class Slow(Exception): pass
def h(*_): raise Slow()
signal.signal(signal.SIGALRM, h)
gen = random.Random(7)
def lin():   # a small polynomial: constant + linear + a z term
    c = [gen.choice([-2, -1, 0, 1, 2]) for _ in range(4)]
    return c, f"(({c[0]})+({c[1]})*x+({c[2]})*y+({c[3]})*z)"
res = {"predicted whole, engine whole": 0, "predicted whole, engine NOT whole": 0, "clean but unproven/slow": 0,
       "not clean (no prediction)": 0}
examples = []
for i in range(int(sys.argv[1]) if len(sys.argv) > 1 else 40):
    (ca, a), (cb, b), (cc, c), (cd, d) = lin(), lin(), lin(), lin()
    # det on the z-axis: (a0 + a3 z)(d0 + d3 z) - (b0 + b3 z)(c0 + c3 z)
    det = [ca[0]*cd[0]-cb[0]*cc[0], ca[0]*cd[3]+ca[3]*cd[0]-cb[0]*cc[3]-cb[3]*cc[0], ca[3]*cd[3]-cb[3]*cc[3]]
    clean = any(det)
    P, Q = f"x*{a}+y*{b}", f"x*{c}+y*{d}"
    f = f"({P})/({Q})"
    signal.alarm(60)
    try:
        r = profile3(f); st = r["signed_text"] if r.get("complete") and r.get("exact", True) else None
    except Slow: st = None
    except Exception as e: st = None
    finally: signal.alarm(0)
    if not clean: res["not clean (no prediction)"] += 1; examples.append(("not clean", f, st)); continue
    if st is None: res["clean but unproven/slow"] += 1; continue
    whole = st.strip() == "[-oo, +oo]"
    res["predicted whole, engine whole" if whole else "predicted whole, engine NOT whole"] += 1
    if not whole: print("COUNTEREXAMPLE?", f, st, flush=True)
print(res)
for e in examples[:6]: print(*e)
