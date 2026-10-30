"""Explore Part 4: f = P/Q in 3 variables with the z-axis inside {P = Q = 0} but NOT clean (P in (x,y)^2).
Compare FS at the origin (3D engine) with FS at points p = (0,0,t) on the axis (J), t = +-1/1000, +-1/100."""
import sys, random, signal, os
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from fmsl.verify import limit_at, limit_at_origin
class Slow(Exception): pass
def h(*_): raise Slow()
signal.signal(signal.SIGALRM, h)
gen = random.Random(5)
def c(): return gen.choice([-2, -1, 1, 2])
def quad():   # an element of (x, y)^2, with z-dependence
    return "+".join(f"({c()})*{m}*{gen.choice(['1','z','(1+z)','(2-z)'])}" for m in gen.sample(["x^2", "x*y", "y^2"], gen.randint(1, 3)))
def lin():    # an element of (x, y)
    return "+".join(f"({c()})*{m}*{gen.choice(['1','z','(1+z)'])}" for m in gen.sample(["x", "y"], gen.randint(1, 2)))
def run(fn, *a):
    signal.alarm(40)
    try:
        r = fn(*a); return r["set"] if r["proven"] else None
    except Slow: return None
    except Exception as e: return "ERR"
    finally: signal.alarm(0)
rows = 0
for i in range(int(sys.argv[1]) if len(sys.argv) > 1 else 25):
    P = quad(); Q = gen.choice([quad(), quad() + "+" + f"({c()})*z^2*x"])
    f = f"({P})/({Q})"
    fs0 = run(limit_at_origin, f, ["x", "y", "z"])
    J = [run(limit_at, f, [0, 0, t], ["x", "y", "z"]) for t in ("1/1000", "-1/1000", "1/100")]
    print(i, f, "| FS(0) =", fs0, "| FS on axis at t=1e-3,-1e-3,1e-2:", J, flush=True)
