"""Does learning help? (python benchmarks/learning_benchmark.py [n_train] [n_test] [cap_s])

Three-variable formulas with a distance, where the choice of starting direction for the squaring
map can change the proof time from 0.1 s to never. Each formula gets at most cap_s seconds.

  A  single start (FMSL 0.25 behaviour): the first coordinate direction only
  B  multi-start without learning: a 5 s probe of each start, then the first one without limit
  C  multi-start with learning (fmsl.memory): the same, but starts are ordered by what worked on
     earlier formulas of the same kind; trained on n_train formulas, measured on n_test new ones
  D  memory: C's test formulas asked again in a new session: proven answers come from the saved
     proofs, and formulas not proven within the limit are reported at once (remembered failures)
Reported: how many are proven within the cap, and the total and median time.
"""
import os, random, signal, statistics, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src")); sys.path.insert(0, HERE)
import fmsl.memory as memory
import fmsl.verify as V
from fmsl.distance import reduce_distance


class Cap(BaseException): pass


def capped(fn, cap):
    def h(*_): raise Cap()
    old = signal.signal(signal.SIGALRM, h); signal.alarm(cap)
    t = time.time()
    try:
        r = fn(); ok = bool(r and r["proven"])
    except Cap:
        ok = False
    except Exception:
        ok = False
    finally:
        signal.alarm(0); signal.signal(signal.SIGALRM, old)
    return ok, min(time.time() - t, cap)


def formula(gen):
    """a rational function of x, y, z and a distance, with uneven degrees so that directions differ"""
    Q = gen.choice(["x^2+y^2+z^2", "x^2+y^2+z^2", "x^2+2*y^2+z^2"])
    R = f"sqrt({Q})"
    def mono(): return f"x^{gen.randint(0, 2)}*y^{gen.randint(0, 2)}*z^{gen.randint(0, 3)}"
    num = "+".join(f"({gen.choice([-2, -1, 1, 2])})*{mono()}*{R}^{gen.randint(0, 1)}" for _ in range(gen.randint(1, 2)))
    den = "+".join(f"({gen.choice([-2, -1, 1, 2])})*{mono()}*{R}^{gen.randint(0, 1)}" for _ in range(gen.randint(1, 3)))
    return f"({num})/({den})"


def run(fs, how, cap):
    res = []
    for f in fs:
        if how == "A":
            g = reduce_distance(f, ["x", "y", "z"], None, 0)
            ok, t = capped(lambda: V._limit_at_origin(g, ["x", "y", "z"]), cap)
        elif how == "B":
            ok, t = capped(lambda: V.limit_at_origin(f, ["x", "y", "z"]), cap)
        else:                                           # FMSL's own time limit, remembered by fmsl.memory
            t0 = time.time()
            try:
                r = V.limit_at_origin(f, ["x", "y", "z"], time_limit=cap); ok = bool(r["proven"])
            except Exception:
                ok = False
            t = min(time.time() - t0, cap)
        res.append((ok, t))
    return res


def line(name, res):
    ts = [t for _, t in res]
    return f"| {name} | {sum(ok for ok, _ in res)} of {len(res)} | {sum(ts):.0f} s | {statistics.median(ts):.1f} s |"


def main(n_train=30, n_test=30, cap=60, only_b=0):
    gen = random.Random(41)
    fs = [formula(gen) for _ in range(n_train + n_test)]
    train, test = fs[:n_train], fs[n_train:]
    memory.disable()
    if only_b:
        print(line("B multi-start", run(test, "B", cap))); return
    A = run(test, "A", cap)
    B = run(test, "B", cap)
    path = os.path.join(tempfile.mkdtemp(), "memory.json")
    memory.enable(path)
    memory._STATE["data"]["results"] = {}
    run(train, "C", cap)                                # learning only: these are not scored
    memory._STATE["data"]["results"] = {}               # no answers carried over, only strategy
    C = run(test, "C", cap)
    memory.disable(); memory.enable(path)               # a new session, same memory file
    D = run(test, "C", cap)
    print("| method | proven within %d s | total time | median |" % cap)
    print("|---|---|---|---|")
    for name, res in [("A single start (0.25)", A), ("B multi-start", B), ("C multi-start + learning", C),
                      ("D asked again, new session (memory)", D)]:
        print(line(name, res))
    print("\nlearned:", {b: v for b, v in memory.stats()["strategy"].items()})


if __name__ == "__main__":
    a = [int(x) for x in sys.argv[1:5]]
    main(*a)
