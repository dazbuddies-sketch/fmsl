"""Check of fmsl.indeterminacy (python benchmarks/indeterminacy_crosscheck.py [n]).

Random rational functions of x, y, z vanishing at the origin. Whenever the Part 2 test finds a clean
indeterminacy line, its answer (every value) is compared with FMSL's three-variable zooming engine
(fmsl.space.profile3), which knows nothing of the theorem. Reports agreements, disagreements, and
how many of the Part 2 answers the engine could not finish within the time limit.
"""
import os, random, signal, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from fmsl.space import parse3, profile3
from fmsl.indeterminacy import clean_line

LIMIT = 45


class Slow(Exception): pass


def rand(gen):
    def poly():
        terms = []
        for _ in range(gen.randint(1, 3)):
            e = [gen.randint(0, 2) for _ in range(3)]
            if sum(e) == 0: e[gen.randrange(3)] = 1
            terms.append(f"({gen.choice([-2, -1, 1, 2])})*x^{e[0]}*y^{e[1]}*z^{e[2]}")
        return "+".join(terms)
    return f"({poly()})/({poly()})"


def main(n=300):
    gen = random.Random(99)
    def h(*_): raise Slow()
    signal.signal(signal.SIGALRM, h)
    tally = {"fires, engine agrees": 0, "fires, engine DISAGREES": 0, "fires, engine not finished": 0,
             "does not fire": 0}
    t0 = time.time()
    for _ in range(n):
        f = rand(gen)
        try:
            N, D = parse3(f)
        except Exception:
            continue
        if clean_line(N, D) is None:
            tally["does not fire"] += 1; continue
        signal.alarm(LIMIT)
        try:
            r = profile3(f)
            st = r["signed_text"] if r.get("complete") and r.get("exact", True) else None
        except (Slow, Exception):
            st = None
        finally:
            signal.alarm(0)
        if st is None: tally["fires, engine not finished"] += 1
        elif st.strip() == "[-oo, +oo]": tally["fires, engine agrees"] += 1
        else:
            tally["fires, engine DISAGREES"] += 1; print("DISAGREE", f, st, flush=True)
    print(tally, f"({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 300)
