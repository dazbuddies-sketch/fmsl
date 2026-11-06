"""Search for polynomial maps F: R^3 -> R^3 whose Jelonek set S_F is a nonempty curve
(codimension 2), and classify their Jacobians (Jelonek's question, November 2026).

For each map: proven points of S_F with arcs (inner bound); the exact outer bound V(a_1) u ...
from polynomial relations, in two coordinate systems; a random-line test of the outer bound's
dimension; the same line test on det J_F. A map is a codimension-2 candidate when S_F has a
proven point and no line out of 200 meets the outer bound.
Jacobian classes (line test on det J, then signs at sample points):
  surface   det J vanishes on a surface (some random line crosses its zero set)
  small     det J vanishes somewhere, but no random line met its zero set (a curve or points)
  nonzero   det J never vanished on the lines or at 2,000 sample points (would be a candidate
            counterexample to Jelonek's real Jacobian conjecture; it would need a proof)
Usage: python benchmarks/jelonek3_search.py [n_maps] [seed]
"""
import os, random, signal, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from fmsl.jelonek_n import analyse, _restrict, _roots
from fmsl.witness import _ev

V = ["x", "y", "z"]


class Slow(Exception): pass


def mono(gen, dmax):
    e = [0, 0, 0]
    for _ in range(gen.randint(1, dmax)): e[gen.randrange(3)] += 1
    return "*".join(f"{v}^{k}" for v, k in zip(V, e) if k) or "1"


def poly(gen, terms, dmax):
    return "+".join(f"({gen.choice([-2, -1, 1, 1, 2])})*{mono(gen, dmax)}" for _ in range(terms))


def family(gen, k):
    kind = k % 5
    if kind == 4:     # a triangular codimension-2 map composed with a polynomial automorphism:
                      # S_F is unchanged, but the map is no longer triangular (a test of the tool)
        A = gen.choice(["(x^2+y^2)", "(x^2+2*y^2)"])
        h = f"{A}*z^2+{gen.choice(['z', '-z', 'x+z', 'y-2*z'])}"
        sub = {"x": "x", "y": f"(y+{gen.choice(['x^2', 'x', '2*x^2'])})", "z": f"(z+{gen.choice(['x*y', 'y^2', 'x^2'])})"}
        comp = lambda e: "".join(sub.get(ch, ch) for ch in e)
        return [comp("x"), comp("y"), comp(h)]
    if kind == 0:     # triangular: (x, y, h), leading z-coefficient vanishing only at (0, 0)
        A = gen.choice(["(x^2+y^2)", "(x^2+2*y^2)", "(x^2+x*y+y^2)", "(x^4+y^2)"])
        h = f"{A}*z^{gen.choice([2, 3])}+{poly(gen, gen.randint(1, 3), 2)}"
        return ["x", "y", h]
    if kind == 1:     # (x, g, h)
        return ["x", poly(gen, gen.randint(1, 3), 3), poly(gen, gen.randint(1, 3), 3)]
    if kind == 2:     # quadratic maps with a linear part
        return [f"{v}+{poly(gen, gen.randint(0, 2) or 1, 2)}" for v in V]
    return [poly(gen, gen.randint(1, 3), 3) for _ in V]


def det_class(J, gen, lines=200):
    for _ in range(lines):
        p0 = [Fr(gen.randint(-500, 500), 89) for _ in V]
        v = [Fr(gen.randint(-97, 97), gen.randint(1, 97)) for _ in V]
        u = _restrict(J, p0, v)
        if len(u) > 1 and _roots(u, Fr(-10 ** 6), Fr(10 ** 6)): return "surface"
        if not u: return "surface"
    signs = set()
    for _ in range(2000):
        x = [Fr(gen.randint(-300, 300), 37) for _ in V]
        v = _ev(J, x)
        signs.add((v > 0) - (v < 0))
    return "small" if 0 in signs or len(signs) > 1 else "nonzero"


def main(N=200, seed=7):
    gen = random.Random(seed)
    tally = {}
    def h(*_): raise Slow()
    for k in range(N):
        F = family(gen, k)
        signal.signal(signal.SIGALRM, h); signal.alarm(30)
        try:
            r = analyse(F, V)
        except (Slow, Exception) as e:
            tally["unfinished"] = tally.get("unfinished", 0) + 1; continue
        finally:
            signal.alarm(0)
        if r["outer"] is None: tally["no relation"] = tally.get("no relation", 0) + 1; continue
        hits, tried = r["surface_lines"]
        nonempty = bool(r["inner_points"])
        if not nonempty and hits == 0: cat = "S_F empty or small, no proven point"
        elif not nonempty: cat = "outer bound has a surface, no proven point"
        elif hits > 0: cat = "S_F has a proven point; outer bound has a surface"
        else: cat = "CODIM-2 CANDIDATE"
        if cat == "CODIM-2 CANDIDATE":
            dc = det_class(r["jacobian"], gen)
            cat += f" / det J {dc}"
            print(cat, F, "| points e.g.", [tuple(map(str, p)) for p in r["inner_points"][:3]], flush=True)
        tally[cat] = tally.get(cat, 0) + 1
    for c, v in sorted(tally.items()): print(f"{v:5d}  {c}")


if __name__ == "__main__":
    main(*(int(a) for a in sys.argv[1:3]))
