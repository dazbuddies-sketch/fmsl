"""The 21 examples of Dinh, Guo, Nguyen, Pham, 'Limits of real bivariate rational functions'
(arXiv:2202.04889, 2022, Section 5.3), run through FMSL.  Run:  python benchmarks/dgnp2022_examples.py"""
import sys
sys.path.insert(0, "src")
from fmsl.multivar import profile

F = {1: "x^4 + x^2*y + y^2", 2: "x^4 + 3*x^2*y - x^2 - y^2",
     3: "2*y^5 + x^2*y^2 - 8*x*y^3 - 13*y^4 - 2*x^3 + 6*x^2*y + 28*x*y^2 + 24*y^3 - 4*x^2 - 12*x*y - 9*y^2",
     4: "4*x^2*y^2 - 4*x*y^3 + y^4 - 2*x*y^2 + y^3",
     5: "10*x^2*y^2 + x^3 + 2*x^2*y + 4*x*y^2 + 6*x^2 + 6*x*y + 3*y^2",
     6: "2*x^2*y^2 + x^2*y + 2*x*y^2 + y^3 + x^2 + 2*x*y + 2*y^2",
     8: "10*x^2*y - 26*x^3 + 37*x^2*y^2 - 8*x*y^3 + 2*y^5 - 18*x*y^4 + 3*y^6"}
G = {1: "x^2 + y^2", 2: "x^2 + y^2",
     3: "y^4 - 5*x*y^2 - 4*y^3 + 7*x^2 + 10*x*y + 4*y^2",
     4: "8*x^2*y^2 - 8*x*y^3 + 3*y^4 + 8*x^2 - 8*x*y + 2*y^2",
     5: "3*x^2*y^2 + 2*x^3 + 2*x*y + y^2",
     6: "x^2*y^2 + x^2 + 2*x*y + 2*y^2",
     7: "3*x^2*y^2 + 2*x^2 + 2*x*y + y^2",
     8: "24*x^2 + 3*y^2 - 21*x*y^2 - 5*x*y + 2*y^3 + 5*y^4",
     9: "x^4*y^4 + x^2 + 2*x*y + 2*y^2"}
F[7] = F[5]; F[9] = F[6]
P = lambda *ps: "*".join(f"({p})" for p in ps)
EX = {i: (F[i], G[i]) for i in range(1, 10)}
EX[10] = (P(F[2], F[7]), P(G[2], G[7]))
EX[11] = (P(F[6], F[7]), P(G[6], G[7]))
EX[12] = (P(F[2], F[6]), P(G[2], G[6]))
EX[13] = (P(F[2], F[6]) + " + x^6*y^6", P(G[2], G[6]) + " + x^4*y^4")
EX[14] = (P(F[2], F[6]) + " + x^10*y^10", P(G[2], G[6]) + " + x^8*y^8")
EX[15] = (P(F[2], F[6], F[7]), P(G[2], G[6], G[7]))
EX[16] = (P(F[2], F[6], F[7], F[8]), P(G[2], G[6], G[7], G[8]))
EX[17] = (P(F[2], F[6], F[7], F[9]), P(G[2], G[6], G[7], G[9]))
EX[18] = ("x^2", "x^4 + y^4"); EX[19] = ("x^3", "x^4 + y^4")
EX[20] = ("x^4 + x^2*y + y^2", "x^6 + y^2"); EX[21] = ("x^4 + x^2*y^2 + y^4", "x^6 + y^4")
PAPER = {1: "[0, 1]", 2: "limit -1", 3: "[-19/3, 0]", 4: "[-sqrt2/4, sqrt2/4]",
         5: "(non-isolated zero: not computed)", 6: "limit 1", 7: "limit 3", 8: "limit 0", 9: "limit 1",
         10: "limit -3", 11: "limit 3", 12: "limit -1", 13: "limit -1", 14: "limit -1", 15: "limit -3",
         16: "limit 0", 17: "limit -3", 18: "[0, +oo]", 19: "[-oo, +oo]", 20: "[3/4, +oo]  (Maple 2021: (3/4, 1), wrong)",
         21: "[1, +oo]"}

if __name__ == "__main__":
    import time
    for i in sorted(EX):
        f, g = EX[i]
        t = time.time()
        r = profile(f"({f})/({g})")
        cl = r["cluster"]
        ours = cl.get("signed_text", cl["text"])
        print(f"Ex {i:2d}  paper: {PAPER[i]:45s}  FMSL: {ours}   [{r['verdict'][0]}]  {time.time() - t:.1f}s", flush=True)
