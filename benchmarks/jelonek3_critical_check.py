"""Re-check the critical-value verdicts of the 72 candidates two ways.
(1) Levenberg-Marquardt minimisation of F1^2 + F2^2 + det^2 (works when the solutions form a curve,
    where plain Newton fails).  (2) For plainly triangular maps (x, y, h): exact criterion -- the
    critical values meet the line u1 = u2 = 0 iff h_z(0, 0, z) has a real root or is identically 0
    (Sturm, exact)."""
import ast, random, re, sys
from fractions import Fraction as Fr
sys.path.insert(0, '/home/claude/fmsl-repo/src')
from fmsl.jelonek_n import _polys, jacobian, _deriv
from fmsl.witness import _roots, _utrim
V = ['x', 'y', 'z']
lines = [l for l in open(sys.argv[1]) if l.startswith('CODIM-2')]
def fl(Q): return [(k, float(v)) for k, v in Q.items()]
def ev(T, x):
    s = 0.0
    for k, v in T:
        m = v
        for xi, e in zip(x, k):
            if e: m *= xi ** e
        s += m
    return s
def lm(G, D, gen, starts=200):
    best = 1e300
    for _ in range(starts):
        x = [gen.uniform(-3, 3) for _ in V]; lam = 1e-3
        for _ in range(200):
            g = [ev(T, x) for T in G]; f0 = sum(a*a for a in g)
            Jm = [[ev(T, x) for T in row] for row in D]
            JTJ = [[sum(Jm[k][i]*Jm[k][j] for k in range(3)) + (lam if i == j else 0) for j in range(3)] for i in range(3)]
            JTg = [sum(Jm[k][i]*g[k] for k in range(3)) for i in range(3)]
            M = [r[:] + [b] for r, b in zip(JTJ, JTg)]
            try:
                for c in range(3):
                    p = max(range(c, 3), key=lambda i: abs(M[i][c])); M[c], M[p] = M[p], M[c]
                    for i in range(3):
                        if i != c:
                            f = M[i][c]/M[c][c]; M[i] = [a-f*b for a, b in zip(M[i], M[c])]
                dx = [M[i][3]/M[i][i] for i in range(3)]
            except ZeroDivisionError:
                break
            y = [a-b for a, b in zip(x, dx)]
            if max(abs(a) for a in y) > 1e4: break
            f1 = sum(ev(T, y)**2 for T in G)
            if f1 < f0: x, lam = y, lam/3
            else: lam *= 4
            if f0 < 1e-26: break
        best = min(best, sum(ev(T, x)**2 for T in G))
        if best < 1e-24: return best, x
    return best, None
gen = random.Random(2)
agree = disagree = 0
for l in lines:
    F = ast.literal_eval(re.search(r"(\[.*?\]) \|", l).group(1))
    P = _polys(F, V); J = jacobian(P, 3)
    G = [fl(P[0]), fl(P[1]), fl(J)]; D = [[fl(_deriv(g, k)) for k in range(3)] for g in (P[0], P[1], J)]
    best, x = lm(G, D, gen)
    num = "meets" if best < 1e-20 else "clear"
    exact = None
    if F[0] == 'x' and F[1] == 'y':
        hz = _deriv(P[2], 2)
        uni = [Fr(0)] * (1 + max([k[2] for k in hz] + [0]))
        for k, v in hz.items():
            if k[0] == 0 and k[1] == 0: uni[k[2]] += v
        uni = _utrim(uni)
        exact = "meets" if not uni or (len(uni) > 1 and _roots(uni, Fr(-10**9), Fr(10**9))) else "clear"
    if exact is not None: agree += exact == num; disagree += exact != num
    print(num, "(exact: %s)" % exact if exact else "", F, "min %.2e" % best, flush=True)
print("plainly triangular: LM agrees with the exact criterion", agree, "times, disagrees", disagree)
