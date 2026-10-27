"""Matrix benchmark: FMSL vs textbook spectral formulas vs scipy, against an exact reference.
Run:  python benchmarks/matrix_bench.py"""
import time
from fractions import Fraction as Fr
from math import factorial

import numpy as np
import scipy.linalg as sl

from fmsl.gauss import G
from fmsl import matrix as fm
from fmsl.reference import phi_ref, frechet_ref

# ---------- exact Hermitian test matrices with prescribed eigenvalues ----------
def rot(c, s): return [[G(c), G(-s)], [G(s), G(c)]]
def kron(A, B):
    return [[A[i // len(B)][j // len(B)] * B[i % len(B)][j % len(B)]
             for j in range(len(A) * len(B))] for i in range(len(A) * len(B))]
def matmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), G(0)) for j in range(len(B[0]))] for i in range(len(A))]
def dagger(A): return [[A[j][i].conj() for j in range(len(A))] for i in range(len(A[0]))]

Q = kron(kron(rot(Fr(3, 5), Fr(4, 5)), rot(Fr(5, 13), Fr(12, 13))), rot(Fr(8, 17), Fr(15, 17)))
phases = [G(Fr(3, 5), Fr(4, 5)), G(1), G(Fr(5, 13), Fr(12, 13)), G(1), G(Fr(8, 17), Fr(-15, 17)), G(1), G(0, 1), G(1)]
Q = [[Q[i][j] * phases[j] for j in range(8)] for i in range(8)]          # exact unitary

def hermitian(D):
    Dm = [[G(D[i]) if i == j else G(0) for j in range(8)] for i in range(8)]
    return matmul(matmul(Q, Dm), dagger(Q))

def tonp(M): return np.array([[complex(v) for v in row] for row in M])

rng = np.random.default_rng(7)
def rational_hermitian(n=8):
    M = [[G(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            re = Fr(int(rng.integers(-6, 7)), 7)
            im = Fr(0) if i == j else Fr(int(rng.integers(-6, 7)), 7)
            M[i][j] = G(re, im); M[j][i] = G(re, -im)
    return M

# ---------- textbook spectral implementations (what is commonly written) ----------
def naive_phi(k, H, t, guard=True):
    E, V = np.linalg.eigh(H); a = -1j * E * t
    with np.errstate(all="ignore"):
        f = (np.exp(a) - sum(a ** j / factorial(j) for j in range(k))) / a ** k
    if guard: f = np.where(a == 0, 1 / factorial(k), f)
    return (V * f) @ V.conj().T

def naive_frechet(H, dH, t, tol=0.0):
    E, V = np.linalg.eigh(H); a = -1j * E * t
    D = a[:, None] - a[None, :]
    ea = np.exp(a)
    with np.errstate(all="ignore"):
        Gm = (ea[:, None] - ea[None, :]) / D
    Gm = np.where(np.abs(D) <= tol, (ea[:, None] + ea[None, :]) / 2, Gm)
    Vh = V.conj().T
    return V @ (Gm * (Vh @ (-1j * t * dH) @ V)) @ Vh

def scipy_phi(k, A):
    n = len(A); B = np.zeros(((k + 1) * n, (k + 1) * n), dtype=complex)
    B[:n, :n] = A
    for b in range(k): B[b * n:(b + 1) * n, (b + 1) * n:(b + 2) * n] = np.eye(n)
    return sl.expm(B)[:n, k * n:(k + 1) * n]

def err(X, R):
    if not np.all(np.isfinite(X)): return float("inf")
    return np.linalg.norm(X - R) / np.linalg.norm(R)

CASES = {
    "exact degeneracy + zero energy": [0, 0, 1, 1, 1, 2, 3, 3],
    "near-degenerate (gaps 1e-9, 1e-12)": [0, Fr(1, 10**9), 1, 1 + Fr(1, 10**9), 1 + Fr(2, 10**9), 2, 3, 3 + Fr(1, 10**12)],
    "weak coupling (all |E| <= 1e-6)": [Fr(k, 10**7) for k in (-9, -5, -2, 0, 1, 3, 6, 10)],
    "generic spectrum": [Fr(-13, 10), Fr(-7, 10), Fr(1, 5), Fr(1, 2), Fr(11, 10), Fr(19, 10), Fr(12, 5), 3],
}
t = Fr(1)
dHx = rational_hermitian()
dH = tonp(dHx)

def fmt(e): return "  FAILS (inf/NaN)" if e == float("inf") else f"{e:9.1e}"

print("Relative error (Frobenius) vs exact 200-bit reference; 8x8 Hamiltonians, t = 1\n")
for name, D in CASES.items():
    Hx = hermitian(D); H = tonp(Hx)
    A = [[v * G(0, -1) * t for v in row] for row in Hx]                   # A = -i H t (exact)
    Ex = [[v * G(0, -1) * t for v in row] for row in dHx]
    An = tonp(A); En = tonp(Ex)
    print(f"== {name}")
    R = frechet_ref(A, Ex)
    print(f"  gradient d exp(-iHt)/dH:  FMSL {fmt(err(fm.frechet(H, dH, float(t)), R))}"
          f" | textbook (exact-equality guard) {fmt(err(naive_frechet(H, dH, float(t)), R))}"
          f" | textbook (tol 1e-8) {fmt(err(naive_frechet(H, dH, float(t), 1e-8), R))}"
          f" | scipy expm_frechet {fmt(err(sl.expm_frechet(An, En, compute_expm=False), R))}")
    for k in (1, 2, 3):
        R = phi_ref(k, A)
        print(f"  phi_{k}(-iHt):              FMSL {fmt(err(fm.phi(k, H, float(t)), R))}"
              f" | textbook (zero guard)            {fmt(err(naive_phi(k, H, float(t)), R))}"
              f" |                     "
              f" | scipy augmented expm {fmt(err(scipy_phi(k, An), R))}")
    print()

print("Speed: quantum-control gradient, m control directions sharing one Hamiltonian")
for n, m in ((32, 10), (64, 10), (128, 10)):
    Hr = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); Hr = (Hr + Hr.conj().T) / 2
    dHs = []
    for _ in range(m):
        d = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); dHs.append((d + d.conj().T) / 2)
    tt = 0.3
    for _ in range(2):     # warm-up
        fm.frechet_many(Hr, dHs[:1], tt); sl.expm_frechet(-1j * Hr * tt, -1j * tt * dHs[0], compute_expm=False)
    t0 = time.perf_counter(); Lf = fm.frechet_many(Hr, dHs, tt); tf = time.perf_counter() - t0
    t0 = time.perf_counter()
    Ls = [sl.expm_frechet(-1j * Hr * tt, -1j * tt * d, compute_expm=False) for d in dHs]
    ts = time.perf_counter() - t0
    agree = max(np.linalg.norm(a - b) / np.linalg.norm(b) for a, b in zip(Lf, Ls))
    print(f"  n={n:4d}, m={m}: FMSL {1000*tf:8.1f} ms | scipy {1000*ts:8.1f} ms | speedup {ts/tf:5.1f}x | agreement {agree:.1e}")
