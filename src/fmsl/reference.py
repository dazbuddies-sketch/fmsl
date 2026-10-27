"""
fmsl.reference -- high-precision reference values for matrix exponentials (pure Python).

Matrices with exact rational / Gaussian-rational entries are exponentiated in binary
fixed point with P = 200 fractional bits, using scaling and squaring of a Taylor series.
The result is accurate to roughly 1e-50 for the moderate sizes and norms used in tests,
far below double precision, so it serves as an independent "exact" answer.

  expm_ref(A)                 exp(A)
  phi_ref(k, A)               phi_k(A)   via exp of the augmented block matrix
  frechet_ref(A, E)           d/de exp(A + eE) via exp([[A, E], [0, A]])
Entries may be int, Fraction, or fmsl.gauss.G. Results are returned as complex numpy arrays.
"""
from fractions import Fraction as Fr

import numpy as np

from .gauss import G

P = 200
ONE = 1 << P


def _fx(v):
    v = G(v)
    return ((v.re.numerator << P) // v.re.denominator, (v.im.numerator << P) // v.im.denominator)


def _mul(A, B):
    n, m, q = len(A), len(B), len(B[0])
    out = []
    for i in range(n):
        Ai = A[i]
        row = []
        for k in range(q):
            re = im = 0
            for j in range(m):
                ar, ai = Ai[j]
                br, bi = B[j][k]
                if ar or ai:
                    re += ar * br - ai * bi
                    im += ar * bi + ai * br
            row.append((re >> P, im >> P))
        out.append(row)
    return out


def _norm1(A):
    return max(sum(abs(A[i][j][0]) + abs(A[i][j][1]) for i in range(len(A))) for j in range(len(A[0]))) / ONE


def _expm_fx(A, terms=40):
    n = len(A)
    s = 0
    nrm = _norm1(A)
    while nrm > 0.5:
        nrm /= 2; s += 1
    As = [[(a[0] >> s, a[1] >> s) for a in row] for row in A]   # exact shift: A / 2^s
    S = [[(ONE, 0) if i == j else (0, 0) for j in range(n)] for i in range(n)]
    T = [row[:] for row in S]
    for j in range(1, terms + 1):
        T = _mul(T, As)
        T = [[(a[0] // j, a[1] // j) for a in row] for row in T]
        S = [[(S[i][k][0] + T[i][k][0], S[i][k][1] + T[i][k][1]) for k in range(n)] for i in range(n)]
    for _ in range(s):
        S = _mul(S, S)
    return S


def _to_np(S):
    return np.array([[complex(a[0] / ONE, a[1] / ONE) for a in row] for row in S])


def _block(blocks):
    """assemble a block matrix from a 2D list of fixed-point blocks (None = zero)"""
    sizes = [next(len(b) for b in row if b is not None) for row in blocks]
    cols = [next(len(blocks[i][j][0]) for i in range(len(blocks)) if blocks[i][j] is not None)
            for j in range(len(blocks[0]))]
    out = []
    for bi, row in enumerate(blocks):
        for r in range(sizes[bi]):
            line = []
            for bj, blk in enumerate(row):
                line += blk[r] if blk is not None else [(0, 0)] * cols[bj]
            out.append(line)
    return out


def _fxmat(A):
    return [[_fx(v) for v in row] for row in A]


def expm_ref(A):
    return _to_np(_expm_fx(_fxmat(A)))


def phi_ref(k, A):
    """phi_k(A) = top-right block of exp of [[A, I, 0..], [0, 0, I, ..], ..., [0 ... 0]]"""
    n = len(A)
    Afx = _fxmat(A)
    I = [[(ONE, 0) if i == j else (0, 0) for j in range(n)] for i in range(n)]
    blocks = [[None] * (k + 1) for _ in range(k + 1)]
    blocks[0][0] = Afx
    for b in range(k):
        blocks[b][b + 1] = I
    Z = [[(0, 0)] * n for _ in range(n)]
    for i in range(k + 1):
        for j in range(k + 1):
            if blocks[i][j] is None: blocks[i][j] = Z
    S = _expm_fx(_block(blocks))
    return _to_np([row[k * n:(k + 1) * n] for row in S[:n]])


def frechet_ref(A, E):
    n = len(A)
    Afx, Efx = _fxmat(A), _fxmat(E)
    Z = [[(0, 0)] * n for _ in range(n)]
    S = _expm_fx(_block([[Afx, Efx], [Z, Afx]]))
    return _to_np([row[n:] for row in S[:n]])
