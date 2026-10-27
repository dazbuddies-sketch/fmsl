"""Matrix tests (skipped if numpy is not installed).  python -m unittest discover tests"""
import unittest
from fractions import Fraction as Fr
from math import factorial

try:
    import numpy as np
    from fmsl import matrix as fm
    from fmsl.reference import phi_ref, frechet_ref, expm_ref
    HAVE_NP = True
except ImportError:
    HAVE_NP = False

from fmsl.gauss import G


def exact_hermitian(D):
    """H = Q diag(D) Q^dagger with an exact rational unitary Q (4x4)."""
    r1 = [[G(Fr(3, 5)), G(Fr(-4, 5))], [G(Fr(4, 5)), G(Fr(3, 5))]]
    r2 = [[G(Fr(5, 13)), G(Fr(-12, 13))], [G(Fr(12, 13)), G(Fr(5, 13))]]
    Q = [[r1[i // 2][j // 2] * r2[i % 2][j % 2] for j in range(4)] for i in range(4)]
    ph = [G(Fr(3, 5), Fr(4, 5)), G(1), G(0, 1), G(1)]
    Q = [[Q[i][j] * ph[j] for j in range(4)] for i in range(4)]
    H = [[sum((Q[i][k] * G(D[k]) * Q[j][k].conj() for k in range(4)), G(0)) for j in range(4)] for i in range(4)]
    return H


def to_np(M):
    return np.array([[complex(v) for v in row] for row in M])


@unittest.skipUnless(HAVE_NP, "numpy not installed")
class ScalarPhi(unittest.TestCase):
    def ref(self, k, z):
        Z = G(Fr(z.real), Fr(z.imag)); acc, term = G(0), G(1)
        for j in range(int(60 + 4 * abs(z))):
            acc = acc + term / factorial(j + k); term = term * Z
        return complex(acc)

    def test_accuracy_all_regimes(self):
        rng = np.random.default_rng(3)
        mags = 10.0 ** rng.uniform(-12, 1.3, 300)
        zs = np.concatenate([mags * np.exp(1j * rng.uniform(0, 6.3, 300)), 1j * mags[:100], -mags[:100]])
        for k in (1, 2, 3, 4):
            f = fm.phi_scalar(k, zs)
            r = np.array([self.ref(k, z) for z in zs])
            self.assertLess((np.abs(f - r) / np.abs(r)).max(), 2e-15, f"phi_{k}")

    def test_exact_zero(self):
        self.assertEqual(fm.phi_scalar(2, np.array([0j]))[0], 0.5)

    def test_proven_bound_small(self):
        for k in (1, 2, 3, 4):
            b = fm.proven_phi_bound(k)
            self.assertIsNotNone(b)
            self.assertLess(float(b), 1e-12)


@unittest.skipUnless(HAVE_NP, "numpy not installed")
class MatrixFunctions(unittest.TestCase):
    def check(self, D, tol=5e-15):
        Hx = exact_hermitian(D); H = to_np(Hx)
        dHx = [[G(Fr(i + j, 7), Fr(i - j, 5)) for j in range(4)] for i in range(4)]   # Hermitian
        dH = to_np(dHx)
        A = [[v * G(0, -1) for v in row] for row in Hx]
        E = [[v * G(0, -1) for v in row] for row in dHx]
        def rel(X, R): return np.linalg.norm(X - R) / np.linalg.norm(R)
        self.assertLess(rel(fm.propagator(H, 1.0), expm_ref(A)), tol)
        self.assertLess(rel(fm.frechet(H, dH, 1.0), frechet_ref(A, E)), tol)
        for k in (1, 2, 3):
            self.assertLess(rel(fm.phi(k, H, 1.0), phi_ref(k, A)), tol)

    def test_degenerate_with_zero(self):
        self.check([0, 0, 1, 1])

    def test_near_degenerate(self):
        self.check([0, Fr(1, 10**10), 1, 1 + Fr(1, 10**12)])

    def test_weak_coupling(self):
        self.check([Fr(-3, 10**7), Fr(1, 10**8), Fr(2, 10**7), Fr(5, 10**7)])

    def test_frechet_many_matches_single(self):
        rng = np.random.default_rng(0)
        H = rng.normal(size=(6, 6)); H = H + H.T
        ds = [(lambda d: d + d.T)(rng.normal(size=(6, 6))) for _ in range(3)]
        many = fm.frechet_many(H, ds, 0.7)
        for d, L in zip(ds, many):
            self.assertLess(np.abs(L - fm.frechet(H, d, 0.7)).max(), 1e-14)


if __name__ == "__main__":
    unittest.main()
