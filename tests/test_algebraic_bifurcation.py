"""Exact algebraic numbers and bifurcations.  python -m unittest discover tests"""
import unittest
from fractions import Fraction as Fr

from fmsl.algebraic import AlgField, A, isolate_real_roots
from fmsl.bifurcation import (iterate, divide_exact, shift, branches, period_polynomial,
                              fixed_polynomial, _add, branch_stability)
from fmsl.multivar import profile


def sqrt2():
    iv = isolate_real_roots([Fr(-2), 0, Fr(1)])[-1]
    return A.root(AlgField([Fr(-2), 0, Fr(1)], iv[1], iv[2]))


class Algebraic(unittest.TestCase):
    def test_arithmetic_is_exact(self):
        s = sqrt2()
        self.assertEqual(s * s, 2)
        self.assertEqual((1 / s) * s, 1)
        self.assertTrue(s > Fr(1414, 1000) and s < Fr(1415, 1000))
        self.assertNotEqual(s, Fr(3, 2))

    def test_reducible_definition_shrinks(self):
        p = [Fr(6), 0, Fr(-5), 0, Fr(1)]                      # (c^2-2)(c^2-3)
        iv = isolate_real_roots(p)[-1]
        F = AlgField(p, iv[1], iv[2]); s3 = A.root(F)
        self.assertEqual(s3 * s3, 3)
        self.assertNotEqual(s3 * s3, 2)
        self.assertEqual(len(F.p), 3)                            # now c^2 - 3

    def test_irrational_balance_constant_is_resolved(self):
        r = profile("x^3/(y^2 - 2*x^2)")
        self.assertEqual(r["verdict"][0], "NO LIMIT (proven)")
        self.assertTrue(any("algebraic" in e["path"] for e in r["entries"]))


class Bifurcations(unittest.TestCase):
    logistic = {(1, 1): Fr(1), (2, 1): Fr(-1)}

    def test_period_doubling_at_3(self):
        P = shift(period_polynomial(self.logistic, 2), Fr(2, 3), 3)
        self.assertEqual(P.get((0, 0), 0), 0)                     # 0/0 point of (f^2 - x)/(f - x)
        bs = [b for b in branches(P) if "a" in b]
        self.assertTrue(bs and all(b["a"] == 2 and b["c"] == 9 for b in bs))   # e = 9 u^2

    def test_classic_normal_forms(self):
        pitch = fixed_polynomial({(1, 0): Fr(1), (1, 1): Fr(1), (3, 0): Fr(-1)})
        self.assertTrue(any(b.get("a") == 2 and b.get("c") == 1 for b in branches(pitch)))
        trans = shift(fixed_polynomial(self.logistic), 0, 1)
        self.assertTrue(any(b.get("a") == 1 for b in branches(trans)))

    def test_second_period_doubling_exact_coefficient(self):
        p = [Fr(1), 0, Fr(-10), 0, Fr(1)]                        # alpha = sqrt2 + sqrt3
        iv = isolate_real_roots(p)[-1]
        al = A.root(AlgField(p, iv[1], iv[2]))
        s6 = (al ** 2 - 5) / 2; s2 = (al ** 3 - 9 * al) / 2; s3 = (11 * al - al ** 3) / 2
        r0 = 1 + s6
        x0 = (2 + s6 + s2) / (2 * (1 + s6))
        Q4 = divide_exact(_add(iterate(self.logistic, 4), {(1, 0): Fr(-1)}),
                          _add(iterate(self.logistic, 2), {(1, 0): Fr(-1)}))
        bs = [b for b in branches(shift(Q4, x0, r0)) if "a" in b]
        expected = 20 + 10 * s2 + Fr(35, 3) * s3 + Fr(35, 3) * s6
        self.assertTrue(bs and all(b["a"] == 2 and b["c"] == expected for b in bs))
        # the new period-4 orbit is attracting; its multiplier is 1 - 8 sqrt6 * e + ...
        lam0, k, m, verdict = branch_stability(self.logistic, 4, x0, r0, bs[0])
        self.assertEqual((lam0, verdict), (1, "stable"))
        self.assertTrue(k / bs[0]["c"] == -8 * s6)

    def test_stability_of_classic_branches(self):
        sn = {(1, 0): Fr(1), (0, 1): Fr(1), (2, 0): Fr(-1)}             # x + r - x^2
        res = {b["sigma"]: branch_stability(sn, 1, 0, 0, b)[3] for b in branches(fixed_polynomial(sn)) if "a" in b}
        self.assertEqual(res, {1: "stable", -1: "unstable"})
        b = [b for b in branches(shift(period_polynomial(self.logistic, 2), Fr(2, 3), 3)) if "a" in b][0]
        lam0, k, m, verdict = branch_stability(self.logistic, 2, Fr(2, 3), 3, b)
        self.assertEqual((lam0, k, m, verdict), (1, -36, 2, "stable"))     # 1 - 36 u^2 = 1 - 4e


if __name__ == "__main__":
    unittest.main()


class TestAlgebraicCoefficients(unittest.TestCase):
    def test_two_fields(self):
        from fmsl.algebraic import A, AlgField
        r2 = A.root(AlgField([-2, 0, 1], 1, 2)); r3 = A.root(AlgField([-3, 0, 1], 1, 2))
        s = r2 + r3
        self.assertEqual(s * s, 5 + 2 * (r2 * r3))          # (sqrt2 + sqrt3)^2
        self.assertEqual((r2 * r3) ** 2, 6)

    def test_roots_with_algebraic_coefficients(self):
        from fractions import Fraction as Fr
        from fmsl.algebraic import A, AlgField, roots_alg
        r2 = A.root(AlgField([-2, 0, 1], 1, 2))
        rts = roots_alg([Fr(-1), -r2, Fr(1)])                 # t^2 - sqrt2 t - 1
        self.assertEqual(len(rts), 2)
        for t in rts:
            self.assertEqual(t * t - r2 * t - 1, 0)           # exact
        self.assertAlmostEqual(float(rts[1]), (2 ** .5 + 6 ** .5) / 2, places=12)

