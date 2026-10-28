import unittest
from fractions import Fraction as Fr
from fmsl.rate import rate, rate_1d


class TestRate(unittest.TestCase):
    def test_simple(self):
        self.assertEqual(rate("x^2*y/(x^2+y^2)")["exponent"], 1)
        self.assertEqual(rate("x^4/(x^2+y^2)")["exponent"], 2)
        self.assertEqual(rate("(x^2*y^2)/(x^2+y^4)")["exponent"], 2)       # along x = c y^2
        self.assertEqual(rate("(2*x^4+2*y^2+x^3*y)/(x^4+y^2)")["exponent"], 1)

    def test_attained_after_zooming(self):
        # x^m/((y - x^k)^2 + x^(2n)): theta = m - 2n, attained along y = x^k + c x^n
        for f, th in [("x^7/((y-x^2)^2+x^6)", 1), ("x^9/((y-x^2-x^3)^2+x^8)", 1), ("y^13/((x-y^3)^2+y^10)", 3)]:
            r = rate(f)
            self.assertEqual(r["exponent"], th, f); self.assertTrue(r["proven"])

    def test_no_limit(self):
        self.assertIsNone(rate("x^3*y/(x^6+y^2)")["exponent"])

    def test_one_variable(self):
        self.assertEqual(rate_1d("sin(x)/x")["exponent"], 2)
        self.assertEqual(rate_1d("(exp(x)-1)/x")["exponent"], 1)


if __name__ == "__main__":
    unittest.main()
